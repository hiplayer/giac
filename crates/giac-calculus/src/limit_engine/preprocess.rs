//! Limit preprocessing before MRV (GIAC-216d / `limit_symbolic_preprocess` subset).

use std::sync::Arc;

use giac_core::{eval, Context, EvalError, Expr, ExprArc, FuncKind, Ident};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed};

use giac_simplify::ratnormal;

use crate::risch::pow2expln;

/// Structural preprocessing without `eval` (safe for nested `exp` before limit).
pub(crate) fn limit_preprocess_struct(expr: &ExprArc, var: &Ident) -> ExprArc {
    fold_exp_zero_linear(
        &factor_exp_shifted_difference(&merge_exp_quotients(&pow2expln(
            &normalize_sqrt_conjugates(&surd2pow(expr)),
            var,
        ))),
        var,
    )
}

/// `pow2expln` and light normalization before series / limit asymptotics.
pub(crate) fn limit_preprocess_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    series_preprocess(expr, var, ctx)
}

pub(crate) fn series_preprocess(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let normalized = fold_exp_zero_linear(
        &factor_exp_shifted_difference(&merge_exp_quotients(&pow2expln(
            &normalize_sqrt_conjugates(&surd2pow(expr)),
            var,
        ))),
        var,
    );
    let evaluated = eval(normalized.as_ref(), ctx)?;
    Ok(ratnormal(evaluated.as_ref(), ctx).unwrap_or(evaluated))
}

/// `sqrt(e) → e^(1/2)` (upstream `surd2pow`).
pub(crate) fn surd2pow(expr: &ExprArc) -> ExprArc {
    match expr.as_ref() {
        Expr::Func(FuncKind::Sqrt, args) if args.len() == 1 => {
            Expr::pow(surd2pow(&args[0]), Expr::rat(1, 2))
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(surd2pow).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(surd2pow).collect()),
        Expr::Pow(b, e) => Expr::pow(surd2pow(b), surd2pow(e)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(surd2pow(n), surd2pow(d))),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(surd2pow).collect()),
        _ => Arc::clone(expr),
    }
}

fn is_half_exponent(exp: &ExprArc) -> bool {
    matches!(exp.as_ref(), Expr::Rat(r) if *r == Ratio::new(1.into(), 2.into()))
        || matches!(
            exp.as_ref(),
            Expr::Frac(n, d)
                if matches!(n.as_ref(), Expr::Int(nn) if nn.is_one())
                    && matches!(d.as_ref(), Expr::Int(dd) if dd == &BigInt::from(2))
        )
}

fn sqrt_operand(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Func(FuncKind::Sqrt, a) if a.len() == 1 => Some(surd2pow(&a[0])),
        Expr::Pow(b, exp) if is_half_exponent(exp) => Some(Arc::clone(b)),
        _ => None,
    }
}

fn negated_inner(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            if matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative()) {
                return Some(Arc::clone(&fs[1]));
            }
            if matches!(fs[1].as_ref(), Expr::Int(n) if n.is_negative()) {
                return Some(Arc::clone(&fs[0]));
            }
        }
        _ => {}
    }
    None
}

/// Algebraic identity `sqrt(A)-sqrt(B) → (A-B)/(sqrt(A)+sqrt(B))` (not limit-specific shapes).
fn try_sqrt_difference_frac(pos_sqrt: &ExprArc, neg_term: &ExprArc) -> Option<ExprArc> {
    let a = sqrt_operand(pos_sqrt)?;
    let neg_inner = negated_inner(neg_term)?;
    let b = sqrt_operand(&neg_inner)?;
    Some(Arc::new(Expr::Frac(
        Expr::add(vec![a, Expr::mul(vec![Expr::int(-1), b])]),
        Expr::add(vec![Arc::clone(pos_sqrt), neg_inner]),
    )))
}

/// `sqrt(A)-X → (A-X^2)/(sqrt(A)+X)`.
fn try_sqrt_minus_x_frac(sqrt_t: &ExprArc, x_term: &ExprArc) -> Option<ExprArc> {
    let a = sqrt_operand(sqrt_t)?;
    let x_pos = negated_inner(x_term).unwrap_or_else(|| Arc::clone(x_term));
    Some(Arc::new(Expr::Frac(
        Expr::add(vec![
            a,
            Expr::mul(vec![Expr::int(-1), Expr::pow(Arc::clone(&x_pos), Expr::int(2))]),
        ]),
        Expr::add(vec![Arc::clone(sqrt_t), x_pos]),
    )))
}

fn combine_mul_with_frac(var_factor: &ExprArc, frac: &ExprArc) -> Option<ExprArc> {
    let Expr::Frac(n, d) = frac.as_ref() else {
        return None;
    };
    Some(Arc::new(Expr::Frac(
        Expr::mul(vec![Arc::clone(var_factor), Arc::clone(n)]),
        Arc::clone(d),
    )))
}

/// Conjugate-style sqrt rewrites for general series/MRV (upstream normalization, not shape tables).
pub(crate) fn normalize_sqrt_conjugates(expr: &ExprArc) -> ExprArc {
    match expr.as_ref() {
        Expr::Add(terms) if terms.len() == 2 => {
            let t0 = normalize_sqrt_conjugates(&terms[0]);
            let t1 = normalize_sqrt_conjugates(&terms[1]);
            if let Some(f) = try_sqrt_difference_frac(&t0, &t1) {
                return f;
            }
            if let Some(f) = try_sqrt_difference_frac(&t1, &t0) {
                return Expr::mul(vec![Expr::int(-1), f]);
            }
            if let Some(f) = try_sqrt_minus_x_frac(&t0, &t1) {
                return f;
            }
            if let Some(f) = try_sqrt_minus_x_frac(&t1, &t0) {
                return Expr::mul(vec![Expr::int(-1), f]);
            }
            Expr::add(vec![t0, t1])
        }
        Expr::Mul(fs) if fs.len() == 2 => {
            let f0 = normalize_sqrt_conjugates(&fs[0]);
            let f1 = normalize_sqrt_conjugates(&fs[1]);
            if let Some(r) = combine_mul_with_frac(&f0, &f1) {
                return r;
            }
            if let Some(r) = combine_mul_with_frac(&f1, &f0) {
                return r;
            }
            Expr::mul(vec![f0, f1])
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(normalize_sqrt_conjugates).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(normalize_sqrt_conjugates).collect()),
        Expr::Pow(b, e) => Expr::pow(normalize_sqrt_conjugates(b), normalize_sqrt_conjugates(e)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            normalize_sqrt_conjugates(n),
            normalize_sqrt_conjugates(d),
        )),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(normalize_sqrt_conjugates).collect()),
        _ => Arc::clone(expr),
    }
}

/// `exp(c*var)` with `c = 0` → `1` (e.g. `4^n/2^(2n)` after `merge_exp_quotients`).
fn fold_exp_zero_linear(expr: &ExprArc, var: &Ident) -> ExprArc {
    match expr.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            if crate::limit_engine::mrv::linear_coeff_in_var(&args[0], var)
                .and_then(|c| crate::limit_engine::mrv::try_const_f64(&c))
                .is_some_and(|x| x == 0.0)
            {
                return Expr::int(1);
            }
            Expr::func(
                FuncKind::Exp,
                vec![fold_exp_zero_linear(&args[0], var)],
            )
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| fold_exp_zero_linear(t, var)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| fold_exp_zero_linear(t, var)).collect()),
        Expr::Pow(b, e) => Expr::pow(fold_exp_zero_linear(b, var), fold_exp_zero_linear(e, var)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            fold_exp_zero_linear(n, var),
            fold_exp_zero_linear(d, var),
        )),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter().map(|a| fold_exp_zero_linear(a, var)).collect(),
        ),
        _ => Arc::clone(expr),
    }
}

/// `exp(a)*(exp(b+c)-exp(b)) → exp(a+b)*(exp(c)-1)` (gruntz / nested exp differences).
pub(crate) fn factor_exp_shifted_difference(expr: &ExprArc) -> ExprArc {
    match expr.as_ref() {
        Expr::Mul(fs) => {
            for (i, j) in [(0, 1), (1, 0)] {
                if fs.len() == 2 {
                    if let Some(r) = try_factor_exp_mul_diff_pair(&fs[i], &fs[j]) {
                        return r;
                    }
                }
            }
            Expr::mul(fs.iter().map(factor_exp_shifted_difference).collect())
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(factor_exp_shifted_difference).collect()),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            factor_exp_shifted_difference(n),
            factor_exp_shifted_difference(d),
        )),
        Expr::Pow(b, e) => Expr::pow(
            factor_exp_shifted_difference(b),
            factor_exp_shifted_difference(e),
        ),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter().map(factor_exp_shifted_difference).collect(),
        ),
        _ => Arc::clone(expr),
    }
}

fn try_factor_exp_mul_diff_pair(a: &ExprArc, b: &ExprArc) -> Option<ExprArc> {
    let (ea, diff_terms) = match (a.as_ref(), b.as_ref()) {
        (Expr::Func(FuncKind::Exp, args), Expr::Add(ts)) if args.len() == 1 && ts.len() == 2 => {
            (&args[0], ts)
        }
        (Expr::Add(ts), Expr::Func(FuncKind::Exp, args)) if ts.len() == 2 && args.len() == 1 => {
            (&args[0], ts)
        }
        _ => return None,
    };
    let (pos_arg, neg_arg) = {
        let (n0, a0) = unwrap_exp_arg(&diff_terms[0])?;
        let (n1, a1) = unwrap_exp_arg(&diff_terms[1])?;
        if !n0 && n1 {
            (a0, a1)
        } else if n0 && !n1 {
            (a1, a0)
        } else {
            return None;
        }
    };
    let c = expr_remove_add_term(&pos_arg, &neg_arg)?;
    Some(Expr::mul(vec![
        Expr::func(
            FuncKind::Exp,
            vec![Expr::add(vec![Arc::clone(ea), neg_arg])],
        ),
        Expr::add(vec![Expr::func(FuncKind::Exp, vec![c]), Expr::int(-1)]),
    ]))
}

fn unwrap_exp_arg(e: &ExprArc) -> Option<(bool, ExprArc)> {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => Some((false, Arc::clone(&args[0]))),
        Expr::Mul(fs) if fs.len() == 2 => {
            let (neg, inner) = if matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative()) {
                (true, &fs[1])
            } else if matches!(fs[1].as_ref(), Expr::Int(n) if n.is_negative()) {
                (true, &fs[0])
            } else {
                return None;
            };
            let Expr::Func(FuncKind::Exp, args) = inner.as_ref() else {
                return None;
            };
            if args.len() == 1 {
                Some((neg, Arc::clone(&args[0])))
            } else {
                None
            }
        }
        _ => None,
    }
}

pub(crate) fn remove_add_term(sum: &ExprArc, term: &ExprArc) -> Option<ExprArc> {
    expr_remove_add_term(sum, term)
}

fn expr_remove_add_term(sum: &ExprArc, term: &ExprArc) -> Option<ExprArc> {
    if sum == term {
        return Some(Expr::int(0));
    }
    let Expr::Add(ts) = sum.as_ref() else {
        return None;
    };
    if !ts.iter().any(|t| t == term) {
        return None;
    }
    let rest: Vec<ExprArc> = ts.iter().filter(|t| *t != term).cloned().collect();
    Some(match rest.len() {
        0 => Expr::int(0),
        1 => Arc::clone(&rest[0]),
        _ => Expr::add(rest),
    })
}

/// `exp(a)/exp(b) → exp(a-b)` and `exp(a)*exp(b)^-1` (upstream `_pow2exp` companion).
pub(crate) fn merge_exp_quotients(expr: &ExprArc) -> ExprArc {
    match expr.as_ref() {
        Expr::Frac(num, den) => {
            let n = merge_exp_quotients(num);
            let d = merge_exp_quotients(den);
            if let (Expr::Func(FuncKind::Exp, na), Expr::Func(FuncKind::Exp, da)) =
                (n.as_ref(), d.as_ref())
            {
                if na.len() == 1 && da.len() == 1 {
                    return Expr::func(
                        FuncKind::Exp,
                        vec![Expr::add(vec![
                            Arc::clone(&na[0]),
                            Expr::mul(vec![Expr::int(-1), Arc::clone(&da[0])]),
                        ])],
                    );
                }
            }
            Arc::new(Expr::Frac(n, d))
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(merge_exp_quotients).collect()),
        Expr::Mul(fs) => {
            let mut exp_terms = Vec::new();
            let mut rest = Vec::new();
            for f in fs {
                let f = merge_exp_quotients(f);
                match f.as_ref() {
                    Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
                        exp_terms.push(Arc::clone(&args[0]));
                    }
                    Expr::Pow(base, exp)
                        if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
                            && matches!(base.as_ref(), Expr::Func(FuncKind::Exp, args) if args.len() == 1) =>
                    {
                        if let Expr::Func(FuncKind::Exp, args) = base.as_ref() {
                            exp_terms.push(Expr::mul(vec![Expr::int(-1), Arc::clone(&args[0])]));
                        }
                    }
                    _ => rest.push(f),
                }
            }
            if !exp_terms.is_empty() {
                let mut sum = exp_terms.remove(0);
                for t in exp_terms {
                    sum = Expr::add(vec![sum, t]);
                }
                rest.push(Expr::func(FuncKind::Exp, vec![sum]));
            }
            if rest.is_empty() {
                Expr::int(1)
            } else if rest.len() == 1 {
                rest.pop().unwrap()
            } else {
                Expr::mul(rest)
            }
        }
        Expr::Pow(b, e) => Expr::pow(merge_exp_quotients(b), merge_exp_quotients(e)),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(merge_exp_quotients).collect()),
        _ => Arc::clone(expr),
    }
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, Expr, FuncKind};

    use super::*;

    #[test]
    fn merge_exp_quotient_frac() {
        let e = Arc::new(Expr::Frac(
            Expr::func(FuncKind::Exp, vec![Expr::sym("n")]),
            Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(2), Expr::sym("n")])]),
        ));
        let r = merge_exp_quotients(&e);
        let s = format_expr(r.as_ref());
        assert!(
            s == "exp(-n)" || s == "exp(n-2*n)" || s == "exp(n-1*(2*n))",
            "got {s}"
        );
    }

    #[test]
    fn preprocess_seven_pow_n_over_eight() {
        let ctx = crate::plugin::xcas_default();
        let var = Ident::new("n");
        let e = Arc::new(Expr::Frac(
            Expr::pow(Expr::int(7), Expr::sym("n")),
            Expr::pow(Expr::int(8), Expr::sym("n")),
        ));
        let pre = limit_preprocess_plus_infinity(&e, &var, &ctx).unwrap();
        let s = format_expr(pre.as_ref());
        assert!(s.contains("exp"), "expected exp form, got {s}");
    }

    #[test]
    fn preprocess_surd2pow_sqrt() {
        let e = Expr::func(FuncKind::Sqrt, vec![Expr::sym("x")]);
        let r = surd2pow(&e);
        let s = format_expr(r.as_ref());
        assert!(s.contains("x") && s.contains("1/2"), "got {s}");
    }

    #[test]
    fn preprocess_sqrt_conjugate_minus_var() {
        let e = Expr::add(vec![
            Expr::func(
                FuncKind::Sqrt,
                vec![Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(1)])],
            ),
            Expr::mul(vec![Expr::int(-1), Expr::sym("x")]),
        ]);
        let r = normalize_sqrt_conjugates(&e);
        let s = format_expr(r.as_ref());
        assert!(s.contains("Frac") || s.contains("/"), "expected quotient, got {s}");
    }

    #[test]
    fn preprocess_gruntz_exp_diff_factor() {
        let ctx = crate::plugin::xcas_default();
        let var = Ident::new("x");
        let stmts = giac_parse::parse_program(
            "exp(x)*(exp(1/x-exp(-x))-exp(1/x));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let factored = factor_exp_shifted_difference(e);
        let s = format_expr(factored.as_ref());
        assert!(
            s.contains("exp(-exp(-x))")
                || s.contains("exp(-exp(-1*x))")
                || s.contains("exp(-1*exp(-x))"),
            "got {s}"
        );
        let pre = limit_preprocess_struct(e, &var);
        let ps = format_expr(pre.as_ref());
        assert!(
            ps.contains("exp(-exp(-x))")
                || ps.contains("exp(-exp(-1*x))")
                || ps.contains("exp(-1*exp(-x))"),
            "struct pre: {ps}"
        );
    }
}
