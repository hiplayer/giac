//! Shared `exp` difference normalization (GIAC-limit-exp-diff P1).
//!
//! Upstream `remove_lnexp` (`series.cc`) applies `ln_expand` / `exp_series` on coefficients;
//! the identity `exp(f) - w^{-1} = w^{-1}(exp(f + ln(w)) - 1)` is the MRV (`w`) layer.
//!
//! `fold_exp_shifted_difference` is the same ε→0 pattern on the original variable:
//! `exp(a)*(exp(b+c)-exp(b)) → exp(a+b)*(exp(c)-1)`.

use std::sync::Arc;

use giac_core::{eval, Context, Expr, ExprArc, FuncKind, Ident};
use num_bigint::BigInt;
use num_traits::{One, Signed};

use crate::expr_util::depends_on_var;

/// Bottom-up `exp` difference factorization (x-layer; used by preprocess).
pub(crate) fn fold_exp_shifted_difference(expr: &ExprArc) -> ExprArc {
    match expr.as_ref() {
        Expr::Mul(fs) => {
            if fs.len() == 2 {
                for (i, j) in [(0, 1), (1, 0)] {
                    if let Some(r) = try_rewrite_exp_mul_shifted_pair(&fs[i], &fs[j]) {
                        return fold_exp_shifted_difference(&r);
                    }
                }
            }
            Expr::mul(fs.iter().map(fold_exp_shifted_difference).collect())
        }
        Expr::Add(ts) => {
            if ts.len() == 2 {
                if let Some(r) = try_rewrite_exp_add_difference(&ts[0], &ts[1]) {
                    return fold_exp_shifted_difference(&r);
                }
                if let Some(r) = try_rewrite_exp_add_difference(&ts[1], &ts[0]) {
                    return fold_exp_shifted_difference(&Expr::mul(vec![Expr::int(-1), r]));
                }
            }
            Expr::add(ts.iter().map(fold_exp_shifted_difference).collect())
        }
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            fold_exp_shifted_difference(n),
            fold_exp_shifted_difference(d),
        )),
        Expr::Pow(b, e) => Expr::pow(
            fold_exp_shifted_difference(b),
            fold_exp_shifted_difference(e),
        ),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter().map(fold_exp_shifted_difference).collect(),
        ),
        _ => Arc::clone(expr),
    }
}

/// `scale * (exp(ε) - 1)` — common target form (P1).
pub(crate) fn exp_scale_times_exp_minus_one(scale: ExprArc, epsilon: ExprArc) -> ExprArc {
    Expr::mul(vec![
        scale,
        Expr::add(vec![Expr::func(FuncKind::Exp, vec![epsilon]), Expr::int(-1)]),
    ])
}

/// `exp(a)*(exp(b+c)-exp(b)) → exp(a+b)*(exp(c)-1)`.
pub(crate) fn try_rewrite_exp_mul_shifted_pair(a: &ExprArc, b: &ExprArc) -> Option<ExprArc> {
    let (ea, diff_terms) = match (a.as_ref(), b.as_ref()) {
        (Expr::Func(FuncKind::Exp, args), Expr::Add(ts)) if args.len() == 1 && ts.len() == 2 => {
            (&args[0], ts)
        }
        (Expr::Add(ts), Expr::Func(FuncKind::Exp, args)) if ts.len() == 2 && args.len() == 1 => {
            (&args[0], ts)
        }
        _ => return None,
    };
    let (pos_arg, neg_arg) = exp_diff_pair_from_add(diff_terms)?;
    let c = remove_add_term(&pos_arg, &neg_arg)?;
    Some(exp_scale_times_exp_minus_one(
        Expr::func(FuncKind::Exp, vec![Expr::add(vec![Arc::clone(ea), neg_arg])]),
        c,
    ))
}

/// `exp(B) - exp(B') → exp(B')*(exp(B-B')-1)` (x-layer add; same rule as mul factor).
pub(crate) fn try_rewrite_exp_add_difference(a: &ExprArc, b: &ExprArc) -> Option<ExprArc> {
    let (pos_arg, neg_arg) = match (unwrap_exp_arg(a), unwrap_exp_arg(b)) {
        (Some((false, p)), Some((true, n))) => (p, n),
        (Some((true, p)), Some((false, n))) => (n, p),
        _ => return None,
    };
    let c = remove_add_term(&pos_arg, &neg_arg)?;
    Some(exp_scale_times_exp_minus_one(
        Expr::func(FuncKind::Exp, vec![neg_arg]),
        c,
    ))
}

/// `exp(f) - scale^{-1} → scale^{-1} * (exp(f + ln(scale)) - 1)` (w-layer / MRV).
pub(crate) fn rewrite_exp_minus_scale_inv(
    f: &ExprArc,
    scale_inv: ExprArc,
    ln_scale: ExprArc,
) -> ExprArc {
    let shifted = Expr::add(vec![Arc::clone(f), ln_scale]);
    exp_scale_times_exp_minus_one(scale_inv, shifted)
}

/// After [`super::preprocess::limit_preprocess_struct`]: `exp(L)*(exp(S)-1)` @ `+∞` (P2).
pub(crate) fn try_limit_exp_times_exp_minus_one_preprocessed(
    pre: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let (outer_arg, small_shift) = match_exp_times_exp_minus_one(pre)?;
    let rest = exp_rest_after_unit_var(&outer_arg, var)?;
    if !is_neg_exp_of_neg_var(&small_shift, var) {
        return None;
    }
    limit_neg_exp_rest_at_plus_infinity(&rest, var, ctx)
}

pub(crate) fn match_exp_times_exp_minus_one(expr: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    let Expr::Mul(fs) = expr.as_ref() else {
        return None;
    };
    if fs.len() != 2 {
        return None;
    }
    for (i, j) in [(0, 1), (1, 0)] {
        if let Some(outer) = exp_func_arg(&fs[i]) {
            if let Some(shift) = exp_minus_one_inner_arg(&fs[j]) {
                return Some((outer, shift));
            }
        }
    }
    None
}

fn exp_diff_pair_from_add(diff_terms: &[ExprArc]) -> Option<(ExprArc, ExprArc)> {
    let (n0, a0) = unwrap_exp_arg(&diff_terms[0])?;
    let (n1, a1) = unwrap_exp_arg(&diff_terms[1])?;
    if !n0 && n1 {
        Some((a0, a1))
    } else if n0 && !n1 {
        Some((a1, a0))
    } else {
        None
    }
}

pub(crate) fn unwrap_exp_arg(e: &ExprArc) -> Option<(bool, ExprArc)> {
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

fn exp_func_arg(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => Some(Arc::clone(&args[0])),
        _ => None,
    }
}

fn exp_minus_one_inner_arg(e: &ExprArc) -> Option<ExprArc> {
    let Expr::Add(ts) = e.as_ref() else {
        return None;
    };
    if ts.len() != 2 {
        return None;
    }
    let (pos, neg) = if matches!(ts[1].as_ref(), Expr::Int(n) if n.is_negative()) {
        (&ts[0], &ts[1])
    } else if matches!(ts[0].as_ref(), Expr::Int(n) if n.is_negative()) {
        (&ts[1], &ts[0])
    } else {
        return None;
    };
    if !matches!(neg.as_ref(), Expr::Int(n) if n.is_negative()) {
        return None;
    }
    exp_func_arg(pos)
}

fn remove_add_term(sum: &ExprArc, term: &ExprArc) -> Option<ExprArc> {
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

fn exp_rest_after_unit_var(arg: &ExprArc, var: &Ident) -> Option<ExprArc> {
    if is_unit_var_term(arg, var) {
        return Some(Expr::int(0));
    }
    let Expr::Add(ts) = arg.as_ref() else {
        return None;
    };
    let mut rest = Vec::new();
    let mut saw_unit_var = false;
    for t in ts {
        if is_unit_var_term(t, var) {
            if saw_unit_var {
                return None;
            }
            saw_unit_var = true;
        } else {
            rest.push(Arc::clone(t));
        }
    }
    if !saw_unit_var {
        return None;
    }
    Some(match rest.len() {
        0 => Expr::int(0),
        1 => Arc::clone(&rest[0]),
        _ => Expr::add(rest),
    })
}

fn is_unit_var_term(e: &ExprArc, var: &Ident) -> bool {
    if is_var(e, var) {
        return true;
    }
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            (is_var(&fs[0], var) && is_one(&fs[1])) || (is_var(&fs[1], var) && is_one(&fs[0]))
        }
        _ => false,
    }
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_one(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n == &BigInt::from(1))
}

fn limit_neg_exp_rest_at_plus_infinity(rest: &ExprArc, var: &Ident, ctx: &Context) -> Option<ExprArc> {
    if !depends_on_var(rest, var) {
        return eval(
            Expr::mul(vec![
                Expr::int(-1),
                Expr::func(FuncKind::Exp, vec![Arc::clone(rest)]),
            ])
            .as_ref(),
            ctx,
        )
        .ok();
    }
    if vanishes_at_plus_infinity(rest, var) {
        return Some(Expr::int(-1));
    }
    None
}

fn vanishes_at_plus_infinity(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Pow(b, exp) if is_var(b, var) => {
            matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
        }
        Expr::Frac(n, d) => {
            matches!(n.as_ref(), Expr::Int(nn) if nn.is_one()) && is_var(d, var)
        }
        Expr::Mul(fs) => fs
            .iter()
            .all(|f| vanishes_at_plus_infinity(f, var) || !depends_on_var(f, var)),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => super::mrv::linear_coeff_in_var(
            &args[0],
            var,
        )
        .and_then(|c| super::mrv::try_const_f64(&c))
        .is_some_and(|x| x < 0.0),
        _ => false,
    }
}

fn is_neg_exp_of_neg_var(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            let neg = fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()));
            let exp_inner = fs.iter().find_map(|f| {
                if let Expr::Func(FuncKind::Exp, args) = f.as_ref() {
                    if args.len() == 1 {
                        return Some(&args[0]);
                    }
                }
                None
            });
            neg && exp_inner.is_some_and(|a| is_neg_var_exp(a, var))
        }
        _ => is_neg_var_exp(e, var),
    }
}

fn is_neg_var_exp(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Mul(fs) if fs.len() == 2
            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
            && fs.iter().any(|f| is_var(f, var))
    ) || matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp, args)
            if args.len() == 1
                && matches!(
                    args[0].as_ref(),
                    Expr::Mul(fs) if fs.len() == 2
                        && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
                        && fs.iter().any(|f| is_var(f, var))
                )
    )
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn fold_mul_shifted_difference_gruntz() {
        let ctx = xcas_default();
        let stmts = giac_parse::parse_program(
            "exp(x)*(exp(1/x-exp(-x))-exp(1/x));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let r = fold_exp_shifted_difference(e);
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("exp(-exp(-x))") || s.contains("exp(-exp(-1*x))"),
            "got {s}"
        );
    }

    #[test]
    fn rewrite_exp_minus_w_inv_matches_remove_lnexp() {
        use super::super::mrv_w::{mrv_ln_w_expr, mrv_w_expr};
        let w = mrv_w_expr();
        let inner = Expr::add(vec![
            Expr::mul(vec![Expr::int(-1), mrv_ln_w_expr()]),
            Expr::sym("eps"),
        ]);
        let _e = Expr::add(vec![
            Expr::func(FuncKind::Exp, vec![Arc::clone(&inner)]),
            Expr::mul(vec![Expr::int(-1), Expr::pow(Arc::clone(&w), Expr::int(-1))]),
        ]);
        let r = exp_scale_times_exp_minus_one(Expr::pow(w, Expr::int(-1)), inner);
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("exp(eps)") || s.contains("exp(-ln(_mrv_w)+eps)"),
            "expected w^-1*(exp(eps)-1) style, got {s}"
        );
        assert!(s.contains("_mrv_w"), "got {s}");
    }
}
