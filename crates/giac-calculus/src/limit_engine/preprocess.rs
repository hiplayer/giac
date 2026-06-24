//! Limit preprocessing before MRV (GIAC-216d / `limit_symbolic_preprocess` subset).
//!
//! **API 分层：** [`giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md)
//! **专项契约：** [`limit-engine-expr-api.md`](../../../../../.doc/limit-engine-expr-api.md)、[`exp-diff-expr-api.md`](../../../../../.doc/exp-diff-expr-api.md)
//!
//! | 层级 | 内容 |
//! |------|------|
//! | **Stable** | `factor_exp_shifted_difference`（→ `canonical_exp_diff`） |
//! | **Pipeline** | `limit_preprocess_*`、`merge_exp_quotients`、`surd2pow` 等 |
//! | **Pipeline private** | `sqrt` 共轭与分式辅助 |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{Context, EvalError, Expr, ExprArc, FuncKind, Ident};


use num_traits::{Signed};

use giac_simplify::{normal, ratnormal};

use crate::risch::pow2expln;
use super::util::is_half_exponent;

use super::exp_diff::{
    algebraize_exp_vanishing_products, balance_exp_arguments_frac_var, first_order_exp_vanishing_epsilon,
    canonical_exp_diff, is_exp_minus_one_factor, simplify_exp_argument_adds,
};
use super::mrv_series_lead::{normalize_inverse_sums, unify_top_quotient};

/// **Pipeline** — MRV 路径预处理（无 Gruntz ε 改写）
pub(crate) fn limit_preprocess_mrv(expr: &ExprArc, var: &Ident) -> ExprArc {
    let folded = canonical_exp_diff(expr);
    let normalized = fold_exp_zero_linear(
        &merge_exp_quotients(&pow2expln(
            &normalize_sqrt_conjugates(&surd2pow(&folded)),
            var,
        )),
        var,
    );
    let refolded = canonical_exp_diff(&normalized);
    simplify_exp_argument_adds(&balance_exp_arguments_frac_var(&refolded, var))
}

/// **Pipeline** — 结构预处理（含 ε 展开；勿用于 MRV peel 前）
pub(crate) fn limit_preprocess_struct(expr: &ExprArc, var: &Ident) -> ExprArc {
    let expr = unify_top_quotient(&normalize_inverse_sums(expr));
    // Fold `exp(A)-exp(B)` while `1/x` is still `Frac(1,x)`; `pow2expln` rewrites to `x^-1`
    // and breaks shared-subterm detection in `detect_exp_difference_add`.
    let folded = canonical_exp_diff(&expr);
    let normalized = fold_exp_zero_linear(
        &merge_exp_quotients(&pow2expln(
            &normalize_sqrt_conjugates(&surd2pow(&folded)),
            var,
        )),
        var,
    );
    // `pow2expln` / `merge_exp_quotients` can expand `exp(A)-exp(B)`; refold before ε rewrite.
    let refolded = canonical_exp_diff(&normalized);
    let epsilon_expanded = first_order_exp_vanishing_epsilon(&refolded, var);
    let simplified = simplify_exp_argument_adds(&epsilon_expanded);
    let algebraized = algebraize_exp_vanishing_products(&simplified, var);
    simplify_exp_argument_adds(&balance_exp_arguments_frac_var(&algebraized, var))
}

/// **Pipeline** — `+∞` 预处理编排
pub(crate) fn limit_preprocess_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    series_preprocess(expr, var, ctx)
}

/// **Pipeline** — 级数路径预处理
pub(crate) fn series_preprocess(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let folded = canonical_exp_diff(expr);
    let normalized = fold_exp_zero_linear(
        &merge_exp_quotients(&pow2expln(
            &normalize_sqrt_conjugates(&surd2pow(&folded)),
            var,
        )),
        var,
    );
    let rat = ratnormal(normalized.as_ref(), ctx).unwrap_or_else(|_| Arc::clone(&normalized));
    Ok(normal(rat.as_ref(), ctx).unwrap_or(rat))
}

/// **Pipeline** — `sqrt` → 有理指数
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

// **Pipeline private** — sqrt operand
fn sqrt_operand(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Func(FuncKind::Sqrt, a) if a.len() == 1 => Some(surd2pow(&a[0])),
        Expr::Pow(b, exp) if is_half_exponent(exp) => Some(Arc::clone(b)),
        _ => None,
    }
}

// **Pipeline private** — negated inner
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

// **Pipeline private** — try sqrt difference frac
fn try_sqrt_difference_frac(pos_sqrt: &ExprArc, neg_term: &ExprArc) -> Option<ExprArc> {
    let a = sqrt_operand(pos_sqrt)?;
    let neg_inner = negated_inner(neg_term)?;
    let b = sqrt_operand(&neg_inner)?;
    Some(Arc::new(Expr::Frac(
        Expr::add(vec![a, Expr::mul(vec![Expr::int(-1), b])]),
        Expr::add(vec![Arc::clone(pos_sqrt), neg_inner]),
    )))
}

// **Pipeline private** — try sqrt minus x frac
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

// **Pipeline private** — combine mul with frac
fn combine_mul_with_frac(var_factor: &ExprArc, frac: &ExprArc) -> Option<ExprArc> {
    let Expr::Frac(n, d) = frac.as_ref() else {
        return None;
    };
    Some(Arc::new(Expr::Frac(
        Expr::mul(vec![Arc::clone(var_factor), Arc::clone(n)]),
        Arc::clone(d),
    )))
}

/// **Pipeline** — 共轭有理化 `sqrt` 差分
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

/// **Pipeline** — `exp(0·x+…)` 线性折叠
pub(crate) fn fold_exp_zero_linear(expr: &ExprArc, var: &Ident) -> ExprArc {
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

/// **Stable** — 转发 `canonical_exp_diff`
pub(crate) fn factor_exp_shifted_difference(expr: &ExprArc) -> ExprArc {
    canonical_exp_diff(expr)
}

/// **Pipeline** — 合并 exp 商式
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
            let merged_fs: Vec<ExprArc> = fs.iter().map(merge_exp_quotients).collect();
            if merged_fs.iter().any(is_exp_minus_one_factor) {
                return if merged_fs.len() == 1 {
                    merged_fs.into_iter().next().unwrap()
                } else {
                    Expr::mul(merged_fs)
                };
            }
            let mut exp_terms = Vec::new();
            let mut rest = Vec::new();
            for f in merged_fs {
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
    //! Test tiers: **A′** / **B** / **C** — `.doc/test-writing-spec.md`
    //! Audit: `.doc/issues/GIAC-expr-api-test-audit.md` §4

    use giac_core::{format_expr, Expr, FuncKind, Ident};
    use giac_simplify::assert_equiv;

    use super::*;

    // **C** — merge_exp_quotients; multi-alternative display (TODO: B assert_equiv → exp(-n)).
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

    // **A′** — limit_preprocess_plus_infinity(7^n/8^n) → limit == 0.
    #[test]
    fn preprocess_seven_pow_n_over_eight() {
        let ctx = crate::plugin::xcas_default();
        let var = Ident::new("n");
        let e = Arc::new(Expr::Frac(
            Expr::pow(Expr::int(7), Expr::sym("n")),
            Expr::pow(Expr::int(8), Expr::sym("n")),
        ));
        let pre = limit_preprocess_plus_infinity(&e, &var, &ctx).unwrap();
        let r = crate::limit_engine::asymptotic::limit_at_plus_infinity(&pre, &var, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "0");
    }

    // **B** — surd2pow; assert_equiv to x^(1/2).
    #[test]
    fn preprocess_surd2pow_sqrt() {
        let ctx = crate::plugin::xcas_default();
        let e = Expr::func(FuncKind::Sqrt, vec![Expr::sym("x")]);
        let r = surd2pow(&e);
        let expected = Expr::pow(Expr::sym("x"), Expr::rat(1, 2));
        assert!(
            assert_equiv(r.as_ref(), expected.as_ref(), &ctx).unwrap(),
            "got {}",
            format_expr(r.as_ref())
        );
    }

    // **B** — normalize_sqrt_conjugates; output is Frac.
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
        assert!(
            matches!(r.as_ref(), Expr::Frac(_, _)),
            "expected quotient, got {}",
            format_expr(r.as_ref())
        );
    }

    // **A′** — nested Gruntz preprocess chain → limit == 1 (includes B checks on fold).
    #[test]
    fn preprocess_nested_gruntz_exp_diff_factor() {
        use crate::limit_engine::asymptotic::limit_at_plus_infinity;
        use crate::limit_engine::exp_diff::{canonical_exp_diff, match_exp_times_exp_minus_one};
        let ctx = crate::plugin::xcas_default();
        let var = Ident::new("x");
        let stmts = giac_parse::parse_program(
            "exp(x)*(exp(1/x+exp(-x)+exp(-(x^2)))-exp(1/x-exp(-exp(x))));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let folded = canonical_exp_diff(e);
        assert!(
            match_exp_times_exp_minus_one(&folded).is_some(),
            "fold step: {}",
            format_expr(folded.as_ref())
        );
        let after_pow = merge_exp_quotients(&pow2expln(
            &normalize_sqrt_conjugates(&surd2pow(&folded)),
            &var,
        ));
        let ps = format_expr(after_pow.as_ref());
        assert!(
            match_exp_times_exp_minus_one(&after_pow).is_some(),
            "expected exp(L)*(exp(S)-1) after pow2expln, got {ps} (fold was {})",
            format_expr(folded.as_ref())
        );
        let pre = limit_preprocess_struct(e, &var);
        let r = limit_at_plus_infinity(&pre, &var, &ctx)
            .unwrap_or_else(|e| panic!("limit failed for {}: {e:?}", format_expr(pre.as_ref())));
        assert_eq!(format_expr(r.as_ref()), "1");
    }

    // **A′** — factor_exp_shifted_difference + struct pre → limit == -1.
    #[test]
    fn preprocess_gruntz_exp_diff_factor() {
        use crate::limit_engine::exp_diff::match_exp_times_exp_minus_one;
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
        let (_, epsilon) = match_exp_times_exp_minus_one(&factored)
            .unwrap_or_else(|| panic!("expected exp(L)*(exp(S)-1), got {}", format_expr(factored.as_ref())));
        let neg_exp_neg_x = Expr::mul(vec![
            Expr::int(-1),
            Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]),
        ]);
        assert!(
            assert_equiv(epsilon.as_ref(), neg_exp_neg_x.as_ref(), &ctx).unwrap(),
            "factored epsilon got {}",
            format_expr(epsilon.as_ref())
        );
        let pre = limit_preprocess_struct(e, &var);
        let r = crate::limit_engine::asymptotic::limit_at_plus_infinity(&pre, &var, &ctx)
            .unwrap_or_else(|e| panic!("limit failed for {}: {e:?}", format_expr(pre.as_ref())));
        assert_eq!(format_expr(r.as_ref()), "-1");
    }
}
