//! `series` / `taylor` eval hook.
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../.doc/giac-calculus-api-stability.md) §7.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `eval_series` |
//! | **Pipeline private** | `taylor_series*`, `parse_series_*`, `series_*`, `fold_elementary`, `factorial` |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, bigint_to_i64, ident_from_expr, var_to_expr, Context, EvalError, Expr,
    ExprArc, FuncKind, Ident, RelOp,
};
use num_traits::Zero;

use crate::diff::diff;
use crate::limit_engine::{asymptotic_series_at_infinity, series_at_center};
use crate::limit_engine::preprocess::series_preprocess;

/// **Stable** — `series(f,var,center,order)` / `taylor(f,var,center,order)` (GIAC-216).
pub fn eval_series(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() < 3 {
        return Err(EvalError::TooFewArgs("series"));
    }
    let f = eval(args[0].as_ref(), ctx)?;
    let (var, center, order) = parse_series_location(&args[1..], ctx)?;
    let f = series_preprocess(&f, &var, ctx)?;
    if is_plus_infinity(&center) {
        return asymptotic_series_at_infinity(&f, &var, order, ctx);
    }
    taylor_series(&f, &var, &center, order, ctx)
}

// **Pipeline private** — detect `+infinity` center symbol.
fn is_plus_infinity(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Symbol(id) if id.as_str() == "+infinity" || id.as_str() == "infinity"
    )
}

// **Pipeline private** — parse `(var, center, order)` from 2- or 3-arg tails.
fn parse_series_location(
    args: &[ExprArc],
    ctx: &Context,
) -> Result<(Ident, ExprArc, usize), EvalError> {
    match args.len() {
        2 => {
            let (var, center) = parse_var_center(&args[0])?;
            let order = series_order_arg(&args[1], ctx)?;
            Ok((var, center, order))
        }
        3 => {
            let var = ident_from_expr(&args[0])?;
            let center = eval(args[1].as_ref(), ctx)?;
            let order = series_order_arg(&args[2], ctx)?;
            Ok((var, center, order))
        }
        _ => Err(EvalError::TooFewArgs("series")),
    }
}

// **Pipeline private** — `x=0` or `x` alone defaults center to 0.
fn parse_var_center(e: &ExprArc) -> Result<(Ident, ExprArc), EvalError> {
    match e.as_ref() {
        Expr::Relation(RelOp::Eq, lhs, rhs) => {
            Ok((ident_from_expr(lhs)?, Arc::clone(rhs)))
        }
        Expr::Symbol(id) => Ok((id.clone(), Expr::int(0))),
        _ => Err(EvalError::TypeError("series center expected")),
    }
}

// **Pipeline private** — evaluate and coerce series order to `usize`.
fn series_order_arg(e: &ExprArc, ctx: &Context) -> Result<usize, EvalError> {
    let ev = eval(e.as_ref(), ctx)?;
    match ev.as_ref() {
        Expr::Int(n) => usize::try_from(bigint_to_i64(n)?).map_err(|_| EvalError::TypeError("series order")),
        _ => Err(EvalError::TypeError("series order")),
    }
}

// **Pipeline private** — Taylor at finite center (diff path then `series_at_center`).
fn taylor_series(
    f: &ExprArc,
    var: &Ident,
    center: &ExprArc,
    order: usize,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if order == 0 {
        return Ok(Expr::int(0));
    }
    if let Ok(r) = taylor_series_diff(f, var, center, order, ctx) {
        return Ok(r);
    }
    if let Ok(r) = series_at_center(f, var, center, order, ctx) {
        return Ok(r);
    }
    Err(EvalError::NotImplemented("series"))
}

// **Pipeline private** — Taylor via repeated differentiation and substitution.
fn taylor_series_diff(
    f: &ExprArc,
    var: &Ident,
    center: &ExprArc,
    order: usize,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let mut terms = Vec::new();
    let mut fk = Arc::clone(f);
    for k in 0..order {
        let coeff = eval_at(&fk, var, center, ctx)?;
        if !is_zero(&coeff) {
            terms.push(series_term(&coeff, var, center, k)?);
        }
        if k + 1 < order {
            fk = diff(&fk, var)?;
            fk = eval(fk.as_ref(), ctx)?;
        }
    }
    if terms.is_empty() {
        return Ok(Expr::int(0));
    }
    let sum = Expr::add(terms);
    let sum = eval(sum.as_ref(), ctx)?;
    giac_simplify::normal(sum.as_ref(), ctx)
}

// **Pipeline private** — single Taylor term `(coeff/k!) * (x - center)^k`.
fn series_term(coeff: &ExprArc, var: &Ident, center: &ExprArc, k: usize) -> Result<ExprArc, EvalError> {
    let scaled = if k == 0 {
        Arc::clone(coeff)
    } else {
        let fact = factorial(k)?;
        Expr::mul(vec![coeff.clone(), Expr::rat(1, fact)])
    };
    if k == 0 {
        return Ok(scaled);
    }
    let delta = if is_zero(center) {
        var_to_expr(var)
    } else {
        Expr::add(vec![var_to_expr(var), Expr::mul(vec![Expr::int(-1), Arc::clone(center)])])
    };
    if k == 1 {
        return Ok(Expr::mul(vec![scaled, delta]));
    }
    Ok(Expr::mul(vec![scaled, Expr::pow(delta, Expr::int(i64::try_from(k).unwrap_or(0)))]))
}

// **Pipeline private** — substitute `var ↦ center` and fold elementary values.
fn eval_at(expr: &ExprArc, var: &Ident, center: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let mut subs = HashMap::new();
    subs.insert(var.clone(), Arc::clone(center));
    let subbed = eval_subst_map(expr, &subs)?;
    let v = eval(subbed.as_ref(), ctx)?;
    fold_elementary(&v, ctx)
}

// **Pipeline private** — zero-test for substitution argument.
fn is_zero_arg(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

// **Pipeline private** — fold `sin(0)`, `cos(0)`, `exp(0)` at series coefficients.
fn fold_elementary(e: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let folded = match e.as_ref() {
        Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_zero_arg(&args[0]) => {
            Expr::int(0)
        }
        Expr::Func(FuncKind::Cos, args) if args.len() == 1 && is_zero_arg(&args[0]) => {
            Expr::int(1)
        }
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 && is_zero_arg(&args[0]) => {
            Expr::int(1)
        }
        Expr::Add(ts) => {
            let mut out = Vec::with_capacity(ts.len());
            for t in ts {
                out.push(fold_elementary(t, ctx)?);
            }
            Expr::add(out)
        }
        Expr::Mul(fs) => {
            let mut out = Vec::with_capacity(fs.len());
            for f in fs {
                out.push(fold_elementary(f, ctx)?);
            }
            Expr::mul(out)
        }
        Expr::Frac(n, d) => Expr::mul(vec![
            fold_elementary(n, ctx)?,
            Expr::pow(fold_elementary(d, ctx)?, Expr::int(-1)),
        ]),
        _ => return Ok(Arc::clone(e)),
    };
    eval(folded.as_ref(), ctx)
}

// **Pipeline private** — zero-test for series coefficient.
fn is_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

// **Pipeline private** — integer factorial for Taylor denominators.
fn factorial(n: usize) -> Result<i64, EvalError> {
    let mut acc = 1_i64;
    for i in 2..=n {
        acc = acc
            .checked_mul(i as i64)
            .ok_or(EvalError::TypeError("factorial overflow"))?;
    }
    Ok(acc)
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn series_preprocess_pow2expln() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let e = Expr::pow(
            Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
            Expr::int(2),
        );
        let r = crate::limit_engine::preprocess::series_preprocess(&e, &var, &ctx).unwrap();
        let text = format_expr(r.as_ref());
        assert!(text.contains("exp"), "got {text}");
    }

    #[test]
    fn series_exp_at_zero() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Series,
            vec![
                Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
                Expr::sym("x"),
                Expr::int(0),
                Expr::int(4),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("1"), "got {s}");
        assert!(s.contains("x"), "got {s}");
    }

    #[test]
    fn taylor_sin_x_equals_zero() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Taylor,
            vec![
                Expr::func(FuncKind::Sin, vec![Expr::sym("x")]),
                Arc::new(Expr::Relation(
                    RelOp::Eq,
                    Expr::sym("x"),
                    Expr::int(0),
                )),
                Expr::int(5),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("x"), "got {s}");
    }

    #[test]
    fn series_at_plus_infinity_exp() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Series,
            vec![
                Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
                Expr::sym("x"),
                Expr::sym("+infinity"),
                Expr::int(4),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(!s.is_empty(), "got {s}");
    }

    #[test]
    fn taylor_at_one() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Taylor,
            vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::sym("x"),
                Expr::int(1),
                Expr::int(3),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("x"), "got {s}");
    }

    #[test]
    fn series_too_few_args() {
        let ctx = xcas_default();
        let e = Expr::func(FuncKind::Series, vec![Expr::sym("x")]);
        assert!(eval(e.as_ref(), &ctx).is_err());
    }

    #[test]
    fn series_order_zero() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Series,
            vec![
                Expr::sym("x"),
                Expr::sym("x"),
                Expr::int(0),
                Expr::int(0),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "0");
    }

    #[test]
    fn series_two_arg_var_symbol_defaults_center_zero() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Series,
            vec![Expr::sym("x"), Expr::sym("x"), Expr::int(2)],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("x"), "got {s}");
    }

    #[test]
    fn series_ln_one_plus_x() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Series,
            vec![
                Expr::func(
                    FuncKind::Ln,
                    vec![Expr::add(vec![Expr::int(1), Expr::sym("x")])],
                ),
                Arc::new(Expr::Relation(
                    RelOp::Eq,
                    Expr::sym("x"),
                    Expr::int(0),
                )),
                Expr::int(4),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("x"), "got {s}");
    }

    #[test]
    fn series_ck_int_65_cos_exp() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Series,
            vec![
                Expr::mul(vec![
                    Expr::func(FuncKind::Cos, vec![Expr::sym("x")]),
                    Expr::func(
                        FuncKind::Exp,
                        vec![Expr::add(vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")]), Expr::int(1)])],
                    ),
                ]),
                Expr::sym("x"),
                Expr::int(0),
                Expr::int(4),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("exp(1)"),
            "expected exp(1) in Taylor series, got {s}"
        );
    }

    #[test]
    fn series_cos_taylor_diff_fallback() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Taylor,
            vec![
                Expr::func(FuncKind::Cos, vec![Expr::sym("x")]),
                Expr::sym("x"),
                Expr::int(0),
                Expr::int(4),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("1") || s.contains("x"), "got {s}");
    }

    #[test]
    fn series_taylor_at_center_two() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Taylor,
            vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::sym("x"),
                Expr::int(2),
                Expr::int(3),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("x"), "got {s}");
    }

    #[test]
    fn eval_series_direct_bad_order() {
        let ctx = xcas_default();
        let args = vec![
            Expr::sym("x"),
            Expr::sym("x"),
            Expr::sym("y"),
            Expr::sym("z"),
        ];
        assert!(eval_series(&args, &ctx).is_err());
    }
}
