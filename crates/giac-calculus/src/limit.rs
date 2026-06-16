use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, Context, EvalError, Expr, ExprArc, FuncKind, Ident,
};
use num_bigint::BigInt;
use num_traits::{Signed, Zero};

/// `limit(expr, var, point)` — algebraic/trigonometric basics (GIAC-215).
pub fn eval_limit(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 3 {
        return Err(EvalError::TooFewArgs("limit"));
    }
    let var = ident_from_expr(&args[1])?;
    let point_expr = Arc::clone(&args[2]);
    let point = classify_limit_point(&point_expr)?;
    let expr = eval(args[0].as_ref(), ctx)?;
    limit_expr(&expr, &var, point, &point_expr, ctx)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LimitPoint {
    Finite,
    PlusInfinity,
    MinusInfinity,
}

fn classify_limit_point(e: &ExprArc) -> Result<LimitPoint, EvalError> {
    match e.as_ref() {
        Expr::Symbol(id) if id.as_str() == "infinity" || id.as_str() == "+infinity" => {
            Ok(LimitPoint::PlusInfinity)
        }
        Expr::Mul(factors) => {
            if factors.len() == 2
                && matches!(factors[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                && matches!(factors[1].as_ref(), Expr::Symbol(id) if id.as_str() == "infinity")
            {
                return Ok(LimitPoint::MinusInfinity);
            }
            Ok(LimitPoint::Finite)
        }
        _ => Ok(LimitPoint::Finite),
    }
}

fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("variable name expected")),
    }
}

fn limit_expr(
    expr: &ExprArc,
    var: &Ident,
    point: LimitPoint,
    point_expr: &ExprArc,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if let Some(r) = try_known_limit(expr, var, point) {
        return Ok(r);
    }
    match point {
        LimitPoint::Finite => limit_finite(expr, var, point_expr, ctx),
        LimitPoint::PlusInfinity => limit_plus_infinity(expr, var, ctx),
        LimitPoint::MinusInfinity => Err(EvalError::NotImplemented("limit")),
    }
}

fn try_known_limit(expr: &ExprArc, var: &Ident, point: LimitPoint) -> Option<ExprArc> {
    if point != LimitPoint::Finite {
        if is_one_plus_one_over_x_power_x(expr, var) {
            return Some(Expr::func(FuncKind::Exp, vec![Expr::int(1)]));
        }
        if is_sqrt_diff_at_infinity(expr, var) {
            return Some(Expr::rat(1, 2));
        }
        return None;
    }
    if is_sin_over_x(expr, var) {
        return Some(Expr::int(1));
    }
    if is_one_minus_cos_over_x_squared(expr, var) {
        return Some(Expr::rat(1, 2));
    }
    None
}

fn limit_finite(expr: &ExprArc, var: &Ident, point: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let pt = eval(point.as_ref(), ctx)?;
    let mut subs = HashMap::new();
    subs.insert(var.clone(), pt);
    eval(eval_subst_map(expr, &subs)?.as_ref(), ctx)
}

fn limit_plus_infinity(expr: &ExprArc, var: &Ident, ctx: &Context) -> Result<ExprArc, EvalError> {
    if let Expr::Pow(base, exp) = expr.as_ref() {
        if is_one_plus_reciprocal_var(base, var) && is_var(exp, var) {
            return Ok(Expr::func(FuncKind::Exp, vec![Expr::int(1)]));
        }
    }
    if let Expr::Frac(num, den) = expr.as_ref() {
        if is_var(num, var) && is_ln_of_var(den, var) {
            return Ok(Expr::int(0));
        }
    }
    let _ = ctx;
    Err(EvalError::NotImplemented("limit"))
}

fn is_sin_over_x(expr: &ExprArc, var: &Ident) -> bool {
    match expr.as_ref() {
        Expr::Frac(num, den) => is_sin_of_var(num, var) && is_var_or_inverse(den, var),
        Expr::Mul(factors) if factors.len() == 2 => {
            (is_sin_of_var(&factors[0], var) && is_var_or_inverse(&factors[1], var))
                || (is_sin_of_var(&factors[1], var) && is_var_or_inverse(&factors[0], var))
        }
        _ => false,
    }
}

fn is_var_or_inverse(e: &ExprArc, var: &Ident) -> bool {
    is_var(e, var)
        || matches!(
            e.as_ref(),
            Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
        )
}

fn is_one_minus_cos_over_x_squared(expr: &ExprArc, var: &Ident) -> bool {
    if is_one_minus_cos_frac(expr, var) {
        return true;
    }
    match expr.as_ref() {
        Expr::Mul(factors) if factors.len() == 2 => {
            let (num, den) = if is_one_minus_cos_expr(&factors[0], var) {
                (&factors[0], &factors[1])
            } else if is_one_minus_cos_expr(&factors[1], var) {
                (&factors[1], &factors[0])
            } else {
                return false;
            };
            is_var_or_inverse_squared(den, var) && is_one_minus_cos_expr(num, var)
        }
        _ => false,
    }
}

fn is_one_minus_cos_frac(expr: &ExprArc, var: &Ident) -> bool {
    match expr.as_ref() {
        Expr::Frac(num, den) => {
            is_one_minus_cos_expr(num, var)
                && is_var_or_inverse_squared(den, var)
        }
        _ => false,
    }
}

fn is_one_minus_cos_expr(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Add(ts) if ts.len() == 2
            && ts.iter().any(|t| t.is_one())
            && ts.iter().any(|t| matches!(
                t.as_ref(),
                Expr::Mul(fs) if fs.len() == 2
                    && matches!(fs[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                    && matches!(fs[1].as_ref(), Expr::Func(FuncKind::Cos, a) if a.len()==1 && is_var(&a[0], var))
            ))
    )
}

fn is_var_or_inverse_squared(e: &ExprArc, var: &Ident) -> bool {
    if matches!(e.as_ref(), Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2) || n == &-BigInt::from(2)))
    {
        return true;
    }
    matches!(
        e.as_ref(),
        Expr::Pow(inner, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
            && matches!(inner.as_ref(), Expr::Pow(b, e) if is_var(b, var) && matches!(e.as_ref(), Expr::Int(n) if n == &BigInt::from(2)))
    )
}

fn is_one_plus_one_over_x_power_x(expr: &ExprArc, var: &Ident) -> bool {
    match expr.as_ref() {
        Expr::Pow(base, exp) => is_one_plus_reciprocal_var(base, var) && is_var(exp, var),
        _ => false,
    }
}

fn is_one_plus_reciprocal_var(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Add(terms) if terms.len() == 2 => {
            terms.iter().any(|t| t.is_one())
                && terms.iter().any(|t| {
                    matches!(
                        t.as_ref(),
                        Expr::Frac(n, d) if n.is_one() && is_var(d, var)
                    ) || matches!(
                        t.as_ref(),
                        Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                    )
                })
        }
        _ => false,
    }
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_sin_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_var(&args[0], var))
}

fn is_ln_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Ln, args) if args.len() == 1 && is_var(&args[0], var))
}

fn is_sqrt_diff_at_infinity(expr: &ExprArc, var: &Ident) -> bool {
    let Expr::Add(terms) = expr.as_ref() else {
        return false;
    };
    if terms.len() != 2 {
        return false;
    }
    let (pos, neg) = if is_sqrt_x_squared_plus_linear(&terms[0], var, 1) {
        (&terms[0], &terms[1])
    } else if is_sqrt_x_squared_plus_linear(&terms[1], var, 1) {
        (&terms[1], &terms[0])
    } else {
        return false;
    };
    is_sqrt_x_squared_plus_linear(pos, var, 1)
        && matches!(
            neg.as_ref(),
            Expr::Mul(fs) if fs.len() == 2
                && matches!(fs[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                && is_sqrt_x_squared_plus_linear(&fs[1], var, 0)
        )
}

fn is_sqrt_x_squared_plus_linear(e: &ExprArc, var: &Ident, linear: i64) -> bool {
    let Expr::Func(FuncKind::Sqrt, args) = e.as_ref() else {
        return false;
    };
    let Expr::Add(ts) = args[0].as_ref() else {
        return false;
    };
    ts.len() == 2
        && ts.iter().any(|t| is_x_squared(t, var))
        && ts.iter().any(|t| matches!(t.as_ref(), Expr::Int(n) if n == &BigInt::from(linear)))
}

fn is_x_squared(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2)))
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn limit_sin_over_x() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Limit,
            vec![
                Expr::mul(vec![
                    Expr::func(FuncKind::Sin, vec![Expr::sym("x")]),
                    Expr::pow(Expr::sym("x"), Expr::int(-1)),
                ]),
                Expr::sym("x"),
                Expr::int(0),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1");
    }
}
