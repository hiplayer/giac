use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, bigint_to_i64, Context, EvalError, Expr, ExprArc, Ident,
    RelOp,
};
use num_traits::Zero;

use crate::diff::diff;

/// `series(f,var,center,order)` / `taylor(f,var,center,order)` (GIAC-216).
pub fn eval_series(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() < 3 {
        return Err(EvalError::TooFewArgs("series"));
    }
    let f = eval(args[0].as_ref(), ctx)?;
    let (var, center, order) = parse_series_location(&args[1..], ctx)?;
    taylor_series(&f, &var, &center, order, ctx)
}

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

fn parse_var_center(e: &ExprArc) -> Result<(Ident, ExprArc), EvalError> {
    match e.as_ref() {
        Expr::Relation(RelOp::Eq, lhs, rhs) => {
            Ok((ident_from_expr(lhs)?, Arc::clone(rhs)))
        }
        Expr::Symbol(id) => Ok((id.clone(), Expr::int(0))),
        _ => Err(EvalError::TypeError("series center expected")),
    }
}

fn series_order_arg(e: &ExprArc, ctx: &Context) -> Result<usize, EvalError> {
    let ev = eval(e.as_ref(), ctx)?;
    match ev.as_ref() {
        Expr::Int(n) => usize::try_from(bigint_to_i64(n)?).map_err(|_| EvalError::TypeError("series order")),
        _ => Err(EvalError::TypeError("series order")),
    }
}

fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("variable name expected")),
    }
}

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
    eval(sum.as_ref(), ctx)
}

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

fn var_to_expr(var: &Ident) -> ExprArc {
    Expr::sym(var.as_str())
}

fn eval_at(expr: &ExprArc, var: &Ident, center: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let mut subs = HashMap::new();
    subs.insert(var.clone(), Arc::clone(center));
    let subbed = eval_subst_map(expr, &subs)?;
    eval(subbed.as_ref(), ctx)
}

fn is_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

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
}
