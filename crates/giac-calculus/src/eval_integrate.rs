use std::collections::HashMap;

use giac_core::{eval, eval_subst_map, Context, EvalError, Expr, ExprArc};

use crate::integrate::integrate;

pub fn eval_integrate(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 && args.len() != 4 {
        return Err(EvalError::TooFewArgs("integrate"));
    }
    let var = match args[1].as_ref() {
        Expr::Symbol(id) => id.clone(),
        _ => return Err(EvalError::TypeError("integration variable")),
    };
    let antideriv = integrate(&args[0], &var)?;
    if args.len() == 2 {
        return Ok(antideriv);
    }
    let lo_bound = eval(&args[2], ctx)?;
    let hi_bound = eval(&args[3], ctx)?;
    let mut hi_sub = HashMap::new();
    hi_sub.insert(var.clone(), hi_bound);
    let hi = eval(eval_subst_map(&antideriv, &hi_sub)?.as_ref(), ctx)?;
    let mut lo_sub = HashMap::new();
    lo_sub.insert(var, lo_bound);
    let lo = eval(eval_subst_map(&antideriv, &lo_sub)?.as_ref(), ctx)?;
    eval(
        Expr::add(vec![hi, Expr::mul(vec![Expr::int(-1), lo])]).as_ref(),
        ctx,
    )
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use crate::plugin::xcas_default;

    #[test]
    fn eval_integrate_definite_one() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![Expr::int(1), Expr::sym("x"), Expr::int(-1), Expr::int(1)],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2");
    }

    #[test]
    fn eval_integrate_inv_x() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![Expr::pow(Expr::sym("x"), Expr::int(-1)), Expr::sym("x")],
        );
        assert!(eval(e.as_ref(), &ctx).is_ok());
    }
}
