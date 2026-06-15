use giac_core::{Context, EvalError, ExprArc};

use crate::diff::diff;

pub fn eval_diff(args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("diff"));
    }
    let var = match args[1].as_ref() {
        giac_core::Expr::Symbol(id) => id.clone(),
        _ => return Err(EvalError::TypeError("differentiation variable")),
    };
    diff(&args[0], &var)
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use crate::plugin::xcas_default;

    #[test]
    fn eval_diff_x_squared() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Diff,
            vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::sym("x")],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("2") && s.contains("x"));
    }
}
