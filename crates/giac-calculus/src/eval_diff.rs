use giac_core::{Context, EvalError, Expr, ExprArc};
use std::sync::Arc;

use crate::diff::diff;

pub fn eval_diff(args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("diff"));
    }
    match args[1].as_ref() {
        Expr::Symbol(id) => diff(&args[0], id),
        Expr::Seq(vars) | Expr::List(vars) => {
            let mut parts = Vec::with_capacity(vars.len());
            for v in vars {
                let id = match v.as_ref() {
                    Expr::Symbol(id) => id.clone(),
                    _ => return Err(EvalError::TypeError("differentiation variable")),
                };
                parts.push(diff(&args[0], &id)?);
            }
            Ok(Arc::new(Expr::List(parts)))
        }
        _ => Err(EvalError::TypeError("differentiation variable")),
    }
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use crate::plugin::xcas_default;

    use super::*;

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

    #[test]
    fn eval_derive_multivariate() {
        let ctx = xcas_default();
        let f = Expr::add(vec![
            Expr::mul(vec![
                Expr::int(2),
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::sym("y"),
            ]),
            Expr::mul(vec![Expr::int(-1), Expr::sym("x"), Expr::pow(Expr::sym("z"), Expr::int(3))]),
        ]);
        let e = Expr::func(
            FuncKind::Derive,
            vec![
                f,
                Arc::new(Expr::Seq(vec![
                    Expr::sym("x"),
                    Expr::sym("y"),
                    Expr::sym("z"),
                ])),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.starts_with('[') && s.ends_with(']'), "expected List, got {s}");
        assert!(s.contains("4*y") || s.contains("2*y"));
        assert!(s.contains("2*x^2"));
        assert!(s.contains("-z^3") || s.contains("3*z^2"));
    }

    #[test]
    fn eval_diff_too_few_args() {
        let ctx = xcas_default();
        let e = Expr::func(FuncKind::Diff, vec![Expr::sym("x")]);
        assert!(eval(e.as_ref(), &ctx).is_err());
    }

    #[test]
    fn eval_diff_list_multivariate() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Diff,
            vec![
                Expr::mul(vec![Expr::sym("x"), Expr::sym("y")]),
                Arc::new(Expr::List(vec![Expr::sym("x"), Expr::sym("y")])),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.starts_with('['), "got {s}");
    }

    #[test]
    fn eval_diff_bad_variable() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Diff,
            vec![Expr::sym("x"), Expr::int(1)],
        );
        assert!(eval(e.as_ref(), &ctx).is_err());
    }
}
