use std::sync::Arc;

use giac_core::{
    eval, expr_to_poly, poly_to_expr, Context, EvalError, Expr, ExprArc, Ident, RelOp,
};
use giac_linalg::eval_linsolve;
use giac_poly::{roots, PolyError, Var};

use crate::rootof::quadratic_rootof_roots;

/// `solve(equation, var)` or `solve([equations], [vars])`.
pub fn eval_solve(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("solve"));
    }
    if matches!(args[0].as_ref(), Expr::List(_) | Expr::Seq(_)) {
        return eval_linsolve(&args[0], &args[1], ctx);
    }
    let var = ident_from_expr(&args[1])?;
    let poly = equation_to_poly(args[0].as_ref(), ctx)?;
    let v = Var::from(var.as_str());
    let items: Vec<ExprArc> = match roots(&poly, &v) {
        Ok(rs) => rs.into_iter().map(|p| poly_to_expr(&p)).collect(),
        Err(PolyError::NotImplemented(_)) => quadratic_rootof_roots(&poly, &v)?,
        Err(e) => return Err(poly_err(e)),
    };
    eval(Arc::new(Expr::List(items)).as_ref(), ctx)
}

fn equation_to_poly(eq: &Expr, ctx: &Context) -> Result<giac_poly::Poly, EvalError> {
    let diff = match eq {
        Expr::Relation(RelOp::Eq, lhs, rhs) => Expr::add(vec![
            Arc::clone(lhs),
            Expr::mul(vec![Expr::int(-1), Arc::clone(rhs)]),
        ]),
        other => Arc::new(other.clone()),
    };
    let diff = eval(diff.as_ref(), ctx)?;
    expr_to_poly(diff.as_ref())
}

fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("variable name expected")),
    }
}

fn poly_err(e: giac_poly::PolyError) -> EvalError {
    match e {
        giac_poly::PolyError::NotImplemented(s) => EvalError::NotImplemented(s),
        giac_poly::PolyError::TypeError(s) => EvalError::TypeError(s),
        giac_poly::PolyError::DivisionByZero => EvalError::TypeError("division by zero"),
    }
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn solve_quadratic_double_root() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Solve,
            vec![
                Arc::new(Expr::Relation(
                    RelOp::Eq,
                    Expr::add(vec![
                        Expr::pow(Expr::sym("x"), Expr::int(2)),
                        Expr::mul(vec![Expr::int(-2), Expr::sym("x")]),
                        Expr::int(1),
                    ]),
                    Expr::int(0),
                )),
                Expr::sym("x"),
            ],
        );
        let r = giac_core::eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "[1]");
    }

    #[test]
    fn solve_linear_system_delegates_to_linsolve() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Solve,
            vec![
                Arc::new(Expr::List(vec![
                    Arc::new(Expr::Relation(
                        RelOp::Eq,
                        Expr::add(vec![Expr::sym("x"), Expr::sym("y")]),
                        Expr::int(3),
                    )),
                    Arc::new(Expr::Relation(
                        RelOp::Eq,
                        Expr::add(vec![
                            Expr::sym("x"),
                            Expr::mul(vec![Expr::int(-1), Expr::sym("y")]),
                        ]),
                        Expr::int(1),
                    )),
                ])),
                Arc::new(Expr::List(vec![Expr::sym("x"), Expr::sym("y")])),
            ],
        );
        let r = giac_core::eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("2") && s.contains("1"), "got {s}");
    }
}
