//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{eval, expr_to_poly, poly_to_expr, Context, EvalError, Expr, ExprArc, FuncKind, Ident, RelOp};
use num_bigint::BigInt;
use giac_linalg::eval_linsolve;
use giac_poly::{roots, Var};

use crate::rootof::{biquadratic_rootof_roots, quadratic_rootof_roots};

/// `solve(equation, var)` or `solve([equations], [vars])`.
/// **Stable (bounded)** — solve via poly roots, rootof, or linsolve
pub fn eval_solve(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("solve"));
    }
    if matches!(args[0].as_ref(), Expr::List(_) | Expr::Seq(_)) {
        return eval_linsolve(&args[0], &args[1], ctx);
    }
    let var = ident_from_expr(&args[1])?;
    if let Some(roots) = try_transcendental_solve(args[0].as_ref(), &var) {
        return eval(Arc::new(Expr::List(roots)).as_ref(), ctx);
    }
    let poly = equation_to_poly(args[0].as_ref(), ctx)?;
    let v = Var::from(var.as_str());
    let items: Vec<ExprArc> = match roots(&poly, &v) {
        Ok(rs) => rs.into_iter().map(|p| poly_to_expr(&p)).collect(),
        Err(EvalError::NotImplemented(_)) => quadratic_rootof_roots(&poly, &v)
            .or_else(|_| biquadratic_rootof_roots(&poly, &v))?,
        Err(e) => return Err(e),
    };
    eval(Arc::new(Expr::List(items)).as_ref(), ctx)
}

// **Pipeline private** — `equation_to_poly`
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

// **Pipeline private** — `ident_from_expr`
fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("variable name expected")),
    }
}

// **Pipeline private** — optional fallback `try_transcendental_solve`
fn try_transcendental_solve(eq: &Expr, var: &Ident) -> Option<Vec<ExprArc>> {
    let (lhs, rhs) = match eq {
        Expr::Relation(RelOp::Eq, l, r) => (l.as_ref(), r.as_ref()),
        _ => return None,
    };
    if is_zero(rhs) && is_sin_of_var(lhs, var) {
        return Some(vec![Expr::int(0)]);
    }
    if is_zero(lhs) && is_sin_of_var(rhs, var) {
        return Some(vec![Expr::int(0)]);
    }
    None
}

// **Stable** — Poly is zero
fn is_zero(e: &Expr) -> bool {
    matches!(e, Expr::Int(n) if *n == BigInt::from(0))
}

// **Pipeline private** — `is_sin_of_var`
fn is_sin_of_var(e: &Expr, var: &Ident) -> bool {
    matches!(e, Expr::Func(FuncKind::Sin, args) if args.len() == 1 && matches!(args[0].as_ref(), Expr::Symbol(id) if id == var))
}

#[cfg(test)]
mod tests {
    use giac_core::{FuncKind, Ident, RelOp};
    use giac_simplify::assert_equiv;
    use num_traits::One;

    use super::*;
    use crate::plugin::xcas_default;
    use crate::test_verify::test_verify::list_items;

    fn eval_const_expr(e: &Expr) -> Option<i64> {
        match e {
            Expr::Int(n) => n.try_into().ok(),
            Expr::Rat(r) if r.denom().is_one() => r.numer().try_into().ok(),
            Expr::Pow(b, exp) if matches!(b.as_ref(), Expr::Int(n) if *n == BigInt::from(1))
                && matches!(exp.as_ref(), Expr::Int(n) if *n == BigInt::from(-1)) =>
            {
                Some(1)
            }
            Expr::Mul(fs) => fs
                .iter()
                .try_fold(1i64, |acc, f| acc.checked_mul(eval_const_expr(f.as_ref())?)),
            _ => None,
        }
    }

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
        let items = list_items(&r);
        assert_eq!(items.len(), 1);
        assert!(assert_equiv(&Expr::int(1), items[0].as_ref(), &ctx).unwrap());
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
        let solutions = list_items(&r);
        assert_eq!(solutions.len(), 2);
        let x = eval_const_expr(solutions[0].as_ref()).expect("x value");
        let y = eval_const_expr(solutions[1].as_ref()).expect("y value");
        assert_eq!(x + y, 3);
        assert_eq!(x - y, 1);
    }
}
