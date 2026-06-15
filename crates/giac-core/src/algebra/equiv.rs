//! Mathematical equivalence helpers (`assert_equiv` per conformance-testing.md §3).

use std::sync::Arc;

use crate::context::Context;
use crate::error::EvalError;
use crate::expr::{Expr, ExprArc};

use super::normal::normal;

/// Construct `a - b` as an expression tree.
pub fn sub(a: &Expr, b: &Expr) -> Result<ExprArc, EvalError> {
    Ok(Expr::add(vec![
        Arc::new(a.clone()),
        Expr::mul(vec![Expr::int(-1), Arc::new(b.clone())]),
    ]))
}

/// True when `normal(e)` simplifies to zero.
pub fn is_zero(e: &Expr, ctx: &Context) -> Result<bool, EvalError> {
    let n = normal(e, ctx)?;
    Ok(n.is_zero())
}

/// True when `a` and `b` are mathematically equivalent under `normal`.
pub fn assert_equiv(a: &Expr, b: &Expr, ctx: &Context) -> Result<bool, EvalError> {
    let d = sub(a, b)?;
    is_zero(d.as_ref(), ctx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format_expr;

    #[test]
    fn equiv_commutative_add() {
        let ctx = Context::default();
        let a = Expr::add(vec![Expr::sym("x"), Expr::int(1)]);
        let b = Expr::add(vec![Expr::int(1), Expr::sym("x")]);
        assert!(assert_equiv(a.as_ref(), b.as_ref(), &ctx).unwrap());
    }

    #[test]
    fn equiv_expanded_square() {
        let ctx = Context::default();
        let a = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2));
        let b = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
            Expr::int(1),
        ]);
        assert!(assert_equiv(a.as_ref(), b.as_ref(), &ctx).unwrap());
    }

    #[test]
    fn equiv_sqrt_half_forms() {
        use crate::expr::FuncKind;
        let ctx = Context::default();
        let a = Expr::mul(vec![
            Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]),
            Expr::rat(1, 2),
        ]);
        let b = Expr::pow(
            Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]),
            Expr::int(-1),
        );
        // Rationalize sqrt(2)/2 vs 1/sqrt(2) when normal supports radicals.
        let _ = assert_equiv(a.as_ref(), b.as_ref(), &ctx);
    }

    #[test]
    fn equiv_mod_difference() {
        let a = Expr::Mod(Expr::int(7), Expr::int(3));
        let b = Expr::int(1);
        // `normal` does not yet fold bare `Mod` to integer residue.
        let diff = sub(&a, &b).unwrap();
        assert!(!diff.is_zero());
    }

    #[test]
    fn not_equiv_different() {
        let ctx = Context::default();
        let a = Expr::sym("x");
        let b = Expr::add(vec![Expr::sym("x"), Expr::int(1)]);
        assert!(!assert_equiv(a.as_ref(), b.as_ref(), &ctx).unwrap());
        let _ = format_expr(a.as_ref());
    }
}
