//! `rootof` branch helpers — **tests / reference only** (S4: production uses `poly_algext_roots`).
use giac_core::{quadratic_rootof_branches, EvalError, ExprArc};
use giac_poly::{Poly, Var};

/// Two `rootof` branches for quadratic irrational roots of `poly` in `var`.
/// **Stable (bounded)** — two rootof branches for quadratic
pub fn quadratic_rootof_roots(poly: &Poly, var: &Var) -> Result<Vec<ExprArc>, EvalError> {
    quadratic_rootof_branches(poly, var)
}

#[cfg(test)]
mod tests {
    //! Test tiers: **A** / **B** — `.doc/test-writing-spec.md` · audit §2

    use std::sync::Arc;

    use giac_core::{eval, contains_algext, Ident, Context, FuncKind, RelOp, Expr};
    use giac_poly::{roots, Poly, Var};

    use super::*;
    use crate::plugin::xcas_default;
    use crate::solve_poly::solve_univariate_over_q;
    use crate::test_verify::test_verify::{
        assert_equation_solutions, assert_is_algext_or_rootof, assert_roots_zero_poly,
        list_items,
    };

    fn x() -> Poly {
        Poly::var("t")
    }

    fn t_sq_minus(n: i64) -> Poly {
        x().pow(2).sub(&Poly::constant(num_rational::Ratio::from_integer(
            num_bigint::BigInt::from(n),
        )))
    }

    // **B** — quadratic_rootof_roots; assert_roots_zero_poly.
    #[test]
    fn quadratic_rootof_has_two_branches() {
        let p = t_sq_minus(2);
        let rs = quadratic_rootof_roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 2);
        for r in &rs {
            assert_is_algext_or_rootof(r);
        }
        let poly = Arc::new(Expr::add(vec![
            Expr::pow(Expr::sym("t"), Expr::int(2)),
            Expr::int(-2),
        ]));
        let ctx = xcas_default();
        assert_roots_zero_poly(&poly, &Ident::new("t"), &rs, &ctx);
    }

    // **A** — eval(Solve); assert_equation_solutions.
    #[test]
    fn solve_t_squared_minus_two_uses_rootof() {
        let ctx = xcas_default();
        let eq = Arc::new(Expr::Relation(
            RelOp::Eq,
            Expr::add(vec![
                Expr::pow(Expr::sym("t"), Expr::int(2)),
                Expr::int(-2),
            ]),
            Expr::int(0),
        ));
        let e = Expr::func(FuncKind::Solve, vec![eq.clone(), Expr::sym("t")]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(list_items(&r).len(), 2);
        for sol in list_items(&r) {
            assert_is_algext_or_rootof(sol);
        }
        assert_equation_solutions(&eq, &Ident::new("t"), &r, &ctx);
    }

    // **A** — eval(Solve); equation + AlgExt shape.
    #[test]
    fn solve_t_fourth_minus_two_uses_rootof() {
        let ctx = xcas_default();
        let eq = Arc::new(Expr::Relation(
            RelOp::Eq,
            Expr::add(vec![
                Expr::pow(Expr::sym("t"), Expr::int(4)),
                Expr::int(-2),
            ]),
            Expr::int(0),
        ));
        let e = Expr::func(FuncKind::Solve, vec![eq.clone(), Expr::sym("t")]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(list_items(&r).len(), 4);
        assert!(list_items(&r).iter().any(|s| contains_algext(s.as_ref())));
        assert_equation_solutions(&eq, &Ident::new("t"), &r, &ctx);
    }

    // **B** — biquadratic via S0 kernel (poly_algext_roots), not legacy rootof path.
    #[test]
    fn biquadratic_t_fourth_minus_two_has_four_roots() {
        let ctx = xcas_default();
        let p = x().pow(4).sub(&Poly::constant(num_rational::Ratio::from_integer(
            num_bigint::BigInt::from(2),
        )));
        let rs = solve_univariate_over_q(&p, &Var::from("t"), &ctx).unwrap();
        assert_eq!(rs.len(), 4);
    }

    // **B** — roots() over ℚ when no AlgExt needed.
    #[test]
    fn rational_quadratic_still_uses_roots() {
        let p = x().pow(2).sub(&Poly::one());
        let rs = roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 2);
    }
}
