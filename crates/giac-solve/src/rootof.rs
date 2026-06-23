//! `rootof` branch helpers — **tests / reference only** (S4: production uses `poly_algext_roots`).
use giac_core::{
    algext_sqrt_branches, quadratic_rootof_branches, EvalError, Expr, ExprArc,
};
use giac_poly::{coeff_at, univariate_degree, Poly, Var};
use num_rational::Ratio;
use num_traits::{One, Zero};

/// Two `rootof` branches for quadratic irrational roots of `poly` in `var`.
/// **Stable (bounded)** — two rootof branches for quadratic
pub fn quadratic_rootof_roots(poly: &Poly, var: &Var) -> Result<Vec<ExprArc>, EvalError> {
    quadratic_rootof_branches(poly, var)
}

/// Four `rootof` branches for biquadratic `a·t⁴ + b·t² + c` (odd terms zero).
/// **Partial** — biquadratic rootof; general quartic NotImplemented
pub fn biquadratic_rootof_roots(poly: &Poly, var: &Var) -> Result<Vec<ExprArc>, EvalError> {
    if univariate_degree(poly, var) != 4 {
        return Err(EvalError::TypeError("expected quartic"));
    }
    if !coeff_at(poly, var, 3).is_zero() || !coeff_at(poly, var, 1).is_zero() {
        return Err(EvalError::NotImplemented("general quartic rootof"));
    }
    let lc = coeff_at(poly, var, 4);
    if lc.is_zero() {
        return Err(EvalError::TypeError("leading coefficient zero"));
    }
    let scale = Ratio::one() / lc;
    let b = coeff_at(poly, var, 2) * scale.clone();
    let c = coeff_at(poly, var, 0) * scale;
    let u_var = Poly::var(var.clone());
    let u_poly = u_var
        .pow(2)
        .mul_scalar(&Ratio::one())
        .add(&u_var.mul_scalar(&b))
        .add(&Poly::constant(c));
    let u_roots = quadratic_rootof_roots(&u_poly, var)?;
    let mut out = Vec::new();
    for u in u_roots {
        let u_data = match u.as_ref() {
            Expr::AlgExt(a) => (**a).clone(),
            _ => return Err(EvalError::TypeError("rootof expected")),
        };
        for t in algext_sqrt_branches(&u_data)? {
            out.push(t);
        }
    }
    if out.is_empty() {
        return Err(EvalError::NotImplemented("biquadratic rootof"));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    //! Test tiers: **A** / **B** — `.doc/test-writing-spec.md` · audit §2

    use std::sync::Arc;

    use giac_core::{eval, contains_algext, Ident, Context, FuncKind, RelOp};
    use giac_poly::{roots, Poly, Var};

    use super::*;
    use crate::plugin::xcas_default;
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

    // **B** — biquadratic_rootof_roots branch count.
    #[test]
    fn biquadratic_t_fourth_minus_two_has_four_roots() {
        let p = x().pow(4).sub(&Poly::constant(num_rational::Ratio::from_integer(
            num_bigint::BigInt::from(2),
        )));
        let rs = biquadratic_rootof_roots(&p, &Var::from("t")).unwrap();
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
