//! Shared univariate solve kernel over ℚ (S0).
//!
//! **Pipeline private** — `solve_univariate_over_q`, `solve_irreducible_factor`.

use std::sync::Arc;

use giac_core::{
    factor_into_algext, format_expr, poly_algext_from_poly, poly_algext_roots_for_ctx, poly_to_expr,
    rootof_from_minpoly, univariate_poly_to_poly1_expr, Context, EvalError, Expr, ExprArc,
    PolyAlgExt,
};
use giac_poly::{factor_univariate_pairs, roots, univariate_degree, Poly, Var};

/// All roots of `poly` w.r.t. `var` over ℚ (with algebraic / rootof branches).
// **Pipeline private** — sqff × factor → per-factor roots
pub(crate) fn solve_univariate_over_q(
    poly: &Poly,
    var: &Var,
    ctx: &Context,
) -> Result<Vec<ExprArc>, EvalError> {
    if poly.is_zero() {
        return Err(EvalError::TypeError("zero polynomial"));
    }
    let deg = univariate_degree(poly, var);
    if deg == 0 {
        return Err(EvalError::TypeError("constant polynomial"));
    }

    let mut roots = Vec::new();
    for (factor, mult) in factor_univariate_pairs(poly, var)? {
        if factor.is_zero() || univariate_degree(&factor, var) == 0 {
            continue;
        }
        let factor_roots = solve_factor_with_k_split(&factor, var, ctx)?;
        for _ in 0..mult {
            roots.extend(factor_roots.iter().cloned());
        }
    }
    Ok(dedup_expr_roots(roots))
}

// **Pipeline private** — optional K factor split then roots (T2-3).
fn solve_factor_with_k_split(
    factor: &Poly,
    var: &Var,
    ctx: &Context,
) -> Result<Vec<ExprArc>, EvalError> {
    // Only split quadratics irreducible over ℚ (e.g. x²+1). Cubics+ via
    // try_factor_by_roots yield K-linear factors whose roots fail eval substitution.
    if univariate_degree(factor, var) == 2 {
        let p_alg = poly_algext_from_poly(factor)?;
        if let Some(splits) = factor_into_algext(&p_alg)? {
            let mut out = Vec::new();
            for f in splits {
                out.extend(solve_irreducible_factor_algext(&f, var, ctx)?);
            }
            return Ok(out);
        }
    }
    solve_irreducible_factor(factor, var, ctx)
}

// **Pipeline private** — roots of one K[var] factor after split
fn solve_irreducible_factor_algext(
    factor: &PolyAlgExt,
    var: &Var,
    ctx: &Context,
) -> Result<Vec<ExprArc>, EvalError> {
    let d = factor.degree_wrt(var);
    match d {
        0 => Ok(vec![]),
        n if n >= 5 => Ok(vec![irreducible_rootof_branch_algext(factor, var)?]),
        _ => {
            let rs = poly_algext_roots_for_ctx(factor, var, ctx)?;
            Ok(rs
                .into_iter()
                .map(|r| Arc::new(r.as_inner().to_expr()))
                .collect())
        }
    }
}

/// Roots of one irreducible (or low-degree) factor.
// **Pipeline private** — deg≤4 algext; deg≥5 single rootof branch
pub(crate) fn solve_irreducible_factor(
    factor: &Poly,
    var: &Var,
    ctx: &Context,
) -> Result<Vec<ExprArc>, EvalError> {
    let d = univariate_degree(factor, var);
    match d {
        0 => Ok(vec![]),
        n if n <= 3 => {
            if let Ok(rs) = roots(factor, var) {
                return Ok(rs.into_iter().map(|p| poly_to_expr(&p)).collect());
            }
            let p_alg = poly_algext_from_poly(factor)?;
            let rs = poly_algext_roots_for_ctx(&p_alg, var, ctx)?;
            Ok(rs
                .into_iter()
                .map(|r| Arc::new(r.as_inner().to_expr()))
                .collect())
        }
        n if n <= 4 => {
            let p_alg = poly_algext_from_poly(factor)?;
            let rs = poly_algext_roots_for_ctx(&p_alg, var, ctx)?;
            Ok(rs
                .into_iter()
                .map(|r| Arc::new(r.as_inner().to_expr()))
                .collect())
        }
        _ => Ok(vec![irreducible_rootof_branch(factor, var)?]),
    }
}

// **Pipeline private** — S0 deg≥5 rootof branch in K[var]
fn irreducible_rootof_branch_algext(factor: &PolyAlgExt, var: &Var) -> Result<ExprArc, EvalError> {
    if factor.degree_wrt(var) < 5 {
        return Err(EvalError::TypeError("expected factor of degree >= 5"));
    }
    let minpoly = giac_core::algext_poly_to_expr(factor)?;
    rootof_from_minpoly(&[1, 0], &minpoly)
}

/// One `rootof([1,0], minpoly)` branch for degree ≥ 5 irreducible factors.
// **Pipeline private** — S0 deg≥5 rootof (not biquadratic fallback)
pub(crate) fn irreducible_rootof_branch(poly: &Poly, var: &Var) -> Result<ExprArc, EvalError> {
    if univariate_degree(poly, var) < 5 {
        return Err(EvalError::TypeError("expected factor of degree >= 5"));
    }
    let minpoly = univariate_poly_to_poly1_expr(poly, var);
    rootof_from_minpoly(&[1, 0], &minpoly)
}

// **Pipeline private** — collapse repeated roots (solve lists unique roots)
fn dedup_expr_roots(roots: Vec<ExprArc>) -> Vec<ExprArc> {
    let mut out = Vec::new();
    for r in roots {
        let key = format_expr(r.as_ref());
        if !out
            .iter()
            .any(|u| Arc::ptr_eq(u, &r) || format_expr(u.as_ref()) == key)
        {
            out.push(r);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use giac_core::{contains_algext, FuncKind, RelOp};
    use giac_poly::Poly;

    use super::*;
    use crate::plugin::xcas_default;
    use crate::test_verify::test_verify::{
        assert_equation_solutions, assert_is_algext_or_rootof, list_items,
    };

    fn poly_x5_minus_x_plus_1() -> Poly {
        Poly::var("x")
            .pow(5)
            .sub(&Poly::var("x"))
            .add(&Poly::one())
    }

    // **B** — S0: deg≥5 → one rootof branch, not biquadratic TypeError.
    #[test]
    fn solve_irreducible_deg5_one_branch() {
        let ctx = xcas_default();
        let p = poly_x5_minus_x_plus_1();
        let rs = solve_univariate_over_q(&p, &Var::from("x"), &ctx).unwrap();
        assert_eq!(rs.len(), 1);
        assert_is_algext_or_rootof(&rs[0]);
    }

    #[test]
    fn solve_x2_plus_1_has_two_roots() {
        use crate::test_verify::test_verify::assert_roots_zero_poly;
        let ctx = xcas_default();
        let x = Poly::var("x");
        let quad = x.pow(2).add(&Poly::one());
        let rs = solve_univariate_over_q(&quad, &Var::from("x"), &ctx).unwrap();
        assert_eq!(rs.len(), 2);
        let quad_expr = giac_core::Expr::add(vec![
            giac_core::Expr::pow(giac_core::Expr::sym("x"), giac_core::Expr::int(2)),
            giac_core::Expr::int(1),
        ]);
        assert_roots_zero_poly(&quad_expr, &giac_core::Ident::new("x"), &rs, &ctx);
    }

    #[test]
    fn solve_x3_minus_x_plus_1_has_three_roots() {
        use crate::test_verify::test_verify::assert_roots_zero_poly;
        let ctx = xcas_default();
        let x = Poly::var("x");
        let cubic = x.pow(3).sub(&x).add(&Poly::one());
        let rs = solve_univariate_over_q(&cubic, &Var::from("x"), &ctx).unwrap();
        assert_eq!(rs.len(), 3);
        let cubic_expr = giac_core::Expr::add(vec![
            giac_core::Expr::pow(giac_core::Expr::sym("x"), giac_core::Expr::int(3)),
            giac_core::Expr::mul(vec![giac_core::Expr::int(-1), giac_core::Expr::sym("x")]),
            giac_core::Expr::int(1),
        ]);
        assert_roots_zero_poly(&cubic_expr, &giac_core::Ident::new("x"), &rs, &ctx);
    }

    // **A** — S0: factor descent yields 5 roots (2 from x²+1 + 3 from x³−x+1).
    #[test]
    fn solve_reducible_deg5_product() {
        let ctx = xcas_default();
        let x = Poly::var("x");
        let p = x
            .pow(2)
            .add(&Poly::one())
            .mul(&x.pow(3).sub(&x).add(&Poly::one()));
        let rs = solve_univariate_over_q(&p, &Var::from("x"), &ctx).unwrap();
        assert_eq!(rs.len(), 5);
    }

    #[test]
    fn solve_reducible_deg5_product_zeros_poly() {
        use crate::test_verify::test_verify::assert_roots_zero_poly;
        let ctx = xcas_default();
        let x = Poly::var("x");
        let p = x
            .pow(2)
            .add(&Poly::one())
            .mul(&x.pow(3).sub(&x).add(&Poly::one()));
        let rs = solve_univariate_over_q(&p, &Var::from("x"), &ctx).unwrap();
        let product_expr = giac_core::Expr::mul(vec![
            giac_core::Expr::add(vec![
                giac_core::Expr::pow(giac_core::Expr::sym("x"), giac_core::Expr::int(2)),
                giac_core::Expr::int(1),
            ]),
            giac_core::Expr::add(vec![
                giac_core::Expr::pow(giac_core::Expr::sym("x"), giac_core::Expr::int(3)),
                giac_core::Expr::mul(vec![giac_core::Expr::int(-1), giac_core::Expr::sym("x")]),
                giac_core::Expr::int(1),
            ]),
        ]);
        assert_roots_zero_poly(&product_expr, &giac_core::Ident::new("x"), &rs, &ctx);
    }

    // **A** — S0: factor descends degree (x⁴−1 → 4 roots, not whole-quartic slam).
    #[test]
    fn solve_x4_minus_1_factor_descent() {
        let ctx = xcas_default();
        let x = Poly::var("x");
        let p = x.pow(4).sub(&Poly::one());
        let rs = solve_univariate_over_q(&p, &Var::from("x"), &ctx).unwrap();
        assert_eq!(rs.len(), 4);
    }

    // **A** — eval(Solve) end-to-end for irreducible quintic.
    #[test]
    fn eval_solve_irreducible_deg5() {
        let ctx = xcas_default();
        let eq = Arc::new(Expr::Relation(
            RelOp::Eq,
            giac_core::Expr::add(vec![
                giac_core::Expr::pow(giac_core::Expr::sym("x"), giac_core::Expr::int(5)),
                giac_core::Expr::mul(vec![giac_core::Expr::int(-1), giac_core::Expr::sym("x")]),
                giac_core::Expr::int(1),
            ]),
            giac_core::Expr::int(0),
        ));
        let e = giac_core::Expr::func(FuncKind::Solve, vec![eq.clone(), giac_core::Expr::sym("x")]);
        let r = giac_core::eval(e.as_ref(), &ctx).unwrap();
        let items = list_items(&r);
        assert_eq!(items.len(), 1);
        assert!(contains_algext(items[0].as_ref()) || matches!(items[0].as_ref(), Expr::Func(_, _)));
        assert_equation_solutions(&eq, &giac_core::Ident::new("x"), &r, &ctx);
    }
}
