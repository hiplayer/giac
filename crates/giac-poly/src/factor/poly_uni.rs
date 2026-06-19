//! Nested-ring polynomial ops (ℚ[others][var]) and sqff-over-coeff-ring factor chain.
//!
//! **Stable:** `coeff_wrt_poly`, `content_wrt`, `primitive_part_wrt`, `substitute_poly`, …
//! **Partial:** `factor_sqff_over_coeff_ring` (upstream `do_factor_hensel` slice).
//! **Ring (crate-internal):** `subresultant::{quo_exact_wrt, quo_exact_coeff, univariate_div_rem_wrt}`.
//! **Temporary (retired):** `try_factor_bivariate_eval`, `try_kronecker_bivariate` — removed (FAC-G3 covered by sparse→Hensel).

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::nested::{MainVar, UnivariateIn};
use crate::poly::Poly;
use crate::resultant::univariate_degree;

use super::ctx::{FactorSet, SqffFactorRecFn, SqffRingCtx};
use super::univariate::factor_univariate_flat;

/// **Stable** — Coefficient of `var^exp` as a polynomial in the remaining variables.
pub fn coeff_wrt_poly(p: &Poly, var: &Var, exp: u64) -> Poly {
    UnivariateIn::new(p, MainVar::new(var.clone())).coeff_at(exp)
}

/// **Stable** — Content of `p` w.r.t. `var`: gcd of all x-coefficients in ℚ[others].
pub fn content_wrt(p: &Poly, var: &Var) -> Poly {
    crate::subresultant::content_wrt_impl(p, var)
}

/// **Stable** — `p / content_wrt(p, var)` in ℚ[others][var].
pub fn primitive_part_wrt(p: &Poly, var: &Var) -> PolyResult<Poly> {
    let content = content_wrt(p, var);
    if content.is_one() {
        return Ok(p.clone());
    }
    let d = univariate_degree(p, var);
    let mut pp = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt_poly(p, var, e);
        if c.is_zero() {
            continue;
        }
        let q = crate::nested::CoeffRingPoly::new(&c)
            .exact_quo(&crate::nested::CoeffRingPoly::new(&content))
            .ok_or(PolyError::NotImplemented("poly division"))?;
        pp = pp.add(&term_with_var(&q, var, e));
    }
    Ok(pp)
}

/// **Stable** — coeff * var^exp as Poly
pub fn term_with_var(coeff: &Poly, var: &Var, exp: u64) -> Poly {
    if exp == 0 {
        return coeff.clone();
    }
    coeff.mul(&Poly::var(var.clone()).pow(exp))
}

/// **Stable** — ∂p/∂var treating coefficients in ℚ[others].
pub fn derivative_wrt(p: &Poly, var: &Var) -> Poly {
    let d = univariate_degree(p, var);
    let mut out = Poly::zero();
    for e in 1..=d {
        let c = coeff_wrt_poly(p, var, e);
        if c.is_zero() {
            continue;
        }
        let scaled = c.mul_scalar(&Ratio::from_integer(BigInt::from(e)));
        out = out.add(&term_with_var(&scaled, var, e - 1));
    }
    out
}

/// **Stable** — Square-free factorization w.r.t. `var` over ℚ[others] (Yun-style via gcd).
pub fn square_free_wrt(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    match square_free_wrt_impl(p, var) {
        Ok(f) => Ok(f),
        Err(crate::error::PolyError::NotImplemented("poly division")) => Ok(vec![(p.clone(), 1)]),
        Err(e) => Err(e),
    }
}

// **Pipeline private** — `square_free_wrt_impl`
fn square_free_wrt_impl(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    if p.is_zero() {
        return Err(PolyError::TypeError("zero polynomial"));
    }
    let mut w = p.clone();
    let mut y = derivative_wrt(&w, var);
    let g0 = w.gcd(&y);
    if !g0.is_one() {
        w = crate::subresultant::quo_exact_wrt(&w, &g0, var)?;
        y = crate::subresultant::quo_exact_wrt(&y, &g0, var)?;
    }
    y = y.sub(&derivative_wrt(&w, var));

    let mut factors = Vec::new();
    let mut k = 1usize;
    let max_k = univariate_degree(p, var) as usize + 2;
    while !y.is_zero() && k <= max_k {
        let g = w.gcd(&y);
        if !g.is_one() {
            factors.push((g.clone(), k));
            w = crate::subresultant::quo_exact_wrt(&w, &g, var)?;
        }
        y = y.sub(&derivative_wrt(&w, var));
        k += 1;
    }
    if !w.is_one() {
        factors.push((w, k));
    }
    Ok(factors)
}

/// **Stable** — Substitute `sub_var -> sub_poly` in `p`.
pub fn substitute_poly(p: &Poly, sub_var: &Var, sub_poly: &Poly) -> Poly {
    let d = univariate_degree(p, sub_var);
    let mut out = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt_poly(p, sub_var, e);
        if c.is_zero() {
            continue;
        }
        out = out.add(&c.mul(&sub_poly.pow(e)));
    }
    out
}

/// **Partial** — Factor square-free `g` in ℚ[others][var] recursively.
pub fn factor_sqff_over_coeff_ring(
    g: &Poly,
    var: &Var,
    others: &[Var],
    factor_rec: SqffFactorRecFn,
) -> PolyResult<Vec<Poly>> {
    let ctx = SqffRingCtx::new(g, var.clone(), others);
    factor_sqff_over_coeff_ring_ctx(ctx, factor_rec).map(|fs| fs.into_polys())
}

// **Pipeline private** — typed sqff factor chain
fn factor_sqff_over_coeff_ring_ctx(
    ctx: SqffRingCtx<'_>,
    factor_rec: SqffFactorRecFn,
) -> PolyResult<FactorSet> {
    let g = ctx.poly;
    let var = ctx.main.as_var();
    let others = ctx.others.as_slice();
    let dx = ctx.main_degree();
    let main = ctx.main.clone();
    if dx == 0 {
        return Ok(FactorSet::from_polys(factor_rec(ctx)?, main));
    }
    if dx == 1 {
        return Ok(FactorSet::irreducible(g.clone(), main));
    }
    if others.is_empty() {
        return Ok(FactorSet::from_polys(
            factor_univariate_flat(g, var)?,
            main,
        ));
    }
    if others.len() == 1 {
        let other = &others[0];
        let aux_refs: Vec<&Var> = vec![other];
        if super::eval::looks_irreducible_by_good_eval(g, var, &aux_refs) {
            return Ok(FactorSet::irreducible(g.clone(), main.clone()));
        }
        for (main_var, aux_var) in [(var, other), (other, var)] {
            if let Some(f) = super::sparse::try_sparse_factor(g, main_var, aux_var) {
                let set = FactorSet::from_polys(f, MainVar::new(main_var.clone()));
                if set.product_equals(g) {
                    return Ok(set);
                }
            }
            if let Some(f) = super::hensel::try_hensel_lift_bivariate(g, main_var, aux_var) {
                let set = FactorSet::from_polys(f, MainVar::new(main_var.clone()));
                if set.product_equals(g) {
                    return Ok(set);
                }
            }
        }
        for (main_var, aux_var) in [(var, other), (other, var)] {
            if let Some(f) = super::unitary::try_unitary_factor_bivariate(g, main_var, aux_var) {
                let set = FactorSet::from_polys(f, MainVar::new(main_var.clone()));
                if set.product_equals(g) {
                    return Ok(set);
                }
            }
        }
    }
    if others.len() >= 2 {
        let tower = ctx.as_poly_factor_tower();
        if tower.is_parametric_tower() {
            let set = tower.factor_sqff_chain();
            if set.product_equals(g) {
                return Ok(set);
            }
        } else {
            let aux_refs: Vec<&Var> = others.iter().collect();
            if super::eval::looks_irreducible_by_good_eval(g, var, &aux_refs) {
                return Ok(FactorSet::irreducible(g.clone(), main.clone()));
            }
        }
        let mut all_vars = vec![var.clone()];
        all_vars.extend(others.iter().cloned());
        if let Some(f) = super::unitary::try_unitary_factor(g, &all_vars) {
            let set = FactorSet::from_polys(f, main.clone());
            if set.product_equals(g) {
                return Ok(set);
            }
        }
    }
    Ok(FactorSet::irreducible(g.clone(), main))
}

trait PolyConstant {
    // **Pipeline private** — `as_constant`
    fn as_constant(&self) -> Option<Ratio<BigInt>>;
}

impl PolyConstant for Poly {
    // **Stable** — `Poly::as_constant`
    fn as_constant(&self) -> Option<Ratio<BigInt>> {
        if self.terms.len() == 1 {
            self.terms.values().next().cloned()
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    use num_traits::One;

    use super::super::multivariate::factor_multivariate_rec;
    use super::super::ctx::SqffRingCtx;

    fn sqff_rec(ctx: SqffRingCtx<'_>) -> PolyResult<Vec<Poly>> {
        factor_multivariate_rec(ctx.poly, ctx.others.as_slice())
    }

    #[test]
    fn sqff_chain_l20_via_factor_sqff_over_coeff_ring() {
        let x = Poly::var("x");
        let b = Poly::var("b");
        let c = Poly::var("c");
        let p = x
            .add(&b)
            .add(&c)
            .mul(
                &x.pow(2)
                    .sub(&x.mul(&b))
                    .sub(&x.mul(&c))
                    .add(&b.pow(2))
                    .sub(&b.mul(&c))
                    .add(&c.pow(2)),
            );
        let others = [Var::from("c"), Var::from("x")];
        let facs = factor_sqff_over_coeff_ring(&p, &Var::from("b"), &others, sqff_rec).unwrap();
        assert!(facs.len() >= 2);
        assert_eq!(facs.iter().fold(Poly::one(), |acc, f| acc.mul(f)), p);
    }

    #[test]
    fn sqff_chain_l21_ternary_linears() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        let p = x
            .sub(&y)
            .sub(&z)
            .mul(&x.sub(&y).add(&z))
            .mul(&x.add(&y).add(&z));
        let others = [Var::from("y"), Var::from("z")];
        let facs = factor_sqff_over_coeff_ring(&p, &Var::from("x"), &others, sqff_rec).unwrap();
        assert!(facs.len() >= 2);
        assert_eq!(facs.iter().fold(Poly::one(), |acc, f| acc.mul(f)), p);
    }

    #[test]
    fn content_wrt_xy_plus_y_squared() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        // xy + y² = y(x + y): content w.r.t. x is y, primitive part is x + y.
        let p = x.mul(&y).add(&y.pow(2));
        assert_eq!(content_wrt(&p, &Var::from("x")), y);
        let pp = primitive_part_wrt(&p, &Var::from("x")).unwrap();
        assert_eq!(pp, x.add(&y));
    }


    #[test]
    fn rational_content_vs_wrt_content() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = Poly::constant(Ratio::from_integer(6.into()))
            .mul(&x.mul(&y).add(&y.pow(2)));
        assert_eq!(p.content(), Ratio::from_integer(6.into()));
        assert_eq!(content_wrt(&p, &Var::from("x")), y);
        let pp = primitive_part_wrt(&p, &Var::from("x")).unwrap();
        assert_eq!(
            pp,
            x.add(&y).mul_scalar(&Ratio::from_integer(6.into()))
        );
    }
}
