//! Multivariate factorization over ℚ: pattern table → univariate → main-variable tower.
//!
//! **Pipeline:** `factor_multivariate` → `factor_multivariate_rec` → `factor_wrt_main_var`.
//! **Stable (bounded):** `factor_into_poly`, `factor_multivariate`.
//! **Pipeline private:** `factor_multivariate_rec`, `factor_wrt_main_var`.

use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;

use super::patterns::try_factor_patterns;
use super::ctx::SqffRingCtx;
use super::poly_uni::{
    content_wrt, factor_sqff_over_coeff_ring, primitive_part_wrt, square_free_wrt,
};
use super::univariate::factor_univariate_flat;
use super::util::{extract_var_power_factors, is_univariate_in, main_var, vars_in};

/// **Stable (bounded)** — factor_multivariate ok→Some
pub fn factor_into_poly(p: &Poly) -> Option<Vec<Poly>> {
    factor_multivariate(p).ok()
}

/// **Stable (bounded)** — multivariate factorization
pub fn factor_multivariate(p: &Poly) -> PolyResult<Vec<Poly>> {
    if p.is_zero() {
        return Err(PolyError::TypeError("zero polynomial"));
    }
    if p.is_one() {
        return Ok(vec![]);
    }
    let vars = vars_in(p);
    factor_multivariate_rec(p, &vars)
}

// **Pipeline private** — multivariate factor recursion (patterns→uni→main var)
pub(crate) fn factor_multivariate_rec(p: &Poly, vars: &[Var]) -> PolyResult<Vec<Poly>> {
    if let Some(f) = try_factor_patterns(p) {
        return Ok(f);
    }
    let (p, mut factors) = extract_var_power_factors(p);
    if p.is_one() {
        return Ok(factors);
    }
    let vars = if vars.is_empty() { vars_in(&p) } else { vars.to_vec() };
    let rest = match vars.len() {
        0 => vec![p],
        1 => factor_univariate_flat(&p, &vars[0])?,
        _ => {
            let x = main_var(&p, &vars);
            let others: Vec<Var> = vars.iter().filter(|v| *v != &x).cloned().collect();
            factor_wrt_main_var(&p, &x, &others)?
        }
    };
    factors.extend(rest);
    Ok(factors)
}

// **Pipeline private** — `factor_wrt_main_var`
fn factor_wrt_main_var(p: &Poly, var: &Var, others: &[Var]) -> PolyResult<Vec<Poly>> {
    if is_univariate_in(p, var) {
        return factor_univariate_flat(p, var);
    }

    let content = content_wrt(p, var);
    let mut factors = if content.is_one() {
        Vec::new()
    } else {
        factor_multivariate_rec(&content, others)?
    };

    let pp = primitive_part_wrt(p, var)?;
    let sqff = square_free_wrt(&pp, var)?;
    for (g, k) in sqff {
        let gf = factor_sqff_over_coeff_ring(&g, var, others, factor_multivariate_rec_sqff)?;
        for _ in 0..k {
            factors.extend(gf.iter().cloned());
        }
    }
    Ok(factors)
}

// **Pipeline private** — [`SqffFactorRecFn`] adapter for multivariate recursion
fn factor_multivariate_rec_sqff(ctx: SqffRingCtx<'_>) -> PolyResult<Vec<Poly>> {
    factor_multivariate_rec(ctx.poly, ctx.others.as_slice())
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    use num_traits::One;
    use crate::poly::Poly;

    #[test]
    fn factor_bivariate_product() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.mul(&y);
        let f = factor_multivariate(&p).unwrap();
        assert_eq!(f.len(), 2);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    fn factor_bivariate_mixed() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x
            .mul(&x.pow(2).add(&x.mul(&y)).add(&y.pow(2)))
            .mul(
                &Poly::constant(Ratio::from_integer(3.into()))
                    .mul(&x)
                    .sub(&Poly::constant(Ratio::from_integer(2.into())).mul(&y)),
            );
        if let Ok(f) = factor_multivariate(&p) {
            assert!(f.len() >= 2);
            let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
            assert_eq!(prod, p);
        }
    }

    #[test]
    fn factor_three_shifted_linears() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x
            .sub(&y)
            .add(&Poly::one())
            .mul(&x.sub(&y))
            .mul(&x.sub(&y).sub(&Poly::one()));
        let f = factor_multivariate(&p).unwrap();
        assert!(f.len() >= 3);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    fn factor_var_power_times_linear() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.pow(2).mul(&y.pow(2)).mul(&x.sub(&Poly::one()));
        let f = factor_multivariate(&p).unwrap();
        assert_eq!(f.len(), 5);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
        assert_eq!(
            f.iter().filter(|q| **q == x).count(),
            2
        );
        assert_eq!(
            f.iter().filter(|q| **q == y).count(),
            2
        );
    }

    #[test]
    fn factor_repeated_linear_pairs() {
        let x = Poly::var("x");
        let p = x.sub(&Poly::one()).pow(2).mul(&x.add(&Poly::constant(Ratio::from_integer(2.into()))).pow(2));
        let f = factor_multivariate(&p).unwrap();
        assert_eq!(f.len(), 4);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }
}
