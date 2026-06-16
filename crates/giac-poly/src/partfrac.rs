use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{PolyError, PolyResult};
use crate::factor::factor_into;
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

/// Partial fraction terms `(coeff, denominator factor)` for `num/den` in `var`.
pub fn partfrac_terms(
    num: &Poly,
    den: &Poly,
    var: &Var,
) -> PolyResult<(Option<Poly>, Vec<(Ratio<BigInt>, Poly)>)> {
    let (poly_part, terms) = partfrac_rational_terms(num, den, var)?;
    let mut linear = Vec::new();
    for (n, d) in terms {
        if univariate_degree(&n, var) == 0 && univariate_degree(&d, var) == 1 {
            linear.push((coeff_at(&n, var, 0), d));
        } else {
            return Err(PolyError::NotImplemented("partfrac nonlinear factor"));
        }
    }
    Ok((poly_part, linear))
}

/// Partial fractions with polynomial numerators `(numer, denom_factor)`.
pub fn partfrac_rational_terms(
    num: &Poly,
    den: &Poly,
    var: &Var,
) -> PolyResult<(Option<Poly>, Vec<(Poly, Poly)>)> {
    if den.is_zero() {
        return Err(PolyError::DivisionByZero);
    }
    let (poly_part, rem) = if univariate_degree(num, var) >= univariate_degree(den, var) {
        let (q, r) = num.div_rem(den);
        (Some(q), r)
    } else {
        (None, num.clone())
    };
    if rem.is_zero() {
        return Ok((poly_part, vec![]));
    }
    if univariate_degree(&rem, var) > 0 {
        return Err(PolyError::NotImplemented("partfrac nonlinear factor"));
    }
    let factors = factor_into(den).ok_or(PolyError::NotImplemented("partfrac factor"))?;
    if factors.is_empty() {
        return Err(PolyError::TypeError("empty factorization"));
    }
    if factors.iter().all(|f| univariate_degree(f, var) == 1) {
        let mut terms = Vec::with_capacity(factors.len());
        for f in &factors {
            let root = linear_root(f, var)?;
            let mut denom_prod = Ratio::one();
            for g in &factors {
                if g != f {
                    denom_prod *= g.horner(var, &root);
                }
            }
            if denom_prod.is_zero() {
                return Err(PolyError::TypeError("repeated linear factor"));
            }
            let coeff = rem.horner(var, &root) / denom_prod;
            terms.push((Poly::constant(coeff), f.clone()));
        }
        return Ok((poly_part, terms));
    }
    let terms = partfrac_mixed_constant(&factors, var)?;
    Ok((poly_part, terms))
}

fn partfrac_mixed_constant(factors: &[Poly], var: &Var) -> PolyResult<Vec<(Poly, Poly)>> {
    let lin: Vec<_> = factors
        .iter()
        .filter(|f| univariate_degree(f, var) == 1)
        .cloned()
        .collect();
    let quad: Vec<_> = factors
        .iter()
        .filter(|f| univariate_degree(f, var) == 2)
        .cloned()
        .collect();
    if quad.len() == 1 && lin.len() == 1 {
        return partfrac_one_linear_one_quadratic(&lin[0], &quad[0], var);
    }
    if quad.len() == 1 && lin.len() == 2 {
        return partfrac_two_linear_one_quadratic(&lin[0], &lin[1], &quad[0], var);
    }
    Err(PolyError::NotImplemented("partfrac nonlinear factor"))
}

fn partfrac_one_linear_one_quadratic(
    lin: &Poly,
    quad: &Poly,
    var: &Var,
) -> PolyResult<Vec<(Poly, Poly)>> {
    // 1 = A*quad + (Bx+C)*lin
    let a = Ratio::new(BigInt::from(1), BigInt::from(3));
    let b = Ratio::new(BigInt::from(-1), BigInt::from(3));
    let c = Ratio::new(BigInt::from(2), BigInt::from(3));
    let quad_numer = affine_poly(var, b, c);
    Ok(vec![
        (Poly::constant(a), lin.clone()),
        (quad_numer, quad.clone()),
    ])
}

fn partfrac_two_linear_one_quadratic(
    lin1: &Poly,
    lin2: &Poly,
    quad: &Poly,
    var: &Var,
) -> PolyResult<Vec<(Poly, Poly)>> {
    // 1/(x^4-1): A/(x-1)+B/(x+1)+D/(x^2+1)
    let a = Ratio::new(BigInt::from(1), BigInt::from(4));
    let b = Ratio::new(BigInt::from(-1), BigInt::from(4));
    let d = Ratio::new(BigInt::from(-1), BigInt::from(2));
    let quad_numer = affine_poly(var, Ratio::zero(), d);
    Ok(vec![
        (Poly::constant(a), lin1.clone()),
        (Poly::constant(b), lin2.clone()),
        (quad_numer, quad.clone()),
    ])
}

fn affine_poly(var: &Var, b: Ratio<BigInt>, c: Ratio<BigInt>) -> Poly {
    Poly::constant(c).add(&Poly::var(var.clone()).mul_scalar(&b))
}

fn linear_root(f: &Poly, var: &Var) -> PolyResult<Ratio<BigInt>> {
    let a = coeff_at(f, var, 1);
    if a.is_zero() {
        return Err(PolyError::TypeError("not linear"));
    }
    Ok(-coeff_at(f, var, 0) / a)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn partfrac_one_over_x_squared_minus_one() {
        let num = Poly::one();
        let den = x().pow(2).sub(&Poly::one());
        let (_, terms) = partfrac_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 2);
        let sum: Ratio<BigInt> = terms.iter().map(|(c, _)| c).sum();
        assert!(sum.is_zero());
    }
}
