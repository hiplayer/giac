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
    let factors = factor_into(den).ok_or(PolyError::NotImplemented("partfrac factor"))?;
    if factors.is_empty() {
        return Err(PolyError::TypeError("empty factorization"));
    }
    for f in &factors {
        if univariate_degree(f, var) != 1 {
            return Err(PolyError::NotImplemented("partfrac nonlinear factor"));
        }
    }
    let mut terms = Vec::with_capacity(factors.len());
    for f in &factors {
        let root = linear_root(f, var)?;
        let mut denom_prod = Ratio::one();
        for g in &factors {
            if g != f {
                let r2 = linear_root(g, var)?;
                denom_prod *= root.clone() - r2;
            }
        }
        if denom_prod.is_zero() {
            return Err(PolyError::TypeError("repeated linear factor"));
        }
        let coeff = rem.horner(var, &root) / denom_prod;
        terms.push((coeff, f.clone()));
    }
    Ok((poly_part, terms))
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
