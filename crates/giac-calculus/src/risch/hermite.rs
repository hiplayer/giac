//! GIAC-227: univariate Hermite reduction for `∫ num / factor^n` (n > 1).

use giac_poly::{
    abcuv, quo, rem, univariate_degree, univariate_derivative, Poly, PolyError, PolyResult, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::One;

/// Extracted rational term `-v / factor^power` (integrates to part of the input).
#[derive(Debug, Clone, PartialEq)]
pub struct HermiteTerm {
    pub numer: Poly,
    pub factor: Poly,
    pub power: usize,
}

/// Hermite-reduce `numer / factor^mult` w.r.t. `var`.
///
/// Returns extracted terms and remaining `(remainder_num, remainder_mult)` with `remainder_mult <= 1`.
pub fn hermite_reduce(
    numer: &Poly,
    factor: &Poly,
    mult: usize,
    var: &Var,
) -> PolyResult<(Vec<HermiteTerm>, Poly, usize)> {
    if mult <= 1 {
        return Ok((vec![], numer.clone(), mult));
    }
    if univariate_degree(factor, var) == 0 {
        return Err(PolyError::TypeError("constant factor"));
    }

    let mut extracted = Vec::new();
    let mut a = numer.clone();
    let g = factor.clone();
    let gp = univariate_derivative(&g, var);
    let mut n = mult;

    while n > 1 {
        let (_u, v_raw) = abcuv(&g, &gp, &a)?;
        let v = rem(&v_raw, &g)?;
        let u = quo(&(a.sub(&v.mul(&gp))), &g)?;
        extracted.push(HermiteTerm {
            numer: v.clone(),
            factor: g.clone(),
            power: n - 1,
        });
        // ∫ P/Q^n = -V/((n-1)Q^{n-1}) + ∫ (U + V'/((n-1)))/Q^{n-1}
        let vp = univariate_derivative(&v, var);
        let scale = Ratio::from_integer(BigInt::from((n - 1) as i64));
        a = u.add(&vp.mul_scalar(&(Ratio::one() / scale)));
        n -= 1;
    }

    Ok((extracted, a, n))
}

#[cfg(test)]
mod tests {
    use giac_poly::coeff_at;
    use num_traits::{One, Zero};

    use super::*;

    fn x_var() -> Var {
        Var::from("x")
    }

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn hermite_x_over_x_squared_plus_one_squared_term_shape() {
        let var = x_var();
        let g = x().pow(2).add(&Poly::one());
        let (terms, rem, mult) = hermite_reduce(&x(), &g, 2, &var).unwrap();
        assert_eq!(terms.len(), 1);
        assert_eq!(univariate_degree(&terms[0].numer, &var), 0);
        assert_eq!(terms[0].power, 1);
        assert!(rem.is_zero());
        assert_eq!(mult, 1);
    }

    #[test]
    fn hermite_x_over_x_squared_plus_one_squared() {
        let var = x_var();
        let g = x().pow(2).add(&Poly::one());
        let (terms, rem, mult) = hermite_reduce(&x(), &g, 2, &var).unwrap();
        assert_eq!(mult, 1);
        assert_eq!(terms.len(), 1);
        assert!(rem.is_zero());
    }

    #[test]
    fn hermite_one_over_x_fourth_plus_one_squared_rem() {
        let var = x_var();
        let g = x().pow(4).add(&Poly::one());
        let (_, rem, _) = hermite_reduce(&Poly::one(), &g, 2, &var).unwrap();
        assert_eq!(coeff_at(&rem, &var, 0), Ratio::new(3.into(), 4.into()));
    }

    #[test]
    fn hermite_one_over_x_fourth_plus_one_squared_terms() {
        let var = x_var();
        let g = x().pow(4).add(&Poly::one());
        let (terms, rem, mult) = hermite_reduce(&Poly::one(), &g, 2, &var).unwrap();
        assert_eq!(mult, 1);
        assert_eq!(terms.len(), 1);
        assert_eq!(univariate_degree(&terms[0].numer, &var), 1);
        assert!(univariate_degree(&rem, &var) <= 0);
    }

    #[test]
    fn hermite_one_over_x_fourth_plus_one_squared_v_numer() {
        let var = x_var();
        let g = x().pow(4).add(&Poly::one());
        let (terms, _, _) = hermite_reduce(&Poly::one(), &g, 2, &var).unwrap();
        assert_eq!(coeff_at(&terms[0].numer, &var, 1), Ratio::new((-1).into(), 4.into()));
    }

    #[test]
    fn hermite_one_over_x_squared_plus_one_squared_rem_half() {
        let var = x_var();
        let g = x().pow(2).add(&Poly::one());
        let (_, rem, _) = hermite_reduce(&Poly::one(), &g, 2, &var).unwrap();
        assert_eq!(coeff_at(&rem, &var, 0), Ratio::new(1.into(), 2.into()));
    }

    #[test]
    fn hermite_mult_one_unchanged() {
        let var = x_var();
        let g = x().pow(2).add(&Poly::one());
        let (terms, rem, mult) = hermite_reduce(&Poly::one(), &g, 1, &var).unwrap();
        assert!(terms.is_empty());
        assert_eq!(mult, 1);
        assert_eq!(rem, Poly::one());
    }

    #[test]
    fn hermite_two_step_cubic_power() {
        let var = x_var();
        let g = x().pow(2).sub(&Poly::one());
        let (terms, rem, mult) = hermite_reduce(&Poly::one(), &g, 3, &var).unwrap();
        assert_eq!(mult, 1);
        assert_eq!(terms.len(), 2);
        let _ = rem;
    }
}
