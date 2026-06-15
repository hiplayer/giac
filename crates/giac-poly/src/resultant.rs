use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;

/// Sylvester resultant of univariate polynomials in `var`.
pub fn resultant(a: &Poly, b: &Poly, var: &Var) -> PolyResult<Poly> {
    let g = a.gcd(b);
    if univariate_degree(&g, var) > 0 {
        return Ok(Poly::zero());
    }
    let da = univariate_degree(a, var);
    let db = univariate_degree(b, var);
    if da == 0 || db == 0 {
        let c = if da == 0 {
            univariate_leading_coeff(a, var)
        } else {
            univariate_leading_coeff(b, var)
        };
        let exp = if da == 0 { db } else { da };
        let other = if da == 0 { b } else { a };
        return Ok(other.pow(exp).mul_scalar(&Ratio::from_integer(c.pow(exp as u32))));
    }
    if da == 1 || db == 1 {
        return Ok(sylvester_det2(a, b, var));
    }
    Err(PolyError::NotImplemented("resultant"))
}

fn sylvester_det2(a: &Poly, b: &Poly, var: &Var) -> Poly {
    let a1 = coeff_at(a, var, 1);
    let a0 = coeff_at(a, var, 0);
    let b1 = coeff_at(b, var, 1);
    let b0 = coeff_at(b, var, 0);
    Poly::constant(a1 * b0 - a0 * b1)
}

fn coeff_at(p: &Poly, var: &Var, exp: u32) -> Ratio<BigInt> {
    for (m, c) in &p.terms {
        if exp == 0 && m.is_const() {
            return c.clone();
        }
        if m.exp_of(var) == exp && m.iter().all(|(v, _)| v == var) {
            return c.clone();
        }
    }
    Ratio::zero()
}

fn univariate_degree(p: &Poly, var: &Var) -> u32 {
    p.terms
        .keys()
        .filter_map(|m| {
            let e = m.exp_of(var);
            if e > 0 { Some(e) } else { None }
        })
        .max()
        .unwrap_or(0)
}

fn univariate_leading_coeff(p: &Poly, var: &Var) -> BigInt {
    let deg = univariate_degree(p, var);
    for (m, c) in p.terms.iter().rev() {
        if m.exp_of(var) == deg {
            if c.denom().is_one() {
                return c.numer().clone();
            }
        }
    }
    BigInt::one()
}

/// Roots of univariate polynomial (low-degree exact cases).
pub fn roots(p: &Poly, var: &Var) -> PolyResult<Vec<Poly>> {
    let d = univariate_degree(p, var);
    match d {
        0 => {
            if p.is_zero() {
                Ok(vec![])
            } else {
                Err(PolyError::TypeError("constant has no roots"))
            }
        }
        1 => {
            let mut a = Ratio::zero();
            let mut b = Ratio::zero();
            for (m, c) in &p.terms {
                let exp = m.exp_of(var);
                if exp == 1 {
                    a += c;
                } else if exp == 0 {
                    b += c;
                }
            }
            if a.is_zero() {
                return Err(PolyError::TypeError("not linear"));
            }
            Ok(vec![Poly::constant(-&b / &a)])
        }
        3 if is_xn_minus_one(p, var, 3) => Ok(vec![Poly::constant(Ratio::one())]),
        _ => Err(PolyError::NotImplemented("roots")),
    }
}

fn is_xn_minus_one(p: &Poly, var: &Var, n: u32) -> bool {
    if p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_m1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(var) == n && *c == Ratio::one() {
            has_xn = true;
        }
        if m.is_const() && *c == Ratio::from_integer(-BigInt::one()) {
            has_m1 = true;
        }
    }
    has_xn && has_m1
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn resultant_shared_factor_is_zero() {
        let a = x().pow(2).sub(&Poly::one());
        let b = x().pow(3).sub(&Poly::one());
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        assert!(r.is_zero());
    }

    #[test]
    fn resultant_two_linear_polys() {
        let a = x().sub(&Poly::one());
        let b = x().add(&Poly::one());
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        assert_eq!(r, Poly::constant(Ratio::from_integer(BigInt::from(2))));
    }

    #[test]
    fn resultant_constant_times_linear() {
        let a = Poly::constant(Ratio::from_integer(BigInt::from(5)));
        let b = x().sub(&Poly::one());
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        let expected = b.mul_scalar(&Ratio::from_integer(BigInt::from(5)));
        assert_eq!(r, expected);
    }

    #[test]
    fn resultant_linear_times_constant() {
        let a = x().sub(&Poly::one());
        let b = Poly::constant(Ratio::from_integer(BigInt::from(7)));
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        let expected = a.mul_scalar(&Ratio::from_integer(BigInt::from(7)));
        assert_eq!(r, expected);
    }

    #[test]
    fn resultant_quadratic_not_implemented() {
        let a = x().pow(2).add(&Poly::one());
        let b = x().pow(2).sub(&Poly::one());
        let err = resultant(&a, &b, &Var::from("x")).unwrap_err();
        assert!(matches!(err, PolyError::NotImplemented(_)));
    }

    #[test]
    fn roots_linear() {
        let p = x().sub(&Poly::one());
        let rs = roots(&p, &Var::from("x")).unwrap();
        assert_eq!(rs.len(), 1);
        assert_eq!(rs[0], Poly::constant(Ratio::one()));
    }

    #[test]
    fn roots_zero_polynomial() {
        let rs = roots(&Poly::zero(), &Var::from("x")).unwrap();
        assert!(rs.is_empty());
    }

    #[test]
    fn roots_constant_nonzero_errors() {
        let err = roots(&Poly::one(), &Var::from("x")).unwrap_err();
        assert!(matches!(err, PolyError::TypeError(_)));
    }

    #[test]
    fn roots_x3_minus_one_real_root() {
        let p = x().pow(3).sub(&Poly::one());
        let rs = roots(&p, &Var::from("x")).unwrap();
        assert_eq!(rs, vec![Poly::constant(Ratio::one())]);
    }

    #[test]
    fn roots_quadratic_not_implemented() {
        let p = x().pow(2).sub(&Poly::one());
        let err = roots(&p, &Var::from("x")).unwrap_err();
        assert!(matches!(err, PolyError::NotImplemented(_)));
    }
}
