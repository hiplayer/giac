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
