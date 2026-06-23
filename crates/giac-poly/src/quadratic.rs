//! Quadratic polynomials over ℚ: coefficient extraction and rational roots.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{Signed, Zero};

use crate::error::{EvalError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;

/// Coefficients `(a, b, c)` of `a·var² + b·var + c`; `None` if not quadratic or `a = 0`.
pub fn quadratic_abc(
    p: &Poly,
    var: &Var,
) -> Option<(Ratio<BigInt>, Ratio<BigInt>, Ratio<BigInt>)> {
    let mut a = Ratio::zero();
    let mut b = Ratio::zero();
    let mut c = Ratio::zero();
    let mut has_quad = false;
    for (m, coeff) in &p.terms {
        if coeff.is_zero() {
            continue;
        }
        match m.exp_of(var) {
            2 => {
                has_quad = true;
                a += coeff;
            }
            1 => b += coeff,
            0 => c += coeff,
            _ => return None,
        }
    }
    if !has_quad || a.is_zero() {
        return None;
    }
    Some((a, b, c))
}

/// **Stable** — `(a, b, c)` for `a·var² + b·var + c`.
pub struct QuadraticCoeffs<C> {
    pub a: C,
    pub b: C,
    pub c: C,
}

/// **Stable** — extract `(a, b, c)` from a quadratic univariate `p`.
pub fn quadratic_coeffs(p: &Poly, var: &Var) -> Option<QuadraticCoeffs<Ratio<BigInt>>> {
    quadratic_abc(p, var).map(|(a, b, c)| QuadraticCoeffs { a, b, c })
}

/// **Stable (bounded)** — exact rational roots when the discriminant is a square in ℚ.
pub fn quadratic_rational_roots(
    q: &QuadraticCoeffs<Ratio<BigInt>>,
) -> PolyResult<Vec<Ratio<BigInt>>> {
    let disc = &q.b * &q.b - Ratio::from_integer(BigInt::from(4)) * &q.a * &q.c;
    if disc.is_zero() {
        let r = -&q.b / (Ratio::from_integer(BigInt::from(2)) * &q.a);
        return Ok(vec![r]);
    }
    let sqrt_d = ratio_perfect_sqrt(&disc).ok_or(EvalError::NotImplemented("roots"))?;
    let two_a = Ratio::from_integer(BigInt::from(2)) * &q.a;
    let r1 = (-&q.b + &sqrt_d) / &two_a;
    let r2 = (-&q.b - &sqrt_d) / &two_a;
    Ok(vec![r1, r2])
}

// ponytail: duplicate of `factor::util::ratio_perfect_sqrt` — breaks factor↔resultant cycle
fn ratio_perfect_sqrt(r: &Ratio<BigInt>) -> Option<Ratio<BigInt>> {
    if r.is_zero() {
        return Some(Ratio::zero());
    }
    let sn = integer_nth_root(r.numer(), 2)?;
    let sd = integer_nth_root(r.denom(), 2)?;
    Some(Ratio::new(sn, sd))
}

fn integer_nth_root(n: &BigInt, exp: u64) -> Option<BigInt> {
    if n.is_negative() && exp % 2 == 0 {
        return None;
    }
    let exp_u32 = u32::try_from(exp).ok()?;
    let mut lo = BigInt::from(0);
    let mut hi = n.abs() + BigInt::from(1);
    while lo < hi {
        let mid = (&lo + &hi) / BigInt::from(2);
        let pow = mid.pow(exp_u32);
        match pow.cmp(n) {
            std::cmp::Ordering::Equal => return Some(if n.is_negative() { -mid } else { mid }),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    None
}
