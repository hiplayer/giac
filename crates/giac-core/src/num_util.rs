use num_bigint::BigInt;
use num_integer::gcd;
use num_traits::{One, Signed, Zero};

use crate::EvalError;
use crate::limits::MAX_POLY_EXPONENT;

/// Exact integer `k`-th root when `n = r^k` for some integer `r`; otherwise `None`.
///
/// Negative `n` is allowed only for odd `k`.
pub fn integer_nth_root(n: &BigInt, exp: u64) -> Option<BigInt> {
    if n.is_negative() && exp % 2 == 0 {
        return None;
    }
    let exp_u32 = u32::try_from(exp).ok()?;
    let mut lo = BigInt::zero();
    let mut hi = n.abs() + BigInt::one();
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

/// Exact integer square root when `n` is a perfect square; otherwise `None`.
pub fn integer_sqrt(n: &BigInt) -> Option<BigInt> {
    integer_nth_root(n, 2)
}

/// Convert a non-negative `BigInt` to `u32`, for bounded exponents etc.
pub fn bigint_to_nonneg_u32(n: &BigInt) -> Result<u32, EvalError> {
    if n.is_negative() {
        return Err(EvalError::TypeError("expected non-negative integer"));
    }
    n.to_string()
        .parse()
        .map_err(|_| EvalError::TypeError("integer out of u32 range"))
}

/// Convert a non-negative `BigInt` to `u64`.
pub fn bigint_to_nonneg_u64(n: &BigInt) -> Result<u64, EvalError> {
    if n.is_negative() {
        return Err(EvalError::TypeError("expected non-negative integer"));
    }
    n.to_string()
        .parse()
        .map_err(|_| EvalError::TypeError("integer out of u64 range"))
}

/// Convert a `BigInt` to `i64` when it fits.
pub fn bigint_to_i64(n: &BigInt) -> Result<i64, EvalError> {
    n.to_string()
        .parse()
        .map_err(|_| EvalError::TypeError("integer out of i64 range"))
}

/// Polynomial power exponent for `expr_to_poly`: non-negative, `<= MAX_POLY_EXPONENT`.
pub fn bigint_to_poly_exponent(n: &BigInt) -> Result<u64, EvalError> {
    let e = bigint_to_nonneg_u64(n)?;
    if e > MAX_POLY_EXPONENT {
        return Err(EvalError::TypeError("polynomial exponent exceeds limit"));
    }
    Ok(e)
}

/// Absolute value as `u32` (for negative exponents on integers).
pub fn bigint_to_u32_abs(n: &BigInt) -> Result<u32, EvalError> {
    bigint_to_nonneg_u32(&n.abs())
}

/// Best rational approximation of `v` via continued fractions.
///
/// `tol` is the absolute error tolerated when reconstructing `v` as `n/d` as `f64`.
#[allow(dead_code)]
pub fn f64_to_rational(v: f64, tol: f64) -> (i64, i64) {
    if !v.is_finite() {
        return (0, 1);
    }
    if v.abs() < tol {
        return (0, 1);
    }
    let sign = if v < 0.0 { -1 } else { 1 };
    let v = v.abs();
    let a0 = v.floor() as i64;
    if (v - a0 as f64).abs() < tol {
        return (sign * a0, 1);
    }

    let mut p0 = 1_i64;
    let mut p1 = a0;
    let mut q0 = 0_i64;
    let mut q1 = 1_i64;
    let mut x = v - a0 as f64;
    const MAX_DEN: i64 = 10_000_000_000_000;

    loop {
        if x.abs() < tol {
            break;
        }
        x = 1.0 / x;
        let a = x.floor() as i64;
        let p = a.saturating_mul(p1).saturating_add(p0);
        let q = a.saturating_mul(q1).saturating_add(q0);
        if q > MAX_DEN || q == 0 {
            break;
        }
        let approx = p as f64 / q as f64;
        if (approx - v).abs() < tol {
            return reduce_rational(sign * p, q);
        }
        p0 = p1;
        p1 = p;
        q0 = q1;
        q1 = q;
        x -= a as f64;
    }
    reduce_rational(sign * p1, q1)
}

fn reduce_rational(n: i64, d: i64) -> (i64, i64) {
    if d == 0 {
        return (0, 1);
    }
    let (n, d) = if d < 0 { (-n, -d) } else { (n, d) };
    let g = gcd(n.abs(), d);
    (n / g, d / g)
}

/// Reduce a rational pair (for use outside this module).
pub fn reduce_rational_pair(n: i64, d: i64) -> (i64, i64) {
    reduce_rational(n, d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MAX_POLY_EXPONENT;
    use num_bigint::BigInt;

    #[test]
    fn integer_sqrt_perfect_and_non_perfect() {
        assert_eq!(integer_sqrt(&BigInt::from(9)).unwrap(), BigInt::from(3));
        assert_eq!(integer_sqrt(&BigInt::from(0)).unwrap(), BigInt::from(0));
        assert!(integer_sqrt(&BigInt::from(2)).is_none());
        assert!(integer_sqrt(&BigInt::from(-4)).is_none());
    }

    #[test]
    fn integer_nth_root_cube() {
        assert_eq!(integer_nth_root(&BigInt::from(27), 3).unwrap(), BigInt::from(3));
        assert!(integer_nth_root(&BigInt::from(-8), 2).is_none());
    }

    #[test]
    fn bigint_to_nonneg_u64_ok() {
        assert_eq!(bigint_to_nonneg_u64(&BigInt::from(1024)).unwrap(), 1024);
    }

    #[test]
    fn bigint_to_poly_exponent_at_limit() {
        assert_eq!(
            bigint_to_poly_exponent(&BigInt::from(MAX_POLY_EXPONENT)).unwrap(),
            MAX_POLY_EXPONENT
        );
    }

    #[test]
    fn bigint_to_poly_exponent_over_limit() {
        let over = BigInt::from(MAX_POLY_EXPONENT) + BigInt::from(1);
        assert!(bigint_to_poly_exponent(&over).is_err());
    }

    #[test]
    fn bigint_to_nonneg_u64_overflow() {
        let huge = BigInt::from(u64::MAX) + BigInt::from(1);
        assert!(matches!(
            bigint_to_nonneg_u64(&huge),
            Err(EvalError::TypeError("integer out of u64 range"))
        ));
    }

    #[test]
    fn bigint_to_nonneg_u32_ok() {
        assert_eq!(bigint_to_nonneg_u32(&BigInt::from(42)).unwrap(), 42);
    }

    #[test]
    fn bigint_to_nonneg_u32_negative() {
        assert_eq!(
            bigint_to_nonneg_u32(&BigInt::from(-1)),
            Err(EvalError::TypeError("expected non-negative integer"))
        );
    }

    #[test]
    fn bigint_to_nonneg_u32_overflow() {
        let huge = BigInt::from(u64::MAX) + BigInt::from(1);
        assert!(bigint_to_nonneg_u32(&huge).is_err());
    }

    #[test]
    fn bigint_to_u32_abs_negative_input() {
        assert_eq!(bigint_to_u32_abs(&BigInt::from(-7)).unwrap(), 7);
    }

    #[test]
    fn f64_to_rational_integer() {
        assert_eq!(f64_to_rational(5.0, 1e-12), (5, 1));
    }

    #[test]
    fn f64_to_rational_half() {
        let (n, d) = f64_to_rational(0.5, 1e-12);
        assert_eq!(n, 1);
        assert_eq!(d, 2);
    }

    #[test]
    fn f64_to_rational_reconstructs_svd_entry() {
        let v = 202_277_f64 / 500_000_f64;
        let (n, d) = f64_to_rational(v, 1e-12);
        assert!((n as f64 / d as f64 - v).abs() < 1e-12);
    }
}
