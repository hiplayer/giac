use num_bigint::BigInt;
use num_traits::Signed;

use crate::EvalError;
use crate::limits::MAX_POLY_EXPONENT;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MAX_POLY_EXPONENT;
    use num_bigint::BigInt;

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
}
