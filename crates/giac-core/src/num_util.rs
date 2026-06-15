use num_bigint::BigInt;
use num_traits::Signed;

use crate::EvalError;

/// Convert a non-negative `BigInt` to `u32`, for bounded exponents etc.
pub fn bigint_to_nonneg_u32(n: &BigInt) -> Result<u32, EvalError> {
    if n.is_negative() {
        return Err(EvalError::TypeError("expected non-negative integer"));
    }
    n.to_string()
        .parse()
        .map_err(|_| EvalError::TypeError("integer out of u32 range"))
}

/// Absolute value as `u32` (for negative exponents on integers).
pub fn bigint_to_u32_abs(n: &BigInt) -> Result<u32, EvalError> {
    bigint_to_nonneg_u32(&n.abs())
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

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
