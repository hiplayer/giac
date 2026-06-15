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
