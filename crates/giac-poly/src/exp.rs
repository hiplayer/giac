//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use num_bigint::BigInt;

/// `num-bigint` exposes `BigInt::pow(u32)` only; bridge from unified `u64` exponents.
/// **Pipeline private** — BigInt pow with overflow check
pub fn bigint_pow(b: &BigInt, exp: u64) -> Option<BigInt> {
    Some(b.pow(u32::try_from(exp).ok()?))
}
