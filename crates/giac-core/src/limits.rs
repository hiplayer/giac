//! MVP bounds for the Rust port (see `.doc/known-divergences.md` §MVP 实施限制).

/// Maximum exponent when lowering `Expr::Pow` → `Poly` ([`expr_to_poly`](crate::algebra::poly::expr_to_poly)).
///
/// Aligned with giac default `GBASISF4_MAX_TOTALDEG` (`global.cc`, value `1024` in stock builds).
///
/// Validation path: `BigInt` → [`bigint_to_nonneg_u64`](crate::num_util::bigint_to_nonneg_u64)
/// → compare to this cap → `u64` for [`giac_poly::Poly::pow`]. Monomial exponents in
/// `giac-poly` are also `u64` end-to-end; only `BigInt::pow` still uses `u32` internally
/// via [`giac_poly::exp::bigint_pow`] when an exponent exceeds `u32::MAX`.
pub const MAX_POLY_EXPONENT: u64 = 1024;
