//! MVP bounds for the Rust port (see `.doc/known-divergences.md` §MVP 实施限制).

/// Maximum exponent when lowering `Expr::Pow` → `Poly` ([`expr_to_poly`](crate::algebra::poly::expr_to_poly)).
///
/// Aligned with giac default `GBASISF4_MAX_TOTALDEG` (`global.cc`, value `1024` in stock builds).
/// That constant guards F4 Groebner total degree; we reuse the same numeric bound at the
/// `Expr` → `Poly` bridge so oversized powers fail with `TypeError` instead of unbounded work.
///
/// Validation path: `BigInt` → `u64` ([`bigint_to_nonneg_u64`]) → compare to this cap → `u32`
/// for [`giac_poly::Poly::pow`]. The `u64` step only parses arbitrary integer literals safely;
/// internal monomial exponents stay `u32` because every allowed exponent fits in `u32` when
/// this cap is `1024` (see `known-divergences.md`).
pub const MAX_POLY_EXPONENT: u64 = 1024;
