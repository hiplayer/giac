//! Dense univariate polynomial arithmetic (giac `poly1` convention).
//!
//! **Tier:** Stable (crate-internal) — see [GIAC-dense-poly1-refactor](.doc/issues/GIAC-dense-poly1-refactor.md) D1.

mod convert;
mod poly1;
mod ratio_ring;

#[cfg(test)]
mod tests;

pub use convert::{
    ascending_to_dense_high_first, dense_high_first_to_ascending, dense_high_first_to_sparse,
    reverse_coeffs, sparse_ascending_to_dense_high_first,
};
pub use poly1::{
    Poly1Order, Poly1RingCtx, add, div_rem, ext_gcd, inv_mod, mul, neg, poly_degree, reduce_mod_monic,
    scale, sub, trim,
};
pub use ratio_ring::{RatioRingCtx, RatioRingOps};
