#![deny(unsafe_code)]
// ponytail: ext tower / field_session WIP — trim allows when API stabilizes
#![allow(dead_code, clippy::unwrap_used, clippy::expect_used)]
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!

pub mod alg_ext;
pub mod alg_ext_c;
pub mod ext_tower;
pub mod field_arith;
pub mod field_session;
pub(crate) mod galois_automorphism;
pub mod poly;
pub mod poly_alg_coeff;
pub mod poly_conv;
pub mod poly_alg_factor;
pub mod poly_alg_partfrac;
pub mod poly_alg_ops;
pub mod poly_alg_resultant;
pub mod poly_alg_sturm;
pub mod poly_roots;

#[cfg(test)]
pub mod test_fixtures;
