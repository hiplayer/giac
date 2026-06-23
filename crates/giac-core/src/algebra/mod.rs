#![deny(unsafe_code)]
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
pub mod poly_alg_ops;
pub mod poly_roots;

#[cfg(test)]
pub mod test_fixtures;
