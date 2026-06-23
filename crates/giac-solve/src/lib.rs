#![deny(unsafe_code)]
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod plugin;
mod solve;
mod solve_poly;
mod rootof;
mod sturm;
mod fsolve;
mod froot;
mod realroot;
mod stubs;

#[cfg(test)]
mod test_verify;

pub use giac_linalg::eval_linsolve;
pub use plugin::{install_solve, xcas_default, DefaultSolvePlugin};
pub use solve::eval_solve;
pub use froot::eval_froot;
pub use realroot::eval_realroot;
pub use stubs::{eval_fsolve, eval_sturm, eval_sturmab};
