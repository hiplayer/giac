#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod plugin;
mod solve;
mod stubs;

pub use giac_linalg::eval_linsolve;
pub use plugin::{install_solve, xcas_default, DefaultSolvePlugin};
pub use solve::eval_solve;
pub use stubs::{eval_fsolve, eval_realroot, eval_sturm};
