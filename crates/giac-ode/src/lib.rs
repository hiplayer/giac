#![deny(unsafe_code)]
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-ode-api-stability.md`.
//!
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod desolve;
mod plugin;
mod stubs;

pub mod test_verify;

pub use desolve::eval_desolve;
pub use plugin::{install_ode, xcas_default, DefaultOdePlugin};
