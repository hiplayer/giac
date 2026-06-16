#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod desolve;
mod plugin;
mod stubs;

pub use desolve::eval_desolve;
pub use plugin::{install_ode, xcas_default, DefaultOdePlugin};
