#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod expand;
mod equiv;
mod ifactor;
mod factor;
mod ratnormal;
mod trig;
mod plugin;

pub use expand::{expand, normal};
pub use equiv::{assert_equiv, is_zero, sub};
pub use factor::factor;
pub use ifactor::ifactor;
pub use ratnormal::ratnormal;
pub use trig::{halftan, lin, texpand};
pub use plugin::{install_simplify, xcas_default, DefaultAlgebraPlugin};
