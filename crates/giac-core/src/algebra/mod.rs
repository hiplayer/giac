#![deny(unsafe_code)]

mod expand;
pub(crate) mod poly;
mod equiv;
mod factor;
mod normal;
mod ratnormal;
mod trig;

pub use expand::expand;
pub use equiv::{assert_equiv, is_zero, sub};
pub use factor::factor;
pub use normal::normal;
pub use ratnormal::ratnormal;
pub use trig::{halftan, lin, texpand};
