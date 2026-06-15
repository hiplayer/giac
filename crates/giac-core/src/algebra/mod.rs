#![deny(unsafe_code)]

mod expand;
pub(crate) mod poly;
mod equiv;
mod factor;
mod normal;
mod ratnormal;

pub use expand::expand;
pub use equiv::{assert_equiv, is_zero, sub};
pub use factor::factor;
pub use normal::normal;
pub use ratnormal::ratnormal;
