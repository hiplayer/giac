#![deny(unsafe_code)]

mod expand;
pub(crate) mod poly;
mod factor;
mod normal;
mod ratnormal;

pub use expand::expand;
pub use factor::factor;
pub use normal::normal;
pub use ratnormal::ratnormal;
