#![deny(unsafe_code)]

mod error;
mod exp;
mod monomial;
mod poly;
mod modint;
mod modular;
mod resultant;
mod factor;
mod chinrem;
mod ops;

pub use error::{PolyError, PolyResult};
pub use monomial::{Monomial, Var};
pub use poly::{Poly, abcuv, egcd, quo, rem, simp2};
pub use modint::{ModInt, irem, smod};
pub use modular::{PolyMod, modp};
pub use resultant::{resultant, roots};
pub use factor::{as_perfect_power, factor_into, factor_poly, factor_poly_mod};
pub use chinrem::{chinrem, chinrem_lists};
pub use ops::{content, gauss};

#[cfg(test)]
mod tests_phase2;
