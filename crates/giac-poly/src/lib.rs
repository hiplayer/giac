#![deny(unsafe_code)]

mod univariate;
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
mod partfrac;

pub use error::{PolyError, PolyResult};
pub use monomial::{Monomial, Var};
pub use poly::{Poly, abcuv, egcd, quo, rem, simp2};
pub use modint::{ModInt, irem, smod};
pub use modular::{PolyMod, modp};
pub use resultant::{coeff_at, resultant, roots, univariate_degree};
pub use factor::{as_perfect_power, factor_into, factor_poly, factor_poly_mod};
pub use partfrac::{partfrac_rational_terms, partfrac_terms};
pub use chinrem::{chinrem, chinrem_lists};
pub use ops::{content, gauss};
pub use univariate::{
    eval_univariate_at, odd_multiplicity_part, sign_variations, sturm_sequence,
    sturm_sign_variations_at, sturmab_count, univariate_derivative,
};

#[cfg(test)]
mod tests_phase2;
