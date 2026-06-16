#![deny(unsafe_code)]

mod tresultant;
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
pub use tresultant::{
    biquartic_conjugate_pairs, eval_param_poly, num_minus_t_derivative, rational_roots_in_t,
    tresultant_eliminate_x, AlgebraicRt, ConjugatePair,
};
pub use univariate::{
    eval_univariate_at, odd_multiplicity_part, sign_variations, square_free_factorization,
    sturm_sequence, sturm_sign_variations_at, sturmab_count, univariate_derivative,
};

#[cfg(test)]
mod tests_phase2;
