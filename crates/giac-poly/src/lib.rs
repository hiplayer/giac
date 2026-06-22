#![deny(unsafe_code)]
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
mod nested;
mod subresultant;
mod tresultant;
mod univariate;
mod error;
mod exp;
mod monomial;
mod poly;
mod poly_coeff;
mod modint;
mod modular;
mod resultant;
mod factor;
mod chinrem;
mod ops;
mod partfrac;

/// Dense poly1 arithmetic (crate-internal; see `dense::poly1`).
pub mod dense;
pub use nested::{
    BivariateEmbed, CoeffRingPoly, FlatUni, FlatUniQ, MainVar, MultivariatePoly, PrimitivePart,
    TnEmbed, UnivariateIn, UnivariateOver, UnivariatePoly, UnivariatePolyQ,
};
pub use factor::GoodEval;
pub use giac_error::EvalError;
pub use error::PolyResult;
pub use monomial::{Monomial, Var};
pub use poly::{Poly, PolyQ, abcuv, egcd, quo, rem, simp2};
pub use poly_coeff::PolyCoeff;
pub use modint::{ModInt, irem, smod};
pub use modular::{PolyMod, modp};
pub use resultant::{coeff_at, quadratic_abc, resultant, roots, univariate_coeffs_ascending, univariate_degree};
pub use factor::{
    as_perfect_power, factor_into, factor_mod_irreducibles, factor_poly, factor_poly_mod,
    quadratic_sqrt_factor_exprs, ratio_perfect_sqrt, try_linear_power, vars_in,
};
pub use partfrac::{partfrac_rational_terms, partfrac_terms};
pub use chinrem::{chinrem, chinrem_lists};
pub use ops::{content, gauss};
pub use tresultant::{
    biquadratic_res_conjugate_pairs, biquartic_conjugate_pairs, eval_param_poly,
    num_minus_t_derivative, rational_roots_in_t, tresultant_eliminate_x, AlgebraicRt,
    ConjugatePair,
};
pub use univariate::{
    eval_univariate_at, odd_multiplicity_part, sign_variations, square_free_factorization,
    square_free_part, substitute_univariate,
    sturm_sequence, sturm_sign_variations_at, sturmab_count, univariate_derivative,
};

#[cfg(test)]
mod tests_phase2;
