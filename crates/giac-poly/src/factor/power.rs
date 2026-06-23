//! Perfect-power detection: `(base)^exp` and `(linear)^n` shapes before general factor.
//!
//! **Stable:** `as_perfect_power`, `try_linear_power`.
//! **Pipeline private:** `try_nth_root`, `try_binomial_square`.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::monomial::{Monomial, Var};
use crate::nested::{FlatUni, MainVar};
use crate::poly::Poly;
use crate::resultant::coeff_at;

use super::util::{integer_nth_root, rational_nth_root};

/// **Stable** — If `p` is a perfect power, return `(base, exponent)`.
pub fn as_perfect_power(p: &Poly) -> Option<(Poly, u64)> {
    if p.degree() <= 1 {
        return None;
    }
    for exp in 2..=p.degree() {
        if p.degree() % exp != 0 {
            continue;
        }
        if let Some(base) = try_nth_root(p, exp) {
            return Some((base, exp));
        }
        if exp == 2 {
            if let Some(base) = try_binomial_square(p) {
                return Some((base, 2));
            }
        }
    }
    None
}

// **Pipeline private** — optional fallback `try_nth_root`
fn try_nth_root(p: &Poly, exp: u64) -> Option<Poly> {
    if p.terms.len() == 1 {
        return None;
    }
    if p.degree() != exp {
        return None;
    }
    let constant = p.terms.get(&Monomial::one())?;
    if constant.denom() != &BigInt::one() {
        return None;
    }
    let c = constant.numer();
    let a = integer_nth_root(c, exp)?;
    let mut base = Poly::constant(Ratio::from_integer(a));
    base = base.add(&Poly::var("x"));
    if base.pow(exp) == *p {
        return Some(base);
    }
    None
}

// **Pipeline private** — optional fallback `try_binomial_square`
fn try_binomial_square(p: &Poly) -> Option<Poly> {
    let d = p.degree();
    if d % 2 != 0 {
        return None;
    }
    let k = d / 2;
    let var = Var::from("x");
    for deg in 1..d {
        if deg != k && !coeff_at(p, &var, deg).is_zero() {
            return None;
        }
    }
    if coeff_at(p, &var, d) != Ratio::one() {
        return None;
    }
    let c2 = coeff_at(p, &var, 0);
    let c = rational_nth_root(&c2, 2)?;
    let mid = coeff_at(p, &var, k);
    if mid != Ratio::from_integer(2.into()) * c.clone() {
        return None;
    }
    let base = Poly::var("x").pow(k).add(&Poly::constant(c));
    if base.pow(2) == *p {
        Some(base)
    } else {
        None
    }
}

/// **Partial** — detect (linear)^n
pub fn try_linear_power(p: &Poly, var: &Var) -> Option<(Poly, u64)> {
    if let Some((base, exp)) = as_perfect_power(p) {
        if crate::resultant::univariate_degree(&base, var) == 1 {
            return Some((base, exp));
        }
    }
    let deg = crate::resultant::univariate_degree(p, var);
    if deg < 2 {
        return None;
    }
    let root = super::univariate::find_rational_root(p, var)?;
    let lin = super::util::linear_poly(var, &root);
    let flat = FlatUni::try_new(p.clone(), MainVar::new(var.clone())).ok()?;
    for exp in (2..=deg).rev() {
        let lin_pow = lin.pow(exp);
        let lin_u = FlatUni::try_new(lin_pow, MainVar::new(var.clone())).ok()?;
        if flat.div_rem(&lin_u).map(|(_, r)| r.is_zero()).unwrap_or(false) {
            return Some((lin, exp));
        }
    }
    None
}
