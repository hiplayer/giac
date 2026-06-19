//! Fast-path pattern factorization (x^n±1, x^n±y^n, cyclotomic hooks) before Hensel.
//!
//! **Partial:** `try_factor_patterns`, `factor_xn_minus_one_display`.
//! **Pipeline private:** `factor_xn_minus_yn`, `is_binomial_diff_power`, …

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::One;

use crate::monomial::Var;
use crate::poly::Poly;

use super::cyclotomic::{
    try_factor_x2n_plus_xn_plus_1, try_factor_xn_minus_one, try_factor_xn_plus_one,
};
use super::util::vars_in;

pub use super::cyclotomic::{factor_xn_minus_one, is_xn_minus_one_poly};

/// **Partial** — cyclotomic/binomial pattern table
pub fn try_factor_patterns(p: &Poly) -> Option<Vec<Poly>> {
    let vars = vars_in(p);
    if vars.len() == 1 {
        let v = &vars[0];
        if let Some(f) = try_factor_xn_minus_one(p, v) {
            return Some(f);
        }
        if let Some(f) = try_factor_xn_plus_one(p, v) {
            return Some(f);
        }
        if let Some(f) = try_factor_x2n_plus_xn_plus_1(p, v) {
            return Some(f);
        }
        return None;
    }
    if vars.len() == 2 {
        if let Some(f) = factor_xn_minus_yn(p, &vars[0], &vars[1]) {
            return Some(f);
        }
    }
    None
}

/// `x^n - y^n` (homogeneous binomial difference).
// **Pipeline private** — `factor_xn_minus_yn`
fn factor_xn_minus_yn(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    let xv = Poly::var(x.clone());
    let yv = Poly::var(y.clone());
    for n in [2u64, 3, 4, 6] {
        if is_binomial_diff_power(p, x, y, n) {
            return Some(factor_xn_minus_yn_explicit(&xv, &yv, n));
        }
    }
    None
}

// **Pipeline private** — `is_binomial_diff_power`
fn is_binomial_diff_power(p: &Poly, x: &Var, y: &Var, n: u64) -> bool {
    if p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_yn = false;
    for (m, c) in &p.terms {
        if m.exp_of(x) == n && m.exp_of(y) == 0 && *c == Ratio::one() {
            has_xn = true;
        }
        if m.exp_of(y) == n && m.exp_of(x) == 0 && *c == Ratio::from_integer(-BigInt::one()) {
            has_yn = true;
        }
    }
    has_xn && has_yn
}

// **Pipeline private** — `factor_xn_minus_yn_explicit`
fn factor_xn_minus_yn_explicit(x: &Poly, y: &Poly, n: u64) -> Vec<Poly> {
    match n {
        2 => vec![x.sub(y), x.add(y)],
        3 => vec![x.sub(y), x.pow(2).add(&x.mul(y)).add(&y.pow(2))],
        4 => vec![
            x.sub(y),
            x.add(y),
            x.pow(2).sub(&x.mul(y)).add(&y.pow(2)),
            x.pow(2).add(&x.mul(y)).add(&y.pow(2)),
        ],
        6 => vec![
            x.sub(y),
            x.add(y),
            x.pow(2).sub(&x.mul(y)).add(&y.pow(2)),
            x.pow(2).add(&x.mul(y)).add(&y.pow(2)),
        ],
        _ => vec![x.pow(n).sub(&y.pow(n))],
    }
}

/// **Partial** — Legacy wrapper used by `factor_poly`.
pub fn factor_xn_minus_one_display(p: &Poly) -> Option<Poly> {
    let vars = vars_in(p);
    if vars.len() != 1 {
        return None;
    }
    let v = &vars[0];
    if is_xn_minus_one_poly(p, v, 4) {
        let x = Poly::var(v.clone());
        return Some(
            x.sub(&Poly::one())
                .mul(&x.add(&Poly::one()))
                .mul(&x.pow(2).add(&Poly::one())),
        );
    }
    None
}
