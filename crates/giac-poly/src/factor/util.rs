//! Shared utilities for factor pipeline: vars, content, primitive part, nth roots, …
//!
//! **Stable:** `vars_in`, `ratio_perfect_sqrt`, `coeff_wrt`, `main_var`, …

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use std::collections::{BTreeMap, BTreeSet};

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::monomial::{Monomial, Var};
use crate::nested::MultivariatePoly;
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

/// **Stable** — All variables appearing in `p`, lexicographically sorted.
pub fn vars_in(p: &Poly) -> Vec<Var> {
    let mut set = BTreeSet::new();
    for m in p.terms.keys() {
        for (v, e) in m.iter() {
            if e > 0 {
                set.insert(v.clone());
            }
        }
    }
    set.into_iter().collect()
}

/// **Stable** — `is_univariate_in`
pub fn is_univariate_in(p: &Poly, var: &Var) -> bool {
    p.terms
        .keys()
        .all(|m| m.iter().all(|(v, _)| v == var))
}

/// **Stable** — Variable of minimum degree (giac `factor_multivar` main var heuristic).
pub fn main_var(p: &Poly, vars: &[Var]) -> Var {
    vars.iter()
        .min_by_key(|v| univariate_degree(p, v))
        .cloned()
        .unwrap_or_else(|| vars[0].clone())
}

/// **Stable** — Minimum exponent of each variable across all terms (missing var counts as 0).
pub fn min_var_exponents(p: &Poly) -> BTreeMap<Var, u64> {
    let all_vars = vars_in(p);
    let mut out = BTreeMap::new();
    for v in all_vars {
        let min_e = p
            .terms
            .keys()
            .map(|m| m.exp_of(&v))
            .min()
            .unwrap_or(0);
        if min_e > 0 {
            out.insert(v, min_e);
        }
    }
    out
}

/// **Stable** — Split `p = (∏ v^{e_v}) * rest` where `e_v` is the minimum exponent of `v` in `p`.
pub fn extract_var_power_factors(p: &Poly) -> (Poly, Vec<Poly>) {
    let powers = min_var_exponents(p);
    if powers.is_empty() {
        return (p.clone(), Vec::new());
    }
    let mut divisor = Poly::one();
    for (v, e) in &powers {
        divisor = divisor.mul(&Poly::var(v.clone()).pow(*e));
    }
    let (rest, rem) = MultivariatePoly::new(p.clone()).div_rem(&divisor);
    if !rem.is_zero() {
        return (p.clone(), Vec::new());
    }
    let mut factors = Vec::new();
    for (v, e) in powers {
        for _ in 0..e {
            factors.push(Poly::var(v.clone()));
        }
    }
    (rest, factors)
}

/// **Stable** — Integer gcd of rational coefficients.
pub fn coeff_gcd(p: &Poly) -> Ratio<BigInt> {
    p.content()
}

/// **Stable** — divide out content
pub fn primitive_part(p: &Poly) -> Poly {
    p.primitive_part()
}

/// **Stable** — `linear_poly`
pub fn linear_poly(var: &Var, root: &Ratio<BigInt>) -> Poly {
    Poly::var(var.clone()).sub(&Poly::constant(root.clone()))
}

/// **Stable** — `monic_quadratic_poly`
pub fn monic_quadratic_poly(var: &Var, u: Ratio<BigInt>, v: Ratio<BigInt>) -> Poly {
    Poly::var(var.clone())
        .pow(2)
        .add(&Poly::var(var.clone()).mul_scalar(&u))
        .add(&Poly::constant(v))
}

/// **Stable** — `integer_divisors`
pub fn integer_divisors(n: &BigInt) -> Vec<BigInt> {
    if n.is_zero() {
        return vec![BigInt::zero()];
    }
    let a = n.abs();
    let mut divs = Vec::new();
    let mut i = BigInt::one();
    while &i * &i <= a {
        if (&a % &i).is_zero() {
            divs.push(i.clone());
            divs.push(&a / &i);
        }
        i += BigInt::one();
    }
    divs.sort();
    divs.dedup();
    divs
}

/// **Stable** — `integer_nth_root`
pub fn integer_nth_root(n: &BigInt, exp: u64) -> Option<BigInt> {
    if n.is_negative() && exp % 2 == 0 {
        return None;
    }
    let exp_u32 = u32::try_from(exp).ok()?;
    let mut lo = BigInt::zero();
    let mut hi = n.abs() + BigInt::one();
    while lo < hi {
        let mid = (&lo + &hi) / BigInt::from(2);
        let pow = mid.pow(exp_u32);
        match pow.cmp(n) {
            std::cmp::Ordering::Equal => return Some(if n.is_negative() { -mid } else { mid }),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    None
}

/// **Stable** — `rational_nth_root`
pub fn rational_nth_root(r: &Ratio<BigInt>, exp: u64) -> Option<Ratio<BigInt>> {
    let num = integer_nth_root(r.numer(), exp)?;
    let den = integer_nth_root(r.denom(), exp)?;
    Some(Ratio::new(num, den))
}

/// **Stable** — detect perfect square Ratio
pub fn ratio_perfect_sqrt(r: &Ratio<BigInt>) -> Option<Ratio<BigInt>> {
    if r.is_zero() {
        return Some(Ratio::zero());
    }
    let sn = integer_nth_root(r.numer(), 2)?;
    let sd = integer_nth_root(r.denom(), 2)?;
    Some(Ratio::new(sn, sd))
}

/// **Stable** — `rational_factor_pairs`
pub fn rational_factor_pairs(a0: &Ratio<BigInt>) -> Vec<(Ratio<BigInt>, Ratio<BigInt>)> {
    if a0.is_zero() {
        return vec![(Ratio::zero(), Ratio::one())];
    }
    let mut pairs = Vec::new();
    for p in integer_divisors(a0.numer()) {
        for q in integer_divisors(a0.denom()) {
            if q.is_zero() {
                continue;
            }
            let qq = Ratio::new(p.clone(), q.clone());
            let ss = a0 / qq.clone();
            pairs.push((qq.clone(), ss.clone()));
            if qq != ss {
                pairs.push((ss, qq));
            }
        }
    }
    pairs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    pairs.dedup();
    pairs
}

/// **Stable** — Coefficient of `var^exp` (quotient by `var^exp` on matching terms).
pub use crate::subresultant::coeff_wrt;
pub fn is_monic_univariate(p: &Poly, var: &Var) -> bool {
    let d = univariate_degree(p, var);
    coeff_at(p, var, d) == Ratio::one()
}
