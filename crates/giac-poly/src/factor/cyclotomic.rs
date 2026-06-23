//! Cyclotomic polynomials Φ_n and specialized x^n ± 1 factorization.
//!
//! **Stable:** `cyclotomic_poly`.
//! **Partial:** `factor_xn_minus_one`, `try_factor_xn_minus_one`, `try_factor_xn_plus_one`, …

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{EvalError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::nested::{FlatUni, MainVar};
use crate::resultant::{coeff_at, univariate_degree};

use super::util::is_univariate_in;

/// **Pipeline private** — Positive divisors of `n`, sorted ascending.
pub fn divisors_u64(n: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut i = 1u64;
    while i * i <= n {
        if n % i == 0 {
            out.push(i);
            if i != n / i {
                out.push(n / i);
            }
        }
        i += 1;
    }
    out.sort();
    out
}

/// **Stable** — n-th cyclotomic polynomial Φ_n(x) over ℚ.
pub fn cyclotomic_poly(n: u64, var: &Var) -> PolyResult<Poly> {
    if n == 0 {
        return Err(EvalError::TypeError("cyclotomic n=0"));
    }
    if n == 1 {
        return Ok(Poly::var(var.clone()).sub(&Poly::one()));
    }
    let x = Poly::var(var.clone());
    let mut phi = FlatUni::new(x.pow(n).sub(&Poly::one()), MainVar::new(var.clone()));
    for d in divisors_u64(n) {
        if d == n {
            continue;
        }
        let sub = cyclotomic_poly(d, var)?;
        let sub_u = FlatUni::new(sub, MainVar::new(var.clone()));
        let (_, r) = phi.div_rem(&sub_u).expect("cyclotomic div_rem");
        if !r.is_zero() {
            return Err(EvalError::NotImplemented("cyclotomic division"));
        }
        phi = FlatUni::new(phi.exact_quo(&sub_u).expect("quotient"), MainVar::new(var.clone()));
    }
    Ok(phi.into_poly())
}

/// **Partial** — `x^n - 1 = ∏_{d|n} Φ_d(x)`.
pub fn factor_xn_minus_one(var: &Var, n: u64) -> PolyResult<Vec<Poly>> {
    if n == 0 {
        return Ok(vec![]);
    }
    divisors_u64(n)
        .into_iter()
        .map(|d| cyclotomic_poly(d, var))
        .collect()
}

/// **Partial** — Factors of `x^(2n)+x^n+1 = (x^(3n)-1)/(x^n-1)` via cyclotomic selection.
pub fn factor_x2n_plus_xn_plus_1(var: &Var, n: u64) -> PolyResult<Vec<Poly>> {
    let mut out = Vec::new();
    for d in divisors_u64(3 * n) {
        if n % d == 0 {
            continue;
        }
        out.push(cyclotomic_poly(d, var)?);
    }
    Ok(out)
}

/// **Partial** — detect and factor x^2n+x^n+1
pub fn try_factor_x2n_plus_xn_plus_1(p: &Poly, var: &Var) -> Option<Vec<Poly>> {
    if !is_univariate_in(p, var) {
        return None;
    }
    let d = univariate_degree(p, var);
    if d == 0 || d % 2 != 0 {
        return None;
    }
    let n = d / 2;
    if !is_x2n_plus_xn_plus_1_sparse(p, var, n) {
        return None;
    }
    let f = factor_x2n_plus_xn_plus_1(var, n).ok()?;
    let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
    if prod == *p || prod.neg() == *p {
        Some(f)
    } else {
        None
    }
}

// **Pipeline private** — `is_x2n_plus_xn_plus_1_sparse`
fn is_x2n_plus_xn_plus_1_sparse(p: &Poly, var: &Var, n: u64) -> bool {
    let d = 2 * n;
    if coeff_at(p, var, d) != num_rational::Ratio::one()
        || coeff_at(p, var, n) != num_rational::Ratio::one()
        || coeff_at(p, var, 0) != num_rational::Ratio::one()
    {
        return false;
    }
    for e in 1..d {
        if e != n && !coeff_at(p, var, e).is_zero() {
            return false;
        }
    }
    true
}

/// **Partial** — detect x^n-1
pub fn try_factor_xn_minus_one(p: &Poly, var: &Var) -> Option<Vec<Poly>> {
    if !is_univariate_in(p, var) {
        return None;
    }
    for (m, c) in &p.terms {
        if !m.is_const() && m.exp_of(var) > 0 && *c != num_rational::Ratio::one() {
            return None;
        }
        if m.is_const() && *c != Ratio::from_integer(-BigInt::one()) {
            return None;
        }
    }
    let n = univariate_degree(p, var);
    if n == 0 {
        return None;
    }
    if !is_xn_minus_one_poly(p, var, n) {
        return None;
    }
    factor_xn_minus_one(var, n).ok()
}

/// **Pipeline private** — shape test x^n-1
pub fn is_xn_minus_one_poly(p: &Poly, var: &Var, n: u64) -> bool {
    if !is_univariate_in(p, var) || p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_m1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(var) == n && *c == num_rational::Ratio::one() {
            has_xn = true;
        }
        if m.is_const() && *c == Ratio::from_integer(-BigInt::one()) {
            has_m1 = true;
        }
    }
    has_xn && has_m1
}

/// **Partial** — detect x^n+1
pub fn try_factor_xn_plus_one(p: &Poly, var: &Var) -> Option<Vec<Poly>> {
    if !is_univariate_in(p, var) || p.terms.len() != 2 {
        return None;
    }
    let n = univariate_degree(p, var);
    let mut has_xn = false;
    let mut has_p1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(var) == n && *c == num_rational::Ratio::one() {
            has_xn = true;
        }
        if m.is_const() && *c == num_rational::Ratio::one() {
            has_p1 = true;
        }
    }
    if !(has_xn && has_p1) {
        return None;
    }
    factor_xn_plus_one(var, n).ok()
}

/// `x^n + 1 = ∏_{d|2n, d∤n} Φ_d(x)`.
// **Pipeline private** — `factor_xn_plus_one`
fn factor_xn_plus_one(var: &Var, n: u64) -> PolyResult<Vec<Poly>> {
    let mut out = Vec::new();
    for d in divisors_u64(2 * n) {
        if n % d == 0 {
            continue;
        }
        out.push(cyclotomic_poly(d, var)?);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::Poly;

    #[test]
    fn cyclotomic_phi3() {
        let x = Poly::var("x");
        let phi = cyclotomic_poly(3, &Var::from("x")).unwrap();
        assert_eq!(phi, x.pow(2).add(&x).add(&Poly::one()));
    }

    #[test]
    fn factor_x100_plus_x50_plus_1() {
        let x = Poly::var("x");
        let p = x.pow(100).add(&x.pow(50)).add(&Poly::one());
        let f = try_factor_x2n_plus_xn_plus_1(&p, &Var::from("x")).unwrap();
        assert_eq!(f.len(), 6);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    fn factor_x10_minus_1() {
        let x = Poly::var("x");
        let p = x.pow(10).sub(&Poly::one());
        let f = try_factor_xn_minus_one(&p, &Var::from("x")).unwrap();
        assert!(f.len() >= 4);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }
}
