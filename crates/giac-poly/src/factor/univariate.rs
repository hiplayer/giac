//! Univariate factorization over ℚ: rational roots, quadratics, Zassenhaus, low-degree patterns.
//!
//! **Stable (bounded):** `factor_univariate_flat`, `factor_univariate_pairs`, `factor_power_pairs`.
//! **Partial:** `try_factor_biquadratic`, `try_factor_two_cubics`.
//! **Pipeline private:** `find_rational_root`, `factor_square_free`, `factor_quadratic`, …

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};
use crate::univariate::{eval_univariate_at, square_free_factorization};

use super::power::as_perfect_power;
use super::util::{
    integer_divisors, is_univariate_in, linear_poly, monic_quadratic_poly, rational_factor_pairs,
    ratio_perfect_sqrt,
};

/// **Partial** — Flat irreducible (or fully split) factor list.
pub fn factor_univariate_flat(p: &Poly, var: &Var) -> PolyResult<Vec<Poly>> {
    let pairs = factor_univariate_pairs(p, var)?;
    let mut out = Vec::new();
    for (f, m) in pairs {
        for _ in 0..m {
            out.push(f.clone());
        }
    }
    Ok(out)
}

/// **Partial** — pairs with multiplicity
pub fn factor_univariate_pairs(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    if p.is_zero() {
        return Err(PolyError::TypeError("zero polynomial"));
    }
    if p.is_one() {
        return Ok(vec![]);
    }
    if !is_univariate_in(p, var) {
        return Err(PolyError::NotImplemented("factor"));
    }
    let mut out = Vec::new();
    let content = p.content();
    if !content.is_one() && !content.is_zero() {
        out.push((Poly::constant(content), 1));
    }
    let mut pp = p.primitive_part();
    if pp.is_one() {
        return Ok(out);
    }
    if univariate_degree(&pp, var) == 1 {
        out.push((pp, 1));
        return Ok(out);
    }
    let sqff = square_free_factorization(&pp, var)?;
    for (g, k) in sqff {
        let mut factors = factor_square_free(&g, var)?;
        for f in factors {
            if let Some((last, m)) = out.last_mut() {
                if last == &f {
                    *m += k;
                    continue;
                }
            }
            out.push((f, k));
        }
    }
    Ok(out)
}

// **Pipeline private** — `factor_square_free`
fn factor_square_free(g: &Poly, var: &Var) -> PolyResult<Vec<Poly>> {
    if g.is_one() {
        return Ok(vec![]);
    }
    let d = univariate_degree(g, var);
    if d == 0 {
        return Ok(vec![g.clone()]);
    }
    if d == 1 {
        return Ok(vec![g.clone()]);
    }
    if let Some(factors) = factor_by_rational_roots(g, var) {
        return Ok(factors);
    }
    if d == 2 {
        return factor_quadratic(g, var);
    }
    if d == 4 {
        if let Some(biq) = try_factor_biquadratic(g, var) {
            let mut out = Vec::new();
            for (f, m) in biq {
                for _ in 0..m {
                    out.push(f.clone());
                }
            }
            return Ok(out);
        }
    }
    if d >= 3 {
        if let Some(facs) = super::zassenhaus::try_zassenhaus_factor(g, var) {
            return Ok(facs);
        }
    }
    if d == 6 {
        if let Some(facs) = try_factor_two_cubics(g, var) {
            return Ok(facs);
        }
    }
    if let Some((base, exp)) = as_perfect_power(g) {
        if crate::resultant::univariate_degree(&base, var) <= 2 {
            let inner = factor_square_free(&base, var)?;
            let mut out = Vec::new();
            for f in inner {
                for _ in 0..exp as usize {
                    out.push(f.clone());
                }
            }
            return Ok(out);
        }
    }
    Ok(vec![g.clone()])
}

// **Pipeline private** — `factor_by_rational_roots`
fn factor_by_rational_roots(p: &Poly, var: &Var) -> Option<Vec<Poly>> {
    let pairs = factor_power_pairs_core(p, var).ok()?;
    let mut out = Vec::new();
    for (f, m) in pairs {
        for _ in 0..m {
            out.push(f.clone());
        }
    }
    Some(out)
}

/// **Stable (bounded)** — factors with multiplicities
pub fn factor_power_pairs(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    if p.is_zero() {
        return Err(PolyError::TypeError("zero polynomial"));
    }
    if p.is_one() {
        return Ok(vec![]);
    }
    if !is_univariate_in(p, var) {
        return Err(PolyError::NotImplemented("factor"));
    }
    let mut out = Vec::new();
    let content = p.content();
    if !content.is_one() && !content.is_zero() {
        out.push((Poly::constant(content), 1));
    }
    let pp = p.primitive_part();
    if pp.is_one() {
        return Ok(out);
    }
    out.extend(factor_power_pairs_core(&pp, var)?);
    Ok(out)
}

// **Pipeline private** — `factor_power_pairs_core`
fn factor_power_pairs_core(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    let mut rest = p.clone();
    let mut factors = Vec::new();
    while univariate_degree(&rest, var) > 0 {
        let root = match find_rational_root(&rest, var) {
            Some(r) => r,
            None => {
                if univariate_degree(&rest, var) == 4 {
                    if let Some(biq) = try_factor_biquadratic(&rest, var) {
                        factors.extend(biq);
                        rest = Poly::one();
                        break;
                    }
                }
                if univariate_degree(&rest, var) <= 2 {
                    break;
                }
                return Err(PolyError::NotImplemented("factor"));
            }
        };
        let lin = linear_poly(var, &root);
        let mut mult = 0usize;
        loop {
            let (_, r) = rest.div_rem(&lin);
            if !r.is_zero() {
                break;
            }
            mult += 1;
            rest = rest.div_rem(&lin).0;
        }
        if mult == 0 {
            return Err(PolyError::NotImplemented("factor"));
        }
        factors.push((lin, mult));
    }
    if rest.is_one() || rest.is_zero() {
        return Ok(factors);
    }
    if let Some((base, exp)) = as_perfect_power(&rest) {
        if univariate_degree(&base, var) <= 2 {
            factors.push((base, exp as usize));
            return Ok(factors);
        }
    }
    if univariate_degree(&rest, var) <= 2 {
        factors.push((rest, 1));
        return Ok(factors);
    }
    Err(PolyError::NotImplemented("factor"))
}

// **Pipeline private** — rational root via rational root theorem
pub(crate) fn find_rational_root(p: &Poly, var: &Var) -> Option<Ratio<BigInt>> {
    let deg = univariate_degree(p, var);
    if deg == 0 {
        return None;
    }
    let a0 = coeff_at(p, var, 0);
    let an = coeff_at(p, var, deg);
    for p_cand in integer_divisors(a0.numer()) {
        for q_cand in integer_divisors(an.numer()) {
            if q_cand.is_zero() {
                continue;
            }
            for &(pn, qn) in &[(1, 1), (1, -1), (-1, 1), (-1, -1)] {
                let r = Ratio::new(&p_cand * pn, &q_cand * qn);
                if eval_univariate_at(p, var, &r).is_zero() {
                    let lin = linear_poly(var, &r);
                    if p.div_rem(&lin).1.is_zero() {
                        return Some(r);
                    }
                }
            }
        }
    }
    None
}

// **Pipeline private** — `factor_quadratic`
fn factor_quadratic(p: &Poly, var: &Var) -> PolyResult<Vec<Poly>> {
    let mut a = Ratio::zero();
    let mut b = Ratio::zero();
    let mut c = Ratio::zero();
    for (m, coeff) in &p.terms {
        match m.exp_of(var) {
            2 => a += coeff,
            1 => b += coeff,
            0 => c += coeff,
            _ => return Err(PolyError::TypeError("not quadratic")),
        }
    }
    if a.is_zero() {
        return Ok(vec![p.clone()]);
    }
    let disc = &b * &b - Ratio::from_integer(BigInt::from(4)) * &a * &c;
    if disc.is_zero() {
        let r = -&b / (Ratio::from_integer(BigInt::from(2)) * &a);
        return Ok(vec![linear_poly(var, &r)]);
    }
    if let Some(sqrt_d) = ratio_perfect_sqrt(&disc) {
        let two_a = Ratio::from_integer(BigInt::from(2)) * &a;
        let r1 = (-&b + &sqrt_d) / &two_a;
        let r2 = (-&b - &sqrt_d) / &two_a;
        return Ok(vec![linear_poly(var, &r1), linear_poly(var, &r2)]);
    }
    Ok(vec![p.clone()])
}

// **Pipeline private** — optional fallback `try_factor_biquadratic`
fn try_factor_biquadratic(p: &Poly, var: &Var) -> Option<Vec<(Poly, usize)>> {
    if univariate_degree(p, var) != 4 {
        return None;
    }
    let lc = coeff_at(p, var, 4);
    if lc.is_zero() {
        return None;
    }
    let scale = Ratio::one() / lc.clone();
    let a3 = coeff_at(p, var, 3) * scale.clone();
    let a2 = coeff_at(p, var, 2) * scale.clone();
    let a1 = coeff_at(p, var, 1) * scale.clone();
    let a0 = coeff_at(p, var, 0) * scale;
    for (q, s) in rational_factor_pairs(&a0) {
        let sum_pr = a2.clone() - q.clone() - s.clone();
        let disc = a3.clone() * a3.clone()
            - Ratio::from_integer(BigInt::from(4)) * sum_pr.clone();
        if disc < Ratio::zero() {
            continue;
        }
        let sqrt_d = ratio_perfect_sqrt(&disc)?;
        let two = Ratio::from_integer(BigInt::from(2));
        let p_coef = (a3.clone() + sqrt_d.clone()) / two.clone();
        let r_coef = (a3.clone() - sqrt_d) / two;
        if p_coef.clone() * s.clone() + q.clone() * r_coef.clone() != a1 {
            continue;
        }
        let f1 = monic_quadratic_poly(var, p_coef, q);
        let f2 = monic_quadratic_poly(var, r_coef, s);
        let prod = f1.clone().mul(&f2);
        if prod == *p {
            return Some(vec![(f1, 1), (f2, 1)]);
        }
        if prod.neg() == *p {
            return Some(vec![(f1.neg(), 1), (f2, 1)]);
        }
    }
    None
}

// **Pipeline private** — `monic_cubic_poly`
fn monic_cubic_poly(var: &Var, a2: i64, a1: i64, a0: i64) -> Poly {
    let x = Poly::var(var.clone());
    let mut out = x.pow(3);
    if a2 != 0 {
        out = out.add(&x.pow(2).mul_scalar(&Ratio::from_integer(BigInt::from(a2))));
    }
    if a1 != 0 {
        out = out.add(&x.mul_scalar(&Ratio::from_integer(BigInt::from(a1))));
    }
    if a0 != 0 {
        out = out.add(&Poly::constant(Ratio::from_integer(BigInt::from(a0))));
    }
    out
}

/// Split a degree-6 square-free polynomial into two monic cubics (Issue 3.1 MVP).
// **Pipeline private** — optional fallback `try_factor_two_cubics`
fn try_factor_two_cubics(p: &Poly, var: &Var) -> Option<Vec<Poly>> {
    if univariate_degree(p, var) != 6 {
        return None;
    }
    let lc = coeff_at(p, var, 6);
    if lc.is_zero() {
        return None;
    }
    let scale = Ratio::one() / lc.clone();
    let mut p_m = Poly::zero();
    for e in 0..=6 {
        let c = coeff_at(p, var, e) * scale.clone();
        if !c.is_zero() {
            p_m = p_m.add(&term_with_var(&Poly::constant(c), var, e));
        }
    }

    let bound = 8i64;
    for a2 in -bound..=bound {
        for a1 in -bound..=bound {
            for a0 in -bound..=bound {
                let f = monic_cubic_poly(var, a2, a1, a0);
                let (_, r) = p_m.div_rem(&f);
                if !r.is_zero() {
                    continue;
                }
                let g = p_m.div_rem(&f).0;
                if univariate_degree(&g, var) != 3 {
                    continue;
                }
                let f = f.mul_scalar(&lc.clone());
                let g = g.mul_scalar(&lc);
                if f.mul(&g) == *p {
                    return Some(vec![f, g]);
                }
            }
        }
    }
    None
}

// **Stable** — coeff * var^exp as Poly
fn term_with_var(coeff: &Poly, var: &Var, exp: u64) -> Poly {
    if exp == 0 {
        return coeff.clone();
    }
    coeff.mul(&Poly::var(var.clone()).pow(exp))
}
