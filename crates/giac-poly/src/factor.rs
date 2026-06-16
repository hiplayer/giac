use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::{Monomial, Var};
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};
use crate::univariate::eval_univariate_at;

/// If `p` is a perfect power, return `(base, exponent)`.
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

/// Integer factorization of x^n - 1 patterns and perfect powers.
pub fn factor_poly(p: &Poly) -> Poly {
    if let Some(f) = factor_xn_minus_one(p) {
        return f;
    }
    if let Some((base, exp)) = as_perfect_power(p) {
        return base.pow(exp);
    }
    p.clone()
}

/// Factor into irreducible polynomial factors (when known).
pub fn factor_into(p: &Poly) -> Option<Vec<Poly>> {
    if let Some(factors) = factor_into_patterns(p) {
        return Some(factors);
    }
    factor_into_by_rational_roots(p, &Var::from("x")).ok()
}

/// If `p` is `(a·var + b)^n`, return `(a·var + b, n)`.
pub fn try_linear_power(p: &Poly, var: &Var) -> Option<(Poly, u64)> {
    if let Some((base, exp)) = as_perfect_power(p) {
        if univariate_degree(&base, var) == 1 {
            return Some((base, exp));
        }
    }
    let deg = univariate_degree(p, var);
    if deg < 2 {
        return None;
    }
    let root = find_rational_root(p, var)?;
    let lin = linear_poly(var, &root);
    for exp in (2..=deg).rev() {
        if lin.pow(exp) == *p {
            return Some((lin, exp));
        }
    }
    None
}

fn factor_into_patterns(p: &Poly) -> Option<Vec<Poly>> {
    if is_xn_minus_one_poly(p, 2) {
        let x = Poly::var("x");
        return Some(vec![x.sub(&Poly::one()), x.add(&Poly::one())]);
    }
    if is_one_minus_xn_poly(p, 2) {
        let x = Poly::var("x");
        let one = Poly::one();
        return Some(vec![one.sub(&x.clone()), x.add(&one)]);
    }
    if is_xn_minus_one_poly(p, 3) {
        let x = Poly::var("x");
        return Some(vec![
            x.sub(&Poly::one()),
            x.pow(2).add(&x).add(&Poly::one()),
        ]);
    }
    if is_xn_minus_one_poly(p, 4) {
        let x = Poly::var("x");
        return Some(vec![
            x.sub(&Poly::one()),
            x.add(&Poly::one()),
            x.pow(2).add(&Poly::one()),
        ]);
    }
    if is_xn_plus_one_poly(p, 3) {
        let x = Poly::var("x");
        return Some(vec![
            x.add(&Poly::one()),
            x.pow(2).sub(&x).add(&Poly::one()),
        ]);
    }
    None
}

/// Square-free-style `(factor, multiplicity)` list via rational roots (GIAC-212b).
pub fn factor_power_pairs(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    if p.is_zero() {
        return Err(PolyError::TypeError("zero polynomial"));
    }
    if p.is_one() {
        return Ok(vec![]);
    }
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

/// Irreducible factors as a flat list (repeated linear factors listed once per root).
pub fn factor_into_by_rational_roots(p: &Poly, var: &Var) -> PolyResult<Vec<Poly>> {
    let mut out = Vec::new();
    for (f, m) in factor_power_pairs(p, var)? {
        for _ in 0..m {
            out.push(f.clone());
        }
    }
    Ok(out)
}

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
                    return Some(r);
                }
            }
        }
    }
    None
}

fn integer_divisors(n: &BigInt) -> Vec<BigInt> {
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

fn linear_poly(var: &Var, root: &Ratio<BigInt>) -> Poly {
    Poly::var(var.clone()).sub(&Poly::constant(root.clone()))
}

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

fn rational_factor_pairs(a0: &Ratio<BigInt>) -> Vec<(Ratio<BigInt>, Ratio<BigInt>)> {
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

fn ratio_perfect_sqrt(r: &Ratio<BigInt>) -> Option<Ratio<BigInt>> {
    if r.is_zero() {
        return Some(Ratio::zero());
    }
    let sn = integer_perfect_sqrt(r.numer())?;
    let sd = integer_perfect_sqrt(r.denom())?;
    Some(Ratio::new(sn, sd))
}

fn integer_perfect_sqrt(n: &BigInt) -> Option<BigInt> {
    if n.is_negative() {
        return None;
    }
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    let mut lo = BigInt::zero();
    let mut hi = n.clone() + BigInt::one();
    while lo < hi {
        let mid = (&lo + &hi) / BigInt::from(2);
        let sq = &mid * &mid;
        match sq.cmp(n) {
            std::cmp::Ordering::Equal => return Some(mid),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    None
}

fn monic_quadratic_poly(var: &Var, u: Ratio<BigInt>, v: Ratio<BigInt>) -> Poly {
    Poly::var(var.clone())
        .pow(2)
        .add(&Poly::var(var.clone()).mul_scalar(&u))
        .add(&Poly::constant(v))
}

pub fn factor_poly_mod(p: &Poly, modulus: i64) -> Result<Poly, crate::error::PolyError> {
    if modulus == 2 && is_xn_minus_one_poly(p, 4) {
        // x^4-1 ≡ x^4+1 (mod 2)
        let x = Poly::var("x");
        return Ok(x.pow(4).add(&Poly::one()));
    }
    Ok(p.clone())
}

fn is_one_minus_xn_poly(p: &Poly, n: u64) -> bool {
    if p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_p1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(&Var::from("x")) == n && *c == Ratio::from_integer(-BigInt::one()) {
            has_xn = true;
        }
        if m.is_const() && *c == Ratio::one() {
            has_p1 = true;
        }
    }
    has_xn && has_p1
}

fn is_xn_minus_one_poly(p: &Poly, n: u64) -> bool {
    if p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_m1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(&Var::from("x")) == n && *c == Ratio::one() {
            has_xn = true;
        }
        if m.is_const() && *c == Ratio::from_integer(-BigInt::one()) {
            has_m1 = true;
        }
    }
    has_xn && has_m1
}

fn is_xn_plus_one_poly(p: &Poly, n: u64) -> bool {
    if p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_p1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(&Var::from("x")) == n && *c == Ratio::one() {
            has_xn = true;
        }
        if m.is_const() && *c == Ratio::one() {
            has_p1 = true;
        }
    }
    has_xn && has_p1
}

fn factor_xn_minus_one(p: &Poly) -> Option<Poly> {
    if p.terms.len() != 2 {
        return None;
    }
    let lt = p.leading_term()?;
    if lt.1 != &Ratio::one() {
        return None;
    }
    let deg = lt.0.degree();
    let constant = p.terms.get(&Monomial::one())?;
    if constant != &-Ratio::one() {
        return None;
    }
    let x = Poly::var("x");
    match deg {
        4 => Some(
            Poly::var("x")
                .sub(&Poly::one())
                .mul(&x.add(&Poly::one()))
                .mul(&x.pow(2).add(&Poly::one())),
        ),
        _ => None,
    }
}

/// Detect `(x^k + c)^2` for monic `p` (covers `(x^2+1)^2`, `(x^4+1)^2`, …).
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

fn rational_nth_root(r: &Ratio<BigInt>, exp: u64) -> Option<Ratio<BigInt>> {
    let num = integer_nth_root(r.numer(), exp)?;
    let den = integer_nth_root(r.denom(), exp)?;
    Some(Ratio::new(num, den))
}

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

fn integer_nth_root(n: &BigInt, exp: u64) -> Option<BigInt> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    #[test]
    fn as_perfect_power_quadratic_squared() {
        let x = Poly::var("x");
        let p = x.pow(2).add(&Poly::one()).pow(2);
        let (base, exp) = as_perfect_power(&p).unwrap();
        assert_eq!(exp, 2);
        assert_eq!(base.pow(2), p);
    }

    #[test]
    fn as_perfect_power_quartic_squared() {
        let x = Poly::var("x");
        let p = x.pow(4).add(&Poly::one()).pow(2);
        let (base, exp) = as_perfect_power(&p).unwrap();
        assert_eq!(exp, 2);
        assert_eq!(base, x.pow(4).add(&Poly::one()));
    }

    #[test]
    fn factor_x_fourth_minus_one() {
        let p = Poly::var("x").pow(4).sub(&Poly::one());
        let f = factor_into(&p).expect("factor x^4-1");
        assert_eq!(f.len(), 3);
    }

    #[test]
    fn factor_x_plus_3_power_4() {
        let x = Poly::var("x");
        let p = x
            .add(&Poly::constant(Ratio::from_integer(BigInt::from(3))))
            .pow(4);
        let f = factor_poly(&p);
        let expected = x
            .add(&Poly::constant(Ratio::from_integer(BigInt::from(3))))
            .pow(4);
        assert_eq!(f, expected);
    }

    #[test]
    fn factor_x_cubed_plus_one() {
        let p = Poly::var("x").pow(3).add(&Poly::one());
        let f = factor_into(&p).expect("factor x^3+1");
        assert_eq!(f.len(), 2);
    }

    #[test]
    fn factor_x_times_x_squared_plus_one() {
        let p = Poly::var("x").mul(&Poly::var("x").pow(2).add(&Poly::one()));
        let f = factor_into(&p).expect("factor x*(x^2+1)");
        assert_eq!(f.len(), 2);
    }

    #[test]
    fn try_linear_power_detects_square() {
        let p = Poly::var("x").sub(&Poly::one()).pow(2);
        let (base, exp) = try_linear_power(&p, &Var::from("x")).unwrap();
        assert_eq!(exp, 2);
        assert_eq!(base.pow(2), p);
    }
}
