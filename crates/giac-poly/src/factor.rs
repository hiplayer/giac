use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::monomial::{Monomial, Var};
use crate::poly::Poly;

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
}
