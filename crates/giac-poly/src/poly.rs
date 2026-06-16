use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::{Monomial, Var};

/// Sparse multivariate polynomial over ℚ.
#[derive(Clone, Debug, PartialEq)]
pub struct Poly {
    pub terms: BTreeMap<Monomial, Ratio<BigInt>>,
}

impl Poly {
    pub fn zero() -> Self {
        Self {
            terms: BTreeMap::new(),
        }
    }

    pub fn one() -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(Monomial::one(), Ratio::one());
        Self { terms }
    }

    pub fn constant(c: Ratio<BigInt>) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        let mut terms = BTreeMap::new();
        terms.insert(Monomial::one(), c);
        Self { terms }
    }

    pub fn var(name: impl Into<Var>) -> Self {
        Self {
            terms: [(Monomial::var(name), Ratio::one())].into(),
        }
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn is_one(&self) -> bool {
        self.terms.len() == 1
            && self
                .terms
                .get(&Monomial::one())
                .is_some_and(|c| c.is_one())
    }

    pub fn leading_term(&self) -> Option<(&Monomial, &Ratio<BigInt>)> {
        self.terms.iter().next_back()
    }

    /// Leading term in lex order induced by `var_order` (first variable is greatest).
    pub fn leading_term_lex(&self, var_order: &[Var]) -> Option<(&Monomial, &Ratio<BigInt>)> {
        self.terms.iter().max_by(|(m1, _), (m2, _)| m1.cmp_lex(m2, var_order))
    }

    pub fn term(monom: Monomial, coeff: Ratio<BigInt>) -> Self {
        if coeff.is_zero() {
            return Self::zero();
        }
        Self {
            terms: [(monom, coeff)].into(),
        }
    }

    pub fn degree(&self) -> u64 {
        self.terms.keys().map(Monomial::degree).max().unwrap_or(0)
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut terms = self.terms.clone();
        for (m, c) in &other.terms {
            let entry = terms.entry(m.clone()).or_insert_with(Ratio::zero);
            *entry += c;
            if entry.is_zero() {
                terms.remove(m);
            }
        }
        Self { terms }
    }

    pub fn sub(&self, other: &Self) -> Self {
        let mut terms = self.terms.clone();
        for (m, c) in &other.terms {
            let entry = terms.entry(m.clone()).or_insert_with(Ratio::zero);
            *entry -= c;
            if entry.is_zero() {
                terms.remove(m);
            }
        }
        Self { terms }
    }

    pub fn neg(&self) -> Self {
        let terms = self
            .terms
            .iter()
            .map(|(m, c)| (m.clone(), -c.clone()))
            .collect();
        Self { terms }
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut out = BTreeMap::new();
        for (m1, c1) in &self.terms {
            for (m2, c2) in &other.terms {
                let m = m1.mul(m2);
                *out.entry(m).or_insert_with(Ratio::zero) += c1 * c2;
            }
        }
        out.retain(|_, c| !c.is_zero());
        Self { terms: out }
    }

    pub fn mul_scalar(&self, s: &Ratio<BigInt>) -> Self {
        if s.is_zero() {
            return Self::zero();
        }
        let terms = self
            .terms
            .iter()
            .map(|(m, c)| (m.clone(), c * s))
            .collect();
        Self { terms }
    }

    pub fn pow(&self, exp: u64) -> Self {
        if exp == 0 {
            return Self::one();
        }
        let mut result = Self::one();
        let mut base = self.clone();
        let mut e = exp;
        while e > 0 {
            if e % 2 == 1 {
                result = result.mul(&base);
            }
            base = base.mul(&base);
            e /= 2;
        }
        result
    }

    pub fn content(&self) -> Ratio<BigInt> {
        let mut g = Ratio::zero();
        for c in self.terms.values() {
            if g.is_zero() {
                g = c.clone();
            } else {
                g = integer_content_gcd(&g, c);
            }
        }
        g
    }

    pub fn primitive_part(&self) -> Self {
        let c = self.content();
        if c.is_zero() || c.is_one() {
            return self.clone();
        }
        self.mul_scalar(&(&Ratio::one() / &c))
    }

    pub fn monic(&self) -> Self {
        if let Some((_, lc)) = self.leading_term() {
            if lc.is_one() {
                return self.clone();
            }
            return self.mul_scalar(&(&Ratio::one() / lc));
        }
        self.clone()
    }

    pub fn div_rem(&self, divisor: &Self) -> (Self, Self) {
        if divisor.is_zero() {
            return (Self::zero(), self.clone());
        }
        let mut remainder = self.clone();
        let mut quotient = Self::zero();
        let Some((div_lt, div_lc)) = divisor.leading_term() else {
            return (Self::zero(), self.clone());
        };

        loop {
            let Some((r_lt, r_lc)) = remainder.leading_term() else {
                break;
            };
            if r_lt.degree() < div_lt.degree() {
                break;
            }
            if !r_lt.is_dividing(div_lt) {
                break;
            }
            let Some(q_m) = r_lt.div_exact(div_lt) else {
                break;
            };
            let q_c = r_lc.clone() / div_lc.clone();
            let q_term = Poly {
                terms: [(q_m, q_c)].into(),
            };
            quotient = quotient.add(&q_term);
            remainder = remainder.sub(&q_term.mul(divisor));
        }
        (quotient, remainder)
    }

    pub fn div_exact(&self, divisor: &Self) -> Option<Self> {
        let (q, r) = self.div_rem(divisor);
        if r.is_zero() {
            Some(q)
        } else {
            None
        }
    }

    pub fn gcd(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }
        let mut a = self.primitive_part();
        let mut b = other.primitive_part();
        let mut steps = 0usize;
        const MAX_GCD_STEPS: usize = 512;
        loop {
            if b.is_zero() {
                return a.monic();
            }
            if steps >= MAX_GCD_STEPS {
                return a.monic();
            }
            steps += 1;
            let (_, r) = a.div_rem(&b);
            a = b;
            b = r;
        }
    }

    pub fn lcm(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }
        if self.is_one() {
            return other.clone();
        }
        if other.is_one() {
            return self.clone();
        }
        let g = self.gcd(other);
        self.mul(other).div_rem(&g).0.monic()
    }

    pub fn horner(&self, var: &Var, x: &Ratio<BigInt>) -> Ratio<BigInt> {
        let mut coeffs: BTreeMap<u64, Ratio<BigInt>> = BTreeMap::new();
        for (m, c) in &self.terms {
            if m.is_const() {
                *coeffs.entry(0).or_insert_with(Ratio::zero) += c;
            } else if m.exp_of(var) == m.degree() && m.iter().count() == 1 {
                let exp = m.exp_of(var);
                *coeffs.entry(exp).or_insert_with(Ratio::zero) += c;
            }
        }
        if coeffs.is_empty() {
            // fallback: evaluate by substitution in general case
            return Ratio::zero();
        }
        let max = *coeffs.keys().next_back().unwrap_or(&0);
        let mut acc = Ratio::zero();
        for e in (0..=max).rev() {
            acc = acc * x + coeffs.get(&e).cloned().unwrap_or_else(Ratio::zero);
        }
        acc
    }
}

fn integer_content_gcd(a: &Ratio<BigInt>, b: &Ratio<BigInt>) -> Ratio<BigInt> {
    let na = a.numer().abs();
    let da = a.denom().abs();
    let nb = b.numer().abs();
    let db = b.denom().abs();
    let g_num = na.gcd(&nb);
    let g_den = da.lcm(&db);
    Ratio::new(g_num, g_den)
}

pub fn quo(a: &Poly, b: &Poly) -> PolyResult<Poly> {
    if b.is_zero() {
        return Err(PolyError::DivisionByZero);
    }
    Ok(a.div_rem(b).0)
}

pub fn rem(a: &Poly, b: &Poly) -> PolyResult<Poly> {
    if b.is_zero() {
        return Err(PolyError::DivisionByZero);
    }
    Ok(a.div_rem(b).1)
}

pub fn egcd(a: &Poly, b: &Poly) -> (Poly, Poly, Poly) {
    let mut old_r = a.clone();
    let mut r = b.clone();
    let mut old_s = Poly::one();
    let mut s = Poly::zero();
    let mut old_t = Poly::zero();
    let mut t = Poly::one();

    while !r.is_zero() {
        let (q, new_r) = old_r.div_rem(&r);
        old_r = r;
        r = new_r;
        let new_s = old_s.sub(&q.mul(&s));
        old_s = s;
        s = new_s;
        let new_t = old_t.sub(&q.mul(&t));
        old_t = t;
        t = new_t;
    }
    (old_r.monic(), old_s, old_t)
}

pub fn simp2(num: &Poly, den: &Poly) -> (Poly, Poly) {
    let g = num.gcd(den);
    (num.div_rem(&g).0, den.div_rem(&g).0)
}

/// Extended gcd solving a*u + b*v = c when c divides gcd(a,b).
pub fn abcuv(a: &Poly, b: &Poly, c: &Poly) -> PolyResult<(Poly, Poly)> {
    let (g, u, v) = egcd(a, b);
    let (q, r) = c.div_rem(&g);
    if !r.is_zero() {
        return Err(PolyError::TypeError("c not divisible by gcd(a,b)"));
    }
    Ok((u.mul(&q), v.mul(&q)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn gcd_x3_x2() {
        let p = x().pow(3).sub(&Poly::one());
        let q = x().pow(2).sub(&Poly::one());
        let g = p.gcd(&q);
        assert_eq!(g, x().sub(&Poly::one()));
    }

    #[test]
    fn quo_rem() {
        let p = x().pow(3).sub(&Poly::one());
        let d = x().sub(&Poly::one());
        assert_eq!(quo(&p, &d).unwrap(), x().pow(2).add(&x()).add(&Poly::one()));
        assert!(rem(&p, &d).unwrap().is_zero());
    }

    #[test]
    fn abcuv_linear_one() {
        let g = x().sub(&Poly::one());
        let gp = Poly::one();
        let a = Poly::one();
        let (u, v) = abcuv(&g, &gp, &a).unwrap();
        assert_eq!(u.mul(&g).add(&v.mul(&gp)), a);
    }
}
