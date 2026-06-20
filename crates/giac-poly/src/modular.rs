//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::error::{EvalError, PolyResult};
use crate::modint::ModInt;
use crate::monomial::Monomial;
use crate::poly::Poly;

/// Polynomial over ℤ/pℤ.
#[derive(Clone, Debug, PartialEq)]
pub struct PolyMod {
    pub terms: BTreeMap<Monomial, ModInt>,
    pub modulus: BigInt,
}

impl PolyMod {
    /// **Stable** — Poly zero
    pub fn zero(modulus: BigInt) -> Self {
        Self {
            terms: BTreeMap::new(),
            modulus,
        }
    }

    /// **Stable** — Poly one
    pub fn one(modulus: BigInt) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(
            Monomial::one(),
            ModInt::new(BigInt::one(), modulus.clone()).unwrap(),
        );
        Self { terms, modulus }
    }

    /// **Stable** — `from_poly`
    pub fn from_poly(p: &Poly, modulus: BigInt) -> PolyResult<Self> {
        let mut terms = BTreeMap::new();
        for (m, c) in &p.terms {
            if !c.denom().is_one() {
                return Err(EvalError::TypeError("non-integer coeff for mod poly"));
            }
            let mi = ModInt::new(c.numer().clone(), modulus.clone())?;
            if !mi.is_zero() {
                terms.insert(m.clone(), mi);
            }
        }
        Ok(Self { terms, modulus })
    }

    /// **Stable** — Poly is zero
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// **Stable** — leading term by total degree
    pub fn leading_term(&self) -> Option<(&Monomial, &ModInt)> {
        self.terms.iter().next_back()
    }

    /// **Stable** — Poly addition
    pub fn add(&self, other: &Self) -> PolyResult<Self> {
        let mut terms = self.terms.clone();
        for (m, c) in &other.terms {
            let entry = if let Some(e) = terms.get_mut(m) {
                e.add(c)?
            } else {
                c.clone()
            };
            if entry.is_zero() {
                terms.remove(m);
            } else {
                terms.insert(m.clone(), entry);
            }
        }
        Ok(Self {
            terms,
            modulus: self.modulus.clone(),
        })
    }

    /// **Stable** — Poly subtraction
    pub fn sub(&self, other: &Self) -> PolyResult<Self> {
        let mut terms = self.terms.clone();
        for (m, c) in &other.terms {
            let neg = ModInt::new(-&c.val, self.modulus.clone())?;
            let entry = if let Some(e) = terms.get_mut(m) {
                e.add(&neg)?
            } else {
                neg
            };
            if entry.is_zero() {
                terms.remove(m);
            } else {
                terms.insert(m.clone(), entry);
            }
        }
        Ok(Self {
            terms,
            modulus: self.modulus.clone(),
        })
    }

    /// **Stable** — Poly multiplication
    pub fn mul(&self, other: &Self) -> PolyResult<Self> {
        let mut out = BTreeMap::new();
        for (m1, c1) in &self.terms {
            for (m2, c2) in &other.terms {
                let m = m1.mul(m2);
                let c = c1.mul(c2)?;
                let entry = out
                    .remove(&m)
                    .unwrap_or_else(|| ModInt::new(BigInt::zero(), self.modulus.clone()).unwrap())
                    .add(&c)?;
                if !entry.is_zero() {
                    out.insert(m, entry);
                }
            }
        }
        Ok(Self {
            terms: out,
            modulus: self.modulus.clone(),
        })
    }

    /// **Stable** — multivariate division with remainder
    pub fn div_rem(&self, divisor: &Self) -> PolyResult<(Self, Self)> {
        if divisor.is_zero() {
            return Err(EvalError::DivisionByZero);
        }
        let mut remainder = self.clone();
        let mut quotient = Self::zero(self.modulus.clone());
        let Some((div_lt, div_lc)) = divisor.leading_term() else {
            return Ok((quotient, remainder));
        };
        let div_lc_inv = div_lc.inv()?;

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
            let q_c = r_lc.mul(&div_lc_inv)?;
            let q_term = Self {
                terms: [(q_m, q_c)].into(),
                modulus: self.modulus.clone(),
            };
            quotient = quotient.add(&q_term)?;
            remainder = remainder.sub(&q_term.mul(divisor)?)?;
        }
        Ok((quotient, remainder))
    }

    /// **Stable** — Poly gcd via subresultant
    pub fn gcd(&self, other: &Self) -> PolyResult<Self> {
        if self.is_zero() {
            return Ok(other.clone());
        }
        if other.is_zero() {
            return Ok(self.clone());
        }
        let mut a = self.clone();
        let mut b = other.clone();
        loop {
            if b.is_zero() {
                return Ok(a);
            }
            let (_, r) = a.div_rem(&b)?;
            a = b;
            b = r;
        }
    }
}

/// **Stable** — Poly → PolyMod mod p
pub fn modp(poly: &Poly, p: i64) -> PolyResult<PolyMod> {
    PolyMod::from_poly(poly, BigInt::from(p))
}
