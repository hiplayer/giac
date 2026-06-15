use std::collections::BTreeMap;
use std::sync::Arc;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::{EvalError, Expr, ExprArc, Ident};

/// Exponent vector sorted by variable name.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Monomial(BTreeMap<Ident, u32>);

impl Monomial {
    pub fn one() -> Self {
        Self(BTreeMap::new())
    }

    pub fn var(id: Ident) -> Self {
        let mut m = BTreeMap::new();
        m.insert(id, 1);
        Self(m)
    }

    pub fn degree(&self) -> u32 {
        self.0.values().sum()
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut out = self.0.clone();
        for (v, e) in &other.0 {
            *out.entry(v.clone()).or_insert(0) += e;
        }
        Self(out)
    }

    pub fn div_exact(&self, other: &Self) -> Option<Self> {
        let mut out = self.0.clone();
        for (v, e) in &other.0 {
            let entry = out.get_mut(v)?;
            if *entry < *e {
                return None;
            }
            *entry -= e;
            if *entry == 0 {
                out.remove(v);
            }
        }
        Some(Self(out))
    }

    pub fn is_dividing(&self, other: &Self) -> bool {
        other.0.iter().all(|(v, e)| self.0.get(v).copied().unwrap_or(0) >= *e)
    }
}

/// Sparse multivariate polynomial with rational coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct Poly {
    pub(crate) terms: BTreeMap<Monomial, Ratio<BigInt>>,
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

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn leading_term(&self) -> Option<(&Monomial, &Ratio<BigInt>)> {
        self.terms.iter().next_back()
    }

    pub fn degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|m| m.degree())
            .max()
            .unwrap_or(0)
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

    pub fn pow(&self, exp: u32) -> Self {
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

    pub fn gcd(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }
        let mut a = self.clone();
        let mut b = other.clone();
        loop {
            if b.is_zero() {
                return a.monic();
            }
            let (_, r) = a.div_rem(&b);
            a = b;
            b = r;
        }
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
        let div_lt = divisor.leading_term().unwrap().0.clone();
        let div_lc = divisor.leading_term().unwrap().1.clone();

        loop {
            let Some((r_lt, r_lc)) = remainder.leading_term() else {
                break;
            };
            if r_lt.degree() < div_lt.degree() {
                break;
            }
            if !r_lt.is_dividing(&div_lt) {
                break;
            }
            let q_m = r_lt.div_exact(&div_lt).unwrap();
            let q_c = r_lc.clone() / div_lc.clone();
            let q_term = Poly {
                terms: [(q_m, q_c)].into(),
            };
            quotient = quotient.add(&q_term);
            remainder = remainder.sub(&q_term.mul(divisor));
        }
        (quotient, remainder)
    }

    /// Least common multiple for polynomial denominators.
    pub fn lcm(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }
        if *self == Self::one() {
            return other.clone();
        }
        if *other == Self::one() {
            return self.clone();
        }
        let g = self.gcd(other);
        self.mul(other).div_rem(&g).0
    }

    /// Divide exactly when remainder is zero.
    pub fn div_exact(&self, divisor: &Self) -> Option<Self> {
        let (q, r) = self.div_rem(divisor);
        if r.is_zero() {
            Some(q)
        } else {
            None
        }
    }
}

/// Try to convert an expression to a polynomial.
pub fn expr_to_poly(expr: &Expr) -> Result<Poly, EvalError> {
    match expr {
        Expr::Int(n) => Ok(Poly::constant(Ratio::from_integer(n.clone()))),
        Expr::Rat(r) => Ok(Poly::constant(r.clone())),
        Expr::Symbol(id) => Ok(Poly {
            terms: [(Monomial::var(id.clone()), Ratio::one())].into(),
        }),
        Expr::Add(terms) => terms
            .iter()
            .map(|t| expr_to_poly(t))
            .try_fold(Poly::zero(), |acc, p| Ok(acc.add(&p?))),
        Expr::Mul(factors) => factors
            .iter()
            .map(|f| expr_to_poly(f))
            .try_fold(Poly::one(), |acc, p| Ok(acc.mul(&p?))),
        Expr::Pow(base, exp) => {
            let base_p = expr_to_poly(base)?;
            if let Expr::Int(e) = exp.as_ref() {
                if e >= &BigInt::zero() && e <= &BigInt::from(50) {
                    return Ok(base_p.pow(e.to_string().parse().unwrap()));
                }
            }
            Err(EvalError::TypeError("non-polynomial power"))
        }
        _ => Err(EvalError::TypeError("not a polynomial expression")),
    }
}

pub fn poly_to_expr(poly: &Poly) -> ExprArc {
    if poly.is_zero() {
        return Expr::int(0);
    }
    let mut terms: Vec<(u32, ExprArc)> = Vec::new();
    for (m, c) in &poly.terms {
        let coeff = ratio_to_expr(c);
        let term = monomial_to_expr(m, coeff);
        terms.push((m.degree(), term));
    }
    terms.sort_by(|a, b| b.0.cmp(&a.0));
    Expr::add(terms.into_iter().map(|(_, t)| t).collect())
}

fn ratio_to_expr(r: &Ratio<BigInt>) -> ExprArc {
    if r.is_zero() {
        Expr::int(0)
    } else if r.denom() == &BigInt::one() {
        let n = r.numer();
        if n >= &BigInt::zero() {
            Expr::int(n.to_string().parse().unwrap_or(0))
        } else {
            Expr::int(n.to_string().parse().unwrap_or(0))
        }
    } else {
        Arc::new(Expr::Rat(r.clone()))
    }
}

fn monomial_to_expr(m: &Monomial, coeff: ExprArc) -> ExprArc {
    if m.0.is_empty() {
        return coeff;
    }
    let mut factors: Vec<ExprArc> = Vec::new();
    if !matches!(coeff.as_ref(), Expr::Int(n) if n.is_one()) {
        factors.push(coeff);
    }
    for (v, e) in &m.0 {
        let base = Expr::sym(v.as_str());
        if *e == 1 {
            factors.push(base);
        } else {
            factors.push(Expr::pow(base, Expr::int(*e as i64)));
        }
    }
    if factors.is_empty() {
        Expr::int(1)
    } else if factors.len() == 1 {
        Arc::clone(&factors[0])
    } else {
        Expr::mul(factors)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Expr;

    #[test]
    fn poly_add_mul() {
        let x = Ident::new("x");
        let p1 = Poly {
            terms: [(Monomial::var(x.clone()), Ratio::one())].into(),
        };
        let p2 = Poly::constant(Ratio::from_integer(BigInt::from(3)));
        let sum = p1.add(&p2);
        let e = poly_to_expr(&sum);
        match e.as_ref() {
            Expr::Add(v) => assert_eq!(v.len(), 2),
            _ => panic!("expected add"),
        }
    }

    #[test]
    fn gcd_linear() {
        let x = Ident::new("x");
        let one = |c: i64| Poly::constant(Ratio::from_integer(BigInt::from(c)));
        let x_m = |c: i64, e: u32| {
            let mut m = BTreeMap::new();
            if e > 0 {
                m.insert(x.clone(), e);
            }
            Poly {
                terms: [(Monomial(m), Ratio::from_integer(BigInt::from(c)))].into(),
            }
        };
        // x^2 - 2x + 1 = (x-1)^2
        let p = x_m(1, 2).sub(&x_m(2, 1)).add(&one(1));
        // x^3 - 1
        let q = x_m(1, 3).sub(&one(1));
        let g = p.gcd(&q);
        let g_expr = poly_to_expr(&g);
        // gcd should be x-1
        let expected = Expr::add(vec![Expr::sym("x"), Expr::int(-1)]);
        assert_eq!(g_expr, expected);
    }

    #[test]
    fn poly_lcm_and_div_exact() {
        let x = Ident::new("x");
        let p = Poly {
            terms: [(Monomial::var(x.clone()), Ratio::one())].into(),
        };
        let one = Poly::one();
        let lcm = p.lcm(&one);
        assert_eq!(lcm, p);
        let sq = p.mul(&p);
        let half = sq.div_exact(&p).unwrap();
        assert_eq!(half, p);
    }
}
