//! Generic Yun square-free factorization over univariate polynomial rings.
//!
//! **Stable (bounded):** [`square_free_yun`] (char 0), [`square_free_yun_mod`] (finite field).
//! **Pipeline private:** [`SquareFreeRing`], [`UniVarRing`].

use crate::error::{EvalError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::univariate_degree;
use crate::univariate::{univariate_derivative, univariate_div_exact, univariate_gcd};

/// Ring operations for Yun square-free factorization.
pub(crate) trait SquareFreeRing {
    type Poly: Clone;

    fn is_zero(&self, p: &Self::Poly) -> bool;
    fn is_one(&self, p: &Self::Poly) -> bool;
    fn derivative(&self, p: &Self::Poly) -> PolyResult<Self::Poly>;
    fn gcd(&self, a: &Self::Poly, b: &Self::Poly) -> PolyResult<Self::Poly>;
    fn div_exact(&self, a: &Self::Poly, b: &Self::Poly) -> PolyResult<Self::Poly>;
    fn sub(&self, a: &Self::Poly, b: &Self::Poly) -> PolyResult<Self::Poly>;
    fn max_exponent(&self, p: &Self::Poly) -> usize;
}

// **Pipeline private** — divide `w,y` by gcd(w,y); return that gcd
pub(crate) fn gcd_reduce<R: SquareFreeRing>(
    ring: &R,
    w: &mut R::Poly,
    y: &mut R::Poly,
) -> PolyResult<R::Poly> {
    let g = ring.gcd(w, y)?;
    if !ring.is_one(&g) {
        *w = ring.div_exact(w, &g)?;
        *y = ring.div_exact(y, &g)?;
    }
    Ok(g)
}

/// Char-0 Yun: `p = ∏ f_k^k` (giac `Tsqff_char0`).
pub(crate) fn square_free_yun<R: SquareFreeRing>(
    ring: &R,
    p: &R::Poly,
) -> PolyResult<Vec<(R::Poly, usize)>> {
    if ring.is_zero(p) {
        return Err(EvalError::TypeError("zero polynomial"));
    }
    let mut w = p.clone();
    let mut y = ring.derivative(p)?;
    gcd_reduce(ring, &mut w, &mut y)?;
    y = ring.sub(&y, &ring.derivative(&w)?)?;

    let mut factors = Vec::new();
    let mut k = 1usize;
    let max_k = ring.max_exponent(p) + 2;
    while !ring.is_zero(&y) && k <= max_k {
        let g = gcd_reduce(ring, &mut w, &mut y)?;
        if !ring.is_one(&g) {
            factors.push((g, k));
        }
        y = ring.sub(&y, &ring.derivative(&w)?)?;
        k += 1;
    }
    if !ring.is_one(&w) {
        factors.push((w, k));
    }
    Ok(factors)
}

/// Finite-field Yun (char *p* loop; used by `factor/fpx`).
pub(crate) fn square_free_yun_mod<R: SquareFreeRing>(
    ring: &R,
    p: &R::Poly,
    degree: impl Fn(&R::Poly) -> u64,
) -> PolyResult<Vec<(R::Poly, usize)>> {
    if ring.is_zero(p) {
        return Err(EvalError::TypeError("zero polynomial"));
    }
    let mut factors = Vec::new();
    let mut w = p.clone();
    let mut y = ring.derivative(p)?;
    gcd_reduce(ring, &mut w, &mut y)?;
    y = ring.sub(&y, &ring.derivative(&w)?)?;
    let mut k = 1usize;
    loop {
        if degree(&w) == 0 {
            break;
        }
        let g = ring.gcd(&w, &y)?;
        if !ring.is_one(&g) {
            factors.push((g.clone(), k));
            w = ring.div_exact(&w, &g)?;
        }
        if degree(&w) == 0 {
            break;
        }
        k += 1;
        y = ring.div_exact(&y, &g)?;
    }
    if !ring.is_one(&w) {
        factors.push((w, k));
    }
    Ok(factors)
}

/// ℚ[x] viewed univariately in `var` (coefficients in ℚ).
pub(crate) struct UniVarRing<'a> {
    pub var: &'a Var,
}

impl SquareFreeRing for UniVarRing<'_> {
    type Poly = Poly;

    fn is_zero(&self, p: &Poly) -> bool {
        p.is_zero()
    }

    fn is_one(&self, p: &Poly) -> bool {
        p.is_one()
    }

    fn derivative(&self, p: &Poly) -> PolyResult<Poly> {
        Ok(univariate_derivative(p, self.var))
    }

    fn gcd(&self, a: &Poly, b: &Poly) -> PolyResult<Poly> {
        Ok(univariate_gcd(a, b, self.var))
    }

    fn div_exact(&self, a: &Poly, b: &Poly) -> PolyResult<Poly> {
        univariate_div_exact(a, b, self.var).ok_or(EvalError::NotImplemented("poly division"))
    }

    fn sub(&self, a: &Poly, b: &Poly) -> PolyResult<Poly> {
        Ok(a.sub(b))
    }

    fn max_exponent(&self, p: &Poly) -> usize {
        univariate_degree(p, self.var) as usize
    }
}
