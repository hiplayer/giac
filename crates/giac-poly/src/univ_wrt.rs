//! Univariate division in **K[var]** for generic coefficient ring [`PolyCoeff`].
//!
//! **Ring context:** flat univariate over a **field** K (coefficients are scalars in K, not
//! nested polynomials). Requires exact Euclidean division: `deg r < deg b` or `r = 0`.
//!
//! **Upstream:** `_EXT` coefficient `quo`/`rem` in `gausspol.cc` (flat univariate over K).
//!
//! **Not** for nested ℚ[others][main] — use [`crate::subresultant::univariate_div_rem_wrt`]
//! there (leading-term quotient; `deg(d)=0` usually **does not** divide).
//! See `.doc/issues/GIAC-poly-flat-field-division-layering.md` §4.2.

use crate::error::{EvalError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::poly_coeff::{FieldCoeff, PolyCoeff};

/// **Stable** — whether `p` uses only powers of `var`.
pub fn is_univariate_in<C: PolyCoeff>(p: &Poly<C>, var: &Var) -> bool {
    p.terms
        .keys()
        .all(|m| m.iter().all(|(v, _)| v == var))
}

/// **Stable** — coefficient of `var^exp` as a scalar in K (0 if absent).
pub fn scalar_coeff_wrt<C: PolyCoeff>(p: &Poly<C>, var: &Var, exp: u64) -> C {
    for (m, c) in &p.terms {
        if exp == 0 && m.is_const() {
            return c.clone();
        }
        if m.exp_of(var) == exp && m.iter().all(|(v, _)| v == var) {
            return c.clone();
        }
    }
    C::coeff_zero()
}

// **Pipeline private** — `term_with_var`
fn term_with_var<C: PolyCoeff>(coeff: &C, var: &Var, exp: u64) -> PolyResult<Poly<C>> {
    if exp == 0 {
        return Ok(Poly::ring_constant(coeff.clone()));
    }
    Poly::ring_var(var.clone())
        .try_pow(exp)?
        .try_mul(&Poly::ring_constant(coeff.clone()))
}

// **Pipeline private** — debug-only Euclidean postcondition
#[cfg(debug_assertions)]
fn debug_assert_euclidean_post<C: PolyCoeff>(
    divisor: &Poly<C>,
    remainder: &Poly<C>,
    var: &Var,
    deg_b: u64,
) {
    if divisor.is_zero() || deg_b == 0 && !divisor.is_zero() {
        return;
    }
    let lc_b = scalar_coeff_wrt(divisor, var, deg_b);
    if lc_b.coeff_is_zero() {
        return;
    }
    debug_assert!(
        remainder.is_zero() || remainder.degree_wrt(var) < deg_b,
        "Euclidean invariant violated: deg(r) >= deg(b)"
    );
}

/// **Stable** — Euclidean `(q, r)` with `a = q*b + r` in K[var].
///
/// **Pre:** K is a field (`coeff_div` exact on nonzero divisors). Prefer [`crate::nested::FlatUni`]
/// for typed flat context. Not for nested ℚ[others][var] — see [`crate::subresultant::univariate_div_rem_wrt`].
///
/// **Errors:** `TypeError` if `b = 0`, or `lc(b) = 0` with `deg(b) > 0`, or nonzero constant
/// divisor with zero scalar coefficient.
///
/// **Post (`b` valid):** `a = q*b + r`; `r = 0` or `deg(r) < deg(b)`; if `deg(b)=0` and `b ≠ 0`
/// then `r = 0`.
///
/// | Condition | flat (this fn) | nested [`subresultant::univariate_div_rem_wrt`] |
/// |-----------|----------------|--------------------------------------------------|
/// | `b = 0` | `Err` | caller must avoid |
/// | `deg(b)=0`, b≠0 | `r=0`, exact `q=a/b` | usually `(0, a)` |
/// | `lc(b)=0`, deg>0 | `Err` | pseudo / break |
pub fn univariate_div_rem_wrt<C: PolyCoeff>(
    a: &Poly<C>,
    b: &Poly<C>,
    var: &Var,
) -> PolyResult<(Poly<C>, Poly<C>)> {
    let mut remainder = a.clone();
    let mut quotient = Poly::ring_zero();
    let db = b.degree_wrt(var);
    if db == 0 {
        if b.is_zero() {
            return Err(EvalError::TypeError("division by zero polynomial"));
        }
        if remainder.is_zero() {
            return Ok((quotient, Poly::ring_zero()));
        }
        // K is a field: nonzero constant divisor divides exactly.
        let bc = scalar_coeff_wrt(b, var, 0);
        if bc.coeff_is_zero() {
            return Err(EvalError::TypeError("zero constant divisor"));
        }
        let inv = C::coeff_one().coeff_div(&bc)?;
        let deg = remainder.degree_wrt(var);
        for exp in 0..=deg {
            let c = scalar_coeff_wrt(&remainder, var, exp);
            if c.coeff_is_zero() {
                continue;
            }
            let scaled = c.coeff_mul(&inv)?;
            quotient = quotient.try_add(&term_with_var(&scaled, var, exp)?)?;
        }
        remainder = Poly::ring_zero();
        debug_assert_euclidean_post(b, &remainder, var, db);
        return Ok((quotient, remainder));
    }
    let lc_b = scalar_coeff_wrt(b, var, db);
    if lc_b.coeff_is_zero() {
        return Err(EvalError::TypeError("leading coefficient zero"));
    }

    loop {
        let dr = remainder.degree_wrt(var);
        if dr < db || remainder.is_zero() {
            break;
        }
        let lc_r = scalar_coeff_wrt(&remainder, var, dr);
        let q_coeff = lc_r.coeff_div(&lc_b)?;
        let d = dr - db;
        let q_term = term_with_var(&q_coeff, var, d)?;
        quotient = quotient.try_add(&q_term)?;
        remainder = remainder.try_sub(&q_term.try_mul(b)?)?;
    }
    debug_assert_euclidean_post(b, &remainder, var, db);
    Ok((quotient, remainder))
}

/// **Stable** — exact quotient in K[var]; fails if remainder is nonzero.
pub fn quo_exact_wrt<C: PolyCoeff>(
    num: &Poly<C>,
    den: &Poly<C>,
    var: &Var,
) -> PolyResult<Poly<C>> {
    let (q, r) = univariate_div_rem_wrt(num, den, var)?;
    if r.is_zero() {
        Ok(q)
    } else {
        Err(EvalError::NotImplemented("poly division"))
    }
}

/// **Stable** — formal derivative w.r.t. `var` in K[var].
pub fn derivative_wrt<C: PolyCoeff>(p: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    if !is_univariate_in(p, var) {
        return Err(EvalError::TypeError("not univariate"));
    }
    let deg = p.degree_wrt(var);
    let mut out = Poly::ring_zero();
    for exp in 1..=deg {
        let c = scalar_coeff_wrt(p, var, exp);
        if c.coeff_is_zero() {
            continue;
        }
        let coeff = c.coeff_mul(&scalar_coeff_from_u64(exp as u64)?)?;
        out = out.try_add(&term_with_var(&coeff, var, exp - 1)?)?;
    }
    Ok(out)
}

// **Pipeline private** — embed small integer into C via 1+1+…
fn scalar_coeff_from_u64<C: PolyCoeff>(n: u64) -> PolyResult<C> {
    let mut acc = C::coeff_zero();
    for _ in 0..n {
        acc = acc.coeff_add(&C::coeff_one())?;
    }
    Ok(acc)
}

/// **Stable** — extended gcd in K[var]; returns `(g, s, t)` with `g` monic.
pub fn egcd_wrt<C: PolyCoeff>(
    a: &Poly<C>,
    b: &Poly<C>,
    var: &Var,
) -> PolyResult<(Poly<C>, Poly<C>, Poly<C>)> {
    let mut old_r = a.clone();
    let mut r = b.clone();
    let mut old_s = Poly::ring_one();
    let mut s = Poly::ring_zero();
    let mut old_t = Poly::ring_zero();
    let mut t = Poly::ring_one();

    while !r.is_zero() {
        let (q, new_r) = univariate_div_rem_wrt(&old_r, &r, var)?;
        old_r = r;
        r = new_r;
        let new_s = old_s.try_sub(&q.try_mul(&s)?)?;
        old_s = s;
        s = new_s;
        let new_t = old_t.try_sub(&q.try_mul(&t)?)?;
        old_t = t;
        t = new_t;
    }
    Ok((monic_wrt(&old_r, var)?, old_s, old_t))
}

/// **Stable** — gcd in K[var] (monic).
pub fn gcd_wrt<C: PolyCoeff>(a: &Poly<C>, b: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    Ok(egcd_wrt(a, b, var)?.0)
}

/// **Stable** — gcd of scalar coefficients (content in K).
pub fn content_scalars<C: PolyCoeff>(p: &Poly<C>) -> PolyResult<C> {
    let mut g: Option<C> = None;
    for c in p.terms.values() {
        if c.coeff_is_zero() {
            continue;
        }
        g = Some(match g {
            None => c.clone(),
            Some(prev) => coeff_gcd(&prev, c)?,
        });
    }
    Ok(g.unwrap_or_else(C::coeff_zero))
}

/// **Stable** — content w.r.t. `var` (univariate; equals scalar content in K).
pub fn content_wrt<C: PolyCoeff>(p: &Poly<C>, var: &Var) -> PolyResult<C> {
    if !is_univariate_in(p, var) {
        return Err(EvalError::TypeError("not univariate"));
    }
    content_scalars(p)
}

/// **Stable** — primitive part w.r.t. `var` in K[var].
pub fn primitive_part_wrt<C: PolyCoeff>(p: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    let c = content_wrt(p, var)?;
    if c.coeff_is_zero() {
        return Ok(p.clone());
    }
    if c.coeff_is_one() {
        return Ok(p.clone());
    }
    let inv = C::coeff_one().coeff_div(&c)?;
    let deg = p.degree_wrt(var);
    let mut out = Poly::ring_zero();
    for exp in 0..=deg {
        let coeff = scalar_coeff_wrt(p, var, exp);
        if coeff.coeff_is_zero() {
            continue;
        }
        let scaled = coeff.coeff_mul(&inv)?;
        out = out.try_add(&term_with_var(&scaled, var, exp)?)?;
    }
    Ok(out)
}

/// **Stable** — square-free factorization `p = ∏ g_k^k` in K[var] (Yun).
pub fn square_free_factorization_wrt<C: FieldCoeff>(
    p: &Poly<C>,
    var: &Var,
) -> PolyResult<Vec<(Poly<C>, usize)>> {
    if !is_univariate_in(p, var) {
        return Err(EvalError::TypeError("not univariate"));
    }
    let ring = crate::square_free::FlatCoeffVarRing::new(var);
    crate::square_free::square_free_yun(&ring, p)
}

/// **Stable** — square-free part w.r.t. `var` in K[var].
pub fn square_free_part_wrt<C: FieldCoeff>(p: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    let ring = crate::square_free::FlatCoeffVarRing::new(var);
    let mut prod = Poly::ring_one();
    for (g, _) in crate::square_free::square_free_yun(&ring, p)? {
        prod = prod.try_mul(&g)?;
    }
    Ok(prod)
}

/// **Stable** — `(a, b, c)` for univariate quadratic `a·var² + b·var + c`.
pub fn quadratic_coeffs_wrt<C: PolyCoeff>(
    p: &Poly<C>,
    var: &Var,
) -> Option<(C, C, C)> {
    if p.degree_wrt(var) != 2 {
        return None;
    }
    if !is_univariate_in(p, var) {
        return None;
    }
    let a = scalar_coeff_wrt(p, var, 2);
    if a.coeff_is_zero() {
        return None;
    }
    Some((a, scalar_coeff_wrt(p, var, 1), scalar_coeff_wrt(p, var, 0)))
}

// **Pipeline private** — gcd in coefficient ring K
fn coeff_gcd<C: PolyCoeff>(a: &C, b: &C) -> PolyResult<C> {
    if a.coeff_is_zero() {
        return Ok(b.clone());
    }
    if b.coeff_is_zero() {
        return Ok(a.clone());
    }
    let mut x = a.clone();
    let mut y = b.clone();
    loop {
        if y.coeff_is_zero() {
            break;
        }
        match x.coeff_div(&y) {
            Ok(q) => {
                let qy = q.coeff_mul(&y)?;
                x = x.coeff_sub(&qy)?;
                std::mem::swap(&mut x, &mut y);
            }
            Err(_) => {
                // K is a field: coprime scalars → unit content
                return Ok(C::coeff_one());
            }
        }
    }
    if x.coeff_is_zero() {
        Ok(C::coeff_zero())
    } else {
        Ok(C::coeff_one())
    }
}

/// **Stable** — divide by leading coefficient w.r.t. `var` (monic in K).
pub fn monic_wrt<C: PolyCoeff>(p: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    let deg = p.degree_wrt(var);
    let lc = scalar_coeff_wrt(p, var, deg);
    if lc.coeff_is_zero() {
        return Err(EvalError::TypeError("leading coefficient zero"));
    }
    if !is_univariate_in(p, var) {
        return Err(EvalError::TypeError("not univariate"));
    }
    let inv = C::coeff_one().coeff_div(&lc)?;
    let mut out = Poly::ring_zero();
    for exp in 0..=deg {
        let c = scalar_coeff_wrt(p, var, exp);
        if c.coeff_is_zero() {
            continue;
        }
        let scaled = c.coeff_mul(&inv)?;
        out = out.try_add(&term_with_var(&scaled, var, exp)?)?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use num_rational::Ratio;
    use num_traits::One;

    use super::*;
    use crate::monomial::Var;
    use crate::poly::Poly;

    fn x() -> Var {
        Var::from("x")
    }

    #[test]
    fn div_rem_rational_x_squared_minus_one() {
        let p = Poly::var(x()).try_pow(2).unwrap().try_sub(&Poly::one()).unwrap();
        let d = Poly::var(x()).try_sub(&Poly::one()).unwrap();
        let (q, r) = univariate_div_rem_wrt(&p, &d, &x()).unwrap();
        assert!(r.is_zero());
        assert_eq!(q, Poly::var(x()).try_add(&Poly::one()).unwrap());
    }

    #[test]
    fn monic_wrt_makes_leading_one() {
        let two = Ratio::from_integer(2.into());
        let four = Ratio::from_integer(4.into());
        let p = Poly::constant(two)
            .try_mul(&Poly::var(x()).try_pow(2).unwrap())
            .unwrap()
            .try_sub(&Poly::constant(four))
            .unwrap();
        let m = monic_wrt(&p, &x()).unwrap();
        let lc = scalar_coeff_wrt(&m, &x(), 2);
        assert!(lc.coeff_is_one());
    }

    #[test]
    fn gcd_wrt_rational_coprimality() {
        let p = Poly::var(x()).try_pow(2).unwrap().try_add(&Poly::one()).unwrap();
        let d = Poly::var(x()).try_sub(&Poly::one()).unwrap();
        assert!(gcd_wrt(&p, &d, &x()).unwrap().is_one());
        let p2 = Poly::var(x()).try_pow(2).unwrap().try_sub(&Poly::one()).unwrap();
        let d2 = Poly::var(x()).try_sub(&Poly::one()).unwrap();
        assert_eq!(gcd_wrt(&p2, &d2, &x()).unwrap(), d2);
    }

    /// Constant-divisor Euclidean step (P1 regression): rem(x²−2, x+1) = −1, then gcd → 1.
    #[test]
    fn gcd_wrt_rational_x_squared_minus_2_and_x_plus_one() {
        let two = Ratio::from_integer(2.into());
        let p = Poly::var(x())
            .try_pow(2)
            .unwrap()
            .try_sub(&Poly::constant(two))
            .unwrap();
        let d = Poly::var(x()).try_add(&Poly::one()).unwrap();
        assert!(gcd_wrt(&p, &d, &x()).unwrap().is_one());
    }

    #[test]
    fn div_rem_by_zero_is_error() {
        let p = Poly::var(x()).try_pow(2).unwrap().try_sub(&Poly::one()).unwrap();
        let zero = Poly::zero();
        assert!(univariate_div_rem_wrt(&p, &zero, &x()).is_err());
    }

    #[test]
    fn div_rem_constant_divisor_exact() {
        let p = Poly::var(x()).try_add(&Poly::one()).unwrap();
        let two = Poly::constant(Ratio::from_integer(2.into()));
        let (q, r) = univariate_div_rem_wrt(&p, &two, &x()).unwrap();
        assert!(r.is_zero());
        let half = Ratio::from_integer(1.into()) / Ratio::from_integer(2.into());
        assert_eq!(scalar_coeff_wrt(&q, &x(), 1), half);
    }
}
