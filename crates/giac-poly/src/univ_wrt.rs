//! Univariate division in K[var] for generic coefficient ring [`PolyCoeff`].
//!
//! **Upstream:** `_EXT` coefficient `quo`/`rem` in `gausspol.cc` (flat univariate over K).
//! **Not** for nested ℚ[others][main] — use [`crate::subresultant::univariate_div_rem_wrt`] there.

use crate::error::{EvalError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::poly_coeff::PolyCoeff;

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

/// **Stable** — Euclidean `(q, r)` with `a = q*b + r` in K[var].
pub fn univariate_div_rem_wrt<C: PolyCoeff>(
    a: &Poly<C>,
    b: &Poly<C>,
    var: &Var,
) -> PolyResult<(Poly<C>, Poly<C>)> {
    let mut remainder = a.clone();
    let mut quotient = Poly::ring_zero();
    let db = b.degree_wrt(var);
    if db == 0 {
        return Ok((quotient, remainder));
    }
    let lc_b = scalar_coeff_wrt(b, var, db);
    if lc_b.coeff_is_zero() {
        return Ok((quotient, remainder));
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
}
