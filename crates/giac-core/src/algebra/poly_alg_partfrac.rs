//! Partial fractions over K for `Poly<AlgExtC>` (T2-4).
//!
//! **Upstream:** `sym2poly.cc` partfrac after `ext_factor`.

use giac_poly::{
    coeff_at, factor_into, partfrac_rational_terms, square_free_factorization, FlatUni, MainVar,
    Poly, FieldCoeff, PolyCoeff, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Zero;

use crate::error::EvalError;

use super::field_session::FieldSession;
use super::poly::{poly_algext_from_poly, PolyAlgExt};
use super::poly_alg_coeff::AlgExtCPolyCoeff;
use super::poly_alg_factor::factor_univariate_over_k;
use super::poly_alg_ops::{infer_ambient_field, normalize_algext_poly};

// **Pipeline private** — irreducible quadratic over ℚ with disc > 0 (splits in K).
fn quadratic_factor_needs_k_in_q(factor: &Poly, var: &Var) -> bool {
    if giac_poly::univariate_degree(factor, var) != 2 {
        return false;
    }
    if factor_into(factor).is_some_and(|fs| fs.len() > 1) {
        return false;
    }
    let a = coeff_at(factor, var, 2);
    let b = coeff_at(factor, var, 1);
    let c = coeff_at(factor, var, 0);
    if a.is_zero() {
        return false;
    }
    let disc = b.clone() * b - Ratio::from_integer(BigInt::from(4)) * a.clone() * c;
    disc > Ratio::zero()
}

/// **Stable** — ℚ partfrac cannot split an irreducible quadratic (disc>0) → use K[var].
///
/// Covers constant numerator (`1/(x²−2)`) and general proper rationals (`x/(x²−2)`).
pub fn partfrac_needs_k_split(num: &Poly, den: &Poly, var: &Var) -> bool {
    if den.is_zero() {
        return false;
    }
    if let Ok(sqff) = square_free_factorization(den, var) {
        if sqff
            .iter()
            .any(|(f, _)| quadratic_factor_needs_k_in_q(f, var))
        {
            return true;
        }
    }
    match partfrac_rational_terms(num, den, var) {
        Err(EvalError::NotImplemented(msg)) if msg == "partfrac real quadratic split" => true,
        Ok((_, terms)) => terms
            .iter()
            .any(|(_, d)| quadratic_factor_needs_k_in_q(d, var)),
        _ => false,
    }
}

/// **Stable** — sqff quadratic factor irreducible in ℚ but splits in K (disc > 0).
pub fn partfrac_sqff_factor_needs_k(factor: &Poly, var: &Var) -> bool {
    quadratic_factor_needs_k_in_q(factor, var)
}

// ponytail: was partfrac_rational_terms Err probe; disc>0 check is O(1) and covers non-constant numerators.

/// **Stable (bounded)** — `(poly_part, terms)` with denominators factored in K[var].
pub fn partfrac_rational_terms_over_k(
    num: &Poly,
    den: &Poly,
    var: &Var,
) -> Result<(Option<Poly>, Vec<(PolyAlgExt, PolyAlgExt)>), EvalError> {
    if den.is_zero() {
        return Err(EvalError::DivisionByZero);
    }
    let (poly_part, rem) = if giac_poly::univariate_degree(num, var)
        >= giac_poly::univariate_degree(den, var)
    {
        let (q, r) = num.div_rem(den);
        (Some(q), r)
    } else {
        (None, num.clone())
    };
    if rem.is_zero() {
        return Ok((poly_part, vec![]));
    }
    if giac_poly::univariate_degree(&rem, var) >= giac_poly::univariate_degree(den, var) {
        return Err(EvalError::TypeError("improper rational remainder"));
    }

    let num_a = poly_algext_from_poly(&rem)?;
    let den_a = poly_algext_from_poly(den)?;
    let field = infer_ambient_field(&den_a)?;
    let session = FieldSession::new(field);
    let num_n = normalize_algext_poly(&num_a, &session)?;
    let den_n = normalize_algext_poly(&den_a, &session)?;
    let flat = FlatUni::try_new(den_n, MainVar::new(var.clone())).map_err(Into::into)?;
    let pairs = factor_univariate_over_k(&session, &flat)?;

    if pairs.is_empty() {
        return Err(EvalError::TypeError("empty factorization"));
    }
    for (f, k) in &pairs {
        if *k != 1 || f.degree_wrt(var) != 1 {
            return Err(EvalError::NotImplemented("partfrac over K: nonlinear/repeated"));
        }
    }
    if pairs.len() <= 1 && giac_poly::univariate_degree(den, var) > 1 {
        return Err(EvalError::NotImplemented("partfrac over K: unsplit denominator"));
    }

    let factors: Vec<_> = pairs.into_iter().map(|(f, _)| f).collect();
    let mut terms = Vec::with_capacity(factors.len());
    for f in &factors {
        let root = linear_root_coeff(f, var)?;
        let mut prod = AlgExtCPolyCoeff::coeff_one();
        for g in &factors {
            if g != f {
                let val = eval_wrt(g, var, &root)?;
                prod = prod.coeff_mul(&val)?;
            }
        }
        if prod.coeff_is_zero() {
            return Err(EvalError::TypeError("repeated linear factor"));
        }
        let numer_val = eval_wrt(&num_n, var, &root)?;
        let coeff = numer_val.field_div(&prod)?;
        if !coeff.coeff_is_zero() {
            terms.push((PolyAlgExt::ring_constant(coeff), f.clone()));
        }
    }
    Ok((poly_part, terms))
}

// **Pipeline private** — `-c0/lc` for monic-linear `lc*x + c0`.
fn linear_root_coeff(f: &PolyAlgExt, var: &Var) -> Result<AlgExtCPolyCoeff, EvalError> {
    let lc = giac_poly::scalar_coeff_wrt(f, var, 1);
    if lc.coeff_is_zero() {
        return Err(EvalError::TypeError("not linear"));
    }
    let c0 = giac_poly::scalar_coeff_wrt(f, var, 0);
    let root = c0.field_div(&lc)?;
    root.coeff_neg().map_err(Into::into)
}

// **Pipeline private** — Horner evaluation in K.
fn eval_wrt(p: &PolyAlgExt, var: &Var, x: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
    let deg = p.degree_wrt(var);
    let mut acc = AlgExtCPolyCoeff::coeff_zero();
    for exp in (0..=deg).rev() {
        acc = acc.coeff_mul(x)?;
        acc = acc.coeff_add(&giac_poly::scalar_coeff_wrt(p, var, exp))?;
    }
    Ok(acc)
}

#[cfg(test)]
mod tests {
    

    use super::*;

    use num_bigint::BigInt;
    use num_rational::Ratio;

    fn x() -> Var {
        Var::from("x")
    }

    #[test]
    fn partfrac_one_over_x_squared_minus_two() {
        // Algorithm slice (no Expr): 1/(x²−2) → two sqff linear factors in K[x].
        let num = Poly::one();
        let den = Poly::var("x")
            .pow(2)
            .sub(&Poly::constant(Ratio::from_integer(BigInt::from(2))));
        assert!(partfrac_needs_k_split(&num, &den, &x()));
        let (_, terms) = partfrac_rational_terms_over_k(&num, &den, &x()).unwrap();
        assert_eq!(terms.len(), 2);
        assert!(terms.iter().all(|(_, d)| d.degree_wrt(&x()) == 1));
    }

    #[test]
    fn partfrac_x_over_x_squared_minus_two() {
        let num = Poly::var("x");
        let den = Poly::var("x")
            .pow(2)
            .sub(&Poly::constant(Ratio::from_integer(BigInt::from(2))));
        assert!(partfrac_needs_k_split(&num, &den, &x()));
        let (poly_part, terms) = partfrac_rational_terms_over_k(&num, &den, &x()).unwrap();
        assert!(poly_part.is_none());
        assert_eq!(terms.len(), 2);
        assert!(terms.iter().all(|(_, d)| d.degree_wrt(&x()) == 1));
        for (n, d) in &terms {
            super::super::poly::algext_poly_to_expr(n).expect("numer to expr");
            super::super::poly::algext_poly_to_expr(d).expect("denom to expr");
        }
    }

    #[test]
    fn partfrac_needs_k_split_false_for_x_squared_minus_one() {
        let num = Poly::one();
        let den = Poly::var("x").pow(2).sub(&Poly::one());
        assert!(!partfrac_needs_k_split(&num, &den, &x()));
    }
}
