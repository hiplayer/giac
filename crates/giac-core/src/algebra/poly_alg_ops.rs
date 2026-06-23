//! `Poly<AlgExtC>` field alignment and univariate ops over K (T0-1 / T0-2 bridge).
//!
//! **Upstream:** `common_EXT` + `ext_reduce` before `_EXT` polynomial arithmetic.

use std::sync::Arc;

use giac_poly::{
    scalar_coeff_wrt, univariate_div_rem_wrt, FlatUni, Poly, Var, monic_wrt, quo_exact_wrt,
};

use crate::error::EvalError;

use super::ext_tower::ExtensionField;
use super::field_session::FieldSession;
use super::poly::PolyAlgExt;
use super::poly_alg_coeff::AlgExtCPolyCoeff;

/// **Stable** — infer ambient **K** from polynomial coefficients (max-dimension field).
pub fn infer_ambient_field(p: &PolyAlgExt) -> Result<Arc<ExtensionField>, EvalError> {
    let mut field: Option<Arc<ExtensionField>> = None;
    for c in p.terms.values() {
        let f = Arc::clone(&c.as_inner().field);
        field = Some(match field {
            None => f,
            Some(prev) if f.dimension() > prev.dimension() => f,
            Some(prev) if prev.dimension() > f.dimension() => prev,
            Some(prev) if prev.dimension() == 1 && f.dimension() > 1 => f,
            Some(prev) => prev,
        });
    }
    Ok(field.unwrap_or_else(ExtensionField::rational))
}

/// **Stable** — lift all coefficients into session working field **L**.
pub fn normalize_algext_poly(
    p: &PolyAlgExt,
    session: &FieldSession,
) -> Result<PolyAlgExt, EvalError> {
    let mut out = PolyAlgExt::ring_zero();
    for (m, c) in &p.terms {
        let c = session.lift(c)?;
        out = out.try_add(
            &PolyAlgExt::ring_constant(c)
                .try_mul(&monomial_to_poly(m)?)?,
        )?;
    }
    Ok(out)
}

/// **Stable** — align two polynomials to a common working field **L** (monotone bump).
pub fn align_algext_polys(
    session: &FieldSession,
    a: &PolyAlgExt,
    b: &PolyAlgExt,
) -> Result<(PolyAlgExt, PolyAlgExt), EvalError> {
    ensure_common_field_for_polys(session, &[a, b])?;
    Ok((
        normalize_algext_poly(a, session)?,
        normalize_algext_poly(b, session)?,
    ))
}

/// **Stable** — bump session **L** to contain all coefficient fields in `polys`.
pub fn ensure_common_field_for_polys(
    session: &FieldSession,
    polys: &[&PolyAlgExt],
) -> Result<(), EvalError> {
    let mut anchor = session.zero();
    for p in polys {
        for c in p.terms.values() {
            let (_, c_aligned) = session.align(&anchor, c)?;
            anchor = c_aligned;
        }
    }
    Ok(())
}

/// **Stable** — `(q, r)` with `a = q*b + r` in K[var] after coefficient alignment.
pub fn div_rem_wrt_algext(
    session: &FieldSession,
    a: &PolyAlgExt,
    b: &PolyAlgExt,
    var: &Var,
) -> Result<(PolyAlgExt, PolyAlgExt), EvalError> {
    let (a, b) = align_algext_polys(session, a, b)?;
    univariate_div_rem_wrt(&a, &b, var).map_err(Into::into)
}

/// **Stable** — exact quotient in K[var]; `Err` if remainder nonzero.
pub fn quo_exact_wrt_algext(
    session: &FieldSession,
    num: &PolyAlgExt,
    den: &PolyAlgExt,
    var: &Var,
) -> Result<PolyAlgExt, EvalError> {
    let (a, b) = align_algext_polys(session, num, den)?;
    quo_exact_wrt(&a, &b, var).map_err(Into::into)
}

/// **Stable** — monic normalize w.r.t. `var` in K[var].
pub fn monic_wrt_algext(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<PolyAlgExt, EvalError> {
    let p = normalize_algext_poly(p, session)?;
    monic_wrt(&p, var).map_err(Into::into)
}

// **Pipeline private** — `monomial_to_poly`
fn monomial_to_poly(m: &giac_poly::Monomial) -> Result<PolyAlgExt, EvalError> {
    let mut out = PolyAlgExt::ring_one();
    for (v, e) in m.iter() {
        out = out.try_mul(&PolyAlgExt::ring_var(v.clone()).try_pow(e)?)?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use giac_poly::MainVar;

    use giac_poly::PolyCoeff;

    use super::*;
    use crate::algebra::test_fixtures::{k1_adjoin_sqrt2, sqrt2_algext};
    use crate::algebra::poly_alg_coeff::AlgExtCPolyCoeff;

    fn x_var() -> Var {
        Var::from("x")
    }

    fn sqrt2_coeff() -> AlgExtCPolyCoeff {
        AlgExtCPolyCoeff::from(
            crate::algebra::alg_ext_c::AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap(),
        )
    }

    fn x_squared_minus_2(session: &FieldSession) -> PolyAlgExt {
        let x = PolyAlgExt::ring_var(x_var());
        let two = session.int(2).unwrap();
        x.try_mul(&x).unwrap().try_sub(&PolyAlgExt::ring_constant(two)).unwrap()
    }

    fn x_minus_sqrt2(session: &FieldSession) -> PolyAlgExt {
        let alpha = session.lift(&sqrt2_coeff()).unwrap();
        PolyAlgExt::ring_var(x_var())
            .try_sub(&PolyAlgExt::ring_constant(alpha))
            .unwrap()
    }

    #[test]
    fn rem_x_squared_minus_2_mod_x_minus_sqrt2() {
        let k = k1_adjoin_sqrt2();
        let session = FieldSession::new(k);
        let p = x_squared_minus_2(&session);
        let d = x_minus_sqrt2(&session);
        let (_, r) = div_rem_wrt_algext(&session, &p, &d, &x_var()).unwrap();
        assert!(r.is_zero());
    }

    #[test]
    fn monic_wrt_leading_one() {
        let k = k1_adjoin_sqrt2();
        let session = FieldSession::new(k);
        let p = x_squared_minus_2(&session);
        let m = monic_wrt_algext(&session, &p, &x_var()).unwrap();
        let lc = scalar_coeff_wrt(&m, &x_var(), 2);
        assert!(lc.coeff_is_one());
    }

    #[test]
    fn align_scales_do_not_change_divisibility() {
        let k = k1_adjoin_sqrt2();
        let session = FieldSession::new(k);
        let p = x_squared_minus_2(&session);
        let d = x_minus_sqrt2(&session);
        let alpha = session.lift(&sqrt2_coeff()).unwrap();
        let scaled_p = p.try_mul(&PolyAlgExt::ring_constant(alpha)).unwrap();
        let (_, r1) = div_rem_wrt_algext(&session, &p, &d, &x_var()).unwrap();
        let (_, r2) = div_rem_wrt_algext(&session, &scaled_p, &d, &x_var()).unwrap();
        assert!(r1.is_zero());
        assert!(r2.is_zero());
    }

    #[test]
    fn flat_uni_quo_exact_matches_div_rem() {
        let k = k1_adjoin_sqrt2();
        let session = FieldSession::new(k);
        let p = x_squared_minus_2(&session);
        let d = x_minus_sqrt2(&session);
        let (a, b) = align_algext_polys(&session, &p, &d).unwrap();
        let flat_p = FlatUni::new(a, MainVar::new(x_var()));
        let flat_d = FlatUni::new(b, MainVar::new(x_var()));
        let (q, r) = flat_p.div_rem(&flat_d).unwrap();
        assert!(r.is_zero());
        assert_eq!(q, flat_p.exact_quo(&flat_d).unwrap());
    }
}
