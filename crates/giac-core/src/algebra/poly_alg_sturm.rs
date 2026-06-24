//! Sturm sequences over K = `AlgExtCPolyCoeff` (P3-4 / adoption B-04).
//!
//! **Upstream:** giac `sturm` / `sturmab` on `_EXT` coefficients.

use std::sync::Arc;

use giac_poly::{derivative_wrt, scalar_coeff_wrt, Poly, PolyCoeff, Var};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Zero;

use crate::error::EvalError;

use super::alg_ext::AlgExtData;
use super::alg_ext_c::AlgExtCData;
use super::field_arith::coords_to_expr;
use super::field_session::{is_negative_rational, FieldSession};
use super::poly::PolyAlgExt;
use super::poly_alg_coeff::AlgExtCPolyCoeff;
use super::poly::poly_algext_from_poly;
use super::poly_alg_ops::{
    div_rem_wrt_algext, ensure_common_field_for_polys, infer_ambient_field,
    normalize_algext_poly, square_free_part_wrt_algext,
};

/// **Stable** — classical Sturm chain in K[var]: P₀ = sqff(p), P₁ = P₀′, Pᵢ₊₁ = −rem(Pᵢ₋₁, Pᵢ).
pub fn sturm_sequence_wrt_algext(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<PolyAlgExt>, EvalError> {
    if p.degree_wrt(var) == 0 {
        return Err(EvalError::TypeError("constant polynomial"));
    }
    ensure_common_field_for_polys(session, &[p])?;
    let p_n = normalize_algext_poly(p, session)?;
    let sq = square_free_part_wrt_algext(session, &p_n, var)?;
    let dp = derivative_wrt(&sq, var).map_err(Into::<EvalError>::into)?;
    let mut seq = vec![sq, dp];
    while seq.len() >= 2 {
        let n = seq.len();
        let (_, r) = div_rem_wrt_algext(session, &seq[n - 2], &seq[n - 1], var)?;
        if r.is_zero() {
            break;
        }
        seq.push(negate_algext_poly(session, &r)?);
    }
    Ok(seq)
}

/// **Stable** — sign-variation count V(a) for a Sturm sequence at rational `a` ∈ ℚ ⊂ K.
pub fn sturm_sign_variations_at_algext(
    session: &FieldSession,
    seq: &[PolyAlgExt],
    var: &Var,
    a: &Ratio<BigInt>,
) -> Result<usize, EvalError> {
    let mut signs = Vec::new();
    for p in seq {
        let v = eval_algext_at_rational(session, p, var, a)?;
        let s = sign_of_algext_value(&v);
        if s != 0 {
            signs.push(s);
        }
    }
    Ok(sign_variations_i8(&signs))
}

/// **Stable** — root count in `(a, b]` for `p` ∈ K[var] (odd-multiplicity sqff convention on ℚ embed).
pub fn sturmab_count_wrt_algext(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
    a: &Ratio<BigInt>,
    b: &Ratio<BigInt>,
) -> Result<usize, EvalError> {
    let seq = sturm_sequence_wrt_algext(session, p, var)?;
    let va = sturm_sign_variations_at_algext(session, &seq, var, a)?;
    let b_eval = if b.is_zero() {
        Ratio::new(BigInt::from(-1), BigInt::from(1_000_000))
    } else {
        b.clone()
    };
    let vb = sturm_sign_variations_at_algext(session, &seq, var, &b_eval)?;
    Ok(va.saturating_sub(vb))
}

/// **Stable** — `sturmab_count_wrt_algext` for ℚ[var] via `poly_algext_from_poly` lift.
pub fn sturmab_count_rational_poly(
    p: &Poly,
    var: &Var,
    a: &Ratio<BigInt>,
    b: &Ratio<BigInt>,
) -> Result<usize, EvalError> {
    let p = poly_algext_from_poly(p)?;
    let field = infer_ambient_field(&p)?;
    let session = FieldSession::new(field);
    sturmab_count_wrt_algext(&session, &p, var, a, b)
}

// **Pipeline private** — Horner evaluation at x ∈ ℚ ⊂ L.
fn eval_algext_at_rational(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
    x: &Ratio<BigInt>,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let x_c = rational_coeff(session, x)?;
    let deg = p.degree_wrt(var);
    let mut acc = session.zero();
    for exp in (0..=deg).rev() {
        acc = acc.coeff_mul(&x_c)?;
        let c = session.lift(&scalar_coeff_wrt(p, var, exp))?;
        acc = acc.coeff_add(&c)?;
    }
    Ok(acc)
}

// **Pipeline private** — embed ℚ constant into current L.
fn rational_coeff(session: &FieldSession, r: &Ratio<BigInt>) -> Result<AlgExtCPolyCoeff, EvalError> {
    let field = session.working();
    let coords = field.embed_rational(r);
    let a = AlgExtData::from_field_coords(Arc::clone(&field), coords_to_expr(&coords)?)?;
    Ok(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&a)?))
}

// **Pipeline private** — −p in K[var].
fn negate_algext_poly(session: &FieldSession, p: &PolyAlgExt) -> Result<PolyAlgExt, EvalError> {
    let neg_one = session.int(-1)?;
    let mut out = PolyAlgExt::ring_zero();
    for (m, c) in &p.terms {
        let term = PolyAlgExt::ring_constant(session.lift(c)?.coeff_mul(&neg_one)?);
        let term = if m.is_const() {
            term
        } else {
            let mut t = term;
            for (v, e) in m.iter() {
                t = t.try_mul(&PolyAlgExt::ring_var(v.clone()).try_pow(e)?)?;
            }
            t
        };
        out = out.try_add(&term)?;
    }
    Ok(out)
}

// **Pipeline private** — sign for Sturm (rational or positive real constant in ℚ ⊂ K).
fn sign_of_algext_value(c: &AlgExtCPolyCoeff) -> i8 {
    if c.coeff_is_zero() {
        0
    } else if is_negative_rational(c) {
        -1
    } else {
        1
    }
}

// **Pipeline private** — sign changes in Sturm sequence values.
fn sign_variations_i8(signs: &[i8]) -> usize {
    if signs.len() < 2 {
        return 0;
    }
    signs.windows(2).filter(|w| w[0] != w[1]).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_poly::Poly;

    use crate::algebra::poly::poly_algext_from_poly;

    fn x_var() -> Var {
        Var::from("x")
    }

    fn session_for(p: &PolyAlgExt) -> FieldSession {
        let field = super::super::poly_alg_ops::infer_ambient_field(p).unwrap();
        FieldSession::new(field)
    }

    // **B** — Sturm chain length for x²−2 over K.
    #[test]
    fn sturm_sequence_x_squared_minus_two() {
        let p_q = Poly::var("x")
            .pow(2)
            .sub(&Poly::constant(Ratio::from_integer(2.into())));
        let p = poly_algext_from_poly(&p_q).unwrap();
        let session = session_for(&p);
        let seq = sturm_sequence_wrt_algext(&session, &p, &x_var()).unwrap();
        assert_eq!(seq.len(), 3);
    }

    // **B** — two real roots of x²−2 in (−10, 10).
    #[test]
    fn sturmab_count_x_squared_minus_two() {
        let p_q = Poly::var("x")
            .pow(2)
            .sub(&Poly::constant(Ratio::from_integer(2.into())));
        let p = poly_algext_from_poly(&p_q).unwrap();
        let session = session_for(&p);
        let count = sturmab_count_wrt_algext(
            &session,
            &p,
            &x_var(),
            &Ratio::from_integer((-10).into()),
            &Ratio::from_integer(10.into()),
        )
        .unwrap();
        assert_eq!(count, 2);
    }
}
