//! Univariate resultant over K = `AlgExtCPolyCoeff` (P3-4 / adoption B-04).
//!
//! **Upstream:** giac `resultant` on `_EXT` coefficients.

use giac_poly::{scalar_coeff_wrt, FieldCoeff, PolyCoeff, Var};

use crate::error::EvalError;

use super::field_session::FieldSession;
use super::poly::PolyAlgExt;
use super::poly_alg_coeff::AlgExtCPolyCoeff;
use super::poly_alg_ops::{
    ensure_common_field_for_polys, gcd_wrt_algext, normalize_algext_poly,
};

/// **Stable** — Sylvester resultant of univariate `a`, `b` ∈ K[var]; scalar in K.
pub fn resultant_wrt_algext(
    session: &FieldSession,
    a: &PolyAlgExt,
    b: &PolyAlgExt,
    var: &Var,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    ensure_common_field_for_polys(session, &[a, b])?;
    let a = normalize_algext_poly(a, session)?;
    let b = normalize_algext_poly(b, session)?;
    let g = gcd_wrt_algext(session, &a, &b, var)?;
    if g.degree_wrt(var) > 0 {
        return Ok(session.zero());
    }
    let da = a.degree_wrt(var);
    let db = b.degree_wrt(var);
    if da == 0 || db == 0 {
        let lc = if da == 0 {
            leading_coeff_wrt(session, &a, var)?
        } else {
            leading_coeff_wrt(session, &b, var)?
        };
        let exp = if da == 0 { db } else { da };
        return coeff_pow(&lc, exp);
    }
    if da == 1 && db == 1 {
        return sylvester_det2_wrt(session, &a, &b, var);
    }
    let ac = univariate_coeffs_asc(session, &a, var, da)?;
    let bc = univariate_coeffs_asc(session, &b, var, db)?;
    sylvester_det_field(&ac, da as usize, &bc, db as usize)
}

// **Pipeline private** — leading coefficient in K[var].
fn leading_coeff_wrt(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let deg = p.degree_wrt(var);
    for exp in (0..=deg).rev() {
        let c = session.lift(&scalar_coeff_wrt(p, var, exp))?;
        if !c.coeff_is_zero() {
            return Ok(c);
        }
    }
    Ok(session.zero())
}

// **Pipeline private** — ascending coeffs c₀ + c₁ var + … in K.
fn univariate_coeffs_asc(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
    degree: u64,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    (0..=degree)
        .map(|e| session.lift(&scalar_coeff_wrt(p, var, e)))
        .collect()
}

// **Pipeline private** — Res(ax+b, cx+d) = ad − bc.
fn sylvester_det2_wrt(
    session: &FieldSession,
    a: &PolyAlgExt,
    b: &PolyAlgExt,
    var: &Var,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let a1 = session.lift(&scalar_coeff_wrt(a, var, 1))?;
    let a0 = session.lift(&scalar_coeff_wrt(a, var, 0))?;
    let b1 = session.lift(&scalar_coeff_wrt(b, var, 1))?;
    let b0 = session.lift(&scalar_coeff_wrt(b, var, 0))?;
    a1.coeff_mul(&b0)?.coeff_sub(&a0.coeff_mul(&b1)?)
}

// **Pipeline private** — Sylvester matrix determinant in K.
fn sylvester_det_field(
    a: &[AlgExtCPolyCoeff],
    deg_a: usize,
    b: &[AlgExtCPolyCoeff],
    deg_b: usize,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let n = deg_b;
    let m = deg_a;
    let size = m + n;
    let mut mat = vec![vec![AlgExtCPolyCoeff::coeff_zero(); size]; size];
    for i in 0..n {
        for j in 0..=m {
            mat[i][i + j] = a[m - j].clone();
        }
    }
    for j in 0..m {
        for k in 0..=n {
            mat[n + j][j + k] = b[n - k].clone();
        }
    }
    det_field(&mut mat)
}

// **Pipeline private** — Gaussian elimination determinant in K.
fn det_field(mat: &mut [Vec<AlgExtCPolyCoeff>]) -> Result<AlgExtCPolyCoeff, EvalError> {
    let n = mat.len();
    let mut det = AlgExtCPolyCoeff::coeff_one();
    for k in 0..n {
        let mut pivot_row = k;
        while pivot_row < n && mat[pivot_row][k].coeff_is_zero() {
            pivot_row += 1;
        }
        if pivot_row == n {
            return Ok(AlgExtCPolyCoeff::coeff_zero());
        }
        if pivot_row != k {
            mat.swap(k, pivot_row);
            det = det.coeff_neg()?;
        }
        let pivot = mat[k][k].clone();
        det = det.coeff_mul(&pivot)?;
        for i in (k + 1)..n {
            if mat[i][k].coeff_is_zero() {
                continue;
            }
            let factor = mat[i][k].field_div(&pivot)?;
            for j in (k + 1)..n {
                let sub = mat[k][j].coeff_mul(&factor)?;
                mat[i][j] = mat[i][j].coeff_sub(&sub)?;
            }
        }
    }
    Ok(det)
}

// **Pipeline private** — c^exp in K.
fn coeff_pow(c: &AlgExtCPolyCoeff, exp: u64) -> Result<AlgExtCPolyCoeff, EvalError> {
    let mut out = AlgExtCPolyCoeff::coeff_one();
    let mut base = c.clone();
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            out = out.coeff_mul(&base)?;
        }
        if e > 1 {
            base = base.coeff_mul(&base)?;
        }
        e >>= 1;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use giac_poly::Poly;
    use num_rational::Ratio;

    use super::*;
    use crate::algebra::poly::poly_algext_from_poly;
    use crate::algebra::test_fixtures::{k1_adjoin_sqrt2, sqrt2_algext};
    use crate::algebra::poly_alg_coeff::AlgExtCPolyCoeff;

    fn x_var() -> Var {
        Var::from("x")
    }

    fn session_for(p: &PolyAlgExt) -> FieldSession {
        let field = super::super::poly_alg_ops::infer_ambient_field(p).unwrap();
        FieldSession::new(field)
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

    fn from_q(p: &Poly) -> PolyAlgExt {
        poly_algext_from_poly(p).unwrap()
    }

    // **B** — Res(x−1, x+1) = 2 over ℚ ⊂ K.
    #[test]
    fn resultant_two_linear_polys() {
        let a = Poly::var("x").sub(&Poly::one());
        let b = Poly::var("x").add(&Poly::one());
        let pa = from_q(&a);
        let pb = from_q(&b);
        let session = session_for(&pa);
        let r = resultant_wrt_algext(&session, &pa, &pb, &x_var()).unwrap();
        assert_eq!(r, session.int(2).unwrap());
    }

    // **B** — common factor → zero.
    #[test]
    fn resultant_shared_factor_is_zero() {
        let a = Poly::var("x").pow(2).sub(&Poly::constant(Ratio::from_integer(2.into())));
        let pa = from_q(&a);
        let session = session_for(&pa);
        let r = resultant_wrt_algext(&session, &pa, &pa, &x_var()).unwrap();
        assert!(r.coeff_is_zero());
    }

    // **B** — Res(x²−2, x−√2) = 0 over K.
    #[test]
    fn resultant_x_squared_minus_2_and_x_minus_sqrt2() {
        let k = k1_adjoin_sqrt2();
        let session = FieldSession::new(k);
        let p = x_squared_minus_2(&session);
        let d = x_minus_sqrt2(&session);
        let r = resultant_wrt_algext(&session, &p, &d, &x_var()).unwrap();
        assert!(r.coeff_is_zero());
    }

    // **B** — coprime quadratics: Res(x²−2, x+1) = −1.
    #[test]
    fn resultant_x_squared_minus_2_and_x_plus_one() {
        let a = Poly::var("x").pow(2).sub(&Poly::constant(Ratio::from_integer(2.into())));
        let b = Poly::var("x").add(&Poly::one());
        let pa = from_q(&a);
        let pb = from_q(&b);
        let session = session_for(&pa);
        let r = resultant_wrt_algext(&session, &pa, &pb, &x_var()).unwrap();
        assert_eq!(r, session.int(-1).unwrap());
    }

    // **B** — quadratic × quadratic over ℚ.
    #[test]
    fn resultant_quadratic() {
        let a = Poly::var("x").pow(2).add(&Poly::one());
        let b = Poly::var("x").pow(2).sub(&Poly::one());
        let pa = from_q(&a);
        let pb = from_q(&b);
        let session = session_for(&pa);
        let r = resultant_wrt_algext(&session, &pa, &pb, &x_var()).unwrap();
        assert_eq!(r, session.int(4).unwrap());
    }
}
