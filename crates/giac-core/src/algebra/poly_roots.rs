//! Univariate exact roots of `Poly<AlgExtCPolyCoeff>` (P2-1/6, P3-6).
//!
//! ```text
//! infer K from coefficients → normalize → monic → FieldSession::new(K) → deg dispatch
//!   1: linear
//!   2: quadratic + sqrt(Δ)
//!   3: one Cardano/cbrt root → deflate → quadratic
//!   4: biquadratic or resolvent cubic → quadratics
//! ```

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::sync::Arc;

use giac_poly::{PolyCoeff, Var};

use crate::error::EvalError;

use super::alg_ext::{algext_cube_root, algext_square_roots, AlgExtData};
use super::alg_ext_c::AlgExtCData;
use super::ext_tower::ExtensionField;
use super::field_arith::coords_to_expr;
use super::field_session::FieldSession;
use super::poly::PolyAlgExt;
use super::poly_alg_coeff::AlgExtCPolyCoeff;

/// Exact roots of univariate `p` w.r.t. `var` over the coefficient field of `p`.
/// **Stable (bounded)** — exact AlgExtC roots deg 1–4; quartic resolvent gap
pub fn poly_algext_roots(p: &PolyAlgExt, var: &Var) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let ambient = infer_field(p)?;
    let mut session = FieldSession::new(Arc::clone(&ambient));
    let p = normalize_coeffs(p, &mut session)?;
    let monic = monic_univariate(&p, var, &session)?;
    let d = monic.degree_wrt(var);
    match d {
        0 => {
            if monic.is_zero() {
                Ok(vec![])
            } else {
                Err(EvalError::TypeError("constant has no roots"))
            }
        }
        1 => linear_root(&mut session, &monic, var),
        2 => quadratic_roots(&mut session, &monic, var),
        3 => cubic_roots(&mut session, &monic, var),
        4 => quartic_roots(&mut session, &monic, var),
        _ => Err(EvalError::NotImplemented("PolyAlgExt::roots")),
    }
}

// **Pipeline private** — infer ambient K from PolyAlgExt coefficients
fn infer_field(p: &PolyAlgExt) -> Result<Arc<ExtensionField>, EvalError> {
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

// **Stable** — univariate coefficient at exponent
fn coeff_at(p: &PolyAlgExt, var: &Var, exp: u64, session: &FieldSession) -> AlgExtCPolyCoeff {
    for (m, c) in &p.terms {
        if exp == 0 && m.is_const() {
            return c.clone();
        }
        if m.exp_of(var) == exp && m.iter().all(|(v, _)| v == var) {
            return c.clone();
        }
    }
    session.zero()
}

// **Pipeline private** — lift all coeffs to ambient K
fn normalize_coeffs(
    p: &PolyAlgExt,
    session: &mut FieldSession,
) -> Result<PolyAlgExt, EvalError> {
    let mut out = PolyAlgExt::ring_zero();
    for (m, c) in &p.terms {
        let c = session.lift(c)?;
        out = out.try_add(&PolyAlgExt::ring_constant(c).try_mul(&monomial_to_poly(m)?)?)?;
    }
    Ok(out)
}

// **Pipeline private** — `monomial_to_poly`
fn monomial_to_poly(m: &giac_poly::Monomial) -> Result<PolyAlgExt, EvalError> {
    let mut out = PolyAlgExt::ring_one();
    for (v, e) in m.iter() {
        out = out.try_mul(&PolyAlgExt::ring_var(v.clone()).try_pow(e)?)?;
    }
    Ok(out)
}

// **Pipeline private** — `monic_univariate`
fn monic_univariate(
    p: &PolyAlgExt,
    var: &Var,
    session: &FieldSession,
) -> Result<PolyAlgExt, EvalError> {
    let deg = p.degree_wrt(var);
    let lc = coeff_at(p, var, deg, session);
    if lc.coeff_is_zero() {
        return Err(EvalError::TypeError("leading coefficient zero"));
    }
    for (m, _) in &p.terms {
        if !m.iter().all(|(v, _)| v == var) {
            return Err(EvalError::TypeError("not univariate"));
        }
    }
    let inv = lc.coeff_inv()?;
    let mut out = PolyAlgExt::ring_zero();
    for exp in 0..=deg {
        let c = coeff_at(p, var, exp, session);
        if c.coeff_is_zero() {
            continue;
        }
        let scaled = c.coeff_mul(&inv)?;
        let mut term = PolyAlgExt::ring_constant(scaled);
        if exp > 0 {
            term = term.try_mul(&PolyAlgExt::ring_var(var.clone()).try_pow(exp)?)?;
        }
        out = out.try_add(&term)?;
    }
    Ok(out)
}

// **Pipeline private** — `linear_root`
fn linear_root(
    session: &mut FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let a = coeff_at(p, var, 1, session);
    let b = coeff_at(p, var, 0, session);
    if a.coeff_is_zero() {
        return Err(EvalError::TypeError("not linear"));
    }
    Ok(vec![session.div(&b, &a)?])
}

// **Pipeline private** — `quadratic_roots`
fn quadratic_roots(
    session: &mut FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let b = coeff_at(p, var, 1, session);
    let c = coeff_at(p, var, 0, session);
    let four = session.int(4)?;
    let bb = session.mul(&b, &b)?;
    let c4 = session.mul(&c, &four)?;
    let disc = bb.coeff_sub(&c4)?;
    if disc.coeff_is_zero() {
        let two = session.int(2)?;
        let nb = session.neg(&b)?;
        return Ok(vec![session.div(&nb, &two)?]);
    }
    let sqrt_d = sqrt_disc(&disc)?;
    session.bump_to(&sqrt_d.as_inner().field);
    let b_lift = session.lift(&b)?;
    let mut nb = session.neg(&b_lift)?;
    let mut sqrt_d = session.lift(&sqrt_d)?;
    let two0 = session.int(2)?;
    let mut two = session.lift(&two0)?;
    let (nb, sqrt_d) = session.align(&nb, &sqrt_d)?;
    two = session.lift(&two)?;
    let (two, sqrt_d) = session.align(&two, &sqrt_d)?;
    let num_plus = session.add(&nb, &sqrt_d)?;
    let r_plus = session.div(&num_plus, &two)?;
    let neg_sqrt = session.neg(&sqrt_d)?;
    let num_minus = session.add(&nb, &neg_sqrt)?;
    let r_minus = session.div(&num_minus, &two)?;
    Ok(vec![r_plus, r_minus])
}

// **Pipeline private** — `sqrt_disc`
fn sqrt_disc(disc: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
    if let Ok(mut roots) = algext_c_sqrt(disc) {
        if let Some(r) = roots.pop() {
            if r.as_inner().im.iter().all(|e| e.is_zero()) {
                return Ok(r);
            }
        }
    }
    let r = algext_c_sqrt(&disc.coeff_neg()?)?
        .into_iter()
        .next()
        .ok_or(EvalError::NotImplemented("algext sqrt"))?;
    mul_i(&r)
}

// **Pipeline private** — `mul_i`
fn mul_i(z: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
    let inner = z.as_inner();
    Ok(AlgExtCPolyCoeff::from(AlgExtCData {
        field: Arc::clone(&inner.field),
        re: coords_to_expr(&inner.field.zero_coords())?,
        im: inner.re.clone(),
        root_index: None,
    }))
}

// **Pipeline private** — `cubic_roots`
fn cubic_roots(
    session: &mut FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let a2 = coeff_at(p, var, 2, session);
    let a1 = coeff_at(p, var, 1, session);
    if a2.coeff_is_zero() && a1.coeff_is_zero() {
        let a0 = coeff_at(p, var, 0, session);
        let r = algext_c_cube_root(&a0.coeff_neg()?)?;
        session.bump_to(&r.as_inner().field);
        return Ok(vec![session.lift(&r)?]);
    }
    let r0 = one_cubic_root(session, p, var)?;
    let mut roots = vec![r0.clone()];
    let quad = deflate_monic(session, p, var, &r0)?;
    roots.extend(quadratic_roots(session, &quad, var)?);
    Ok(roots)
}

// **Pipeline private** — `one_cubic_root`
fn one_cubic_root(
    session: &mut FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let a2 = coeff_at(p, var, 2, session);
    let a1 = coeff_at(p, var, 1, session);
    let a0 = coeff_at(p, var, 0, session);
    let three = session.int(3)?;
    let neg_a2 = session.neg(&a2)?;
    let shift = session.div(&neg_a2, &three)?;
    let a2_shift = session.mul(&a2, &shift)?;
    let neg_a2_shift = session.neg(&a2_shift)?;
    let p_num = session.add(&a1, &neg_a2_shift)?;
    let p_dep = session.div(&p_num, &three)?;
    let t1 = session.mul(&a1, &shift)?;
    let a2s = session.mul(&a2, &shift)?;
    let t2 = session.mul(&a2s, &shift)?;
    let ss = session.mul(&shift, &shift)?;
    let t3 = session.mul(&ss, &shift)?;
    let neg_t1 = session.neg(&t1)?;
    let s0 = session.add(&a0, &neg_t1)?;
    let s1 = session.add(&s0, &t2)?;
    let q_dep = session.add(&s1, &t3)?;
    if p_dep.coeff_is_zero() {
        let neg_q = session.neg(&q_dep)?;
        let r = algext_c_cube_root(&neg_q)?;
        session.bump_to(&r.as_inner().field);
        let r_lift = session.lift(&r)?;
        let shift_lift = session.lift(&shift)?;
        return session.add(&r_lift, &shift_lift);
    }
    let half = session.half()?;
    let twenty_seven = session.int(27)?;
    let qh = session.mul(&q_dep, &half)?;
    let qh2 = session.mul(&qh, &qh)?;
    let p2 = session.mul(&p_dep, &p_dep)?;
    let p3 = session.mul(&p2, &p_dep)?;
    let term = session.div(&p3, &twenty_seven)?;
    let delta = session.add(&qh2, &term)?;
    let sqrt_delta = algext_c_sqrt(&delta)?
        .into_iter()
        .next()
        .ok_or(EvalError::NotImplemented("algext sqrt"))?;
    session.bump_to(&sqrt_delta.as_inner().field);
    let q_half = session.mul(&q_dep, &half)?;
    let neg_q_half = session.neg(&q_half)?;
    let sqrt_lift = session.lift(&sqrt_delta)?;
    let u_arg = session.add(&neg_q_half, &sqrt_lift)?;
    let u = algext_c_cube_root(&u_arg)?;
    session.bump_to(&u.as_inner().field);
    let neg_sqrt = session.neg(&sqrt_lift)?;
    let v_arg = session.add(&neg_q_half, &neg_sqrt)?;
    let v = algext_c_cube_root(&v_arg)?;
    let (u, v) = session.align(&u, &v)?;
    let uv = session.add(&u, &v)?;
    let shift_lift = session.lift(&shift)?;
    session.add(&uv, &shift_lift)
}

// **Pipeline private** — `quartic_roots`
fn quartic_roots(
    session: &mut FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let a3 = coeff_at(p, var, 3, session);
    let a1 = coeff_at(p, var, 1, session);
    if a3.coeff_is_zero() && a1.coeff_is_zero() {
        return biquadratic_roots(session, p, var);
    }
    let four = session.int(4)?;
    let neg_a3 = session.neg(&a3)?;
    let shift = session.div(&neg_a3, &four)?;
    let dep = depress_quartic(session, p, var, &shift)?;
    let p2 = coeff_at(&dep, var, 2, session);
    let p1 = coeff_at(&dep, var, 1, session);
    let p0 = coeff_at(&dep, var, 0, session);
    let res = build_resolvent_cubic(session, &p2, &p1, &p0)?;
    let z_roots = cubic_roots(session, &res, &Var::from("_z"))?;
    let mut all = Vec::new();
    for z in z_roots {
        let mut rs = split_depressed_quartic(session, &dep, var, &z)?;
        for r in rs.drain(..) {
            if !all.iter().any(|x: &AlgExtCPolyCoeff| x.eq_mod(&r).unwrap_or(false)) {
                all.push(r);
            }
        }
    }
    if all.len() != 4 {
        return Err(EvalError::NotImplemented("quartic roots"));
    }
    Ok(all)
}

// **Pipeline private** — `biquadratic_roots`
fn biquadratic_roots(
    session: &mut FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let b = coeff_at(p, var, 2, session);
    let c = coeff_at(p, var, 0, session);
    let u = PolyAlgExt::ring_var(Var::from("u"))
        .try_pow(2)?
        .try_add(&PolyAlgExt::ring_constant(b))?
        .try_add(&PolyAlgExt::ring_constant(c))?;
    let u_roots = quadratic_roots(session, &u, &Var::from("u"))?;
    let mut out = Vec::new();
    for ur in u_roots {
        out.extend(algext_c_sqrt(&ur)?);
    }
    let _ = var;
    Ok(out)
}

// **Pipeline private** — `depress_quartic`
fn depress_quartic(
    session: &mut FieldSession,
    p: &PolyAlgExt,
    var: &Var,
    shift: &AlgExtCPolyCoeff,
) -> Result<PolyAlgExt, EvalError> {
    let x = PolyAlgExt::ring_var(var.clone());
    let t = x.try_add(&PolyAlgExt::ring_constant(shift.clone()))?;
    let mut sum = PolyAlgExt::ring_zero();
    for e in 0..=4 {
        let c = coeff_at(p, var, e, session);
        if c.coeff_is_zero() {
            continue;
        }
        sum = sum.try_add(&PolyAlgExt::ring_constant(c).try_mul(&t.try_pow(e)?)?)?;
    }
    monic_univariate(&sum, var, session)
}

// **Pipeline private** — `build_resolvent_cubic`
fn build_resolvent_cubic(
    session: &mut FieldSession,
    p: &AlgExtCPolyCoeff,
    q: &AlgExtCPolyCoeff,
    r: &AlgExtCPolyCoeff,
) -> Result<PolyAlgExt, EvalError> {
    let four = session.int(4)?;
    let z = PolyAlgExt::ring_var(Var::from("_z"));
    let r4 = session.mul(r, &four)?;
    let neg_r4 = session.neg(&r4)?;
    let rq = session.mul(r, p)?;
    let rqp4 = session.mul(&rq, &four)?;
    let qq = session.mul(q, q)?;
    let neg_qq = session.neg(&qq)?;
    let const_term = session.add(&rqp4, &neg_qq)?;
    Ok(z.try_pow(3)?
        .try_sub(&PolyAlgExt::ring_constant(p.clone()).try_mul(&z.try_pow(2)?)?)?
        .try_sub(&PolyAlgExt::ring_constant(neg_r4).try_mul(&z)?)?
        .try_add(&PolyAlgExt::ring_constant(const_term))?)
}

// **Pipeline private** — `split_depressed_quartic`
fn split_depressed_quartic(
    session: &mut FieldSession,
    dep: &PolyAlgExt,
    var: &Var,
    z: &AlgExtCPolyCoeff,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let q = coeff_at(dep, var, 1, session);
    let r = coeff_at(dep, var, 0, session);
    let four = session.int(4)?;
    let zz = session.mul(z, z)?;
    let r4 = session.mul(&r, &four)?;
    let neg_r4 = session.neg(&r4)?;
    let disc = session.add(&zz, &neg_r4)?;
    let m = algext_c_sqrt(&disc)?
        .into_iter()
        .next()
        .ok_or(EvalError::NotImplemented("algext sqrt"))?;
    session.bump_to(&m.as_inner().field);
    let half = session.half()?;
    let m_lift = session.lift(&m)?;
    let neg_m = session.neg(&m_lift)?;
    let mut out = Vec::new();
    for mp in [session.add(z, &m_lift)?, session.add(z, &neg_m)?] {
        let c0 = session.mul(&mp, &half)?;
        let c1 = if q.coeff_is_zero() {
            session.zero()
        } else {
            let q_over_m = session.div(&q, &m_lift)?;
            session.mul(&q_over_m, &half)?
        };
        let quad = PolyAlgExt::ring_var(var.clone())
            .try_pow(2)?
            .try_add(&PolyAlgExt::ring_constant(c1).try_mul(&PolyAlgExt::ring_var(var.clone()))?)?
            .try_add(&PolyAlgExt::ring_constant(c0))?;
        out.extend(quadratic_roots(session, &quad, var)?);
    }
    Ok(out)
}

// **Pipeline private** — `deflate_monic`
fn deflate_monic(
    session: &mut FieldSession,
    p: &PolyAlgExt,
    var: &Var,
    root: &AlgExtCPolyCoeff,
) -> Result<PolyAlgExt, EvalError> {
    session.bump_to(&root.as_inner().field);
    let n = p.degree_wrt(var) as usize;
    let mut qs = vec![session.zero(); n];
    qs[n - 1] = session.one();
    for k in (0..n - 1).rev() {
        let ak = session.lift(&coeff_at(p, var, k as u64, session))?;
        let term = session.mul(root, &qs[k + 1])?;
        qs[k] = session.add(&ak, &term)?;
    }
    let mut out = PolyAlgExt::ring_zero();
    for (exp, c) in qs.iter().enumerate() {
        if c.coeff_is_zero() {
            continue;
        }
        let mut term = PolyAlgExt::ring_constant(c.clone());
        if exp > 0 {
            term = term.try_mul(&PolyAlgExt::ring_var(var.clone()).try_pow(exp as u64)?)?;
        }
        out = out.try_add(&term)?;
    }
    Ok(out)
}

// **Pipeline private** — `algext_c_sqrt`
fn algext_c_sqrt(z: &AlgExtCPolyCoeff) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let inner = z.as_inner();
    if inner.im.iter().all(|e| e.is_zero()) {
        let re = AlgExtData::from_field_coords(Arc::clone(&inner.field), inner.re.clone())?;
        let mut out = Vec::new();
        for a in algext_square_roots(&re)? {
            out.push(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&a)?));
        }
        return Ok(out);
    }
    Err(EvalError::NotImplemented("algext c sqrt complex"))
}

// **Pipeline private** — `algext_c_cube_root`
fn algext_c_cube_root(z: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
    let inner = z.as_inner();
    if inner.im.iter().all(|e| e.is_zero()) {
        let re = AlgExtData::from_field_coords(Arc::clone(&inner.field), inner.re.clone())?;
        let a = algext_cube_root(&re)?;
        return Ok(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&a)?));
    }
    Err(EvalError::NotImplemented("algext c cbrt complex"))
}

trait CoeffInv {
    // **Pipeline private** — `coeff_inv`
    fn coeff_inv(&self) -> Result<Self, EvalError>
    where
        Self: Sized;
    // **Pipeline private** — `eq_mod`
    fn eq_mod(&self, other: &Self) -> Result<bool, EvalError>;
}

impl CoeffInv for AlgExtCPolyCoeff {
    // **Stable** — `Poly::coeff_inv`
    fn coeff_inv(&self) -> Result<Self, EvalError> {
        Ok(Self::from(self.as_inner().inv()?))
    }
    // **Stable** — `Poly::eq_mod`
    fn eq_mod(&self, other: &Self) -> Result<bool, EvalError> {
        self.as_inner().eq_mod(other.as_inner())
    }
}

// **Pipeline private** — verify root vanishes mod minpoly
fn verify_root(p: &PolyAlgExt, var: &Var, root: &AlgExtCPolyCoeff) {
    let ambient = infer_field(p).expect("infer field");
    let mut session = FieldSession::new(ambient);
    let p = monic_univariate(
        &normalize_coeffs(p, &mut session).expect("normalize"),
        var,
        &session,
    )
    .expect("monic");
    session.bump_to(&root.as_inner().field);
    let mut val = session.zero();
    for (m, c) in &p.terms {
        let exp = m.exp_of(var);
        let mut pow = session.one();
        for _ in 0..exp {
            pow = session.mul(root, &pow).unwrap();
        }
        let c = session.lift(c).unwrap();
        let (c, pow) = session.align(&c, &pow).unwrap();
        let term = session.mul(&c, &pow).unwrap();
        val = session.add(&val, &term).unwrap();
    }
    assert!(val.coeff_is_zero(), "root does not vanish");
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_poly::{PolyCoeff, Var};
    use num_rational::Ratio;
    use serial_test::serial;

    use crate::algebra::test_fixtures::sqrt2_algext;
    use crate::expr::Expr;

    use super::*;
    use super::super::alg_ext_c::AlgExtCData;
    use super::super::poly::poly_alg_from_expr;

    fn q_session() -> FieldSession {
        FieldSession::new(ExtensionField::rational())
    }

    fn rat_coeff(session: &FieldSession, n: i64) -> AlgExtCPolyCoeff {
        session.int(n).unwrap()
    }

    #[serial]
    #[test]
    fn cubic_one_root_vanishes() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(3)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(&session, 2)))
            .unwrap();
        let r = one_cubic_root(&mut session, &p, &Var::from("t")).unwrap();
        verify_root(&p, &Var::from("t"), &r);
    }

    #[serial]
    #[test]
    fn quadratic_sqrt_four_times_sqrt2_over_k1() {
        let k1 = Arc::clone(
            &AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap())
                .as_inner()
                .field,
        );
        let mut session = FieldSession::new(k1);
        let sqrt2 = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap());
        let four = session.int(4).unwrap();
        let disc = session.mul(&sqrt2, &four).unwrap();
        let rs = algext_c_sqrt(&disc).unwrap();
        assert_eq!(rs.len(), 2);
        for beta in &rs {
            session.bump_to(&beta.as_inner().field);
            let sq = session.mul(beta, beta).unwrap();
            let (sq_a, d_a) = session.align(&sq, &disc).unwrap();
            assert!(sq_a.coeff_sub(&d_a).unwrap().coeff_is_zero());
            let two = session.int(2).unwrap();
            let half = session.half().unwrap();
            let root_mul = session.mul(beta, &half).unwrap();
            let root_div = session.div(beta, &two).unwrap();
            let (rm, rd) = session.align(&root_mul, &root_div).unwrap();
            assert!(
                rm.coeff_sub(&rd).unwrap().coeff_is_zero(),
                "mul/2 and div/2 should agree"
            );
            let sqrt2 = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap());
            let root_sq = session.mul(&root_div, &root_div).unwrap();
            let (rsq, s2) = session.align(&root_sq, &sqrt2).unwrap();
            assert!(
                rsq.coeff_sub(&s2).unwrap().coeff_is_zero(),
                "root^2 should equal sqrt2"
            );
        }
    }

    #[serial]
    #[test]
    fn quadratic_x2_minus_sqrt2_roots_vanish() {
        let mut session = q_session();
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(-1), sqrt2_algext().into_expr()]),
        ]);
        let p = poly_alg_from_expr(&e).unwrap();
        session = FieldSession::new(infer_field(&p).unwrap());
        let rs = quadratic_roots(&mut session, &p, &Var::from("x")).unwrap();
        assert_eq!(rs.len(), 2);
        for r in &rs {
            verify_root(&p, &Var::from("x"), r);
        }
    }

    #[serial]
    #[test]
    fn roots_quadratic_x2_minus_2() {
        let x = PolyAlgExt::ring_var("x");
        let mut session = q_session();
        let p = x
            .try_pow(2)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(&session, 2)))
            .unwrap();
        let rs = poly_algext_roots(&p, &Var::from("x")).unwrap();
        assert_eq!(rs.len(), 2);
        for r in &rs {
            verify_root(&p, &Var::from("x"), r);
        }
    }

    #[serial]
    #[test]
    fn roots_quadratic_x2_minus_sqrt2_over_k() {
        let mut session = q_session();
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(-1), sqrt2_algext().into_expr()]),
        ]);
        let p = poly_alg_from_expr(&e).unwrap();
        let rs = poly_algext_roots(&p, &Var::from("x")).unwrap();
        assert_eq!(rs.len(), 2);
        let sqrt2 = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap());
        for r in &rs {
            let sq = r.coeff_mul(r).unwrap();
            session = FieldSession::new(infer_field(&p).unwrap());
            let (sq, s2) = session.align(&sq, &sqrt2).unwrap();
            assert!(
                sq.coeff_sub(&s2).unwrap().coeff_is_zero(),
                "root^2 should equal sqrt2"
            );
            verify_root(&p, &Var::from("x"), r);
        }
    }

    #[serial]
    #[test]
    fn roots_cubic_t3_minus_2() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(3)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(&session, 2)))
            .unwrap();
        let rs = poly_algext_roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 1);
        verify_root(&p, &Var::from("t"), &rs[0]);
    }

    #[serial]
    #[test]
    #[ignore = "quartic resolvent Cardano exceeds 10s; needs common-field optimization"]
    fn roots_quartic_t4_plus_t_plus_1() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)).try_mul(&t).unwrap())
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)))
            .unwrap();
        let rs = poly_algext_roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 4);
        for r in &rs {
            verify_root(&p, &Var::from("t"), r);
        }
    }

    #[serial]
    #[test]
    fn roots_biquadratic_t4_minus_2() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(&session, 2)))
            .unwrap();
        let rs = poly_algext_roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 4);
        for r in &rs {
            verify_root(&p, &Var::from("t"), r);
        }
    }
}
