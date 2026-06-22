//! Univariate exact roots of `Poly<AlgExtCPolyCoeff>` (P2-1/6, P3-6).
//!
//! ```text
//! infer K from coefficients → normalize → monic → FieldSession::new(K) → deg dispatch
//!   1: linear
//!   2: quadratic + sqrt(Δ)
//!   3: pure t³+a₀ → β·ω^k (PR-C′) or Cardano+deflate
//!   4: biquadratic or resolvent R(z)=z³−pz²−4rz+(4pr−q²) → quadratics (PR-D′)
//! ```

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::sync::Arc;

use giac_poly::{PolyCoeff, Var};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;

use super::alg_ext::{algext_cube_root, algext_square_roots, AlgExtData};
use super::alg_ext_c::AlgExtCData;
use super::ext_tower::ExtensionField;
use super::field_arith::{coords_to_expr, pad_to_len, rationalize_poly1, CoordsQ};
use super::field_session::FieldSession;
use super::poly::PolyAlgExt;
use super::poly_alg_coeff::AlgExtCPolyCoeff;
use crate::context::Context;

/// Exact roots of univariate `p` w.r.t. `var` over the coefficient field of `p`.
/// **Stable (bounded)** — exact AlgExtC roots deg 1–4; quartic resolvent gap
pub fn poly_algext_roots(p: &PolyAlgExt, var: &Var) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let ambient = infer_field(p)?;
    let session = FieldSession::new(Arc::clone(&ambient));
    poly_algext_roots_in_session(p, var, &session)
}

/// Like [`poly_algext_roots`] using `ctx` session caches (R5 solve / eval).
/// **Stable (bounded)** — roots with Context session
pub fn poly_algext_roots_for_ctx(
    p: &PolyAlgExt,
    var: &Var,
    ctx: &Context,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let ambient = infer_field(p)?;
    let session = ctx.session().fork_ambient(ambient);
    poly_algext_roots_in_session(p, var, &session)
}

// **Pipeline private** — shared roots dispatch
fn poly_algext_roots_in_session(
    p: &PolyAlgExt,
    var: &Var,
    session: &FieldSession,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let p = normalize_coeffs(p, session)?;
    let monic = monic_univariate(&p, var, session)?;
    let d = monic.degree_wrt(var);
    match d {
        0 => {
            if monic.is_zero() {
                Ok(vec![])
            } else {
                Err(EvalError::TypeError("constant has no roots"))
            }
        }
        1 => linear_root(session, &monic, var),
        2 => quadratic_roots(session, &monic, var),
        3 => cubic_roots(session, &monic, var),
        4 => quartic_roots(session, &monic, var),
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
    session: &FieldSession,
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
    session: &FieldSession,
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
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let b = session.lift(&coeff_at(p, var, 1, session))?;
    let c = session.lift(&coeff_at(p, var, 0, session))?;
    let parent = session.working();
    let one = parent.one_coords();
    let b_c = embed_real_coords_for_parent(session, &parent, &b)?;
    let c_c = embed_real_coords_for_parent(session, &parent, &c)?;
    let field = if parent.dimension() == 1 {
        ExtensionField::adjoin_irreducible_over_q(vec![
            Ratio::one(),
            pad_to_len(&b_c, 1)[0].clone(),
            pad_to_len(&c_c, 1)[0].clone(),
        ])?
    } else {
        ExtensionField::adjoin_irreducible_parent_coeffs(
            &parent,
            vec![one, b_c, c_c],
        )?
    };
    let r1 = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&AlgExtData::from_field_coords(
        Arc::clone(&field),
        coords_to_expr(&field.generator_coords())?,
    )?)?);
    session.bump_to(&field);
    let b_lift = session.lift(&b)?;
    let r1_lift = session.lift(&r1)?;
    let r2 = session.add(&session.neg(&b_lift)?, &session.neg(&r1_lift)?)?;
    let _ = var;
    Ok(vec![r1_lift, r2])
}

// **Pipeline private** — quadratic roots via √Δ (no x²+bx+c adjoin layer)
fn quadratic_roots_formula(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let b = session.lift(&coeff_at(p, var, 1, session))?;
    let c = session.lift(&coeff_at(p, var, 0, session))?;
    let half = session.half()?;
    let b2 = session.mul(&b, &b)?;
    let four = session.int(4)?;
    let four_c = session.mul(&four, &c)?;
    let disc = session.add(&b2, &session.neg(&four_c)?)?;
    let sqrt_d = sqrt_disc(session, &disc)?;
    let neg_b = session.neg(&b)?;
    let neg_b_half = session.mul(&neg_b, &half)?;
    let sh = session.mul(&sqrt_d, &half)?;
    let r1 = session.add(&neg_b_half, &sh)?;
    let neg_sh = session.neg(&sh)?;
    let r2 = session.add(&neg_b_half, &neg_sh)?;
    let _ = var;
    Ok(vec![r1, r2])
}

// **Pipeline private** — negative constant in ℚ ⊂ K (for Δ<0 guard)
fn is_negative_rational(c: &AlgExtCPolyCoeff) -> bool {
    let inner = c.as_inner();
    if !inner.im.iter().all(|e| e.is_zero()) {
        return false;
    }
    let Ok(re) = rationalize_poly1(&inner.re) else {
        return false;
    };
    if inner.field.dimension() != 1 {
        return false;
    }
    pad_to_len(&re, 1)[0] < Ratio::zero()
}

// **Pipeline private** — sqrt(Δ) via session adjoin; imaginary branch when Δ<0 in K
fn sqrt_disc(
    session: &FieldSession,
    disc: &AlgExtCPolyCoeff,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    if is_negative_rational(disc) {
        let abs = session.neg(disc)?;
        let beta = session.adjoin_sqrt(&abs)?;
        return mul_i(session, &beta);
    }
    if disc.as_inner().im.iter().all(|e| e.is_zero()) {
        if let Ok(beta) = session.adjoin_sqrt(disc) {
            let sq = session.mul(&beta, &beta)?;
            let (sq_a, d_a) = session.align(&sq, disc)?;
            if sq_a.coeff_sub(&d_a).unwrap().coeff_is_zero()
                && beta.as_inner().im.iter().all(|e| e.is_zero())
            {
                return Ok(beta);
            }
        }
    }
    let neg_disc = session.neg(disc)?;
    let beta = session.adjoin_sqrt(&neg_disc)?;
    mul_i(session, &beta)
}

// **Pipeline private** — formal i times real z on session working field
fn mul_i(session: &FieldSession, z: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
    session.mul_formal_i(z)
}

// **Pipeline private** — Euler sqrt: real adjoin for u>0, i·√(−u) for u<0 (approx sign in K)
fn pick_sqrt_euler(
    session: &FieldSession,
    u: &AlgExtCPolyCoeff,
    force_imaginary: bool,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let u = session.lift(u)?;
    if u.coeff_is_zero() {
        return Ok(session.zero());
    }
    let use_imag = force_imaginary
        || is_negative_rational(&u)
        || approx_real_sign(session, &u).is_some_and(|s| s < 0.0);
    if use_imag {
        let abs = session.neg(&u)?;
        let beta = session.adjoin_sqrt(&abs)?;
        return mul_i(session, &beta);
    }
    session.sqrt_principal(&u)
}

// **Pipeline private** — rough real embedding sign (Euler sqrt branch only; deg ≤ 6)
fn approx_real_sign(session: &FieldSession, c: &AlgExtCPolyCoeff) -> Option<f64> {
    use num_traits::ToPrimitive;

    let inner = c.as_inner();
    if !inner.im.iter().all(|e| e.is_zero()) {
        return None;
    }
    let dim = inner.field.dimension();
    if dim > 6 {
        return None;
    }
    let coords = rationalize_poly1(&inner.re).ok()?;
    let coords = pad_to_len(&coords, dim);
    let mp = session.flatten_min_poly_over_q(&inner.field).ok()?;
    let roots = all_real_roots_minpoly(&mp);
    if roots.is_empty() {
        return None;
    }
    let vals: Vec<f64> = roots
        .iter()
        .map(|theta| {
            let mut v = 0.0;
            for c in &coords {
                v = v * theta + c.to_f64().unwrap_or(0.0);
            }
            v
        })
        .collect();
    if vals.iter().all(|v| *v < -1e-8) {
        return Some(-1.0);
    }
    if vals.iter().all(|v| *v > 1e-8) {
        return Some(1.0);
    }
    vals.into_iter()
        .max_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap())
}

// **Pipeline private** — real roots of monic minpoly (high-first, deg 1–3)
fn all_real_roots_minpoly(mp: &[Ratio<BigInt>]) -> Vec<f64> {
    use num_traits::ToPrimitive;

    let d = mp.len().checked_sub(1).unwrap_or(0);
    let coeffs: Vec<f64> = mp.iter().map(|r| r.to_f64()).collect::<Option<_>>().unwrap_or_default();
    match d {
        1 => {
            let a = coeffs[0];
            let b = coeffs[1];
            if a.abs() < 1e-15 {
                vec![]
            } else {
                vec![-b / a]
            }
        }
        2 => {
            let (a, b, c) = (coeffs[0], coeffs[1], coeffs[2]);
            let disc = b * b - 4.0 * a * c;
            if disc < 0.0 {
                vec![]
            } else {
                let s = disc.sqrt();
                vec![(-b + s) / (2.0 * a), (-b - s) / (2.0 * a)]
            }
        }
        3 => {
            let mut found = Vec::new();
            for start in [-3.0f64, 0.0, 3.0] {
                let mut x = start;
                for _ in 0..64 {
                    let f = coeffs[0] * x.powi(3) + coeffs[1] * x.powi(2) + coeffs[2] * x + coeffs[3];
                    let fp = 3.0 * coeffs[0] * x.powi(2) + 2.0 * coeffs[1] * x + coeffs[2];
                    if fp.abs() < 1e-15 {
                        break;
                    }
                    x -= f / fp;
                }
                if found.iter().all(|r: &f64| (r - x).abs() > 1e-5) {
                    found.push(x);
                }
            }
            found
        }
        _ => vec![],
    }
}

// **Pipeline private** — `cubic_roots`
fn cubic_roots(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let a2 = coeff_at(p, var, 2, session);
    let a1 = coeff_at(p, var, 1, session);
    if a2.coeff_is_zero() && a1.coeff_is_zero() {
        return pure_cubic_roots(session, &coeff_at(p, var, 0, session));
    }
    let r0 = one_cubic_root(session, p, var)?;
    let mut roots = vec![r0.clone()];
    let quad = deflate_monic(session, p, var, &r0)?;
    roots.extend(quadratic_roots(session, &quad, var)?);
    dedup_roots(roots)
}

// **Pipeline private** — monic t³+a₀=0: β=∛(−a₀), roots β·ω^k
fn pure_cubic_roots(
    session: &FieldSession,
    a0: &AlgExtCPolyCoeff,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let neg_a0 = session.neg(a0)?;
    let beta = session.adjoin_cbrt(&neg_a0)?;
    let omega = session.adjoin_primitive_cube_root_of_unity()?;
    let beta = session.lift(&beta)?;
    let omega2 = session.mul(&omega, &omega)?;
    let r1 = session.mul(&omega, &beta)?;
    let r2 = session.mul(&omega2, &beta)?;
    dedup_roots(vec![beta, r1, r2])
}

// **Pipeline private** — merge roots modulo minpoly
fn dedup_roots(roots: Vec<AlgExtCPolyCoeff>) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let mut out = Vec::new();
    for r in roots {
        if !out
            .iter()
            .any(|x: &AlgExtCPolyCoeff| x.eq_mod(&r).unwrap_or(false))
        {
            out.push(r);
        }
    }
    Ok(out)
}

// **Pipeline private** — `one_cubic_root`
fn one_cubic_root(
    session: &FieldSession,
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
    let p_num = session.add(&a1, &a2_shift)?;
    let p_dep = p_num;
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
        let r = session.adjoin_cbrt(&neg_q)?;
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
    if is_negative_rational(&delta) {
        return casus_adjoin_cubic_root(session, &p_dep, &q_dep, &shift);
    }
    let sqrt_delta = sqrt_disc(session, &delta)?;
    if !sqrt_delta
        .as_inner()
        .im
        .iter()
        .all(|e| e.is_zero())
    {
        return casus_adjoin_cubic_root(session, &p_dep, &q_dep, &shift);
    }
    let q_half = session.mul(&q_dep, &half)?;
    let neg_q_half = session.neg(&q_half)?;
    let sqrt_lift = session.lift(&sqrt_delta)?;
    let u_arg = session.add(&neg_q_half, &sqrt_lift)?;
    let u = session.adjoin_cbrt(&u_arg)?;
    let neg_sqrt = session.neg(&sqrt_lift)?;
    let v_arg = session.add(&neg_q_half, &neg_sqrt)?;
    let v = session.adjoin_cbrt(&v_arg)?;
    let (u, v) = session.align(&u, &v)?;
    let uv = session.add(&u, &v)?;
    let shift_lift = session.lift(&shift)?;
    session.add(&uv, &shift_lift)
}

// **Pipeline private** — casus: adjoin one root of monic t³+pt+q (Δ<0) directly
fn casus_adjoin_cubic_root(
    session: &FieldSession,
    p: &AlgExtCPolyCoeff,
    q: &AlgExtCPolyCoeff,
    shift: &AlgExtCPolyCoeff,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let parent = session.working();
    let zero = parent.zero_coords();
    let one = parent.one_coords();
    let p_c = embed_real_coords_for_parent(session, &parent, p)?;
    let q_c = embed_real_coords_for_parent(session, &parent, q)?;
    let field = if parent.dimension() == 1 {
        ExtensionField::adjoin_irreducible_over_q(vec![
            Ratio::one(),
            Ratio::zero(),
            pad_to_len(&p_c, 1)[0].clone(),
            pad_to_len(&q_c, 1)[0].clone(),
        ])?
    } else {
        ExtensionField::adjoin_irreducible_parent_coeffs(
            &parent,
            vec![one, zero.clone(), p_c, q_c],
        )?
    };
    let root = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&AlgExtData::from_field_coords(
        Arc::clone(&field),
        coords_to_expr(&field.generator_coords())?,
    )?)?);
    session.bump_to(&field);
    let shift_lift = session.lift(shift)?;
    session.add(&root, &shift_lift)
}

fn embed_real_coords_for_parent(
    session: &FieldSession,
    parent: &Arc<ExtensionField>,
    c: &AlgExtCPolyCoeff,
) -> Result<CoordsQ, EvalError> {
    let inner = c.as_inner();
    if !inner.im.iter().all(|e| e.is_zero()) {
        return Err(EvalError::TypeError("expected real coefficient"));
    }
    let re = rationalize_poly1(&inner.re)?;
    let aligned = session.align_elements(
        &inner.field,
        &re,
        parent,
        &parent.zero_coords(),
    )?;
    Ok(pad_to_len(&aligned.left, parent.dimension()))
}

// **Pipeline private** — `quartic_roots`
fn quartic_roots(
    session: &FieldSession,
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
    let p2 = session.lift(&coeff_at(&dep, var, 2, session))?;
    let p1 = session.lift(&coeff_at(&dep, var, 1, session))?;
    let p0 = session.lift(&coeff_at(&dep, var, 0, session))?;
    let res = build_resolvent_cubic(session, &p2, &p1, &p0)?;
    // Resolvent roots: one Cardano/casus root + deflate quadratic (formula, not
    // `quadratic_roots` adjoin). Full `cubic_roots` over-splits the tower and
    // breaks subsequent `adjoin_sqrt` (ε²≠u); see adjoin_sqrt_after_* tests.
    let z_var = Var::from("_z");
    let z0 = one_cubic_root(session, &res, &z_var)?;
    let quad = deflate_monic(session, &res, &z_var, &z0)?;
    let mut z_roots = vec![z0];
    z_roots.extend(quadratic_roots_formula(session, &quad, &z_var)?);
    let z_roots = dedup_roots(z_roots)?;
    if z_roots.len() < 3 {
        return Err(EvalError::NotImplemented("quartic resolvent"));
    }
    let all = euler_depressed_quartic_roots(session, &dep, var, &z_roots)?;
    let all = dedup_roots(all)?;
    if all.len() != 4 {
        return Err(EvalError::NotImplemented("quartic roots"));
    }
    Ok(all)
}

// **Pipeline private** — `biquadratic_roots`
fn biquadratic_roots(
    session: &FieldSession,
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
    session: &FieldSession,
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

// **Pipeline private** — Ferrari resolvent R(z)=z³−pz²−4rz+(4pr−q²) on session
fn build_resolvent_cubic(
    session: &FieldSession,
    p: &AlgExtCPolyCoeff,
    q: &AlgExtCPolyCoeff,
    r: &AlgExtCPolyCoeff,
) -> Result<PolyAlgExt, EvalError> {
    let four = session.int(4)?;
    let z = PolyAlgExt::ring_var(Var::from("_z"));
    let r4 = session.mul(r, &four)?;
    let rq = session.mul(r, p)?;
    let rqp4 = session.mul(&rq, &four)?;
    let qq = session.mul(q, q)?;
    let neg_qq = session.neg(&qq)?;
    let const_term = session.add(&rqp4, &neg_qq)?;
    Ok(z.try_pow(3)?
        .try_sub(&PolyAlgExt::ring_constant(p.clone()).try_mul(&z.try_pow(2)?)?)?
        .try_sub(&PolyAlgExt::ring_constant(r4).try_mul(&z)?)?
        .try_add(&PolyAlgExt::ring_constant(const_term))?)
}

// **Pipeline private** — Euler resolvent: four roots from three resolvent zeros α,β,γ
fn euler_depressed_quartic_roots(
    session: &FieldSession,
    dep: &PolyAlgExt,
    var: &Var,
    z_roots: &[AlgExtCPolyCoeff],
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    if z_roots.len() < 3 {
        return Err(EvalError::NotImplemented("quartic euler"));
    }
    let q = session.lift(&coeff_at(dep, var, 1, session))?;
    let checkpoint = session.working();
    let permutations = [
        [0usize, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for perm in permutations {
        session.set_working(&checkpoint);
        if let Ok(roots) = euler_depressed_quartic_roots_perm(session, dep, var, &q, z_roots, perm) {
            return Ok(roots);
        }
    }
    Err(EvalError::NotImplemented("quartic euler"))
}

fn euler_depressed_quartic_roots_perm(
    session: &FieldSession,
    dep: &PolyAlgExt,
    var: &Var,
    q: &AlgExtCPolyCoeff,
    z_roots: &[AlgExtCPolyCoeff],
    perm: [usize; 3],
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let alpha = session.lift(&z_roots[perm[0]])?;
    let beta = session.lift(&z_roots[perm[1]])?;
    let gamma = session.lift(&z_roots[perm[2]])?;
    let checkpoint = session.working();
    let flip_g_opts: &[bool] = if q.coeff_is_zero() { &[false, true] } else { &[false] };
    for flip_a in [false, true] {
        for flip_b in [false, true] {
            for &flip_g in flip_g_opts {
                session.set_working(&checkpoint);
                let sa = pick_sqrt_euler(session, &alpha, flip_a)?;
                let sb = pick_sqrt_euler(session, &beta, flip_b)?;
                let roots = if q.coeff_is_zero() {
                    let sg = pick_sqrt_euler(session, &gamma, flip_g)?;
                    euler_four_roots_from_triple(session, q, &sa, &sb, &sg)?
                } else {
                    let sqrt_g = euler_derived_sqrt_gamma(session, q, &sa, &sb)?;
                    euler_four_roots_from_triple(session, q, &sa, &sb, &sqrt_g)?
                };
                if roots_all_vanish(dep, var, &roots) {
                    return Ok(roots);
                }
            }
        }
    }
    Err(EvalError::NotImplemented("quartic euler perm"))
}

// **Pipeline private** — √γ = −q / (√α·√β) for depressed x⁴+px²+qx+r (q≠0)
fn euler_derived_sqrt_gamma(
    session: &FieldSession,
    q: &AlgExtCPolyCoeff,
    sqrt_a: &AlgExtCPolyCoeff,
    sqrt_b: &AlgExtCPolyCoeff,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let num = session.neg(q)?;
    let denom = session.mul(sqrt_a, sqrt_b)?;
    session.div(&num, &denom)
}

fn euler_four_roots_from_triple(
    session: &FieldSession,
    _q: &AlgExtCPolyCoeff,
    sqrt_a: &AlgExtCPolyCoeff,
    sqrt_b: &AlgExtCPolyCoeff,
    sqrt_g: &AlgExtCPolyCoeff,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let half = session.half()?;
    let sa = session.lift(sqrt_a)?;
    let sb = session.lift(sqrt_b)?;
    let sg = session.lift(sqrt_g)?;
    let neg_sa = session.neg(&sa)?;
    let neg_sb = session.neg(&sb)?;
    let neg_sg = session.neg(&sg)?;
    let mut out = Vec::new();
    for (ca, cb, cg) in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)] {
        let ta = if ca > 0 { &sa } else { &neg_sa };
        let tb = if cb > 0 { &sb } else { &neg_sb };
        let tg = if cg > 0 { &sg } else { &neg_sg };
        let mid = session.add(tb, tg)?;
        let sum = session.add(ta, &mid)?;
        out.push(session.mul(&sum, &half)?);
    }
    Ok(out)
}

fn roots_all_vanish(dep: &PolyAlgExt, var: &Var, roots: &[AlgExtCPolyCoeff]) -> bool {
    roots.iter().all(|r| root_vanishes(dep, var, r))
}

// **Pipeline private** — verify root vanishes (fresh session, like verify_root)
fn root_vanishes(p: &PolyAlgExt, var: &Var, root: &AlgExtCPolyCoeff) -> bool {
    let ambient = match infer_field(p) {
        Ok(f) => f,
        Err(_) => return false,
    };
    let mut session = FieldSession::new(ambient);
    let p = match normalize_coeffs(p, &session) {
        Ok(p) => p,
        Err(_) => return false,
    };
    let p = match monic_univariate(&p, var, &session) {
        Ok(p) => p,
        Err(_) => return false,
    };
    eval_vanishes(&session, &p, var, root)
}

// **Pipeline private** — quick vanishing check on working session
fn eval_vanishes(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
    root: &AlgExtCPolyCoeff,
) -> bool {
    session.bump_to(&root.as_inner().field);
    let mut val = session.zero();
    for (m, c) in &p.terms {
        let exp = m.exp_of(var);
        let mut pow = session.one();
        for _ in 0..exp {
            pow = match session.mul(root, &pow) {
                Ok(p) => p,
                Err(_) => return false,
            };
        }
        let c = match session.lift(c) {
            Ok(c) => c,
            Err(_) => return false,
        };
        let (c, pow) = match session.align(&c, &pow) {
            Ok(x) => x,
            Err(_) => return false,
        };
        let term = match session.mul(&c, &pow) {
            Ok(t) => t,
            Err(_) => return false,
        };
        val = match session.add(&val, &term) {
            Ok(v) => v,
            Err(_) => return false,
        };
    }
    val.coeff_is_zero()
}

// **Pipeline private** — `split_depressed_quartic` (legacy Ferrari; kept for tests)
fn split_depressed_quartic(
    session: &FieldSession,
    dep: &PolyAlgExt,
    var: &Var,
    z: &AlgExtCPolyCoeff,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let q = session.lift(&coeff_at(dep, var, 1, session))?;
    let r = session.lift(&coeff_at(dep, var, 0, session))?;
    let z_lift = session.lift(z)?;
    session.bump_to(&z_lift.as_inner().field);
    let four = session.int(4)?;
    let zz = session.mul(&z_lift, &z_lift)?;
    let r4 = session.mul(&r, &four)?;
    let neg_r4 = session.neg(&r4)?;
    let disc = session.add(&zz, &neg_r4)?;
    let m = session.adjoin_sqrt(&disc)?;
    let half = session.half()?;
    let m_lift = session.lift(&m)?;
    let neg_m = session.neg(&m_lift)?;
    let mut out = Vec::new();
    for mp in [
        session.add(&z_lift, &m_lift)?,
        session.add(&z_lift, &neg_m)?,
    ] {
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
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
    root: &AlgExtCPolyCoeff,
) -> Result<PolyAlgExt, EvalError> {
    session.bump_to(&root.as_inner().field);
    let n = p.degree_wrt(var) as usize;
    let mut qs = vec![session.zero(); n];
    qs[n - 1] = session.one();
    for k in (0..n - 1).rev() {
        let ak = session.lift(&coeff_at(p, var, (k + 1) as u64, session))?;
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
        &normalize_coeffs(p, &session).expect("normalize"),
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
    //! Test tiers: **B** — `.doc/test-writing-spec.md` §6.1
    //! Inventory: `.doc/giac-core-algebra-api-stability.md` (`poly_roots.rs` tests)
    use std::sync::Arc;

    use giac_poly::{PolyCoeff, Var};
    use num_rational::Ratio;

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

    // **B** — one Cardano root vanishes on t³−2.
    #[test]
    fn cubic_one_root_vanishes() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(3)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(&session, 2)))
            .unwrap();
        let r = one_cubic_root(&session, &p, &Var::from("t")).unwrap();
        verify_root(&p, &Var::from("t"), &r);
    }

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

    #[test]
    fn quadratic_x2_minus_sqrt2_roots_vanish() {
        let mut session = q_session();
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(-1), sqrt2_algext().into_expr()]),
        ]);
        let p = poly_alg_from_expr(&e).unwrap();
        session = FieldSession::new(infer_field(&p).unwrap());
        let rs = quadratic_roots(&session, &p, &Var::from("x")).unwrap();
        assert_eq!(rs.len(), 2);
        for r in &rs {
            verify_root(&p, &Var::from("x"), r);
        }
    }

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

    // **B** — roots_cubic_t3_minus_2: three roots β, ωβ, ω²β (PR-C′).
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
        assert_eq!(rs.len(), 3);
        for r in &rs {
            verify_root(&p, &Var::from("t"), r);
        }
    }

    // **B** — PR-E′: all three resolvent roots of z³−4z−1.
    // FIX-Z3M4Z1
    #[test]
    fn resolvent_cubic_all_roots_vanish() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let rs = cubic_roots(&session, &res, &Var::from("_z")).unwrap();
        assert_eq!(rs.len(), 3);
        for r in &rs {
            verify_root(&res, &Var::from("_z"), r);
        }
    }

    // **B** — PR-E′: one Cardano/casus root of resolvent z³−4z−1.
    // FIX-Z3M4Z1
    #[test]
    fn resolvent_one_cubic_root_z3_minus_4z_minus_1() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let z = one_cubic_root(&session, &res, &Var::from("_z")).unwrap();
        verify_root(&res, &Var::from("_z"), &z);
    }

    // **B** — PR-D′ golden: depressed (p,q,r)=(0,1,1) → R(z)=z³−4z−1.
    // FIX-T4P1 · Lean: Giac.Examples.T4PlusTPlus1.t4_plus_t_plus_1_resolvent
    #[test]
    fn resolvent_golden_t4_plus_t_plus_1() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let z = Var::from("_z");
        let c3 = coeff_at(&res, &z, 3, &session);
        let c2 = coeff_at(&res, &z, 2, &session);
        let c1 = coeff_at(&res, &z, 1, &session);
        let c0 = coeff_at(&res, &z, 0, &session);
        let one = session.one();
        let neg_four = session.int(-4).unwrap();
        let neg_one = session.int(-1).unwrap();
        assert!(c3.eq_mod(&one).unwrap());
        assert!(c2.eq_mod(&session.zero()).unwrap());
        assert!(c1.eq_mod(&neg_four).unwrap());
        assert!(c0.eq_mod(&neg_one).unwrap());
    }

    // **B** — PR-D′: t⁴+t+1 depression yields q=1, r=1, p=0.
    // FIX-T4P1 · Lean: Giac.Examples.T4PlusTPlus1.t4_plus_t_plus_1_depressed_coeffs
    #[test]
    fn depressed_t4_plus_t_plus_1_coeffs() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)).try_mul(&t).unwrap())
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)))
            .unwrap();
        let shift = session.zero();
        let dep = depress_quartic(&session, &p, &Var::from("t"), &shift).unwrap();
        let var = Var::from("t");
        assert!(coeff_at(&dep, &var, 3, &session).coeff_is_zero());
        assert!(coeff_at(&dep, &var, 2, &session).coeff_is_zero());
        assert!(coeff_at(&dep, &var, 1, &session).eq_mod(&session.int(1).unwrap()).unwrap());
        assert!(coeff_at(&dep, &var, 0, &session).eq_mod(&session.int(1).unwrap()).unwrap());
    }

    // **B** — PR-E′: adjoin_sqrt preserves ε²=u on one resolvent root (L₃ only).
    // FIX-SQRT-RESOLVENT · Lean: Giac.Tower.SqrtInField (future)
    #[test]
    fn adjoin_sqrt_squares_one_resolvent_root() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let z = one_cubic_root(&session, &res, &Var::from("_z")).unwrap();
        let z = session.lift(&z).unwrap();
        let sa = session.adjoin_sqrt(&z).unwrap();
        let sq = session.mul(&sa, &sa).unwrap();
        let (sq, z_a) = session.align(&sq, &z).unwrap();
        assert!(
            sq.coeff_sub(&z_a).unwrap().coeff_is_zero(),
            "adjoin_sqrt generator squared should equal input"
        );
    }

    // **B** — PR-E′: deflate alone does not break adjoin_sqrt (regression).
    // FIX-SQRT-RESOLVENT
    #[test]
    fn adjoin_sqrt_after_resolvent_deflate_only() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let z_var = Var::from("_z");
        let z0 = one_cubic_root(&session, &res, &z_var).unwrap();
        let _quad = deflate_monic(&session, &res, &z_var, &z0).unwrap();
        let z = session.lift(&z0).unwrap();
        let sa = session.adjoin_sqrt(&z).unwrap();
        let sq = session.mul(&sa, &sa).unwrap();
        let (sq, z_a) = session.align(&sq, &z).unwrap();
        assert!(
            sq.coeff_sub(&z_a).unwrap().coeff_is_zero(),
            "sqrt square after deflate only"
        );
    }

    // **B** — PR-E′: Euler four roots vanish on depressed t⁴+t+1.
    // FIX-T4P1 · Lean: Giac.Examples.T4PlusTPlus1.t4_plus_t_plus_1_euler_vanishes (partial)
    #[test]
    #[ignore = "tower sqrt embedding: adjoin_sqrt after resolvent split breaks ε²=u"]
    fn euler_four_roots_vanish() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)).try_mul(&t).unwrap())
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)))
            .unwrap();
        let shift = session.zero();
        let dep = depress_quartic(&session, &p, &Var::from("t"), &shift).unwrap();
        let p2 = session.lift(&coeff_at(&dep, &Var::from("t"), 2, &session)).unwrap();
        let p1 = session.lift(&coeff_at(&dep, &Var::from("t"), 1, &session)).unwrap();
        let p0 = session.lift(&coeff_at(&dep, &Var::from("t"), 0, &session)).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let z_roots = cubic_roots(&session, &res, &Var::from("_z")).unwrap();
        let rs = euler_depressed_quartic_roots(&session, &dep, &Var::from("t"), &z_roots).unwrap();
        assert_eq!(rs.len(), 4);
        for r in &rs {
            verify_root(&dep, &Var::from("t"), r);
        }
    }

    // FIX-T4P1 · Lean: Giac.Examples.T4PlusTPlus1.t4_plus_t_plus_1_four_roots (partial)
    #[test]
    #[ignore = "tower sqrt embedding: adjoin_sqrt after resolvent split breaks ε²=u"]
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
