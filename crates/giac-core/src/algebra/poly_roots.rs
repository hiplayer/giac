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

use super::alg_ext::AlgExtData;
use super::alg_ext_c::AlgExtCData;
use super::ext_tower::{self, ExtensionField};
use super::field_arith::{coords_to_expr, pad_to_len, rationalize_poly1, CoordsQ};
use super::field_session::{
    coeff_from_coords, coords_in_field, is_negative_rational, FieldSession, POLY_ROOTS_DIM_HARD,
    POLY_ROOTS_DIM_QUARTIC_OUT, POLY_ROOTS_DIM_RESOLVENT,
};
use super::galois_automorphism;
use super::poly::PolyAlgExt;
use super::poly_alg_coeff::AlgExtCPolyCoeff;
use crate::context::Context;

/// Normalized monic univariate input over fixed ambient **K** (algorithm-expr-api §6.3 / 1B).
struct PolyInK {
    poly: PolyAlgExt,
    ambient: Arc<ExtensionField>,
}

impl PolyInK {
    // **Pipeline private** — infer K → normalize → monic; sole entry to roots_dispatch
    fn prepare(
        p: &PolyAlgExt,
        var: &Var,
        session: &FieldSession,
    ) -> Result<Self, EvalError> {
        let p = normalize_coeffs(p, session)?;
        let poly = monic_univariate(&p, var, session)?;
        Ok(Self {
            poly,
            ambient: Arc::clone(session.ambient()),
        })
    }

    // **Pipeline private** — prepare with fresh session from coefficient field of p
    fn prepare_with_session(p: &PolyAlgExt, var: &Var) -> Result<(Self, FieldSession), EvalError> {
        let ambient = infer_field(p)?;
        let session = FieldSession::new(ambient);
        let input = Self::prepare(p, var, &session)?;
        Ok((input, session))
    }
}

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
    let input = PolyInK::prepare(p, var, session)?;
    roots_dispatch(&input, var, session)
}

// **Pipeline private** — degree dispatch on prepared monic input
fn roots_dispatch(
    input: &PolyInK,
    var: &Var,
    session: &FieldSession,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    debug_assert!(Arc::ptr_eq(&input.ambient, session.ambient()));
    let monic = &input.poly;
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
        2 => quadratic_roots_formula(session, &monic, var),
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

// **Pipeline private** — sqrt(Δ) via session; imaginary branch when Δ<0 in ℚ ⊂ K
fn sqrt_disc(
    session: &FieldSession,
    disc: &AlgExtCPolyCoeff,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    if is_negative_rational(disc) {
        let abs = session.neg(disc)?;
        let beta = session.sqrt_in_field(&abs)?;
        return mul_i(session, &beta);
    }
    session.sqrt_in_field(disc)
}

// **Pipeline private** — formal i times real z on session working field
fn mul_i(session: &FieldSession, z: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
    session.mul_formal_i(z)
}

// **Pipeline private** — √γ = −q / (√α·√β) for depressed x⁴+px²+qx+r (q≠0)
fn cubic_depressed_parts(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<(AlgExtCPolyCoeff, AlgExtCPolyCoeff, AlgExtCPolyCoeff, AlgExtCPolyCoeff), EvalError> {
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
    let half = session.half()?;
    let twenty_seven = session.int(27)?;
    let qh = session.mul(&q_dep, &half)?;
    let qh2 = session.mul(&qh, &qh)?;
    let p2 = session.mul(&p_dep, &p_dep)?;
    let p3 = session.mul(&p2, &p_dep)?;
    let term = session.div(&p3, &twenty_seven)?;
    let delta = session.add(&qh2, &term)?;
    Ok((p_dep, q_dep, shift, delta))
}

fn finish_cubic_roots(
    session: &FieldSession,
    roots: Vec<AlgExtCPolyCoeff>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    if session.working().dimension() > POLY_ROOTS_DIM_HARD {
        return Err(EvalError::NotImplemented("poly roots dim bound"));
    }
    canonical_sort_roots(session, dedup_roots(roots)?)
}

// **Pipeline private** — F2 A′: monic t³+a₀ → β·ω^k (ω adjoin once).
fn f2_pure_cubic_roots_in_session(
    session: &FieldSession,
    a0: &AlgExtCPolyCoeff,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    pure_cubic_roots(session, a0)
}

// **Pipeline private** — F2 deflate path: one root via irreducible t³+p_dep·t+q (no Cardano stack).
fn f2_one_cubic_root_for_deflate(
    session: &FieldSession,
    p_dep: &AlgExtCPolyCoeff,
    q_dep: &AlgExtCPolyCoeff,
    shift: &AlgExtCPolyCoeff,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    casus_adjoin_cubic_root(session, p_dep, q_dep, shift)
}

// **Pipeline private** — F2: z₀ + deflate quadratic; √Δ only via sqrt_in_field (F1).
fn f2_split_cubic_via_deflate_in_session(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
    p_dep: &AlgExtCPolyCoeff,
    q_dep: &AlgExtCPolyCoeff,
    shift: &AlgExtCPolyCoeff,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let z0 = f2_one_cubic_root_for_deflate(session, p_dep, q_dep, shift)?;
    // ponytail: casus p≠0 (e.g. z³−4z−1, x³−x+1): ω·z₀ is not a root; deflate + F1 √Δ.
    let quad = deflate_monic(session, p, var, &z0)?;
    let mut roots = vec![session.lift(&z0)?];
    roots.extend(quadratic_roots_formula(session, &quad, var)?);
    finish_cubic_roots(session, roots)
}

// **Pipeline private** — unified cubic split (P3-6 §3 / F2).
fn split_monic_cubic_roots_in_session(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let a2 = coeff_at(p, var, 2, session);
    let a1 = coeff_at(p, var, 1, session);
    if a2.coeff_is_zero() && a1.coeff_is_zero() {
        return finish_cubic_roots(session, f2_pure_cubic_roots_in_session(session, &coeff_at(p, var, 0, session))?);
    }

    let (p_dep, q_dep, shift, _delta) = cubic_depressed_parts(session, p, var)?;

    if p_dep.coeff_is_zero() {
        let shift_lift = session.lift(&shift)?;
        let mut rs = f2_pure_cubic_roots_in_session(session, &q_dep)?;
        rs = rs
            .into_iter()
            .map(|r| session.add(&r, &shift_lift))
            .collect::<Result<Vec<_>, EvalError>>()?;
        return finish_cubic_roots(session, rs);
    }

    f2_split_cubic_via_deflate_in_session(session, p, var, &p_dep, &q_dep, &shift)
}

// **Pipeline private** — resolvent three roots in L (split_monic_cubic_roots_in_session)
fn resolvent_cubic_roots_in_session(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    split_monic_cubic_roots_in_session(session, p, var)
}

// **Pipeline private** — `cubic_roots` (delegates to F2 resolvent API)
fn cubic_roots(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    resolvent_cubic_roots_in_session(session, p, var)
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

// **Pipeline private** — lex sort roots by coords in working field (F3 canonical α,β,γ)
fn root_coords_lex_key(
    session: &FieldSession,
    r: &AlgExtCPolyCoeff,
) -> Result<Vec<Ratio<BigInt>>, EvalError> {
    let lifted = session.lift(r)?;
    embed_real_coords_for_parent(session, &session.working(), &lifted)
}

fn canonical_sort_roots(
    session: &FieldSession,
    roots: Vec<AlgExtCPolyCoeff>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let mut keyed: Vec<_> = roots
        .into_iter()
        .map(|r| {
            let key = root_coords_lex_key(session, &r)?;
            Ok((key, r))
        })
        .collect::<Result<Vec<_>, EvalError>>()?;
    keyed.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(keyed.into_iter().map(|(_, r)| r).collect())
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
    let (p_dep, q_dep, shift, delta) = cubic_depressed_parts(session, p, var)?;
    if p_dep.coeff_is_zero() {
        let neg_q = session.neg(&q_dep)?;
        let r = session.adjoin_cbrt(&neg_q)?;
        let r_lift = session.lift(&r)?;
        let shift_lift = session.lift(&shift)?;
        return session.add(&r_lift, &shift_lift);
    }
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
    let half = session.half()?;
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
        session.adjoin_irreducible(
            &parent,
            vec![
                Ratio::one(),
                Ratio::zero(),
                pad_to_len(&p_c, 1)[0].clone(),
                pad_to_len(&q_c, 1)[0].clone(),
            ],
        )?
    } else {
        session.adjoin_irreducible_parent_coeffs(
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
    let z_var = Var::from("_z");
    let z_roots = resolvent_cubic_roots_in_session(session, &res, &z_var)?;
    debug_assert!(
        session.working().dimension() <= POLY_ROOTS_DIM_RESOLVENT,
        "F2 resolvent stage dim bound"
    );
    if z_roots.len() < 3 {
        return Err(EvalError::NotImplemented("quartic resolvent"));
    }
    let all = euler_depressed_quartic_roots(session, &dep, var, &z_roots)?;
    let all = dedup_roots(all)?;
    if all.len() != 4 {
        return Err(EvalError::NotImplemented("quartic roots"));
    }
    if session.working().dimension() > POLY_ROOTS_DIM_HARD {
        // F4 exit gate — `.doc/issues/GIAC-poly-p3-6-roots-algorithm-spec.md` §5.4
        return Err(EvalError::NotImplemented("poly roots dim bound"));
    }
    Ok(all)
}

// **Pipeline private** — fallback: adjoin one quartic root, then solve deflated cubic.
fn quartic_roots_by_adjoin_deflate(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let root = adjoin_one_root_of_monic(session, p, var)?;
    let cubic = literal_monic_univariate(
        session,
        &monic_univariate(&deflate_monic(session, p, var, &root)?, var, session)?,
        var,
    )?;
    let root2 = adjoin_one_root_of_monic(session, &cubic, var)?;
    let quad = literal_monic_univariate(
        session,
        &monic_univariate(&deflate_monic(session, &cubic, var, &root2)?, var, session)?,
        var,
    )?;
    let root3 = adjoin_one_root_of_monic(session, &quad, var)?;
    let b = session.lift(&coeff_at(&quad, var, 1, session))?;
    let root3_lift = session.lift(&root3)?;
    let root4 = session.neg(&session.add(&b, &root3_lift)?)?;
    let out = vec![root, root2, root3, root4];
    let out = dedup_roots(out)?;
    if out.len() == 4 && roots_all_vanish(p, var, &out) {
        Ok(out)
    } else {
        Err(EvalError::NotImplemented("quartic adjoin-deflate fallback"))
    }
}

fn literal_monic_univariate(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<PolyAlgExt, EvalError> {
    let deg = p.degree_wrt(var);
    let mut out = PolyAlgExt::ring_var(var.clone()).try_pow(deg)?;
    for exp in 0..deg {
        let c = session.lift(&coeff_at(p, var, exp, session))?;
        if c.coeff_is_zero() {
            continue;
        }
        let mut term = PolyAlgExt::ring_constant(c);
        if exp > 0 {
            term = term.try_mul(&PolyAlgExt::ring_var(var.clone()).try_pow(exp)?)?;
        }
        out = out.try_add(&term)?;
    }
    Ok(out)
}

fn adjoin_one_root_of_monic(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let parent = session.working();
    let degree = p.degree_wrt(var);
    let field = if parent.is_base() {
        let mut min_poly = Vec::new();
        for exp in (0..=degree).rev() {
            if exp == degree {
                min_poly.push(Ratio::one());
                continue;
            }
            let c = session.lift(&coeff_at(p, var, exp, session))?;
            let re = rationalize_poly1(&c.as_inner().re)?;
            min_poly.push(pad_to_len(&re, 1)[0].clone());
        }
        session.adjoin_irreducible(&parent, min_poly)?
    } else {
        let mut blocks = Vec::new();
        for exp in (0..=degree).rev() {
            if exp == degree {
                blocks.push(parent.one_coords());
                continue;
            }
            let c = session.lift(&coeff_at(p, var, exp, session))?;
            let inner = c.as_inner();
            if !inner.im.iter().all(|e| e.is_zero()) {
                return Err(EvalError::TypeError("expected real coefficient"));
            }
            if !(Arc::ptr_eq(&inner.field, &parent) || *inner.field == *parent) {
                return Err(EvalError::TypeError("coefficient not lifted to parent"));
            }
            let re = rationalize_poly1(&inner.re)?;
            blocks.push(pad_to_len(&re, parent.dimension()));
        }
        ext_tower::build_adjoin_parent_coeffs(&parent, blocks)?
    };
    session.bump_to(&field);
    let root_data = AlgExtData::from_field_coords(
        Arc::clone(&field),
        coords_to_expr(&field.generator_coords())?,
    )?;
    let root = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&root_data)?);
    Ok(root)
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
    let u_roots = quadratic_roots_formula(session, &u, &Var::from("u"))?;
    let mut out = Vec::new();
    for ur in u_roots {
        out.extend(sqrt_branches(session, &ur)?);
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

// **Pipeline private** — F4′: pick β among conjugates with √(β/α) ∈ L(√α); else Galois σ; else blind adjoin.
fn sqrt_in_field_euler_second(
    session: &FieldSession,
    u: &AlgExtCPolyCoeff,
    candidates: &[AlgExtCPolyCoeff],
    sqrt_u: &AlgExtCPolyCoeff,
    z_roots: &[AlgExtCPolyCoeff],
) -> Result<AlgExtCPolyCoeff, EvalError> {
    for v in candidates {
        if let Some(sv) = session.try_sqrt_in_field(v)? {
            return Ok(sv);
        }
        let ratio = session.div(v, u)?;
        if let Some(sr) = session.try_sqrt_in_field(&ratio)? {
            return session.mul(&sr, sqrt_u);
        }
        if let Some(sr) = session.try_sqrt_in_field_shallow(&ratio)? {
            return session.mul(&sr, sqrt_u);
        }
    }
    if let Some(sv) = try_galois_sqrt_second_coeff(session, u, candidates, sqrt_u, z_roots)? {
        return Ok(sv);
    }
    if session.working().dimension() >= POLY_ROOTS_DIM_HARD {
        return Err(EvalError::NotImplemented("poly roots dim bound"));
    }
    // ponytail: first candidate still needs adjoin (ratio not square in L(√α))
    session.adjoin_sqrt_new(candidates.first().ok_or(EvalError::NotImplemented(
        "euler second sqrt",
    ))?)
}

// **Pipeline private** — F4′: embed coeffs in parent P and call `galois_automorphism`.
fn try_galois_sqrt_second_coeff(
    session: &FieldSession,
    u: &AlgExtCPolyCoeff,
    candidates: &[AlgExtCPolyCoeff],
    sqrt_u: &AlgExtCPolyCoeff,
    z_roots: &[AlgExtCPolyCoeff],
) -> Result<Option<AlgExtCPolyCoeff>, EvalError> {
    let field = session.working();
    let parent = match field.parent_field() {
        Some(p) => Arc::clone(p),
        None => return Ok(None),
    };
    let u_c = coords_in_field(session, u, &parent)?;
    let sqrt_u_c = coords_in_field(session, sqrt_u, &field)?;
    let z_coords: Result<Vec<_>, _> = z_roots
        .iter()
        .map(|z| coords_in_field(session, z, &parent))
        .collect();
    let z_coords = z_coords?;
    for v in candidates {
        let v_c = coords_in_field(session, v, &parent)?;
        if let Some(h) = galois_automorphism::try_galois_sqrt_second(
            &field, &u_c, &v_c, &sqrt_u_c, &z_coords,
        )? {
            return Ok(Some(coeff_from_coords(&field, &h)?));
        }
    }
    Ok(None)
}

// **Pipeline private** — F4′: one Euler branch with α = sorted[alpha_idx].
fn try_euler_one_alpha(
    session: &FieldSession,
    dep: &PolyAlgExt,
    var: &Var,
    sorted: &[AlgExtCPolyCoeff],
    alpha_idx: usize,
    q: &AlgExtCPolyCoeff,
    z_roots: &[AlgExtCPolyCoeff],
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let alpha = session.lift(&sorted[alpha_idx])?;
    let others: Vec<_> = (0..sorted.len())
        .filter(|&i| i != alpha_idx)
        .map(|i| session.lift(&sorted[i]))
        .collect::<Result<Vec<_>, _>>()?;
    if others.len() < 2 {
        return Err(EvalError::NotImplemented("quartic euler"));
    }
    let sa = session.sqrt_in_field(&alpha)?;
    let sb = sqrt_in_field_euler_second(session, &alpha, &others, &sa, z_roots)?;
    if q.coeff_is_zero() {
        let sg = sqrt_in_field_euler_second(
            session,
            &alpha,
            &[others[1].clone(), others[0].clone()],
            &sa,
            z_roots,
        )?;
        euler_four_roots_from_triple(session, q, &sa, &sb, &sg)
    } else {
        let sqrt_g = euler_derived_sqrt_gamma(session, q, &sa, &sb)?;
        euler_four_roots_from_triple(session, q, &sa, &sb, &sqrt_g)
    }
}

// **Pipeline private** — Euler resolvent: four roots from three resolvent zeros α,β,γ (F3 single path)
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
    let sorted = canonical_sort_roots(session, dedup_roots(z_roots.to_vec())?)?;
    if sorted.len() < 3 {
        return Err(EvalError::NotImplemented("quartic euler"));
    }
    let roots = try_euler_one_alpha(session, dep, var, &sorted, 0, &q, z_roots)?;
    if !roots_all_vanish(dep, var, &roots) {
        return Err(EvalError::NotImplemented("quartic euler"));
    }
    Ok(roots)
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
    val.as_inner()
        .eq_mod(session.zero().as_inner())
        .unwrap_or(false)
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

// **Pipeline private** — ±√z via session (replaces blind `algext_square_roots`)
fn sqrt_branches(
    session: &FieldSession,
    z: &AlgExtCPolyCoeff,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let inner = z.as_inner();
    if !inner.im.iter().all(|e| e.is_zero()) {
        return Err(EvalError::NotImplemented("sqrt branches complex"));
    }
    let beta = session.sqrt_in_field(z)?;
    let neg = session.neg(&beta)?;
    Ok(vec![beta, neg])
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

// **Pipeline private** — verify root vanishes mod minpoly (same prepare path as roots)
fn verify_root(p: &PolyAlgExt, var: &Var, root: &AlgExtCPolyCoeff) {
    let (input, mut session) = PolyInK::prepare_with_session(p, var).expect("prepare");
    session.bump_to(&root.as_inner().field);
    let mut val = session.zero();
    for (m, c) in &input.poly.terms {
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
    assert!(
        val.as_inner()
            .eq_mod(session.zero().as_inner())
            .unwrap_or(false),
        "root does not vanish (dim={})",
        root.as_inner().field.dimension()
    );
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
        let rs = sqrt_branches(&session, &disc).unwrap();
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
        let rs = quadratic_roots_formula(&session, &p, &Var::from("x")).unwrap();
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
        let rs = resolvent_cubic_roots_in_session(&session, &res, &Var::from("_z")).unwrap();
        assert!(
            session.working().dimension() <= 6,
            "resolvent splitting field dim bound"
        );
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

    // **B** — F1: quadratic_roots_formula must not break subsequent adjoin_sqrt (ε²=u).
    // FIX-SQRT-RESOLVENT
    #[test]
    fn adjoin_sqrt_after_quadratic_roots_formula() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let z_var = Var::from("_z");
        let z0 = one_cubic_root(&session, &res, &z_var).unwrap();
        let quad = deflate_monic(&session, &res, &z_var, &z0).unwrap();
        quadratic_roots_formula(&session, &quad, &z_var).unwrap();
        let z = session.lift(&z0).unwrap();
        let sa = session.sqrt_in_field(&z).unwrap_or_else(|e| {
            panic!("sqrt_in_field failed after formula: {e:?}, dim={}", session.working().dimension())
        });
        let sq = session.mul(&sa, &sa).unwrap();
        let (sq, z_a) = session.align(&sq, &z).unwrap();
        assert!(
            sq.coeff_sub(&z_a).unwrap().coeff_is_zero(),
            "adjoin_sqrt after quadratic_roots_formula"
        );
    }

    // **B** — F1: sqrt_disc on resolvent deflate quadratic, then adjoin_sqrt on z0.
    // FIX-SQRT-RESOLVENT
    #[test]
    fn adjoin_sqrt_after_sqrt_disc_on_resolvent() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let z_var = Var::from("_z");
        let z0 = one_cubic_root(&session, &res, &z_var).unwrap();
        let quad = deflate_monic(&session, &res, &z_var, &z0).unwrap();
        let b = session.lift(&coeff_at(&quad, &z_var, 1, &session)).unwrap();
        let c = session.lift(&coeff_at(&quad, &z_var, 0, &session)).unwrap();
        let four = session.int(4).unwrap();
        let b2 = session.mul(&b, &b).unwrap();
        let four_c = session.mul(&four, &c).unwrap();
        let disc = session
            .add(&b2, &session.neg(&four_c).unwrap())
            .unwrap();
        let _sqrt_d = sqrt_disc(&session, &disc).unwrap();
        let z = session.lift(&z0).unwrap();
        let sa = session.adjoin_sqrt(&z).unwrap();
        let sq = session.mul(&sa, &sa).unwrap();
        let (sq, z_a) = session.align(&sq, &z).unwrap();
        assert!(
            sq.coeff_sub(&z_a).unwrap().coeff_is_zero(),
            "adjoin_sqrt after sqrt_disc on resolvent deflate"
        );
    }

    // **B** — F2/M2: t⁴+t+1 resolvent stage dim ≤ 6.
    #[test]
    fn resolvent_dim_bound_t4_plus_t_plus_1() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let _ = resolvent_cubic_roots_in_session(&session, &res, &Var::from("_z")).unwrap();
        assert!(
            session.working().dimension() <= 6,
            "resolvent splitting field for t^4+t+1"
        );
    }

    // **B** — F3: eval at adjoined root of t⁴+t+1 (dim-4 smoke).
    #[test]
    fn eval_vanishes_direct_root_dim4() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)).try_mul(&t).unwrap())
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)))
            .unwrap();
        let dep = monic_univariate(&p, &Var::from("t"), &session).unwrap();
        let minpoly = vec![
            Ratio::from_integer(1.into()),
            Ratio::from_integer(0.into()),
            Ratio::from_integer(0.into()),
            Ratio::from_integer(1.into()),
            Ratio::from_integer(1.into()),
        ];
        let k = ExtensionField::adjoin_irreducible_over_q(minpoly).unwrap();
        session.bump_to(&k);
        let data = AlgExtData::from_field_coords(
            Arc::clone(&k),
            coords_to_expr(&k.generator_coords()).unwrap(),
        )
        .unwrap();
        let root = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&data).unwrap());
        assert!(eval_vanishes(&session, &dep, &Var::from("t"), &root));
    }

    // **B** — F3: γ relation √α√β√γ + q = 0 on t⁴+t+1 resolvent.
    #[test]
    fn euler_gamma_relation_holds() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let z_roots = resolvent_cubic_roots_in_session(&session, &res, &Var::from("_z")).unwrap();
        let sorted = canonical_sort_roots(&session, z_roots).unwrap();
        let q = session.int(1).unwrap();
        let alpha = session.lift(&sorted[0]).unwrap();
        let beta = session.lift(&sorted[1]).unwrap();
        let gamma = session.lift(&sorted[2]).unwrap();
        let alpha_beta = session.add(&alpha, &beta).unwrap();
        let sum_roots = session.add(&alpha_beta, &gamma).unwrap();
        assert!(sum_roots.coeff_is_zero(), "alpha+beta+gamma");
        let ab = session.mul(&alpha, &beta).unwrap();
        let ag = session.mul(&alpha, &gamma).unwrap();
        let bg = session.mul(&beta, &gamma).unwrap();
        let pair_sum = session.add(&session.add(&ab, &ag).unwrap(), &bg).unwrap();
        let minus_four = session.int(-4).unwrap();
        let (pair_sum, minus_four) = session.align(&pair_sum, &minus_four).unwrap();
        assert!(
            pair_sum
                .coeff_sub(&minus_four)
                .unwrap()
                .coeff_is_zero(),
            "alpha beta + alpha gamma + beta gamma"
        );
        let product = session.mul(&ab, &gamma).unwrap();
        let one = session.one();
        let (product, one) = session.align(&product, &one).unwrap();
        assert!(product.coeff_sub(&one).unwrap().coeff_is_zero(), "alpha beta gamma");
        let sa = session.sqrt_principal(&alpha).unwrap();
        let sb = session.sqrt_principal(&beta).unwrap();
        let sg = euler_derived_sqrt_gamma(&session, &q, &sa, &sb).unwrap();
        let prod = session.mul(&session.mul(&sa, &sb).unwrap(), &sg).unwrap();
        let sum = session.add(&prod, &q).unwrap();
        assert!(sum.coeff_is_zero(), "sqrt_alpha*sqrt_beta*sqrt_gamma + q");
        let sg2 = session.mul(&sg, &sg).unwrap();
        let (sg2, gamma) = session.align(&sg2, &gamma).unwrap();
        assert!(sg2.coeff_sub(&gamma).unwrap().coeff_is_zero(), "sg^2 = gamma");
    }

    // **B** — PR-E′: Euler four roots vanish on depressed t⁴+t+1.
    // FIX-T4P1 · Lean: Giac.Examples.T4PlusTPlus1.t4_plus_t_plus_1_euler_vanishes (partial)
    #[test]
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
        let z_roots = resolvent_cubic_roots_in_session(&session, &res, &Var::from("_z")).unwrap();
        let rs = euler_depressed_quartic_roots(&session, &dep, &Var::from("t"), &z_roots).unwrap();
        assert_eq!(rs.len(), 4);
        for r in &rs {
            verify_root(&dep, &Var::from("t"), r);
        }
    }

    #[test]
    fn field_session_dimension_bound_quartic() {
        let session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)).try_mul(&t).unwrap())
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)))
            .unwrap();
        let rs = poly_algext_roots_in_session(&p, &Var::from("t"), &session).unwrap();
        assert_eq!(rs.len(), 4);
        let dim = session.working().dimension();
        assert!(
            dim <= 24,
            "t^4+t+1 splitting field [K:Q] <= 24 (S4), got {dim}"
        );
        for r in &rs {
            verify_root(&p, &Var::from("t"), r);
        }
    }

    // **B** — F4: t⁴+t+1 baseline d_L ≤ D_HARD (24); D_QUARTIC_OUT=12 is A₄ reference only.
    #[test]
    fn field_session_dimension_bound_quartic_tight() {
        let session = q_session();
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)).try_mul(&t).unwrap())
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)))
            .unwrap();
        let rs = poly_algext_roots_in_session(&p, &Var::from("t"), &session).unwrap();
        assert_eq!(rs.len(), 4);
        let dim = session.working().dimension();
        assert!(
            dim <= POLY_ROOTS_DIM_HARD,
            "t^4+t+1 splitting field [K:Q] <= {POLY_ROOTS_DIM_HARD}, got {dim}"
        );
        // ponytail: informational only — not asserted; see algorithm spec §5.4.
        let _ = dim <= POLY_ROOTS_DIM_QUARTIC_OUT;
        for r in &rs {
            verify_root(&p, &Var::from("t"), r);
        }
    }

    // FIX-T4P1 · Lean: Giac.Examples.T4PlusTPlus1.t4_plus_t_plus_1_four_roots (partial)
    #[test]
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

    // **B** — P3-6-C1: general cubic x³−x+1 (S0 solve 暴露).
    #[test]
    fn roots_x3_minus_x_plus_1_vanish() {
        let mut session = q_session();
        let t = PolyAlgExt::ring_var("x");
        let p = t
            .try_pow(3)
            .unwrap()
            .try_sub(&t)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)))
            .unwrap();
        let rs = poly_algext_roots(&p, &Var::from("x")).unwrap();
        assert_eq!(rs.len(), 3);
        for r in &rs {
            verify_root(&p, &Var::from("x"), r);
        }
    }

    // **B** — P3-6-C3: S3 splitting field dim bound for x³−x+1.
    #[test]
    fn field_session_dimension_bound_cubic_x3_minus_x_plus_1() {
        let session = q_session();
        let t = PolyAlgExt::ring_var("x");
        let p = t
            .try_pow(3)
            .unwrap()
            .try_sub(&t)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(&session, 1)))
            .unwrap();
        let rs = poly_algext_roots_in_session(&p, &Var::from("x"), &session).unwrap();
        assert_eq!(rs.len(), 3);
        let dim = session.working().dimension();
        assert!(
            dim <= 6,
            "x^3-x+1 splitting field [K:Q] <= 6 (S3), got {dim}"
        );
    }

    // **B** — F2: resolvent stage dim ≤ 6 for t⁴+t+1 (d_K=1).
    #[test]
    fn f2_resolvent_split_no_cardano_stack() {
        let mut session = q_session();
        let p2 = session.zero();
        let p1 = session.int(1).unwrap();
        let p0 = session.int(1).unwrap();
        let res = build_resolvent_cubic(&session, &p2, &p1, &p0).unwrap();
        let rs = resolvent_cubic_roots_in_session(&session, &res, &Var::from("_z")).unwrap();
        assert_eq!(rs.len(), 3);
        assert!(
            session.working().dimension() <= POLY_ROOTS_DIM_RESOLVENT,
            "F2 resolvent dim bound"
        );
        for r in &rs {
            verify_root(&res, &Var::from("_z"), r);
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
