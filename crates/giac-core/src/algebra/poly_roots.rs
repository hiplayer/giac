//! Univariate exact roots of `Poly<AlgExtCPolyCoeff>` (P2-1/6, P3-6).
//!
//! ```text
//! infer K from coefficients → monic → deg dispatch
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
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;

use super::alg_ext::{algext_cube_root, algext_square_roots, AlgExtData};
use super::alg_ext_c::AlgExtCData;
use super::field_arith::coords_to_expr;
use super::ext_tower::ExtensionField;
use super::poly::PolyAlgExt;
use super::poly_alg_coeff::AlgExtCPolyCoeff;

/// Exact roots of univariate `p` w.r.t. `var` over the coefficient field of `p`.
/// **Stable (bounded)** — exact AlgExtC roots deg 1–4; quartic resolvent gap
pub fn poly_algext_roots(p: &PolyAlgExt, var: &Var) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let field = infer_field(p)?;
    let p = normalize_coeffs(p, &field)?;
    let monic = monic_univariate(&p, var, &field)?;
    let d = monic.degree_wrt(var);
    match d {
        0 => {
            if monic.is_zero() {
                Ok(vec![])
            } else {
                Err(EvalError::TypeError("constant has no roots"))
            }
        }
        1 => linear_root(&monic, var, &field),
        2 => quadratic_roots(&monic, var, &field),
        3 => cubic_roots(&monic, var, &field),
        4 => quartic_roots(&monic, var, &field),
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
fn coeff_at(p: &PolyAlgExt, var: &Var, exp: u64, field: &Arc<ExtensionField>) -> AlgExtCPolyCoeff {
    for (m, c) in &p.terms {
        if exp == 0 && m.is_const() {
            return c.clone();
        }
        if m.exp_of(var) == exp && m.iter().all(|(v, _)| v == var) {
            return c.clone();
        }
    }
    ring_zero(field)
}

// **Pipeline private** — `ring_zero`
fn ring_zero(field: &Arc<ExtensionField>) -> AlgExtCPolyCoeff {
    AlgExtCPolyCoeff::from(AlgExtCData::zero(Arc::clone(field)).expect("zero"))
}

// **Pipeline private** — `ring_one`
fn ring_one(field: &Arc<ExtensionField>) -> AlgExtCPolyCoeff {
    AlgExtCPolyCoeff::from(AlgExtCData::one(Arc::clone(field)).expect("one"))
}

// **Pipeline private** — `rat_in_field`
fn rat_in_field(field: &Arc<ExtensionField>, r: Ratio<BigInt>) -> Result<AlgExtCPolyCoeff, EvalError> {
    let coords = field.embed_rational(&r);
    let a = AlgExtData::from_field_coords(Arc::clone(field), coords_to_expr(&coords)?)?;
    Ok(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&a)?))
}

// **Pipeline private** — `ring_int`
fn ring_int(field: &Arc<ExtensionField>, n: i64) -> Result<AlgExtCPolyCoeff, EvalError> {
    rat_in_field(field, Ratio::from_integer(BigInt::from(n)))
}

// **Pipeline private** — `ring_half`
fn ring_half(field: &Arc<ExtensionField>) -> Result<AlgExtCPolyCoeff, EvalError> {
    rat_in_field(field, Ratio::new(1.into(), 2.into()))
}

// **Pipeline private** — lift all coeffs to ambient K
fn normalize_coeffs(p: &PolyAlgExt, field: &Arc<ExtensionField>) -> Result<PolyAlgExt, EvalError> {
    let mut out = PolyAlgExt::ring_zero();
    for (m, c) in &p.terms {
        let c = lift_to_field(c, field)?;
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
    field: &Arc<ExtensionField>,
) -> Result<PolyAlgExt, EvalError> {
    let deg = p.degree_wrt(var);
    let lc = coeff_at(p, var, deg, field);
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
        let c = coeff_at(p, var, exp, field);
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
    let _ = field;
    Ok(out)
}

// **Pipeline private** — `linear_root`
fn linear_root(
    p: &PolyAlgExt,
    var: &Var,
    field: &Arc<ExtensionField>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let a = coeff_at(p, var, 1, field);
    let b = coeff_at(p, var, 0, field);
    if a.coeff_is_zero() {
        return Err(EvalError::TypeError("not linear"));
    }
    Ok(vec![b.coeff_neg()?.coeff_div(&a)?])
}

// **Pipeline private** — embed coeff into target ExtensionField
fn lift_to_field(
    c: &AlgExtCPolyCoeff,
    field: &Arc<ExtensionField>,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let (a, _) = align_coeff(c, &ring_zero(field))?;
    Ok(a)
}

// **Pipeline private** — `quadratic_roots`
fn quadratic_roots(
    p: &PolyAlgExt,
    var: &Var,
    field: &Arc<ExtensionField>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let b = coeff_at(p, var, 1, field);
    let c = coeff_at(p, var, 0, field);
    let four = ring_int(field, 4)?;
    let disc = b.coeff_mul(&b)?.coeff_sub(&c.coeff_mul(&four)?)?;
    if disc.coeff_is_zero() {
        let ext = Arc::clone(&b.as_inner().field);
        let two = ring_int(&ext, 2)?;
        let nb = lift_to_field(&b, &ext)?.coeff_neg()?;
        return Ok(vec![nb.coeff_div(&two)?]);
    }
    let sqrt_d = sqrt_disc(&disc)?;
    let ext = Arc::clone(&sqrt_d.as_inner().field);
    let mut nb = lift_to_field(&b, &ext)?.coeff_neg()?;
    let mut sqrt_d = sqrt_d;
    let mut two = ring_int(&ext, 2)?;
    let (nb_a, sqrt_a) = align_coeff(&nb, &sqrt_d)?;
    let (two_a, nb_a) = align_coeff(&two, &nb_a)?;
    nb = nb_a;
    sqrt_d = sqrt_a;
    two = two_a;
    let (two_b, sqrt_b) = align_coeff(&two, &sqrt_d)?;
    two = two_b;
    Ok(vec![
        nb.clone().coeff_add(&sqrt_b)?.coeff_div(&two)?,
        nb.coeff_sub(&sqrt_b)?.coeff_div(&two)?,
    ])
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
    p: &PolyAlgExt,
    var: &Var,
    field: &Arc<ExtensionField>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let a2 = coeff_at(p, var, 2, field);
    let a1 = coeff_at(p, var, 1, field);
    if a2.coeff_is_zero() && a1.coeff_is_zero() {
        let a0 = coeff_at(p, var, 0, field);
        let r = algext_c_cube_root(&a0.coeff_neg()?)?;
        return Ok(vec![r]);
    }
    let r0 = one_cubic_root(p, var, field)?;
    let mut roots = vec![r0.clone()];
    let quad = deflate_monic(p, var, &r0)?;
    let ext_field = Arc::clone(&r0.as_inner().field);
    roots.extend(quadratic_roots(&quad, var, &ext_field)?);
    Ok(roots)
}

// **Pipeline private** — `one_cubic_root`
fn one_cubic_root(
    p: &PolyAlgExt,
    var: &Var,
    field: &Arc<ExtensionField>,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let a2 = coeff_at(p, var, 2, field);
    let a1 = coeff_at(p, var, 1, field);
    let a0 = coeff_at(p, var, 0, field);
    let three = ring_int(field, 3)?;
    let shift = a2.coeff_div(&three)?.coeff_neg()?;
    let p_dep = a1
        .coeff_sub(&a2.coeff_mul(&shift)?)?
        .coeff_div(&three)?;
    let q_dep = a0
        .coeff_sub(&a1.coeff_mul(&shift)?)?
        .coeff_add(&a2.coeff_mul(&shift)?.coeff_mul(&shift)?)?
        .coeff_add(&shift.coeff_mul(&shift)?.coeff_mul(&shift)?)?;
    if p_dep.coeff_is_zero() {
        let r = algext_c_cube_root(&q_dep.coeff_neg()?)?;
        return r.coeff_add(&shift);
    }
    let half = ring_half(field)?;
    let twenty_seven = ring_int(field, 27)?;
    let delta = q_dep
        .coeff_mul(&half)?
        .coeff_mul(&half)?
        .coeff_add(&p_dep.coeff_mul(&p_dep)?.coeff_mul(&p_dep)?.coeff_div(&twenty_seven)?)?;
    let sqrt_delta = algext_c_sqrt(&delta)?
        .into_iter()
        .next()
        .ok_or(EvalError::NotImplemented("algext sqrt"))?;
    let u = algext_c_cube_root(
        &q_dep.coeff_mul(&half)?.coeff_neg()?.coeff_add(&sqrt_delta)?,
    )?;
    let v = algext_c_cube_root(
        &q_dep.coeff_mul(&half)?.coeff_neg()?.coeff_sub(&sqrt_delta)?,
    )?;
    let (u, v) = align_coeff(&u, &v)?;
    u.coeff_add(&v)?.coeff_add(&shift)
}

// **Pipeline private** — `quartic_roots`
fn quartic_roots(
    p: &PolyAlgExt,
    var: &Var,
    field: &Arc<ExtensionField>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let a3 = coeff_at(p, var, 3, field);
    let a1 = coeff_at(p, var, 1, field);
    if a3.coeff_is_zero() && a1.coeff_is_zero() {
        return biquadratic_roots(p, var, field);
    }
    let four = ring_int(field, 4)?;
    let shift = a3.coeff_div(&four)?.coeff_neg()?;
    let dep = depress_quartic(p, var, &shift, field)?;
    let p2 = coeff_at(&dep, var, 2, field);
    let p1 = coeff_at(&dep, var, 1, field);
    let p0 = coeff_at(&dep, var, 0, field);
    let res = build_resolvent_cubic(&p2, &p1, &p0, field)?;
    let z_roots = cubic_roots(&res, &Var::from("_z"), field)?;
    let mut all = Vec::new();
    for z in z_roots {
        let mut rs = split_depressed_quartic(&dep, var, &z, field)?;
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
    p: &PolyAlgExt,
    var: &Var,
    field: &Arc<ExtensionField>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let b = coeff_at(p, var, 2, field);
    let c = coeff_at(p, var, 0, field);
    let u = PolyAlgExt::ring_var(Var::from("u"))
        .try_pow(2)?
        .try_add(&PolyAlgExt::ring_constant(b))?
        .try_add(&PolyAlgExt::ring_constant(c))?;
    let u_roots = quadratic_roots(&u, &Var::from("u"), field)?;
    let mut out = Vec::new();
    for ur in u_roots {
        out.extend(algext_c_sqrt(&ur)?);
    }
    let _ = var;
    Ok(out)
}

// **Pipeline private** — `depress_quartic`
fn depress_quartic(
    p: &PolyAlgExt,
    var: &Var,
    shift: &AlgExtCPolyCoeff,
    field: &Arc<ExtensionField>,
) -> Result<PolyAlgExt, EvalError> {
    let x = PolyAlgExt::ring_var(var.clone());
    let t = x.try_add(&PolyAlgExt::ring_constant(shift.clone()))?;
    let mut sum = PolyAlgExt::ring_zero();
    for e in 0..=4 {
        let c = coeff_at(p, var, e, field);
        if c.coeff_is_zero() {
            continue;
        }
        sum = sum.try_add(&PolyAlgExt::ring_constant(c).try_mul(&t.try_pow(e)?)?)?;
    }
    monic_univariate(&sum, var, field)
}

// **Pipeline private** — `build_resolvent_cubic`
fn build_resolvent_cubic(
    p: &AlgExtCPolyCoeff,
    q: &AlgExtCPolyCoeff,
    r: &AlgExtCPolyCoeff,
    field: &Arc<ExtensionField>,
) -> Result<PolyAlgExt, EvalError> {
    let four = ring_int(field, 4)?;
    let z = PolyAlgExt::ring_var(Var::from("_z"));
    Ok(z.try_pow(3)?
        .try_sub(&PolyAlgExt::ring_constant(p.clone()).try_mul(&z.try_pow(2)?)?)?
        .try_sub(&PolyAlgExt::ring_constant(r.coeff_mul(&four)?.coeff_neg()?).try_mul(&z)?)?
        .try_add(&PolyAlgExt::ring_constant(
            r.coeff_mul(p)?.coeff_mul(&four)?.coeff_sub(&q.coeff_mul(q)?)?,
        ))?)
}

// **Pipeline private** — `split_depressed_quartic`
fn split_depressed_quartic(
    dep: &PolyAlgExt,
    var: &Var,
    z: &AlgExtCPolyCoeff,
    field: &Arc<ExtensionField>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError> {
    let q = coeff_at(dep, var, 1, field);
    let r = coeff_at(dep, var, 0, field);
    let four = ring_int(field, 4)?;
    let two = ring_int(field, 2)?;
    let m = algext_c_sqrt(&z.coeff_mul(z)?.coeff_sub(&r.coeff_mul(&four)?)?)?
        .into_iter()
        .next()
        .ok_or(EvalError::NotImplemented("algext sqrt"))?;
    let half = ring_half(field)?;
    let mut out = Vec::new();
    for mp in [z.coeff_add(&m.clone())?, z.coeff_sub(&m)?] {
        let c0 = mp.coeff_mul(&half)?;
        let c1 = if q.coeff_is_zero() {
            ring_zero(field)
        } else {
            q.coeff_div(&m)?.coeff_mul(&half)?
        };
        let quad = PolyAlgExt::ring_var(var.clone())
            .try_pow(2)?
            .try_add(&PolyAlgExt::ring_constant(c1).try_mul(&PolyAlgExt::ring_var(var.clone()))?)?
            .try_add(&PolyAlgExt::ring_constant(c0))?;
        out.extend(quadratic_roots(&quad, var, field)?);
    }
    Ok(out)
}

// **Pipeline private** — align two AlgExtCPolyCoeff to common field; retire FieldSession
fn align_coeff(a: &AlgExtCPolyCoeff, b: &AlgExtCPolyCoeff) -> Result<(AlgExtCPolyCoeff, AlgExtCPolyCoeff), EvalError> {
    let (aa, bb) = AlgExtCData::align_pair(a.as_inner(), b.as_inner())?;
    Ok((AlgExtCPolyCoeff::from(aa), AlgExtCPolyCoeff::from(bb)))
}

// **Pipeline private** — `deflate_monic`
fn deflate_monic(
    p: &PolyAlgExt,
    var: &Var,
    root: &AlgExtCPolyCoeff,
) -> Result<PolyAlgExt, EvalError> {
    let n = p.degree_wrt(var) as usize;
    let field = Arc::clone(&root.as_inner().field);
    let mut qs = vec![ring_zero(&field); n];
    qs[n - 1] = ring_one(&field);
    for k in (0..n - 1).rev() {
        let ak = coeff_at(p, var, k as u64, &field);
        qs[k] = ak.coeff_add(&root.coeff_mul(&qs[k + 1])?)?;
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
    if inner.im.iter().all(|e| matches!(e.as_ref(), crate::Expr::Int(n) if n.is_zero())) {
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
    if inner.im.iter().all(|e| matches!(e.as_ref(), crate::Expr::Int(n) if n.is_zero())) {
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

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_poly::Var;
    use num_rational::Ratio;
    use serial_test::serial;

    use crate::algebra::test_fixtures::sqrt2_algext;
    use crate::expr::Expr;

    use super::*;
    use super::super::alg_ext_c::AlgExtCData;
    use super::super::poly::poly_alg_from_expr;

    fn q_field() -> Arc<ExtensionField> {
        ExtensionField::rational()
    }

    fn rat_coeff(n: i64) -> AlgExtCPolyCoeff {
        rat_in_field(&q_field(), Ratio::from_integer(BigInt::from(n))).unwrap()
    }

    fn verify_root(p: &PolyAlgExt, var: &Var, root: &AlgExtCPolyCoeff) {
        let poly_field = infer_field(p).expect("infer field");
        let p = monic_univariate(
            &normalize_coeffs(p, &poly_field).expect("normalize"),
            var,
            &poly_field,
        )
        .expect("monic");
        let eval_field = Arc::clone(&root.as_inner().field);
        let mut val = ring_zero(&eval_field);
        for (m, c) in &p.terms {
            let exp = m.exp_of(var);
            let mut pow = ring_one(&eval_field);
            for _ in 0..exp {
                pow = pow.coeff_mul(root).unwrap();
            }
            let c = lift_to_field(c, &eval_field).unwrap();
            let (c, pow) = align_coeff(&c, &pow).unwrap();
            val = val.coeff_add(&c.coeff_mul(&pow).unwrap()).unwrap();
        }
        assert!(val.coeff_is_zero(), "root does not vanish");
    }

    #[serial]
    #[test]
    fn cubic_one_root_vanishes() {
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(3)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(2)))
            .unwrap();
        let field = q_field();
        let r = one_cubic_root(&p, &Var::from("t"), &field).unwrap();
        verify_root(&p, &Var::from("t"), &r);
    }

    #[serial]
    #[test]
    fn quadratic_sqrt_four_times_sqrt2_over_k1() {
        let sqrt2 = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap());
        let four = ring_int(&Arc::clone(&sqrt2.as_inner().field), 4).unwrap();
        let disc = sqrt2.coeff_mul(&four).unwrap();
        let rs = algext_c_sqrt(&disc).unwrap();
        assert_eq!(rs.len(), 2);
        for beta in &rs {
            let sq = beta.coeff_mul(beta).unwrap();
            let (sq_a, d_a) = align_coeff(&sq, &disc).unwrap();
            assert!(sq_a.coeff_sub(&d_a).unwrap().coeff_is_zero());
            let field = Arc::clone(&beta.as_inner().field);
            let two = ring_int(&field, 2).unwrap();
            let half = ring_half(&field).unwrap();
            let root_mul = beta.coeff_mul(&half).unwrap();
            let root_div = beta.coeff_div(&two).unwrap();
            let (rm, rd) = align_coeff(&root_mul, &root_div).unwrap();
            assert!(
                rm.coeff_sub(&rd).unwrap().coeff_is_zero(),
                "mul/2 and div/2 should agree"
            );
            let sqrt2 = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap());
            let root_sq = root_div.coeff_mul(&root_div).unwrap();
            let (rsq, s2) = align_coeff(&root_sq, &sqrt2).unwrap();
            assert!(
                rsq.coeff_sub(&s2).unwrap().coeff_is_zero(),
                "root^2 should equal sqrt2"
            );
        }
    }

    #[serial]
    #[test]
    fn quadratic_x2_minus_sqrt2_roots_vanish() {
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(-1), sqrt2_algext().into_expr()]),
        ]);
        let p = poly_alg_from_expr(&e).unwrap();
        let field = infer_field(&p).unwrap();
        let rs = quadratic_roots(&p, &Var::from("x"), &field).unwrap();
        assert_eq!(rs.len(), 2);
        for r in &rs {
            verify_root(&p, &Var::from("x"), r);
        }
    }

    #[serial]
    #[test]
    fn roots_quadratic_x2_minus_2() {
        let x = PolyAlgExt::ring_var("x");
        let p = x
            .try_pow(2)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(2)))
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
            let (sq, s2) = align_coeff(&sq, &sqrt2).unwrap();
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
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(3)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(2)))
            .unwrap();
        let rs = poly_algext_roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 1);
        verify_root(&p, &Var::from("t"), &rs[0]);
    }

    #[serial]
    #[test]
    #[ignore = "quartic resolvent Cardano exceeds 10s; needs common-field optimization"]
    fn roots_quartic_t4_plus_t_plus_1() {
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(1)).try_mul(&t).unwrap())
            .unwrap()
            .try_add(&PolyAlgExt::ring_constant(rat_coeff(1)))
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
        let t = PolyAlgExt::ring_var("t");
        let p = t
            .try_pow(4)
            .unwrap()
            .try_sub(&PolyAlgExt::ring_constant(rat_coeff(2)))
            .unwrap();
        let rs = poly_algext_roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 4);
        for r in &rs {
            verify_root(&p, &Var::from("t"), r);
        }
    }
}
