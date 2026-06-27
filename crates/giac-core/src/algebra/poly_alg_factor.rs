//! Univariate factorization over K for `Poly<AlgExtC>` (T2-1 / L3-1).
//!
//! **Upstream:** `ext_factor` / `ext_factor_nodegck` — sqff → linear / quadratic split / roots / witness.

use std::sync::Arc;

use giac_poly::{factor_into, FlatUni, MainVar, Poly, PolyCoeff, Var};
use num_traits::Zero;

use crate::error::EvalError;

use super::ext_tower::ExtensionField;
use super::field_arith::{expr_to_ratio, poly_degree, CoordsQ};
use super::field_session::{coeff_from_coords, FieldSession};
use super::poly::{poly_algext_from_poly, PolyAlgExt};
use super::poly_alg_coeff::AlgExtCPolyCoeff;
use super::poly_alg_ops::{
    div_rem_wrt_algext, eval_univariate_algext_at, infer_ambient_field, normalize_algext_poly,
    split_quadratic_factor,
};
use super::poly_roots::poly_algext_roots;

/// Formal layer variable for minpoly factorization over **K**.
pub(crate) fn common_minpoly_var() -> Var {
    Var::from("__u")
}

/// **Stable (bounded)** — irreducible factors over K for univariate `p`; `None` if not split.
pub fn factor_into_algext(p: &PolyAlgExt) -> Result<Option<Vec<PolyAlgExt>>, EvalError> {
    let Some(var) = univariate_main_var(p) else {
        return Ok(None);
    };
    let field = infer_ambient_field(p)?;
    let session = FieldSession::new(field);
    let flat = FlatUni::try_new(
        normalize_algext_poly(p, &session)?,
        MainVar::new(var.clone()),
    )
    .map_err(Into::into)?;
    let pairs = factor_univariate_over_k(&session, &flat)?;
    let factors = expand_factor_pairs(pairs);
    if factors.len() <= 1 {
        return Ok(None);
    }
    Ok(Some(factors))
}

// **Pipeline private** — single main variable when `p` is univariate in it.
fn univariate_main_var(p: &PolyAlgExt) -> Option<Var> {
    use std::collections::BTreeSet;

    let mut vars = BTreeSet::new();
    for m in p.terms.keys() {
        for (v, _) in m.iter() {
            vars.insert(v.clone());
        }
    }
    match vars.len() {
        0 => None,
        1 => vars.into_iter().next(),
        _ => {
            let v = vars.iter().next()?.clone();
            if giac_poly::is_univariate_in(p, &v) {
                Some(v)
            } else {
                None
            }
        }
    }
}

fn expand_factor_pairs(pairs: Vec<(PolyAlgExt, usize)>) -> Vec<PolyAlgExt> {
    let mut out = Vec::new();
    for (f, k) in pairs {
        for _ in 0..k {
            out.push(f.clone());
        }
    }
    out
}

/// **Stable** — factor `Poly` over ℚ via K[var] pipeline when univariate split exists.
pub fn factor_into_via_algext(p: &Poly) -> Result<Option<Vec<PolyAlgExt>>, EvalError> {
    factor_into_algext(&super::poly::poly_algext_from_poly(p)?)
}

/// **Stable (bounded)** — sqff factor pairs over K for univariate `Poly` over ℚ (T2-3).
pub fn factor_univariate_pairs_over_k(
    p: &Poly,
    var: &Var,
) -> Result<Vec<(PolyAlgExt, usize)>, EvalError> {
    let p_alg = super::poly::poly_algext_from_poly(p)?;
    let field = infer_ambient_field(&p_alg)?;
    let session = FieldSession::new(field);
    let flat = FlatUni::try_new(
        normalize_algext_poly(&p_alg, &session)?,
        MainVar::new(var.clone()),
    )
    .map_err(Into::into)?;
    factor_univariate_over_k(&session, &flat)
}

/// **Stable** — factor `flat` in K[var] as `(factor, multiplicity)` pairs.
pub fn factor_univariate_over_k(
    session: &FieldSession,
    flat: &FlatUni<AlgExtCPolyCoeff>,
) -> Result<Vec<(PolyAlgExt, usize)>, EvalError> {
    let _var = flat.var().as_var().clone();
    let flat = aligned_flat(session, flat)?;
    let mut out = Vec::new();
    let c = flat.content()?;
    if !c.coeff_is_one() && !c.coeff_is_zero() {
        push_factor(&mut out, PolyAlgExt::ring_constant(c), 1);
    }
    let pp = flat.primitive_part()?;
    if pp.as_poly().is_one() {
        return Ok(out);
    }
    let deg = pp.degree();
    if deg <= 1 {
        push_factor(&mut out, pp.into_poly(), 1);
        return Ok(out);
    }
    for (g, k) in pp.square_free_factorization()? {
        for f in factor_square_free_over_k(session, &g)? {
            push_factor(&mut out, f, k);
        }
    }
    Ok(out)
}

/// **Stable** — flat list of factors (with repetition for multiplicity).
pub fn factor_univariate_flat_over_k(
    session: &FieldSession,
    flat: &FlatUni<AlgExtCPolyCoeff>,
) -> Result<Vec<PolyAlgExt>, EvalError> {
    let mut out = Vec::new();
    for (f, k) in factor_univariate_over_k(session, flat)? {
        for _ in 0..k {
            out.push(f.clone());
        }
    }
    Ok(out)
}

// **Pipeline private** — normalize and re-wrap as [`FlatUni`].
fn aligned_flat(
    session: &FieldSession,
    flat: &FlatUni<AlgExtCPolyCoeff>,
) -> Result<FlatUni<AlgExtCPolyCoeff>, EvalError> {
    FlatUni::try_new(
        normalize_algext_poly(flat.as_poly(), session)?,
        flat.var().clone(),
    )
    .map_err(Into::into)
}

// **Pipeline private** — embed `FlatUni` as ℚ `Poly` when all coeffs lie in ℚ.
fn try_as_rational_poly(g: &FlatUni<AlgExtCPolyCoeff>) -> Option<Poly> {
    let p = g.as_poly();
    let mut out = Poly::ring_zero();
    for (m, c) in &p.terms {
        let alg_c = c.as_inner();
        if !alg_c.im.iter().all(|e| e.is_zero()) || !alg_c.field.is_base() {
            return None;
        }
        let r = expr_to_ratio(alg_c.re.first()?.as_ref()).ok()?;
        out = out.add(&Poly::term(m.clone(), r));
    }
    Some(out)
}

// **Pipeline private** — square-free factor in K[var].
fn factor_square_free_over_k(
    session: &FieldSession,
    g: &FlatUni<AlgExtCPolyCoeff>,
) -> Result<Vec<PolyAlgExt>, EvalError> {
    let var = g.var().as_var();
    let p = g.as_poly();
    let deg = g.degree();
    if deg <= 1 {
        return Ok(vec![p.clone()]);
    }
    if deg == 2 {
        let lin = split_quadratic_factor(session, p, var)?;
        if lin.len() == 2 {
            return Ok(lin);
        }
        return Ok(vec![p.clone()]);
    }
    if deg >= 3 {
        if let Some(qpoly) = try_as_rational_poly(g) {
            if let Some(qfactors) = factor_into(&qpoly) {
                if qfactors.len() > 1 {
                    let mut out = Vec::new();
                    for qf in qfactors {
                        let lifted = poly_algext_from_poly(&qf)?;
                        let flat = FlatUni::try_new(lifted, MainVar::new(var.clone()))
                            .map_err(Into::into)?;
                        out.extend(factor_square_free_over_k(session, &flat)?);
                    }
                    return Ok(out);
                }
            }
        }
    }
    if deg == 4 && is_even_only(p, var) {
        if let Some(facs) = try_factor_via_x_squared(session, p, var)? {
            return Ok(facs);
        }
    }
    if deg == 3 || deg == 4 {
        if let Some(facs) = try_factor_by_roots(session, p, var)? {
            return Ok(facs);
        }
    }
    Ok(vec![p.clone()])
}

// **Pipeline private** — full linear split when all roots lie in K.
fn try_factor_by_roots(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Option<Vec<PolyAlgExt>>, EvalError> {
    let deg = p.degree_wrt(var);
    let roots = poly_algext_roots(p, var)?;
    if roots.len() != deg as usize {
        return Ok(None);
    }
    let mut out = Vec::new();
    for r in roots {
        let r = session.lift(&r)?;
        out.push(
            PolyAlgExt::ring_var(var.clone())
                .try_sub(&PolyAlgExt::ring_constant(r))?,
        );
    }
    if !product_divides(session, p, &out, var)? {
        return Ok(None);
    }
    Ok(Some(out))
}

// **Pipeline private** — p(x)=h(x²); factor h quadratically and lift to x²−r.
fn try_factor_via_x_squared(
    session: &FieldSession,
    p: &PolyAlgExt,
    var: &Var,
) -> Result<Option<Vec<PolyAlgExt>>, EvalError> {
    let h = halve_exponents(p, var);
    if h.degree_wrt(var) != 2 {
        return Ok(None);
    }
    let lin = split_quadratic_factor(session, &h, var)?;
    if lin.len() != 2 {
        return Ok(None);
    }
    let mut out = Vec::new();
    for f in lin {
        let r = linear_root_coeff(session, &f, var)?;
        let x2 = PolyAlgExt::ring_var(var.clone()).try_pow(2)?;
        out.push(x2.try_sub(&PolyAlgExt::ring_constant(r))?);
    }
    if !product_divides(session, p, &out, var)? {
        return Ok(None);
    }
    Ok(Some(out))
}

fn is_even_only(p: &PolyAlgExt, var: &Var) -> bool {
    p.terms
        .keys()
        .all(|m| m.exp_of(var) % 2 == 0)
}

fn halve_exponents(p: &PolyAlgExt, var: &Var) -> PolyAlgExt {
    let deg = p.degree_wrt(var);
    let mut out = PolyAlgExt::ring_zero();
    for exp in (0..=deg).step_by(2) {
        let c = giac_poly::scalar_coeff_wrt(p, var, exp);
        if c.coeff_is_zero() {
            continue;
        }
        let half = exp / 2;
        let term = if half == 0 {
            PolyAlgExt::ring_constant(c)
        } else {
            PolyAlgExt::ring_constant(c)
                .try_mul(&PolyAlgExt::ring_var(var.clone()).try_pow(half).expect("pow"))
                .expect("mul")
        };
        out = out.try_add(&term).expect("add");
    }
    out
}

fn linear_root_coeff(
    session: &FieldSession,
    lin: &PolyAlgExt,
    var: &Var,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let deg = lin.degree_wrt(var);
    if deg != 1 {
        return Err(EvalError::TypeError("not linear"));
    }
    let lc = giac_poly::scalar_coeff_wrt(lin, var, 1);
    if !lc.coeff_is_one() {
        return Err(EvalError::TypeError("not monic linear"));
    }
    let c0 = giac_poly::scalar_coeff_wrt(lin, var, 0);
    session.neg(&c0)
}

fn product_divides(
    session: &FieldSession,
    p: &PolyAlgExt,
    factors: &[PolyAlgExt],
    var: &Var,
) -> Result<bool, EvalError> {
    let prod = factors
        .iter()
        .try_fold(PolyAlgExt::ring_one(), |acc, f| acc.try_mul(f))?;
    let (_, r) = div_rem_wrt_algext(session, p, &prod, var)?;
    Ok(r.is_zero())
}

// **Pipeline private** — embed monic `poly1` / ℚ into `PolyAlgExt` over `base`.
pub(crate) fn minpoly_q_to_poly_over_field(
    minpoly: &CoordsQ,
    base: &Arc<ExtensionField>,
    var: &Var,
) -> Result<PolyAlgExt, EvalError> {
    let deg = poly_degree(minpoly);
    let u = PolyAlgExt::ring_var(var.clone());
    let mut p = PolyAlgExt::ring_zero();
    for (i, c) in minpoly.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let exp = deg - i;
        let coords = base.embed_rational(c);
        let coeff = coeff_from_coords(base, &coords)?;
        let mut term = PolyAlgExt::ring_constant(coeff);
        if exp > 0 {
            term = term.try_mul(&u.try_pow(exp as u64)?)?;
        }
        p = p.try_add(&term)?;
    }
    Ok(p)
}

// **Pipeline private** — `poly1` minpoly over ℚ → `Poly` in `var`.
fn minpoly_coords_to_poly(minpoly: &CoordsQ, var: &Var) -> Poly {
    let deg = poly_degree(minpoly);
    let u = Poly::var(&**var);
    let mut p = Poly::ring_zero();
    for (i, c) in minpoly.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let exp = deg - i;
        let term = if exp == 0 {
            Poly::constant(c.clone())
        } else {
            Poly::constant(c.clone()).mul(&u.pow(exp as u64))
        };
        p = p.add(&term);
    }
    p
}

/// Upstream `ext_factor` on `b` minpoly over **K** = `base` (U2).
// **Pipeline private** — `factor_minpoly_over_field`
pub(crate) fn factor_minpoly_over_field(
    minpoly: &CoordsQ,
    _base: &Arc<ExtensionField>,
) -> Result<Vec<PolyAlgExt>, EvalError> {
    let var = common_minpoly_var();
    let session = FieldSession::new(Arc::clone(_base));
    let q_poly = minpoly_coords_to_poly(minpoly, &var);
    let q_chunks = match factor_into(&q_poly) {
        Some(facs) if !facs.is_empty() => facs,
        _ => vec![q_poly],
    };
    q_chunks
        .into_iter()
        .map(|qf| {
            let p_alg = poly_algext_from_poly(&qf)?;
            normalize_algext_poly(&p_alg, &session)
        })
        .collect()
}

/// Upstream `common_EXT` factor choice: vanishing at `b` gen, else minimum degree (U2).
// **Pipeline private** — `select_factor_for_common`
pub(crate) fn select_factor_for_common(
    factors: &[PolyAlgExt],
    b: &Arc<ExtensionField>,
    _base: &Arc<ExtensionField>,
) -> Result<PolyAlgExt, EvalError> {
    if factors.is_empty() {
        return Err(EvalError::TypeError("select_factor: empty"));
    }
    if factors.len() == 1 {
        return Ok(factors[0].clone());
    }
    let var = common_minpoly_var();
    let eval_session = FieldSession::new(Arc::clone(_base));
    let can_eval_at_b = ExtensionField::is_subfield_of(b, _base)
        || ExtensionField::is_subfield_of(_base, b)
        || Arc::ptr_eq(b, _base);
    if can_eval_at_b {
        let b_gen = if ExtensionField::is_subfield_of(b, _base) {
            let emb = ExtensionField::try_subfield_embedding(b, _base)?
                .ok_or(EvalError::TypeError("select_factor: embed b"))?;
            coeff_from_coords(_base, &emb.apply(&b.generator_coords()))?
        } else {
            coeff_from_coords(b, &b.generator_coords())?
        };
        for f in factors {
            let val = eval_univariate_algext_at(&eval_session, f, &var, &b_gen)?;
            if val.coeff_is_zero() {
                return Ok(f.clone());
            }
        }
    }
    factors
        .iter()
        .min_by_key(|f| f.degree_wrt(&var))
        .cloned()
        .ok_or(EvalError::TypeError("select_factor: no factor"))
}

fn push_factor(out: &mut Vec<(PolyAlgExt, usize)>, f: PolyAlgExt, k: usize) {
    if f.is_one() {
        return;
    }
    if let Some((last, m)) = out.last_mut() {
        if last == &f {
            *m += k;
            return;
        }
    }
    out.push((f, k));
}

#[cfg(test)]
mod tests {
    

    use num_rational::Ratio;
    use num_traits::Zero;

    use super::*;
    use crate::algebra::ext_tower::ExtensionField;
    use crate::algebra::test_fixtures::k1_adjoin_sqrt2;

    fn x_var() -> Var {
        Var::from("x")
    }

    fn flat(_session: &FieldSession, p: PolyAlgExt) -> FlatUni<AlgExtCPolyCoeff> {
        FlatUni::try_new(p, MainVar::new(x_var())).unwrap()
    }

    fn x_squared_minus_2(session: &FieldSession) -> PolyAlgExt {
        let x = PolyAlgExt::ring_var(x_var());
        let two = session.int(2).unwrap();
        x.try_mul(&x).unwrap().try_sub(&PolyAlgExt::ring_constant(two)).unwrap()
    }

    fn x_fourth_minus_4(session: &FieldSession) -> PolyAlgExt {
        let x2 = PolyAlgExt::ring_var(x_var()).try_pow(2).unwrap();
        let four = session.int(4).unwrap();
        x2.try_mul(&x2).unwrap().try_sub(&PolyAlgExt::ring_constant(four)).unwrap()
    }

    #[test]
    fn factor_into_algext_x_squared_minus_2() {
        let session = FieldSession::new(ExtensionField::rational());
        let p = x_squared_minus_2(&session);
        let factors = factor_into_algext(&p).unwrap().expect("split");
        assert_eq!(factors.len(), 2);
        let prod = factors
            .iter()
            .try_fold(PolyAlgExt::ring_one(), |acc, f| acc.try_mul(f))
            .unwrap();
        let (_, r) = div_rem_wrt_algext(&session, &p, &prod, &x_var()).unwrap();
        assert!(r.is_zero());
    }

    #[test]
    fn factor_into_via_algext_from_rational_poly() {
        use giac_poly::Poly;

        let x = Poly::var("x");
        let p = x.pow(2).sub(&Poly::constant(Ratio::from_integer(2.into())));
        let factors = factor_into_via_algext(&p).unwrap().expect("split");
        assert_eq!(factors.len(), 2);
    }

    #[test]
    fn factor_x_squared_minus_2_over_k1() {
        let k = k1_adjoin_sqrt2();
        let session = FieldSession::new(k);
        let p = x_squared_minus_2(&session);
        let factors = factor_univariate_flat_over_k(&session, &flat(&session, p.clone())).unwrap();
        assert_eq!(factors.len(), 2);
        let prod = factors
            .iter()
            .try_fold(PolyAlgExt::ring_one(), |acc, f| acc.try_mul(f))
            .unwrap();
        let (_, r) = div_rem_wrt_algext(&session, &p, &prod, &x_var()).unwrap();
        assert!(r.is_zero());
    }

    #[test]
    fn factor_irreducible_deg5_witness_over_k() {
        use giac_poly::Poly;

        let x = Poly::var("x");
        let p_rational = x.pow(5).sub(&x).add(&Poly::one());
        let p_alg = super::super::poly::poly_algext_from_poly(&p_rational).unwrap();
        let session = FieldSession::new(ExtensionField::rational());
        let flat = flat(&session, p_alg.clone());
        let pairs = factor_univariate_over_k(&session, &flat).unwrap();
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].0, p_alg);
        assert_eq!(pairs[0].1, 1);
    }

    #[test]
    fn factor_quartic_reducible_over_q_yields_four_linear() {
        use giac_poly::Poly;

        let t = Var::from("__t");
        let tv = Poly::var("__t");
        let p = tv
            .pow(4)
            .mul_scalar(&Ratio::from_integer((-1).into()))
            .add(&tv.pow(3).mul_scalar(&Ratio::from_integer(4.into())))
            .add(&tv.pow(2).mul_scalar(&Ratio::from_integer((-2).into())))
            .add(&tv.mul_scalar(&Ratio::from_integer(4.into())))
            .add(&Poly::constant(Ratio::from_integer((-1).into())));
        let p_alg = super::super::poly::poly_algext_from_poly(&p).unwrap();
        let session = FieldSession::new(ExtensionField::rational());
        let flat = FlatUni::try_new(p_alg.clone(), MainVar::new(t.clone())).unwrap();
        let factors = factor_univariate_flat_over_k(&session, &flat).unwrap();
        assert_eq!(factors.len(), 4);
        assert!(factors.iter().all(|f| f.degree_wrt(&t) == 1));
        assert!(product_divides(&session, &p_alg, &factors, &t).unwrap());
    }

    #[test]
    fn factor_minpoly_x4_minus_4_over_sqrt2() {
        let k1 = k1_adjoin_sqrt2();
        let mb = vec![
            Ratio::from_integer(1.into()),
            Ratio::zero(),
            Ratio::zero(),
            Ratio::zero(),
            Ratio::from_integer((-4).into()),
        ];
        let factors = factor_minpoly_over_field(&mb, &k1).expect("factor");
        assert!(factors.len() >= 2);
        assert!(factors.iter().all(|f| f.degree_wrt(&common_minpoly_var()) == 2));
    }

    #[test]
    fn factor_x_fourth_minus_4_over_q() {
        let session = FieldSession::new(ExtensionField::rational());
        let p = x_fourth_minus_4(&session);
        let pairs = factor_univariate_over_k(&session, &flat(&session, p.clone())).unwrap();
        assert_eq!(pairs.len(), 4);
        assert!(pairs.iter().all(|(_, k)| *k == 1));
        let prod = pairs
            .iter()
            .map(|(f, k)| {
                (0..*k).try_fold(PolyAlgExt::ring_one(), |acc, _| acc.try_mul(f))
            })
            .try_fold(PolyAlgExt::ring_one(), |acc, t| acc.try_mul(&t?))
            .unwrap();
        let (_, r) = div_rem_wrt_algext(&session, &p, &prod, &x_var()).unwrap();
        assert!(r.is_zero());
        assert!(pairs.iter().all(|(f, _)| f.degree_wrt(&x_var()) == 1));
    }
}
