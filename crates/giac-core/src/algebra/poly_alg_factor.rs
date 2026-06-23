//! Univariate factorization over K for `Poly<AlgExtC>` (T2-1 / L3-1).
//!
//! **Upstream:** `ext_factor` / `ext_factor_nodegck` — sqff → linear / quadratic split / roots / witness.

use giac_poly::{FlatUni, MainVar, PolyCoeff, Var};

use crate::error::EvalError;

use super::field_session::FieldSession;
use super::poly::PolyAlgExt;
use super::poly_alg_coeff::AlgExtCPolyCoeff;
use super::poly_alg_ops::{
    div_rem_wrt_algext, normalize_algext_poly, split_quadratic_factor,
};
use super::poly_roots::poly_algext_roots;

/// **Stable** — factor `flat` in K[var] as `(factor, multiplicity)` pairs.
pub fn factor_univariate_over_k(
    session: &FieldSession,
    flat: &FlatUni<AlgExtCPolyCoeff>,
) -> Result<Vec<(PolyAlgExt, usize)>, EvalError> {
    let var = flat.var().as_var().clone();
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
    use giac_poly::PolyCoeff;

    use super::*;
    use crate::algebra::ext_tower::ExtensionField;
    use crate::algebra::test_fixtures::k1_adjoin_sqrt2;

    fn x_var() -> Var {
        Var::from("x")
    }

    fn flat(session: &FieldSession, p: PolyAlgExt) -> FlatUni<AlgExtCPolyCoeff> {
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
    fn factor_x_fourth_minus_4_over_q() {
        let session = FieldSession::new(ExtensionField::rational());
        let p = x_fourth_minus_4(&session);
        let pairs = factor_univariate_over_k(&session, &flat(&session, p.clone())).unwrap();
        assert_eq!(pairs.len(), 2);
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
        assert!(pairs.iter().all(|(f, _)| f.degree_wrt(&x_var()) == 2));
    }
}
