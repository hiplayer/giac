//! Multivariate gcd via subresultant PRS for [`Poly<C>`](crate::poly::Poly) with `C: FieldCoeff`.
//!
//! **Upstream:** `gcd_ext` / `gcdheu` coefficient-ring recursion over K.
//! **Ring context:** coefficients of nested polynomials are scalars in field K (`FieldCoeff`).

use std::collections::BTreeSet;

use crate::error::{EvalError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::poly_coeff::FieldCoeff;
use crate::univ_wrt::is_univariate_in;

// **Pipeline private** — `monomial_pow`
fn monomial_pow(var: &Var, exp: u64) -> crate::monomial::Monomial {
    let mut m = crate::monomial::Monomial::one();
    for _ in 0..exp {
        m = m.mul(&crate::monomial::Monomial::var(var.clone()));
    }
    m
}

/// **Stable** — coefficient of `var^exp` as an element of K[others].
pub fn coeff_wrt_field<C: FieldCoeff>(p: &Poly<C>, var: &Var, exp: u64) -> PolyResult<Poly<C>> {
    let mut out = Poly::ring_zero();
    let div = monomial_pow(var, exp);
    for (m, c) in &p.terms {
        if m.exp_of(var) == exp {
            if let Some(rest_m) = m.div_exact(&div) {
                out = out.try_add(&Poly::term(rest_m.clone(), c.clone()))?;
            }
        }
    }
    Ok(out)
}

// **Pipeline private** — `term_with_var_field`
fn term_with_var_field<C: FieldCoeff>(coeff: &Poly<C>, var: &Var, exp: u64) -> PolyResult<Poly<C>> {
    if exp == 0 {
        return Ok(coeff.clone());
    }
    Ok(coeff.try_mul(&Poly::ring_var(var.clone()).try_pow(exp)?)?)
}

// **Pipeline private** — sorted variables in `p`
fn vars_in_field<C: FieldCoeff>(p: &Poly<C>) -> Vec<Var> {
    let mut set = BTreeSet::new();
    for m in p.terms.keys() {
        for (v, e) in m.iter() {
            if e > 0 {
                set.insert(v.clone());
            }
        }
    }
    set.into_iter().collect()
}

// **Pipeline private** — `vars_union_field`
fn vars_union_field<C: FieldCoeff>(a: &Poly<C>, b: &Poly<C>) -> Vec<Var> {
    let mut set = BTreeSet::new();
    for p in [a, b] {
        for v in vars_in_field(p) {
            set.insert(v);
        }
    }
    set.into_iter().collect()
}

// **Pipeline private** — `main_var_for_gcd_field`
fn main_var_for_gcd_field<C: FieldCoeff>(a: &Poly<C>, b: &Poly<C>) -> Option<Var> {
    let vars = vars_union_field(a, b);
    vars.into_iter().min_by_key(|v| {
        a.degree_wrt(v)
            .max(b.degree_wrt(v))
            .max(1)
    })
}

/// **Stable (crate-internal)** — exact quotient `a/b` in K[others] when `b | a`.
pub(crate) fn div_exact_coeff_field<C: FieldCoeff>(
    a: &Poly<C>,
    b: &Poly<C>,
) -> Option<Poly<C>> {
    let (q, r) = div_rem_leading_field(a, b).ok()?;
    if r.is_zero() { Some(q) } else { None }
}

/// **Stable (crate-internal)** — multivariate leading-term division for `Poly<C>`.
pub(crate) fn div_rem_leading_field<C: FieldCoeff>(
    dividend: &Poly<C>,
    divisor: &Poly<C>,
) -> PolyResult<(Poly<C>, Poly<C>)> {
    if divisor.is_zero() {
        return Ok((Poly::ring_zero(), dividend.clone()));
    }
    let mut remainder = dividend.clone();
    let mut quotient = Poly::ring_zero();
    let Some((div_lt, div_lc)) = divisor.leading_term() else {
        return Ok((Poly::ring_zero(), remainder));
    };

    loop {
        let Some((r_lt, r_lc)) = remainder.leading_term() else {
            break;
        };
        if r_lt.degree() < div_lt.degree() {
            break;
        }
        if !r_lt.is_dividing(div_lt) {
            break;
        }
        let Some(q_m) = r_lt.div_exact(div_lt) else {
            break;
        };
        let q_c = r_lc.field_div(div_lc)?;
        let q_term = Poly::term(q_m, q_c);
        quotient = quotient.try_add(&q_term)?;
        remainder = remainder.try_sub(&q_term.try_mul(divisor)?)?;
    }
    Ok((quotient, remainder))
}

/// **Stable (crate-internal)** — division in K[others][var] (nested; not flat K[var]).
pub(crate) fn nested_div_rem_wrt_in_field<C: FieldCoeff>(
    a: &Poly<C>,
    b: &Poly<C>,
    var: &Var,
) -> PolyResult<(Poly<C>, Poly<C>)> {
    let mut remainder = a.clone();
    let mut quotient = Poly::ring_zero();
    let db = b.degree_wrt(var);
    if db == 0 {
        return Ok((Poly::ring_zero(), remainder));
    }
    let lc_b = coeff_wrt_field(b, var, db)?;
    if lc_b.is_zero() {
        return Ok((Poly::ring_zero(), remainder));
    }

    loop {
        let dr = remainder.degree_wrt(var);
        if dr < db || remainder.is_zero() {
            break;
        }
        let lc_r = coeff_wrt_field(&remainder, var, dr)?;
        let Some(q_coeff) = div_exact_coeff_field(&lc_r, &lc_b) else {
            break;
        };
        let d = dr - db;
        let q_term = term_with_var_field(&q_coeff, var, d)?;
        quotient = quotient.try_add(&q_term)?;
        remainder = remainder.try_sub(&q_term.try_mul(b)?)?;
    }
    Ok((quotient, remainder))
}

// **Pipeline private** — pseudo-remainder in K[others][var]
fn pseudo_rem_wrt_field<C: FieldCoeff>(a: &Poly<C>, b: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    if b.is_zero() {
        return Ok(a.clone());
    }
    let da = a.degree_wrt(var);
    let db = b.degree_wrt(var);
    if da < db {
        return Ok(a.clone());
    }
    let lc_b = coeff_wrt_field(b, var, db)?;
    if lc_b.is_zero() {
        return Ok(a.clone());
    }
    let exp = da - db + 1;
    let scaled = a.try_mul(&lc_b.try_pow(exp)?)?;
    let (_, r) = nested_div_rem_wrt_in_field(&scaled, b, var)?;
    Ok(primitive_part_wrt_field(&r, var)?)
}

// **Pipeline private** — monic normalize w.r.t. lex leading coefficient in K
fn scalar_monic_poly<C: FieldCoeff>(p: &Poly<C>) -> PolyResult<Poly<C>> {
    if p.is_zero() || p.is_one() {
        return Ok(p.clone());
    }
    let Some((_, lc)) = p.leading_term() else {
        return Ok(p.clone());
    };
    if lc.coeff_is_one() {
        return Ok(p.clone());
    }
    let inv = C::coeff_one().coeff_div(lc)?;
    let mut out = Poly::ring_zero();
    for (mon, c) in &p.terms {
        let scaled = c.coeff_mul(&inv)?;
        out = out.try_add(&Poly::term(mon.clone(), scaled))?;
    }
    Ok(out)
}

// **Pipeline private** — content in K[others] w.r.t. `var`
fn content_wrt_field<C: FieldCoeff>(p: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    let d = p.degree_wrt(var);
    let mut g = Poly::ring_zero();
    for e in 0..=d {
        let c = coeff_wrt_field(p, var, e)?;
        if c.is_zero() {
            continue;
        }
        g = if g.is_zero() {
            c
        } else {
            subresultant_gcd_field(&g, &c)?
        };
    }
    if g.is_zero() {
        Ok(Poly::ring_one())
    } else {
        scalar_monic_poly(&g)
    }
}

// **Pipeline private** — primitive part in K[others][var]
fn primitive_part_wrt_field<C: FieldCoeff>(p: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    let content = content_wrt_field(p, var)?;
    if content.is_one() {
        return Ok(p.clone());
    }
    let d = p.degree_wrt(var);
    let mut pp = Poly::ring_zero();
    for e in 0..=d {
        let c = coeff_wrt_field(p, var, e)?;
        if c.is_zero() {
            continue;
        }
        let q = nested_exact_quo_coeff_field(&c, &content, var)?;
        pp = pp.try_add(&term_with_var_field(&q, var, e)?)?;
    }
    Ok(pp)
}

// **Pipeline private** — exact quotient in coefficient ring K[others]
fn nested_exact_quo_coeff_field<C: FieldCoeff>(
    num: &Poly<C>,
    den: &Poly<C>,
    pivot: &Var,
) -> PolyResult<Poly<C>> {
    if den.is_zero() {
        return Err(EvalError::DivisionByZero);
    }
    if den.is_one() {
        return Ok(num.clone());
    }
    if is_univariate_in(num, pivot) && is_univariate_in(den, pivot) {
        return crate::univ_wrt::quo_exact_wrt(num, den, pivot);
    }
    let (q, r) = div_rem_leading_field(num, den)?;
    if r.is_zero() {
        Ok(q)
    } else {
        Err(EvalError::NotImplemented("poly division"))
    }
}

// **Pipeline private** — gcd when `var` side is constant in main
fn gcd_constant_wrt_field<C: FieldCoeff>(
    var: &Poly<C>,
    other: &Poly<C>,
    x: &Var,
) -> PolyResult<Poly<C>> {
    let d = other.degree_wrt(x);
    let mut g = var.clone();
    for e in 0..=d {
        let c = coeff_wrt_field(other, x, e)?;
        if !c.is_zero() {
            g = subresultant_gcd_field(&g, &c)?;
        }
    }
    scalar_monic_poly(&g)
}

// **Pipeline private** — subresultant gcd w.r.t. main variable
fn subresultant_gcd_wrt_field<C: FieldCoeff>(
    a: &Poly<C>,
    b: &Poly<C>,
    var: &Var,
) -> PolyResult<Poly<C>> {
    let da = a.degree_wrt(var);
    let db = b.degree_wrt(var);

    if da == 0 && db == 0 {
        let mut vars = vars_union_field(a, b);
        vars.retain(|v| v != var);
        if vars.is_empty() {
            return scalar_monic_poly(a);
        }
        return subresultant_gcd_wrt_field(a, b, &vars[0]);
    }
    if da == 0 {
        return gcd_constant_wrt_field(a, b, var);
    }
    if db == 0 {
        return gcd_constant_wrt_field(b, a, var);
    }

    let mut a = primitive_part_wrt_field(a, var)?;
    let mut b = primitive_part_wrt_field(b, var)?;
    if a.degree_wrt(var) < b.degree_wrt(var) {
        std::mem::swap(&mut a, &mut b);
    }
    while !b.is_zero() {
        let db = b.degree_wrt(var);
        if db == 0 {
            return gcd_constant_wrt_field(&b, &a, var);
        }
        let r = pseudo_rem_wrt_field(&a, &b, var)?;
        if r.is_zero() {
            return scalar_monic_poly(&b);
        }
        let dr = r.degree_wrt(var);
        if dr >= a.degree_wrt(var) && !r.is_zero() {
            return scalar_monic_poly(&a);
        }
        a = b;
        b = r;
    }
    scalar_monic_poly(&a)
}

// **Pipeline private** — scalar content gcd in K
fn field_scalar_content_gcd<C: FieldCoeff>(a: &C, b: &C) -> PolyResult<C> {
    if a.coeff_is_zero() {
        return Ok(b.clone());
    }
    if b.coeff_is_zero() {
        return Ok(a.clone());
    }
    let mut x = a.clone();
    let mut y = b.clone();
    loop {
        if y.coeff_is_zero() {
            break;
        }
        match x.coeff_div(&y) {
            Ok(q) => {
                let qy = q.coeff_mul(&y)?;
                x = x.coeff_sub(&qy)?;
                std::mem::swap(&mut x, &mut y);
            }
            Err(_) => return Ok(C::coeff_one()),
        }
    }
    if x.coeff_is_zero() {
        Ok(C::coeff_zero())
    } else {
        Ok(C::coeff_one())
    }
}

// **Pipeline private** — gcd univariate in one variable via flat K[var]
fn gcd_univariate_field<C: FieldCoeff>(a: &Poly<C>, b: &Poly<C>, var: &Var) -> PolyResult<Poly<C>> {
    crate::univ_wrt::gcd_wrt(a, b, var)
}

/// **Stable** — multivariate gcd over field coefficients K via subresultant PRS.
pub fn subresultant_gcd_field<C: FieldCoeff>(a: &Poly<C>, b: &Poly<C>) -> PolyResult<Poly<C>> {
    if a.is_zero() {
        return Ok(if b.is_zero() {
            Poly::ring_zero()
        } else {
            scalar_monic_poly(b)?
        });
    }
    if b.is_zero() {
        return scalar_monic_poly(a);
    }

    let mut g_const = C::coeff_zero();
    for c in a.terms.values().chain(b.terms.values()) {
        g_const = if g_const.coeff_is_zero() {
            c.clone()
        } else {
            field_scalar_content_gcd(&g_const, c)?
        };
    }

    let a = primitive_part_scalar_field(a)?;
    let b = primitive_part_scalar_field(b)?;

    let vars = vars_union_field(&a, &b);
    if vars.len() == 1 {
        let v = &vars[0];
        let g = gcd_univariate_field(&a, &b, v)?;
        return scale_gcd_by_scalar_content(g, &g_const);
    }

    let Some(main) = main_var_for_gcd_field(&a, &b) else {
        return Ok(Poly::ring_constant(if g_const.coeff_is_zero() {
            C::coeff_one()
        } else {
            g_const
        }));
    };

    if is_univariate_in(&a, &main) && is_univariate_in(&b, &main) {
        let g = gcd_univariate_field(&a, &b, &main)?;
        return scale_gcd_by_scalar_content(g, &g_const);
    }

    let g_poly = subresultant_gcd_wrt_field(&a, &b, &main)?;
    scale_gcd_by_scalar_content(g_poly, &g_const)
}

// **Pipeline private** — divide out scalar content from K
fn primitive_part_scalar_field<C: FieldCoeff>(p: &Poly<C>) -> PolyResult<Poly<C>> {
    let mut g: Option<C> = None;
    for c in p.terms.values() {
        if c.coeff_is_zero() {
            continue;
        }
        g = Some(match g {
            None => c.clone(),
            Some(prev) => field_scalar_content_gcd(&prev, c)?,
        });
    }
    let Some(g) = g else {
        return Ok(p.clone());
    };
    if g.coeff_is_one() {
        return Ok(p.clone());
    }
    let inv = C::coeff_one().coeff_div(&g)?;
    let mut out = Poly::ring_zero();
    for (m, c) in &p.terms {
        out = out.try_add(&Poly::term(m.clone(), c.coeff_mul(&inv)?))?;
    }
    Ok(out)
}

// **Pipeline private** — attach scalar content to polynomial gcd
fn scale_gcd_by_scalar_content<C: FieldCoeff>(g: Poly<C>, c: &C) -> PolyResult<Poly<C>> {
    if c.coeff_is_zero() || c.coeff_is_one() || g.is_zero() {
        return scalar_monic_poly(&g);
    }
    if g.is_one() {
        return Ok(Poly::ring_constant(c.clone()));
    }
    let mut out = Poly::ring_zero();
    for (m, coeff) in &g.terms {
        out = out.try_add(&Poly::term(m.clone(), coeff.coeff_mul(c)?))?;
    }
    scalar_monic_poly(&out)
}

#[cfg(test)]
mod tests {
    use num_rational::Ratio;

    use super::*;
    use crate::monomial::Var;
    use crate::poly::Poly;

    #[test]
    fn gcd_field_xy_and_y_matches_rational() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let g = subresultant_gcd_field(&x.try_mul(&y).unwrap(), &y).unwrap();
        assert_eq!(g, y);
    }

    #[test]
    fn gcd_field_univariate_matches_univ_wrt() {
        let x = Var::from("x");
        let p = Poly::var(x.clone())
            .try_pow(3)
            .unwrap()
            .try_sub(&Poly::one())
            .unwrap();
        let q = Poly::var(x.clone())
            .try_pow(2)
            .unwrap()
            .try_sub(&Poly::one())
            .unwrap();
        assert_eq!(
            subresultant_gcd_field(&p, &q).unwrap(),
            Poly::var(x).try_sub(&Poly::one()).unwrap()
        );
    }

    #[test]
    fn gcd_field_bivariate_linear_and_quadratic_coprime() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let three = Poly::constant(Ratio::from_integer(3.into()));
        let two = Poly::constant(Ratio::from_integer(2.into()));
        let a = three
            .try_mul(&x)
            .unwrap()
            .try_sub(&two.try_mul(&y).unwrap())
            .unwrap();
        let b = x.try_pow(2)
            .unwrap()
            .try_add(&x.try_mul(&y).unwrap())
            .unwrap()
            .try_add(&y.try_pow(2).unwrap())
            .unwrap();
        assert!(subresultant_gcd_field(&a, &b).unwrap().is_one());
    }
}
