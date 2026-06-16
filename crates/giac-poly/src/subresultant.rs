//! Multivariate gcd via subresultant PRS in ℚ[y₁,…][x].

use std::collections::BTreeSet;

use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::univariate_degree;
use crate::univariate::gcd_univariate;

/// Coefficient of `var^exp` (quotient by `var^exp` on matching terms).
pub fn coeff_wrt(p: &Poly, var: &Var, exp: u64) -> Poly {
    let mut out = Poly::zero();
    let div = monomial_pow(var, exp);
    for (m, c) in &p.terms {
        if m.exp_of(var) == exp {
            if let Some(rest_m) = m.div_exact(&div) {
                out = out.add(&Poly::term(rest_m, c.clone()));
            }
        }
    }
    out
}

fn monomial_pow(var: &Var, exp: u64) -> crate::monomial::Monomial {
    let mut m = crate::monomial::Monomial::one();
    for _ in 0..exp {
        m = m.mul(&crate::monomial::Monomial::var(var.clone()));
    }
    m
}

fn term_with_var(coeff: &Poly, var: &Var, exp: u64) -> Poly {
    if exp == 0 {
        return coeff.clone();
    }
    coeff.mul(&Poly::var(var.clone()).pow(exp))
}

fn vars_in(p: &Poly) -> Vec<Var> {
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

fn vars_union(a: &Poly, b: &Poly) -> Vec<Var> {
    let mut set = BTreeSet::new();
    for p in [a, b] {
        for v in vars_in(p) {
            set.insert(v);
        }
    }
    set.into_iter().collect()
}

fn is_univariate_in(p: &Poly, var: &Var) -> bool {
    p.terms
        .keys()
        .all(|m| m.iter().all(|(v, _)| v == var))
}

/// Main variable for gcd: minimum of max(deg_x a, deg_x b) over variables.
fn main_var_for_gcd(a: &Poly, b: &Poly) -> Option<Var> {
    let vars = vars_union(a, b);
    vars.into_iter().min_by_key(|v| {
        univariate_degree(a, v)
            .max(univariate_degree(b, v))
            .max(1)
    })
}

fn rational_primitive(p: &Poly) -> Poly {
    p.primitive_part()
}

/// Exact division `a / b` in the coefficient ring when `b | a`.
pub(crate) fn div_exact_coeff(a: &Poly, b: &Poly) -> Option<Poly> {
    let (q, r) = a.div_rem(b);
    if r.is_zero() {
        Some(q)
    } else {
        None
    }
}

/// Division in ℚ[others][var]: divide `a` by `b` treating them as univariate in `var`.
pub(crate) fn univariate_div_rem_wrt(a: &Poly, b: &Poly, var: &Var) -> (Poly, Poly) {
    let mut remainder = a.clone();
    let mut quotient = Poly::zero();
    let db = univariate_degree(b, var);
    if db == 0 {
        return (Poly::zero(), remainder);
    }
    let lc_b = coeff_wrt(b, var, db);
    if lc_b.is_zero() {
        return (Poly::zero(), remainder);
    }

    loop {
        let dr = univariate_degree(&remainder, var);
        if dr < db || remainder.is_zero() {
            break;
        }
        let lc_r = coeff_wrt(&remainder, var, dr);
        let Some(q_coeff) = div_exact_coeff(&lc_r, &lc_b) else {
            break;
        };
        let d = dr - db;
        let q_term = term_with_var(&q_coeff, var, d);
        quotient = quotient.add(&q_term);
        remainder = remainder.sub(&q_term.mul(b));
    }
    (quotient, remainder)
}

/// Pseudo-remainder of `a` modulo `b` w.r.t. `var` (coefficients in ℚ[others]).
fn pseudo_rem_wrt(a: &Poly, b: &Poly, var: &Var) -> Poly {
    if b.is_zero() {
        return a.clone();
    }
    let da = univariate_degree(a, var);
    let db = univariate_degree(b, var);
    if da < db {
        return a.clone();
    }
    let lc_b = coeff_wrt(b, var, db);
    if lc_b.is_zero() {
        return a.clone();
    }
    let exp = da - db + 1;
    let scaled = a.mul(&lc_b.pow(exp));
    let (_, r) = univariate_div_rem_wrt(&scaled, b, var);
    r.primitive_part()
}

fn content_wrt(p: &Poly, var: &Var) -> Poly {
    content_wrt_impl(p, var)
}

fn primitive_part_wrt(p: &Poly, var: &Var) -> Poly {
    primitive_part_wrt_impl(p, var)
}

/// Gcd of coefficient polynomials w.r.t. `var` (used by factor and gcd).
pub(crate) fn content_wrt_impl(p: &Poly, var: &Var) -> Poly {
    let d = univariate_degree(p, var);
    let mut g = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt(p, var, e);
        if c.is_zero() {
            continue;
        }
        g = if g.is_zero() {
            c
        } else {
            subresultant_gcd(&g, &c)
        };
    }
    if g.is_zero() {
        Poly::one()
    } else {
        g.monic()
    }
}

pub(crate) fn primitive_part_wrt_impl(p: &Poly, var: &Var) -> Poly {
    let content = content_wrt_impl(p, var);
    if content.is_one() {
        return p.clone();
    }
    let d = univariate_degree(p, var);
    let mut pp = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt(p, var, e);
        if c.is_zero() {
            continue;
        }
        let Some(q) = div_exact_coeff(&c, &content) else {
            return p.clone();
        };
        pp = pp.add(&term_with_var(&q, var, e));
    }
    pp
}

fn gcd_constant_wrt(var: &Poly, other: &Poly, x: &Var) -> Poly {
    let d = univariate_degree(other, x);
    let mut g = var.clone();
    for e in 0..=d {
        let c = coeff_wrt(other, x, e);
        if !c.is_zero() {
            g = subresultant_gcd(&g, &c);
        }
    }
    g.monic()
}

fn subresultant_gcd_wrt(a: &Poly, b: &Poly, var: &Var) -> Poly {
    let da = univariate_degree(a, var);
    let db = univariate_degree(b, var);

    if da == 0 && db == 0 {
        let mut vars = vars_union(a, b);
        vars.retain(|v| v != var);
        if vars.is_empty() {
            return a.monic();
        }
        return subresultant_gcd_wrt(a, b, &vars[0]);
    }
    if da == 0 {
        return gcd_constant_wrt(a, b, var);
    }
    if db == 0 {
        return gcd_constant_wrt(b, a, var);
    }

    let mut a = primitive_part_wrt(a, var);
    let mut b = primitive_part_wrt(b, var);
    if univariate_degree(&a, var) < univariate_degree(&b, var) {
        std::mem::swap(&mut a, &mut b);
    }
    while !b.is_zero() {
        let db = univariate_degree(&b, var);
        if db == 0 {
            return gcd_constant_wrt(&b, &a, var);
        }
        let r = pseudo_rem_wrt(&a, &b, var);
        if r.is_zero() {
            return b.monic();
        }
        let dr = univariate_degree(&r, var);
        if dr >= univariate_degree(&a, var) && !r.is_zero() {
            return a.monic();
        }
        a = b;
        b = r;
    }
    a.monic()
}

/// Gcd of multivariate polynomials over ℚ via subresultant PRS.
pub fn subresultant_gcd(a: &Poly, b: &Poly) -> Poly {
    if a.is_zero() {
        return b.clone();
    }
    if b.is_zero() {
        return a.clone();
    }

    let ca = a.content();
    let cb = b.content();
    let mut g_const = if ca.is_zero() {
        cb.clone()
    } else if cb.is_zero() {
        ca.clone()
    } else {
        crate::poly::integer_content_gcd(&ca, &cb)
    };

    let mut a = rational_primitive(a);
    let mut b = rational_primitive(b);

    let vars = vars_union(&a, &b);
    if vars.len() == 1 {
        let g = gcd_univariate(&a, &b, &vars[0]);
        return scale_gcd_by_content(g, &g_const);
    }

    let Some(main) = main_var_for_gcd(&a, &b) else {
        return Poly::constant(g_const).monic();
    };

    if is_univariate_in(&a, &main) && is_univariate_in(&b, &main) {
        let g = gcd_univariate(&a, &b, &main);
        return scale_gcd_by_content(g, &g_const);
    }

    let g_poly = subresultant_gcd_wrt(&a, &b, &main);
    scale_gcd_by_content(g_poly, &g_const)
}

fn scale_gcd_by_content(mut g: Poly, c: &Ratio<num_bigint::BigInt>) -> Poly {
    if c.is_one() || g.is_zero() {
        return g.monic();
    }
    if g.is_one() && g.terms.len() == 1 {
        return Poly::constant(c.clone()).monic();
    }
    g = g.mul_scalar(c);
    g.monic()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_rational::Ratio;

    #[test]
    fn gcd_xy_and_y() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let g = subresultant_gcd(&x.mul(&y), &y);
        assert_eq!(g, y);
    }

    #[test]
    fn gcd_bivariate_linear_and_quadratic() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let a = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x)
            .sub(&Poly::constant(Ratio::from_integer(2.into())).mul(&y));
        let b = x.pow(2).add(&x.mul(&y)).add(&y.pow(2));
        let g = subresultant_gcd(&a, &b);
        assert!(g.is_one());
    }

    #[test]
    fn gcd_univariate_matches_subresultant() {
        let x = Poly::var("x");
        let p = x.pow(3).sub(&Poly::one());
        let q = x.pow(2).sub(&Poly::one());
        assert_eq!(subresultant_gcd(&p, &q), x.sub(&Poly::one()));
    }

    #[test]
    fn content_wrt_y_of_xy() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.mul(&y);
        assert_eq!(content_wrt(&p, &Var::from("x")), y);
    }
}
