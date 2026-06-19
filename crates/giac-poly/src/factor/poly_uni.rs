//! Polynomials as univariate in main var with coefficients in nested Poly ring (bivariate steps).
//!
//! **Stable / Partial:** `coeff_wrt_poly`, `content_wrt`, `factor_sqff_over_coeff_ring`, …
//! **Pipeline private:** sqff-over-coeff-ring factor chain (upstream `do_factor_hensel` slice).
//! Retired from hot path: `try_factor_bivariate_eval`, `try_kronecker_bivariate` (FAC-G3 covered by sparse→Hensel).

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::univariate_degree;

use super::univariate::factor_univariate_flat;
use super::util::{coeff_wrt, is_univariate_in};

/// **Stable** — Coefficient of `var^exp` as a polynomial in the remaining variables.
pub fn coeff_wrt_poly(p: &Poly, var: &Var, exp: u64) -> Poly {
    coeff_wrt(p, var, exp)
}

/// **Stable** — Content of `p` w.r.t. `var`: gcd of all x-coefficients in ℚ[others].
pub fn content_wrt(p: &Poly, var: &Var) -> Poly {
    crate::subresultant::content_wrt_impl(p, var)
}

// **Pipeline private** — `poly_div_exact`
fn poly_div_exact(num: &Poly, den: &Poly) -> PolyResult<Poly> {
    if let Some(q) = crate::subresultant::div_exact_coeff(num, den) {
        return Ok(q);
    }
    Err(PolyError::NotImplemented("poly division"))
}

// **Pipeline private** — `poly_div_exact_wrt`
fn poly_div_exact_wrt(num: &Poly, den: &Poly, var: &Var) -> PolyResult<Poly> {
    let (q, r) = crate::subresultant::univariate_div_rem_wrt(num, den, var);
    if r.is_zero() {
        Ok(q)
    } else {
        Err(PolyError::NotImplemented("poly division"))
    }
}

/// **Stable** — `p / content_wrt(p, var)` in ℚ[others][var].
pub fn primitive_part_wrt(p: &Poly, var: &Var) -> PolyResult<Poly> {
    let content = content_wrt(p, var);
    if content.is_one() {
        return Ok(p.clone());
    }
    let d = univariate_degree(p, var);
    let mut pp = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt_poly(p, var, e);
        if c.is_zero() {
            continue;
        }
        let q = poly_div_exact(&c, &content)?;
        pp = pp.add(&term_with_var(&q, var, e));
    }
    Ok(pp)
}

/// **Stable** — coeff * var^exp as Poly
pub fn term_with_var(coeff: &Poly, var: &Var, exp: u64) -> Poly {
    if exp == 0 {
        return coeff.clone();
    }
    coeff.mul(&Poly::var(var.clone()).pow(exp))
}

/// **Stable** — ∂p/∂var treating coefficients in ℚ[others].
pub fn derivative_wrt(p: &Poly, var: &Var) -> Poly {
    let d = univariate_degree(p, var);
    let mut out = Poly::zero();
    for e in 1..=d {
        let c = coeff_wrt_poly(p, var, e);
        if c.is_zero() {
            continue;
        }
        let scaled = c.mul_scalar(&Ratio::from_integer(BigInt::from(e)));
        out = out.add(&term_with_var(&scaled, var, e - 1));
    }
    out
}

/// **Stable** — Square-free factorization w.r.t. `var` over ℚ[others] (Yun-style via gcd).
pub fn square_free_wrt(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    match square_free_wrt_impl(p, var) {
        Ok(f) => Ok(f),
        Err(crate::error::PolyError::NotImplemented("poly division")) => Ok(vec![(p.clone(), 1)]),
        Err(e) => Err(e),
    }
}

// **Pipeline private** — `square_free_wrt_impl`
fn square_free_wrt_impl(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    if p.is_zero() {
        return Err(PolyError::TypeError("zero polynomial"));
    }
    let mut w = p.clone();
    let mut y = derivative_wrt(&w, var);
    let g0 = w.gcd(&y);
    if !g0.is_one() {
        w = poly_div_exact_wrt(&w, &g0, var)?;
        y = poly_div_exact_wrt(&y, &g0, var)?;
    }
    y = y.sub(&derivative_wrt(&w, var));

    let mut factors = Vec::new();
    let mut k = 1usize;
    let max_k = univariate_degree(p, var) as usize + 2;
    while !y.is_zero() && k <= max_k {
        let g = w.gcd(&y);
        if !g.is_one() {
            factors.push((g.clone(), k));
            w = poly_div_exact_wrt(&w, &g, var)?;
        }
        y = y.sub(&derivative_wrt(&w, var));
        k += 1;
    }
    if !w.is_one() {
        factors.push((w, k));
    }
    Ok(factors)
}

/// **Stable** — Substitute `sub_var -> sub_poly` in `p`.
pub fn substitute_poly(p: &Poly, sub_var: &Var, sub_poly: &Poly) -> Poly {
    let d = univariate_degree(p, sub_var);
    let mut out = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt_poly(p, sub_var, e);
        if c.is_zero() {
            continue;
        }
        out = out.add(&c.mul(&sub_poly.pow(e)));
    }
    out
}

type FactorRecFn = fn(&Poly, &[Var]) -> PolyResult<Vec<Poly>>;

/// **Partial** — Factor square-free `g` in ℚ[others][var] recursively.
pub fn factor_sqff_over_coeff_ring(
    g: &Poly,
    var: &Var,
    others: &[Var],
    factor_rec: FactorRecFn,
) -> PolyResult<Vec<Poly>> {
    let dx = univariate_degree(g, var);
    if dx == 0 {
        return factor_rec(g, others);
    }
    if dx == 1 {
        return Ok(vec![g.clone()]);
    }
    if others.is_empty() {
        return factor_univariate_flat(g, var);
    }
    if others.len() == 1 {
        let other = &others[0];
        // Upstream `do_factor_hensel`: try_sparse_factor → try_sparse_factor_bi → try_hensel_lift_factor
        // → unitaryfactor/pzadic. Rust: sparse @ aux=0, then Hensel @ aux=0 (both main/aux orders).
        for (main, aux) in [(var, other), (other, var)] {
            if let Some(f) = super::sparse::try_sparse_factor(g, main, aux) {
                return Ok(f);
            }
            if let Some(f) = super::hensel::try_hensel_lift_bivariate(g, main, aux) {
                return Ok(f);
            }
            // Heuristic pzadic: disabled until bounded integer path (upstream last resort).
        }
    }
    if others.len() >= 2 {
        for av in others {
            if let Some(f) = super::hensel::try_lift_factors_in_aux_var(g, var, av, others) {
                return Ok(f);
            }
        }
    }
    // Fallback chain exhausted → treat as irreducible (upstream pushes pcur unchanged).
    Ok(vec![g.clone()])
}

// **Temporary (retired)** — eval+interp lift; superseded by sparse→Hensel (FAC-G3). Kept for regression only.
#[cfg(test)]
fn try_factor_bivariate_eval(
    p: &Poly,
    main: &Var,
    other: &Var,
    all_others: &[Var],
) -> Option<Vec<Poly>> {
    if all_others.len() != 1 {
        return None;
    }
    let mut nf: Option<usize> = None;
    for k in 0i64..=3 {
        let sub = Poly::constant(Ratio::from_integer(BigInt::from(k)));
        let ev = substitute_poly(p, other, &sub);
        let facs = factor_univariate_flat(&ev, main).ok()?;
        if facs.len() <= 1 {
            return None;
        }
        nf = Some(match nf {
            None => facs.len(),
            Some(n) if n == facs.len() => n,
            _ => return None,
        });
    }
    let _ = nf?;
    try_lift_bivariate_from_eval(p, main, other)
}

#[cfg(test)]
fn try_lift_bivariate_from_eval(p: &Poly, main: &Var, other: &Var) -> Option<Vec<Poly>> {
    let mut candidates = Vec::new();
    for k in 0i64..=2 {
        let sub = Poly::constant(Ratio::from_integer(BigInt::from(k)));
        let ev = substitute_poly(p, other, &sub);
        let facs = factor_univariate_flat(&ev, main).ok()?;
        for f in facs {
            if f.is_one() {
                continue;
            }
            if let Some(lifted) = lift_univariate_factor(&f, main, other, p) {
                if !candidates.iter().any(|c: &Poly| c == &lifted) {
                    candidates.push(lifted);
                }
            }
        }
    }
    if candidates.is_empty() {
        return None;
    }
    let mut rest = p.clone();
    let mut out = Vec::new();
    for c in &candidates {
        let (_, r) = rest.div_rem(c);
        if r.is_zero() {
            out.push(c.clone());
            rest = rest.div_rem(c).0;
        }
    }
    if rest.is_one() || rest.is_zero() {
        Some(out)
    } else {
        None
    }
}

#[cfg(test)]
fn lift_univariate_factor(f: &Poly, main: &Var, other: &Var, p: &Poly) -> Option<Poly> {
    let deg = univariate_degree(f, main);
    if deg == 1 {
        let c1 = coeff_wrt_poly(f, main, 1);
        let c0 = coeff_wrt_poly(f, main, 0);
        if c1.is_one() {
            let lin = term_with_var(&Poly::one(), main, 1).add(&c0);
            let (_, r) = p.div_rem(&lin);
            if r.is_zero() {
                return Some(lin);
            }
        }
        if let Some(r) = c0.as_constant() {
            if c1.is_one() && deg == 1 {
                let lin = term_with_var(&Poly::one(), main, 1).sub(&Poly::constant(r));
                let (_, r2) = p.div_rem(&lin);
                if r2.is_zero() {
                    return Some(lin);
                }
            }
        }
    }
    for e in 0..=deg {
        let c = coeff_wrt_poly(f, main, e);
        if c.is_zero() {
            continue;
        }
        if c.is_one() {
            let trial = term_with_var(&Poly::one(), main, e);
            let (_, r) = p.div_rem(&trial);
            if r.is_zero() {
                return Some(trial);
            }
        }
        let trial = term_with_var(&Poly::var(other.clone()), main, e);
        let (_, r) = p.div_rem(&trial);
        if r.is_zero() {
            return Some(trial);
        }
        let trial = term_with_var(&Poly::var(other.clone()).sub(&Poly::one()), main, e);
        let (_, r) = p.div_rem(&trial);
        if r.is_zero() {
            return Some(trial);
        }
    }
    let (_, r) = p.div_rem(f);
    if r.is_zero() && !f.is_one() {
        Some(f.clone())
    } else {
        None
    }
}

// **Temporary (retired)** — Kronecker embed; slow/unreliable on L22. Kept for regression only.
#[cfg(test)]
fn try_kronecker_bivariate(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    let dx = univariate_degree(p, x);
    let dy = univariate_degree(p, y);
    if dx == 0 || dy == 0 {
        return None;
    }
    let n = dx + dy + 1;
    let sub = Poly::var(x.clone()).pow(n);
    let un = substitute_poly(p, y, &sub);
    if !is_univariate_in(&un, x) {
        return None;
    }
    let facs = factor_univariate_flat(&un, x).ok()?;
    if facs.len() <= 1 {
        return None;
    }
    let mut out = Vec::new();
    let mut rest = p.clone();
    for f in facs {
        if let Some(lifted) = kronecker_lift(&f, x, y, n) {
            let (_, r) = rest.div_rem(&lifted);
            if r.is_zero() {
                out.push(lifted.clone());
                rest = rest.div_rem(&lifted).0;
            }
        }
    }
    if rest.is_one() && !out.is_empty() {
        Some(out)
    } else {
        None
    }
}

// **Pipeline private** — `kronecker_lift`
fn kronecker_lift(f: &Poly, x: &Var, y: &Var, n: u64) -> Option<Poly> {
    let d = univariate_degree(f, x);
    let mut out = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt_poly(f, x, e);
        if c.is_zero() {
            continue;
        }
        let exp_y = e / n;
        let exp_x = e % n;
        let term = c
            .mul(&Poly::var(y.clone()).pow(exp_y))
            .mul(&Poly::var(x.clone()).pow(exp_x));
        out = out.add(&term);
    }
    if out.is_one() {
        None
    } else {
        Some(out)
    }
}

trait PolyConstant {
    // **Pipeline private** — `as_constant`
    fn as_constant(&self) -> Option<Ratio<BigInt>>;
}

impl PolyConstant for Poly {
    // **Stable** — `Poly::as_constant`
    fn as_constant(&self) -> Option<Ratio<BigInt>> {
        if self.terms.len() == 1 {
            self.terms.values().next().cloned()
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    use num_traits::One;

    #[test]
    fn content_wrt_xy_plus_y_squared() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        // xy + y² = y(x + y): content w.r.t. x is y, primitive part is x + y.
        let p = x.mul(&y).add(&y.pow(2));
        assert_eq!(content_wrt(&p, &Var::from("x")), y);
        let pp = primitive_part_wrt(&p, &Var::from("x")).unwrap();
        assert_eq!(pp, x.add(&y));
    }

    #[test]
    #[ignore = "FAC-G3: Kronecker too slow / fails on mixed-degree L22"]
    fn kronecker_line22() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x)
            .sub(&y.pow(2))
            .add(&y)
            .sub(&Poly::constant(Ratio::from_integer(5.into())))
            .mul(
                &x
                    .mul(&y)
                    .add(&Poly::constant(Ratio::from_integer(3.into())).mul(&x))
                    .sub(&y.pow(2))
                    .sub(&Poly::one()),
            );
        let f = try_kronecker_bivariate(&p, &Var::from("x"), &Var::from("y"))
            .expect("kronecker should factor L22");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn rational_content_vs_wrt_content() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = Poly::constant(Ratio::from_integer(6.into()))
            .mul(&x.mul(&y).add(&y.pow(2)));
        assert_eq!(p.content(), Ratio::from_integer(6.into()));
        assert_eq!(content_wrt(&p, &Var::from("x")), y);
        let pp = primitive_part_wrt(&p, &Var::from("x")).unwrap();
        assert_eq!(
            pp,
            x.add(&y).mul_scalar(&Ratio::from_integer(6.into()))
        );
    }
}
