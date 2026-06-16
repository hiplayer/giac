//! Bivariate factorization via evaluation + coefficient interpolation (Hensel MVP).

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::univariate_degree;

use super::poly_uni::{coeff_wrt_poly, substitute_poly, term_with_var};
use super::univariate::factor_univariate_flat;

fn as_rational_constant(p: &Poly) -> Option<Ratio<BigInt>> {
    if p.is_zero() {
        return Some(Ratio::zero());
    }
    if p.terms.len() == 1 {
        let m = p.terms.keys().next()?;
        if m.is_const() {
            return p.terms.values().next().cloned();
        }
    }
    None
}

fn linear_root(f: &Poly, x: &Var) -> Option<Ratio<BigInt>> {
    if univariate_degree(f, x) != 1 {
        return None;
    }
    let lc = as_rational_constant(&coeff_wrt_poly(f, x, 1))?;
    if !lc.is_one() {
        return None;
    }
    let c0 = as_rational_constant(&coeff_wrt_poly(f, x, 0))?;
    Some(-c0)
}

/// Sort key for matching univariate factors across `y = k` evaluations.
fn factor_match_key(f: &Poly, x: &Var) -> (u64, Ratio<BigInt>, Ratio<BigInt>) {
    if let Some(r) = linear_root(f, x) {
        return (1, r, Ratio::zero());
    }
    let d = univariate_degree(f, x);
    if d == 2 {
        let a = as_rational_constant(&coeff_wrt_poly(f, x, 1)).unwrap_or_else(Ratio::zero);
        let b = as_rational_constant(&coeff_wrt_poly(f, x, 0)).unwrap_or_else(Ratio::zero);
        return (2, a, b);
    }
    let c0 = as_rational_constant(&coeff_wrt_poly(f, x, 0)).unwrap_or_else(Ratio::zero);
    (d, c0, Ratio::zero())
}

fn lagrange_interpolate_y(y: &Var, points: &[(i64, Ratio<BigInt>)]) -> Poly {
    let mut out = Poly::zero();
    for (i, (yi, fi)) in points.iter().enumerate() {
        let mut basis = Poly::constant(fi.clone());
        for (j, (yj, _)) in points.iter().enumerate() {
            if i == j {
                continue;
            }
            let num = Poly::var(y.clone())
                .sub(&Poly::constant(Ratio::from_integer(BigInt::from(*yj))));
            let den = Ratio::from_integer(BigInt::from(*yi - *yj));
            basis = basis.mul(&num.mul_scalar(&(Ratio::one() / den)));
        }
        out = out.add(&basis);
    }
    out
}

fn lift_factor_from_evals(samples: &[Poly], x: &Var, y: &Var) -> Option<Poly> {
    let deg = univariate_degree(&samples[0], x);
    let mut out = Poly::zero();
    for e in 0..=deg {
        let mut points = Vec::with_capacity(samples.len());
        for (k, f) in samples.iter().enumerate() {
            let c = coeff_wrt_poly(f, x, e);
            let val = as_rational_constant(&c)?;
            points.push((k as i64, val));
        }
        let cy = lagrange_interpolate_y(y, &points);
        out = out.add(&term_with_var(&cy, x, e));
    }
    Some(out)
}

/// Factor `p(x,y)` in ℚ[y][x] by evaluating at `y = 0..dy` and lifting coefficients.
pub fn try_hensel_lift_bivariate(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    let dx = univariate_degree(p, x);
    let dy = univariate_degree(p, y);
    if dx == 0 || dy == 0 {
        return None;
    }

    let npts = (dy + 1).max(3) as i64;
    let mut evals: Vec<Vec<Poly>> = Vec::new();
    for k in 0..npts {
        let pk = substitute_poly(p, y, &Poly::constant(Ratio::from_integer(BigInt::from(k))));
        let facs = factor_univariate_flat(&pk, x).ok()?;
        if facs.len() <= 1 {
            return None;
        }
        if !evals.is_empty() && facs.len() != evals[0].len() {
            return None;
        }
        evals.push(facs);
    }

    let nf = evals[0].len();
    let mut tracks: Vec<Vec<Poly>> = vec![Vec::new(); nf];

    // Sort factors at y=0 by match key; track same slot at each evaluation.
    let mut order: Vec<usize> = (0..nf).collect();
    order.sort_by(|&a, &b| {
        factor_match_key(&evals[0][a], x).cmp(&factor_match_key(&evals[0][b], x))
    });
    for slot in 0..nf {
        tracks[slot].push(evals[0][order[slot]].clone());
    }

    for k in 1..evals.len() {
        let mut perm: Vec<usize> = (0..nf).collect();
        perm.sort_by(|&a, &b| {
            factor_match_key(&evals[k][a], x).cmp(&factor_match_key(&evals[k][b], x))
        });
        for slot in 0..nf {
            tracks[slot].push(evals[k][perm[slot]].clone());
        }
    }

    let mut factors = Vec::with_capacity(nf);
    for track in &tracks {
        factors.push(lift_factor_from_evals(track, x, y)?);
    }

    let prod = factors.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod != *p {
        return None;
    }
    Some(factors)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    #[test]
    fn hensel_three_linear_shifted() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x
            .sub(&y)
            .add(&Poly::one())
            .mul(&x.sub(&y))
            .mul(&x.sub(&y).sub(&Poly::one()));
        let f = try_hensel_lift_bivariate(&p, &Var::from("x"), &Var::from("y"))
            .expect("hensel lift failed");
        assert_eq!(f.len(), 3);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    fn hensel_two_bilinear_factors() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.add(&y).sub(&Poly::one()).mul(&x.add(&y).add(&Poly::one()));
        let f = try_hensel_lift_bivariate(&p, &Var::from("x"), &Var::from("y")).unwrap();
        assert_eq!(f.len(), 2);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }
}
