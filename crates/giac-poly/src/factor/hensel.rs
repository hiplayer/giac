//! Bivariate factorization in ℚ[y][x] via Hensel lifting (GIAC `try_hensel_lift_factor` / `hensel_lift`).

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

use super::poly_uni::{coeff_wrt_poly, substitute_poly, term_with_var};
use super::univariate::factor_univariate_flat;

// ---------------------------------------------------------------------------
// Univariate ℚ[x] helpers (coefficient vectors)
// ---------------------------------------------------------------------------

type RatVec = Vec<Ratio<BigInt>>;

fn trim_rat(c: &[Ratio<BigInt>]) -> RatVec {
    let mut out = c.to_vec();
    while out.len() > 1 && out.last().is_some_and(|v| v.is_zero()) {
        out.pop();
    }
    if out.is_empty() {
        vec![Ratio::zero()]
    } else {
        out
    }
}

fn rat_mul(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> RatVec {
    if a.is_empty() || b.is_empty() {
        return vec![Ratio::zero()];
    }
    let mut out = vec![Ratio::zero(); a.len() + b.len() - 1];
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            out[i + j] += ai.clone() * bj.clone();
        }
    }
    trim_rat(&out)
}

fn rat_div_rem(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> (RatVec, RatVec) {
    let mut r = trim_rat(a);
    let b = trim_rat(b);
    if b.len() == 1 && b[0].is_zero() {
        return (vec![Ratio::zero()], r);
    }
    if b.len() == 1 {
        let q: RatVec = r.iter().map(|c| c / &b[0]).collect();
        return (trim_rat(&q), vec![Ratio::zero()]);
    }
    let db = b.len() - 1;
    let mut q = vec![Ratio::zero(); r.len().saturating_sub(db).max(1)];
    while r.len() > db {
        let da = r.len() - 1;
        if r[da].is_zero() {
            r.pop();
            if r.is_empty() {
                r = vec![Ratio::zero()];
            }
            continue;
        }
        let coeff = r[da].clone() / b[db].clone();
        let shift = da - db;
        if shift >= q.len() {
            q.resize(shift + 1, Ratio::zero());
        }
        q[shift] += coeff.clone();
        for j in 0..=db {
            r[j + shift] -= coeff.clone() * b[j].clone();
        }
        r = trim_rat(&r);
    }
    (trim_rat(&q), r)
}

fn rat_egcd(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> (RatVec, RatVec, RatVec) {
    let mut old_r = trim_rat(a);
    let mut r = trim_rat(b);
    let mut old_s = vec![Ratio::one()];
    let mut s = vec![Ratio::zero()];
    let mut old_t = vec![Ratio::zero()];
    let mut t = vec![Ratio::one()];
    while !(r.len() == 1 && r[0].is_zero()) {
        let (q, new_r) = rat_div_rem(&old_r, &r);
        old_r = r;
        r = new_r;
        let new_s = rat_sub(&old_s, &rat_mul(&q, &s));
        old_s = s;
        s = new_s;
        let new_t = rat_sub(&old_t, &rat_mul(&q, &t));
        old_t = t;
        t = new_t;
    }
    let lc = old_r.last().cloned().unwrap_or_else(Ratio::one);
    if !lc.is_zero() && !lc.is_one() {
        let inv = Ratio::one() / lc.clone();
        old_r = old_r.iter().map(|c| c * &inv).collect();
        old_s = old_s.iter().map(|c| c * &inv).collect();
        old_t = old_t.iter().map(|c| c * &inv).collect();
    }
    (old_r, old_s, old_t)
}

fn rat_add(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> RatVec {
    let n = a.len().max(b.len());
    let mut out = vec![Ratio::zero(); n];
    for (i, c) in a.iter().enumerate() {
        out[i] += c.clone();
    }
    for (i, c) in b.iter().enumerate() {
        out[i] += c.clone();
    }
    trim_rat(&out)
}

fn rat_sub(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> RatVec {
    let n = a.len().max(b.len());
    let mut out = vec![Ratio::zero(); n];
    for (i, c) in a.iter().enumerate() {
        out[i] += c.clone();
    }
    for (i, c) in b.iter().enumerate() {
        out[i] -= c.clone();
    }
    trim_rat(&out)
}

fn poly_univariate_rat(p: &Poly, x: &Var) -> RatVec {
    let d = univariate_degree(p, x);
    (0..=d).map(|e| coeff_at(p, x, e)).collect()
}

fn poly_from_rat(x: &Var, c: &[Ratio<BigInt>]) -> Poly {
    let mut out = Poly::zero();
    for (e, coeff) in c.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        out = out.add(&term_with_var(&Poly::constant(coeff.clone()), x, e as u64));
    }
    out
}

/// `Σ u[i] * Π_{j≠i} f[j] = 1` for pairwise coprime univariate `f[i]` (GIAC `modpoly::egcd`).
fn egcd_factor_list(factors: &[Poly], x: &Var) -> Option<Vec<Poly>> {
    let n = factors.len();
    if n == 0 {
        return None;
    }
    if n == 1 {
        return Some(vec![Poly::one()]);
    }

    let f: Vec<RatVec> = factors.iter().map(|p| poly_univariate_rat(p, x)).collect();
    let mut pi: Vec<RatVec> = vec![f[n - 1].clone()];
    for k in 1..n - 1 {
        pi.push(rat_mul(&pi[k - 1], &f[n - k - 1]));
    }

    let mut u: Vec<RatVec> = Vec::with_capacity(n);
    let mut c = vec![Ratio::one()];
    for k in 0..n - 1 {
        // modpoly: a[k]*v + pi*U = g; u.push = U*c mod a[k]; c = v*c mod pi.
        let (g, v, big_u) = rat_egcd(&f[k], &pi[n - k - 2]);
        if !(g.len() == 1 && g[0].is_one()) {
            return None;
        }
        let (_, r) = rat_div_rem(&rat_mul(&big_u, &c), &f[k]);
        u.push(r);
        let (_, new_c) = rat_div_rem(&rat_mul(&v, &c), &pi[n - k - 2]);
        c = new_c;
    }
    u.push(c);
    Some(u.into_iter().map(|v| poly_from_rat(x, &v)).collect())
}

// ---------------------------------------------------------------------------
// ℚ[y][x] helpers
// ---------------------------------------------------------------------------

/// Terms whose exponent of `y` is at most `max` (GIAC `reduce_poly` with b=0: degree < deg+1).
fn truncate_y(p: &Poly, y: &Var, max: u64) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        if m.exp_of(y) <= max {
            out = out.add(&Poly::term(m.clone(), c.clone()));
        }
    }
    out
}

fn is_independent_of_y(p: &Poly, y: &Var) -> bool {
    p.terms.keys().all(|m| m.exp_of(y) == 0)
}

fn leading_coeff_x(p: &Poly, x: &Var) -> Poly {
    let d = univariate_degree(p, x);
    coeff_wrt_poly(p, x, d)
}

fn scale_univariate_x(p: &Poly, x: &Var, scale: &Ratio<BigInt>) -> Poly {
    let d = univariate_degree(p, x);
    let mut out = Poly::zero();
    for e in 0..=d {
        let c = coeff_at(p, x, e) * scale;
        if !c.is_zero() {
            out = out.add(&term_with_var(&Poly::constant(c), x, e));
        }
    }
    out
}

fn div_rem_x_over_qy(rem: &Poly, div: &Poly, x: &Var, y: &Var) -> Option<(Poly, Poly)> {
    if !is_independent_of_y(div, y) {
        return None;
    }
    let dd = univariate_degree(div, x);
    if dd == 0 {
        return None;
    }
    let lc = coeff_at(div, x, dd);
    if lc.is_zero() {
        return None;
    }

    let mut r = rem.clone();
    let mut q = Poly::zero();
    loop {
        let dr = univariate_degree(&r, x);
        if r.is_zero() || dr < dd {
            break;
        }
        let lc_r = coeff_wrt_poly(&r, x, dr);
        let quo_coeff = lc_r.mul_scalar(&(Ratio::one() / lc.clone()));
        let shift = dr - dd;
        let q_term = term_with_var(&quo_coeff, x, shift);
        q = q.add(&q_term);
        r = r.sub(&q_term.mul(div));
    }
    Some((q, r))
}

fn hensel_lift_two_at_zero(
    p: &Poly,
    x: &Var,
    y: &Var,
    f0: &Poly,
    g0: &Poly,
) -> Option<Vec<Poly>> {
    let fa = poly_univariate_rat(f0, x);
    let ga = poly_univariate_rat(g0, x);
    let (g, qu, ru) = rat_egcd(&fa, &ga);
    if !(g.len() == 1 && g[0].is_one()) {
        return None;
    }
    let qu_p = poly_from_rat(x, &qu);
    let ru_p = poly_from_rat(x, &ru);
    let dy = univariate_degree(p, y);
    let mut f = f0.clone();
    let mut g = g0.clone();
    // GIAC `hensel_lift` (linear_lift, b=0): keep aux vars with total degree < deg.
    for deg in 1..=dy {
        let prod = truncate_y(&f.mul(&g), y, deg.saturating_sub(1));
        let err = truncate_y(&p.sub(&prod), y, deg.saturating_sub(1));
        if err.is_zero() {
            if f.mul(&g) == *p {
                break;
            }
            continue;
        }
        // r' = qu*err, q' = ru*err (swapped vs naive Bezout pairing).
        let rprime = truncate_y(&err.mul(&qu_p), y, deg.saturating_sub(1));
        let qprime = truncate_y(&err.mul(&ru_p), y, deg.saturating_sub(1));
        let (_fq, frem) = div_rem_x_over_qy(&qprime, f0, x, y)?;
        let (_gq, grem) = div_rem_x_over_qy(&rprime, g0, x, y)?;
        f = f.add(&truncate_y(&frem, y, deg.saturating_sub(1)));
        g = g.add(&truncate_y(&grem, y, deg.saturating_sub(1)));
    }
    if f.mul(&g) == *p {
        Some(vec![f, g])
    } else {
        None
    }
}

/// Full Hensel lift at `y = 0` (GIAC `try_hensel_lift_factor`, `b = 0` branch).
fn hensel_lift_at_zero(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    let dy = univariate_degree(p, y);
    if dy == 0 {
        return None;
    }

    let p0 = substitute_poly(p, y, &Poly::zero());
    let mut f0 = factor_univariate_flat(&p0, x).ok()?;
    let s = f0.len();
    if s <= 1 {
        return None;
    }

    if s == 2 {
        let mut a = f0[0].clone();
        let mut b = f0[1].clone();
        if univariate_degree(&a, x) < univariate_degree(&b, x) {
            std::mem::swap(&mut a, &mut b);
        }
        return hensel_lift_two_at_zero(p, x, y, &a, &b);
    }

    let lcp = leading_coeff_x(p, x);
    if lcp.is_zero() {
        return None;
    }

    for fi in &mut f0 {
        let d = univariate_degree(fi, x);
        let lc_i = coeff_at(fi, x, d);
        if lc_i.is_zero() {
            return None;
        }
        let lcv = if is_independent_of_y(&lcp, y) {
            coeff_at(&lcp, x, univariate_degree(&lcp, x))
        } else {
            let lcp0 = substitute_poly(&lcp, y, &Poly::zero());
            coeff_at(&lcp0, x, univariate_degree(&lcp0, x))
        };
        *fi = scale_univariate_x(fi, x, &(lcv / lc_i));
    }

    let mut p_adj = p.clone();
    if !lcp.is_one() {
        let mut extra = Poly::one();
        for _ in 0..s - 1 {
            extra = extra.mul(&lcp);
        }
        p_adj = p_adj.mul(&extra);
    }

    let u = egcd_factor_list(&f0, x)?;
    let mut p_lift: Vec<Poly> = f0.clone();
    let p0_ref: Vec<Poly> = f0.clone();

    for deg in 1..=dy {
        let mut prod = p_lift[0].clone();
        for i in 1..s {
            prod = truncate_y(&prod.mul(&p_lift[i]), y, deg.saturating_sub(1));
        }
        let err = truncate_y(
            &truncate_y(&p_adj, y, deg.saturating_sub(1)).sub(&prod),
            y,
            deg.saturating_sub(1),
        );
        if err.is_zero() {
            let full = p_lift.iter().fold(Poly::one(), |acc, f| acc.mul(f));
            if full == p_adj {
                break;
            }
            continue;
        }

        for i in 0..s {
            let rem = truncate_y(&err.mul(&u[i]), y, deg.saturating_sub(1));
            let (_quo, lift) = div_rem_x_over_qy(&rem, &p0_ref[i], x, y)?;
            p_lift[i] = p_lift[i].add(&truncate_y(&lift, y, deg.saturating_sub(1)));
        }
    }

    let prod = p_lift.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod == p_adj || prod == *p {
        Some(p_lift)
    } else if !is_independent_of_y(&lcp, y) {
        None
    } else {
        let mut extra = Poly::one();
        for _ in 0..s - 1 {
            extra = extra.mul(&lcp);
        }
        if prod == p.mul(&extra) {
            Some(p_lift)
        } else {
            None
        }
    }
}

// ---------------------------------------------------------------------------
// Fallback: evaluation + Lagrange interpolation
// ---------------------------------------------------------------------------

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

fn try_hensel_lift_interp(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    let dy = univariate_degree(p, y);
    if dy == 0 {
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

/// Factor `p(x,y)` in ℚ[y][x]: Hensel lift at `y=0`, then interpolation fallback.
pub fn try_hensel_lift_bivariate(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    if univariate_degree(p, x) == 0 || univariate_degree(p, y) == 0 {
        return None;
    }
    hensel_lift_at_zero(p, x, y).or_else(|| try_hensel_lift_interp(p, x, y))
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

    #[test]
    fn hensel_egcd_two_factors() {
        let x = Poly::var("x");
        let f1 = x.sub(&Poly::one());
        let f2 = x.add(&Poly::one());
        let a = poly_univariate_rat(&f1, &Var::from("x"));
        let b = poly_univariate_rat(&f2, &Var::from("x"));
        let (g, s, t) = rat_egcd(&a, &b);
        let direct = rat_add(&rat_mul(&s, &a), &rat_mul(&t, &b));
        assert_eq!(direct, g);
        let u = egcd_factor_list(&[f1.clone(), f2.clone()], &Var::from("x")).unwrap();
        let sum = u[0].mul(&f2).add(&u[1].mul(&f1));
        assert_eq!(sum, Poly::one());
    }

    #[test]
    fn hensel_egcd_factor_list() {
        let x = Poly::var("x");
        let f1 = x.sub(&Poly::one());
        let f2 = x.add(&Poly::one());
        let f3 = x.add(&Poly::constant(Ratio::from_integer(2.into())));
        let u = egcd_factor_list(&[f1.clone(), f2.clone(), f3.clone()], &Var::from("x")).unwrap();
        assert_eq!(u.len(), 3);
        let mut sum = Poly::zero();
        let pairs = [
            (0, f2.mul(&f3)),
            (1, f1.mul(&f3)),
            (2, f1.mul(&f2)),
        ];
        for (i, cof) in pairs {
            sum = sum.add(&u[i].mul(&cof));
        }
        assert_eq!(sum, Poly::one());
    }

    #[test]
    fn div_rem_xy_by_x_minus_one() {
        let x = Var::from("x");
        let y = Var::from("y");
        let qprime = Poly::var("x").mul(&Poly::var("y"));
        let f0 = Poly::var("x").sub(&Poly::one());
        let (q, r) = div_rem_x_over_qy(&qprime, &f0, &x, &y).unwrap();
        assert_eq!(q, Poly::var("y"));
        assert_eq!(r, Poly::var("y"));
    }

    #[test]
    fn hensel_two_bilinear_at_zero() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.add(&y).sub(&Poly::one()).mul(&x.add(&y).add(&Poly::one()));
        let f = hensel_lift_at_zero(&p, &Var::from("x"), &Var::from("y")).expect("2-factor hensel");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn hensel_uses_lift_at_zero_for_linears() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x
            .sub(&y)
            .add(&Poly::one())
            .mul(&x.sub(&y))
            .mul(&x.sub(&y).sub(&Poly::one()));
        let f = hensel_lift_at_zero(&p, &Var::from("x"), &Var::from("y")).expect("full hensel at y=0");
        assert_eq!(f.len(), 3);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }
}
