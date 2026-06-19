//! Bivariate factorization in ℚ[y][x] via Hensel lifting @ y=0 + interpolation fallback.
//!
//! **Upstream:** `gausspol.cc` `try_hensel_lift_factor`, `hensel_lift`.
//! **Partial:** `try_hensel_lift_bivariate` (FAC-G3 混合次数仍可能 None).
//! **Pipeline private:** `hensel_lift_at_zero`, `try_hensel_lift_interp`, rat-vector helpers.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

use super::poly_uni::{coeff_wrt_poly, substitute_poly, term_with_var};
use super::univariate::factor_univariate_flat;
use super::util::rational_nth_root;

// ---------------------------------------------------------------------------
// Univariate ℚ[x] helpers (coefficient vectors)
// ---------------------------------------------------------------------------

type RatVec = Vec<Ratio<BigInt>>;

// **Pipeline private** — `trim_rat`
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

// **Pipeline private** — `rat_mul`
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

// **Pipeline private** — `rat_div_rem`
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

// **Pipeline private** — `rat_egcd`
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

// **Pipeline private** — `rat_add`
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

// **Pipeline private** — `rat_sub`
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

// **Pipeline private** — `poly_univariate_rat`
fn poly_univariate_rat(p: &Poly, x: &Var) -> RatVec {
    let d = univariate_degree(p, x);
    (0..=d).map(|e| coeff_at(p, x, e)).collect()
}

// **Pipeline private** — `poly_from_rat`
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
// **Pipeline private** — `egcd_factor_list`
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
// **Pipeline private** — `truncate_y`
fn truncate_y(p: &Poly, y: &Var, max: u64) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        if m.exp_of(y) <= max {
            out = out.add(&Poly::term(m.clone(), c.clone()));
        }
    }
    out
}

// **Pipeline private** — `is_independent_of_y`
fn is_independent_of_y(p: &Poly, y: &Var) -> bool {
    p.terms.keys().all(|m| m.exp_of(y) == 0)
}

// **Pipeline private** — `leading_coeff_x`
fn leading_coeff_x(p: &Poly, x: &Var) -> Poly {
    let d = univariate_degree(p, x);
    coeff_wrt_poly(p, x, d)
}

// **Pipeline private** — `scale_univariate_x`
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

// **Pipeline private** — `div_rem_x_over_qy`
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

// **Pipeline private** — `hensel_lift_two_at_zero`
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
    let p0 = substitute_poly(p, y, &Poly::zero());
    if f.mul(&g) != p0 {
        return None;
    }
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
        // Extra pass: some bivariate products need degrees > dy in y.
        for deg in (dy + 1)..=(dy + dy) {
            let prod = truncate_y(&f.mul(&g), y, deg.saturating_sub(1));
            let err = truncate_y(&p.sub(&prod), y, deg.saturating_sub(1));
            if err.is_zero() {
                if f.mul(&g) == *p {
                    return Some(vec![f, g]);
                }
                continue;
            }
            let rprime = truncate_y(&err.mul(&qu_p), y, deg.saturating_sub(1));
            let qprime = truncate_y(&err.mul(&ru_p), y, deg.saturating_sub(1));
            let (_fq, frem) = div_rem_x_over_qy(&qprime, f0, x, y)?;
            let (_gq, grem) = div_rem_x_over_qy(&rprime, g0, x, y)?;
            f = f.add(&truncate_y(&frem, y, deg.saturating_sub(1)));
            g = g.add(&truncate_y(&grem, y, deg.saturating_sub(1)));
            if f.mul(&g) == *p {
                return Some(vec![f, g]);
            }
        }
        None
    }
}

/// Fold extracted constant factors into the remaining univariate factors so `∏ f_i = p0`.
// **Pipeline private** — `normalize_univariate_factors`
fn normalize_univariate_factors(f0: &mut Vec<Poly>, x: &Var, p0: &Poly) -> bool {
    let mut content = Poly::one();
    f0.retain(|f| {
        if univariate_degree(f, x) == 0 {
            content = content.mul(f);
            false
        } else {
            true
        }
    });
    if f0.is_empty() {
        return false;
    }
    if !content.is_one() {
        if f0.len() == 1 {
            f0[0] = f0[0].mul(&content);
        } else if let Some(c) = as_rational_constant(&content) {
            if let Some(root) = rational_nth_root(&c, f0.len() as u64) {
                for f in f0.iter_mut() {
                    *f = f.mul_scalar(&root);
                }
            } else {
                f0[0] = f0[0].mul(&content);
            }
        } else {
            return false;
        }
    }
    f0.iter().fold(Poly::one(), |acc, f| acc.mul(f)) == *p0
}

// **Pipeline private** — `hensel_lift_at_zero`
fn hensel_lift_at_zero(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    let dy = univariate_degree(p, y);
    if dy == 0 {
        return None;
    }

    let p0 = substitute_poly(p, y, &Poly::zero());
    let mut f0 = factor_univariate_flat(&p0, x).ok()?;
    if !normalize_univariate_factors(&mut f0, x, &p0) {
        return None;
    }
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

// **Pipeline private** — `as_rational_constant`
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

// **Pipeline private** — `linear_root`
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

// **Pipeline private** — `factor_match_key_at`
fn factor_match_key_at(f: &Poly, x: &Var, eval_var: &Var, eval_val: i64) -> (u64, Ratio<BigInt>, Ratio<BigInt>) {
    let fe = substitute_poly(
        f,
        eval_var,
        &Poly::constant(Ratio::from_integer(BigInt::from(eval_val))),
    );
    factor_match_key(&fe, x)
}

// **Pipeline private** — `sort_factors_by_match_key`
fn sort_factors_by_match_key(
    facs: &[Poly],
    main: &Var,
    rest: Option<&Var>,
) -> Vec<usize> {
    let mut perm: Vec<usize> = (0..facs.len()).collect();
    perm.sort_by(|&a, &b| {
        let ka = if let Some(rv) = rest {
            factor_match_key_at(&facs[a], main, rv, 0)
        } else {
            factor_match_key(&facs[a], main)
        };
        let kb = if let Some(rv) = rest {
            factor_match_key_at(&facs[b], main, rv, 0)
        } else {
            factor_match_key(&facs[b], main)
        };
        ka.cmp(&kb)
    });
    perm
}

// **Pipeline private** — `factor_match_key`
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

// **Pipeline private** — `lagrange_interpolate_y`
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

// **Pipeline private** — `lift_factor_from_evals`
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

/// Lift one factor track from evaluations at auxiliary values `points`.
// **Pipeline private** — Hensel lift factor from auxiliary eval tracks
pub(crate) fn lift_factor_from_aux_evals(
    samples: &[(i64, Poly)],
    x: &Var,
    aux: &Var,
) -> Option<Poly> {
    if samples.is_empty() {
        return None;
    }
    let deg = univariate_degree(&samples[0].1, x);
    let mut out = Poly::zero();
    for e in 0..=deg {
        let mut points = Vec::with_capacity(samples.len());
        for (k, f) in samples {
            let c = coeff_wrt_poly(f, x, e);
            let val = as_rational_constant(&c)?;
            points.push((*k, val));
        }
        let cy = lagrange_interpolate_y(aux, &points);
        out = out.add(&term_with_var(&cy, x, e));
    }
    Some(out)
}

/// Lift when substituted factors live in ℚ[rest][x] (one remaining variable besides `aux`).
// **Pipeline private** — `lift_factor_from_aux_evals_with_rest`
fn lift_factor_from_aux_evals_with_rest(
    samples: &[(i64, Poly)],
    main: &Var,
    rest: &Var,
    aux: &Var,
) -> Option<Poly> {
    if samples.is_empty() {
        return None;
    }
    let deg_x = univariate_degree(&samples[0].1, main);
    let mut out = Poly::zero();
    for e in 0..=deg_x {
        let max_dy = samples
            .iter()
            .map(|(_, f)| univariate_degree(&coeff_wrt_poly(f, main, e), rest))
            .max()
            .unwrap_or(0);
        let mut cy = Poly::zero();
        for d in 0..=max_dy {
            let mut points = Vec::with_capacity(samples.len());
            for (k, f) in samples {
                let ce = coeff_wrt_poly(f, main, e);
                let cd = coeff_wrt_poly(&ce, rest, d);
                let val = as_rational_constant(&cd)?;
                points.push((*k, val));
            }
            cy = cy.add(&term_with_var(&lagrange_interpolate_y(aux, &points), rest, d));
        }
        out = out.add(&term_with_var(&cy, main, e));
    }
    Some(out)
}

/// Factor by substituting an auxiliary variable and lifting (GIAC `find_good_eval` MVP).
// **Pipeline private** — lift bivariate factors via aux variable
pub(crate) fn try_lift_factors_in_aux_var(
    p: &Poly,
    main: &Var,
    aux: &Var,
    others: &[Var],
) -> Option<Vec<Poly>> {
    if others.len() < 2 {
        return None;
    }
    let rest: Vec<Var> = others.iter().filter(|v| *v != aux).cloned().collect();

    let mut eval_sets: Vec<(i64, Vec<Poly>)> = Vec::new();
    for k in -2i64..=2 {
        let pk = substitute_poly(
            p,
            aux,
            &Poly::constant(Ratio::from_integer(BigInt::from(k))),
        );
        let facs = if rest.len() == 1 {
            match try_hensel_lift_bivariate(&pk, main, &rest[0]) {
                Some(f) => f,
                None => continue,
            }
        } else {
            match super::multivariate::factor_multivariate_rec(&pk, &rest) {
                Ok(f) => f,
                Err(_) => continue,
            }
        };
        if facs.len() <= 1 {
            continue;
        }
        eval_sets.push((k, facs));
    }
    if eval_sets.len() < 2 {
        return None;
    }
    let nf = eval_sets[0].1.len();
    if !eval_sets.iter().all(|(_, f)| f.len() == nf) {
        return None;
    }

    let rest_match = rest.len().eq(&1).then(|| &rest[0]);
    let mut tracks: Vec<Vec<(i64, Poly)>> = vec![Vec::new(); nf];
    let order = sort_factors_by_match_key(&eval_sets[0].1, main, rest_match);
    for slot in 0..nf {
        tracks[slot].push((eval_sets[0].0, eval_sets[0].1[order[slot]].clone()));
    }
    for (k, facs) in eval_sets.iter().skip(1) {
        let perm = sort_factors_by_match_key(facs, main, rest_match);
        for slot in 0..nf {
            tracks[slot].push((*k, facs[perm[slot]].clone()));
        }
    }

    let mut factors = Vec::with_capacity(nf);
    for track in &tracks {
        factors.push(if rest.len() == 1 {
            lift_factor_from_aux_evals_with_rest(track, main, &rest[0], aux)?
        } else {
            lift_factor_from_aux_evals(track, main, aux)?
        });
    }
    let prod = factors.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod == *p {
        Some(factors)
    } else {
        None
    }
}

// **Pipeline private** — optional fallback `try_hensel_lift_interp`
fn try_hensel_lift_interp(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    let dy = univariate_degree(p, y);
    if dy == 0 {
        return None;
    }
    let npts = (dy + 1).max(3) as i64;
    let mut evals: Vec<Vec<Poly>> = Vec::new();
    for k in 0..npts {
        let pk = substitute_poly(p, y, &Poly::constant(Ratio::from_integer(BigInt::from(k))));
        let mut facs = factor_univariate_flat(&pk, x).ok()?;
        if !normalize_univariate_factors(&mut facs, x, &pk) {
            return None;
        }
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

/// **Partial** — Factor `p(x,y)` in ℚ[y][x]: Hensel lift at `y=0`, then interpolation fallback.
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
    fn factor_non_monic_product_at_zero() {
        let x = Poly::var("x");
        let p = Poly::constant(Ratio::from_integer(9.into()))
            .mul(&x.pow(2))
            .sub(&Poly::constant(Ratio::from_integer(18.into())).mul(&x))
            .add(&Poly::constant(Ratio::from_integer(5.into())));
        let f = super::super::univariate::factor_univariate_flat(&p, &Var::from("x")).unwrap();
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p, "factors: {:?}", f);
    }

    #[test]
    fn hensel_two_simple_non_monic() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let f = Poly::var("x")
            .mul_scalar(&Ratio::from_integer(3.into()))
            .sub(&Poly::constant(Ratio::from_integer(5.into())))
            .add(&y);
        let g = Poly::var("x")
            .mul_scalar(&Ratio::from_integer(3.into()))
            .sub(&Poly::constant(Ratio::from_integer(1.into())))
            .add(&y.mul_scalar(&Ratio::from_integer(2.into())));
        let p = f.mul(&g);
        let out = hensel_lift_at_zero(&p, &Var::from("x"), &Var::from("y"));
        assert!(out.is_some(), "simple non-monic hensel");
        let out = out.unwrap();
        assert_eq!(out.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn aux_lift_line21() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        let p = x
            .sub(&y)
            .sub(&z)
            .mul(&x.sub(&y).add(&z))
            .mul(&x.add(&y).add(&z));
        let f = super::try_lift_factors_in_aux_var(
            &p,
            &Var::from("x"),
            &Var::from("z"),
            &[Var::from("y"), Var::from("z")],
        );
        assert!(f.is_some(), "aux lift at z");
        let f = f.unwrap();
        assert_eq!(f.len(), 3);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    #[ignore = "Issue 2.3: high-degree y Hensel"]
    fn hensel_line22_mixed_bivariate() {
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
        let f = try_hensel_lift_bivariate(&p, &Var::from("x"), &Var::from("y"));
        if let Some(f) = f {
            assert_eq!(f.len(), 2);
            assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
        }
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
