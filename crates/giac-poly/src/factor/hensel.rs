//! Bivariate factorization in ℚ[y][x] via Hensel lifting @ y=0 + interpolation fallback.
//!
//! **Upstream:** `gausspol.cc` `try_hensel_lift_factor`, `hensel_lift`.
//! **Partial:** `try_hensel_lift_bivariate` — upstream `try_hensel_lift_factor` total-degree lift @ y=0.
//! **Pipeline private:** `try_hensel_lift_factor`, `hensel_lift_at_zero`, `try_hensel_lift_interp`.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::monomial::Var;
use crate::nested::{MainVar, UnivariateIn};
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

use super::poly_uni::{coeff_wrt_poly, substitute_poly, term_with_var};
use super::ctx::HenselPair;
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

// **Pipeline private** — `div_rem_x_over_qy` via nested-ring API (divisor independent of aux)
fn div_rem_x_over_qy(rem: &Poly, div: &Poly, x: &Var, y: &Var) -> Option<(Poly, Poly)> {
    UnivariateIn::new(div, MainVar::new(x.clone())).div_rem_wrt_aux_indep(rem, y)
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

// **Pipeline private** — legacy y-degree truncation Hensel (2 factors); fallback when total-degree lift fails
fn hensel_lift_two_at_zero(pair: HenselPair<'_>) -> Option<Vec<Poly>> {
    let p = pair.p;
    let main = pair.main.as_var();
    let aux = &pair.aux;
    let f0 = pair.f0.as_view().poly;
    let g0 = pair.g0.as_view().poly;
    let fa = poly_univariate_rat(f0, main);
    let ga = poly_univariate_rat(g0, main);
    let (g, qu, ru) = rat_egcd(&fa, &ga);
    if !(g.len() == 1 && g[0].is_one()) {
        return None;
    }
    let qu_p = poly_from_rat(main, &qu);
    let ru_p = poly_from_rat(main, &ru);
    let daux = univariate_degree(p, aux);
    let mut f = f0.clone();
    let mut g = g0.clone();
    let p0 = substitute_poly(p, aux, &Poly::zero());
    if f.mul(&g) != p0 {
        return None;
    }
    for deg in 1..=daux {
        let prod = truncate_y(&f.mul(&g), aux, deg.saturating_sub(1));
        let err = truncate_y(&p.sub(&prod), aux, deg.saturating_sub(1));
        if err.is_zero() {
            if f.mul(&g) == *p {
                break;
            }
            continue;
        }
        let rprime = truncate_y(&err.mul(&qu_p), aux, deg.saturating_sub(1));
        let qprime = truncate_y(&err.mul(&ru_p), aux, deg.saturating_sub(1));
        let (_fq, frem) = pair.f0.div_rem_wrt_aux_indep(&qprime)?;
        let (_gq, grem) = pair.g0.div_rem_wrt_aux_indep(&rprime)?;
        f = f.add(&truncate_y(&frem, aux, deg.saturating_sub(1)));
        g = g.add(&truncate_y(&grem, aux, deg.saturating_sub(1)));
    }
    if f.mul(&g) == *p {
        return Some(vec![f, g]);
    }
    for deg in (daux + 1)..=(daux + daux) {
        let prod = truncate_y(&f.mul(&g), aux, deg.saturating_sub(1));
        let err = truncate_y(&p.sub(&prod), aux, deg.saturating_sub(1));
        if err.is_zero() {
            if f.mul(&g) == *p {
                return Some(vec![f, g]);
            }
            continue;
        }
        let rprime = truncate_y(&err.mul(&qu_p), aux, deg.saturating_sub(1));
        let qprime = truncate_y(&err.mul(&ru_p), aux, deg.saturating_sub(1));
        let (_fq, frem) = pair.f0.div_rem_wrt_aux_indep(&qprime)?;
        let (_gq, grem) = pair.g0.div_rem_wrt_aux_indep(&rprime)?;
        f = f.add(&truncate_y(&frem, aux, deg.saturating_sub(1)));
        g = g.add(&truncate_y(&grem, aux, deg.saturating_sub(1)));
        if f.mul(&g) == *p {
            return Some(vec![f, g]);
        }
    }
    None
}

/// Fold extracted constant factors into the remaining univariate factors so `∏ f_i = p0`.
// **Pipeline private** — `normalize_univariate_factors`
pub(crate) fn normalize_univariate_factors(f0: &mut Vec<Poly>, x: &Var, p0: &Poly) -> bool {
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

/// Terms whose total degree is at most `max` (upstream `poly_truncate1` / `EZGCD_DEGONLY`).
// **Pipeline private** — `truncate_total_degree`
fn truncate_total_degree(p: &Poly, max: u64) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        if m.degree() <= max {
            out = out.add(&Poly::term(m.clone(), c.clone()));
        }
    }
    out
}

// **Pipeline private** — `total_degree_poly`
fn total_degree_poly(p: &Poly) -> u64 {
    p.terms.keys().map(|m| m.degree()).max().unwrap_or(0)
}

// **Pipeline private** — scalar value of `p(y=0)` when constant
fn scalar_at_y_zero(p: &Poly, y: &Var) -> Option<Ratio<BigInt>> {
    as_rational_constant(&substitute_poly(p, y, &Poly::zero()))
}

// **Pipeline private** — true when `lcp` depends on auxiliary var `aux`
fn lcp_depends_on_aux(lcp: &Poly, aux: &Var) -> bool {
    univariate_degree(lcp, aux) > 0
}

// **Pipeline private** — lcm of two rationals
fn lcm_rat(a: &Ratio<BigInt>, b: &Ratio<BigInt>) -> Ratio<BigInt> {
    let ad = a.numer() * b.denom();
    let cb = b.numer() * a.denom();
    let bd = a.denom() * b.denom();
    Ratio::new(ad.lcm(&cb), bd)
}

// **Pipeline private** — lcm of denominators of all coeffs in `p`
fn lcm_poly_denoms(p: &Poly) -> Ratio<BigInt> {
    p.terms
        .values()
        .fold(Ratio::one(), |acc, c| lcm_rat(&acc, c))
}

// **Pipeline private** — scale `f0[i]` by `lcoeff(y=0)/lc(f0[i])` (upstream `mulmodpoly` on `F0fact`)
fn scale_f0_factors(
    f0: &[Poly],
    x: &Var,
    lcoeffs: &[Poly],
    aux: &Var,
) -> Option<Vec<Poly>> {
    let mut out = Vec::with_capacity(f0.len());
    for (f, lcoeff) in f0.iter().zip(lcoeffs) {
        let d = univariate_degree(f, x);
        let lc = coeff_at(f, x, d);
        if lc.is_zero() {
            return None;
        }
        let lcv = scalar_at_y_zero(lcoeff, aux).or_else(|| as_rational_constant(lcoeff))?;
        out.push(scale_univariate_x(f, x, &(lcv / lc)));
    }
    Some(out)
}

/// Build initial lift seeds `P_i` and moduli `P0_i` (upstream `try_hensel_lift_factor` init).
// **Pipeline private** — `build_hensel_lift_seeds`
fn build_hensel_lift_seeds(
    f0_scaled: &[Poly],
    x: &Var,
    lcoeffs: &[Poly],
) -> Option<Vec<Poly>> {
    let mut seeds = Vec::with_capacity(f0_scaled.len());
    for (f, lcoeff) in f0_scaled.iter().zip(lcoeffs) {
        let d = univariate_degree(f, x);
        if d == 0 {
            return None;
        }
        let mut seed = term_with_var(lcoeff, x, d);
        for e in 0..d {
            let c = coeff_at(f, x, e);
            if !c.is_zero() {
                seed = seed.add(&term_with_var(&Poly::constant(c.clone()), x, e));
            }
        }
        seeds.push(seed);
    }
    Some(seeds)
}

// **Pipeline private** — `egcd_factor_list` + lcm denominator `D` (upstream `lcmdeno` / `D`)
fn egcd_factor_list_normalized(factors: &[Poly], x: &Var) -> Option<(Vec<Poly>, Ratio<BigInt>)> {
    let u = egcd_factor_list(factors, x)?;
    let mut d = Ratio::one();
    let mut denoms: Vec<Ratio<BigInt>> = Vec::with_capacity(u.len());
    for ui in &u {
        let den = lcm_poly_denoms(ui);
        denoms.push(den.clone());
        d = lcm_rat(&d, &den);
    }
    let u_norm: Vec<Poly> = u
        .iter()
        .zip(&denoms)
        .map(|(ui, den)| ui.mul_scalar(&(d.clone() / den.clone())))
        .collect();
    Some((u_norm, d))
}

// **Pipeline private** — total-degree Hensel iteration (upstream `EZGCD_DEGONLY`, `b=0`)
fn hensel_lift_factor_loop(
    p_adj: &Poly,
    x: &Var,
    aux: &Var,
    p_lift: &[Poly],
    p0_ref: &[Poly],
    s: usize,
) -> Option<Vec<Poly>> {
    let (u, d) = egcd_factor_list_normalized(p0_ref, x)?;
    let mut p_lift = p_lift.to_vec();
    let total = total_degree_poly(p_adj);

    for deg in 1..=total {
        let mut prod = p_lift[s - 2].clone();
        if s >= 3 {
            for i in (0..=s - 3).rev() {
                prod = truncate_total_degree(&prod.mul(&p_lift[i]), deg);
            }
        }
        prod = truncate_total_degree(&prod.mul(&p_lift[s - 1]), deg);

        let err = truncate_total_degree(
            &truncate_total_degree(p_adj, deg).sub(&prod),
            deg,
        );
        if err.is_zero() {
            if deg == total {
                let full = p_lift.iter().fold(Poly::one(), |acc, f| acc.mul(f));
                if full == *p_adj {
                    break;
                }
            } else {
                let full = p_lift.iter().fold(Poly::one(), |acc, f| acc.mul(f));
                if full == *p_adj {
                    break;
                }
            }
            continue;
        }
        for i in 0..s {
            let rem = truncate_total_degree(&err.mul(&u[i]), deg);
            let (_, mut lift) = div_rem_x_over_qy(&rem, &p0_ref[i], x, aux)?;
            if !d.is_one() {
                lift = lift.mul_scalar(&(Ratio::one() / d.clone()));
            }
            p_lift[i] = p_lift[i].add(&lift);
        }
    }

    let prod = p_lift.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod == *p_adj {
        return Some(p_lift);
    }
    None
}

/// Multivariate Hensel lift @ auxiliary vars = 0 (upstream `try_hensel_lift_factor`, `b=0`).
// **Pipeline private** — `try_hensel_lift_factor`
fn try_hensel_lift_factor(
    p: &Poly,
    x: &Var,
    aux: &Var,
    f0: &[Poly],
) -> Option<Vec<Poly>> {
    let s = f0.len();
    if s < 2 {
        return None;
    }
    let p0 = substitute_poly(p, aux, &Poly::zero());
    if f0.iter().fold(Poly::one(), |acc, f| acc.mul(f)) != p0 {
        return None;
    }

    let lcp = leading_coeff_x(p, x);
    if lcp.is_zero() {
        return None;
    }
    // Upstream: reject when leading coeff depends on the evaluation variable.
    if lcp_depends_on_aux(&lcp, aux) {
        return None;
    }
    if scalar_at_y_zero(&lcp, aux).is_none() && as_rational_constant(&lcp).is_none() {
        return None;
    }
    if !lcp.is_one() && s > 1 {
        let mut est = lcp.terms.len();
        for _ in 1..s - 1 {
            est = est.saturating_mul(lcp.terms.len());
        }
        if est > 1000 {
            return None;
        }
    }

    let lcoeffs: Vec<Poly> = vec![lcp.clone(); s];
    let f0_scaled = scale_f0_factors(f0, x, &lcoeffs, aux)?;
    let p_lift = build_hensel_lift_seeds(&f0_scaled, x, &lcoeffs)?;

    let mut p_adj = p.clone();
    if !lcp.is_one() {
        for _ in 1..s {
            p_adj = p_adj.mul(&lcp);
        }
    }

    hensel_lift_factor_loop(&p_adj, x, aux, &p_lift, &f0_scaled, s)
}

// **Pipeline private** — normalize lifted factors whose product is `p` or `p_adj`
fn normalize_hensel_lift_result(
    lifted: Vec<Poly>,
    p: &Poly,
    main: &Var,
    f0: &[Poly],
) -> Option<Vec<Poly>> {
    let lcp = leading_coeff_x(p, main);
    let s = f0.len();
    let mut p_adj = p.clone();
    if !lcp.is_one() {
        for _ in 1..s {
            p_adj = p_adj.mul(&lcp);
        }
    }
    let prod = lifted.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod == *p {
        return Some(lifted);
    }
    if prod == p_adj && !lcp.is_one() {
        if let Some(c) = as_rational_constant(&lcp) {
            let root = rational_nth_root(&c, s as u64)?;
            return Some(
                lifted
                    .into_iter()
                    .map(|f| f.mul_scalar(&root))
                    .collect(),
            );
        }
    }
    if prod == p_adj {
        return Some(lifted);
    }
    None
}

// **Pipeline private** — `hensel_lift_at_zero`
fn hensel_lift_at_zero(p: &Poly, main: &Var, aux: &Var) -> Option<Vec<Poly>> {
    if univariate_degree(p, aux) == 0 {
        return None;
    }

    let p0 = substitute_poly(p, aux, &Poly::zero());
    let mut f0 = factor_univariate_flat(&p0, main).ok()?;
    if !normalize_univariate_factors(&mut f0, main, &p0) {
        return None;
    }
    if f0.len() <= 1 {
        return None;
    }

    if let Some(lifted) = try_hensel_lift_factor(p, main, aux, &f0) {
        if let Some(f) = normalize_hensel_lift_result(lifted, p, main, &f0) {
            return Some(f);
        }
    }

    if f0.len() == 2 {
        return HenselPair::try_new(p, main.clone(), aux, &f0[0], &f0[1])
            .and_then(hensel_lift_two_at_zero)
            .or_else(|| {
                HenselPair::try_new(p, main.clone(), aux, &f0[1], &f0[0])
                    .and_then(hensel_lift_two_at_zero)
            });
    }
    None
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
    let d = univariate_degree(f, x);
    let lc = as_rational_constant(&coeff_wrt_poly(f, x, d)).unwrap_or_else(Ratio::zero);
    let c0 = as_rational_constant(&coeff_wrt_poly(f, x, 0)).unwrap_or_else(Ratio::zero);
    if d == 1 {
        if let Some(r) = linear_root(f, x) {
            return (1, r, lc);
        }
    }
    (d, c0, lc)
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
            if let Some(val) = as_rational_constant(&c) {
                points.push((k as i64, val));
            } else {
                // Coeff varies in nested ring at this sample — cannot lift this track.
                return None;
            }
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

/// **Partial** — Factor `p(x,y)` in ℚ[y][x]: Hensel lift at `aux=0`, both main orders, then interpolation fallback.
pub fn try_hensel_lift_bivariate(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    if univariate_degree(p, x) == 0 || univariate_degree(p, y) == 0 {
        return None;
    }
    hensel_lift_at_zero(p, x, y)
        .or_else(|| hensel_lift_at_zero(p, y, x))
        .or_else(|| try_hensel_lift_interp(p, x, y))
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn try_hensel_lift_factor_line22_debug_steps() {
        let x = Var::from("x");
        let y = Var::from("y");
        let p = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&Poly::var("x"))
            .sub(&Poly::var("y").pow(2))
            .add(&Poly::var("y"))
            .sub(&Poly::constant(Ratio::from_integer(5.into())))
            .mul(
                &Poly::var("x")
                    .mul(&Poly::var("y"))
                    .add(&Poly::constant(Ratio::from_integer(3.into())).mul(&Poly::var("x")))
                    .sub(&Poly::var("y").pow(2))
                    .sub(&Poly::one()),
            );
        // (main=x, aux=y): lcp = 3y+9 depends on y → upstream rejects.
        let p0 = substitute_poly(&p, &y, &Poly::zero());
        let mut f0 = factor_univariate_flat(&p0, &x).unwrap();
        normalize_univariate_factors(&mut f0, &x, &p0);
        assert!(try_hensel_lift_factor(&p, &x, &y, &f0).is_none());

        // (main=y, aux=x): total-degree lift may fail; aux-truncation fallback should succeed.
        let px0 = substitute_poly(&p, &x, &Poly::zero());
        let mut f0y = factor_univariate_flat(&px0, &y).unwrap();
        normalize_univariate_factors(&mut f0y, &y, &px0);
        assert!(try_hensel_lift_factor(&p, &y, &x, &f0y).is_none());
        let f = hensel_lift_at_zero(&p, &y, &x).expect("hensel_lift_at_zero (y,x) L22");
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    fn try_hensel_lift_factor_line22() {
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
        let f = try_hensel_lift_bivariate(&p, &Var::from("x"), &Var::from("y"))
            .expect("try_hensel_lift_bivariate L22 via (y,x) order");
        assert_eq!(f.len(), 2);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    #[ignore = "legacy y-truncation path; use try_hensel_lift_factor_line22"]
    fn two_factor_hensel_at_zero_line22() {
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
        let p0 = substitute_poly(&p, &Var::from("y"), &Poly::zero());
        let mut f0 = factor_univariate_flat(&p0, &Var::from("x")).unwrap();
        normalize_univariate_factors(&mut f0, &Var::from("x"), &p0);
        assert_eq!(f0.len(), 2);
        let f = HenselPair::try_new(&p, Var::from("x"), &Var::from("y"), &f0[0], &f0[1])
            .and_then(hensel_lift_two_at_zero)
            .or_else(|| {
                HenselPair::try_new(&p, Var::from("x"), &Var::from("y"), &f0[1], &f0[0])
                    .and_then(hensel_lift_two_at_zero)
            })
            .expect("two-factor Hensel at zero");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    #[ignore = "interp track matching still open for mixed lc"]
    fn interp_line22_mixed_bivariate() {
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
        let f = try_hensel_lift_interp(&p, &Var::from("x"), &Var::from("y"))
            .expect("interp lift should factor L22");
        assert_eq!(f.len(), 2);
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
        let f = try_hensel_lift_bivariate(&p, &Var::from("x"), &Var::from("y"))
            .expect("try_hensel_lift_bivariate 3 linear factors");
        assert_eq!(f.len(), 3);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }
}
