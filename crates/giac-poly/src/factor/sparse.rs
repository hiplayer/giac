//! Sparse / heuristic multivariate factorization fallbacks.
//!
//! **Upstream:** `ezgcd.cc` `try_sparse_factor`, `gausspol.cc` `unitaryfactor` / `pzadic`.
//! **Partial:** `try_sparse_factor` (FAC-G1), `try_sparse_factor_bi` (FAC-G1 sum-coeff + dilation),
//! `try_heuristic_factor_bivariate` (FAC-G1/G3 last resort).

use std::collections::HashMap;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::monomial::Var;
use crate::nested::{
    CoeffRingPoly, EmbedFactorDraft, EmbedMonomial, MainVar, TnEmbed, UnivariateIn,
};
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

use super::hensel::normalize_univariate_factors;
use super::poly_uni::{coeff_wrt_poly, primitive_part_wrt, substitute_poly, term_with_var};
use super::univariate::factor_univariate_flat;

/// **Partial** — Sparse reconstruction from univariate factors at `other = 0` (FAC-G1).
///
/// Upstream `try_sparse_factor`: `lcp^(s-1)*p = ∏ P_i` with `lcp = Tfirstcoeff(p)` (Poly in
/// `other`), unknown lower coeffs matching eval factor pattern. Supports general `s≥2` via
/// iterative linear solve; optimized 2-factor path with bilinear `A·B` completion.
pub fn try_sparse_factor(p: &Poly, main: &Var, other: &Var) -> Option<Vec<Poly>> {
    try_sparse_factor_at(p, main, other, None)
        .or_else(|| {
            super::eval::find_good_eval(p, main, &[other], &[0])
                .and_then(|ge| ge.first_value().cloned())
                .and_then(|v| try_sparse_factor_at(p, main, other, Some(v)))
        })
        .or_else(|| try_sparse_factor_at(p, other, main, None))
        .or_else(|| {
            super::eval::find_good_eval(p, other, &[main], &[0])
                .and_then(|ge| ge.first_value().cloned())
                .and_then(|v| try_sparse_factor_at(p, other, main, Some(v)))
        })
}

/// **Partial** — sparse factor with optional auxiliary evaluation point (upstream `b0`).
pub fn try_sparse_factor_at(
    p: &Poly,
    main: &Var,
    other: &Var,
    at: Option<Ratio<BigInt>>,
) -> Option<Vec<Poly>> {
    try_sparse_factor_impl(p, main, other, at)
}

// **Pipeline private** — single `(main, other)` attempt via typed sparse stages
fn try_sparse_factor_impl(
    p: &Poly,
    main: &Var,
    other: &Var,
    at: Option<Ratio<BigInt>>,
) -> Option<Vec<Poly>> {
    SparseAtZero::prepare(p, main, other, at)?
        .into_system()?
        .solve()
}

/// Stage 1: univariate image @ `aux = eval` (upstream sparse @0 / good eval).
// **Pipeline private**
struct SparseAtZero<'a> {
    p: &'a Poly,
    main: MainVar,
    aux: Var,
    eval: Ratio<BigInt>,
    p0: Poly,
    facs: Vec<Poly>,
    s: usize,
    lcp: Poly,
    dy: usize,
    n_la: usize,
    tdeg: usize,
}

impl<'a> SparseAtZero<'a> {
    fn prepare(
        p: &'a Poly,
        main: &Var,
        aux: &Var,
        at: Option<Ratio<BigInt>>,
    ) -> Option<Self> {
        let dx = univariate_degree(p, main);
        if dx == 0 || univariate_degree(p, aux) == 0 {
            return None;
        }

        let eval = at.unwrap_or_else(Ratio::zero);
        let p0 = substitute_poly(p, aux, &Poly::constant(eval.clone()));
        let mut facs = factor_univariate_flat(&p0, main).ok()?;
        if !normalize_univariate_factors(&mut facs, main, &p0) {
            return None;
        }
        let s = facs.len();
        if s < 2 {
            return None;
        }

        let n_la: usize = facs
            .iter()
            .map(|f| factor_unknown_count(f, main))
            .sum();
        let tdeg = p.terms.keys().map(|m| m.degree()).max().unwrap_or(0) as usize;
        if n_la == 0 || n_la >= 5usize.max(tdeg / 2) {
            return None;
        }

        let dy = univariate_degree(p, aux) as usize;
        let lcp = leading_coeff_main(p, main);
        if lcp.is_zero() {
            return None;
        }

        Some(Self {
            p,
            main: MainVar::new(main.clone()),
            aux: aux.clone(),
            eval,
            p0,
            facs,
            s,
            lcp,
            dy,
            n_la,
            tdeg,
        })
    }

    fn into_system(self) -> Option<SparseSystem<'a>> {
        let main = self.main.as_var();
        let templates = build_factor_templates(&self.facs, main);
        let n_vars = self.n_la * (self.dy + 1);
        let target = scale_by_lcp_power(self.p, &self.lcp, self.s - 1);
        let mut eqs = build_sparse_equations(
            &templates,
            &self.lcp,
            &target,
            main,
            &self.aux,
            self.dy,
        );
        if eqs.is_empty() {
            return None;
        }
        Some(SparseSystem {
            p: self.p,
            main: self.main,
            aux: self.aux,
            lcp: self.lcp,
            s: self.s,
            templates,
            eqs,
            n_vars,
            dy: self.dy,
            n_la: self.n_la,
            target,
        })
    }
}

/// Stage 2: linear system for unknown coeffs (no further `Poly` round-trips until solve).
// **Pipeline private**
struct SparseSystem<'a> {
    p: &'a Poly,
    main: MainVar,
    aux: Var,
    lcp: Poly,
    s: usize,
    templates: Vec<FactorTemplate>,
    eqs: Vec<SparseEquation>,
    n_vars: usize,
    dy: usize,
    n_la: usize,
    target: Poly,
}

impl<'a> SparseSystem<'a> {
    fn solve(mut self) -> Option<Vec<Poly>> {
        let main = self.main.as_var();
        let sol = if self.s == 2 && self.n_la == 2 {
            solve_sparse_two_factor(
                &mut self.eqs,
                &self.lcp,
                &self.target,
                main,
                &self.aux,
                self.dy,
            )?
        } else {
            solve_sparse_system(&mut self.eqs, self.n_vars)?
        };
        let factors = build_factors(
            &self.templates,
            &sol,
            &self.lcp,
            main,
            &self.aux,
            self.dy,
        );
        let mut factors: Vec<Poly> = factors
            .into_iter()
            .filter_map(|f| primitive_part_wrt(&f, main).ok())
            .collect();
        if factors.len() != self.s {
            return None;
        }
        try_adjust_sparse_scale(&mut factors, self.p, &self.lcp, self.s, main);
        verify_sparse_factors(&factors, self.p, &self.lcp, self.s)
    }
}

// **Pipeline private** — leading coefficient of `p` in `main` (poly in remaining vars)
fn leading_coeff_main(p: &Poly, main: &Var) -> Poly {
    let main_var = MainVar::new(main.clone());
    UnivariateIn::new(p, main_var).leading_coeff()
}

// **Pipeline private** — non-leading term count in univariate factor `f`
fn factor_unknown_count(f: &Poly, main: &Var) -> usize {
    let d = univariate_degree(f, main);
    if d == 0 {
        return 0;
    }
    (0..d).filter(|e| !coeff_at(f, main, *e).is_zero()).count()
}

/// One factor template: `(main exponent, la index)`; leading term uses `lcp`.
// **Pipeline private**
struct FactorTemplate {
    terms: Vec<(u64, Option<usize>)>,
}

// **Pipeline private** — assign global `la` indices to non-leading terms
fn build_factor_templates(facs: &[Poly], main: &Var) -> Vec<FactorTemplate> {
    let mut la_next = 0usize;
    let mut out = Vec::with_capacity(facs.len());
    for f in facs {
        let d = univariate_degree(f, main);
        let mut terms = vec![(d, None)];
        for e in (0..d).rev() {
            if !coeff_at(f, main, e).is_zero() {
                terms.push((e, Some(la_next)));
                la_next += 1;
            }
        }
        out.push(FactorTemplate { terms });
    }
    out
}

// **Pipeline private** — `lcp^(pow) * p`
fn scale_by_lcp_power(p: &Poly, lcp: &Poly, pow: usize) -> Poly {
    if pow == 0 {
        return p.clone();
    }
    let mut scale = lcp.clone();
    for _ in 1..pow {
        scale = scale.mul(lcp);
    }
    p.mul(&scale)
}

/// Sparse coefficient: known part + linear part in flat unknowns.
// **Pipeline private**
#[derive(Clone, Debug, Default)]
struct SparseCoeff {
    known: Ratio<BigInt>,
    linear: Vec<(usize, Ratio<BigInt>)>,
}

impl SparseCoeff {
    fn zero() -> Self {
        Self::default()
    }

    fn is_zero(&self) -> bool {
        self.known.is_zero() && self.linear.iter().all(|(_, c)| c.is_zero())
    }

    fn add(&self, other: &Self) -> Self {
        let mut out = Self {
            known: self.known.clone() + other.known.clone(),
            linear: self.linear.clone(),
        };
        for (i, c) in &other.linear {
            merge_linear(&mut out.linear, *i, c.clone());
        }
        out
    }

    fn mul(&self, other: &Self) -> (Self, Vec<BilinearTerm>) {
        let mut out = Self::zero();
        let mut bilinear = Vec::new();

        if !self.known.is_zero() {
            if !other.known.is_zero() {
                out.known += self.known.clone() * other.known.clone();
            }
            for (i, c) in &other.linear {
                merge_linear(&mut out.linear, *i, self.known.clone() * c.clone());
            }
        }
        if !other.known.is_zero() {
            for (i, c) in &self.linear {
                merge_linear(&mut out.linear, *i, c.clone() * other.known.clone());
            }
        }
        for (i, ci) in &self.linear {
            for (j, cj) in &other.linear {
                if !ci.is_zero() && !cj.is_zero() {
                    bilinear.push(BilinearTerm {
                        i: *i,
                        j: *j,
                        coeff: ci.clone() * cj.clone(),
                    });
                }
            }
        }
        (out, bilinear)
    }
}

// **Pipeline private**
#[derive(Clone, Debug)]
struct BilinearTerm {
    i: usize,
    j: usize,
    coeff: Ratio<BigInt>,
}

// **Pipeline private**
#[derive(Clone, Debug)]
struct SparseEquation {
    known: Ratio<BigInt>,
    linear: Vec<(usize, Ratio<BigInt>)>,
    bilinear: Vec<BilinearTerm>,
}

impl SparseEquation {
    fn is_zero(&self) -> bool {
        self.known.is_zero()
            && self.linear.iter().all(|(_, c)| c.is_zero())
            && self.bilinear.iter().all(|t| t.coeff.is_zero())
    }

    fn from_parts(c: SparseCoeff, bilinear: Vec<BilinearTerm>) -> Self {
        Self {
            known: c.known,
            linear: c.linear,
            bilinear,
        }
    }

    fn is_linear(&self) -> bool {
        self.bilinear.is_empty()
    }

    fn substitute(&mut self, idx: usize, value: &Ratio<BigInt>) {
        if value.is_zero() {
            return;
        }
        for (_, c) in self.linear.iter().filter(|(i, _)| *i == idx) {
            self.known += c.clone() * value.clone();
        }
        self.linear.retain(|(i, _)| *i != idx);

        let mut new_bilinear = Vec::new();
        for t in self.bilinear.drain(..) {
            if t.i == idx {
                merge_linear(&mut self.linear, t.j, t.coeff.clone() * value.clone());
            } else if t.j == idx {
                merge_linear(&mut self.linear, t.i, t.coeff.clone() * value.clone());
            } else {
                new_bilinear.push(t);
            }
        }
        self.bilinear = new_bilinear;
    }
}

// **Pipeline private**
type SparseMap = HashMap<(u64, u64), SparseCoeff>;

// **Pipeline private**
type BilinearMap = HashMap<(u64, u64), Vec<BilinearTerm>>;

// **Pipeline private**
fn merge_linear(v: &mut Vec<(usize, Ratio<BigInt>)>, idx: usize, add: Ratio<BigInt>) {
    if add.is_zero() {
        return;
    }
    if let Some((_, c)) = v.iter_mut().find(|(i, _)| *i == idx) {
        *c += add;
        if c.is_zero() {
            v.retain(|(_, c)| !c.is_zero());
        }
    } else {
        v.push((idx, add));
    }
}

// **Pipeline private**
fn merge_bilinear(v: &mut Vec<BilinearTerm>, add: BilinearTerm) {
    if add.coeff.is_zero() {
        return;
    }
    if let Some(t) = v
        .iter_mut()
        .find(|t| t.i == add.i && t.j == add.j)
    {
        t.coeff += add.coeff.clone();
        if t.coeff.is_zero() {
            v.retain(|t| !t.coeff.is_zero());
        }
    } else {
        v.push(add);
    }
}

// **Pipeline private** — known `(main, other)` exponent map
fn poly_to_sparse_map(p: &Poly, main: &Var, other: &Var) -> SparseMap {
    let mut out = SparseMap::new();
    for (m, c) in &p.terms {
        let xe = m.exp_of(main);
        let ye = m.exp_of(other);
        let entry = out.entry((xe, ye)).or_insert_with(SparseCoeff::zero);
        entry.known += c.clone();
    }
    out
}

// **Pipeline private** — one factor template as sparse map with `la` unknowns
fn template_to_sparse_map(
    tmpl: &FactorTemplate,
    lcp: &Poly,
    other: &Var,
    dy: usize,
) -> SparseMap {
    let mut out = SparseMap::new();
    for (xe, la_idx) in &tmpl.terms {
        match la_idx {
            None => {
                for (m, c) in &lcp.terms {
                    let ye = m.exp_of(other);
                    let entry = out.entry((*xe, ye)).or_insert_with(SparseCoeff::zero);
                    entry.known += c.clone();
                }
            }
            Some(i) => {
                for k in 0..=dy {
                    let entry = out.entry((*xe, k as u64)).or_insert_with(SparseCoeff::zero);
                    let var = i * (dy + 1) + k;
                    merge_linear(&mut entry.linear, var, Ratio::one());
                }
            }
        }
    }
    out
}

// **Pipeline private** — multiply sparse maps, tracking bilinear per monomial
fn sparse_map_mul(left: &SparseMap, right: &SparseMap) -> (SparseMap, BilinearMap) {
    let mut out = SparseMap::new();
    let mut bil_map = BilinearMap::new();
    for (&(x1, y1), a) in left {
        for (&(x2, y2), b) in right {
            let (c, bil) = a.mul(b);
            let key = (x1 + x2, y1 + y2);
            if !c.is_zero() {
                let entry = out.entry(key).or_insert_with(SparseCoeff::zero);
                *entry = entry.add(&c);
            }
            if !bil.is_empty() {
                let entry = bil_map.entry(key).or_insert_with(Vec::new);
                for t in bil {
                    merge_bilinear(entry, t);
                }
            }
        }
    }
    (out, bil_map)
}

// **Pipeline private**
fn sparse_map_sub(left: &SparseMap, right: &SparseMap) -> SparseMap {
    let mut out = left.clone();
    for (k, b) in right {
        let entry = out.entry(*k).or_insert_with(SparseCoeff::zero);
        *entry = entry.add(&SparseCoeff {
            known: -b.known.clone(),
            linear: b
                .linear
                .iter()
                .map(|(i, c)| (*i, -c.clone()))
                .collect(),
        });
    }
    out.retain(|_, v| !v.is_zero());
    out
}

// **Pipeline private**
fn build_sparse_equations(
    templates: &[FactorTemplate],
    lcp: &Poly,
    target: &Poly,
    main: &Var,
    other: &Var,
    dy: usize,
) -> Vec<SparseEquation> {
    let mut prod = template_to_sparse_map(&templates[0], lcp, other, dy);
    let mut bil_map = BilinearMap::new();
    for tmpl in templates.iter().skip(1) {
        let right = template_to_sparse_map(tmpl, lcp, other, dy);
        let (p, b) = sparse_map_mul(&prod, &right);
        prod = p;
        for (k, v) in b {
            let entry = bil_map.entry(k).or_insert_with(Vec::new);
            for t in v {
                merge_bilinear(entry, t);
            }
        }
    }

    let tgt = poly_to_sparse_map(target, main, other);
    let residual = sparse_map_sub(&prod, &tgt);

    let mut keys: Vec<_> = residual.keys().chain(bil_map.keys()).copied().collect();
    keys.sort_unstable();
    keys.dedup();

    let mut eqs = Vec::new();
    for k in keys {
        let c = residual.get(&k).cloned().unwrap_or_else(SparseCoeff::zero);
        let bil = bil_map.get(&k).cloned().unwrap_or_default();
        if !c.is_zero() || !bil.is_empty() {
            eqs.push(SparseEquation::from_parts(c, bil));
        }
    }
    eqs
}

// **Pipeline private** — `x^0` coefficients of `target` by `other`-degree
fn c_from_target_x0(target: &Poly, main: &Var, other: &Var, dy: usize) -> Vec<Ratio<BigInt>> {
    let mut c = vec![Ratio::zero(); 2 * dy + 1];
    for m in 0..=2 * dy {
        let px0 = coeff_wrt_poly(target, main, 0);
        c[m] = coeff_at(&px0, other, m as u64);
    }
    c
}

// **Pipeline private** — `lcp` as univariate poly in `other`: `[lcp_0, lcp_1, …]`
fn lcp_coeffs_in_other(lcp: &Poly, other: &Var) -> Vec<Ratio<BigInt>> {
    let d = univariate_degree(lcp, other);
    (0..=d).map(|i| coeff_at(lcp, other, i)).collect()
}

// **Pipeline private** — `x^1` coefficients of `target` by `other`-degree
fn x1_target_coeffs(
    target: &Poly,
    main: &Var,
    other: &Var,
    max_m: usize,
) -> Vec<Ratio<BigInt>> {
    let px1 = coeff_wrt_poly(target, main, 1);
    (0..=max_m)
        .map(|m| coeff_at(&px1, other, m as u64))
        .collect()
}

/// Solve `∑_{i+j=m} lcp_i · s_j = t_m` for `s_j` (forward substitution).
///
/// Upstream `x^1` block: `lcp·(A+B)` convolved in `other`; then `b_j = s_j - a_j`.
// **Pipeline private**
fn solve_lcp_convolution_x1(
    lcp: &Poly,
    target: &Poly,
    main: &Var,
    other: &Var,
    dy: usize,
) -> Option<Vec<Ratio<BigInt>>> {
    let lcp_c = lcp_coeffs_in_other(lcp, other);
    if lcp_c.is_empty() || lcp_c[0].is_zero() {
        return None;
    }
    let d_lcp = lcp_c.len() - 1;
    let max_m = dy + d_lcp;
    let t = x1_target_coeffs(target, main, other, max_m);
    let mut s = vec![Ratio::zero(); dy + 1];
    for m in 0..=dy {
        let mut rhs = t.get(m).cloned().unwrap_or_else(Ratio::zero);
        for i in 1..=d_lcp {
            if m >= i {
                rhs -= lcp_c[i].clone() * s[m - i].clone();
            }
        }
        s[m] = rhs / lcp_c[0].clone();
    }
    // consistency: degrees `dy+1 … dy+d_lcp`
    for m in (dy + 1)..=max_m {
        let mut lhs = Ratio::zero();
        for i in 0..=d_lcp {
            let j = m - i;
            if j <= dy {
                lhs += lcp_c[i].clone() * s[j].clone();
            }
        }
        let rhs = t.get(m).cloned().unwrap_or_else(Ratio::zero);
        if lhs != rhs {
            return None;
        }
    }
    Some(s)
}

// **Pipeline private** — 2-factor path: `b_k = s_k - a_k`, bilinear in `a_k`
fn solve_sparse_two_factor(
    eqs: &mut [SparseEquation],
    lcp: &Poly,
    target: &Poly,
    main: &Var,
    other: &Var,
    dy: usize,
) -> Option<Vec<Ratio<BigInt>>> {
    let n = 2 * (dy + 1);
    let lk = solve_lcp_convolution_x1(lcp, target, main, other, dy)?;
    let c = c_from_target_x0(target, main, other, dy);
    let _ = eqs;

    let (r1, r2) = a0_quadratic_roots(&lk, &c)?;
    for a0 in [r1, r2] {
        let mut a = vec![Ratio::zero(); dy + 1];
        a[0] = a0;
        if bilinear_at_m(&a, &lk, 0) != c[0] {
            continue;
        }
        if complete_two_factor_a(&mut a, &lk, &c, dy).is_some() {
            let mut sol = vec![Ratio::zero(); n];
            for k in 0..=dy {
                sol[k] = a[k].clone();
                sol[dy + 1 + k] = lk[k].clone() - a[k].clone();
            }
            return Some(sol);
        }
    }
    None
}

fn a0_quadratic_roots(
    lk: &[Ratio<BigInt>],
    c: &[Ratio<BigInt>],
) -> Option<(Ratio<BigInt>, Ratio<BigInt>)> {
    let l0 = lk[0].clone();
    let c0 = c[0].clone();
    let disc = l0.clone() * l0.clone() - Ratio::from_integer(4.into()) * c0.clone();
    let r1 = if disc.is_zero() {
        l0.clone() / Ratio::from_integer(2.into())
    } else {
        let s = rational_sqrt(&disc)?;
        (l0.clone() + s.clone()) / Ratio::from_integer(2.into())
    };
    let r2 = if disc.is_zero() {
        r1.clone()
    } else {
        let s = rational_sqrt(&disc)?;
        (l0 - s) / Ratio::from_integer(2.into())
    };
    Some((r1, r2))
}

fn complete_two_factor_a(
    a: &mut [Ratio<BigInt>],
    lk: &[Ratio<BigInt>],
    c: &[Ratio<BigInt>],
    dy: usize,
) -> Option<()> {
    for m in 1..=2 * dy {
        if bilinear_at_m(a, lk, m) == c[m] {
            continue;
        }
        let mut solved = false;
        for j in 0..=dy {
            if !a[j].is_zero() {
                continue;
            }
            if let Some(v) = solve_a_j_at_m(a, lk, c, dy, m, j) {
                a[j] = v;
                solved = true;
                break;
            }
        }
        if !solved && bilinear_at_m(a, lk, m) != c[m] {
            return None;
        }
    }
    Some(())
}

fn bilinear_at_m(a: &[Ratio<BigInt>], lk: &[Ratio<BigInt>], m: usize) -> Ratio<BigInt> {
    let dy = a.len() - 1;
    let mut sum = Ratio::zero();
    for i in 0..=dy.min(m) {
        let j = m - i;
        if j > dy {
            continue;
        }
        sum += a[i].clone() * (lk[j].clone() - a[j].clone());
    }
    sum
}

fn solve_a_j_at_m(
    a: &[Ratio<BigInt>],
    lk: &[Ratio<BigInt>],
    c: &[Ratio<BigInt>],
    dy: usize,
    m: usize,
    j: usize,
) -> Option<Ratio<BigInt>> {
    // sum_{i+k=m} a_i*(l_k - a_k) with a_j unknown, others known
    let mut partial = Ratio::zero();
    let mut coeff = Ratio::zero();
    for i in 0..=dy.min(m) {
        let k = m - i;
        if k > dy {
            continue;
        }
        if k == j && i != j {
            if a[i].is_zero() {
                continue;
            }
            // a_i * (l_j - a_j): contributes -a_i * a_j
            coeff -= a[i].clone();
            partial += a[i].clone() * lk[j].clone();
        } else if i == j && k != j {
            if a[k].is_zero() && k != j {
                continue;
            }
            // a_j * (l_k - a_k): contributes a_j * (l_k - a_k)
            coeff += lk[k].clone() - a[k].clone();
        } else if i != j && k != j {
            partial += a[i].clone() * (lk[k].clone() - a[k].clone());
        } else if i == j && k == j {
            // a_j*(l_j - a_j): nonlinear, skip
        }
    }
    if coeff.is_zero() {
        return None;
    }
    Some((c[m].clone() - partial) / coeff)
}

fn rational_sqrt(d: &Ratio<BigInt>) -> Option<Ratio<BigInt>> {
    if d.is_zero() {
        return Some(Ratio::zero());
    }
    let sn = integer_sqrt_bigint(d.numer())?;
    let sd = integer_sqrt_bigint(d.denom())?;
    Some(Ratio::new(sn, sd))
}

fn integer_sqrt_bigint(n: &BigInt) -> Option<BigInt> {
    if n.is_negative() {
        return None;
    }
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    let mut lo = BigInt::zero();
    let mut hi = n.clone() + BigInt::one();
    while lo < hi {
        let mid = (&lo + &hi) / BigInt::from(2);
        let sq = &mid * &mid;
        match sq.cmp(n) {
            std::cmp::Ordering::Equal => return Some(mid),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    if (&(&lo - 1) * &(&lo - 1)) == *n {
        Some(&lo - 1)
    } else {
        None
    }
}

// **Pipeline private** — upstream iterative linear extraction + bilinear finish
fn solve_sparse_system(eqs: &mut Vec<SparseEquation>, n_vars: usize) -> Option<Vec<Ratio<BigInt>>> {
    let mut sol = vec![Ratio::zero(); n_vars];

    loop {
        let linear: Vec<SparseEquation> = eqs.iter().filter(|e| e.is_linear()).cloned().collect();
        if linear.is_empty() {
            break;
        }
        let step = solve_linear_equations(&linear, n_vars)?;
        for (idx, val) in step {
            sol[idx] = val.clone();
            for eq in eqs.iter_mut() {
                eq.substitute(idx, &val);
            }
        }
    }

    eqs.retain(|e| !e.is_zero());
    if !eqs.is_empty() {
        solve_bilinear_remaining(eqs, &mut sol, n_vars)?;
    }

    Some(sol)
}

// **Pipeline private** — extract one linear equation and solve via Gaussian elimination
fn solve_linear_equations(
    eqs: &[SparseEquation],
    n_vars: usize,
) -> Option<HashMap<usize, Ratio<BigInt>>> {
    if eqs.is_empty() {
        return None;
    }
    let mut matrix: Vec<Vec<Ratio<BigInt>>> = Vec::new();
    let mut rhs: Vec<Ratio<BigInt>> = Vec::new();
    for eq in eqs {
        let mut row = vec![Ratio::zero(); n_vars];
        for (i, c) in &eq.linear {
            if *i < n_vars {
                row[*i] += c.clone();
            }
        }
        if row.iter().all(|c| c.is_zero()) {
            if !eq.known.is_zero() {
                return None;
            }
            continue;
        }
        matrix.push(row);
        rhs.push(-eq.known.clone());
    }
    if matrix.is_empty() {
        return None;
    }
    let x = gaussian_elimination(&matrix, &rhs)?;
    let mut out = HashMap::new();
    for (i, v) in x.into_iter().enumerate() {
        if !v.is_zero() {
            out.insert(i, v);
        }
    }
    if out.is_empty() { None } else { Some(out) }
}

// **Pipeline private** — handle remaining bilinear equations (2-factor `A*B` terms)
fn solve_bilinear_remaining(
    eqs: &mut Vec<SparseEquation>,
    sol: &mut [Ratio<BigInt>],
    n_vars: usize,
) -> Option<()> {
    for eq in eqs.iter_mut() {
        if eq.bilinear.len() == 1 && eq.linear.is_empty() {
            let t = &eq.bilinear[0];
            if eq.known.is_zero() && t.coeff.is_one() {
                // v_i * v_j + known = 0  with known in substituted eq
                let vi = sol[t.i].clone();
                if t.j < n_vars && sol[t.j].is_zero() && !vi.is_zero() {
                    sol[t.j] = -eq.known.clone() / vi.clone();
                    eq.known = Ratio::zero();
                }
            }
        }
    }
    eqs.retain(|e| !e.is_zero());
    if eqs.is_empty() {
        return Some(());
    }

    // Pair bilinear: substitute from linear relations in same equation
    for eq in eqs.iter_mut() {
        if eq.bilinear.is_empty() {
            if !eq.known.is_zero() {
                return None;
            }
            continue;
        }
        if eq.bilinear.len() == 1 && eq.linear.len() == 1 {
            let t = eq.bilinear[0].clone();
            let (li, lc) = eq.linear[0].clone();
            // lc * v_li + t.coeff * v_i * v_j + known = 0
            if li == t.i && sol[t.j].is_zero() {
                // lc * v_li + t.coeff * v_li * v_j = -known
                // v_j = (-known/lc - v_li) / v_li  when v_li != 0 — skip, use quadratic
            }
        }
    }

    // Full bilinear system: expand and solve via substitution from partial sol
    let mut matrix: Vec<Vec<Ratio<BigInt>>> = Vec::new();
    let mut rhs: Vec<Ratio<BigInt>> = Vec::new();
    for eq in eqs.iter() {
        if !eq.bilinear.is_empty() {
            // For each bilinear term v_i v_j, expand using known v if one is set
            let mut row = vec![Ratio::zero(); n_vars];
            for t in &eq.bilinear {
                if !sol[t.i].is_zero() {
                    row[t.j] += t.coeff.clone() * sol[t.i].clone();
                } else if !sol[t.j].is_zero() {
                    row[t.i] += t.coeff.clone() * sol[t.j].clone();
                } else {
                    return None;
                }
            }
            for (i, v) in eq.linear.iter() {
                row[*i] += v.clone();
            }
            if row.iter().any(|r| !r.is_zero()) {
                matrix.push(row);
                rhs.push(-eq.known.clone());
            } else if !eq.known.is_zero() {
                return None;
            }
            continue;
        }
        let mut row = vec![Ratio::zero(); n_vars];
        for (i, v) in &eq.linear {
            row[*i] += v.clone();
        }
        if row.iter().any(|r| !r.is_zero()) {
            matrix.push(row);
            rhs.push(-eq.known.clone());
        } else if !eq.known.is_zero() {
            return None;
        }
    }
    if matrix.is_empty() {
        return Some(());
    }
    let x = gaussian_elimination(&matrix, &rhs)?;
    for (i, v) in x.into_iter().enumerate() {
        if !v.is_zero() {
            sol[i] = v;
        }
    }
    Some(())
}

// **Pipeline private** — Gaussian elimination over ℚ
fn gaussian_elimination(
    a: &[Vec<Ratio<BigInt>>],
    b: &[Ratio<BigInt>],
) -> Option<Vec<Ratio<BigInt>>> {
    if a.is_empty() {
        return Some(Vec::new());
    }
    let n = a[0].len();
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();
    let mut pivot_row = 0usize;
    let mut pivot_col = 0usize;
    let mut pivot_where: Vec<Option<usize>> = vec![None; n];

    while pivot_row < mat.len() && pivot_col < n {
        let mut sel = None;
        for r in pivot_row..mat.len() {
            if !mat[r][pivot_col].is_zero() {
                sel = Some(r);
                break;
            }
        }
        let Some(r) = sel else {
            pivot_col += 1;
            continue;
        };
        mat.swap(pivot_row, r);
        rhs.swap(pivot_row, r);
        let pivot = mat[pivot_row][pivot_col].clone();
        for c in 0..n {
            mat[pivot_row][c] = mat[pivot_row][c].clone() / pivot.clone();
        }
        rhs[pivot_row] = rhs[pivot_row].clone() / pivot;
        for r in 0..mat.len() {
            if r == pivot_row {
                continue;
            }
            let f = mat[r][pivot_col].clone();
            if f.is_zero() {
                continue;
            }
            for c in 0..n {
                let sub = mat[pivot_row][c].clone() * f.clone();
                mat[r][c] = mat[r][c].clone() - sub;
            }
            rhs[r] = rhs[r].clone() - rhs[pivot_row].clone() * f;
        }
        pivot_where[pivot_col] = Some(pivot_row);
        pivot_row += 1;
        pivot_col += 1;
    }

    let mut x = vec![Ratio::zero(); n];
    for (col, prow) in pivot_where.iter().enumerate() {
        if let Some(r) = prow {
            x[col] = rhs[*r].clone();
        }
    }
    Some(x)
}

// **Pipeline private**
fn la_poly_from_sol(la_idx: usize, dy: usize, sol: &[Ratio<BigInt>], other: &Var) -> Poly {
    let mut out = Poly::zero();
    for k in 0..=dy {
        let v = sol[la_idx * (dy + 1) + k].clone();
        if !v.is_zero() {
            out = out.add(&term_with_var(&Poly::constant(v), other, k as u64));
        }
    }
    out
}

// **Pipeline private**
fn build_factors(
    templates: &[FactorTemplate],
    sol: &[Ratio<BigInt>],
    lcp: &Poly,
    main: &Var,
    other: &Var,
    dy: usize,
) -> Vec<Poly> {
    templates
        .iter()
        .map(|tmpl| {
            let mut out = Poly::zero();
            for (xe, la_idx) in &tmpl.terms {
                let coeff = match la_idx {
                    None => lcp.clone(),
                    Some(i) => la_poly_from_sol(*i, dy, sol, other),
                };
                out = out.add(&term_with_var(&coeff, main, *xe));
            }
            out
        })
        .collect()
}

// **Pipeline private** — divide each `main`-coefficient by `den` when exact
fn divide_poly_coeffs_by(p: &Poly, den: &Poly, main: &Var) -> Option<Poly> {
    let d = univariate_degree(p, main);
    let mut out = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt_poly(p, main, e);
        if c.is_zero() {
            continue;
        }
        let q = crate::subresultant::div_exact_coeff(&c, den)?;
        out = out.add(&term_with_var(&q, main, e));
    }
    Some(out)
}

/// When template solve yields `∏ f_i = lcp^(s-1)·p`, strip one `lcp` factor (upstream scaling).
// **Pipeline private**
fn try_adjust_sparse_scale(
    factors: &mut Vec<Poly>,
    p: &Poly,
    lcp: &Poly,
    s: usize,
    main: &Var,
) {
    if factors.iter().fold(Poly::one(), |acc, f| acc.mul(f)) == *p {
        return;
    }
    if s < 2 {
        return;
    }
    let mut adj = p.clone();
    for _ in 1..s {
        adj = adj.mul(lcp);
    }
    if factors.iter().fold(Poly::one(), |acc, f| acc.mul(f)) != adj {
        return;
    }
    for idx in 0..factors.len() {
        if let Some(q) = divide_poly_coeffs_by(&factors[idx], lcp, main) {
            let mut trial = factors.clone();
            trial[idx] = q;
            let trial: Vec<Poly> = trial
                .into_iter()
                .filter_map(|f| primitive_part_wrt(&f, main).ok())
                .collect();
            if trial.len() == factors.len()
                && trial.iter().fold(Poly::one(), |acc, f| acc.mul(f)) == *p
            {
                *factors = trial;
                return;
            }
        }
    }
}

// **Pipeline private** — accept only `∏ f_i = p` (upstream `divbylgcd` applied above)
fn verify_sparse_factors(factors: &[Poly], p: &Poly, _lcp: &Poly, _s: usize) -> Option<Vec<Poly>> {
    let prod = factors.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod == *p {
        return Some(factors.to_vec());
    }
    None
}

// ---------------------------------------------------------------------------
// Multivariate sparse embedding (upstream `ezgcd.cc` `try_sparse_factor_bi`)
// ---------------------------------------------------------------------------

/// **Partial** — sparse factor via bivariate `eval_tn` embedding (FAC-G1, 3+ vars).
pub fn try_sparse_factor_bi(p: &Poly, main: &Var, auxes: &[&Var]) -> Option<Vec<Poly>> {
    if auxes.len() < 2 {
        return None;
    }
    if auxes.len() == 2 {
        return try_sparse_factor_bi_two_aux(p, main, auxes[0], auxes[1]);
    }
    None
}

// **Pipeline private** — build `Poly` from sorted embed monomials + aux exponents
fn embed_sorted_monomials(p: &Poly, main: &Var, t: &Var) -> Vec<EmbedMonomial> {
    EmbedMonomial::sorted_from_poly(p, main, t)
}

// **Pipeline private** — distinct `main`-degrees with pairwise distinct coeffs (upstream `x_degrees`).
fn bivariate_x_degrees_ok(p: &Poly, main: &Var) -> Option<Vec<u64>> {
    let main_var = MainVar::new(main.clone());
    let view = UnivariateIn::new(p, main_var);
    let mut degs = Vec::new();
    let mut prev: Option<u64> = None;
    let mut coeffs: Vec<Poly> = Vec::new();
    let d = view.degree();
    for e in (0..=d).rev() {
        let c = view.coeff_at(e);
        if c.is_zero() {
            continue;
        }
        if prev == Some(e) {
            return None;
        }
        if coeffs.iter().any(|x| x == &c) {
            return None;
        }
        degs.push(e);
        coeffs.push(c);
        prev = Some(e);
    }
    if degs.is_empty() { None } else { Some(degs) }
}

// **Pipeline private** — pick sparsest bivariate factor candidate.
fn select_bivariate_factor(
    facs: &[Poly],
    main: &Var,
    lcpt: &Poly,
) -> Option<(Poly, Vec<u64>)> {
    let mut best: Option<(Poly, Vec<u64>, usize)> = None;
    for f in facs {
        if f.is_one() {
            continue;
        }
        let degs = bivariate_x_degrees_ok(f, main).unwrap_or_default();
        let lc = leading_coeff_main(f, main);
        if lc.is_zero() {
            continue;
        }
        let multby = CoeffRingPoly::new(&lcpt)
            .exact_quo(&CoeffRingPoly::new(&lc))?;
        let score = multby.terms.len();
        if best.as_ref().map(|(_, _, s)| score < *s).unwrap_or(true) {
            best = Some((multby.mul(f), degs, score));
        }
    }
    best.map(|(p, d, _)| (p, d))
}

// **Pipeline private** — factor of `eval_tn(p)` with matching `main`-degree pattern.
fn matching_embed_factor(
    p: &Poly,
    lcp: &Poly,
    embed: &TnEmbed,
    n1: [usize; 2],
    seldegs: &[u64],
) -> Option<Poly> {
    let emb = embed.clone().with_n(n1);
    let pt_emb = emb.embed(p);
    let pt = primitive_part_wrt(pt_emb.as_poly(), emb.main.as_var()).ok()?;
    let vars = [emb.main.as_var().clone(), emb.t.as_var().clone()];
    let facs = super::multivariate::factor_multivariate_rec(&pt, &vars).ok()?;
    let lcpt = emb.embed(lcp).as_poly().clone();
    let main = emb.main.as_var();
    for f in &facs {
        if f.is_one() {
            continue;
        }
        if !seldegs.is_empty() {
            let degs = bivariate_x_degrees_ok(f, main).unwrap_or_default();
            if degs != seldegs {
                continue;
            }
        }
        let lc = leading_coeff_main(f, main);
        if lc.is_zero() {
            continue;
        }
        let Some(multby) = CoeffRingPoly::new(&lcpt).exact_quo(&CoeffRingPoly::new(&lc)) else {
            continue;
        };
        return Some(multby.mul(f));
    }
    None
}

/// Monomial-by-monomial aux exponent recovery (upstream `try_sparse_factor_bi` inner loop).
// **Pipeline private**
fn reconstruct_sparse_bi_monomials(
    draft: &mut EmbedFactorDraft,
    p: &Poly,
    lcp: &Poly,
    n: [usize; 2],
) -> Option<Poly> {
    const N_AUX: usize = 2;
    let embed = &draft.embed;
    let main = embed.main.as_var();
    let t = embed.t.as_var();
    let seldegs = draft.seldegs.clone();
    let mut i = 0usize;
    let mut increment = 1usize;

    while i < N_AUX {
        let mut n1 = n;
        n1[i] += increment;
        let curp = matching_embed_factor(p, lcp, embed, n1, &seldegs)?;
        let curp_monos = embed_sorted_monomials(&curp, main, t);

        if curp_monos.len() < draft.monos.len() {
            increment += 1;
            if increment > 3 {
                return None;
            }
            continue;
        }
        if curp_monos.len() > draft.monos.len() {
            draft.monos = curp_monos;
            draft.aux_exps = vec![(0, 0); draft.monos.len()];
            return reconstruct_sparse_bi_monomials(draft, p, lcp, n1);
        }

        let ni = n[i] as u64;
        let n1i = n1[i] as u64;
        if n1i == ni {
            return None;
        }
        for ((st, ct), exp) in draft
            .monos
            .iter()
            .zip(curp_monos.iter())
            .zip(draft.aux_exps.iter_mut())
        {
            if st.main_e != ct.main_e {
                return None;
            }
            let delta = (ct.t_e.saturating_sub(st.t_e)) / (n1i - ni);
            let mut ae = exp.0;
            let mut be = exp.1;
            if i == 0 {
                ae = delta;
            } else {
                be = delta;
            }
            if i == N_AUX - 2 {
                let mut rem = ct.t_e;
                rem = rem.saturating_sub(ae * n1[0] as u64);
                be = rem / n1[1].max(1) as u64;
            }
            *exp = (ae, be);
        }

        increment = 1;
        if i == N_AUX - 2 {
            i += 1;
        }
        i += 1;
    }

    if i < N_AUX {
        return None;
    }
    draft.materialize().map(|up| up.poly)
}

/// Dual-embedding reconstruction for monomial-only coeffs (sum-coeff uses monomial loop).
// **Pipeline private**
fn reconstruct_factor_dual_embed(
    selp: &Poly,
    p: &Poly,
    lcp: &Poly,
    embed: &TnEmbed,
    seldegs: &[u64],
) -> Option<Poly> {
    let main = embed.main.as_var();
    let t = embed.t.as_var();
    let aux_a = &embed.aux[0];
    let aux_b = &embed.aux[1];
    let curp_a = matching_embed_factor(p, lcp, embed, [2, 1], seldegs)?;
    let curp_b = matching_embed_factor(p, lcp, embed, [1, 2], seldegs)?;
    let selp_view = UnivariateIn::new(selp, embed.main.clone());
    let mut recon = Poly::zero();
    let deg_iter: Vec<u64> = if seldegs.is_empty() {
        (0..=selp_view.degree()).rev().collect()
    } else {
        seldegs.to_vec()
    };
    for e in deg_iter {
        let c0 = selp_view.coeff_at(e);
        if c0.is_zero() {
            continue;
        }
        let ca = UnivariateIn::new(&curp_a, embed.main.clone()).coeff_at(e);
        let cb = UnivariateIn::new(&curp_b, embed.main.clone()).coeff_at(e);
        let t0 = embed.view_t(&c0).degree();
        let ta = embed.view_t(&ca).degree();
        let tb = embed.view_t(&cb).degree();
        if t0 == 0 && ta == 0 && tb == 0 {
            recon = recon.add(&term_with_var(&c0, main, e));
            continue;
        }
        let b_exp = ta.saturating_sub(t0);
        let c_exp = tb.saturating_sub(t0);
        let c0t = if t0 == 0 {
            coeff_at(&c0, t, 0)
        } else {
            coeff_at(&c0, t, t0)
        };
        let mut mon = term_with_var(&Poly::constant(c0t), main, e);
        if b_exp > 0 {
            mon = mon.mul(&Poly::var(aux_a.clone()).pow(b_exp));
        }
        if c_exp > 0 {
            mon = mon.mul(&Poly::var(aux_b.clone()).pow(c_exp));
        }
        recon = recon.add(&mon);
    }
    primitive_part_wrt(&recon, main).ok().filter(|pp| !pp.is_zero())
}

// **Pipeline private**
fn reconstruct_sparse_bi_factor(
    selp: &Poly,
    seldegs: &[u64],
    p: &Poly,
    lcp: &Poly,
    embed: &TnEmbed,
    n: [usize; 2],
) -> Option<Poly> {
    if let Some(pp) = reconstruct_factor_dual_embed(selp, p, lcp, embed, seldegs) {
        if UnivariateIn::new(&pp, embed.main.clone()).divides(p) {
            return Some(pp);
        }
    }
    let mut draft = EmbedFactorDraft::from_poly(selp, embed, seldegs);
    reconstruct_sparse_bi_monomials(&mut draft, p, lcp, n)
}

// **Pipeline private** — exact quotient `p / factor` in ℚ[others][main].
fn exact_quo_wrt(p: &Poly, factor: &Poly, main: MainVar) -> Option<Poly> {
    UnivariateIn::new(factor, main).exact_quo_dividing(p).ok()
}

// **Pipeline private** — substitute `aux -> factor * aux`
fn dilate_poly(p: &Poly, aux: &Var, factor: i64) -> Poly {
    if factor == 1 {
        return p.clone();
    }
    let sub = if factor == -1 {
        Poly::var(aux.clone()).neg()
    } else {
        Poly::var(aux.clone()).mul_scalar(&Ratio::from_integer(BigInt::from(factor)))
    };
    substitute_poly(p, aux, &sub)
}

// **Pipeline private** — undo dilation: `factor*aux -> aux`
fn undilate_poly(p: &Poly, aux: &Var, factor: i64) -> Poly {
    if factor == 1 {
        return p.clone();
    }
    let inv = Ratio::from_integer(BigInt::from(1)) / Ratio::from_integer(BigInt::from(factor.abs()));
    let sub = if factor == -1 {
        Poly::var(aux.clone()).neg()
    } else {
        Poly::var(aux.clone()).mul_scalar(&inv)
    };
    substitute_poly(p, aux, &sub)
}

// **Pipeline private** — upstream random dilation fallback (deterministic seeds)
fn try_dilation_sparse_bi(
    p: &Poly,
    main: &Var,
    aux_a: &Var,
    aux_b: &Var,
) -> Option<Vec<Poly>> {
    const DILATIONS: [(i64, i64); 4] = [(2, -1), (-1, 2), (2, 2), (-1, -1)];
    for (da, db) in DILATIONS {
        let mut dilated = dilate_poly(p, aux_a, da);
        dilated = dilate_poly(&dilated, aux_b, db);
        let Some(facs) = try_sparse_factor_bi_two_aux_inner(&dilated, main, aux_a, aux_b, false) else {
            continue;
        };
        let undilated: Vec<Poly> = facs
            .iter()
            .map(|f| {
                let mut u = undilate_poly(f, aux_b, db);
                u = undilate_poly(&u, aux_a, da);
                u
            })
            .collect();
        if undilated.iter().fold(Poly::one(), |acc, q| acc.mul(q)) == *p {
            return Some(undilated);
        }
    }
    None
}

// **Pipeline private** — two auxiliary variables; upstream `n` sweep then dilation fallback.
fn try_sparse_factor_bi_two_aux(p: &Poly, main: &Var, aux_a: &Var, aux_b: &Var) -> Option<Vec<Poly>> {
    try_sparse_factor_bi_two_aux_inner(p, main, aux_a, aux_b, true)
}

fn try_sparse_factor_bi_two_aux_inner(
    p: &Poly,
    main: &Var,
    aux_a: &Var,
    aux_b: &Var,
    allow_dilation: bool,
) -> Option<Vec<Poly>> {
    let embed = TnEmbed::unit(main.clone(), aux_a.clone(), aux_b.clone());
    let lcp = leading_coeff_main(p, main);
    if lcp.is_zero() {
        return None;
    }

    let mut n = [1usize, 1usize];
    while n[0] < 4 {
        if let Some(f) = try_sparse_factor_bi_single_n(p, &embed, &lcp, n, allow_dilation) {
            return Some(f);
        }
        n[0] += 1;
    }
    if allow_dilation {
        try_dilation_sparse_bi(p, main, aux_a, aux_b)
    } else {
        None
    }
}

// **Pipeline private** — one embed exponent vector `n`
fn try_sparse_factor_bi_single_n(
    p: &Poly,
    embed: &TnEmbed,
    lcp: &Poly,
    n: [usize; 2],
    allow_dilation: bool,
) -> Option<Vec<Poly>> {
    let emb = embed.clone().with_n(n);
    let pt_emb = emb.embed(p);
    let main = emb.main.as_var();
    let pt = primitive_part_wrt(pt_emb.as_poly(), main).ok()?;
    if pt.is_one() || emb.view_t(&pt).degree() == 0 {
        return None;
    }
    let vars = [main.clone(), emb.t.as_var().clone()];
    let facs = super::multivariate::factor_multivariate_rec(&pt, &vars).ok()?;
    if facs.len() < 2 {
        return None;
    }
    let lcpt = emb.embed(lcp).as_poly().clone();

    let (selp, degs) = select_bivariate_factor(&facs, main, &lcpt)?;

    let recon = reconstruct_sparse_bi_factor(&selp, &degs, p, lcp, &emb, n)?;
    let factor = UnivariateIn::new(&recon, emb.main.clone());
    if recon.is_one() || !factor.divides(p) {
        if allow_dilation && facs.iter().any(|f| !f.is_one() && bivariate_x_degrees_ok(f, main).is_none())
        {
            return try_dilation_sparse_bi(p, main, &emb.aux[0], &emb.aux[1]);
        }
        return None;
    }
    let q = exact_quo_wrt(p, &recon, emb.main.clone())?;
    let mut out = vec![recon];
    if !q.is_one() {
        let rest = match try_sparse_factor_bi_two_aux_inner(&q, main, &emb.aux[0], &emb.aux[1], allow_dilation) {
            Some(r) if r.iter().fold(Poly::one(), |acc, f| acc.mul(f)) == q => r,
            _ => vec![q],
        };
        out.extend(rest);
    }
    if out.iter().fold(Poly::one(), |acc, f| acc.mul(f)) == *p {
        Some(out)
    } else {
        None
    }
}

/// **Partial** — Heuristic factorization via large eval + `pzadic` lift (FAC-G1/G3).
pub fn try_heuristic_factor_bivariate(p: &Poly, main: &Var, other: &Var) -> Option<Vec<Poly>> {
    let dy = univariate_degree(p, other);
    if dy == 0 || univariate_degree(p, main) == 0 {
        return None;
    }
    let norm = linfnorm(p);
    let mut base = BigInt::from(2) * norm.numer().abs() + BigInt::from(2);
    if !norm.denom().is_one() {
        base += norm.denom().abs();
    }

    for _try in 0..12 {
        let ev = substitute_poly(
            p,
            other,
            &Poly::constant(Ratio::from_integer(base.clone())),
        );
        let mut facs = factor_univariate_flat(&ev, main).ok()?;
        if !normalize_univariate_factors(&mut facs, main, &ev) {
            base += BigInt::one();
            continue;
        }
        if facs.len() <= 1 {
            base += BigInt::one();
            continue;
        }
        let lifted: Vec<Poly> = facs
            .iter()
            .map(|f| pzadic_lift(f, main, other, &base))
            .collect();
        let main_var = MainVar::new(main.clone());
        let mut rest = p.clone();
        let mut out = Vec::new();
        for f in &lifted {
            let factor = UnivariateIn::new(f, main_var.clone());
            if factor.divides(&rest) {
                let q = factor.exact_quo_dividing(&rest).ok()?;
                out.push(f.clone());
                rest = q;
            }
        }
        if (rest.is_one() || rest.is_zero()) && out.len() >= 2 {
            return Some(out);
        }
        base = base * BigInt::from(73794) / BigInt::from(27011) + BigInt::one();
    }
    None
}

// **Pipeline private** — max |coeff| of `p`
fn linfnorm(p: &Poly) -> Ratio<BigInt> {
    p.terms
        .values()
        .map(|c| c.abs())
        .max()
        .unwrap_or_else(Ratio::zero)
}

// **Pipeline private** — upstream `pzadic`: base-`n` digit expansion of coeffs → powers of `other`
fn pzadic_lift(f: &Poly, main: &Var, other: &Var, base: &BigInt) -> Poly {
    let b = base.abs();
    if b.is_zero() {
        return f.clone();
    }
    let mut out = Poly::zero();
    let d = univariate_degree(f, main);
    for e in 0..=d {
        let c = coeff_at(f, main, e);
        if c.is_zero() {
            continue;
        }
        let mut num = c.numer().clone();
        let den = c.denom().clone();
        let mut j = 0u64;
        loop {
            let denom_step = den.clone() * b.clone();
            let r = (&num % &denom_step + &denom_step) % &denom_step;
            let digit = &r / den.clone();
            if !digit.is_zero() {
                let rc = Ratio::new(digit, den.clone());
                let term = term_with_var(&Poly::constant(rc), main, e)
                    .mul(&Poly::var(other.clone()).pow(j));
                out = out.add(&term);
            }
            num = (num - r) / b.clone();
            if num.is_zero() {
                break;
            }
            j += 1;
            if j > 128 {
                break;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    fn l22_poly() -> Poly {
        let x = Poly::var("x");
        let y = Poly::var("y");
        Poly::constant(Ratio::from_integer(3.into()))
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
            )
    }

    /// Poly-`lcp` x^1 deconvolution for L22 target `lcp·p`.
    #[test]
    fn lcp_convolution_deconv_line22() {
        let p = l22_poly();
        let main = Var::from("x");
        let other = Var::from("y");
        let dy = univariate_degree(&p, &other) as usize;
        let lcp = leading_coeff_main(&p, &main);
        let target = scale_by_lcp_power(&p, &lcp, 1);
        let s = solve_lcp_convolution_x1(&lcp, &target, &main, &other, dy)
            .expect("poly-lcp deconv");
        assert_eq!(s.len(), dy + 1);
        let px1 = coeff_wrt_poly(&target, &main, 1);
        let t0 = coeff_at(&px1, &other, 0);
        assert_eq!(s[0], t0 / Ratio::from_integer(9.into()));
        // full sparse still returns None: template at `(main=x)` yields ∏f=lcp·p spurious branch
        assert!(try_sparse_factor(&p, &main, &other).is_none());
    }

    #[test]
    fn sparse_factor_tri_var() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        let p = x.add(&y.pow(2).mul(&z.pow(3))).mul(&x.add(&Poly::one()));
        let main = Var::from("x");
        let ya = Var::from("y");
        let za = Var::from("z");
        let f = try_sparse_factor_bi(&p, &main, &[&ya, &za]).expect("sparse_bi");
        assert!(f.len() >= 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn sparse_factor_tri_var_sum_coeff() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        let p = x.add(&y)
            .add(&z)
            .sub(&Poly::one())
            .mul(&x.add(&y).sub(&z).add(&Poly::one()));
        let f = try_sparse_factor_bi(
            &p,
            &Var::from("x"),
            &[&Var::from("y"), &Var::from("z")],
        )
        .expect("sparse_bi should factor sum-coeff trivariate");
        assert!(f.len() >= 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn sparse_factor_bilinear() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.add(&y).sub(&Poly::one()).mul(&x.add(&y).add(&Poly::one()));
        let f = try_sparse_factor(&p, &Var::from("x"), &Var::from("y"))
            .expect("sparse should factor bilinear");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn sparse_factor_non_monic_quadratic() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let f1 = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x)
            .sub(&Poly::constant(Ratio::from_integer(5.into())))
            .add(&y);
        let g1 = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x)
            .sub(&Poly::constant(Ratio::from_integer(1.into())))
            .add(&y.mul_scalar(&Ratio::from_integer(2.into())));
        let p = f1.mul(&g1);
        let f = try_sparse_factor(&p, &Var::from("x"), &Var::from("y"))
            .expect("sparse should factor non-monic pair");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    #[ignore = "FAC-G1: pzadic heuristic needs integer-only path + division guard"]
    fn heuristic_factors_line22() {
        let p = l22_poly();
        let f = try_heuristic_factor_bivariate(&p, &Var::from("x"), &Var::from("y"))
            .expect("heuristic should factor L22");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }
}
