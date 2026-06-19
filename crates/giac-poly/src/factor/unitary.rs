//! Upstream `gausspol.cc` `unitaryfactor` / `pzadic` (FAC-G1 last-resort fallback).
//!
//! **数学原理：** `.doc/giac-poly-factor-unitary-principles.md`（pzadic / P2a 局部窗 + monic + Lagrange）
//!
//! ## Ring / stage typing
//!
//! | Stage | Mathematical object | Rust type | Division API |
//! |-------|---------------------|-----------|--------------|
//! | Input sqff block | `p ∈ ℚ[vars]` primitive | `&Poly` + `vars_rev` | — |
//! | Reversed poly | `p̃` after `reverse()` | `Poly` | — |
//! | Eval image | `p₀ ∈ ℚ[vars′][main]` after `eval ↦ base` | `Poly` | — |
//! | Univariate factor | `f ∈ ℚ[main]` at eval point | `factor_univariate_flat` | `FlatUni::div_rem` |
//! | pzadic draft | `f` + lift metadata (`dim+1`) | `PzadicDraft` | — |
//! | Lifted factor | `f̂ ∈ ℚ[eval_var][main]` | `LiftedFactor` | `.as_univariate_in().divides` |
//! | Output | verified factor list | `FactorSet` via caller | `product_equals` |
//!
//! **Single entry:** [`try_unitary_factor`] → [`unitary_factor_rev`];
//! non-unitary remainder → `unitarize` / `ununitarize` (upstream `do_factor_hensel` L7044–7077).
//!
//! **Upstream convention:** variables in **reversed** order; `main = vars_rev.last()`,
//! `eval_var = vars_rev[0]`; substitute `eval_var ↦ base`, factor w.r.t. `main`, faithful `pzadic`,
//! peel with exact nested-ring division (never [`Poly::div_rem`] on lifted factors).
//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::monomial::Var;
use crate::nested::{LiftedFactor, MainVar, PzadicDraft, UnivariateIn};
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};
use crate::univariate::square_free_part;

use super::hensel::normalize_univariate_factors;
use super::poly_uni::{coeff_wrt_poly, substitute_poly, term_with_var};
use super::univariate::factor_univariate_flat;

/// Upstream `GCDHEU_MAXTRY` bound for eval-point search.
pub(crate) const UNITARY_MAX_TRY: usize = 48;

/// Large-integer evaluation point (`x0 = 2·‖p‖∞ + 2`, then `x0 ← x0·73794/27011`).
#[derive(Clone, Debug)]
pub(crate) struct UnitaryEvalPoint {
    base: BigInt,
}

impl UnitaryEvalPoint {
    // **Pipeline private** — GCDHEU eval base 2·‖p‖∞+2
    pub(crate) fn initial(p: &Poly) -> Self {
        let norm = linfnorm(p);
        let mut base = BigInt::from(2) * norm.numer().abs() + BigInt::from(2);
        if !norm.denom().is_one() {
            base += norm.denom().abs();
        }
        Self { base }
    }

    // **Pipeline private** — eval base accessor
    pub(crate) fn base(&self) -> &BigInt {
        &self.base
    }

    // **Pipeline private** — set eval base
    pub(crate) fn set_base(&mut self, base: BigInt) {
        self.base = base;
    }

    // **Pipeline private** — sqff micro-bump base += 1
    pub(crate) fn bump_sqff(&mut self) {
        self.base += BigInt::one();
    }

    // **Pipeline private** — upstream eval step ⌊base·73794/27011⌋+1
    pub(crate) fn advance(&mut self) {
        self.base = &self.base * BigInt::from(73794) / BigInt::from(27011) + BigInt::one();
    }

    // **Pipeline private** — eval base as Ratio<BigInt>
    fn as_ratio(&self) -> Ratio<BigInt> {
        Ratio::from_integer(self.base.clone())
    }
}

/// Upstream eval trajectory only: `x0 = 2·‖p‖∞+2`, then `x0 ← x0·73794/27011+1`.
/// Sqff micro-bumps (`x0 += 1`) happen per-try in [`unitary_factor_rev`], not here.
struct EvalBaseStream {
    point: UnitaryEvalPoint,
    started: bool,
}

impl EvalBaseStream {
    // **Pipeline private** — EvalBaseStream from poly norm
    fn new(p: &Poly) -> Self {
        Self {
            point: UnitaryEvalPoint::initial(p),
            started: false,
        }
    }

    // **Pipeline private** — current EvalBaseStream point
    fn current(&self) -> &UnitaryEvalPoint {
        &self.point
    }

    // **Pipeline private** — mutable current EvalBaseStream point
    fn current_mut(&mut self) -> &mut UnitaryEvalPoint {
        &mut self.point
    }

    /// Next outer eval base. First call keeps `initial(p)`; later calls `advance()`.
    // **Pipeline private** — advance outer eval base or stop when bits > 256
    fn next(&mut self) -> bool {
        if self.started {
            self.point.advance();
        } else {
            self.started = true;
        }
        self.point.base().bits() as usize <= 256
    }

    /// First `limit` bases on the upstream trajectory (for tests).
    #[cfg(test)]
    // **Pipeline private** — first N bases on EvalBaseStream (tests)
    fn upstream_bases(p: &Poly, limit: usize) -> Vec<BigInt> {
        let mut stream = Self::new(p);
        let mut out = Vec::with_capacity(limit);
        while out.len() < limit {
            if !stream.next() {
                break;
            }
            out.push(stream.current().base().clone());
        }
        out
    }
}

/// Upstream `tensor::reverse()` — swap variable index `i` ↔ `n-1-i` for `order`.
// **Pipeline private** — upstream tensor reverse on variable indices
pub(crate) fn reverse_var_order(p: &Poly, order: &[Var]) -> Poly {
    if order.len() < 2 {
        return p.clone();
    }
    let n = order.len();
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        let mut acc = Poly::constant(c.clone());
        for (var, exp) in m.iter() {
            if exp == 0 {
                continue;
            }
            if let Some(pos) = order.iter().position(|v| v == var) {
                let mapped = &order[n - 1 - pos];
                acc = acc.mul(&Poly::var(mapped.clone()).pow(exp));
            } else {
                acc = acc.mul(&Poly::var(var.clone()).pow(exp));
            }
        }
        out = out.add(&acc);
    }
    out
}

/// Upstream `trunc1()` — drop the first variable dimension (`eval_var` exponents).
// **Pipeline private** — drop eval_var tail (upstream trunc1)
pub(crate) fn trunc1_drop_var(p: &Poly, drop_var: &Var) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        let mut acc = Poly::constant(c.clone());
        for (var, exp) in m.iter() {
            if var == drop_var || exp == 0 {
                continue;
            }
            acc = acc.mul(&Poly::var(var.clone()).pow(exp));
        }
        out = out.add(&acc);
    }
    out
}

/// Upstream `untrunc1(j)` — prepend `drop_var^j` to each term.
// **Pipeline private** — reinsert eval_var with zero exp (upstream untrunc1)
pub(crate) fn untrunc1_insert_var(p: &Poly, var: &Var, j: u64) -> Poly {
    if j == 0 {
        return p.clone();
    }
    p.mul(&Poly::var(var.clone()).pow(j))
}

// **Pipeline private** — group terms by eval_var exponent
fn eval_coeff_groups(p: &Poly, eval_var: &Var) -> Vec<(u64, Poly)> {
    let max = univariate_degree(p, eval_var);
    let mut out = Vec::new();
    for e in 0..=max {
        let c = coeff_wrt_poly(p, eval_var, e);
        if !c.is_zero() {
            out.push((e, c));
        }
    }
    out
}

/// Upstream `unitarize` w.r.t. `eval_var` (first index after `reverse`).
// **Pipeline private** — scale to unitary leading coeff w.r.t. eval_var
pub(crate) fn unitarize(p: &Poly, eval_var: &Var) -> (Poly, Poly) {
    let groups = eval_coeff_groups(p, eval_var);
    if groups.is_empty() {
        return (p.clone(), Poly::one());
    }
    let (max_ev, an) = groups.last().unwrap();
    let max_ev = *max_ev;
    let an = an.clone();
    if an.is_one() {
        return (p.clone(), Poly::one());
    }

    let mut unitaryp = Poly::var(eval_var.clone()).pow(max_ev);
    let mut curanpow = Poly::one();
    let mut savpow = max_ev;
    for (newpow, an_1) in groups.iter().rev().skip(1) {
        let gap = savpow - *newpow;
        if gap > 0 {
            curanpow = curanpow.mul(&pow_poly(&an, gap as usize));
        }
        let piece = term_with_var(&an_1.mul(&curanpow), eval_var, *newpow);
        unitaryp = unitaryp.add(&piece);
        savpow = *newpow;
    }
    (unitaryp, an)
}

/// Upstream `ununitarize`.
// **Pipeline private** — undo unitarize scaling factor
pub(crate) fn ununitarize(unitaryp: &Poly, an: &Poly, eval_var: &Var) -> Poly {
    if an.is_one() {
        return unitaryp.clone();
    }
    let mut ppush = Poly::zero();
    for (curpow, an_1) in eval_coeff_groups(unitaryp, eval_var) {
        let scaled = an_1.mul(&pow_poly(an, curpow as usize));
        ppush = ppush.add(&term_with_var(&scaled, eval_var, curpow));
    }
    ppush.primitive_part()
}

// **Pipeline private** — integer exponentiation in Poly ring
fn pow_poly(base: &Poly, exp: usize) -> Poly {
    if exp == 0 {
        return Poly::one();
    }
    let mut acc = Poly::one();
    let mut b = base.clone();
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            acc = acc.mul(&b);
        }
        b = b.mul(&b);
        e >>= 1;
    }
    acc
}

/// `pzadic`: faithful base-`base` digit lift (`dim+1` via `eval_var` exponents).
#[derive(Clone, Debug)]
pub(crate) struct PzadicLift<'a> {
    pub main: MainVar,
    pub eval_var: &'a Var,
    pub base: BigInt,
}

impl<'a> PzadicLift<'a> {
    // **Pipeline private** — PzadicLift builder
    pub(crate) fn new(main: impl Into<MainVar>, eval_var: &'a Var, base: BigInt) -> Self {
        Self {
            main: main.into(),
            eval_var,
            base,
        }
    }

    // **Pipeline private** — build PzadicDraft from eval factor
    pub(crate) fn draft_from(&self, f: &Poly) -> PzadicDraft {
        PzadicDraft::from_eval_factor(
            f.clone(),
            self.main.clone(),
            self.eval_var.clone(),
            self.base.clone(),
        )
    }

    // **Pipeline private** — pzadic lift candidates from draft
    pub(crate) fn lift_candidates(&self, draft: &PzadicDraft) -> Vec<LiftedFactor> {
        let poly = self.pzadic(&draft.factor_at_eval);
        vec![LiftedFactor::new(poly, draft.main.clone(), 0)]
    }

    /// Upstream `pzadic(p, n)`: expand each coefficient in base `n`, attach `eval_var^j`.
    // **Pipeline private** — faithful base-B digit lift (dim+1 via eval_var)
    pub(crate) fn pzadic(&self, f: &Poly) -> Poly {
        let main = self.main.as_var();
        let b = self.base.abs();
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
                let (digit, new_num) = sym_mod_digit(&num, &den, &b);
                if !digit.is_zero() {
                    let rc = Ratio::new(digit, den.clone());
                    let term = term_with_var(&Poly::constant(rc), main, e)
                        .mul(&Poly::var(self.eval_var.clone()).pow(j));
                    out = out.add(&term);
                }
                num = new_num;
                if num.is_zero() {
                    break;
                }
                j += 1;
                if j > 256 {
                    break;
                }
            }
        }
        out
    }
}

/// Max samples for multi-point coeff interpolation (P2a).
const MULTI_EVAL_MAX_SAMPLES: usize = 8;

// **Pipeline private** — sort key (deg, lc) for eval factors
fn factor_sort_key(f: &Poly, main: &Var) -> (u64, Ratio<BigInt>) {
    let d = univariate_degree(f, main);
    (d, coeff_at(f, main, d))
}

// **Pipeline private** — stable sort eval factor slots
fn sort_eval_factors(fz: &mut [Poly], main: &Var) {
    fz.sort_by(|a, b| factor_sort_key(a, main).cmp(&factor_sort_key(b, main)));
}

/// Lagrange interpolation: points `(x_i, v_i)` → `Poly` in `eval_var`.
// **Pipeline private** — P2a Lagrange coeff in eval_var
fn lagrange_interp_coeff(
    samples: &[(BigInt, Ratio<BigInt>)],
    eval_var: &Var,
) -> Option<Poly> {
    if samples.is_empty() {
        return None;
    }
    let mut out = Poly::zero();
    let n = samples.len();
    for i in 0..n {
        let (xi, vi) = &samples[i];
        let mut basis = Poly::constant(Ratio::one());
        for j in 0..n {
            if i == j {
                continue;
            }
            let (xj, _) = &samples[j];
            let diff = xi - xj;
            if diff.is_zero() {
                return None;
            }
            let num = Poly::var(eval_var.clone())
                .sub(&Poly::constant(Ratio::from_integer(xj.clone())));
            let scale = Ratio::new(BigInt::one(), diff);
            basis = basis.mul(&num.mul_scalar(&scale));
        }
        out = out.add(&basis.mul_scalar(vi));
    }
    Some(out)
}

// **Pipeline private** — normalize factor monic w.r.t. main
fn monic_wrt_main(f: &Poly, main: &Var) -> Option<Poly> {
    let d = univariate_degree(f, main);
    if d == 0 {
        return None;
    }
    let lc = coeff_at(f, main, d);
    if lc.is_zero() {
        return None;
    }
    Some(f.mul_scalar(&lc.recip()))
}

/// Multi-point coeff lift when single-base `pzadic` fails (P2a / line 25).
///
/// Sample bases form a **local window** below `base0` (not a global `2..N` scan):
/// `[base0 - (need-1), …, base0 + tries]` so interpolation can use nearby sqff points
/// while the outer stream stays on upstream `initial` / `advance`.
// **Pipeline private** — P2a local-window multi-point coeff lift
fn lift_factor_multi_eval(
    p: &Poly,
    eval_var: &Var,
    main: &Var,
    factor_slot: usize,
    base0: &BigInt,
    max_y_deg: u64,
) -> Option<Poly> {
    let need = (max_y_deg as usize + 1).min(MULTI_EVAL_MAX_SAMPLES).max(2);
    let back = BigInt::from(need as i64 - 1);
    let mut b = if base0 > &back {
        base0 - &back
    } else {
        BigInt::from(2)
    };
    let mut samples_by_deg: Vec<Vec<(BigInt, Ratio<BigInt>)>> = Vec::new();
    let mut tries = 0usize;
    let max_tries = need + 16;
    while samples_by_deg.first().map(|s| s.len()).unwrap_or(0) < need && tries < max_tries {
        tries += 1;
        let mut ev = substitute_poly(p, eval_var, &Poly::constant(Ratio::from_integer(b.clone())));
        for _ in 0..8 {
            if is_sqff_wrt_main(&ev, main) {
                break;
            }
            b += BigInt::one();
            ev = substitute_poly(p, eval_var, &Poly::constant(Ratio::from_integer(b.clone())));
        }
        if !is_sqff_wrt_main(&ev, main) {
            b += BigInt::one();
            continue;
        }
        let Ok(mut fz) = factor_univariate_flat(&ev, main) else {
            b += BigInt::one();
            continue;
        };
        if !normalize_univariate_factors(&mut fz, main, &ev) || fz.len() < 2 {
            b += BigInt::one();
            continue;
        }
        sort_eval_factors(&mut fz, main);
        if factor_slot >= fz.len() {
            b += BigInt::one();
            continue;
        }
        let f = monic_wrt_main(&fz[factor_slot], main)?;
        let fd = univariate_degree(&f, main);
        if samples_by_deg.is_empty() {
            samples_by_deg.resize(fd as usize + 1, Vec::new());
        } else if fd as usize + 1 != samples_by_deg.len() {
            b += BigInt::one();
            continue;
        }
        for e in 0..fd {
            samples_by_deg[e as usize].push((b.clone(), coeff_at(&f, main, e)));
        }
        b += BigInt::one();
    }
    if samples_by_deg.first().map(|s| s.len()).unwrap_or(0) < need {
        return None;
    }
    let mut lifted = Poly::zero();
    let top = samples_by_deg.len().saturating_sub(1);
    for (e, samples) in samples_by_deg.iter().enumerate() {
        if e == top {
            lifted = lifted.add(&term_with_var(&Poly::one(), main, e as u64));
            continue;
        }
        if samples.iter().all(|(_, v)| v.is_zero()) {
            continue;
        }
        let c_y = lagrange_interp_coeff(samples, eval_var)?;
        if c_y.is_zero() {
            continue;
        }
        lifted = lifted.add(&term_with_var(&c_y, main, e as u64));
    }
    if lifted.is_zero() {
        None
    } else {
        Some(lifted)
    }
}

// **Pipeline private** — pzadic peel then P2a fallback + divides check
fn try_lift_and_peel(
    unitaryp: &Poly,
    p: &Poly,
    eval_var: &Var,
    main: &Var,
    main_tag: &MainVar,
    eval_base: &BigInt,
    f: &Poly,
    factor_slot: usize,
) -> Option<Poly> {
    let lift = PzadicLift::new(main_tag.clone(), eval_var, eval_base.clone());
    let draft = lift.draft_from(f);
    for lifted in lift.lift_candidates(&draft) {
        if lifted.as_univariate_in().divides(unitaryp) {
            return Some(lifted.poly);
        }
    }
    let max_y = univariate_degree(p, eval_var);
    let interp = lift_factor_multi_eval(p, eval_var, main, factor_slot, eval_base, max_y)?;
    let candidate = LiftedFactor::new(interp, main_tag.clone(), 1);
    if candidate.as_univariate_in().divides(unitaryp) {
        Some(candidate.poly)
    } else {
        None
    }
}

/// Centered symmetric digit for upstream `smod` + `iquo((k-r), n)`.
// **Pipeline private** — symmetric mod digit for pzadic expansion
fn sym_mod_digit(num: &BigInt, den: &BigInt, base: &BigInt) -> (BigInt, BigInt) {
    let step = den * base;
    let half = &step / BigInt::from(2);
    let mut r = ((num % &step) + &step) % &step;
    if r > half {
        r -= &step;
    }
    let digit = &r / den;
    let new_num = (num - &r) / base;
    (digit, new_num)
}

/// Map factors from reversed-variable workspace back to `order`.
// **Pipeline private** — inverse of upstream `tensor::reverse()` on factors
fn unreverse_factors(factors: &[Poly], order: &[Var]) -> Vec<Poly> {
    factors
        .iter()
        .map(|f| reverse_var_order(f, order))
        .collect()
}

/// Multivariate `unitaryfactor` with upstream tail (`unitarize` fallback).
// **Partial** — FAC-G1 last-resort; sparse/Hensel fallback; bounded GCDHEU eval stream
pub(crate) fn try_unitary_factor(p: &Poly, vars: &[Var]) -> Option<Vec<Poly>> {
    if vars.is_empty() {
        return None;
    }
    if vars.len() == 1 {
        return factor_univariate_flat(p, &vars[0]).ok();
    }
    let vars_rev: Vec<Var> = vars.iter().rev().cloned().collect();
    let eval_var = &vars_rev[0];
    // U5: upstream `tensor::reverse()` before `unitaryfactor` — required for 3+ vars only.
    // Bivariate already uses `vars_rev`; reversing `p` too misaligns `EvalBaseStream::initial`
    // and forces full UNITARY_MAX_TRY² retries (line25/line26 gate timeout).
    let (work, map_back) = if vars.len() >= 3 {
        (reverse_var_order(p, vars), true)
    } else {
        (p.clone(), false)
    };

    if let Some(f) = unitary_factor_rev(&work, &vars_rev) {
        let mapped = if map_back {
            unreverse_factors(&f, vars)
        } else {
            f
        };
        if let Some(v) = verified_product(mapped, p) {
            return Some(v);
        }
    }

    let (unitaryp, an) = unitarize(&work, eval_var);
    if !an.is_one() {
        if let Some(fz2) = unitary_factor_rev(&unitaryp, &vars_rev) {
            let all: Vec<Poly> = fz2
                .iter()
                .map(|f| ununitarize(f, &an, eval_var))
                .collect();
            let mapped = if map_back {
                unreverse_factors(&all, vars)
            } else {
                all
            };
            if let Some(v) = verified_product(mapped, p) {
                return Some(v);
            }
        }
    }
    None
}

/// Core loop on reversed variable order (upstream `unitaryfactor`).
// **Partial** — core unitaryfactor loop on vars_rev; pzadic peel + P2a fallback
pub(crate) fn unitary_factor_rev(p: &Poly, vars_rev: &[Var]) -> Option<Vec<Poly>> {
    let main = vars_rev.last()?;
    let eval_var = &vars_rev[0];
    let dx = univariate_degree(p, main);
    if dx == 0 {
        return factor_constant_tail(p, vars_rev);
    }
    if dx == 1 {
        return Some(vec![p.clone()]);
    }
    if vars_rev.len() == 1 {
        return factor_univariate_flat(p, main).ok();
    }

    let mut unitaryp = p.primitive_part();
    if unitaryp.is_zero() {
        return None;
    }
    let main_tag = MainVar::new(main.clone());
    let mut factors = Vec::new();
    let mut bases = EvalBaseStream::new(p);
    let mut ntry = 0usize;

    while univariate_degree(&unitaryp, main) > 0 && !unitaryp.is_one() {
        ntry += 1;
        if ntry > UNITARY_MAX_TRY {
            break;
        }
        if !bases.next() {
            break;
        }

        let mut ev = substitute_poly(
            &unitaryp,
            eval_var,
            &Poly::constant(bases.current().as_ratio()),
        );
        for _ in 0..UNITARY_MAX_TRY {
            if is_sqff_wrt_main(&ev, main) {
                break;
            }
            bases.current_mut().bump_sqff();
            ev = substitute_poly(
                &unitaryp,
                eval_var,
                &Poly::constant(bases.current().as_ratio()),
            );
        }
        if !is_sqff_wrt_main(&ev, main) {
            continue;
        }

        let child_rev = &vars_rev[1..];
        let mut fz = factor_at_eval(&ev, main, child_rev)?;
        if !normalize_univariate_factors(&mut fz, main, &ev) {
            continue;
        }
        if fz.is_empty() {
            break;
        }
        sort_eval_factors(&mut fz, main);
        if fz.len() == 1 {
            if factors.is_empty() {
                continue;
            }
            factors.push(unitaryp);
            return verified_product(factors, p);
        }
        if univariate_degree(&unitaryp, main) == 1 {
            factors.push(unitaryp);
            return verified_product(factors, p);
        }

        let eval_base = bases.current().base().clone();
        if let Some(batch) = try_peel_all_at_eval(
            &unitaryp,
            p,
            eval_var,
            main,
            &eval_base,
            &fz,
            main_tag.clone(),
        ) {
            factors.extend(batch);
            unitaryp = Poly::one();
            break;
        }

        let mut peeled = false;
        for (fi, f) in fz.iter().enumerate() {
            if let Some(lifted_poly) =
                try_lift_and_peel(&unitaryp, p, eval_var, main, &main_tag, &eval_base, f, fi)
            {
                let divisor = UnivariateIn::new(&lifted_poly, main_tag.clone());
                let q = divisor.exact_quo_dividing(&unitaryp).ok()?;
                factors.push(lifted_poly);
                unitaryp = q;
                peeled = true;
                break;
            }
        }
        if unitaryp.is_one() {
            return verified_product(factors, p);
        }
        if univariate_degree(&unitaryp, main) == 1 {
            factors.push(unitaryp);
            return verified_product(factors, p);
        }
        if !peeled {
            continue;
        }
    }

    if unitaryp.is_one() {
        return verified_product(factors, p);
    }
    if univariate_degree(&unitaryp, main) == 0 {
        return factor_constant_tail_into(factors, &unitaryp, vars_rev, p);
    }
    if !unitaryp.is_one() {
        factors.push(unitaryp);
    }
    verified_product(factors, p)
}

// **Pipeline private** — recurse constant tail into factor list
fn factor_constant_tail_into(
    mut factors: Vec<Poly>,
    unitaryp: &Poly,
    vars_rev: &[Var],
    orig: &Poly,
) -> Option<Vec<Poly>> {
    if unitaryp.is_one() {
        return verified_product(factors, orig);
    }
    let eval_var = &vars_rev[0];
    let tmp = trunc1_drop_var(unitaryp, eval_var);
    let child_rev = &vars_rev[1..];
    if child_rev.is_empty() {
        factors.push(unitaryp.clone());
        return verified_product(factors, orig);
    }
    if let Some(tail) = unitary_factor_rev(&tmp, child_rev) {
        for f in tail {
            factors.push(untrunc1_insert_var(&f, eval_var, 0));
        }
    } else {
        factors.push(unitaryp.clone());
    }
    verified_product(factors, orig)
}

// **Pipeline private** — factor tail when main degree → 0
fn factor_constant_tail(p: &Poly, vars_rev: &[Var]) -> Option<Vec<Poly>> {
    if p.is_one() {
        return Some(vec![]);
    }
    if vars_rev.len() <= 1 {
        return Some(vec![p.clone()]);
    }
    let eval_var = &vars_rev[0];
    let tmp = trunc1_drop_var(p, eval_var);
    let child_rev = &vars_rev[1..];
    let mut out = unitary_factor_rev(&tmp, child_rev)?;
    for f in &mut out {
        *f = untrunc1_insert_var(f, eval_var, 0);
    }
    if out.is_empty() {
        Some(vec![p.clone()])
    } else {
        verified_product(out.clone(), p).or(Some(out))
    }
}

// **Pipeline private** — batch peel all slots at one eval base
fn try_peel_all_at_eval(
    unitaryp: &Poly,
    orig: &Poly,
    eval_var: &Var,
    main: &Var,
    base: &BigInt,
    fz: &[Poly],
    main_tag: MainVar,
) -> Option<Vec<Poly>> {
    if fz.len() < 2 {
        return None;
    }
    let mut rest = unitaryp.clone();
    let mut out = Vec::new();
    for (fi, f) in fz.iter().enumerate() {
        let lifted_poly = try_lift_and_peel(&rest, orig, eval_var, main, &main_tag, base, f, fi)?;
        let divisor = UnivariateIn::new(&lifted_poly, main_tag.clone());
        let q = divisor.exact_quo_dividing(&rest).ok()?;
        out.push(lifted_poly);
        rest = q;
    }
    if (rest.is_one() || rest.is_zero()) && out.len() >= 2 {
        Some(out)
    } else {
        None
    }
}

// **Pipeline private** — factor eval image w.r.t. main
fn factor_at_eval(ev: &Poly, main: &Var, child_rev: &[Var]) -> Option<Vec<Poly>> {
    if child_rev.len() <= 1 {
        return factor_univariate_flat(ev, main).ok();
    }
    unitary_factor_rev(ev, child_rev)
}

// **Pipeline private** — check factor product equals orig
fn verified_product(factors: Vec<Poly>, orig: &Poly) -> Option<Vec<Poly>> {
    if factors.len() < 2 {
        return None;
    }
    let prod = factors.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod == *orig {
        Some(factors)
    } else {
        None
    }
}

// **Pipeline private** — sqff test w.r.t. main var
fn is_sqff_wrt_main(p: &Poly, main: &Var) -> bool {
    square_free_part(p, main)
        .ok()
        .map(|sf| sf == *p)
        .unwrap_or(false)
}

// **Pipeline private** — L∞ norm of Poly coefficients
fn linfnorm(p: &Poly) -> Ratio<BigInt> {
    p.terms
        .values()
        .map(|c| c.abs())
        .max()
        .unwrap_or_else(Ratio::zero)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nested::{MainVar, UnivariateIn};
    use num_traits::One;

    fn l22_y3_product() -> Poly {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let f1 = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x)
            .sub(&y.pow(2))
            .add(&y)
            .sub(&Poly::constant(Ratio::from_integer(5.into())));
        let f2 = x
            .mul(&y)
            .add(&Poly::constant(Ratio::from_integer(3.into())).mul(&x))
            .sub(&y.pow(2))
            .sub(&Poly::one())
            .add(&y.pow(3));
        f1.mul(&f2)
    }

    fn l22_y3_x2y_product() -> Poly {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let f1 = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x)
            .sub(&y.pow(2))
            .add(&y)
            .sub(&Poly::constant(Ratio::from_integer(5.into())));
        let f2 = x
            .mul(&y)
            .add(&Poly::constant(Ratio::from_integer(3.into())).mul(&x))
            .sub(&y.pow(2))
            .sub(&Poly::one())
            .add(&y.pow(3))
            .sub(&x.pow(2).mul(&y));
        f1.mul(&f2)
    }

    #[test]
    fn pzadic_lift_linear_via_digits() {
        let x = Var::from("x");
        let y = Var::from("y");
        let base = BigInt::from(5);
        let f = Poly::var(x.clone()).sub(&Poly::constant(Ratio::from_integer(6.into())));
        let lift = PzadicLift::new(MainVar::new(x.clone()), &y, base);
        let lifted = lift.pzadic(&f);
        assert_eq!(
            lifted,
            Poly::var(x.clone())
                .sub(&Poly::var(y.clone()))
                .sub(&Poly::one())
        );
    }

    #[test]
    fn reverse_roundtrip() {
        let x = Var::from("x");
        let y = Var::from("y");
        let p = Poly::var(x.clone())
            .mul(&Poly::var(y.clone()))
            .add(&Poly::var(x.clone()).pow(2));
        let order = vec![x.clone(), y.clone()];
        let back = reverse_var_order(&reverse_var_order(&p, &order), &order);
        assert_eq!(back, p);
    }

    #[test]
    fn eval_base_stream_upstream_only() {
        let p = l22_y3_product();
        let initial = UnitaryEvalPoint::initial(&p).base().clone();
        let bases = EvalBaseStream::upstream_bases(&p, 4);
        assert_eq!(bases.len(), 4);
        assert_eq!(bases[0], initial);
        assert!(&bases[0] >= &BigInt::from(2));
        // No small-integer scan: first base is norm-derived, not 2.
        assert!(bases[0] > BigInt::from(31));
        let advanced = &initial * BigInt::from(73794) / BigInt::from(27011) + BigInt::one();
        assert_eq!(bases[1], advanced);
    }

    /// P2a on first 4 upstream bases (partial peel, with sqff bumps).
    #[test]
    #[cfg_attr(debug_assertions, ignore = "slow in debug; run `cargo test --release -p giac-poly p2a_line25_upstream_trajectory`")]
    fn p2a_line25_upstream_trajectory() {
        let p = l22_y3_product();
        let x = Var::from("x");
        let y = Var::from("y");
        let main_tag = MainVar::new(x.clone());
        let mut peel_hits = 0usize;
        for mut point in EvalBaseStream::upstream_bases(&p, 4)
            .into_iter()
            .map(|b| UnitaryEvalPoint { base: b })
        {
            let mut ev = substitute_poly(&p, &y, &Poly::constant(point.as_ratio()));
            for _ in 0..UNITARY_MAX_TRY {
                if is_sqff_wrt_main(&ev, &x) {
                    break;
                }
                point.bump_sqff();
                ev = substitute_poly(&p, &y, &Poly::constant(point.as_ratio()));
            }
            if !is_sqff_wrt_main(&ev, &x) {
                continue;
            }
            let Ok(mut fz) = factor_univariate_flat(&ev, &x) else {
                continue;
            };
            if !normalize_univariate_factors(&mut fz, &x, &ev) || fz.len() < 2 {
                continue;
            }
            sort_eval_factors(&mut fz, &x);
            let base = point.base().clone();
            if fz.iter().enumerate().any(|(fi, f)| {
                try_lift_and_peel(&p, &p, &y, &x, &main_tag, &base, f, fi).is_some()
            }) {
                peel_hits += 1;
            }
        }
        assert!(
            peel_hits >= 1,
            "P2a partial peel on upstream trajectory, got {peel_hits}"
        );
    }

    /// Local sample window below `base0` (not a global `2..N` outer scan).
    #[test]
    fn p2a_sample_window_anchors_below_base0() {
        let p = l22_y3_product();
        let x = Var::from("x");
        let y = Var::from("y");
        let base0 = UnitaryEvalPoint::initial(&p).base().clone();
        let lifted = lift_factor_multi_eval(&p, &y, &x, 1, &base0, univariate_degree(&p, &y))
            .expect("multi-eval local window");
        assert!(!lifted.is_zero());
        assert!(univariate_degree(&lifted, &x) > 0);
        assert!(base0 > BigInt::from(34));
    }

    #[test]
    fn unitary_bilinear_via_single_entry() {
        let x = Var::from("x");
        let y = Var::from("y");
        let xv = Poly::var(x.clone());
        let yv = Poly::var(y.clone());
        let p = xv.add(&yv).sub(&Poly::one()).mul(&xv.sub(&yv).add(&Poly::one()));
        let f = try_unitary_factor(&p, &[x, y]).expect("bilinear unitary");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn unitarize_extracts_leading_coeff() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x)
            .sub(&y.pow(2))
            .add(&y)
            .sub(&Poly::constant(Ratio::from_integer(5.into())));
        let (_up, an) = unitarize(&p, &Var::from("y"));
        assert_eq!(
            an,
            Poly::constant(Ratio::from_integer(BigInt::from(-1)))
        );
    }

    #[test]
    fn unitary_factor_line25_l22_y3() {
        let p = l22_y3_product();
        let f = try_unitary_factor(&p, &[Var::from("x"), Var::from("y")]).expect("unitaryfactor");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn unitary_factor_line26_l22_y3_x2y() {
        let p = l22_y3_x2y_product();
        let f = try_unitary_factor(&p, &[Var::from("x"), Var::from("y")]).expect("unitaryfactor");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }
}
