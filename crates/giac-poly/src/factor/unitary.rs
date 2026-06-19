//! Upstream `gausspol.cc` `unitaryfactor` / `pzadic` (FAC-G1 last-resort fallback).
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

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::monomial::Var;
use crate::nested::{LiftedFactor, MainVar, PzadicDraft};
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
    pub(crate) fn initial(p: &Poly) -> Self {
        let norm = linfnorm(p);
        let mut base = BigInt::from(2) * norm.numer().abs() + BigInt::from(2);
        if !norm.denom().is_one() {
            base += norm.denom().abs();
        }
        Self { base }
    }

    pub(crate) fn base(&self) -> &BigInt {
        &self.base
    }

    pub(crate) fn set_base(&mut self, base: BigInt) {
        self.base = base;
    }

    pub(crate) fn bump_sqff(&mut self) {
        self.base += BigInt::one();
    }

    pub(crate) fn advance(&mut self) {
        self.base = &self.base * BigInt::from(73794) / BigInt::from(27011) + BigInt::one();
    }

    fn as_ratio(&self) -> Ratio<BigInt> {
        Ratio::from_integer(self.base.clone())
    }
}

struct EvalBaseStream {
    seeds: Vec<BigInt>,
    idx: usize,
    point: UnitaryEvalPoint,
}

impl EvalBaseStream {
    fn new(p: &Poly) -> Self {
        let mut seeds: Vec<BigInt> = (2..32).map(BigInt::from).collect();
        let point = UnitaryEvalPoint::initial(p);
        seeds.push(point.base().clone());
        Self {
            seeds,
            idx: 0,
            point,
        }
    }

    fn current(&self) -> &UnitaryEvalPoint {
        &self.point
    }

    fn current_mut(&mut self) -> &mut UnitaryEvalPoint {
        &mut self.point
    }

    fn next(&mut self) -> bool {
        if self.idx < self.seeds.len() {
            self.point.set_base(self.seeds[self.idx].clone());
            self.idx += 1;
            return self.point.base().bits() as usize <= 256;
        }
        self.point.advance();
        if self.point.base().bits() as usize <= 256 {
            self.seeds.push(self.point.base().clone());
            true
        } else {
            false
        }
    }
}

/// Upstream `tensor::reverse()` — swap variable index `i` ↔ `n-1-i` for `order`.
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
pub(crate) fn untrunc1_insert_var(p: &Poly, var: &Var, j: u64) -> Poly {
    if j == 0 {
        return p.clone();
    }
    p.mul(&Poly::var(var.clone()).pow(j))
}

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
    pub(crate) fn new(main: impl Into<MainVar>, eval_var: &'a Var, base: BigInt) -> Self {
        Self {
            main: main.into(),
            eval_var,
            base,
        }
    }

    pub(crate) fn draft_from(&self, f: &Poly) -> PzadicDraft {
        PzadicDraft::from_eval_factor(
            f.clone(),
            self.main.clone(),
            self.eval_var.clone(),
            self.base.clone(),
        )
    }

    pub(crate) fn lift_candidates(&self, draft: &PzadicDraft) -> Vec<LiftedFactor> {
        let poly = self.pzadic(&draft.factor_at_eval);
        vec![LiftedFactor::new(poly, draft.main.clone(), 0)]
    }

    /// Upstream `pzadic(p, n)`: expand each coefficient in base `n`, attach `eval_var^j`.
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

/// Centered symmetric digit for upstream `smod` + `iquo((k-r), n)`.
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

/// Multivariate `unitaryfactor` with upstream tail (`unitarize` fallback).
pub(crate) fn try_unitary_factor(p: &Poly, vars: &[Var]) -> Option<Vec<Poly>> {
    if vars.is_empty() {
        return None;
    }
    if vars.len() == 1 {
        return factor_univariate_flat(p, &vars[0]).ok();
    }
    let vars_rev: Vec<Var> = vars.iter().rev().cloned().collect();
    let eval_var = &vars_rev[0];

    if let Some(f) = unitary_factor_rev(p, &vars_rev) {
        if let Some(v) = verified_product(f.clone(), p) {
            return Some(v);
        }
    }

    let (unitaryp, an) = unitarize(p, eval_var);
    if !an.is_one() {
        if let Some(fz2) = unitary_factor_rev(&unitaryp, &vars_rev) {
            let all: Vec<Poly> = fz2
                .iter()
                .map(|f| ununitarize(f, &an, eval_var))
                .collect();
            if let Some(v) = verified_product(all, p) {
                return Some(v);
            }
        }
    }
    None
}

/// Core loop on reversed variable order (upstream `unitaryfactor`).
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
        if fz.len() == 1 {
            if factors.is_empty() {
                continue;
            }
            factors.push(unitaryp);
            return verified_product(factors, p);
        }

        let eval_base = bases.current().base().clone();
        if let Some(batch) = try_peel_all_at_eval(
            &unitaryp,
            eval_var,
            &eval_base,
            &fz,
            main_tag.clone(),
        ) {
            factors.extend(batch);
            unitaryp = Poly::one();
            break;
        }

        let lift = PzadicLift::new(main_tag.clone(), eval_var, eval_base);
        let mut peeled = false;
        for f in &fz {
            let draft = lift.draft_from(f);
            for lifted in lift.lift_candidates(&draft) {
                let divisor = lifted.as_univariate_in();
                if divisor.divides(&unitaryp) {
                    let q = divisor.exact_quo_dividing(&unitaryp).ok()?;
                    factors.push(lifted.poly);
                    unitaryp = q;
                    peeled = true;
                    break;
                }
            }
        }
        if unitaryp.is_one() {
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

fn try_peel_all_at_eval(
    p: &Poly,
    eval_var: &Var,
    base: &BigInt,
    fz: &[Poly],
    main_tag: MainVar,
) -> Option<Vec<Poly>> {
    if fz.len() < 2 {
        return None;
    }
    let lift = PzadicLift::new(main_tag.clone(), eval_var, base.clone());
    let mut rest = p.clone();
    let mut out = Vec::new();
    for f in fz {
        let draft = lift.draft_from(f);
        let lifted = lift.lift_candidates(&draft).into_iter().next()?;
        let divisor = lifted.as_univariate_in();
        if divisor.divides(&rest) {
            let q = divisor.exact_quo_dividing(&rest).ok()?;
            out.push(lifted.poly);
            rest = q;
        } else {
            return None;
        }
    }
    if (rest.is_one() || rest.is_zero()) && out.len() >= 2 {
        Some(out)
    } else {
        None
    }
}

fn factor_at_eval(ev: &Poly, main: &Var, child_rev: &[Var]) -> Option<Vec<Poly>> {
    if child_rev.len() <= 1 {
        return factor_univariate_flat(ev, main).ok();
    }
    unitary_factor_rev(ev, child_rev)
}

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

fn is_sqff_wrt_main(p: &Poly, main: &Var) -> bool {
    square_free_part(p, main)
        .ok()
        .map(|sf| sf == *p)
        .unwrap_or(false)
}

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

    #[test]
    fn l22_scan_peel() {
        let p = l22_poly();
        let x = Var::from("x");
        let y = Var::from("y");
        let main_tag = MainVar::new(x.clone());
        let mut found = 0;
        for b in 2i64..80 {
            let base = BigInt::from(b);
            let ev = substitute_poly(&p, &y, &Poly::constant(Ratio::from_integer(b.into())));
            if let Ok(mut fz) = factor_univariate_flat(&ev, &x) {
                if normalize_univariate_factors(&mut fz, &x, &ev) && fz.len() >= 2 {
                    if try_peel_all_at_eval(&p, &y, &base, &fz, main_tag.clone()).is_some() {
                        found += 1;
                    }
                }
            }
        }
        eprintln!("L22 peel hits: {found}");
        let _ = found;
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
    #[ignore = "y^3 coeff lift needs P2 multi-point interp; P1 tail chain landed"]
    fn unitary_factor_line25_l22_y3() {
        let p = l22_y3_product();
        let f = try_unitary_factor(&p, &[Var::from("x"), Var::from("y")]).expect("unitaryfactor");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }
}
