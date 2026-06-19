//! Upstream `gausspol.cc` `unitaryfactor` / `pzadic` (FAC-G1 last-resort fallback).
//!
//! ## Ring / stage typing
//!
//! | Stage | Mathematical object | Rust type | Division API |
//! |-------|---------------------|-----------|--------------|
//! | Input sqff block | `p ∈ ℚ[vars]` primitive | `&Poly` + `vars_rev` | — |
//! | Eval image | `p₀ ∈ ℚ[vars′][main]` after `eval ↦ base` | `Poly` | — |
//! | Univariate factor | `f ∈ ℚ[main]` at eval point | `FlatUni` / `factor_univariate_flat` | `FlatUni::div_rem` |
//! | Lifted factor | `f̂ ∈ ℚ[eval_var][main]` | `PzadicLift::lift` → `UnivariateIn` | `.divides` / `.exact_quo_dividing` |
//! | Output | verified factor list | `FactorSet` via caller | `product_equals` |
//!
//! **Upstream convention:** variables in **reversed** order; `main = vars_rev.last()`,
//! `eval_var = vars_rev[0]`; substitute `eval_var ↦ base`, factor w.r.t. `main`, `pzadic` lift,
//! peel with exact nested-ring division (never [`Poly::div_rem`] on lifted factors).

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::monomial::Var;
use crate::nested::{MainVar, UnivariateIn};
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};
use crate::univariate::square_free_part;

use super::hensel::normalize_univariate_factors;
use super::poly_uni::{substitute_poly, term_with_var};
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

    /// Increment until evaluated poly is square-free in `main` (upstream gcd loop).
    pub(crate) fn bump_sqff(&mut self) {
        self.base += BigInt::one();
    }

    /// Next trial after a failed peel (upstream `iquo(x0*73794, 27011)`).
    pub(crate) fn advance(&mut self) {
        self.base = &self.base * BigInt::from(73794) / BigInt::from(27011) + BigInt::one();
    }

    fn as_ratio(&self) -> Ratio<BigInt> {
        Ratio::from_integer(self.base.clone())
    }
}

/// `pzadic`: lift `f ∈ ℚ[main]` (coeffs constant at eval) into `ℚ[eval_var][main]`.
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

    /// Expand rational coefficients in base-`base` digits; attach `eval_var^j`.
    pub(crate) fn lift(&self, f: &Poly) -> Poly {
        self.lift_factor(f)
    }

    /// Lift univariate factor at `eval_var = base` back into `ℚ[eval_var][main]`.
    pub(crate) fn lift_factor(&self, f: &Poly) -> Poly {
        self.lift_factor_candidates(f)
            .into_iter()
            .next()
            .unwrap_or_else(|| self.pzadic_digits(f))
    }

    /// Affine lifts for monic/linear factors at eval (`±eval_var` slopes).
    pub(crate) fn lift_factor_candidates(&self, f: &Poly) -> Vec<Poly> {
        let main = self.main.as_var();
        let d = univariate_degree(f, main);
        if d != 1 {
            return vec![self.pzadic_digits(f)];
        }
        let lc = coeff_at(f, main, 1);
        if lc.is_zero() {
            return vec![self.pzadic_digits(f)];
        }
        let c0 = coeff_at(f, main, 0);
        let base_r = Ratio::from_integer(self.base.clone());
        let head = term_with_var(&Poly::constant(lc.clone()), main, 1);
        vec![
            head.clone()
                .add(
                    &Poly::var(self.eval_var.clone())
                        .mul_scalar(&lc)
                        .add(&Poly::constant(c0.clone() - base_r.clone() * lc.clone())),
                ),
            head.add(
                &Poly::var(self.eval_var.clone())
                    .neg()
                    .mul_scalar(&lc)
                    .add(&Poly::constant(c0 + base_r * lc)),
            ),
            self.pzadic_digits(f),
        ]
    }

    fn pzadic_digits(&self, f: &Poly) -> Poly {
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
                let denom_step = den.clone() * b.clone();
                let r = (&num % &denom_step + &denom_step) % &denom_step;
                let digit = &r / den.clone();
                if !digit.is_zero() {
                    let rc = Ratio::new(digit, den.clone());
                    let term = term_with_var(&Poly::constant(rc), main, e)
                        .mul(&Poly::var(self.eval_var.clone()).pow(j));
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
}

/// Bivariate slice: `p ∈ ℚ[eval_var][main]` (sqff block).
pub(crate) fn try_unitary_factor_bivariate(
    p: &Poly,
    main: &Var,
    eval_var: &Var,
) -> Option<Vec<Poly>> {
    if univariate_degree(p, eval_var) == 0 || univariate_degree(p, main) == 0 {
        return None;
    }
    let mut seeds: Vec<BigInt> = small_eval_seeds(p);
    let mut point = UnitaryEvalPoint::initial(p);
    seeds.push(point.base().clone());

    let mut idx = 0usize;
    while idx < seeds.len() && idx < UNITARY_MAX_TRY {
        point.base = seeds[idx].clone();
        idx += 1;
        if point.base().bits() as usize > 256 {
            continue;
        }
        if let Some(f) = try_unitary_peel_at_base(p, main, eval_var, &point) {
            return Some(f);
        }
        point.advance();
        if point.base().bits() as usize <= 256 {
            seeds.push(point.base().clone());
        }
    }
    None
}

fn small_eval_seeds(p: &Poly) -> Vec<BigInt> {
    let mut out: Vec<BigInt> = (2..32).map(BigInt::from).collect();
    let norm = linfnorm(p);
    let mut b = BigInt::from(2) * norm.numer().abs() + BigInt::from(2);
    if !norm.denom().is_one() {
        b += norm.denom().abs();
    }
    out.push(b);
    out
}

fn try_unitary_peel_at_base(
    p: &Poly,
    main: &Var,
    eval_var: &Var,
    point: &UnitaryEvalPoint,
) -> Option<Vec<Poly>> {
    let main_tag = MainVar::new(main.clone());
    let mut ev = substitute_poly(p, eval_var, &Poly::constant(point.as_ratio()));
    for _ in 0..8 {
        if is_sqff_wrt_main(&ev, main) {
            break;
        }
        let mut bumped = point.clone();
        bumped.bump_sqff();
        ev = substitute_poly(p, eval_var, &Poly::constant(bumped.as_ratio()));
    }
    if !is_sqff_wrt_main(&ev, main) {
        return None;
    }

    let mut fz = factor_univariate_flat(&ev, main).ok()?;
    if !normalize_univariate_factors(&mut fz, main, &ev) || fz.len() < 2 {
        return None;
    }

    let mut rest = p.clone();
    let mut out = Vec::new();
    for f in &fz {
        let lift = PzadicLift::new(main_tag.clone(), eval_var, point.base().clone());
        let mut peeled_one = false;
        for lifted in lift.lift_factor_candidates(f) {
            let divisor = UnivariateIn::new(&lifted, main_tag.clone());
            if divisor.divides(&rest) {
                let q = divisor.exact_quo_dividing(&rest).ok()?;
                out.push(lifted);
                rest = q;
                peeled_one = true;
                break;
            }
        }
        if !peeled_one {
            return None;
        }
    }
    if (rest.is_one() || rest.is_zero()) && out.len() >= 2 {
        Some(out)
    } else {
        None
    }
}

/// Multivariate `unitaryfactor` on `vars` in caller order (internally reversed).
pub(crate) fn try_unitary_factor(p: &Poly, vars: &[Var]) -> Option<Vec<Poly>> {
    if vars.is_empty() {
        return None;
    }
    if vars.len() == 1 {
        return factor_univariate_flat(p, &vars[0]).ok();
    }
    let vars_rev: Vec<Var> = vars.iter().rev().cloned().collect();
    unitary_factor_rev(p, &vars_rev)
}

/// Core loop on reversed variable order (upstream `unitaryfactor`).
fn unitary_factor_rev(p: &Poly, vars_rev: &[Var]) -> Option<Vec<Poly>> {
    let main = vars_rev.last()?;
    let eval_var = &vars_rev[0];
    let dx = univariate_degree(p, main);
    if dx == 0 {
        return Some(vec![]);
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
    let mut point = UnitaryEvalPoint::initial(p);
    let mut ntry = 0usize;

    while univariate_degree(&unitaryp, main) > 0 && !unitaryp.is_one() {
        ntry += 1;
        if ntry > UNITARY_MAX_TRY {
            break;
        }

        let sub = Poly::constant(point.as_ratio());
        let mut ev = substitute_poly(&unitaryp, eval_var, &sub);
        for _ in 0..UNITARY_MAX_TRY {
            if is_sqff_wrt_main(&ev, main) {
                break;
            }
            point.bump_sqff();
            ev = substitute_poly(&unitaryp, eval_var, &Poly::constant(point.as_ratio()));
        }
        if !is_sqff_wrt_main(&ev, main) {
            point.advance();
            continue;
        }

        let child_rev = &vars_rev[1..];
        let mut fz = factor_at_eval(&ev, main, child_rev)?;
        if !normalize_univariate_factors(&mut fz, main, &ev) {
            point.advance();
            continue;
        }
        if fz.is_empty() {
            break;
        }
        if fz.len() == 1 {
            factors.push(unitaryp);
            return verified_product(factors, p);
        }

        let lift = PzadicLift::new(main_tag.clone(), eval_var, point.base.clone());
        let mut peeled = false;
        for f in &fz {
            let lifted = lift.lift(f);
            let divisor = UnivariateIn::new(&lifted, main_tag.clone());
            if divisor.divides(&unitaryp) {
                let q = divisor.exact_quo_dividing(&unitaryp).ok()?;
                factors.push(lifted);
                unitaryp = q;
                peeled = true;
            }
        }
        if unitaryp.is_one() {
            return verified_product(factors, p);
        }
        if !peeled {
            point.advance();
            continue;
        }
        point.advance();
    }

    if !unitaryp.is_one() && univariate_degree(&unitaryp, main) > 0 {
        factors.push(unitaryp);
    }
    verified_product(factors, p)
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
    fn pzadic_lift_affine_linear() {
        let x = Var::from("x");
        let y = Var::from("y");
        let base = BigInt::from(5);
        let f = Poly::var(x.clone()).sub(&Poly::constant(Ratio::from_integer(6.into())));
        let lift = PzadicLift::new(MainVar::new(x.clone()), &y, base);
        let cands = lift.lift_factor_candidates(&f);
        assert!(cands.iter().any(|u| {
            *u == Poly::var(x.clone()).sub(&Poly::var(y.clone())).sub(&Poly::one())
        }));
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
        let mut found = 0;
        for b in 2i64..80 {
            let point = UnitaryEvalPoint {
                base: BigInt::from(b),
            };
            if try_unitary_peel_at_base(&p, &x, &y, &point).is_some() {
                found += 1;
            }
        }
        eprintln!("L22 peel hits: {found}");
    }

    #[test]
    fn unitary_bilinear_linear_lift() {
        let x = Var::from("x");
        let y = Var::from("y");
        let xv = Poly::var(x.clone());
        let yv = Poly::var(y.clone());
        let p = xv.add(&yv).sub(&Poly::one()).mul(&xv.sub(&yv).add(&Poly::one()));
        let f = try_unitary_factor_bivariate(&p, &x, &y).expect("bilinear unitary");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    #[ignore = "L22+y^3 needs higher-degree coeff lift; Hensel-fail gate pending"]
    fn unitary_factor_line25_l22_y3() {
        let p = l22_y3_product();
        let f = try_unitary_factor_bivariate(&p, &Var::from("x"), &Var::from("y"))
            .expect("unitaryfactor L22+y^3");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }
}
