//! Nested-ring polynomial views: encode ℚ[others][main] in the type system.
//!
//! `Poly` is intentionally erased (one sparse representation for all rings).
//! Use [`MainVar`] + [`UnivariateIn`] when the object is univariate in `main`
//! with coefficients in ℚ[remaining vars]. Use [`CoeffRingPoly`] for coefficient-
//! ring exact division. Use [`TnEmbed`] for sparse_bi `eval_tn` maps.
//!
//! **Do not** use [`Poly::div_rem`] to test divisibility in these contexts.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_rational::Ratio;
use num_bigint::BigInt;
use num_traits::Zero;

use crate::error::PolyResult;
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::univariate_degree;
use crate::subresultant::{
    div_exact_coeff, primitive_part_wrt_impl, quo_exact_wrt, univariate_div_rem_wrt,
};

/// Main univariate indeterminate of ℚ[others][main].
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct MainVar(pub Var);

impl MainVar {
    /// **Stable** — wrap a variable name as main indeterminate.
    pub fn new(name: impl Into<Var>) -> Self {
        Self(name.into())
    }

    /// **Stable** — underlying [`Var`].
    pub fn as_var(&self) -> &Var {
        &self.0
    }
}

impl From<Var> for MainVar {
    // **Stable** — `Poly::from`
    fn from(v: Var) -> Self {
        Self(v)
    }
}

/// Element of ℚ[others] (coefficient ring; no distinguished main variable).
#[derive(Clone, Copy, Debug)]
pub struct CoeffRingPoly<'a>(pub &'a Poly);

impl<'a> CoeffRingPoly<'a> {
    /// **Stable** — view `p` as a coefficient-ring element.
    pub fn new(p: &'a Poly) -> Self {
        Self(p)
    }

    /// **Stable** — gcd in ℚ[others] via subresultant (typed coefficient-ring path).
    pub fn gcd(self, other: CoeffRingPoly<'_>) -> Poly {
        use crate::subresultant::subresultant_gcd;
        if self.0.is_zero() {
            return if other.0.is_zero() {
                Poly::one()
            } else {
                other.0.clone().monic()
            };
        }
        if other.0.is_zero() {
            return self.0.clone().monic();
        }
        subresultant_gcd(self.0, other.0).monic()
    }

    /// **Stable** — exact quotient in ℚ[others]; `None` if `den` does not divide `self`.
    pub fn exact_quo(&self, den: &CoeffRingPoly<'_>) -> Option<Poly> {
        div_exact_coeff(self.0, den.0)
    }
}

/// Alias for [`UnivariateIn`] (nested ring ℚ[others][main]).
pub type UnivariateOver<'a> = UnivariateIn<'a>;

/// Borrowed view: `p` ∈ ℚ[others][main].
#[derive(Clone, Debug)]
pub struct UnivariateIn<'a> {
    pub poly: &'a Poly,
    pub main: MainVar,
}

impl<'a> UnivariateIn<'a> {
    /// **Stable** — view `poly` as univariate in `main`.
    pub fn new(poly: &'a Poly, main: MainVar) -> Self {
        Self { poly, main }
    }

    /// **Stable** — main variable.
    pub fn main_var(&self) -> &Var {
        self.main.as_var()
    }

    /// **Stable** — degree w.r.t. `main`.
    pub fn degree(&self) -> u64 {
        univariate_degree(self.poly, self.main.as_var())
    }

    /// **Stable** — coefficient of `main^exp` as ℚ[others].
    pub fn coeff_at(&self, exp: u64) -> Poly {
        coeff_wrt_impl(self.poly, self.main.as_var(), exp)
    }

    /// **Stable** — leading coefficient w.r.t. `main` (element of ℚ[others]).
    pub fn leading_coeff(&self) -> Poly {
        let d = self.degree();
        if d == 0 {
            return self.poly.clone();
        }
        self.coeff_at(d)
    }

    /// **Stable** — whether `self` divides `p` in ℚ[others][main].
    ///
    /// Uses [`quo_exact_wrt`], **not** [`Poly::div_rem`].
    pub fn divides(&self, p: &Poly) -> bool {
        self.exact_quo_dividing(p).is_ok()
    }

    /// **Stable** — exact quotient `p / self` in ℚ[others][main].
    pub fn exact_quo_dividing(&self, p: &Poly) -> PolyResult<Poly> {
        quo_exact_wrt(p, self.poly, self.main.as_var())
    }

    /// **Stable (crate-internal)** — `(q, r)` with `rem = q*self + r` in ℚ[aux][main].
    ///
    /// Divisor [`Self`] must be independent of `aux` (Hensel lift step). Uses
    /// [`univariate_div_rem_wrt`], **not** [`Poly::div_rem`].
    /// **Stable** — `div_rem_wrt_aux_indep`
    pub fn div_rem_wrt_aux_indep(
        &self,
        rem: &Poly,
        aux: &Var,
    ) -> Option<(Poly, Poly)> {
        div_rem_wrt_aux_indep(rem, self.poly, &self.main, aux)
    }

    /// **Stable** — substitute `aux ↦ value` in ℚ[others][main]; `main` unchanged.
    pub fn eval_aux(&self, aux: &Var, value: &CoeffRingPoly<'_>) -> UnivariatePoly {
        let sub = value.0;
        let d = univariate_degree(self.poly, aux);
        let mut out = Poly::zero();
        for e in 0..=d {
            let c = coeff_wrt_impl(self.poly, aux, e);
            if c.is_zero() {
                continue;
            }
            out = out.add(&c.mul(&sub.pow(e)));
        }
        UnivariatePoly::new(out, self.main.clone())
    }
}

/// Owned element of ℚ[others][main].
#[derive(Clone, Debug, PartialEq)]
pub struct UnivariatePoly {
    pub poly: Poly,
    pub main: MainVar,
}

impl UnivariatePoly {
    /// **Stable** — construct owned nested-ring element.
    pub fn new(poly: Poly, main: MainVar) -> Self {
        Self { poly, main }
    }

    /// **Stable** — borrowed view.
    pub fn as_view(&self) -> UnivariateIn<'_> {
        UnivariateIn::new(&self.poly, self.main.clone())
    }
}

/// ℚ[var] flat univariate polynomial — **only** factor subpath that may use classic `div_rem`.
#[derive(Clone, Debug, PartialEq)]
pub struct FlatUni {
    poly: Poly,
    var: MainVar,
}

impl FlatUni {
    /// **Stable** — view `poly` as univariate in `var` over ℚ.
    pub fn new(poly: Poly, var: impl Into<MainVar>) -> Self {
        Self {
            poly,
            var: var.into(),
        }
    }

    /// **Stable** — construct when `poly` is genuinely univariate in `var`.
    pub fn try_new(poly: Poly, var: &Var) -> Option<Self> {
        if !is_flat_univariate_in(&poly, var) {
            return None;
        }
        Some(Self::new(poly, MainVar::new(var.clone())))
    }

    /// **Stable** — underlying polynomial (explicit downgrade).
    pub fn as_poly(&self) -> &Poly {
        &self.poly
    }

    /// **Stable** — consume and return inner [`Poly`].
    pub fn into_poly(self) -> Poly {
        self.poly
    }

    /// **Stable** — main variable.
    pub fn var(&self) -> &MainVar {
        &self.var
    }

    /// **Stable** — Euclidean `(q, r)` in ℚ[var] via [`univariate_div_rem_wrt`].
    pub fn div_rem(&self, divisor: &Self) -> (Poly, Poly) {
        univariate_div_rem_wrt(
            &self.poly,
            &divisor.poly,
            self.var.as_var(),
        )
    }

    /// **Stable** — whether `divisor` divides `self` in ℚ[var].
    pub fn divides(&self, divisor: &Self) -> bool {
        self.div_rem(divisor).1.is_zero()
    }

    /// **Stable** — exact quotient `self / divisor` when remainder is zero.
    pub fn exact_quo(&self, divisor: &Self) -> Option<Poly> {
        let (q, r) = self.div_rem(divisor);
        if r.is_zero() {
            Some(q)
        } else {
            None
        }
    }
}

/// Multivariate polynomial in the erased display ring — explicit home for leading-monomial `div_rem`.
#[derive(Clone, Debug, PartialEq)]
pub struct MultivariatePoly(pub Poly);

impl MultivariatePoly {
    /// **Stable** — wrap a general sparse polynomial.
    pub fn new(poly: Poly) -> Self {
        Self(poly)
    }

    /// **Stable** — underlying [`Poly`] (explicit downgrade).
    pub fn as_poly(&self) -> &Poly {
        &self.0
    }

    /// **Stable** — consume and return inner [`Poly`].
    pub fn into_inner(self) -> Poly {
        self.0
    }

    /// **Stable** — multivariate leading-monomial division (display / heuristic contexts only).
    pub fn div_rem(&self, divisor: &Poly) -> (Poly, Poly) {
        self.0.div_rem(divisor)
    }
}

/// `primitive_part_wrt(p, var)` tagged as ℚ[others][var].
#[derive(Clone, Debug, PartialEq)]
pub struct PrimitivePart {
    poly: Poly,
    wrt: MainVar,
}

impl PrimitivePart {
    /// **Stable** — `p / content_wrt(p, var)` in ℚ[others][var].
    pub fn wrt(p: &Poly, var: &Var) -> PolyResult<Self> {
        let pp = primitive_part_wrt_impl(p, var);
        if pp.is_zero() && !p.is_zero() {
            return Err(crate::error::PolyError::TypeError("zero primitive part"));
        }
        Ok(Self {
            poly: pp,
            wrt: MainVar::new(var.clone()),
        })
    }

    /// **Stable** — primitive part polynomial.
    pub fn as_poly(&self) -> &Poly {
        &self.poly
    }

    /// **Stable** — variable the part is primitive w.r.t.
    pub fn main_var(&self) -> &MainVar {
        &self.wrt
    }

    /// **Stable** — consume into `(poly, main)`.
    pub fn into_inner(self) -> (Poly, MainVar) {
        (self.poly, self.wrt)
    }
}

/// sparse_bi auxiliary dilation `aux ↦ factor * aux` with undo.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DilationMap {
    pub aux_a: Var,
    pub aux_b: Var,
    pub da: i64,
    pub db: i64,
}

impl DilationMap {
    /// **Stable (crate-internal)** — upstream deterministic dilation presets.
    pub const PRESETS: [(i64, i64); 4] = [(2, -1), (-1, 2), (2, 2), (-1, -1)];

    // **Stable** — `Poly::pair`
    pub(crate) fn pair(da: i64, db: i64, aux_a: Var, aux_b: Var) -> Self {
        Self {
            aux_a,
            aux_b,
            da,
            db,
        }
    }

    // **Stable** — `Poly::apply`
    pub(crate) fn apply(&self, p: &Poly) -> Poly {
        let mut out = dilate_one(p, &self.aux_a, self.da);
        dilate_one(&out, &self.aux_b, self.db)
    }

    // **Stable** — `undo`
    pub(crate) fn undo(&self, p: &Poly) -> Poly {
        let mut out = undilate_one(p, &self.aux_b, self.db);
        undilate_one(&out, &self.aux_a, self.da)
    }
}

// **Pipeline private** — `is_flat_univariate_in`
fn is_flat_univariate_in(p: &Poly, var: &Var) -> bool {
    p.terms
        .keys()
        .all(|m| m.iter().all(|(v, _)| v == var))
}

/// **Pipeline private** — substitute `sub_var ↦ sub_poly` (coefficients may depend on other vars).
pub(crate) fn substitute_wrt(p: &Poly, sub_var: &Var, sub_poly: &Poly) -> Poly {
    let d = univariate_degree(p, sub_var);
    let mut out = Poly::zero();
    for e in 0..=d {
        let c = coeff_wrt_impl(p, sub_var, e);
        if c.is_zero() {
            continue;
        }
        out = out.add(&c.mul(&sub_poly.pow(e)));
    }
    out
}

// **Pipeline private** — `dilate_one`
fn dilate_one(p: &Poly, aux: &Var, factor: i64) -> Poly {
    if factor == 1 {
        return p.clone();
    }
    let sub = if factor == -1 {
        Poly::var(aux.clone()).neg()
    } else {
        Poly::var(aux.clone()).mul_scalar(&Ratio::from_integer(BigInt::from(factor)))
    };
    substitute_wrt(p, aux, &sub)
}

// **Pipeline private** — `undilate_one`
fn undilate_one(p: &Poly, aux: &Var, factor: i64) -> Poly {
    if factor == 1 {
        return p.clone();
    }
    let inv = Ratio::from_integer(BigInt::from(1))
        / Ratio::from_integer(BigInt::from(factor.abs()));
    let sub = if factor == -1 {
        Poly::var(aux.clone()).neg()
    } else {
        Poly::var(aux.clone()).mul_scalar(&inv)
    };
    substitute_wrt(p, aux, &sub)
}

/// **Stable (crate-internal)** — single-aux dilation `aux ↦ factor * aux`.
// **Stable** — `dilate_aux`
pub(crate) fn dilate_aux(p: &Poly, aux: &Var, factor: i64) -> Poly {
    dilate_one(p, aux, factor)
}

/// **Stable (crate-internal)** — undo single-aux dilation.
// **Stable** — `undilate_aux`
pub(crate) fn undilate_aux(p: &Poly, aux: &Var, factor: i64) -> Poly {
    undilate_one(p, aux, factor)
}

/// `eval_tn` embedding for sparse_bi: `aux[i] ↦ t^{n[i]}`.
#[derive(Clone, Debug, PartialEq)]
pub struct TnEmbed {
    pub main: MainVar,
    pub aux: [Var; 2],
    pub n: [usize; 2],
    pub t: MainVar,
}

impl TnEmbed {
    /// **Stable** — standard sparse_bi embed tag for auxiliary `t`.
    pub const EMBED_T: &'static str = "__sparse_t__";

    /// **Stable** — `n = [1, 1]` embed with fresh `t` variable.
    pub fn unit(main: impl Into<Var>, aux_a: impl Into<Var>, aux_b: impl Into<Var>) -> Self {
        Self {
            main: MainVar::new(main),
            aux: [aux_a.into(), aux_b.into()],
            n: [1, 1],
            t: MainVar::new(Var::from(Self::EMBED_T)),
        }
    }

    /// **Stable** — copy with new exponent vector.
    pub fn with_n(mut self, n: [usize; 2]) -> Self {
        self.n = n;
        self
    }

    /// **Stable** — map `p ∈ ℚ[main, aux…]` into ℚ[main, t] (typed embed result).
    pub fn embed(&self, p: &Poly) -> BivariateEmbed {
        BivariateEmbed {
            poly: eval_tn_impl(
                p,
                self.main.as_var(),
                &self.aux[0],
                &self.aux[1],
                self.n[0],
                self.n[1],
                self.t.as_var(),
            ),
            embed: self.clone(),
        }
    }

    /// **Stable** — view embedded poly as univariate in `main` (coeffs in ℚ[t]).
    pub fn view_main<'a>(&self, p: &'a Poly) -> UnivariateIn<'a> {
        UnivariateIn::new(p, self.main.clone())
    }

    /// **Stable** — view embedded poly as univariate in `t`.
    pub fn view_t<'a>(&self, p: &'a Poly) -> UnivariateIn<'a> {
        UnivariateIn::new(p, self.t.clone())
    }
}

/// Embedded image `p ∈ ℚ[main, t]` after [`TnEmbed::embed`]; retains `(main, t, n, aux)`.
#[derive(Clone, Debug, PartialEq)]
pub struct BivariateEmbed {
    poly: Poly,
    embed: TnEmbed,
}

impl BivariateEmbed {
    /// **Stable** — construct from polynomial and embed configuration.
    pub fn new(poly: Poly, embed: TnEmbed) -> Self {
        Self { poly, embed }
    }

    /// **Stable** — underlying polynomial (explicit downgrade).
    pub fn as_poly(&self) -> &Poly {
        &self.poly
    }

    /// **Stable** — consume and return the inner [`Poly`].
    pub fn into_poly(self) -> Poly {
        self.poly
    }

    /// **Stable** — embed configuration.
    pub fn config(&self) -> &TnEmbed {
        &self.embed
    }

    /// **Stable** — view as univariate in `main` (coeffs in ℚ[t]).
    pub fn view_main(&self) -> UnivariateIn<'_> {
        self.embed.view_main(&self.poly)
    }

    /// **Stable** — view as univariate in `t`.
    pub fn view_t(&self) -> UnivariateIn<'_> {
        self.embed.view_t(&self.poly)
    }

    /// **Stable (crate-internal)** — sparse_bi reconstruction draft from this candidate.
    // **Stable** — `to_recon_draft`
    pub(crate) fn to_recon_draft(&self, seldegs: &[u64]) -> EmbedFactorDraft {
        EmbedFactorDraft::from_embedded(self, seldegs)
    }
}

/// sparse_bi embed-stage IR: monomials + aux exponents before materializing ℚ[main, aux].
#[derive(Clone, Debug)]
pub(crate) struct EmbedFactorDraft {
    pub monos: Vec<EmbedMonomial>,
    pub aux_exps: Vec<(u64, u64)>,
    pub embed: TnEmbed,
    pub seldegs: Vec<u64>,
}

impl EmbedFactorDraft {
    // **Stable** — `Poly::from_embedded`
    pub(crate) fn from_embedded(selp: &BivariateEmbed, seldegs: &[u64]) -> Self {
        let main = selp.embed.main.as_var();
        let t = selp.embed.t.as_var();
        let monos = EmbedMonomial::sorted_from_poly(&selp.poly, main, t);
        Self {
            aux_exps: vec![(0, 0); monos.len()],
            monos,
            embed: selp.embed.clone(),
            seldegs: seldegs.to_vec(),
        }
    }

    // **Stable** — `Poly::from_poly`
    pub(crate) fn from_poly(selp: &Poly, embed: &TnEmbed, seldegs: &[u64]) -> Self {
        Self::from_embedded(&BivariateEmbed::new(selp.clone(), embed.clone()), seldegs)
    }

    /// Rebuild factor in ℚ[main, aux] from monomial IR (no further embed round-trips).
    // **Stable** — `materialize`
    pub(crate) fn materialize(&self) -> Option<UnivariatePoly> {
        let main = self.embed.main.as_var();
        let aux_a = &self.embed.aux[0];
        let aux_b = &self.embed.aux[1];
        let recon = embed_monomials_to_poly(
            &self.monos,
            main,
            aux_a,
            aux_b,
            &self.aux_exps,
        );
        let pp = PrimitivePart::wrt(&recon, main)
            .ok()
            .map(|t| t.poly)
            .filter(|pp| !pp.is_zero())?;
        Some(UnivariatePoly::new(pp, self.embed.main.clone()))
    }
}

/// `pzadic`-stage IR: univariate factor at eval before lift into `ℚ[eval_var][main]`.
///
/// Upstream `pzadic` expands dimension by one (`dim+1`): digit exponents along `eval_var`
/// encode base-`base` expansions (centered symmetric `smod`). Implemented in
/// [`crate::factor::unitary::PzadicLift::pzadic`].
#[derive(Clone, Debug)]
pub(crate) struct PzadicDraft {
    /// `f ∈ ℚ[main]` after `eval_var ↦ base` (coeffs constant w.r.t. `eval_var`).
    pub factor_at_eval: Poly,
    pub main: MainVar,
    pub eval_var: Var,
    pub base: BigInt,
    /// Upstream `dim+1` sketch: one new digit axis (`eval_var` exponents).
    pub lifted_dim: usize,
}

impl PzadicDraft {
    // **Stable** — `Poly::from_eval_factor`
    pub(crate) fn from_eval_factor(
        factor_at_eval: Poly,
        main: MainVar,
        eval_var: Var,
        base: BigInt,
    ) -> Self {
        Self {
            factor_at_eval,
            main,
            eval_var,
            base,
            lifted_dim: 1,
        }
    }
}

/// Lifted factor candidate in `ℚ[eval_var][main]` (post-pzadic, pre-`divides` check).
///
/// Peel with [`LiftedFactor::as_univariate_in`], **not** [`Poly::div_rem`].
#[derive(Clone, Debug)]
pub(crate) struct LiftedFactor {
    pub poly: Poly,
    pub main: MainVar,
    pub candidate_idx: usize,
}

impl LiftedFactor {
    // **Stable** — `new`
    pub(crate) fn new(poly: Poly, main: MainVar, candidate_idx: usize) -> Self {
        Self {
            poly,
            main,
            candidate_idx,
        }
    }

    /// **Stable (crate-internal)** — nested-ring divisor view for peel.
    // **Stable** — `Poly::as_univariate_in`
    pub(crate) fn as_univariate_in(&self) -> UnivariateIn<'_> {
        UnivariateIn::new(&self.poly, self.main.clone())
    }
}

/// **Stable (crate-internal)** — `rem` has no exponent of `var`.
// **Stable** — `is_independent_of_var`
pub(crate) fn is_independent_of_var(p: &Poly, var: &Var) -> bool {
    p.terms.keys().all(|m| m.exp_of(var) == 0)
}

/// **Stable (crate-internal)** — division in ℚ[aux][main] when `div` is independent of `aux`.
// **Stable** — `div_rem_wrt_aux_indep`
pub(crate) fn div_rem_wrt_aux_indep(
    rem: &Poly,
    div: &Poly,
    main: &MainVar,
    aux: &Var,
) -> Option<(Poly, Poly)> {
    if !is_independent_of_var(div, aux) {
        return None;
    }
    let main_var = main.as_var();
    let dd = univariate_degree(div, main_var);
    if dd == 0 {
        return None;
    }
    let lc = coeff_wrt_impl(div, main_var, dd);
    if lc.is_zero() {
        return None;
    }
    Some(univariate_div_rem_wrt(rem, div, main_var))
}

/// One monomial in `(main, t)` embed space (sparse_bi reconstruction IR).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct EmbedMonomial {
    pub main_e: u64,
    pub t_e: u64,
    pub coeff: Ratio<BigInt>,
}

impl EmbedMonomial {
    // **Stable** — `Poly::sorted_from_poly`
    pub(crate) fn sorted_from_poly(p: &Poly, main: &Var, t: &Var) -> Vec<Self> {
        let mut out = Vec::new();
        for (m, c) in &p.terms {
            out.push(Self {
                main_e: m.exp_of(main),
                t_e: m.exp_of(t),
                coeff: c.clone(),
            });
        }
        out.sort_by(|a, b| {
            b.main_e
                .cmp(&a.main_e)
                .then_with(|| b.t_e.cmp(&a.t_e))
                .then_with(|| b.coeff.cmp(&a.coeff))
        });
        out
    }
}

// **Pipeline private** — build `Poly` from sorted embed monomials + aux exponents
pub(crate) fn embed_monomials_to_poly(
    monos: &[EmbedMonomial],
    main: &Var,
    aux_a: &Var,
    aux_b: &Var,
    aux_exps: &[(u64, u64)],
) -> Poly {
    let mut out = Poly::zero();
    for (mono, &(ae, be)) in monos.iter().zip(aux_exps.iter()) {
        if mono.coeff.is_zero() {
            continue;
        }
        let mut term = term_with_var(&Poly::constant(mono.coeff.clone()), main, mono.main_e);
        if ae > 0 {
            term = term.mul(&Poly::var(aux_a.clone()).pow(ae));
        }
        if be > 0 {
            term = term.mul(&Poly::var(aux_b.clone()).pow(be));
        }
        out = out.add(&term);
    }
    out
}

// **Pipeline private** — shared with `poly_uni::term_with_var`
pub(crate) fn term_with_var(coeff: &Poly, var: &Var, exp: u64) -> Poly {
    if exp == 0 {
        return coeff.clone();
    }
    coeff.mul(&Poly::var(var.clone()).pow(exp))
}

// **Pipeline private**
fn coeff_wrt_impl(p: &Poly, var: &Var, exp: u64) -> Poly {
    crate::subresultant::coeff_wrt(p, var, exp)
}

// **Pipeline private** — `eval_tn`
fn eval_tn_impl(
    p: &Poly,
    main: &Var,
    aux_a: &Var,
    aux_b: &Var,
    n_a: usize,
    n_b: usize,
    t: &Var,
) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        let xe = m.exp_of(main);
        let te = m.exp_of(aux_a) * n_a as u64 + m.exp_of(aux_b) * n_b as u64;
        if te == 0 {
            out = out.add(&term_with_var(&Poly::constant(c.clone()), main, xe));
        } else {
            let tc = Poly::var(t.clone()).pow(te);
            out = out.add(&term_with_var(&tc.mul_scalar(c), main, xe));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::monomial::Var;
    use num_traits::One;

    #[test]
    fn flat_uni_div_rem_vs_multivariate_div_rem() {
        let x = Poly::var("x");
        let p = x.pow(3).sub(&Poly::one());
        let lin = x.sub(&Poly::one());
        let flat = FlatUni::new(p.clone(), MainVar::new("x"));
        let div = FlatUni::new(lin.clone(), MainVar::new("x"));
        let (_, r) = flat.div_rem(&div);
        assert!(r.is_zero());
        let q = flat.exact_quo(&div).expect("quotient");
        assert_eq!(q.mul(&lin), p);
    }

    #[test]
    fn primitive_part_wrt_tags_main() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.mul(&y).add(&x.pow(2));
        let pp = PrimitivePart::wrt(&p, &Var::from("x")).expect("pp");
        assert_eq!(pp.main_var().as_var(), &Var::from("x"));
        assert!(!pp.as_poly().is_zero());
    }

    #[test]
    fn dilation_map_apply_undo() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        let p = x.add(&y).add(&z);
        let map = DilationMap::pair(2, -1, Var::from("y"), Var::from("z"));
        let dil = map.apply(&p);
        let back = map.undo(&dil);
        assert_eq!(back, p);
    }

    #[test]
    fn multivariate_poly_div_rem_extract_powers() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.pow(2).mul(&y);
        let divisor = x.mul(&y);
        let mv = MultivariatePoly::new(p.clone());
        let (rest, rem) = mv.div_rem(&divisor);
        assert!(rem.is_zero());
        assert_eq!(rest, x);
    }

    #[test]
    fn coeff_ring_gcd_divides_both() {
        let y = Poly::var("y");
        let g = CoeffRingPoly::new(&y).gcd(CoeffRingPoly::new(&y.pow(2)));
        assert_eq!(g, y);
        assert!(CoeffRingPoly::new(&y.pow(2)).exact_quo(&CoeffRingPoly::new(&g)).is_some());
    }

    #[test]
    fn univariate_in_divides_vs_div_rem() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        let p = x.add(&y.pow(2).mul(&z.pow(3))).mul(&x.add(&Poly::one()));
        let factor = x.add(&y.pow(2).mul(&z.pow(3)));
        let main = MainVar::new("x");
        let f = UnivariateIn::new(&factor, main);

        assert!(f.divides(&p), "quo_exact_wrt must see factor");
        let (_, rem) = p.div_rem(&factor);
        assert!(
            !rem.is_zero(),
            "multivariate div_rem must not be used for nested-ring divisibility"
        );
        let q = f.exact_quo_dividing(&p).expect("quotient");
        assert_eq!(q.mul(&factor), p);
    }

    #[test]
    fn lifted_factor_peel_vs_div_rem() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.add(&y).sub(&Poly::one()).mul(&x.sub(&y).sub(&Poly::one()));
        let factor = x.add(&y).sub(&Poly::one());
        let main = MainVar::new("x");
        let lifted = LiftedFactor::new(factor.clone(), main.clone(), 0);
        assert!(lifted.as_univariate_in().divides(&p));
        let (_, rem) = p.div_rem(&factor);
        assert!(
            !rem.is_zero(),
            "LiftedFactor peel must not use multivariate div_rem"
        );
    }

    #[test]
    fn coeff_ring_exact_quo() {
        let y = Poly::var("y");
        let z = Poly::var("z");
        let a = y.mul(&z);
        let b = y;
        let q = CoeffRingPoly::new(&a)
            .exact_quo(&CoeffRingPoly::new(&b))
            .expect("z divides yz");
        assert_eq!(q, z);
    }

    #[test]
    fn eval_aux_preserves_main() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.pow(2).add(&x.mul(&y)).add(&Poly::one());
        let main = MainVar::new("x");
        let view = UnivariateIn::new(&p, main);
        let at_y2 = Poly::var("y").pow(2);
        let ev = view.eval_aux(&Var::from("y"), &CoeffRingPoly::new(&at_y2));
        assert_eq!(ev.poly, x.pow(2).add(&x.mul(&at_y2)).add(&Poly::one()));
        assert_eq!(ev.main.as_var(), &Var::from("x"));
    }

    #[test]
    fn tn_embed_maps_aux_to_t() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        let p = x.add(&y).add(&z);
        let emb = TnEmbed::unit("x", "y", "z");
        let pt = emb.embed(&p);
        assert_eq!(pt.view_main().coeff_at(1), Poly::one());
        assert_eq!(
            pt.view_t().coeff_at(1),
            Poly::constant(Ratio::from_integer(2.into()))
        );
    }

    #[test]
    fn embed_factor_draft_materialize_with_aux_exps() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let _z = Poly::var("z");
        let factor = x.add(&y);
        let emb = TnEmbed::unit("x", "y", "z");
        let embedded = emb.embed(&factor);
        let mut draft = embedded.to_recon_draft(&[1]);
        draft.aux_exps = vec![(0, 0), (1, 0)];
        let up = draft.materialize().expect("materialize");
        assert_eq!(up.poly, factor);
    }

    #[test]
    fn div_rem_wrt_aux_indep_matches_univariate() {
        let x = Var::from("x");
        let y = Var::from("y");
        let main = MainVar::new(x.clone());
        let f0 = Poly::var(x.clone()).add(&Poly::one());
        let rem = Poly::var(x.clone())
            .pow(2)
            .add(&Poly::var(x.clone()).mul(&Poly::var(y.clone())))
            .add(&Poly::one());
        let div = UnivariateIn::new(&f0, main.clone());
        let (q, r) = div.div_rem_wrt_aux_indep(&rem, &y).expect("div");
        assert!(r.is_zero() || univariate_degree(&r, &x) < univariate_degree(&f0, &x));
        assert_eq!(q.mul(&f0).add(&r), rem);
    }
}
