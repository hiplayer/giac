//! Factor pipeline context types (Phase 2 nested-ring issue).
//!
//! Encodes ring tower (`SqffRingCtx`), factor lists (`FactorSet`), good eval points
//! (`GoodEval`), and Hensel lift pairs (`HenselPair`).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_rational::Ratio;
use num_bigint::BigInt;

use crate::error::PolyResult;
use crate::monomial::Var;
use crate::nested::{is_independent_of_var, MainVar, UnivariateIn, UnivariatePoly};
use crate::poly::Poly;
use crate::resultant::univariate_degree;

/// Auxiliary variables of ℚ[others][main] (ordered slice).
#[derive(Clone, Debug)]
pub(crate) struct AuxVars<'a> {
    vars: &'a [Var],
}

impl<'a> AuxVars<'a> {
    // **Stable** — `new`
    pub(crate) fn new(vars: &'a [Var]) -> Self {
        Self { vars }
    }

    // **Stable** — `as_slice`
    pub(crate) fn as_slice(&self) -> &[Var] {
        self.vars
    }

    // **Stable** — `len`
    pub(crate) fn len(&self) -> usize {
        self.vars.len()
    }

    // **Stable** — `is_empty`
    pub(crate) fn is_empty(&self) -> bool {
        self.vars.is_empty()
    }

    // **Stable** — `refs`
    pub(crate) fn refs(&self) -> Vec<&'a Var> {
        self.vars.iter().collect()
    }
}

/// Square-free factorization context: `p` ∈ ℚ[others][main].
#[derive(Clone, Debug)]
pub(crate) struct SqffRingCtx<'a> {
    pub poly: &'a Poly,
    pub main: MainVar,
    pub others: AuxVars<'a>,
}

impl<'a> SqffRingCtx<'a> {
    // **Stable** — `new`
    pub(crate) fn new(poly: &'a Poly, main: impl Into<MainVar>, others: &'a [Var]) -> Self {
        Self {
            poly,
            main: main.into(),
            others: AuxVars::new(others),
        }
    }

    // **Stable** — `main_degree`
    pub(crate) fn main_degree(&self) -> u64 {
        univariate_degree(self.poly, self.main.as_var())
    }

    // **Stable** — `with_main_and_others`
    pub(crate) fn with_main_and_others(
        &self,
        main: impl Into<MainVar>,
        others: &'a [Var],
    ) -> Self {
        Self {
            poly: self.poly,
            main: main.into(),
            others: AuxVars::new(others),
        }
    }
}

/// One sqff component with multiplicity in ℚ[others][main].
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct SqffComponent {
    pub factor: UnivariatePoly,
    pub multiplicity: usize,
}

/// Verified-style factor list in ℚ[others][main].
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct FactorSet {
    pub components: Vec<SqffComponent>,
    pub main: MainVar,
}

impl FactorSet {
    // **Stable** — `Poly::irreducible`
    pub(crate) fn irreducible(poly: Poly, main: MainVar) -> Self {
        Self {
            components: vec![SqffComponent {
                factor: UnivariatePoly::new(poly, main.clone()),
                multiplicity: 1,
            }],
            main,
        }
    }

    // **Stable** — `Poly::from_polys`
    pub(crate) fn from_polys(polys: Vec<Poly>, main: MainVar) -> Self {
        Self {
            components: polys
                .into_iter()
                .map(|p| SqffComponent {
                    factor: UnivariatePoly::new(p, main.clone()),
                    multiplicity: 1,
                })
                .collect(),
            main,
        }
    }

    // **Stable** — `into_polys`
    pub(crate) fn into_polys(self) -> Vec<Poly> {
        self.components.into_iter().map(|c| c.factor.poly).collect()
    }

    /// Whether `∏ f_i^{m_i} == p` in ℚ[others][main] (exact product check).
    // **Stable** — `product_equals`
    pub(crate) fn product_equals(&self, p: &Poly) -> bool {
        let mut prod = Poly::one();
        for c in &self.components {
            let mut term = c.factor.poly.clone();
            for _ in 1..c.multiplicity {
                term = term.mul(&c.factor.poly);
            }
            prod = prod.mul(&term);
        }
        prod == *p
    }

    /// Strip to plain factors after nested-ring divisibility checks on each component.
    // **Stable** — `verify_divides_chain`
    pub(crate) fn verify_divides_chain(&self, p: &Poly) -> bool {
        if !self.product_equals(p) {
            return false;
        }
        let mut rest = p.clone();
        for c in &self.components {
            let view = c.factor.as_view();
            if !view.divides(&rest) {
                return false;
            }
            if let Ok(q) = view.exact_quo_dividing(&rest) {
                rest = q;
            } else {
                return false;
            }
        }
        rest.is_one()
    }
}

/// Good evaluation: auxiliary assignment preserving `main`-degree.
#[derive(Clone, Debug, PartialEq)]
pub struct GoodEval {
    pub evaluated: Poly,
    pub values: Vec<Ratio<BigInt>>,
    pub preserved_main_degree: u64,
}

impl GoodEval {
    // **Stable** — `new`
    pub(crate) fn new(evaluated: Poly, values: Vec<Ratio<BigInt>>, preserved_main_degree: u64) -> Self {
        Self {
            evaluated,
            values,
            preserved_main_degree,
        }
    }

    // **Stable** — `Poly::first_value`
    pub(crate) fn first_value(&self) -> Option<&Ratio<BigInt>> {
        self.values.first()
    }
}

/// Hensel lift factor independent of auxiliary variable (ℚ[main] ⊂ ℚ[aux][main]).
#[derive(Clone, Debug)]
pub(crate) struct AuxIndepFactor<'a> {
    view: UnivariateIn<'a>,
    aux: Var,
}

impl<'a> AuxIndepFactor<'a> {
    // **Partial** — optional algorithm path `try_new`
    pub(crate) fn try_new(poly: &'a Poly, main: MainVar, aux: &Var) -> Option<Self> {
        if !is_independent_of_var(poly, aux) {
            return None;
        }
        if univariate_degree(poly, main.as_var()) == 0 {
            return None;
        }
        Some(Self {
            view: UnivariateIn::new(poly, main),
            aux: aux.clone(),
        })
    }

    // **Stable** — `as_view`
    pub(crate) fn as_view(&self) -> &UnivariateIn<'a> {
        &self.view
    }

    // **Stable** — `div_rem_wrt_aux_indep`
    pub(crate) fn div_rem_wrt_aux_indep(
        &self,
        rem: &Poly,
    ) -> Option<(Poly, Poly)> {
        self.view.div_rem_wrt_aux_indep(rem, &self.aux)
    }
}

/// Two-factor Hensel lift @ aux = 0.
#[derive(Clone, Debug)]
pub(crate) struct HenselPair<'a> {
    pub p: &'a Poly,
    pub main: MainVar,
    pub aux: Var,
    pub f0: AuxIndepFactor<'a>,
    pub g0: AuxIndepFactor<'a>,
}

impl<'a> HenselPair<'a> {
    // **Partial** — optional algorithm path `try_new`
    pub(crate) fn try_new(
        p: &'a Poly,
        main: impl Into<MainVar>,
        aux: &Var,
        f0: &'a Poly,
        g0: &'a Poly,
    ) -> Option<Self> {
        let main = main.into();
        Some(Self {
            p,
            main: main.clone(),
            aux: aux.clone(),
            f0: AuxIndepFactor::try_new(f0, main.clone(), aux)?,
            g0: AuxIndepFactor::try_new(g0, main, aux)?,
        })
    }
}

pub(crate) type SqffFactorRecFn =
    for<'a> fn(SqffRingCtx<'a>) -> PolyResult<Vec<Poly>>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factor_set_product_equals() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let f1 = x.sub(&y);
        let f2 = x.add(&y);
        let p = f1.mul(&f2);
        let main = MainVar::new("x");
        let set = FactorSet::from_polys(vec![f1, f2], main);
        assert!(set.product_equals(&p));
        assert!(set.verify_divides_chain(&p));
    }

    #[test]
    fn aux_indep_factor_rejects_aux_dep() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let f = x.add(&y);
        assert!(AuxIndepFactor::try_new(&f, MainVar::new("x"), &Var::from("y")).is_none());
        let g = x.add(&Poly::one());
        assert!(AuxIndepFactor::try_new(&g, MainVar::new("x"), &Var::from("y")).is_some());
    }

    #[test]
    fn hensel_pair_requires_aux_indep() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.sub(&Poly::one()).mul(&x.add(&Poly::one()));
        let f0 = x.sub(&Poly::one());
        let g0 = x.add(&Poly::one());
        assert!(HenselPair::try_new(&p, Var::from("x"), &Var::from("y"), &f0, &g0).is_some());
        let bad = x.add(&y);
        assert!(HenselPair::try_new(&p, Var::from("x"), &Var::from("y"), &bad, &g0).is_none());
    }
}
