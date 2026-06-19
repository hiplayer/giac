//! `poly_factor` coefficient-ring tower (FAC-G2).
//!
//! Models `p ∈ ℚ[coeff_vars][main]` (upstream `gausspol.cc` `poly_factor` /
//! `unsplitmultivarpoly` / `splitmultivarpoly`). giac-rs stores the **unsplit**
//! flat image in a single [`Poly`]; this module encodes the tower and routes
//! factorization (aux-lift before sparse_bi).

use crate::monomial::Var;
use crate::nested::{MainVar, UnivariatePoly};
use crate::poly::Poly;
use crate::resultant::univariate_degree;

use super::ctx::{AuxVars, FactorSet, SqffRingCtx};

/// Coefficient ring ℚ[coeff_vars] of a `poly_factor` tower.
#[derive(Clone, Debug)]
pub(crate) struct CoeffRing<'a> {
    vars: AuxVars<'a>,
}

impl<'a> CoeffRing<'a> {
    pub(crate) fn new(vars: AuxVars<'a>) -> Self {
        Self { vars }
    }

    pub(crate) fn from_slice(vars: &'a [Var]) -> Self {
        Self::new(AuxVars::new(vars))
    }

    /// Upstream `inner_dim`.
    pub(crate) fn inner_dim(&self) -> usize {
        self.vars.len()
    }

    pub(crate) fn as_slice(&self) -> &[Var] {
        self.vars.as_slice()
    }
}

/// `p ∈ ℚ[coeff_ring][main]` — upstream `poly_factor` tower on flat [`Poly`].
#[derive(Clone, Debug)]
pub(crate) struct PolyFactorTower<'a> {
    pub poly: &'a Poly,
    pub main: MainVar,
    pub coeff_ring: CoeffRing<'a>,
}

impl<'a> PolyFactorTower<'a> {
    pub(crate) fn new(poly: &'a Poly, main: impl Into<MainVar>, coeff_vars: &'a [Var]) -> Self {
        Self {
            poly,
            main: main.into(),
            coeff_ring: CoeffRing::from_slice(coeff_vars),
        }
    }

    pub(crate) fn from_sqff_ctx(ctx: &SqffRingCtx<'a>) -> Self {
        Self {
            poly: ctx.poly,
            main: ctx.main.clone(),
            coeff_ring: CoeffRing::new(ctx.others.clone()),
        }
    }

    pub(crate) fn main_degree(&self) -> u64 {
        univariate_degree(self.poly, self.main.as_var())
    }

    /// Total variable count after unsplit (main + coeff ring).
    pub(crate) fn flat_dim(&self) -> usize {
        1 + self.coeff_ring.inner_dim()
    }

    /// All variables in the flat unsplit image (main first).
    pub(crate) fn all_vars(&self) -> Vec<Var> {
        let mut out = vec![self.main.as_var().clone()];
        out.extend(self.coeff_ring.as_slice().iter().cloned());
        out
    }

    /// Whether the FAC-G2 tower path applies (`|coeff| ≥ 2`, positive main degree).
    pub(crate) fn is_parametric_tower(&self) -> bool {
        self.coeff_ring.inner_dim() >= 2 && self.main_degree() > 0
    }

    /// giac `unsplitmultivarpoly`: flat [`Poly`] is already the unsplit image.
    pub(crate) fn unsplit_to_flat(&self) -> Poly {
        self.poly.clone()
    }

    /// giac `splitmultivarpoly`: tag a flat factor as ℚ[coeff_ring][main].
    pub(crate) fn split_factor(&self, flat: Poly) -> UnivariatePoly {
        UnivariatePoly::new(flat, self.main.clone())
    }

    /// Upstream `poly_factor` MVP on flat representation: aux substitution + lift.
    pub(crate) fn try_factor_aux_lift(&self) -> Option<FactorSet> {
        if !self.is_parametric_tower() {
            return None;
        }
        let main = self.main.as_var();
        let coeff = self.coeff_ring.as_slice();
        for av in coeff {
            let Some(f) = super::hensel::try_lift_factors_in_aux_var(self.poly, main, av, coeff) else {
                continue;
            };
            let set = FactorSet::from_polys(f, self.main.clone());
            if set.product_equals(self.poly) {
                return Some(set);
            }
        }
        None
    }

    /// FAC-G2 tower factor chain: aux-lift (poly_factor) then sparse_bi fallback.
    pub(crate) fn try_factor(&self) -> Option<FactorSet> {
        if !self.is_parametric_tower() {
            return None;
        }
        if let Some(set) = self.try_factor_aux_lift() {
            return Some(set);
        }
        let aux_refs: Vec<&Var> = self.coeff_ring.as_slice().iter().collect();
        let f = super::sparse::try_sparse_factor_bi(self.poly, self.main.as_var(), &aux_refs)?;
        let set = FactorSet::from_polys(f, self.main.clone());
        if set.product_equals(self.poly) {
            Some(set)
        } else {
            None
        }
    }
}

impl<'a> SqffRingCtx<'a> {
    pub(crate) fn as_poly_factor_tower(&self) -> PolyFactorTower<'a> {
        PolyFactorTower::from_sqff_ctx(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn l20_poly() -> Poly {
        let x = Poly::var("x");
        let b = Poly::var("b");
        let c = Poly::var("c");
        x.add(&b)
            .add(&c)
            .mul(
                &x.pow(2)
                    .sub(&x.mul(&b))
                    .sub(&x.mul(&c))
                    .add(&b.pow(2))
                    .sub(&b.mul(&c))
                    .add(&c.pow(2)),
            )
    }

    fn l20_others() -> [Var; 2] {
        [Var::from("c"), Var::from("x")]
    }

    #[test]
    fn tower_unsplit_is_flat_poly() {
        let p = l20_poly();
        let others = l20_others();
        let tower = PolyFactorTower::new(&p, Var::from("b"), &others);
        assert!(tower.is_parametric_tower());
        assert_eq!(tower.flat_dim(), 3);
        assert_eq!(tower.unsplit_to_flat(), p);
    }

    #[test]
    fn tower_aux_lift_factors_l20() {
        let p = l20_poly();
        let others = l20_others();
        let tower = PolyFactorTower::new(&p, Var::from("b"), &others);
        let set = tower
            .try_factor_aux_lift()
            .expect("poly_factor aux-lift L20");
        assert_eq!(set.components.len(), 2);
        assert!(set.product_equals(&p));
        assert!(set.verify_divides_chain(&p));
    }

    #[test]
    fn tower_from_sqff_ctx() {
        let p = l20_poly();
        let others = [Var::from("c"), Var::from("x")];
        let ctx = SqffRingCtx::new(&p, Var::from("b"), &others);
        let tower = ctx.as_poly_factor_tower();
        assert_eq!(tower.coeff_ring.inner_dim(), 2);
        assert_eq!(tower.main.as_var(), &Var::from("b"));
    }

    #[test]
    fn split_factor_tags_main() {
        let p = l20_poly();
        let others = l20_others();
        let tower = PolyFactorTower::new(&p, Var::from("b"), &others);
        let f = Poly::var("b").add(&Poly::var("c")).add(&Poly::var("x"));
        let tagged = tower.split_factor(f.clone());
        assert_eq!(tagged.poly, f);
        assert_eq!(tagged.main.as_var(), &Var::from("b"));
    }
}
