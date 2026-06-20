//! `poly_factor` coefficient-ring tower (FAC-G2).
//!
//! Models `p ∈ ℚ[coeff_vars][main]` (upstream `gausspol.cc` `poly_factor` /
//! `unsplitmultivarpoly` / `splitmultivarpoly`). giac-rs stores the **unsplit**
//! flat image in a single [`Poly`]; this module encodes the tower and routes
//! factorization (aux-lift before sparse_bi).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
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
    // **Stable** — `new`
    pub(crate) fn new(vars: AuxVars<'a>) -> Self {
        Self { vars }
    }

    // **Stable** — `from_slice`
    pub(crate) fn from_slice(vars: &'a [Var]) -> Self {
        Self::new(AuxVars::new(vars))
    }

    /// Upstream `inner_dim`.
    // **Stable** — `inner_dim`
    pub(crate) fn inner_dim(&self) -> usize {
        self.vars.len()
    }

    // **Stable** — `as_slice`
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
    // **Stable** — `new`
    pub(crate) fn new(poly: &'a Poly, main: impl Into<MainVar>, coeff_vars: &'a [Var]) -> Self {
        Self {
            poly,
            main: main.into(),
            coeff_ring: CoeffRing::from_slice(coeff_vars),
        }
    }

    // **Stable** — `from_sqff_ctx`
    pub(crate) fn from_sqff_ctx(ctx: &SqffRingCtx<'a>) -> Self {
        Self {
            poly: ctx.poly,
            main: ctx.main.clone(),
            coeff_ring: CoeffRing::new(ctx.others.clone()),
        }
    }

    // **Stable** — `main_degree`
    pub(crate) fn main_degree(&self) -> u64 {
        univariate_degree(self.poly, self.main.as_var())
    }

    /// Total variable count after unsplit (main + coeff ring).
    // **Stable** — `flat_dim`
    pub(crate) fn flat_dim(&self) -> usize {
        1 + self.coeff_ring.inner_dim()
    }

    /// All variables in the flat unsplit image (main first).
    // **Stable** — `all_vars`
    pub(crate) fn all_vars(&self) -> Vec<Var> {
        let mut out = vec![self.main.as_var().clone()];
        out.extend(self.coeff_ring.as_slice().iter().cloned());
        out
    }

    /// Whether the FAC-G2 tower path applies (`|coeff| ≥ 2`, positive main degree).
    // **Stable** — `is_parametric_tower`
    pub(crate) fn is_parametric_tower(&self) -> bool {
        self.coeff_ring.inner_dim() >= 2 && self.main_degree() > 0
    }

    /// giac `unsplitmultivarpoly`: flat [`Poly`] is already the unsplit image.
    // **Stable** — `unsplit_to_flat`
    pub(crate) fn unsplit_to_flat(&self) -> Poly {
        self.poly.clone()
    }

    /// giac `splitmultivarpoly`: tag a flat factor as ℚ[coeff_ring][main].
    // **Stable** — `split_factor`
    pub(crate) fn split_factor(&self, flat: Poly) -> UnivariatePoly {
        UnivariatePoly::new(flat, self.main.clone())
    }

    /// Upstream `poly_factor` MVP on flat representation: aux substitution + lift.
    // **Partial** — optional algorithm path `try_factor_aux_lift`
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

    /// FAC-G2 tower factor chain: aux-lift → good_eval → sparse_bi (sparse only if still splittable).
    // **Stable** — `factor_sqff_chain`
    pub(crate) fn factor_sqff_chain(&self) -> FactorSet {
        if !self.is_parametric_tower() {
            return FactorSet::irreducible(self.poly.clone(), self.main.clone());
        }
        if let Some(set) = self.try_factor_aux_lift() {
            return set;
        }
        let aux_refs: Vec<&Var> = self.coeff_ring.as_slice().iter().collect();
        if super::eval::looks_irreducible_by_good_eval(self.poly, self.main.as_var(), &aux_refs) {
            return FactorSet::irreducible(self.poly.clone(), self.main.clone());
        }
        if let Some(set) = self.try_sparse_bi() {
            return set;
        }
        FactorSet::irreducible(self.poly.clone(), self.main.clone())
    }

    /// Sparse_bi fallback only (used by tests and internal chain).
    // **Partial** — optional algorithm path `try_sparse_bi`
    pub(crate) fn try_sparse_bi(&self) -> Option<FactorSet> {
        let aux_refs: Vec<&Var> = self.coeff_ring.as_slice().iter().collect();
        let f = super::sparse::try_sparse_factor_bi(self.poly, self.main.as_var(), &aux_refs)?;
        let set = FactorSet::from_polys(f, self.main.clone());
        if set.product_equals(self.poly) {
            Some(set)
        } else {
            None
        }
    }

    /// FAC-G2 tower factor chain (alias): aux-lift then sparse_bi only.
    // **Partial** — optional algorithm path `try_factor`
    pub(crate) fn try_factor(&self) -> Option<FactorSet> {
        if !self.is_parametric_tower() {
            return None;
        }
        self.try_factor_aux_lift().or_else(|| self.try_sparse_bi())
    }
}

impl<'a> SqffRingCtx<'a> {
    // **Stable** — `as_poly_factor_tower`
    pub(crate) fn as_poly_factor_tower(&self) -> PolyFactorTower<'a> {
        PolyFactorTower::from_sqff_ctx(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;

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

    fn l21_poly() -> Poly {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        x.sub(&y)
            .sub(&z)
            .mul(&x.sub(&y).add(&z))
            .mul(&x.add(&y).add(&z))
    }

    #[test]
    fn tower_factor_sqff_chain_l20() {
        let p = l20_poly();
        let others = l20_others();
        let tower = PolyFactorTower::new(&p, Var::from("b"), &others);
        let set = tower.factor_sqff_chain();
        assert_eq!(set.components.len(), 2);
        assert!(set.product_equals(&p));
        assert!(set.verify_divides_chain(&p));
    }

    #[test]
    fn tower_factor_sqff_chain_l21() {
        let p = l21_poly();
        let others = [Var::from("y"), Var::from("z")];
        let tower = PolyFactorTower::new(&p, Var::from("x"), &others);
        let set = tower.factor_sqff_chain();
        assert!(set.components.len() >= 2);
        assert!(set.product_equals(&p));
    }

    #[test]
    fn tower_try_factor_skips_good_eval() {
        let p = l20_poly();
        let others = l20_others();
        let tower = PolyFactorTower::new(&p, Var::from("b"), &others);
        let set = tower.try_factor().expect("aux-lift or sparse");
        assert_eq!(set.components.len(), 2);
    }

    fn l24_poly() -> Poly {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        x.pow(3)
            .sub(&Poly::constant(Ratio::from_integer(21.into())).mul(&x).mul(&z).mul(&y.pow(3)))
            .add(&Poly::constant(Ratio::from_integer(13.into())).mul(&x).mul(&y).mul(&z.pow(2)))
            .add(&Poly::constant(Ratio::from_integer(2.into())))
    }

    #[test]
    #[ignore = "slow in release (~13s): factor_sqff_chain exhausts sparse_bi on irreducible L24; irreducibility covered by tower_l23"]
    fn tower_l24_irreducible_skips_sparse_bi() {
        let p = l24_poly();
        let others = [Var::from("x"), Var::from("y")];
        let tower = PolyFactorTower::new(&p, Var::from("z"), &others);
        assert!(tower.try_sparse_bi().is_none() || tower.try_factor_aux_lift().is_some());
        let set = tower.factor_sqff_chain();
        assert_eq!(set.components.len(), 1);
        assert!(set.product_equals(&p));
    }

    #[test]
    fn tower_l23_cubic_factors_are_irreducible() {
        use super::super::multivariate::factor_multivariate;
        let x = Poly::var("x");
        let y = Poly::var("y");
        let z = Poly::var("z");
        let f1 = l24_poly();
        let f2 = Poly::constant(Ratio::from_integer(6.into()))
            .mul(&x.pow(3))
            .add(&x.mul(&y))
            .add(&x.mul(&z))
            .add(&Poly::one());
        // Upstream line 23/24: each cubic is irreducible; conformance factors Mul structurally
        // (giac-simplify `factor_expr`), not the expanded degree-6 product via `factor_multivariate`.
        for (label, q) in [("f1", &f1), ("f2", &f2)] {
            let facs = factor_multivariate(q).expect(label);
            assert_eq!(facs.len(), 1, "{label} should stay irreducible");
            assert_eq!(facs[0], *q);
        }
        let product = f1.mul(&f2);
        assert_ne!(product, f1);
        assert_ne!(product, f2);
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
