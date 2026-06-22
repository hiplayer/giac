//! Shared test fixtures for `giac-core::algebra` (not production API).
//!
//! Tower examples follow [GIAC-lazy-common-tower-plan.md] §11.4 (T1a/T1b).

use std::sync::Arc;

use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::expr::{Expr, ExprArc, FuncKind};

use super::alg_ext::AlgExtData;
use super::ext_tower::ExtensionField;
use super::field_arith::{coords_to_expr, CoordsQ};

/// Second `Arc` with the same tower/minpoly but a distinct field id (for merge tests).
pub fn duplicate_field_arc(f: &Arc<ExtensionField>) -> Arc<ExtensionField> {
    super::ext_tower::duplicate_field_arc_for_test(f)
}

/// Monic `x² + c` in giac `poly1` order (e.g. `c = -2` → x²−2).
pub fn minpoly_x2_plus_c(c: i64) -> CoordsQ {
    vec![
        Ratio::one(),
        Ratio::zero(),
        Ratio::from_integer(c.into()),
    ]
}

/// Alias used in tower tests (`x² + n`).
pub fn minpoly_u2_minus(n: i64) -> CoordsQ {
    minpoly_x2_plus_c(n)
}

/// Monic `x³ − 2` over ℚ (∛2 adjoin).
pub fn minpoly_x3_minus_2() -> CoordsQ {
    vec![
        Ratio::one(),
        Ratio::zero(),
        Ratio::zero(),
        Ratio::from_integer((-2).into()),
    ]
}

/// `Expr` minpoly for x²−2.
pub fn sqrt2_minpoly_expr() -> ExprArc {
    Arc::new(Expr::Func(
        FuncKind::Poly1,
        vec![Arc::new(Expr::Seq(vec![
            Expr::int(1),
            Expr::int(0),
            Expr::int(-2),
        ]))],
    ))
}

/// `Expr` minpoly for x³−2.
pub fn cbrt2_minpoly_expr() -> ExprArc {
    Arc::new(Expr::Func(
        FuncKind::Poly1,
        vec![Arc::new(Expr::Seq(vec![
            Expr::int(1),
            Expr::int(0),
            Expr::int(0),
            Expr::int(-2),
        ]))],
    ))
}

/// Rootof numerator `[1, 0]` = √2 in K₁.
pub fn sqrt2_rootof_num_expr() -> ExprArc {
    Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)]))
}

/// Rootof numerator `[1, 0, 0]` = ∛2 in K₃.
pub fn cbrt2_rootof_num_expr() -> ExprArc {
    Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0), Expr::int(0)]))
}

/// K₁ = ℚ(√2).
pub fn k1_adjoin_sqrt2() -> Arc<ExtensionField> {
    ExtensionField::adjoin_irreducible(
        &ExtensionField::rational(),
        minpoly_x2_plus_c(-2),
    )
    .unwrap()
}

/// ℚ(√3) simple extension (plan: distinct from K₃ = ℚ(∛2) in §11.4).
pub fn k1_adjoin_sqrt3() -> Arc<ExtensionField> {
    ExtensionField::adjoin_irreducible_over_q(minpoly_x2_plus_c(-3)).unwrap()
}

/// K₃ = ℚ(∛2) via legacy `adjoin_irreducible_over_q`.
pub fn k1_adjoin_cbrt2() -> Arc<ExtensionField> {
    ExtensionField::adjoin_irreducible_over_q(minpoly_x3_minus_2()).unwrap()
}

/// [`AlgExtData`] for −√2 in K₁.
pub fn neg_sqrt2_algext() -> AlgExtData {
    AlgExtData::from_rootof(
        &Arc::new(Expr::Seq(vec![Expr::int(-1), Expr::int(0)])),
        &sqrt2_minpoly_expr(),
    )
    .unwrap()
}

/// [`AlgExtData`] for ∛2 in K₃.
pub fn cbrt2_algext() -> AlgExtData {
    AlgExtData::from_rootof(&cbrt2_rootof_num_expr(), &cbrt2_minpoly_expr()).unwrap()
}

/// Rootof numerator `[1, 0]` = √3 in ℚ(√3).
pub fn sqrt3_rootof_num_expr() -> ExprArc {
    Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)]))
}

/// `Expr` minpoly for x²−3.
pub fn sqrt3_minpoly_expr() -> ExprArc {
    Arc::new(Expr::Func(
        FuncKind::Poly1,
        vec![Arc::new(Expr::Seq(vec![
            Expr::int(1),
            Expr::int(0),
            Expr::int(-3),
        ]))],
    ))
}

/// [`AlgExtData`] for √2 in K₁.
pub fn sqrt2_algext() -> AlgExtData {
    AlgExtData::from_rootof(&sqrt2_rootof_num_expr(), &sqrt2_minpoly_expr()).unwrap()
}

/// [`AlgExtData`] for √3 in ℚ(√3).
pub fn sqrt3_algext() -> AlgExtData {
    AlgExtData::from_rootof(&sqrt3_rootof_num_expr(), &sqrt3_minpoly_expr()).unwrap()
}

/// [`Expr::AlgExt`] for √2 in K₁.
pub fn sqrt2_algext_expr() -> ExprArc {
    sqrt2_algext().into_expr()
}

/// T1b: K₂ = K₁(β), β² = 3 (layer minpoly u²−3, constant coeffs in K₁).
pub struct T1bK2 {
    pub k1: Arc<ExtensionField>,
    pub k2: Arc<ExtensionField>,
    /// √2 in K₂ tensor coords.
    pub sqrt2_in_k2: CoordsQ,
    /// β = adjoin generator (√3 in this fixture).
    pub beta_in_k2: CoordsQ,
}

pub fn t1b_k2_adjoin_sqrt3_over_k1() -> T1bK2 {
    let k1 = k1_adjoin_sqrt2();
    let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_x2_plus_c(-3)).unwrap();
    let emb = ExtensionField::try_subfield_embedding(&k1, &k2)
        .unwrap()
        .expect("K1 embeds in K2");
    T1bK2 {
        sqrt2_in_k2: emb.apply(&k1.generator_coords()),
        beta_in_k2: k2.generator_coords(),
        k1,
        k2,
    }
}

/// T3+a: K₂ = K₁(u), u² = √2 (layer minpoly has parent non-constant coefficient).
pub fn t3a_k2_adjoin_u2_minus_sqrt2_over_k1() -> (Arc<ExtensionField>, Arc<ExtensionField>) {
    let k1 = k1_adjoin_sqrt2();
    let zero = k1.zero_coords();
    let one = k1.one_coords();
    let neg_alpha = k1.element_neg(&k1.generator_coords()).unwrap();
    let layer = vec![one, zero, neg_alpha];
    let k2 = ExtensionField::adjoin_irreducible_parent_coeffs(&k1, layer).unwrap();
    (k1, k2)
}

/// Build [`AlgExtData`] from rational coords in the field's operational basis.
pub fn algext_with_coords(field: Arc<ExtensionField>, coords: CoordsQ) -> AlgExtData {
    AlgExtData::from_field_coords(field, coords_to_expr(&coords).unwrap()).unwrap()
}

/// Build [`AlgExtData`] on `fix.k2` from tensor coords.
pub fn algext_on_t1b_k2(fix: &T1bK2, coords: &CoordsQ) -> AlgExtData {
    algext_with_coords(Arc::clone(&fix.k2), coords.clone())
}
