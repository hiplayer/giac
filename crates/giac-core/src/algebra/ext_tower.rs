//! Tower of algebraic extensions K₀=ℚ → K₁ → … → Kₙ and cross-field `common`.
//!
//! Normative model: [GIAC-algext-adoption.md](../../../../.doc/issues/GIAC-algext-adoption.md) §8.2.
//!
//! Phase 0 scope: coefficients in ℚ only; each [`ExtensionField`] is described by a
//! primitive/minimal polynomial over ℚ (single adjoin or `common` composite).
//!
//! T1+: nontrivial adjoins record a true [`ExtensionTower::Adj`] parent chain and
//! `(parent_id, min_poly)` registry keys; see [GIAC-lazy-common-tower-plan.md] T1.

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc};

use super::field_arith::{
    apply_linear_map, char_poly_matrix, coords_all_zero, embed_in_square_extension,
    generator_coords, kron_left, mat_mul, mult_matrix_of_element, pad_to_len, poly_add,
    poly_degree, poly_inv_mod, poly_mul, poly_neg, poly_reduce, poly_sub,
    poly_add_with_coeffs_in_field, poly_inv_mod_with_coeffs_in_field,
    poly_mul_with_coeffs_in_field, poly_neg_with_coeffs_in_field,
    poly_reduce_with_coeffs_in_field, ParentCoeffRing, trim_leading_zero, CoordsQ,
};

static FIELD_ID: AtomicU64 = AtomicU64::new(1);

fn next_field_id() -> u64 {
    FIELD_ID.fetch_add(1, Ordering::Relaxed)
}

/// K₀ = ℚ or K = parent(γ) with γ satisfying `min_poly_q` over the parent.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExtensionTower {
    /// K = ℚ, dimension 1.
    Base,
    /// K = parent(γ); `min_poly_q` is the minimal polynomial of γ over the parent (Phase C
    /// before parent-coeff arithmetic: coefficients are rational constants only).
    ///
    /// For towers built via [`ExtensionField::common_over_q`], `min_poly_q` is the
    /// primitive-element polynomial over ℚ describing the entire field.
    Adj {
        parent: Arc<ExtensionTower>,
        min_poly_q: CoordsQ,
        ext_degree: usize,
        /// When this layer adjoins over a nontrivial parent, `k` in θ = α + k·β used to
        /// build [`ExtensionField::min_poly_over_q`] (T1 compose). `None` for flatten `common`.
        primitive_k: Option<i64>,
    },
}

impl ExtensionTower {
    pub fn dimension(&self) -> usize {
        match self {
            ExtensionTower::Base => 1,
            ExtensionTower::Adj {
                parent,
                ext_degree,
                ..
            } => parent.dimension() * ext_degree,
        }
    }

    /// True when this tower is a single adjoin over ℚ (not a `common` composite).
    pub fn is_simple_over_q(&self) -> bool {
        matches!(
            self,
            ExtensionTower::Adj {
                parent,
                ..
            } if matches!(parent.as_ref(), ExtensionTower::Base)
        )
    }
}

/// Shared handle to an extension field K / ℚ(α₁,…).
#[derive(Clone, Debug)]
pub struct ExtensionField {
    id: u64,
    tower: Arc<ExtensionTower>,
    /// Defining polynomial over ℚ for this field (primitive element when composite).
    min_poly_over_q: CoordsQ,
    /// Immediate parent when registered via [`ExtensionField::adjoin_irreducible`].
    parent_field: Option<Arc<ExtensionField>>,
}

impl PartialEq for ExtensionField {
    fn eq(&self, other: &Self) -> bool {
        self.tower == other.tower && self.min_poly_over_q == other.min_poly_over_q
    }
}

impl Eq for ExtensionField {}

macro_rules! parent_coeff_ring {
    ($parent:ident) => {
        ParentCoeffRing {
            zero: $parent.zero_coords(),
            one: $parent.one_coords(),
            add: &|a, b| $parent.element_add_primitive(a, b),
            sub: &|a, b| $parent.element_sub_primitive(a, b),
            mul: &|a, b| $parent.element_mul_primitive(a, b),
            neg: &|a| $parent.element_neg_primitive(a),
            inv: &|a| $parent.element_inv_primitive(a),
            is_zero: &|a| $parent.element_is_zero(a),
        }
    };
}

impl ExtensionField {
    fn is_base(&self) -> bool {
        matches!(self.tower.as_ref(), ExtensionTower::Base)
    }

    /// ℚ as an extension field (dimension 1).
    pub fn rational() -> Arc<Self> {
        static BASE: OnceLock<Arc<ExtensionField>> = OnceLock::new();
        Arc::clone(BASE.get_or_init(|| {
            Arc::new(Self {
                id: next_field_id(),
                tower: Arc::new(ExtensionTower::Base),
                min_poly_over_q: Vec::new(),
                parent_field: None,
            })
        }))
    }

    pub fn id(&self) -> u64 {
        self.id
    }

    pub fn tower(&self) -> &Arc<ExtensionTower> {
        &self.tower
    }

    pub fn dimension(&self) -> usize {
        self.tower.dimension()
    }

    pub fn min_poly_over_q(&self) -> &[Ratio<BigInt>] {
        &self.min_poly_over_q
    }

    /// Immediate parent field when this handle comes from a tower adjoin (T1+).
    pub fn parent_field(&self) -> Option<&Arc<ExtensionField>> {
        self.parent_field.as_ref()
    }

    /// Adjoin a root of irreducible (or degree ≥ 1) `min_poly_q` over `parent`.
    pub fn adjoin_irreducible(
        parent: &Arc<ExtensionField>,
        min_poly_q: CoordsQ,
    ) -> Result<Arc<Self>, EvalError> {
        let min_poly_q = trim_leading_zero(min_poly_q);
        let d = poly_degree(&min_poly_q);
        if d < 1 {
            return Err(EvalError::TypeError("extension degree >= 1"));
        }
        if !parent.is_base() && !layer_minpoly_rational_constants(&min_poly_q) {
            return Err(EvalError::TypeError(
                "T1 adjoin: layer minpoly must have rational coefficients",
            ));
        }
        field_registry().register_adjoin(parent, &min_poly_q)
    }

    /// Build ℚ(α) from an irreducible (or at least degree ≥ 1) polynomial over ℚ.
    pub fn adjoin_irreducible_over_q(min_poly_q: CoordsQ) -> Result<Arc<Self>, EvalError> {
        Self::adjoin_irreducible(&Self::rational(), min_poly_q)
    }

    /// True when `sub` embeds in `sup` along the registered adjoin parent chain (T1).
    pub fn is_subfield_of(sub: &Arc<ExtensionField>, sup: &Arc<ExtensionField>) -> bool {
        if Arc::ptr_eq(sub, sup) || **sub == **sup {
            return true;
        }
        if sub.is_base() {
            return true;
        }
        if sup.is_base() {
            return false;
        }
        let mut cur = sup.parent_field.as_ref();
        while let Some(p) = cur {
            if Arc::ptr_eq(p, sub) || **p == **sub {
                return true;
            }
            cur = p.parent_field.as_ref();
        }
        false
    }

    /// Linear embedding `sub` ↪ `sup` when [`Self::is_subfield_of`] (T1); `None` if unrelated.
    pub fn try_subfield_embedding(
        sub: &Arc<ExtensionField>,
        sup: &Arc<ExtensionField>,
    ) -> Result<Option<FieldEmbedding>, EvalError> {
        if !Self::is_subfield_of(sub, sup) {
            return Ok(None);
        }
        if Arc::ptr_eq(sub, sup) || **sub == **sup {
            let dim = sub.dimension();
            return Ok(Some(FieldEmbedding {
                source: Arc::clone(sub),
                target: Arc::clone(sup),
                matrix: super::field_arith::identity_matrix(dim),
            }));
        }
        if sub.is_base() {
            return Ok(Some(rational_subfield_embedding(sub, sup)));
        }
        let mut chain = Vec::new();
        let mut cur = Arc::clone(sup);
        loop {
            let parent = cur
                .parent_field
                .as_ref()
                .ok_or(EvalError::TypeError("subfield chain broken"))?;
            chain.push(Arc::clone(&cur));
            if Arc::ptr_eq(parent, sub) || **parent == **sub {
                break;
            }
            if parent.is_base() {
                return Ok(None);
            }
            cur = Arc::clone(parent);
        }
        chain.reverse();
        let mut emb = direct_adjoin_parent_embedding(sub, &chain[0])?;
        for child in &chain[1..] {
            let parent = child.parent_field.as_ref().unwrap();
            let step = direct_adjoin_parent_embedding(parent, child)?;
            emb = compose_field_embeddings(&emb, &step)?;
        }
        Ok(Some(emb))
    }

    /// `min_poly` as `Expr` coefficients (giac `poly1` order).
    pub fn top_min_poly_exprs(&self) -> Result<Vec<ExprArc>, EvalError> {
        if self.is_base() {
            return Ok(vec![Expr::int(0)]);
        }
        super::field_arith::coords_to_expr(&self.min_poly_over_q)
    }

    pub fn zero_coords(&self) -> CoordsQ {
        vec![Ratio::zero(); self.dimension()]
    }

    pub fn one_coords(&self) -> CoordsQ {
        if let Some(parent) = self.parent_field.as_ref() {
            if !parent.is_base() {
                let pd = parent.dimension();
                let mut v = self.zero_coords();
                let parent_one = parent.one_coords();
                for i in 0..pd {
                    if i < parent_one.len() {
                        v[i] = parent_one[i].clone();
                    }
                }
                return v;
            }
        }
        let mut v = self.zero_coords();
        if let Some(c) = v.last_mut() {
            *c = Ratio::one();
        }
        v
    }

    pub(crate) fn generator_coords(&self) -> CoordsQ {
        if let Some(parent) = self.parent_field.as_ref() {
            if !parent.is_base() {
                let pd = parent.dimension();
                let mut v = self.zero_coords();
                let one = parent.one_coords();
                let block1 = pd;
                for i in 0..pd {
                    if block1 + i < v.len() && i < one.len() {
                        v[block1 + i] = one[i].clone();
                    }
                }
                return v;
            }
        }
        let mut v = self.zero_coords();
        if !v.is_empty() {
            v[0] = Ratio::one();
        }
        v
    }

    fn uses_tower_arithmetic(&self) -> bool {
        self.parent_field
            .as_ref()
            .is_some_and(|parent| !parent.is_base())
    }

    /// Embed a rational constant into this field (constant term at layer 0).
    pub fn embed_rational(&self, r: &Ratio<BigInt>) -> CoordsQ {
        if self.is_base() {
            return vec![r.clone()];
        }
        if let Some(parent) = self.parent_field.as_ref() {
            if !parent.is_base() {
                let pd = parent.dimension();
                let block0 = parent.embed_rational(r);
                let mut v = self.zero_coords();
                for i in 0..pd {
                    if i < block0.len() && i < v.len() {
                        v[i] = block0[i].clone();
                    }
                }
                return v;
            }
        }
        let mut v = self.zero_coords();
        if let Some(c) = v.last_mut() {
            *c = r.clone();
        }
        v
    }

    fn layer_ext_degree(&self) -> Result<usize, EvalError> {
        match self.tower.as_ref() {
            ExtensionTower::Adj { ext_degree, .. } => Ok(*ext_degree),
            ExtensionTower::Base => Err(EvalError::TypeError("expected extension field")),
        }
    }

    fn layer_minpoly_parent_coeffs(&self, parent: &ExtensionField) -> Result<Vec<CoordsQ>, EvalError> {
        let min_poly_q = match self.tower.as_ref() {
            ExtensionTower::Adj { min_poly_q, .. } => min_poly_q,
            ExtensionTower::Base => {
                return Err(EvalError::TypeError("expected extension field"));
            }
        };
        Ok(min_poly_q
            .iter()
            .map(|r| parent.embed_rational(r))
            .collect())
    }

    fn unflatten_layer_blocks(&self, v: &CoordsQ) -> Result<Vec<CoordsQ>, EvalError> {
        let parent = self
            .parent_field
            .as_ref()
            .ok_or(EvalError::TypeError("tower field without parent"))?;
        let pd = parent.dimension();
        let e = self.layer_ext_degree()?;
        let v = pad_to_len(v, self.dimension());
        let mut blocks = Vec::with_capacity(e);
        for j in 0..e {
            blocks.push(v[j * pd..(j + 1) * pd].to_vec());
        }
        // Coords are u^0..u^{e-1}; poly ops use leading coeff first.
        blocks.reverse();
        Ok(blocks)
    }

    fn flatten_layer_element(&self, blocks: &[CoordsQ]) -> Result<CoordsQ, EvalError> {
        let parent = self
            .parent_field
            .as_ref()
            .ok_or(EvalError::TypeError("tower field without parent"))?;
        let mut blocks = blocks.to_vec();
        blocks.reverse();
        Ok(pad_to_len(
            &flatten_layer_blocks(&blocks, parent.dimension()),
            self.dimension(),
        ))
    }

    pub fn element_add(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_add_tower(a, b)
        } else {
            self.element_add_primitive(a, b)
        }
    }

    fn element_add_tower(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        let parent = self
            .parent_field
            .as_ref()
            .ok_or(EvalError::TypeError("tower field without parent"))?;
        let ring = parent_coeff_ring!(parent);
        let blocks_a = self.unflatten_layer_blocks(a)?;
        let blocks_b = self.unflatten_layer_blocks(b)?;
        let sum = poly_add_with_coeffs_in_field(&blocks_a, &blocks_b, &ring)?;
        self.flatten_layer_element(&sum)
    }

    fn element_add_primitive(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        self.ensure_same_field_len(b)?;
        if self.is_base() {
            return Ok(vec![pad_to_len(a, 1)[0].clone() + pad_to_len(b, 1)[0].clone()]);
        }
        Ok(pad_to_len(
            &poly_reduce(&poly_add(a, b), &self.min_poly_over_q),
            self.dimension(),
        ))
    }

    pub fn element_sub(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_sub_tower(a, b)
        } else {
            self.element_sub_primitive(a, b)
        }
    }

    fn element_sub_tower(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        let parent = self
            .parent_field
            .as_ref()
            .ok_or(EvalError::TypeError("tower field without parent"))?;
        let ring = parent_coeff_ring!(parent);
        let blocks_a = self.unflatten_layer_blocks(a)?;
        let blocks_b = self.unflatten_layer_blocks(b)?;
        let mut diff = Vec::with_capacity(blocks_a.len().max(blocks_b.len()));
        let d = blocks_a.len().max(blocks_b.len());
        for i in 0..d {
            let ai = blocks_a.get(i).unwrap_or(&ring.zero);
            let bi = blocks_b.get(i).unwrap_or(&ring.zero);
            diff.push((ring.sub)(ai, bi)?);
        }
        self.flatten_layer_element(&diff)
    }

    fn element_sub_primitive(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        self.ensure_same_field_len(b)?;
        if self.is_base() {
            return Ok(vec![pad_to_len(a, 1)[0].clone() - pad_to_len(b, 1)[0].clone()]);
        }
        Ok(pad_to_len(
            &poly_reduce(&poly_sub(a, b), &self.min_poly_over_q),
            self.dimension(),
        ))
    }

    pub fn element_neg(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_neg_tower(a)
        } else {
            self.element_neg_primitive(a)
        }
    }

    fn element_neg_tower(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        let parent = self
            .parent_field
            .as_ref()
            .ok_or(EvalError::TypeError("tower field without parent"))?;
        let ring = parent_coeff_ring!(parent);
        let blocks = self.unflatten_layer_blocks(a)?;
        let neg = poly_neg_with_coeffs_in_field(&blocks, &ring)?;
        self.flatten_layer_element(&neg)
    }

    fn element_neg_primitive(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        if self.is_base() {
            return Ok(vec![-pad_to_len(a, 1)[0].clone()]);
        }
        Ok(pad_to_len(
            &poly_reduce(&poly_neg(a), &self.min_poly_over_q),
            self.dimension(),
        ))
    }

    pub fn element_mul(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_mul_tower(a, b)
        } else {
            self.element_mul_primitive(a, b)
        }
    }

    fn element_mul_tower(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        let parent = self
            .parent_field
            .as_ref()
            .ok_or(EvalError::TypeError("tower field without parent"))?;
        let ring = parent_coeff_ring!(parent);
        let modulus = self.layer_minpoly_parent_coeffs(parent)?;
        let blocks_a = self.unflatten_layer_blocks(a)?;
        let blocks_b = self.unflatten_layer_blocks(b)?;
        let prod = poly_mul_with_coeffs_in_field(&blocks_a, &blocks_b, &ring)?;
        let reduced = poly_reduce_with_coeffs_in_field(&prod, &modulus, &ring)?;
        self.flatten_layer_element(&reduced)
    }

    fn element_mul_primitive(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        self.ensure_same_field_len(b)?;
        if self.is_base() {
            return Ok(vec![pad_to_len(a, 1)[0].clone() * pad_to_len(b, 1)[0].clone()]);
        }
        Ok(pad_to_len(
            &poly_reduce(&poly_mul(a, b), &self.min_poly_over_q),
            self.dimension(),
        ))
    }

    pub fn element_inv(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_inv_tower(a)
        } else {
            self.element_inv_primitive(a)
        }
    }

    fn element_inv_tower(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        let parent = self
            .parent_field
            .as_ref()
            .ok_or(EvalError::TypeError("tower field without parent"))?;
        if coords_all_zero(a) {
            return Err(EvalError::DivisionByZero);
        }
        let ring = parent_coeff_ring!(parent);
        let modulus = self.layer_minpoly_parent_coeffs(parent)?;
        let blocks = self.unflatten_layer_blocks(a)?;
        let inv = poly_inv_mod_with_coeffs_in_field(&blocks, &modulus, &ring)?;
        self.flatten_layer_element(&inv)
    }

    fn element_inv_primitive(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        if self.is_base() {
            let x = pad_to_len(a, 1)[0].clone();
            if x.is_zero() {
                return Err(EvalError::DivisionByZero);
            }
            return Ok(vec![Ratio::one() / x]);
        }
        if coords_all_zero(a) {
            return Err(EvalError::DivisionByZero);
        }
        let inv = poly_inv_mod(a, &self.min_poly_over_q)?;
        Ok(pad_to_len(&inv, self.dimension()))
    }

    pub fn element_eq_mod(&self, a: &CoordsQ, b: &CoordsQ) -> Result<bool, EvalError> {
        let diff = self.element_sub(a, b)?;
        Ok(coords_all_zero(&diff))
    }

    pub fn element_is_zero(&self, a: &CoordsQ) -> bool {
        a.iter().all(|c| c.is_zero())
    }

    pub fn element_is_one(&self, a: &CoordsQ) -> bool {
        pad_to_len(a, self.dimension()) == self.one_coords()
    }

    fn ensure_same_field_len(&self, a: &CoordsQ) -> Result<(), EvalError> {
        if a.len() > self.dimension() {
            return Err(EvalError::TypeError("coords longer than field dimension"));
        }
        Ok(())
    }

    /// Minimal common extension of two ℚ-described fields, with linear embeddings.
    pub fn common_over_q(
        a: &Arc<ExtensionField>,
        b: &Arc<ExtensionField>,
    ) -> Result<Arc<CommonFieldPair>, EvalError> {
        if Arc::ptr_eq(a, b) || **a == **b {
            return Ok(CommonFieldPair::identity(a));
        }
        field_registry().common_cached(a, b)
    }

    /// Pick the embedding whose [`FieldEmbedding::source`] matches `source`.
    ///
    /// `CommonFieldPair::embed_a` / `embed_b` are keyed by sorted field id in the cache,
    /// not by caller operand order — always use this when applying a cached pair.
    pub fn embedding_for<'a>(
        source: &Arc<ExtensionField>,
        pair: &'a CommonFieldPair,
    ) -> Result<&'a FieldEmbedding, EvalError> {
        if Arc::ptr_eq(source, &pair.embed_a.source) || **source == *pair.embed_a.source {
            Ok(&pair.embed_a)
        } else if Arc::ptr_eq(source, &pair.embed_b.source) || **source == *pair.embed_b.source {
            Ok(&pair.embed_b)
        } else {
            Err(EvalError::TypeError("field not in common pair"))
        }
    }

    /// Align two field elements into a common ambient field (S0 unified entry point).
    ///
    /// T2: when one field is a subfield of the other (tower adjoin chain), embed along the
    /// chain and **do not** call [`Self::common_over_q`] / `common_cache`.
    pub fn align_elements(
        a_field: &Arc<ExtensionField>,
        a_coords: &CoordsQ,
        b_field: &Arc<ExtensionField>,
        b_coords: &CoordsQ,
    ) -> Result<AlignedElements, EvalError> {
        if Arc::ptr_eq(a_field, b_field) || **a_field == **b_field {
            return Ok(AlignedElements {
                field: Arc::clone(a_field),
                left: pad_to_len(a_coords, a_field.dimension()),
                right: pad_to_len(b_coords, b_field.dimension()),
            });
        }
        if Self::is_subfield_of(a_field, b_field) {
            let emb_a = Self::try_subfield_embedding(a_field, b_field)?
                .ok_or(EvalError::TypeError("subfield embedding missing"))?;
            return Ok(AlignedElements {
                field: Arc::clone(b_field),
                left: emb_a.apply(a_coords),
                right: pad_to_len(b_coords, b_field.dimension()),
            });
        }
        if Self::is_subfield_of(b_field, a_field) {
            let emb_b = Self::try_subfield_embedding(b_field, a_field)?
                .ok_or(EvalError::TypeError("subfield embedding missing"))?;
            return Ok(AlignedElements {
                field: Arc::clone(a_field),
                left: pad_to_len(a_coords, a_field.dimension()),
                right: emb_b.apply(b_coords),
            });
        }
        let common = Self::common_over_q(a_field, b_field)?;
        let left = Self::embedding_for(a_field, &common)?.apply(a_coords);
        let right = Self::embedding_for(b_field, &common)?.apply(b_coords);
        Ok(AlignedElements {
            field: Arc::clone(&common.field),
            left,
            right,
        })
    }
}

/// Linear embedding from a subfield into a common superfield.
#[derive(Clone, Debug)]
pub struct FieldEmbedding {
    pub source: Arc<ExtensionField>,
    pub target: Arc<ExtensionField>,
    /// `out[i] = Σⱼ matrix[i][j] · in[j]`
    pub matrix: Vec<Vec<Ratio<BigInt>>>,
}

impl FieldEmbedding {
    pub fn apply(&self, coords: &CoordsQ) -> CoordsQ {
        let src_dim = self.source.dimension();
        let tgt_dim = self.target.dimension();
        let v = pad_to_len(coords, src_dim);
        let mut out = vec![Ratio::zero(); tgt_dim];
        for i in 0..tgt_dim {
            if i >= self.matrix.len() {
                continue;
            }
            for j in 0..src_dim {
                if j < self.matrix[i].len() {
                    out[i] += self.matrix[i][j].clone() * v[j].clone();
                }
            }
        }
        out
    }
}

/// Result of merging two extension fields into a common ambient field.
#[derive(Clone, Debug)]
pub struct CommonFieldPair {
    pub field: Arc<ExtensionField>,
    pub embed_a: FieldEmbedding,
    pub embed_b: FieldEmbedding,
}

impl CommonFieldPair {
    pub fn identity(field: &Arc<ExtensionField>) -> Arc<Self> {
        let dim = field.dimension();
        let id = super::field_arith::identity_matrix(dim);
        Arc::new(Self {
            field: Arc::clone(field),
            embed_a: FieldEmbedding {
                source: Arc::clone(field),
                target: Arc::clone(field),
                matrix: id.clone(),
            },
            embed_b: FieldEmbedding {
                source: Arc::clone(field),
                target: Arc::clone(field),
                matrix: id,
            },
        })
    }
}

/// Two extension elements aligned into a common ambient field.
#[derive(Clone, Debug)]
pub struct AlignedElements {
    pub field: Arc<ExtensionField>,
    pub left: CoordsQ,
    pub right: CoordsQ,
}

struct FieldRegistry {
    by_adjoin: Mutex<HashMap<(u64, Vec<u8>), Arc<ExtensionField>>>,
    by_min_poly: Mutex<HashMap<Vec<u8>, Arc<ExtensionField>>>,
    common_cache: Mutex<HashMap<(u64, u64), Arc<CommonFieldPair>>>,
}

impl FieldRegistry {
    fn new() -> Self {
        Self {
            by_adjoin: Mutex::new(HashMap::new()),
            by_min_poly: Mutex::new(HashMap::new()),
            common_cache: Mutex::new(HashMap::new()),
        }
    }

    fn min_poly_key(p: &[Ratio<BigInt>]) -> Vec<u8> {
        let mut key = Vec::new();
        for c in p {
            let n = c.numer().to_signed_bytes_le();
            let d = c.denom().to_signed_bytes_le();
            key.extend_from_slice(&(n.len() as u32).to_le_bytes());
            key.extend_from_slice(&n);
            key.extend_from_slice(&(d.len() as u32).to_le_bytes());
            key.extend_from_slice(&d);
        }
        key
    }

    fn register_adjoin(
        &self,
        parent: &Arc<ExtensionField>,
        min_poly_q: &CoordsQ,
    ) -> Result<Arc<ExtensionField>, EvalError> {
        let adjoin_key = (parent.id(), Self::min_poly_key(min_poly_q));
        if let Some(hit) = self.by_adjoin.lock().unwrap().get(&adjoin_key) {
            return Ok(Arc::clone(hit));
        }
        if parent.is_base() {
            if let Some(hit) = self
                .by_min_poly
                .lock()
                .unwrap()
                .get(&Self::min_poly_key(min_poly_q))
            {
                self.by_adjoin
                    .lock()
                    .unwrap()
                    .insert(adjoin_key, Arc::clone(hit));
                return Ok(Arc::clone(hit));
            }
        }
        let ext_degree = poly_degree(min_poly_q);
        let (min_poly_over_q, primitive_k) = if parent.is_base() {
            (min_poly_q.clone(), None)
        } else {
            compose_min_poly_over_q(parent, min_poly_q)?
        };
        let field = Arc::new(ExtensionField {
            id: next_field_id(),
            tower: Arc::new(ExtensionTower::Adj {
                parent: Arc::clone(parent.tower()),
                min_poly_q: min_poly_q.clone(),
                ext_degree,
                primitive_k,
            }),
            min_poly_over_q,
            parent_field: Some(Arc::clone(parent)),
        });
        self.by_adjoin
            .lock()
            .unwrap()
            .insert(adjoin_key, Arc::clone(&field));
        if parent.is_base() {
            self.by_min_poly
                .lock()
                .unwrap()
                .insert(Self::min_poly_key(min_poly_q), Arc::clone(&field));
        }
        Ok(field)
    }

    fn register_common(
        &self,
        min_poly_q: CoordsQ,
        embed_a: FieldEmbedding,
        embed_b: FieldEmbedding,
    ) -> Arc<CommonFieldPair> {
        let key = Self::min_poly_key(&min_poly_q);
        let field = if let Some(hit) = self.by_min_poly.lock().unwrap().get(&key) {
            Arc::clone(hit)
        } else {
            let d = poly_degree(&min_poly_q);
            let f = Arc::new(ExtensionField {
                id: next_field_id(),
                tower: Arc::new(ExtensionTower::Adj {
                    parent: Arc::new(ExtensionTower::Base),
                    min_poly_q: min_poly_q.clone(),
                    ext_degree: d,
                    primitive_k: None,
                }),
                min_poly_over_q: min_poly_q,
                parent_field: None,
            });
            self.by_min_poly
                .lock()
                .unwrap()
                .insert(Self::min_poly_key(f.min_poly_over_q()), Arc::clone(&f));
            f
        };
        Arc::new(CommonFieldPair {
            field,
            embed_a,
            embed_b,
        })
    }

    fn common_cached(
        &self,
        a: &Arc<ExtensionField>,
        b: &Arc<ExtensionField>,
    ) -> Result<Arc<CommonFieldPair>, EvalError> {
        let key = if a.id() <= b.id() {
            (a.id(), b.id())
        } else {
            (b.id(), a.id())
        };
        if let Some(hit) = self.common_cache.lock().unwrap().get(&key) {
            return Ok(Arc::clone(hit));
        }
        let (fa, fb) = if a.id() <= b.id() {
            (a, b)
        } else {
            (b, a)
        };
        let pair = compute_common_over_q(fa, fb)?;
        self.common_cache.lock().unwrap().insert(key, Arc::clone(&pair));
        Ok(pair)
    }
}

fn field_registry() -> &'static FieldRegistry {
    static REG: OnceLock<FieldRegistry> = OnceLock::new();
    REG.get_or_init(FieldRegistry::new)
}

#[cfg(test)]
pub(crate) fn common_cache_len_for_test() -> usize {
    field_registry().common_cache.lock().unwrap().len()
}

/// Second `Arc` handle with the same tower / minpoly but a distinct id (tests only).
#[cfg(test)]
pub(crate) fn duplicate_field_arc_for_test(f: &Arc<ExtensionField>) -> Arc<ExtensionField> {
    Arc::new(ExtensionField {
        id: next_field_id(),
        tower: Arc::clone(f.tower()),
        min_poly_over_q: f.min_poly_over_q().to_vec(),
        parent_field: f.parent_field.as_ref().map(Arc::clone),
    })
}

fn layer_minpoly_rational_constants(min_poly_q: &[Ratio<BigInt>]) -> bool {
    min_poly_q.iter().all(|c| c.denom().is_one())
}

fn compose_min_poly_over_q(
    parent: &ExtensionField,
    layer_min_poly: &CoordsQ,
) -> Result<(CoordsQ, Option<i64>), EvalError> {
    let na = parent.dimension();
    let nb = poly_degree(layer_min_poly);
    let expected = na * nb;
    for k in 1i64..=12 {
        match common_primitive_sum(
            parent.min_poly_over_q(),
            layer_min_poly,
            na,
            nb,
            k,
        ) {
            Ok((min_g, _, _)) if poly_degree(&min_g) == expected => {
                return Ok((min_g, Some(k)));
            }
            Ok(_) | Err(EvalError::TypeError(_)) => continue,
            Err(e) => return Err(e),
        }
    }
    Err(EvalError::NotImplemented("ExtensionField::adjoin compose minpoly"))
}

fn rational_subfield_embedding(
    sub: &Arc<ExtensionField>,
    sup: &Arc<ExtensionField>,
) -> FieldEmbedding {
    let dim = sup.dimension();
    let mut m = vec![vec![Ratio::zero(); 1]; dim];
    if dim > 0 {
        m[dim - 1][0] = Ratio::one();
    }
    FieldEmbedding {
        source: Arc::clone(sub),
        target: Arc::clone(sup),
        matrix: m,
    }
}

fn flatten_layer_blocks(blocks: &[CoordsQ], parent_dim: usize) -> CoordsQ {
    let mut out = Vec::with_capacity(blocks.len() * parent_dim);
    for block in blocks {
        out.extend(pad_to_len(block, parent_dim).iter().cloned());
    }
    out
}

fn direct_adjoin_parent_embedding(
    parent: &Arc<ExtensionField>,
    child: &Arc<ExtensionField>,
) -> Result<FieldEmbedding, EvalError> {
    if child.parent_field.as_ref().map(|p| p.id()) != Some(parent.id()) {
        return Err(EvalError::TypeError("child parent mismatch"));
    }
    if parent.is_base() {
        return Ok(rational_subfield_embedding(parent, child));
    }
    let pd = parent.dimension();
    let cd = child.dimension();
    let mut matrix = vec![vec![Ratio::zero(); pd]; cd];
    for j in 0..pd {
        matrix[j][j] = Ratio::one();
    }
    Ok(FieldEmbedding {
        source: Arc::clone(parent),
        target: Arc::clone(child),
        matrix,
    })
}

fn compose_field_embeddings(
    inner: &FieldEmbedding,
    outer: &FieldEmbedding,
) -> Result<FieldEmbedding, EvalError> {
    debug_assert!(Arc::ptr_eq(&inner.target, &outer.source) || inner.target == outer.source);
    let src_dim = inner.source.dimension();
    let mid_dim = inner.target.dimension();
    let tgt_dim = outer.target.dimension();
    if inner.matrix.len() != mid_dim
        || outer.matrix.len() != tgt_dim
        || inner
            .matrix
            .iter()
            .any(|row| row.len() != src_dim)
        || outer.matrix.iter().any(|row| row.len() != mid_dim)
    {
        return Err(EvalError::TypeError("embedding compose dimension mismatch"));
    }
    let mut matrix = vec![vec![Ratio::zero(); src_dim]; tgt_dim];
    for i in 0..tgt_dim {
        for j in 0..src_dim {
            for k in 0..mid_dim {
                matrix[i][j] += outer.matrix[i][k].clone() * inner.matrix[k][j].clone();
            }
        }
    }
    Ok(FieldEmbedding {
        source: Arc::clone(&inner.source),
        target: Arc::clone(&outer.target),
        matrix,
    })
}

fn compute_common_over_q(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    if a.is_base() && !b.is_base() {
        return Ok(embed_rationals_into(b, a, b));
    }
    if b.is_base() && !a.is_base() {
        return Ok(embed_rationals_into(a, b, a));
    }
    let ma = a.min_poly_over_q();
    let mb = b.min_poly_over_q();
    let na = poly_degree(ma);
    let nb = poly_degree(mb);
    for k in 1i64..=12 {
        match common_primitive_sum(ma, mb, na, nb, k) {
            Ok((min_g, mat_a, mat_b)) => {
                let pair = field_registry().register_common(
                    min_g,
                    FieldEmbedding {
                        source: Arc::clone(a),
                        target: Arc::clone(a),
                        matrix: mat_a.clone(),
                    },
                    FieldEmbedding {
                        source: Arc::clone(b),
                        target: Arc::clone(b),
                        matrix: mat_b.clone(),
                    },
                );
                let patched = Arc::new(CommonFieldPair {
                    field: Arc::clone(&pair.field),
                    embed_a: FieldEmbedding {
                        source: Arc::clone(a),
                        target: Arc::clone(&pair.field),
                        matrix: mat_a,
                    },
                    embed_b: FieldEmbedding {
                        source: Arc::clone(b),
                        target: Arc::clone(&pair.field),
                        matrix: mat_b,
                    },
                });
                return Ok(patched);
            }
            Err(EvalError::TypeError(_)) => continue,
            Err(e) => return Err(e),
        }
    }
    Err(EvalError::NotImplemented("ExtensionField::common"))
}

/// ℚ ↪ K: rationals embed as the constant term of K.
fn embed_rationals_into(
    ext: &Arc<ExtensionField>,
    source_a: &Arc<ExtensionField>,
    source_b: &Arc<ExtensionField>,
) -> Arc<CommonFieldPair> {
    let dim = ext.dimension();
    let rational_embed = || {
        let mut m = vec![vec![Ratio::zero(); 1]; dim];
        if dim > 0 {
            m[dim - 1][0] = Ratio::one();
        }
        m
    };
    Arc::new(CommonFieldPair {
        field: Arc::clone(ext),
        embed_a: FieldEmbedding {
            source: Arc::clone(source_a),
            target: Arc::clone(ext),
            matrix: if source_a.is_base() {
                rational_embed()
            } else {
                super::field_arith::identity_matrix(dim)
            },
        },
        embed_b: FieldEmbedding {
            source: Arc::clone(source_b),
            target: Arc::clone(ext),
            matrix: if source_b.is_base() {
                rational_embed()
            } else {
                super::field_arith::identity_matrix(dim)
            },
        },
    })
}

fn common_primitive_sum(
    ma: &[Ratio<BigInt>],
    mb: &[Ratio<BigInt>],
    na: usize,
    nb: usize,
    k: i64,
) -> Result<(CoordsQ, Vec<Vec<Ratio<BigInt>>>, Vec<Vec<Ratio<BigInt>>>), EvalError> {
    let gen_a = generator_coords(na);
    let gen_b = generator_coords(nb);
    let mat_a = kron_left(&mult_matrix_of_element(&gen_a, ma), nb);
    let mat_b = kron_left(&mult_matrix_of_element(&gen_b, mb), na);
    let k_rat = Ratio::from_integer(BigInt::from(k));
    let dim = na * nb;
    let mut mat_theta = mat_a.clone();
    for i in 0..dim {
        for j in 0..dim {
            mat_theta[i][j] += k_rat.clone() * mat_b[i][j].clone();
        }
    }
    let min_g = char_poly_matrix(&mat_theta);
    if poly_degree(&min_g) != dim {
        return Err(EvalError::TypeError("common degree mismatch"));
    }
    let embed_a = embedding_matrix_from_theta(&mat_theta, na, dim)?;
    let embed_b = embedding_matrix_from_theta_block(&mat_theta, nb, na, dim)?;
    Ok((min_g, embed_a, embed_b))
}

/// Embed each basis vector of a degree-`source_dim` field into the `dim`-dimensional θ-basis.
fn embedding_matrix_from_theta(
    mat_theta: &[Vec<Ratio<BigInt>>],
    source_dim: usize,
    dim: usize,
) -> Result<Vec<Vec<Ratio<BigInt>>>, EvalError> {
    let mut matrix = vec![vec![Ratio::zero(); source_dim]; dim];
    for j in 0..source_dim {
        let mut raw = vec![Ratio::zero(); dim];
        raw[j] = Ratio::one();
        let col = embed_in_gamma_vector(&raw, mat_theta, dim)?;
        for i in 0..dim {
            matrix[i][j] = col[i].clone();
        }
    }
    Ok(matrix)
}

/// Embed block-placed basis vectors (β^j at offset j·na) into the θ-basis.
fn embedding_matrix_from_theta_block(
    mat_theta: &[Vec<Ratio<BigInt>>],
    source_dim: usize,
    na: usize,
    dim: usize,
) -> Result<Vec<Vec<Ratio<BigInt>>>, EvalError> {
    let mut matrix = vec![vec![Ratio::zero(); source_dim]; dim];
    for j in 0..source_dim {
        let mut raw = vec![Ratio::zero(); dim];
        raw[j * na] = Ratio::one();
        let col = embed_in_gamma_vector(&raw, mat_theta, dim)?;
        for i in 0..dim {
            matrix[i][j] = col[i].clone();
        }
    }
    Ok(matrix)
}

fn embed_in_gamma_vector(
    v: &[Ratio<BigInt>],
    mat_theta: &[Vec<Ratio<BigInt>>],
    dim: usize,
) -> Result<CoordsQ, EvalError> {
    let mut pow = super::field_arith::identity_matrix(dim);
    let mut basis = vec![vec![Ratio::zero(); dim]; dim];
    for i in 0..dim {
        basis[i] = pow[i].clone();
        pow = mat_mul(&pow, mat_theta);
    }
    let mut out = vec![Ratio::zero(); dim];
    for i in 0..dim {
        for j in 0..dim {
            if i < v.len() {
                out[j] += v[i].clone() * basis[i][j].clone();
            }
        }
    }
    Ok(out)
}

/// Embed coords of `source` into `target` via `embedding`.
pub fn embed_coords(embedding: &FieldEmbedding, coords: &CoordsQ) -> CoordsQ {
    embedding.apply(coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::test_fixtures::{
        common_cache_len, k1_adjoin_cbrt2, k1_adjoin_sqrt2, minpoly_u2_minus,
    };

    #[test]
    fn rational_field_dimension_one() {
        let q = ExtensionField::rational();
        assert_eq!(q.dimension(), 1);
        assert!(q.element_is_one(&q.one_coords()));
    }

    #[test]
    fn adjoin_sqrt2_dimension_two() {
        let k = k1_adjoin_sqrt2();
        assert_eq!(k.dimension(), 2);
        let alpha = k.generator_coords();
        let sq = k.element_mul(&alpha, &alpha).unwrap();
        let two = pad_to_len(&[Ratio::from_integer(2.into())], 2);
        assert!(k.element_eq_mod(&sq, &two).unwrap());
    }

    #[test]
    #[ignore = "primitive-element common(√2,∛2) char poly is slow; see ext_tower perf follow-up"]
    fn common_sqrt2_cbrt2_has_degree_six() {
        let a = k1_adjoin_sqrt2();
        let b = k1_adjoin_cbrt2();
        let common = ExtensionField::common_over_q(&a, &b).unwrap();
        assert_eq!(common.field.dimension(), 6);
        let ea = common.embed_a.apply(&a.generator_coords());
        let eb = common.embed_b.apply(&b.generator_coords());
        let sum = common.field.element_add(&ea, &eb).unwrap();
        assert!(!common.field.element_is_zero(&sum));
    }

    #[test]
    fn common_cache_identity_is_fast() {
        let a = k1_adjoin_sqrt2();
        let c1 = ExtensionField::common_over_q(&a, &a).unwrap();
        let c2 = ExtensionField::common_over_q(&a, &a).unwrap();
        assert_eq!(c1.field.id(), c2.field.id());
    }

    #[test]
    #[ignore = "primitive-element common(√2,∛2) char poly is slow; see ext_tower perf follow-up"]
    fn common_cache_hits_same_pair() {
        let a = k1_adjoin_sqrt2();
        let b = k1_adjoin_cbrt2();
        let c1 = ExtensionField::common_over_q(&a, &b).unwrap();
        let c2 = ExtensionField::common_over_q(&b, &a).unwrap();
        assert_eq!(c1.field.id(), c2.field.id());
    }

    #[test]
    fn field_registry_dedup_same_minpoly() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = k1_adjoin_sqrt2();
        assert_eq!(k1.id(), k2.id());
    }

    #[test]
    fn embedding_for_matches_source_not_operand_order() {
        let sqrt2 = k1_adjoin_sqrt2();
        let cbrt2 = k1_adjoin_cbrt2();
        let (low, high) = if sqrt2.id() <= cbrt2.id() {
            (&sqrt2, &cbrt2)
        } else {
            (&cbrt2, &sqrt2)
        };
        let id = crate::algebra::field_arith::identity_matrix;
        let pair = CommonFieldPair {
            field: Arc::clone(high),
            embed_a: FieldEmbedding {
                source: Arc::clone(low),
                target: Arc::clone(high),
                matrix: id(low.dimension()),
            },
            embed_b: FieldEmbedding {
                source: Arc::clone(high),
                target: Arc::clone(high),
                matrix: id(high.dimension()),
            },
        };
        assert_eq!(
            ExtensionField::embedding_for(low, &pair).unwrap().source.id(),
            low.id()
        );
        assert_eq!(
            ExtensionField::embedding_for(high, &pair).unwrap().source.id(),
            high.id()
        );
    }

    #[test]
    fn align_elements_reverse_order_after_cache_warm() {
        let q = ExtensionField::rational();
        let sqrt2 = k1_adjoin_sqrt2();
        let _ = ExtensionField::common_over_q(&q, &sqrt2).unwrap();

        let one_q = q.one_coords();
        let alpha = sqrt2.generator_coords();

        let qr = ExtensionField::align_elements(&q, &one_q, &sqrt2, &alpha).unwrap();
        let rq = ExtensionField::align_elements(&sqrt2, &alpha, &q, &one_q).unwrap();

        assert_eq!(qr.field.id(), rq.field.id());
        assert_eq!(qr.field.id(), sqrt2.id());
        assert!(qr.field.element_eq_mod(&qr.left, &rq.right).unwrap());
        assert!(qr.field.element_eq_mod(&qr.right, &rq.left).unwrap());

        let sum = qr.field.element_add(&qr.left, &qr.right).unwrap();
        assert!(!qr.field.element_is_zero(&sum));
    }

    #[test]
    fn t1a_adjoin_base_sqrt2_matches_legacy() {
        let legacy = ExtensionField::adjoin_irreducible_over_q(minpoly_u2_minus(-2)).unwrap();
        let via_parent =
            ExtensionField::adjoin_irreducible(&ExtensionField::rational(), minpoly_u2_minus(-2))
                .unwrap();
        assert_eq!(*legacy, *via_parent);
        assert_eq!(legacy.dimension(), 2);
        assert!(legacy.parent_field().is_some());
        assert!(legacy.parent_field().unwrap().is_base());
    }

    #[test]
    fn t1b_adjoin_k1_u2_minus_3_has_dimension_four() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        assert_eq!(k2.dimension(), 4);
        assert_eq!(k2.parent_field().map(|p| p.id()), Some(k1.id()));
        assert!(ExtensionField::is_subfield_of(&k1, &k2));
        assert!(!ExtensionField::is_subfield_of(&k2, &k1));
    }

    #[test]
    fn t1_try_subfield_embedding_k1_into_k2() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let emb = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .expect("K1 embeds in K2");
        assert_eq!(emb.source.id(), k1.id());
        assert_eq!(emb.target.id(), k2.id());
        assert_eq!(emb.matrix.len(), k2.dimension());
        assert_eq!(emb.matrix[0].len(), k1.dimension());
        let img = emb.apply(&k1.generator_coords());
        assert_eq!(img.len(), k2.dimension());
    }

    #[test]
    fn t2_align_subfield_does_not_grow_common_cache() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let before = common_cache_len();
        let alpha = k1.generator_coords();
        let alpha_in_k2 = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .unwrap()
            .apply(&alpha);
        let _ = ExtensionField::align_elements(&k1, &alpha, &k2, &alpha_in_k2).unwrap();
        assert_eq!(common_cache_len(), before);
    }

    #[test]
    fn t2_algext_add_same_generator_in_superfield_no_common_cache() {
        use crate::algebra::alg_ext::AlgExtData;
        use crate::algebra::field_arith::coords_to_expr;

        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let alpha_k2 = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .unwrap()
            .apply(&k1.generator_coords());
        let a = AlgExtData::from_field_coords(
            Arc::clone(&k2),
            coords_to_expr(&alpha_k2).unwrap(),
        )
        .unwrap();
        let b = AlgExtData::from_field_coords(Arc::clone(&k2), coords_to_expr(&alpha_k2).unwrap()).unwrap();
        let before = common_cache_len();
        let sum = a.add(&b).unwrap();
        assert_eq!(common_cache_len(), before);
        assert_eq!(sum.field.id(), a.field.id());
        let two_b = b.add(&b).unwrap();
        assert!(sum.eq_mod(&two_b).unwrap());
    }

    #[test]
    fn t3_k1_sqrt2_squared_is_two() {
        let k1 = k1_adjoin_sqrt2();
        let alpha = k1.generator_coords();
        let prod = k1.element_mul(&alpha, &alpha).unwrap();
        let two = k1.embed_rational(&Ratio::from_integer(2.into()));
        assert!(k1.element_eq_mod(&prod, &two).unwrap());
    }

    #[test]
    fn t3_k2_sqrt2_beta_times_beta_is_three_sqrt2() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let emb = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .expect("K1 embeds in K2");
        let sqrt2 = emb.apply(&k1.generator_coords());
        let beta = k2.generator_coords();
        let sqrt2_beta = k2.element_mul(&sqrt2, &beta).unwrap();
        let prod = k2.element_mul(&sqrt2_beta, &beta).unwrap();
        let three = k2.embed_rational(&Ratio::from_integer(3.into()));
        let expected = k2.element_mul(&sqrt2, &three).unwrap();
        assert!(k2.element_eq_mod(&prod, &expected).unwrap());
    }
}
