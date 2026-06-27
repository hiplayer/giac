//! Tower of algebraic extensions K₀=ℚ → K₁ → … → Kₙ and cross-field `common`.
//!
//! Normative model: [GIAC-algext-adoption.md](../../../../.doc/issues/GIAC-algext-adoption.md) §8.2.
//!
//! **Dev note — primitive vs tower + `parent_coeff_ring`:** see
//! [giac-tower-common-math.md](../../../../.doc/giac-tower-common-math.md) §1.1.
//! Summary: `element_*` dispatches primitive (single `Adj{parent:ℚ}`) vs tower
//! (`min_poly_parent_blocks` or nontrivial `parent_field`). Tower layer ops unflatten
//! coords into parent-field coefficients; **`parent_coeff_ring!` must call `element_*`
//! on the parent, never `element_*_primitive`** — otherwise adjoin on a composite parent
//! (e.g. dim-6 after resolvent √Δ) hits `primitive modulus: expected simple over Q`.
//!
//! Phase 0 scope: coefficients in ℚ only; each [`ExtensionField`] is described by a
//! primitive/minimal polynomial over ℚ (single adjoin or `common` composite).
//!
//! T1+: nontrivial adjoins record a true [`ExtensionTower::Adj`] parent chain and
//! `(parent_id, min_poly)` registry keys; see [GIAC-lazy-common-tower-plan.md] T1.
//!
//! ## Square roots in K (F5)
//!
//! | API | Behavior |
//! |-----|----------|
//! | [`ExtensionField::try_square_root_in_field`] | Structural probe S0–S6; **no** adjoin |
//! | [`ExtensionField::try_square_root_in_field_shallow`] | S1–S3 only (Euler hot path) |
//! | [`super::field_session::FieldSession::sqrt_in_field`] | try → adjoin fallback |
//! | [`super::alg_ext::algext_square_roots`] | **Always** adjoin; roots pipeline must not call |
//!
//! **F5 probe order:** S1 layer ±g → S2 basis ±eᵢ → S3 subfield descent → S4 top quadratic
//! closed form → S5 layer-gen products → S6 dim≤6 pairwise ±1 fallback.
//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc, FuncKind};

use super::alg_ext::AlgExtData;
use super::field_arith::{
    char_poly_matrix, coords_all_zero,
    generator_coords, kron_left, mat_mul, mult_matrix_of_adjoin_generator,
    mult_matrix_of_element, pad_to_len, poly_add, poly_degree, poly_inv_mod, poly_mul, poly_neg,
    poly_reduce, poly_sub, poly_add_with_coeffs_in_field, poly_inv_mod_with_coeffs_in_field,
    poly_mul_with_coeffs_in_field, poly_neg_with_coeffs_in_field,
    poly_reduce_with_coeffs_in_field, ParentCoeffRing,
    coords_to_expr, CoordsQ,
};

#[path = "common_minimal.rs"]
mod common_minimal;

static FIELD_ID: AtomicU64 = AtomicU64::new(1);

/// Upper bound for θ = α + k·β primitive-element search (flatten / compose / U1a).
pub(crate) const PRIMITIVE_K_SEARCH_MAX: i64 = 64;

// **Pipeline private** — `next_field_id`
fn next_field_id() -> u64 {
    FIELD_ID.fetch_add(1, Ordering::Relaxed)
}

// **Pipeline private** — one layer-minpoly coefficient (parent-field element) as `Expr`
fn layer_minpoly_block_to_expr(block: &CoordsQ) -> Result<ExprArc, EvalError> {
    use super::field_arith::canonical_poly1_expr;
    let coeffs = canonical_poly1_expr(&coords_to_expr(block)?);
    Ok(match coeffs.as_slice() {
        [c] => Arc::clone(c),
        _ => Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(coeffs))],
        )),
    })
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
        /// Rational-constant layer minpoly (high-first); empty when [`min_poly_parent_blocks`] is set.
        min_poly_q: CoordsQ,
        /// T3+: layer minpoly with parent operational coefficients (high-first blocks).
        min_poly_parent_blocks: Option<Vec<CoordsQ>>,
        ext_degree: usize,
        /// When this layer adjoins over a nontrivial parent, `k` in θ = α + k·β used to
        /// build ℚ-flatten primitive poly (T1 compose). `None` for flatten `common`.
        primitive_k: Option<i64>,
    },
}

impl ExtensionTower {
    /// **Stable** — `Poly::dimension`
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
    /// **Stable** — `Poly::is_simple_over_q`
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

// **Pipeline private** — `min_poly_key` (R1 semantic_key / adjoin cache)
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

// **Pipeline private** — `parent_blocks_key`
fn parent_blocks_key(blocks: &[CoordsQ]) -> Vec<u8> {
    let mut key = b"parent_blocks:".to_vec();
    for b in blocks {
        key.extend_from_slice(&min_poly_key(b));
    }
    key
}

impl ExtensionTower {
    /// Stable byte key for cache / ordering (R1); does not use `field.id()`.
    /// **Pipeline private** — `semantic_key_bytes`
    fn semantic_key_bytes(&self) -> Vec<u8> {
        match self {
            ExtensionTower::Base => Vec::new(),
            ExtensionTower::Adj {
                parent,
                min_poly_q,
                min_poly_parent_blocks,
                ..
            } => {
                let mut key = parent.semantic_key_bytes();
                if let Some(blocks) = min_poly_parent_blocks {
                    key.extend_from_slice(b"pb:");
                    key.extend_from_slice(&parent_blocks_key(blocks));
                } else {
                    key.extend_from_slice(b"rq:");
                    key.extend_from_slice(&min_poly_key(min_poly_q));
                }
                key
            }
        }
    }
}

/// Descriptor for an extension field K / ℚ(α₁,…) (R3).
///
/// [`ExtensionField`] is a stable type alias for this struct.
#[derive(Clone, Debug)]
pub struct ExtensionDesc {
    id: u64,
    tower: Arc<ExtensionTower>,
    /// Immediate parent when registered via adjoin.
    parent_field: Option<Arc<ExtensionDesc>>,
}

/// Shared handle to an extension field K / ℚ(α₁,…).
pub type ExtensionField = ExtensionDesc;

impl PartialEq for ExtensionDesc {
    // **Stable** — `Poly::eq` (R0: tower only; lazy flatten not compared)
    fn eq(&self, other: &Self) -> bool {
        self.tower == other.tower
    }
}

impl Eq for ExtensionDesc {}

/// Diagnostic: how [`ExtensionField::element_*`] reduces coords on this layer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LayerArithMode {
    /// Single `Adj { parent: ℚ, min_poly_q }` with registered parent.
    PrimitiveOverQ,
    /// `min_poly_parent_blocks` or nested adjoin parent chain.
    Tower,
    /// Flatten compositum: `Adj { parent: ℚ }` but `parent_field: None`.
    FlattenOverQ,
}

macro_rules! parent_coeff_ring {
    ($parent:ident) => {
        ParentCoeffRing {
            zero: $parent.zero_coords(),
            one: $parent.one_coords(),
            add: &|a, b| $parent.element_add(a, b),
            sub: &|a, b| $parent.element_sub(a, b),
            mul: &|a, b| $parent.element_mul(a, b),
            neg: &|a| $parent.element_neg(a),
            inv: &|a| $parent.element_inv(a),
            is_zero: &|a| {
                $parent
                    .element_eq_mod(a, &$parent.zero_coords())
                    .unwrap_or(false)
            },
        }
    };
}

impl ExtensionDesc {
    // **Stable** — `Poly::is_base`
    pub(crate) fn is_base(&self) -> bool {
        matches!(self.tower.as_ref(), ExtensionTower::Base)
    }

    /// ℚ as an extension field (dimension 1).
    /// **Stable** — `Poly::rational`
    pub fn rational() -> Arc<Self> {
        static BASE: OnceLock<Arc<ExtensionField>> = OnceLock::new();
        Arc::clone(BASE.get_or_init(|| {
            Arc::new(Self {
                id: next_field_id(),
                tower: Arc::new(ExtensionTower::Base),
                parent_field: None,
            })
        }))
    }

    /// **Stable** — `id`
    pub fn id(&self) -> u64 {
        self.id
    }

    /// **Stable** — `tower`
    pub fn tower(&self) -> &Arc<ExtensionTower> {
        &self.tower
    }

    /// **Stable** — `dimension`
    pub fn dimension(&self) -> usize {
        self.tower.dimension()
    }

    /// Semantic key for cache / canonical ordering (R1); tower + layer minpoly bytes.
    /// **Stable (bounded)** — `semantic_key`
    pub(crate) fn semantic_key(&self) -> Vec<u8> {
        self.tower.semantic_key_bytes()
    }

    /// Modulus for primitive (ℚ-power-basis) arithmetic: layer minpoly over ℚ.
    /// **Pipeline private** — `primitive_modulus`
    fn primitive_modulus(&self) -> Result<&CoordsQ, EvalError> {
        match self.tower.as_ref() {
            ExtensionTower::Base => {
                // ponytail: MSRV 1.75 — process-wide empty modulus for dim-1 ℚ primitive ops;
                // upstream has no `_EXT` minpoly for ℚ; upgrade path: `LazyLock` or inline
                // empty slice if `poly_reduce` accepts it (R6+ may drop this static).
                static Q_MIN: OnceLock<CoordsQ> = OnceLock::new();
                Ok(Q_MIN.get_or_init(Vec::new))
            }
            ExtensionTower::Adj {
                min_poly_parent_blocks: Some(_),
                ..
            } => Err(EvalError::TypeError("primitive modulus: expected simple over Q")),
            ExtensionTower::Adj {
                min_poly_q,
                parent,
                ..
            } if matches!(parent.as_ref(), ExtensionTower::Base) => Ok(min_poly_q),
            ExtensionTower::Adj { .. } => {
                Err(EvalError::TypeError("primitive modulus: expected simple over Q"))
            }
        }
    }

    /// Immediate parent field when this handle comes from a tower adjoin (T1+).
    /// **Stable** — `parent_field`
    pub fn parent_field(&self) -> Option<&Arc<ExtensionField>> {
        self.parent_field.as_ref()
    }

    /// Adjoin a root of irreducible (or degree ≥ 1) `min_poly_q` over `parent`.
    ///
    /// No session cache: each call constructs a fresh field (R5b). Prefer
    /// [`super::field_session::FieldSession::adjoin_irreducible`] when dedup matters.
    /// **Stable** — `adjoin_irreducible`
    pub fn adjoin_irreducible(
        parent: &Arc<ExtensionField>,
        min_poly_q: CoordsQ,
    ) -> Result<Arc<Self>, EvalError> {
        build_adjoin_irreducible(parent, min_poly_q)
    }

    /// Adjoin a root of `layer_blocks` over `parent` when layer minpoly has **parent** coefficients (T3+).
    ///
    /// `layer_blocks` is high-first monic (e.g. `[1, 0, −α]` for `u²−α`).
    /// **Stable** — `adjoin_irreducible_parent_coeffs`
    pub fn adjoin_irreducible_parent_coeffs(
        parent: &Arc<ExtensionField>,
        layer_blocks: Vec<CoordsQ>,
    ) -> Result<Arc<Self>, EvalError> {
        build_adjoin_parent_coeffs(parent, layer_blocks)
    }

    /// Build ℚ(α) from an irreducible (or at least degree ≥ 1) polynomial over ℚ.
    /// **Stable** — `adjoin_irreducible_over_q`
    pub fn adjoin_irreducible_over_q(min_poly_q: CoordsQ) -> Result<Arc<Self>, EvalError> {
        Self::adjoin_irreducible(&Self::rational(), min_poly_q)
    }

    /// True when `sub` embeds in `sup` along the registered adjoin parent chain (T1).
    /// **Stable** — `is_subfield_of`
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
    /// **Partial** — optional algorithm path `try_subfield_embedding`
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

    /// Layer (adjoin) minimal polynomial as `Expr` coefficients (giac `poly1` order).
    ///
    /// Coefficients lie in the parent field for T3+ towers; over ℚ for a single adjoin over ℚ.
    /// **Stable** — `layer_min_poly_exprs`
    pub fn layer_min_poly_exprs(&self) -> Result<Vec<ExprArc>, EvalError> {
        if self.is_base() {
            return Ok(vec![Expr::int(0)]);
        }
        match self.tower.as_ref() {
            ExtensionTower::Adj {
                min_poly_parent_blocks: Some(blocks),
                ..
            } => {
                let mut out = Vec::with_capacity(blocks.len());
                for block in blocks {
                    out.push(layer_minpoly_block_to_expr(block)?);
                }
                Ok(out)
            }
            ExtensionTower::Adj { min_poly_q, .. } => {
                super::field_arith::coords_to_expr(min_poly_q)
            }
            ExtensionTower::Base => Ok(vec![Expr::int(0)]),
        }
    }

    /// Flattened primitive polynomial over ℚ; roots / legacy metadata.
    /// **Stable** — `top_min_poly_exprs`
    pub fn top_min_poly_exprs(&self) -> Result<Vec<ExprArc>, EvalError> {
        if self.is_base() {
            return Ok(vec![Expr::int(0)]);
        }
        super::field_arith::coords_to_expr(&flatten_min_poly_over_q(self, None)?)
    }

    /// **Stable** — `zero_coords`
    pub fn zero_coords(&self) -> CoordsQ {
        vec![Ratio::zero(); self.dimension()]
    }

    /// **Stable** — `one_coords`
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

    // **Stable** — `generator_coords`
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
        if v.len() >= 2 {
            let gen_idx = v.len() - 2;
            v[gen_idx] = Ratio::one();
        } else if !v.is_empty() {
            v[0] = Ratio::one();
        }
        v
    }

    // **Pipeline private** — `uses_tower_arithmetic`
    fn uses_tower_arithmetic(&self) -> bool {
        if matches!(
            self.tower.as_ref(),
            ExtensionTower::Adj {
                min_poly_parent_blocks: Some(_),
                ..
            }
        ) {
            return true;
        }
        self.parent_field
            .as_ref()
            .is_some_and(|parent| !parent.is_base())
    }

    /// **Pipeline private** — `layer_arith_mode`
    pub(crate) fn layer_arith_mode(&self) -> LayerArithMode {
        if self.is_base() {
            return LayerArithMode::PrimitiveOverQ;
        }
        match self.tower.as_ref() {
            ExtensionTower::Adj {
                min_poly_parent_blocks: Some(_),
                ..
            } => LayerArithMode::Tower,
            ExtensionTower::Adj { parent, .. }
                if matches!(parent.as_ref(), ExtensionTower::Base) =>
            {
                if self.parent_field.is_none() {
                    LayerArithMode::FlattenOverQ
                } else {
                    LayerArithMode::PrimitiveOverQ
                }
            }
            ExtensionTower::Adj { .. } => LayerArithMode::Tower,
            ExtensionTower::Base => LayerArithMode::PrimitiveOverQ,
        }
    }

    /// True when [`Self`] may be parent of parent-coeff adjoin (non-trivial `element_*`).
    /// **Pipeline private** — `parent_coeff_ring_capable`
    pub(crate) fn parent_coeff_ring_capable(&self) -> bool {
        !self.is_base()
    }

    /// Embed a rational constant into this field (constant term at layer 0).
    /// **Stable** — `embed_rational`
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

    // **Pipeline private** — `layer_ext_degree`
    pub(crate) fn layer_ext_degree(&self) -> Result<usize, EvalError> {
        match self.tower.as_ref() {
            ExtensionTower::Adj { ext_degree, .. } => Ok(*ext_degree),
            ExtensionTower::Base => Err(EvalError::TypeError("expected extension field")),
        }
    }

    // **Pipeline private** — `layer_minpoly_parent_coeffs`
    fn layer_minpoly_parent_coeffs(&self, parent: &ExtensionField) -> Result<Vec<CoordsQ>, EvalError> {
        match self.tower.as_ref() {
            ExtensionTower::Adj {
                min_poly_parent_blocks: Some(blocks),
                ..
            } => Ok(blocks.clone()),
            ExtensionTower::Adj { min_poly_q, .. } => Ok(min_poly_q
                .iter()
                .map(|r| parent.embed_rational(r))
                .collect()),
            ExtensionTower::Base => Err(EvalError::TypeError("expected extension field")),
        }
    }

    // **Pipeline private** — `unflatten_layer_blocks`
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

    // **Pipeline private** — `flatten_layer_element`
    fn flatten_layer_element(&self, blocks: &[CoordsQ]) -> Result<CoordsQ, EvalError> {
        let parent = self
            .parent_field
            .as_ref()
            .ok_or(EvalError::TypeError("tower field without parent"))?;
        let pd = parent.dimension();
        let e = self.layer_ext_degree()?;
        let mut blocks = blocks.to_vec();
        blocks.reverse();
        // Missing high-u blocks were trimmed; append zero blocks (u^1, u^2, …).
        while blocks.len() < e {
            blocks.push(parent.zero_coords());
        }
        Ok(flatten_layer_blocks(&blocks, pd))
    }

    /// **Stable** — `element_add`
    pub fn element_add(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_add_tower(a, b)
        } else {
            self.element_add_primitive(a, b)
        }
    }

    // **Pipeline private** — `element_add_tower`
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

    // **Pipeline private** — `element_add_primitive`
    fn element_add_primitive(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        self.ensure_same_field_len(b)?;
        if self.is_base() {
            return Ok(vec![pad_to_len(a, 1)[0].clone() + pad_to_len(b, 1)[0].clone()]);
        }
        Ok(pad_to_len(
            &poly_reduce(&poly_add(a, b), self.primitive_modulus()?),
            self.dimension(),
        ))
    }

    /// **Stable** — `element_sub`
    pub fn element_sub(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_sub_tower(a, b)
        } else {
            self.element_sub_primitive(a, b)
        }
    }

    // **Pipeline private** — `element_sub_tower`
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

    // **Pipeline private** — `element_sub_primitive`
    fn element_sub_primitive(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        self.ensure_same_field_len(b)?;
        if self.is_base() {
            return Ok(vec![pad_to_len(a, 1)[0].clone() - pad_to_len(b, 1)[0].clone()]);
        }
        Ok(pad_to_len(
            &poly_reduce(&poly_sub(a, b), self.primitive_modulus()?),
            self.dimension(),
        ))
    }

    /// **Stable** — `element_neg`
    pub fn element_neg(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_neg_tower(a)
        } else {
            self.element_neg_primitive(a)
        }
    }

    // **Pipeline private** — `element_neg_tower`
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

    // **Pipeline private** — `element_neg_primitive`
    fn element_neg_primitive(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        if self.is_base() {
            return Ok(vec![-pad_to_len(a, 1)[0].clone()]);
        }
        Ok(pad_to_len(
            &poly_reduce(&poly_neg(a), self.primitive_modulus()?),
            self.dimension(),
        ))
    }

    /// **Stable** — `element_mul`
    pub fn element_mul(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_mul_tower(a, b)
        } else {
            self.element_mul_primitive(a, b)
        }
    }

    // **Pipeline private** — `element_mul_tower`
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

    // **Pipeline private** — `element_mul_primitive`
    fn element_mul_primitive(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
        self.ensure_same_field_len(a)?;
        self.ensure_same_field_len(b)?;
        if self.is_base() {
            return Ok(vec![pad_to_len(a, 1)[0].clone() * pad_to_len(b, 1)[0].clone()]);
        }
        Ok(pad_to_len(
            &poly_reduce(&poly_mul(a, b), self.primitive_modulus()?),
            self.dimension(),
        ))
    }

    /// **Stable** — `element_inv`
    pub fn element_inv(&self, a: &CoordsQ) -> Result<CoordsQ, EvalError> {
        if self.uses_tower_arithmetic() {
            self.element_inv_tower(a)
        } else {
            self.element_inv_primitive(a)
        }
    }

    // **Pipeline private** — `element_inv_tower`
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

    // **Pipeline private** — `element_inv_primitive`
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
        let inv = poly_inv_mod(a, self.primitive_modulus()?)?;
        Ok(pad_to_len(&inv, self.dimension()))
    }

    /// **Stable** — `element_eq_mod`
    pub fn element_eq_mod(&self, a: &CoordsQ, b: &CoordsQ) -> Result<bool, EvalError> {
        let diff = self.element_sub(a, b)?;
        Ok(coords_all_zero(&diff))
    }

    /// Try to find ε ∈ K with ε² ≡ u (coords). Returns `None` if no probe hits.
    ///
    /// Scan: deg-2 layer generators (±g, k·g), basis ±eᵢ, then bounded ±1 linear
    /// combos when `dim(K) ≤` pairwise/triple ceilings (see `TRY_SQRT_*_DIM_CEILING`).
    ///
    /// Does **not** adjoin. For adjoin fallback use
    /// [`super::alg_ext::algext_square_roots`] or
    /// [`super::field_session::FieldSession::sqrt_in_field`].
    /// **Stable (bounded)** — try square root in extension field
    pub fn try_square_root_in_field(
        field: &Arc<ExtensionField>,
        u: &CoordsQ,
    ) -> Result<Option<CoordsQ>, EvalError> {
        try_square_root_in_field_impl(field, u)
    }

    /// Layer generators + basis squares only (no ±1 combo enum). Hot path for Euler
    /// second sqrt when blind adjoin is expected.
    /// **Pipeline private** — shallow sqrt probe
    pub(crate) fn try_square_root_in_field_shallow(
        field: &Arc<ExtensionField>,
        u: &CoordsQ,
    ) -> Result<Option<CoordsQ>, EvalError> {
        try_square_root_in_field_shallow_impl(field, u)
    }

    /// **Stable** — `element_is_zero`
    pub fn element_is_zero(&self, a: &CoordsQ) -> bool {
        a.iter().all(|c| c.is_zero())
    }

    /// **Stable** — `element_is_one`
    pub fn element_is_one(&self, a: &CoordsQ) -> bool {
        pad_to_len(a, self.dimension()) == self.one_coords()
    }

    // **Pipeline private** — `ensure_same_field_len`
    fn ensure_same_field_len(&self, a: &CoordsQ) -> Result<(), EvalError> {
        if a.len() > self.dimension() {
            return Err(EvalError::TypeError("coords longer than field dimension"));
        }
        Ok(())
    }

    /// Minimal common extension of two ℚ-described fields, with linear embeddings.
    ///
    /// Uses an ephemeral per-call cache (R5b: no cross-call dedup). Prefer
    /// [`super::field_session::FieldSession::common_over_q`] when a session is available.
    /// **Stable** — `common_over_q`
    pub fn common_over_q(
        a: &Arc<ExtensionField>,
        b: &Arc<ExtensionField>,
    ) -> Result<Arc<CommonFieldPair>, EvalError> {
        with_ephemeral_common_cache(|cache| common_over_q_in_cache(a, b, cache, None))
    }

    /// Like [`Self::common_over_q`] with an explicit session cache (R2).
    /// **Stable (bounded)** — common with session cache
    pub fn common_over_q_with_cache(
        a: &Arc<ExtensionField>,
        b: &Arc<ExtensionField>,
        cache: &mut CommonCache,
    ) -> Result<Arc<CommonFieldPair>, EvalError> {
        common_over_q_in_cache(a, b, cache, None)
    }

    /// Pick the embedding whose [`FieldEmbedding::source`] matches `source`.
    ///
    /// `CommonFieldPair::embed_a` / `embed_b` are keyed by sorted semantic keys in the cache,
    /// not by caller operand order — always use this when applying a cached pair.
    /// **Stable** — `embedding_for`
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
    /// Uses an ephemeral per-call cache (R5b). Prefer
    /// [`super::field_session::FieldSession::align_elements`] with a [`Context`](crate::context::Context) session.
    /// **Stable** — `align_elements`
    pub fn align_elements(
        a_field: &Arc<ExtensionField>,
        a_coords: &CoordsQ,
        b_field: &Arc<ExtensionField>,
        b_coords: &CoordsQ,
    ) -> Result<AlignedElements, EvalError> {
        with_ephemeral_common_cache(|cache| {
            align_elements_in_cache(a_field, a_coords, b_field, b_coords, cache, None)
        })
    }

    /// Like [`Self::align_elements`] with an explicit session cache (R2).
    /// **Stable (bounded)** — align with session cache
    pub fn align_elements_with_cache(
        a_field: &Arc<ExtensionField>,
        a_coords: &CoordsQ,
        b_field: &Arc<ExtensionField>,
        b_coords: &CoordsQ,
        cache: &mut CommonCache,
    ) -> Result<AlignedElements, EvalError> {
        align_elements_in_cache(a_field, a_coords, b_field, b_coords, cache, None)
    }
}

/// Session-local `common_over_q` cache (R2).
pub(crate) type CommonCache = HashMap<(Vec<u8>, Vec<u8>), Arc<CommonFieldPair>>;

// **Pipeline private** — ephemeral cache for static API (R5b: no cross-call dedup).
fn with_ephemeral_common_cache<R>(f: impl FnOnce(&mut CommonCache) -> R) -> R {
    f(&mut HashMap::new())
}

// **Pipeline private** — `common_over_q_in_cache`
pub(crate) fn common_over_q_in_cache(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    cache: &mut CommonCache,
    session: Option<&super::field_session::FieldSession>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    if Arc::ptr_eq(a, b) || **a == **b {
        return Ok(CommonFieldPair::identity(a));
    }
    let key = common_cache_key(a, b);
    if let Some(hit) = cache.get(&key) {
        return Ok(Arc::clone(hit));
    }
    let ka = a.semantic_key();
    let kb = b.semantic_key();
    let (fa, fb) = if ka <= kb { (a, b) } else { (b, a) };
    let pair = compute_common_dispatch(fa, fb, session)?;
    cache.insert(key, Arc::clone(&pair));
    Ok(pair)
}

// **Pipeline private** — `align_elements_in_cache`
pub(crate) fn align_elements_in_cache(
    a_field: &Arc<ExtensionField>,
    a_coords: &CoordsQ,
    b_field: &Arc<ExtensionField>,
    b_coords: &CoordsQ,
    cache: &mut CommonCache,
    session: Option<&super::field_session::FieldSession>,
) -> Result<AlignedElements, EvalError> {
    if Arc::ptr_eq(a_field, b_field) || **a_field == **b_field {
        return Ok(AlignedElements {
            field: Arc::clone(a_field),
            left: pad_to_len(a_coords, a_field.dimension()),
            right: pad_to_len(b_coords, b_field.dimension()),
        });
    }
    if ExtensionField::is_subfield_of(a_field, b_field) {
        let emb_a = ExtensionField::try_subfield_embedding(a_field, b_field)?
            .ok_or(EvalError::TypeError("subfield embedding missing"))?;
        return Ok(AlignedElements {
            field: Arc::clone(b_field),
            left: emb_a.apply(a_coords),
            right: pad_to_len(b_coords, b_field.dimension()),
        });
    }
    if ExtensionField::is_subfield_of(b_field, a_field) {
        let emb_b = ExtensionField::try_subfield_embedding(b_field, a_field)?
            .ok_or(EvalError::TypeError("subfield embedding missing"))?;
        return Ok(AlignedElements {
            field: Arc::clone(a_field),
            left: pad_to_len(a_coords, a_field.dimension()),
            right: emb_b.apply(b_coords),
        });
    }
    let common = common_over_q_in_cache(a_field, b_field, cache, session)?;
    let left = ExtensionField::embedding_for(a_field, &common)?.apply(a_coords);
    let right = ExtensionField::embedding_for(b_field, &common)?.apply(b_coords);
    Ok(AlignedElements {
        field: Arc::clone(&common.field),
        left,
        right,
    })
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
    /// **Stable** — `Poly::apply`
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
///
/// Each [`FieldEmbedding`] maps **`source` → `target`** (the common ambient field is
/// [`Self::field`]). Cache entries sort operands by `field.id`; use
/// [`ExtensionField::embedding_for`] to pick the matrix for a given operand, not
/// `embed_a` / `embed_b` position.
#[derive(Clone, Debug)]
pub struct CommonFieldPair {
    pub field: Arc<ExtensionField>,
    pub embed_a: FieldEmbedding,
    pub embed_b: FieldEmbedding,
}

impl CommonFieldPair {
    /// **Stable** — `Poly::identity`
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

/// U0 compositum gate — see `.doc/issues/GIAC-ext-common-unified-path.md` §3.
// **Pipeline private** — `verify_common_pair`
pub(crate) fn verify_common_pair(
    pair: &CommonFieldPair,
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
) -> Result<(), EvalError> {
    let subfield =
        ExtensionField::is_subfield_of(a, b) || ExtensionField::is_subfield_of(b, a);
    if !subfield {
        let expected = a.dimension().checked_mul(b.dimension()).ok_or(
            EvalError::TypeError("common verify: dimension overflow"),
        )?;
        if pair.field.dimension() != expected {
            return Err(EvalError::TypeError("common verify: compositum degree mismatch"));
        }
    }
    verify_embedded_one(pair, a)?;
    verify_embedded_one(pair, b)?;
    verify_embedded_generator(pair, a)?;
    verify_embedded_generator(pair, b)?;
    Ok(())
}

// **Pipeline private** — `element_pow_coords`
fn element_pow_coords(
    field: &ExtensionField,
    base: &CoordsQ,
    exp: u64,
) -> Result<CoordsQ, EvalError> {
    let mut out = field.one_coords();
    if exp == 0 {
        return Ok(out);
    }
    let mut pow = base.clone();
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            out = field.element_mul(&out, &pow)?;
        }
        e >>= 1;
        if e > 0 {
            pow = field.element_mul(&pow, &pow)?;
        }
    }
    Ok(out)
}

// **Pipeline private** — Horner eval of monic poly1 `[c_deg,…,c_0]` at field element.
fn eval_rational_poly1_at_field(
    field: &ExtensionField,
    poly: &[Ratio<BigInt>],
    x: &CoordsQ,
) -> Result<CoordsQ, EvalError> {
    let deg = poly_degree(poly);
    let mut acc = field.zero_coords();
    for (i, c) in poly.iter().enumerate().take(deg + 1) {
        if c.is_zero() {
            continue;
        }
        let exp = (deg - i) as u64;
        let term = field.element_mul(&field.embed_rational(c), &element_pow_coords(field, x, exp)?)?;
        acc = field.element_add(&acc, &term)?;
    }
    Ok(acc)
}

// **Pipeline private** — `verify_embedded_one`
fn verify_embedded_one(
    pair: &CommonFieldPair,
    operand: &Arc<ExtensionField>,
) -> Result<(), EvalError> {
    if operand.is_base() {
        return Ok(());
    }
    let emb = ExtensionField::embedding_for(operand, pair)?;
    let one = emb.apply(&operand.one_coords());
    if !pair
        .field
        .element_eq_mod(&one, &pair.field.one_coords())
        .unwrap_or(false)
    {
        return Err(EvalError::TypeError("common verify: embed(1) != 1"));
    }
    Ok(())
}

// **Pipeline private** — layer generator + minpoly vanishing in common field.
fn verify_embedded_generator(
    pair: &CommonFieldPair,
    operand: &Arc<ExtensionField>,
) -> Result<(), EvalError> {
    if operand.is_base() {
        return Ok(());
    }
    let emb = ExtensionField::embedding_for(operand, pair)?;
    let gen = operand.generator_coords();
    let gen_sq = operand.element_mul(&gen, &gen)?;
    let img = emb.apply(&gen);
    let img_sq = emb.apply(&gen_sq);
    let img_g_sq = pair.field.element_mul(&img, &img)?;
    if !pair
        .field
        .element_eq_mod(&img_g_sq, &img_sq)
        .unwrap_or(false)
    {
        return Err(EvalError::TypeError(
            "common verify: embed(g)^2 != embed(g^2)",
        ));
    }
    match layer_minpoly_coords_for_adjoin(operand)? {
        LayerMinPolyForAdjoin::Rational(m) => {
            let val = eval_rational_poly1_at_field(&pair.field, m, &img)?;
            if !pair.field.element_is_zero(&val) {
                return Err(EvalError::TypeError(
                    "common verify: layer minpoly does not vanish at embed(g)",
                ));
            }
        }
        LayerMinPolyForAdjoin::ParentBlocks(blocks) => {
            let parent = operand
                .parent_field
                .as_ref()
                .ok_or(EvalError::TypeError("common verify: no parent"))?;
            let parent_emb = ExtensionField::embedding_for(parent, pair)?;
            let pd = parent.dimension();
            let deg = blocks.len().checked_sub(1).ok_or(EvalError::TypeError(
                "common verify: layer minpoly",
            ))?;
            let mut acc = pair.field.zero_coords();
            for (j, block) in blocks.iter().enumerate() {
                let c = parent_emb.apply(&pad_to_len(block, pd));
                let exp = (deg - j) as u64;
                let term = pair.field.element_mul(
                    &c,
                    &element_pow_coords(&pair.field, &img, exp)?,
                )?;
                acc = pair.field.element_add(&acc, &term)?;
            }
            if !pair.field.element_is_zero(&acc) {
                return Err(EvalError::TypeError(
                    "common verify: parent-blocks minpoly does not vanish at embed(g)",
                ));
            }
        }
    }
    Ok(())
}

/// Two extension elements aligned into a common ambient field.
#[derive(Clone, Debug)]
pub struct AlignedElements {
    pub field: Arc<ExtensionField>,
    pub left: CoordsQ,
    pub right: CoordsQ,
}

// **Pipeline private** — `common_cache_key`
fn common_cache_key(a: &ExtensionField, b: &ExtensionField) -> (Vec<u8>, Vec<u8>) {
    let ka = a.semantic_key();
    let kb = b.semantic_key();
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

/// **Pipeline private** — byte key for ℚ minpoly dedup (R4 session cache).
pub(crate) fn min_poly_key_bytes(p: &[Ratio<BigInt>]) -> Vec<u8> {
    min_poly_key(p)
}

/// **Pipeline private** — adjoin dedup key (R4).
pub(crate) fn adjoin_cache_key_rational(parent: &ExtensionField, min_poly_q: &CoordsQ) -> Vec<u8> {
    let mut k = parent.semantic_key();
    k.extend_from_slice(b"rq:");
    k.extend_from_slice(&min_poly_key(min_poly_q));
    k
}

/// **Pipeline private** — parent-coeff adjoin dedup key (R4).
pub(crate) fn adjoin_cache_key_parent_blocks(
    parent: &ExtensionField,
    blocks: &[CoordsQ],
) -> Vec<u8> {
    let mut k = parent.semantic_key();
    k.extend_from_slice(&parent_blocks_key(blocks));
    k
}

/// ponytail: S6 pairwise ±1 fallback only (F5); no dim³ enum.
const TRY_SQRT_PAIRWISE_DIM_CEILING: usize = 6;

// **Pipeline private** — collect proper subfields along parent chain (immediate → … → ℚ).
fn proper_subfields_chain(field: &Arc<ExtensionField>) -> Vec<Arc<ExtensionField>> {
    let mut out = Vec::new();
    let mut cur = field.parent_field().map(Arc::clone);
    while let Some(p) = cur {
        out.push(Arc::clone(&p));
        cur = p.parent_field().map(Arc::clone);
    }
    out
}

// **Pipeline private** — solve M·v = u for embedding matrix (tgt × src); None if u ∉ Im(M).
pub(crate) fn try_preimage_under_embedding(
    emb: &FieldEmbedding,
    u: &CoordsQ,
) -> Option<CoordsQ> {
    let rows = emb.matrix.len();
    if rows == 0 {
        return None;
    }
    let cols = emb.matrix[0].len();
    let u = pad_to_len(u, rows);
    let mut aug = vec![vec![Ratio::zero(); cols + 1]; rows];
    for i in 0..rows {
        for j in 0..cols {
            aug[i][j] = emb.matrix[i][j].clone();
        }
        aug[i][cols] = u[i].clone();
    }
    let (pivot_cols, inconsistent) = gauss_elim_rref(&mut aug);
    if inconsistent {
        return None;
    }
    let mut v = vec![Ratio::zero(); cols];
    let mut used = vec![false; cols];
    for (row, &pc) in pivot_cols.iter().enumerate() {
        if pc >= cols {
            continue;
        }
        used[pc] = true;
        v[pc] = aug[row][cols].clone();
    }
    Some(v)
}

// **Pipeline private** — ℚ Gaussian elimination; returns pivot column per row, inconsistent flag.
pub(crate) fn gauss_elim_rref(aug: &mut [Vec<Ratio<BigInt>>]) -> (Vec<usize>, bool) {
    let rows = aug.len();
    let cols = aug[0].len().saturating_sub(1);
    let mut pivot_cols = vec![cols; rows];
    let mut pivot_row = 0;
    for col in 0..cols {
        if pivot_row >= rows {
            break;
        }
        let mut sel = None;
        for r in pivot_row..rows {
            if !aug[r][col].is_zero() {
                sel = Some(r);
                break;
            }
        }
        let Some(r) = sel else { continue };
        aug.swap(pivot_row, r);
        let pivot = aug[pivot_row][col].clone();
        for c in col..=cols {
            aug[pivot_row][c] /= pivot.clone();
        }
        for r in 0..rows {
            if r == pivot_row || aug[r][col].is_zero() {
                continue;
            }
            let factor = aug[r][col].clone();
            let pivot_row_vals: Vec<Ratio<BigInt>> = aug[pivot_row][col..=cols].to_vec();
            for (c_offset, pv) in pivot_row_vals.iter().enumerate() {
                aug[r][col + c_offset] -= factor.clone() * pv.clone();
            }
        }
        pivot_cols[pivot_row] = col;
        pivot_row += 1;
    }
    let inconsistent = (pivot_row..rows).any(|r| {
        !aug[r][cols].is_zero() && (0..cols).all(|c| aug[r][c].is_zero())
    });
    (pivot_cols, inconsistent)
}

// **Pipeline private** — F5 S3: u ∈ F ⊂ L → sqrt in F, embed back.
fn try_sqrt_subfield_descent(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    for sub in proper_subfields_chain(field) {
        let emb = match ExtensionField::try_subfield_embedding(&sub, field)? {
            Some(e) => e,
            None => continue,
        };
        let Some(u_sub) = try_preimage_under_embedding(&emb, u) else {
            continue;
        };
        if let Some(root_sub) = try_square_root_in_field_impl(&sub, &u_sub)? {
            return Ok(Some(emb.apply(&root_sub)));
        }
    }
    Ok(None)
}

// **Pipeline private** — F5 S4: top deg-2 layer u = a + b·g, (a+bg)² = u (standard 2×pd layout).
fn try_sqrt_quadratic_top_layer(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let parent = match field.parent_field() {
        Some(p) => Arc::clone(p),
        None => return Ok(None),
    };
    if field.layer_ext_degree().ok() != Some(2) {
        return Ok(None);
    }
    let pd = parent.dimension();
    let cd = field.dimension();
    if cd != 2 * pd {
        return Ok(None);
    }
    let u = pad_to_len(u, cd);
    let u_a = pad_to_len(&u[0..pd], pd);
    let u_b = pad_to_len(&u[pd..2 * pd], pd);
    let g = field.generator_coords();
    let g_sq = field.element_mul(&g, &g)?;
    let d = pad_to_len(&g_sq[0..pd], pd);

    let try_pair = |a: &CoordsQ, b: &CoordsQ| -> Result<Option<CoordsQ>, EvalError> {
        let mut out = parent.zero_coords();
        out.extend(pad_to_len(a, pd));
        out.extend(pad_to_len(b, pd));
        let out = pad_to_len(&out, cd);
        for cand in [out.clone(), field.element_neg(&out)?] {
            if coords_square_eq_mod(field, &cand, &u)? {
                return Ok(Some(cand));
            }
        }
        Ok(None)
    };

    // u = b²·d (no g-component): √u = ±b·g
    if u_b.iter().all(|c| c.is_zero()) {
        if let Ok(inv_d) = parent.element_inv(&d) {
            let ratio = parent.element_mul(&u_a, &inv_d)?;
            if let Some(b_cand) = try_square_root_in_field_impl(&parent, &ratio)? {
                if let Some(hit) = try_pair(&parent.zero_coords(), &b_cand)? {
                    return Ok(Some(hit));
                }
            }
        }
    } else if !u_a.iter().all(|c| c.is_zero()) {
        // F5 S4 general: (a+bg)² = u_a + u_b·g, g²=d → a² = (u_a ± √(u_a² − u_b²d)) / 2.
        let u_b_sq = parent.element_mul(&u_b, &u_b)?;
        let u_b_sq_d = parent.element_mul(&u_b_sq, &d)?;
        let ua_sq = parent.element_mul(&u_a, &u_a)?;
        let inner = parent.element_sub(&ua_sq, &u_b_sq_d)?;
        let two = parent.embed_rational(&Ratio::from_integer(2.into()));
        let Ok(two_inv) = parent.element_inv(&two) else {
            return Ok(None);
        };
        let Some(sqrt_inner) = try_square_root_in_field_impl(&parent, &inner)? else {
            return Ok(None);
        };
        for sign in [1i32, -1i32] {
            let signed = if sign > 0 {
                sqrt_inner.clone()
            } else {
                parent.element_neg(&sqrt_inner)?
            };
            let a_sq = parent.element_mul(&parent.element_add(&u_a, &signed)?, &two_inv)?;
            if let Some(a_cand) = try_square_root_in_field_impl(&parent, &a_sq)? {
                let two_a = parent.element_mul(&two, &a_cand)?;
                if let Ok(inv_two_a) = parent.element_inv(&two_a) {
                    let b_cand = parent.element_mul(&u_b, &inv_two_a)?;
                    for a_try in [a_cand.clone(), parent.element_neg(&a_cand)?] {
                        if let Some(hit) = try_pair(&a_try, &b_cand)? {
                            return Ok(Some(hit));
                        }
                    }
                }
            }
        }
    }
    Ok(None)
}

// **Pipeline private** — collect deg-2 layer generators embedded in `field`.
fn quadratic_layer_generators(field: &Arc<ExtensionField>) -> Result<Vec<CoordsQ>, EvalError> {
    let mut gens = Vec::new();
    let mut cur: Option<Arc<ExtensionField>> = Some(Arc::clone(field));
    while let Some(layer) = cur {
        if layer.layer_ext_degree().ok() == Some(2) {
            gens.push(embed_coords_in(field, &layer, &layer.generator_coords())?);
        }
        cur = layer.parent_field().map(Arc::clone);
    }
    Ok(gens)
}

// **Pipeline private** — F5 S5: products / ratios of quadratic layer generators.
fn try_sqrt_layer_generator_algebra(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let gens = quadratic_layer_generators(field)?;
    let n = gens.len();
    for i in 0..n {
        for j in i..n {
            let prod = if i == j {
                gens[i].clone()
            } else {
                field.element_mul(&gens[i], &gens[j])?
            };
            for cand in [prod.clone(), field.element_neg(&prod)?] {
                if coords_square_eq_mod(field, &cand, u)? {
                    return Ok(Some(cand));
                }
            }
            if i != j {
                if let Ok(inv_j) = field.element_inv(&gens[j]) {
                    let quot = field.element_mul(&gens[i], &inv_j)?;
                    for cand in [quot.clone(), field.element_neg(&quot)?] {
                        if coords_square_eq_mod(field, &cand, u)? {
                            return Ok(Some(cand));
                        }
                    }
                }
            }
        }
    }
    try_sqrt_small_combo_with_quadratic_gens(field, u)
}

// **Pipeline private** — S6: dim≤6 pairwise ±1 ponytail fallback.
fn try_sqrt_pairwise_fallback(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let dim = field.dimension();
    if dim > TRY_SQRT_PAIRWISE_DIM_CEILING {
        return Ok(None);
    }
    for i in 0..dim {
        for j in (i + 1)..dim {
            for &ci in &[1i32, -1i32] {
                for &cj in &[1i32, -1i32] {
                    let mut v = field.zero_coords();
                    v[i] = Ratio::from_integer(BigInt::from(ci));
                    v[j] = Ratio::from_integer(BigInt::from(cj));
                    if field.element_is_zero(&v) {
                        continue;
                    }
                    for cand in [v.clone(), field.element_neg(&v)?] {
                        if coords_square_eq_mod(field, &cand, u)? {
                            return Ok(Some(cand));
                        }
                    }
                }
            }
        }
    }
    Ok(None)
}

/// ponytail: k·g scan ceiling when dim>8 (F5 S1).
fn try_sqrt_layer_generators(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let k_max = if field.dimension() > 8 { 4 } else { 8 };
    let mut cur: Option<Arc<ExtensionField>> = Some(Arc::clone(field));
    while let Some(layer) = cur {
        if layer.layer_ext_degree().ok() == Some(2) {
            let gen = embed_coords_in(field, &layer, &layer.generator_coords())?;
            for cand in [gen.clone(), field.element_neg(&gen)?] {
                if coords_square_eq_mod(field, &cand, u)? {
                    return Ok(Some(cand));
                }
            }
            for k in 2i32..=k_max {
                let k_rat = field.embed_rational(&Ratio::from_integer(BigInt::from(k)));
                for cand in [
                    field.element_mul(&k_rat, &gen)?,
                    field.element_neg(&field.element_mul(&k_rat, &gen)?)?,
                ] {
                    if coords_square_eq_mod(field, &cand, u)? {
                        return Ok(Some(cand));
                    }
                }
            }
        }
        cur = layer.parent_field().map(Arc::clone);
    }
    Ok(None)
}

// **Pipeline private** — operational basis ±eᵢ.
fn try_sqrt_basis_squares(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let dim = field.dimension();
    for i in 0..dim {
        let mut e = field.zero_coords();
        e[i] = Ratio::one();
        for cand in [e.clone(), field.element_neg(&e)?] {
            if coords_square_eq_mod(field, &cand, u)? {
                return Ok(Some(cand));
            }
        }
    }
    Ok(None)
}

// **Pipeline private** — F5 core: structural probes S0–S6 (no dim³ enum).
fn try_square_root_in_field_impl(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let u = pad_to_len(u, field.dimension());
    if field.element_is_zero(&u) {
        return Ok(Some(field.zero_coords()));
    }
    if let Some(c) = try_sqrt_layer_generators(field, &u)? {
        return Ok(Some(c));
    }
    if let Some(c) = try_sqrt_basis_squares(field, &u)? {
        return Ok(Some(c));
    }
    if let Some(c) = try_sqrt_subfield_descent(field, &u)? {
        return Ok(Some(c));
    }
    if let Some(c) = try_sqrt_quadratic_top_layer(field, &u)? {
        return Ok(Some(c));
    }
    if let Some(c) = try_sqrt_layer_generator_algebra(field, &u)? {
        return Ok(Some(c));
    }
    try_sqrt_pairwise_fallback(field, &u)
}

// **Pipeline private** — shallow: S1–S3 + S5 products (Euler hot path).
fn try_square_root_in_field_shallow_impl(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let u = pad_to_len(u, field.dimension());
    if field.element_is_zero(&u) {
        return Ok(Some(field.zero_coords()));
    }
    if let Some(c) = try_sqrt_layer_generators(field, &u)? {
        return Ok(Some(c));
    }
    if let Some(c) = try_sqrt_basis_squares(field, &u)? {
        return Ok(Some(c));
    }
    if let Some(c) = try_sqrt_subfield_descent(field, &u)? {
        return Ok(Some(c));
    }
    try_sqrt_layer_generator_algebra(field, &u)
}

// **Pipeline private** — F4: e_i ± k·g combinations for last quadratic adjoin layers (dim≤12).
fn try_sqrt_small_combo_with_quadratic_gens(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let dim = field.dimension();
    let mut gens = Vec::new();
    let mut cur: Option<Arc<ExtensionField>> = Some(Arc::clone(field));
    while let Some(layer) = cur {
        if layer.layer_ext_degree().ok() == Some(2) {
            gens.push(embed_coords_in(field, &layer, &layer.generator_coords())?);
        }
        cur = layer.parent_field().map(Arc::clone);
    }
    for i in 0..dim {
        let mut e = field.zero_coords();
        e[i] = Ratio::one();
        for g in &gens {
            for &ci in &[-2i32, -1, 1, 2] {
                for &cj in &[-2i32, -1, 1, 2] {
                    if ci == 0 && cj == 0 {
                        continue;
                    }
                    let ci_r = field.embed_rational(&Ratio::from_integer(BigInt::from(ci)));
                    let cj_r = field.embed_rational(&Ratio::from_integer(BigInt::from(cj)));
                    let mut v = field.zero_coords();
                    if ci != 0 {
                        v = field.element_add(&v, &field.element_mul(&ci_r, &e)?)?;
                    }
                    if cj != 0 {
                        v = field.element_add(&v, &field.element_mul(&cj_r, g)?)?;
                    }
                    if field.element_is_zero(&v) {
                        continue;
                    }
                    for cand in [v.clone(), field.element_neg(&v)?] {
                        if coords_square_eq_mod(field, &cand, u)? {
                            return Ok(Some(cand));
                        }
                    }
                }
            }
            for cand in [
                field.element_mul(&e, g)?,
                field.element_neg(&field.element_mul(&e, g)?)?,
            ] {
                if coords_square_eq_mod(field, &cand, u)? {
                    return Ok(Some(cand));
                }
            }
        }
    }
    Ok(None)
}

fn embed_coords_in(
    sup: &Arc<ExtensionField>,
    sub: &Arc<ExtensionField>,
    coords: &CoordsQ,
) -> Result<CoordsQ, EvalError> {
    let emb = ExtensionField::try_subfield_embedding(sub, sup)?
        .ok_or(EvalError::TypeError("subfield embedding"))?;
    Ok(emb.apply(coords))
}

fn coords_square_eq_mod(
    field: &ExtensionField,
    cand: &CoordsQ,
    u: &CoordsQ,
) -> Result<bool, EvalError> {
    let sq = field.element_mul(cand, cand)?;
    field.element_eq_mod(&sq, u)
}

/// **Pipeline private** — T1 layer minpoly coefficient check.
// **Pipeline private** — `AlgExtData` from operational coords in `field`.
fn algext_from_coords(field: &Arc<ExtensionField>, coords: &CoordsQ) -> Result<AlgExtData, EvalError> {
    AlgExtData::from_field_coords(
        Arc::clone(field),
        coords_to_expr(&pad_to_len(coords, field.dimension()))?,
    )
}

pub(crate) fn layer_minpoly_rational_constants(min_poly_q: &[Ratio<BigInt>]) -> bool {
    min_poly_q.iter().all(|c| c.denom().is_one())
}

/// **Pipeline private** — construct adjoin field without session dedup (R4).
pub(crate) fn build_adjoin_irreducible(
    parent: &Arc<ExtensionField>,
    min_poly_q: CoordsQ,
) -> Result<Arc<ExtensionField>, EvalError> {
    let ext_degree = poly_degree(&min_poly_q);
    Ok(Arc::new(ExtensionDesc {
        id: next_field_id(),
        tower: Arc::new(ExtensionTower::Adj {
            parent: Arc::clone(parent.tower()),
            min_poly_q,
            min_poly_parent_blocks: None,
            ext_degree,
            primitive_k: None,
        }),
        parent_field: Some(Arc::clone(parent)),
    }))
}

/// **Pipeline private** — construct parent-coeff adjoin without session dedup (R4).
pub(crate) fn build_adjoin_parent_coeffs(
    parent: &Arc<ExtensionField>,
    layer_blocks: Vec<CoordsQ>,
) -> Result<Arc<ExtensionField>, EvalError> {
    if !parent.parent_coeff_ring_capable() {
        return Err(EvalError::TypeError(
            "parent not suitable for parent-coeff adjoin",
        ));
    }
    let ext_degree = layer_blocks.len() - 1;
    Ok(Arc::new(ExtensionDesc {
        id: next_field_id(),
        tower: Arc::new(ExtensionTower::Adj {
            parent: Arc::clone(parent.tower()),
            min_poly_q: Vec::new(),
            min_poly_parent_blocks: Some(layer_blocks),
            ext_degree,
            primitive_k: None,
        }),
        parent_field: Some(Arc::clone(parent)),
    }))
}

/// **Pipeline private** — construct base extension without session dedup (R4/R5b).
pub(crate) fn get_or_create_base_by_min_poly(min_poly_q: CoordsQ) -> Arc<ExtensionField> {
    build_base_extension_uncached(min_poly_q)
}

/// **Pipeline private** — construct base extension without dedup (R4).
pub(crate) fn build_base_extension_uncached(min_poly_q: CoordsQ) -> Arc<ExtensionField> {
    let d = poly_degree(&min_poly_q);
    Arc::new(ExtensionDesc {
        id: next_field_id(),
        tower: Arc::new(ExtensionTower::Adj {
            parent: Arc::new(ExtensionTower::Base),
            min_poly_q: min_poly_q.clone(),
            min_poly_parent_blocks: None,
            ext_degree: d,
            primitive_k: None,
        }),
        parent_field: None,
    })
}

#[cfg(test)]
// **Stable** — `compute_common_flatten_for_test`
pub(crate) fn compute_common_flatten_for_test(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    compute_common_flatten(a, b, None)
}

#[cfg(test)]
/// Second `Arc` handle with the same tower / minpoly but a distinct id (tests only).
// **Stable** — `duplicate_field_arc_for_test`
pub(crate) fn duplicate_field_arc_for_test(f: &Arc<ExtensionField>) -> Arc<ExtensionField> {
    Arc::new(ExtensionDesc {
        id: next_field_id(),
        tower: Arc::clone(f.tower()),
        parent_field: f.parent_field.as_ref().map(Arc::clone),
    })
}

/// Layer minpoly for T4a sibling adjoin (rational or T3+ parent blocks).
pub(crate) enum LayerMinPolyForAdjoin<'a> {
    Rational(&'a CoordsQ),
    ParentBlocks(&'a [CoordsQ]),
}

// **Pipeline private** — `layer_minpoly_coords_for_adjoin`
pub(crate) fn layer_minpoly_coords_for_adjoin(
    field: &ExtensionField,
) -> Result<LayerMinPolyForAdjoin<'_>, EvalError> {
    if field.is_base() {
        return Err(EvalError::TypeError("base has no layer minpoly"));
    }
    match field.tower.as_ref() {
        ExtensionTower::Adj {
            min_poly_parent_blocks: Some(blocks),
            ..
        } => Ok(LayerMinPolyForAdjoin::ParentBlocks(blocks)),
        ExtensionTower::Adj { min_poly_q, .. } => Ok(LayerMinPolyForAdjoin::Rational(min_poly_q)),
        ExtensionTower::Base => Err(EvalError::TypeError("base has no layer minpoly")),
    }
}

/// Cold flatten over ℚ (no session memo).
/// **Pipeline private** — `flatten_min_poly_over_q_cold`
pub(crate) fn flatten_min_poly_over_q_cold(
    field: &ExtensionField,
    session: Option<&super::field_session::FieldSession>,
) -> Result<CoordsQ, EvalError> {
    match field.tower.as_ref() {
        ExtensionTower::Base => Ok(Vec::new()),
        ExtensionTower::Adj {
            min_poly_parent_blocks: Some(blocks),
            ..
        } => {
            let parent = field
                .parent_field
                .as_ref()
                .ok_or(EvalError::TypeError("tower field without parent"))?;
            let ring = parent_coeff_ring!(parent);
            Ok(char_poly_matrix(&mult_matrix_of_adjoin_generator(
                parent, blocks, &ring,
            )?))
        }
        ExtensionTower::Adj {
            min_poly_q,
            parent,
            ..
        } => {
            if matches!(parent.as_ref(), ExtensionTower::Base) {
                return Ok(min_poly_q.clone());
            }
            let parent_field = field
                .parent_field
                .as_ref()
                .ok_or(EvalError::TypeError("tower field without parent"))?;
            let (min_g, _k) = compose_min_poly_over_q(parent_field, min_poly_q, session)?;
            Ok(min_g)
        }
    }
}

/// Explicit ℚ-flatten minpoly; optional [`FieldSession`] memo (R6).
/// **Pipeline private** — `flatten_min_poly_over_q`
pub(crate) fn flatten_min_poly_over_q(
    field: &ExtensionField,
    session: Option<&super::field_session::FieldSession>,
) -> Result<CoordsQ, EvalError> {
    if field.is_base() {
        return Ok(Vec::new());
    }
    if field.tower().is_simple_over_q() {
        if let ExtensionTower::Adj { min_poly_q, .. } = field.tower.as_ref() {
            return Ok(min_poly_q.clone());
        }
    }
    if let Some(s) = session {
        let key = field.semantic_key();
        if let Some(hit) = s.flatten_cache.borrow().get(&key) {
            return Ok(hit.clone());
        }
        let poly = flatten_min_poly_over_q_cold(field, Some(s))?;
        s.flatten_cache.borrow_mut().insert(key, poly.clone());
        return Ok(poly);
    }
    flatten_min_poly_over_q_cold(field, None)
}

// **Pipeline private** — compositum minpoly via multiplication-matrix charpoly (no k-search).
fn compose_min_poly_via_charpoly(
    parent: &ExtensionField,
    layer_min_poly: &CoordsQ,
) -> Result<CoordsQ, EvalError> {
    let blocks: Vec<CoordsQ> = layer_min_poly
        .iter()
        .map(|r| parent.embed_rational(r))
        .collect();
    let ring = parent_coeff_ring!(parent);
    Ok(char_poly_matrix(&mult_matrix_of_adjoin_generator(
        parent, &blocks, &ring,
    )?))
}

// **Pipeline private** — `compose_min_poly_over_q`
fn compose_min_poly_over_q(
    parent: &Arc<ExtensionField>,
    layer_min_poly: &CoordsQ,
    session: Option<&super::field_session::FieldSession>,
) -> Result<(CoordsQ, Option<i64>), EvalError> {
    let na = parent.dimension();
    let nb = poly_degree(layer_min_poly);
    let expected = na * nb;
    let parent_mp = flatten_min_poly_over_q(parent.as_ref(), session)?;
    for k in 1i64..=PRIMITIVE_K_SEARCH_MAX {
        match common_primitive_sum(
            &parent_mp,
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
    compose_min_poly_via_charpoly(parent.as_ref(), layer_min_poly).and_then(|min_g| {
        if poly_degree(&min_g) != expected {
            return Err(EvalError::NotImplemented("ExtensionField::adjoin compose minpoly"));
        }
        let field = match session {
            Some(s) => s.get_or_create_base_by_min_poly(min_g.clone()),
            None => get_or_create_base_by_min_poly(min_g.clone()),
        };
        if ExtensionField::try_subfield_embedding(parent, &field)?.is_none() {
            return Err(EvalError::NotImplemented("ExtensionField::adjoin compose minpoly"));
        }
        Ok((min_g, None))
    })
}

// **Pipeline private** — `rational_subfield_embedding`
fn rational_subfield_embedding(
    sub: &Arc<ExtensionField>,
    sup: &Arc<ExtensionField>,
) -> FieldEmbedding {
    let dim = sup.dimension();
    let mut m = vec![vec![Ratio::zero(); 1]; dim];
    for (i, c) in sup.embed_rational(&Ratio::one()).into_iter().enumerate() {
        if i < dim {
            m[i][0] = c;
        }
    }
    FieldEmbedding {
        source: Arc::clone(sub),
        target: Arc::clone(sup),
        matrix: m,
    }
}

// **Pipeline private** — `flatten_layer_blocks`
fn flatten_layer_blocks(blocks: &[CoordsQ], parent_dim: usize) -> CoordsQ {
    let mut out = Vec::with_capacity(blocks.len() * parent_dim);
    for block in blocks {
        out.extend(pad_to_len(block, parent_dim).iter().cloned());
    }
    out
}

// **Pipeline private** — `fields_same_parent`
fn fields_same_parent(a: &Arc<ExtensionField>, b: &Arc<ExtensionField>) -> bool {
    Arc::ptr_eq(a, b) || **a == **b
}

// **Pipeline private** — `direct_adjoin_parent_embedding`
fn direct_adjoin_parent_embedding(
    parent: &Arc<ExtensionField>,
    child: &Arc<ExtensionField>,
) -> Result<FieldEmbedding, EvalError> {
    let parent_ok = child
        .parent_field
        .as_ref()
        .is_some_and(|p| fields_same_parent(p, parent));
    if !parent_ok {
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

// **Pipeline private** — `compose_field_embeddings`
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

/// Multiplication-by-layer-generator matrix for U1b `common_minimal`.
// **Pipeline private** — `mult_matrix_layer_gen`
pub(crate) fn mult_matrix_layer_gen(
    parent: &ExtensionField,
    blocks: &[CoordsQ],
) -> Result<Vec<Vec<Ratio<BigInt>>>, EvalError> {
    let ring = parent_coeff_ring!(parent);
    mult_matrix_of_adjoin_generator(parent, blocks, &ring)
}

// **Pipeline private** — `compute_common_dispatch`
fn compute_common_dispatch(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    session: Option<&super::field_session::FieldSession>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    if ExtensionField::is_subfield_of(a, b) || ExtensionField::is_subfield_of(b, a) {
        return subfield_common_pair(a, b);
    }
    if tower_minimal_eligible(a, b) && !common_minimal::compositum_upstream_active() {
        if let Ok(pair) = common_minimal::compute_common_minimal_pair(a, b, session) {
            return Ok(pair);
        }
    }
    #[cfg(feature = "tower-common")]
    if tower_common_eligible(a, b) {
        return compute_common_tower(a, b);
    }
    compute_common_flatten(a, b, session)
}

/// True when upstream compositum pipeline may apply (not subfield).
// **Pipeline private** — `tower_minimal_eligible`
fn tower_minimal_eligible(a: &Arc<ExtensionField>, b: &Arc<ExtensionField>) -> bool {
    !a.is_base()
        && !b.is_base()
        && !ExtensionField::is_subfield_of(a, b)
        && !ExtensionField::is_subfield_of(b, a)
}

/// T4b: `common(sub, sup)` = inclusion into the superfield (no new compositum, no flatten search).
// **Pipeline private** — `subfield_common_pair`
fn subfield_common_pair(
    fa: &Arc<ExtensionField>,
    fb: &Arc<ExtensionField>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    let (sub, sup) = if ExtensionField::is_subfield_of(fa, fb) {
        (fa, fb)
    } else if ExtensionField::is_subfield_of(fb, fa) {
        (fb, fa)
    } else {
        return Err(EvalError::TypeError("subfield_common_pair: unrelated fields"));
    };
    let embed_sub = ExtensionField::try_subfield_embedding(sub, sup)?
        .ok_or(EvalError::TypeError("subfield common: embed missing"))?;
    let dim = sup.dimension();
    let embed_sup = FieldEmbedding {
        source: Arc::clone(sup),
        target: Arc::clone(sup),
        matrix: super::field_arith::identity_matrix(dim),
    };
    let pick = |op: &Arc<ExtensionField>| {
        if Arc::ptr_eq(op, sub) || **op == **sub {
            embed_sub.clone()
        } else {
            embed_sup.clone()
        }
    };
    Ok(Arc::new(CommonFieldPair {
        field: Arc::clone(sup),
        embed_a: pick(fa),
        embed_b: pick(fb),
    }))
}

/// True when T4a tower adjoin can replace flatten `k=1..12` search.
#[cfg(feature = "tower-common")]
// **Pipeline private** — `tower_common_eligible`
fn tower_common_eligible(a: &Arc<ExtensionField>, b: &Arc<ExtensionField>) -> bool {
    !a.is_base()
        && !b.is_base()
        && a.tower().is_simple_over_q()
        && b.tower().is_simple_over_q()
        && !ExtensionField::is_subfield_of(a, b)
        && !ExtensionField::is_subfield_of(b, a)
}

/// Pick the smaller-degree field as adjoin parent (tie-break: `semantic_key` lex).
#[cfg(feature = "tower-common")]
// **Pipeline private** — `pick_tower_adjoin_parent`
fn pick_tower_adjoin_parent<'a>(
    a: &'a Arc<ExtensionField>,
    b: &'a Arc<ExtensionField>,
) -> (&'a Arc<ExtensionField>, &'a Arc<ExtensionField>) {
    match a.dimension().cmp(&b.dimension()) {
        std::cmp::Ordering::Less => (a, b),
        std::cmp::Ordering::Greater => (b, a),
        std::cmp::Ordering::Equal => {
            if a.semantic_key() <= b.semantic_key() {
                (a, b)
            } else {
                (b, a)
            }
        }
    }
}

#[cfg(feature = "tower-common")]
// **Pipeline private** — `embedding_for_common_operand`
fn embedding_for_common_operand(
    operand: &Arc<ExtensionField>,
    parent: &Arc<ExtensionField>,
    embed_parent: &FieldEmbedding,
    embed_sibling: &FieldEmbedding,
) -> FieldEmbedding {
    if Arc::ptr_eq(operand, parent) || **operand == **parent {
        embed_parent.clone()
    } else {
        embed_sibling.clone()
    }
}

/// T4a: compositum as a tower adjoin `parent(sibling.min_poly)` with subfield + sibling embeddings.
///
/// Math notes: see `.doc/giac-tower-common-math.md` (compositum, primitive element, formalization).
#[cfg(feature = "tower-common")]
// **Pipeline private** — `compute_common_tower`
fn compute_common_tower(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    debug_assert!(tower_common_eligible(a, b));
    let (parent, sibling) = pick_tower_adjoin_parent(a, b);
    let common = match layer_minpoly_coords_for_adjoin(sibling)? {
        LayerMinPolyForAdjoin::Rational(layer) => {
            ExtensionField::adjoin_irreducible(parent, layer.to_vec())?
        }
        LayerMinPolyForAdjoin::ParentBlocks(blocks) => {
            ExtensionField::adjoin_irreducible_parent_coeffs(parent, blocks.to_vec())?
        }
    };
    let expected_dim = parent.dimension() * sibling.dimension();
    if common.dimension() != expected_dim {
        return Err(EvalError::TypeError("tower common: compositum degree mismatch"));
    }
    let embed_parent = ExtensionField::try_subfield_embedding(parent, &common)?
        .ok_or(EvalError::TypeError("tower common: parent embed missing"))?;
    let embed_sibling = simple_over_q_embedding(sibling, &common)?;
    Ok(Arc::new(CommonFieldPair {
        field: Arc::clone(&common),
        embed_a: embedding_for_common_operand(a, parent, &embed_parent, &embed_sibling),
        embed_b: embedding_for_common_operand(b, parent, &embed_parent, &embed_sibling),
    }))
}

/// Embed `source` ≅ ℚ(α) into `target` when `target` adjoins the same layer minpoly over a superfield.
// **Pipeline private** — `simple_over_q_embedding`
fn simple_over_q_embedding(
    source: &Arc<ExtensionField>,
    target: &Arc<ExtensionField>,
) -> Result<FieldEmbedding, EvalError> {
    let src_dim = source.dimension();
    let tgt_dim = target.dimension();
    let mut matrix = vec![vec![Ratio::zero(); src_dim]; tgt_dim];
    for j in 0..src_dim {
        let mut unit = source.zero_coords();
        unit[j] = Ratio::one();
        let image = embed_simple_over_q_coords(source, target, &unit)?;
        for i in 0..tgt_dim {
            matrix[i][j] = image[i].clone();
        }
    }
    Ok(FieldEmbedding {
        source: Arc::clone(source),
        target: Arc::clone(target),
        matrix,
    })
}

/// Evaluate `coords` (poly1 in `source`) at the adjoin generator of `target`.
// **Pipeline private** — `embed_simple_over_q_coords`
fn embed_simple_over_q_coords(
    source: &ExtensionField,
    target: &ExtensionField,
    coords: &CoordsQ,
) -> Result<CoordsQ, EvalError> {
    let src_dim = source.dimension();
    let coords = pad_to_len(coords, src_dim);
    let gen = target.generator_coords();
    let mut powers = vec![target.one_coords()];
    for _ in 1..src_dim {
        powers.push(target.element_mul(powers.last().unwrap(), &gen)?);
    }
    let mut acc = target.zero_coords();
    for i in 0..src_dim {
        let c = coords[i].clone();
        if c.is_zero() {
            continue;
        }
        let exp = src_dim - 1 - i;
        let term = target.element_mul(&target.embed_rational(&c), &powers[exp])?;
        acc = target.element_add(&acc, &term)?;
    }
    Ok(acc)
}

#[cfg(all(test, feature = "tower-common"))]
// **Stable** — `tower_adjoin_parent_for_test`
pub(crate) fn tower_adjoin_parent_for_test<'a>(
    a: &'a Arc<ExtensionField>,
    b: &'a Arc<ExtensionField>,
) -> &'a Arc<ExtensionField> {
    pick_tower_adjoin_parent(a, b).0
}

/// Compositum = adjoin(sibling.layer_minpoly) over parent when sibling is simple-over-ℚ.
// **Pipeline private** — `common_adjoin_sibling_over`
fn common_adjoin_sibling_over(
    parent: &Arc<ExtensionField>,
    sibling: &Arc<ExtensionField>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    if !sibling.tower().is_simple_over_q() {
        return Err(EvalError::TypeError("common adjoin: sibling not simple over Q"));
    }
    let common = match layer_minpoly_coords_for_adjoin(sibling)? {
        LayerMinPolyForAdjoin::Rational(layer) => {
            ExtensionField::adjoin_irreducible(parent, layer.to_vec())?
        }
        LayerMinPolyForAdjoin::ParentBlocks(blocks) => {
            ExtensionField::adjoin_irreducible_parent_coeffs(parent, blocks.to_vec())?
        }
    };
    let expected_dim = parent.dimension() * sibling.dimension();
    if common.dimension() != expected_dim {
        return Err(EvalError::TypeError("common adjoin: degree mismatch"));
    }
    let embed_parent = ExtensionField::try_subfield_embedding(parent, &common)?
        .ok_or(EvalError::TypeError("common adjoin: parent embed missing"))?;
    let embed_sibling = simple_over_q_embedding(sibling, &common)?;
    Ok(Arc::new(CommonFieldPair {
        field: Arc::clone(&common),
        embed_a: embed_parent,
        embed_b: embed_sibling,
    }))
}

// **Pipeline private** — nested + simple-over-ℚ compositum when flatten k-search fails
pub(crate) fn try_common_adjoin_one_simple(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    if b.tower().is_simple_over_q()
        && !ExtensionField::is_subfield_of(b, a)
        && !ExtensionField::is_subfield_of(a, b)
    {
        return common_adjoin_sibling_over(a, b);
    }
    if a.tower().is_simple_over_q()
        && !ExtensionField::is_subfield_of(a, b)
        && !ExtensionField::is_subfield_of(b, a)
    {
        let pair = common_adjoin_sibling_over(b, a)?;
        return Ok(Arc::new(CommonFieldPair {
            field: Arc::clone(&pair.field),
            embed_a: pair.embed_b.clone(),
            embed_b: pair.embed_a.clone(),
        }));
    }
    Err(EvalError::TypeError("common adjoin: no eligible simple operand"))
}

/// Phase 0 flatten compositum: primitive element θ = α + k·β, search `k = 1..PRIMITIVE_K_SEARCH_MAX`.
// **Pipeline private** — `compute_common_flatten`
fn compute_common_flatten(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    session: Option<&super::field_session::FieldSession>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    if a.is_base() && !b.is_base() {
        return Ok(embed_rationals_into(b, a, b));
    }
    if b.is_base() && !a.is_base() {
        return Ok(embed_rationals_into(a, b, a));
    }
    let ma = flatten_min_poly_over_q(a, session)?;
    let mb = flatten_min_poly_over_q(b, session)?;
    let na = poly_degree(&ma);
    let nb = poly_degree(&mb);
    for k in 1i64..=PRIMITIVE_K_SEARCH_MAX {
        match common_primitive_sum(&ma, &mb, na, nb, k) {
            Ok((min_g, mat_a, mat_b)) => {
                let field = match session {
                    Some(s) => s.get_or_create_base_by_min_poly(min_g.clone()),
                    None => get_or_create_base_by_min_poly(min_g.clone()),
                };
                let patched = Arc::new(CommonFieldPair {
                    field: Arc::clone(&field),
                    embed_a: FieldEmbedding {
                        source: Arc::clone(a),
                        target: Arc::clone(&field),
                        matrix: mat_a,
                    },
                    embed_b: FieldEmbedding {
                        source: Arc::clone(b),
                        target: Arc::clone(&field),
                        matrix: mat_b,
                    },
                });
                return Ok(patched);
            }
            Err(EvalError::TypeError(_)) => continue,
            Err(e) => return Err(e),
        }
    }
    if let Ok(pair) = try_common_adjoin_one_simple(a, b) {
        return Ok(pair);
    }
    Err(EvalError::NotImplemented("ExtensionField::common"))
}

/// ℚ ↪ K: rationals embed as the constant term of K.
// **Pipeline private** — `embed_rationals_into`
fn embed_rationals_into(
    ext: &Arc<ExtensionField>,
    source_a: &Arc<ExtensionField>,
    source_b: &Arc<ExtensionField>,
) -> Arc<CommonFieldPair> {
    let dim = ext.dimension();
    let rational_embed = || {
        let mut m = vec![vec![Ratio::zero(); 1]; dim];
        for (i, c) in ext.embed_rational(&Ratio::one()).into_iter().enumerate() {
            if i < dim {
                m[i][0] = c;
            }
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

// **Pipeline private** — `common_primitive_sum`
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
// **Pipeline private** — `embedding_matrix_from_theta`
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
// **Pipeline private** — `embedding_matrix_from_theta_block`
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

// **Pipeline private** — `embed_in_gamma_vector`
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
/// **Stable** — `embed_coords`
pub fn embed_coords(embedding: &FieldEmbedding, coords: &CoordsQ) -> CoordsQ {
    embedding.apply(coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::test_fixtures::{
        k1_adjoin_cbrt2, k1_adjoin_sqrt2, k1_adjoin_sqrt3, minpoly_u2_minus,
    };

    #[derive(Clone)]
    struct CoordsInFieldForTest {
        field: Arc<ExtensionField>,
        coords: CoordsQ,
    }

    impl CoordsInFieldForTest {
        fn new(field: Arc<ExtensionField>, coords: CoordsQ) -> Self {
            let coords = pad_to_len(&coords, field.dimension());
            Self { field, coords }
        }

        fn zero(field: &Arc<ExtensionField>) -> Self {
            Self::new(Arc::clone(field), field.zero_coords())
        }

        fn one(field: &Arc<ExtensionField>) -> Self {
            Self::new(Arc::clone(field), field.one_coords())
        }

        fn rational(field: &Arc<ExtensionField>, n: i64) -> Self {
            Self::new(
                Arc::clone(field),
                field.embed_rational(&Ratio::from_integer(n.into())),
            )
        }

        fn generator(field: &Arc<ExtensionField>) -> Self {
            Self::new(Arc::clone(field), field.generator_coords())
        }

        fn add(&self, rhs: &Self) -> Self {
            assert_eq!(*self.field, *rhs.field);
            Self::new(
                Arc::clone(&self.field),
                self.field.element_add(&self.coords, &rhs.coords).unwrap(),
            )
        }

        fn mul(&self, rhs: &Self) -> Self {
            assert_eq!(*self.field, *rhs.field);
            Self::new(
                Arc::clone(&self.field),
                self.field.element_mul(&self.coords, &rhs.coords).unwrap(),
            )
        }

        fn embedded_by(&self, emb: &FieldEmbedding) -> Self {
            assert_eq!(*self.field, *emb.source);
            Self::new(Arc::clone(&emb.target), emb.apply(&self.coords))
        }

        fn eq_mod(&self, rhs: &Self) -> bool {
            assert_eq!(*self.field, *rhs.field);
            self.field.element_eq_mod(&self.coords, &rhs.coords).unwrap()
        }
    }

    fn assert_embedding_ring_hom(
        emb: &FieldEmbedding,
        a: &CoordsInFieldForTest,
        b: &CoordsInFieldForTest,
    ) {
        let target = &emb.target;
        assert!(CoordsInFieldForTest::zero(&emb.source)
            .embedded_by(emb)
            .eq_mod(&CoordsInFieldForTest::zero(target)));
        assert!(CoordsInFieldForTest::one(&emb.source)
            .embedded_by(emb)
            .eq_mod(&CoordsInFieldForTest::one(target)));
        assert!(a
            .add(b)
            .embedded_by(emb)
            .eq_mod(&a.embedded_by(emb).add(&b.embedded_by(emb))));
        assert!(a
            .mul(b)
            .embedded_by(emb)
            .eq_mod(&a.embedded_by(emb).mul(&b.embedded_by(emb))));
    }

    #[test]
    fn r6_nested_adjoin_layer_two_flatten_explicit_four() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let layer = match layer_minpoly_coords_for_adjoin(&k2).unwrap() {
            LayerMinPolyForAdjoin::Rational(p) => p,
            LayerMinPolyForAdjoin::ParentBlocks(_) => panic!("expected rational layer"),
        };
        assert_eq!(poly_degree(layer), 2);
        let flat = flatten_min_poly_over_q_cold(&k2, None).unwrap();
        assert_eq!(poly_degree(&flat), 4);
    }

    #[test]
    fn compose_minpoly_charpoly_matches_nested_sqrt2_sqrt3() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let layer = match layer_minpoly_coords_for_adjoin(&k2).unwrap() {
            LayerMinPolyForAdjoin::Rational(p) => p.to_vec(),
            LayerMinPolyForAdjoin::ParentBlocks(_) => panic!("expected rational layer"),
        };
        let (via_k, _) = compose_min_poly_over_q(&k1, &layer, None).unwrap();
        let via_char = compose_min_poly_via_charpoly(&k1, &layer).unwrap();
        assert_eq!(poly_degree(&via_k), 4);
        assert_eq!(poly_degree(&via_char), 4);
        // ponytail: charpoly minpoly may differ from θ=α+kβ primitive poly; same degree suffices for fallback.
    }

    #[test]
    fn common_adjoin_nested_sqrt2_sqrt3_with_sqrt5() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let k5 = ExtensionField::adjoin_irreducible(
            &ExtensionField::rational(),
            minpoly_u2_minus(-5),
        )
        .unwrap();
        let pair = common_adjoin_sibling_over(&k2, &k5).unwrap();
        assert_eq!(pair.field.dimension(), 8);
        let emb_k5 = ExtensionField::embedding_for(&k5, &pair).unwrap();
        let b = emb_k5.apply(&k5.generator_coords());
        let sq_b = pair.field.element_mul(&b, &b).unwrap();
        let five = pair.field.embed_rational(&Ratio::from_integer(5.into()));
        assert!(pair.field.element_eq_mod(&sq_b, &five).unwrap());
        let emb_k2 = ExtensionField::embedding_for(&k2, &pair).unwrap();
        let a = emb_k2.apply(&k2.generator_coords());
        let sum = pair.field.element_add(&a, &b).unwrap();
        assert!(!pair.field.element_is_zero(&sum));
    }

    // U0 — compositum verify gate (see `.doc/issues/GIAC-ext-common-unified-path.md`).
    mod u0_verify_common_pair {
        use super::*;

        #[test]
        fn t4a_sqrt2_sqrt3_passes() {
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let pair = ExtensionField::common_over_q(&k1, &k3).unwrap();
            super::super::verify_common_pair(&pair, &k1, &k3).expect("legacy T4a");
        }

        #[test]
        fn subfield_k1_in_k2_passes() {
            let k1 = k1_adjoin_sqrt2();
            let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
            let pair = ExtensionField::common_over_q(&k1, &k2).unwrap();
            super::super::verify_common_pair(&pair, &k1, &k2).expect("T4b subfield");
        }

        #[test]
        fn adjoin_nested_sqrt2_sqrt3_with_sqrt5_passes() {
            let k1 = k1_adjoin_sqrt2();
            let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
            let k5 = ExtensionField::adjoin_irreducible(
                &ExtensionField::rational(),
                minpoly_u2_minus(-5),
            )
            .unwrap();
            let pair = common_adjoin_sibling_over(&k2, &k5).unwrap();
            super::super::verify_common_pair(&pair, &k2, &k5).expect("common adjoin");
        }

        #[test]
        fn flatten_k_search_nested_sqrt5_fails_verify() {
            let k1 = k1_adjoin_sqrt2();
            let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
            let k5 = ExtensionField::adjoin_irreducible(
                &ExtensionField::rational(),
                minpoly_u2_minus(-5),
            )
            .unwrap();
            let pair = compute_common_flatten_for_test(&k2, &k5).unwrap();
            assert!(
                super::super::verify_common_pair(&pair, &k2, &k5).is_err(),
                "flatten k-search must not pass U0 gate on Q(sqrt2,sqrt3) x Q(sqrt5)"
            );
        }
    }

    // U1a — upstream `common_minimal_POLY` simple×simple (see GIAC-ext-common-unified-path.md).
    mod u1a_common_minimal {
        use super::*;
        use super::common_minimal::compute_common_minimal_pair_for_test;

        #[test]
        fn sqrt2_sqrt3_passes_verify() {
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let pair = compute_common_minimal_pair_for_test(&k1, &k3).expect("U1a compositum");
            verify_common_pair(&pair, &k1, &k3).expect("U0 gate");
            assert_eq!(pair.field.dimension(), 4);
        }

        #[test]
        fn sqrt2_cbrt2_passes_verify() {
            let k2 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_cbrt2();
            let pair = compute_common_minimal_pair_for_test(&k2, &k3).expect("U1a compositum");
            verify_common_pair(&pair, &k2, &k3).expect("U0 gate");
            assert_eq!(pair.field.dimension(), 6);
        }

        #[test]
        fn dispatch_prefers_minimal_over_t4a() {
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let pair = ExtensionField::common_over_q(&k1, &k3).unwrap();
            verify_common_pair(&pair, &k1, &k3).expect("dispatch compositum");
            assert!(pair.field.parent_field().is_none(), "U1a flat compositum");
        }
    }

    mod u1b_common_minimal {
        use super::*;
        use super::common_minimal::{
            common_minimal_poly_over_parent_for_test,
        };

        #[test]
        fn over_parent_sqrt2_u2_minus_alpha_degree_four() {
            let k1 = k1_adjoin_sqrt2();
            let one = k1.one_coords();
            let zero = k1.zero_coords();
            let neg_alpha = k1.element_neg(&k1.generator_coords()).unwrap();
            let layer = vec![one, zero, neg_alpha];
            let ma = simple_layer_minpoly_from_field(&k1);
            let (min_g, _k, _w_a, _w_b) =
                common_minimal_poly_over_parent_for_test(&k1, &ma, &layer).expect("U1b core");
            assert_eq!(poly_degree(&min_g), 4);
        }

        fn simple_layer_minpoly_from_field(f: &ExtensionField) -> CoordsQ {
            match layer_minpoly_coords_for_adjoin(f).unwrap() {
                LayerMinPolyForAdjoin::Rational(m) => m.to_vec(),
                LayerMinPolyForAdjoin::ParentBlocks(_) => panic!("expected rational"),
            }
        }

        #[test]
        fn mrref_spec_tensor_dims_consistent() {
            use super::super::common_minimal::{
                compositum_session_active, compute_compositum_upstream_for_test,
            };
            assert!(!compositum_session_active());
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let spec = compute_compositum_upstream_for_test(&k1, &k3).expect("U1a spec");
            assert_eq!(spec.dim(), 4);
            assert_eq!(spec.na() * spec.nb(), spec.dim());
            assert_eq!(spec.mat_theta().len(), spec.dim());
            assert!(!compositum_session_active());
        }

        #[test]
        fn verify_k5_flat_k2_core_embeddings_pass_v2() {
            use super::super::common_minimal::compute_compositum_upstream_for_test;
            let k1 = k1_adjoin_sqrt2();
            let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
            let k5 = ExtensionField::adjoin_irreducible(
                &ExtensionField::rational(),
                minpoly_u2_minus(-5),
            )
            .unwrap();
            let spec = compute_compositum_upstream_for_test(&k5, &k2).expect("U2 minimal core");
            assert_eq!(poly_degree(spec.min_g()), 8);
            assert_eq!(spec.dim(), 8);
            assert_eq!(spec.na() * spec.nb(), spec.dim());
            let pair = super::super::common_minimal::compute_common_minimal_pair_for_test(
                &k2, &k5,
            )
            .expect("compositum with adjoin fallback");
            verify_common_pair(&pair, &k2, &k5).expect("U0 gate");
        }

        #[test]
        fn flatten_k2_sqrt2_sqrt3_is_fast() {
            use super::super::flatten_min_poly_over_q;
            let k1 = k1_adjoin_sqrt2();
            let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
            let flat = flatten_min_poly_over_q(k2.as_ref(), None).unwrap();
            assert_eq!(poly_degree(&flat), 4);
        }

        #[test]
        fn nested_sqrt5_upstream_passes_verify() {
            let k1 = k1_adjoin_sqrt2();
            let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
            let k5 = ExtensionField::adjoin_irreducible(
                &ExtensionField::rational(),
                minpoly_u2_minus(-5),
            )
            .unwrap();
            let pair = super::super::common_minimal::compute_common_minimal_pair_for_test(&k2, &k5)
                .expect("U2 upstream");
            verify_common_pair(&pair, &k2, &k5).expect("U0 gate");
            assert_eq!(pair.field.dimension(), 8);
        }

        #[test]
        fn dispatch_nested_sqrt5_prefers_upstream() {
            let k1 = k1_adjoin_sqrt2();
            let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
            let k5 = ExtensionField::adjoin_irreducible(
                &ExtensionField::rational(),
                minpoly_u2_minus(-5),
            )
            .unwrap();
            let pair = ExtensionField::common_over_q(&k2, &k5).expect("dispatch common");
            assert_eq!(pair.field.dimension(), 8);
            verify_common_pair(&pair, &k2, &k5).expect("U0 gate");
        }
    }

    mod u2_factor_select {
        use super::*;
        use crate::algebra::poly_alg_factor::{
            factor_minpoly_over_field, select_factor_for_common,
        };

        #[test]
        fn factor_x4_minus_4_over_sqrt2_selects_quadratic() {
            let k1 = k1_adjoin_sqrt2();
            // (x²-2)² = x⁴ - 4x² + 4; over Q(√2): x⁴-4 = (x²-2)(x²+2)
            let mb = vec![
                Ratio::one(),
                Ratio::zero(),
                Ratio::zero(),
                Ratio::zero(),
                Ratio::from_integer((-4).into()),
            ];
            let factors = factor_minpoly_over_field(&mb, &k1).expect("factor");
            assert!(factors.len() >= 2);
            let k3 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-2)).unwrap();
            let picked = select_factor_for_common(&factors, &k3, &k1).expect("select");
            assert_eq!(
                picked.degree_wrt(&crate::algebra::poly_alg_factor::common_minpoly_var()),
                2,
            );
        }
    }

    #[test]
    fn r6_parent_coeff_adjoin_flatten_explicit_four() {
        let k1 = k1_adjoin_sqrt2();
        let zero = k1.zero_coords();
        let one = k1.one_coords();
        let neg_alpha = k1.element_neg(&k1.generator_coords()).unwrap();
        let layer = vec![one, zero, neg_alpha];
        let k2 = ExtensionField::adjoin_irreducible_parent_coeffs(&k1, layer).unwrap();
        assert_eq!(k2.dimension(), 4);
        let flat = flatten_min_poly_over_q_cold(&k2, None).unwrap();
        assert_eq!(poly_degree(&flat), 4);
    }

    #[test]
    fn rational_field_dimension_one() {
        let q = ExtensionField::rational();
        assert_eq!(q.dimension(), 1);
        assert!(q.element_is_one(&q.one_coords()));
    }

    #[test]
    fn t3a_adjoin_k1_u2_minus_sqrt2_has_dimension_four() {
        let (k1, k2) = crate::algebra::test_fixtures::t3a_k2_adjoin_u2_minus_sqrt2_over_k1();
        assert_eq!(k2.dimension(), 4);
        assert_eq!(k2.parent_field(), Some(&k1));
        let u = k2.generator_coords();
        let u_sq = k2.element_mul(&u, &u).unwrap();
        let alpha = k1.generator_coords();
        let alpha_in_k2 = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .expect("K1 embeds in K2")
            .apply(&alpha);
        assert!(k2.element_eq_mod(&u_sq, &alpha_in_k2).unwrap());
        let u4 = k2.element_mul(&u_sq, &u_sq).unwrap();
        let two = k2.embed_rational(&Ratio::from_integer(2.into()));
        assert!(k2.element_eq_mod(&u4, &two).unwrap());
        assert_eq!(poly_degree(&flatten_min_poly_over_q_cold(&k2, None).unwrap()), 4);
    }

    #[test]
    fn rational_embedding_into_tower_uses_constant_block() {
        let (_k1, k2) = crate::algebra::test_fixtures::t3a_k2_adjoin_u2_minus_sqrt2_over_k1();
        let q = ExtensionField::rational();
        let emb = ExtensionField::try_subfield_embedding(&q, &k2)
            .unwrap()
            .expect("Q embeds in tower");
        let one = emb.apply(&q.one_coords());
        assert_eq!(one, k2.embed_rational(&Ratio::one()));
    }

    #[test]
    fn embedding_ring_hom_rational_to_tower() {
        let (_k1, k2) = crate::algebra::test_fixtures::t3a_k2_adjoin_u2_minus_sqrt2_over_k1();
        let q = ExtensionField::rational();
        let emb = ExtensionField::try_subfield_embedding(&q, &k2)
            .unwrap()
            .expect("Q embeds in tower");
        let two = CoordsInFieldForTest::rational(&q, 2);
        let three = CoordsInFieldForTest::rational(&q, 3);
        assert_embedding_ring_hom(&emb, &two, &three);
    }

    #[test]
    fn common_rational_to_tower_embedding_uses_constant_block() {
        let (_k1, k2) = crate::algebra::test_fixtures::t3a_k2_adjoin_u2_minus_sqrt2_over_k1();
        let q = ExtensionField::rational();
        let pair = ExtensionField::common_over_q(&q, &k2).unwrap();
        let emb = ExtensionField::embedding_for(&q, &pair).unwrap();
        let one = CoordsInFieldForTest::one(&q).embedded_by(emb);
        assert!(one.eq_mod(&CoordsInFieldForTest::one(&k2)));
    }

    #[test]
    fn embedding_ring_hom_parent_to_child() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let emb = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .expect("K1 embeds in child");
        let alpha = CoordsInFieldForTest::generator(&k1);
        let one = CoordsInFieldForTest::one(&k1);
        let alpha_plus_one = alpha.add(&one);
        assert_embedding_ring_hom(&emb, &alpha, &alpha_plus_one);
    }

    #[test]
    fn embedding_ring_hom_composite_chain() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let k3 = ExtensionField::adjoin_irreducible(&k2, minpoly_u2_minus(-5)).unwrap();
        let emb = ExtensionField::try_subfield_embedding(&k1, &k3)
            .unwrap()
            .expect("K1 embeds through composite parent chain");
        let alpha = CoordsInFieldForTest::generator(&k1);
        let two = CoordsInFieldForTest::rational(&k1, 2);
        let alpha_plus_two = alpha.add(&two);
        assert_embedding_ring_hom(&emb, &alpha, &alpha_plus_two);
    }

    #[test]
    fn adjoin_cbrt2_generator_cubes_to_two() {
        let k = k1_adjoin_cbrt2();
        assert_eq!(k.dimension(), 3);
        let alpha = k.generator_coords();
        let a2 = k.element_mul(&alpha, &alpha).unwrap();
        let a3 = k.element_mul(&a2, &alpha).unwrap();
        let two = k.embed_rational(&Ratio::from_integer(2.into()));
        assert!(k.element_eq_mod(&a3, &two).unwrap());
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

    // **B** — F5 S3: √2 in ℚ(√2,√3) via subfield descent (no new adjoin).
    #[test]
    fn try_square_root_sqrt2_in_q_sqrt2_sqrt3() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        assert_eq!(k2.dimension(), 4);
        let two = k2.embed_rational(&Ratio::from_integer(2.into()));
        let root = ExtensionField::try_square_root_in_field(&k2, &two)
            .unwrap()
            .expect("sqrt(2) in Q(sqrt2,sqrt3) via S3");
        let sq = k2.element_mul(&root, &root).unwrap();
        assert!(k2.element_eq_mod(&sq, &two).unwrap());
    }

    // **B** — F5 S5: √6 = √2·√3 in ℚ(√2,√3).
    #[test]
    fn try_square_root_sqrt6_in_q_sqrt2_sqrt3() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let six = k2.embed_rational(&Ratio::from_integer(6.into()));
        let root = ExtensionField::try_square_root_in_field(&k2, &six)
            .unwrap()
            .expect("sqrt(6) in Q(sqrt2,sqrt3) via S5");
        let sq = k2.element_mul(&root, &root).unwrap();
        assert!(k2.element_eq_mod(&sq, &six).unwrap());
    }

    // **B** — F5: √8 in ℚ(√2) without new adjoin.
    #[test]
    fn try_square_root_sqrt8_in_q_sqrt2() {
        let k = k1_adjoin_sqrt2();
        let eight = k.embed_rational(&Ratio::from_integer(8.into()));
        let root = ExtensionField::try_square_root_in_field(&k, &eight)
            .unwrap()
            .expect("sqrt(8) in Q(sqrt2)");
        let sq = k.element_mul(&root, &root).unwrap();
        assert!(k.element_eq_mod(&sq, &eight).unwrap());
    }

    // **B** — F5: √2 in ℚ(√2) via layer generator scan.
    #[test]
    fn try_square_root_sqrt2_in_q_sqrt2() {
        let k = k1_adjoin_sqrt2();
        let two = k.embed_rational(&Ratio::from_integer(2.into()));
        let root = ExtensionField::try_square_root_in_field(&k, &two)
            .unwrap()
            .expect("sqrt(2) in Q(sqrt2)");
        let sq = k.element_mul(&root, &root).unwrap();
        assert!(k.element_eq_mod(&sq, &two).unwrap());
    }

    // **B** — F5: parent-coeff tower; subfield √2 recognized after embed.
    #[test]
    fn try_square_root_in_parent_coeff_tower() {
        let k1 = k1_adjoin_sqrt2();
        let one = k1.one_coords();
        let sqrt2 = k1.generator_coords();
        let three = k1.embed_rational(&Ratio::from_integer(3.into()));
        let field =
            ExtensionField::adjoin_irreducible_parent_coeffs(&k1, vec![one, sqrt2, three])
                .unwrap();
        let two = field.embed_rational(&Ratio::from_integer(2.into()));
        let root = ExtensionField::try_square_root_in_field(&field, &two)
            .unwrap()
            .expect("sqrt(2) in parent-coeff tower");
        let sq = field.element_mul(&root, &root).unwrap();
        assert!(field.element_eq_mod(&sq, &two).unwrap());
    }

    // **B** — F5: align then probe √ in common superfield (tower compositum; flat U1a S-path TBD).
    #[test]
    #[cfg(feature = "tower-common")]
    fn try_square_root_after_align_in_common() {
        let k1 = k1_adjoin_sqrt2();
        let k3 = k1_adjoin_sqrt3();
        let common = compute_common_tower(&k1, &k3).unwrap();
        let emb = ExtensionField::embedding_for(&k1, &common).unwrap();
        let alpha = emb.apply(&k1.generator_coords());
        let sq_alpha = common.field.element_mul(&alpha, &alpha).unwrap();
        let root = ExtensionField::try_square_root_in_field(&common.field, &sq_alpha)
            .unwrap()
            .expect("sqrt(alpha^2) in common field");
        let sq = common.field.element_mul(&root, &root).unwrap();
        assert!(common.field.element_eq_mod(&sq, &sq_alpha).unwrap());
    }

    #[test]
    #[cfg_attr(
        not(feature = "tower-common"),
        ignore = "flatten common(√2,∛2) is slow; default enables tower-common (T4a)"
    )]
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
        assert_eq!(*c1.field, *c2.field);
    }

    #[test]
    #[cfg_attr(
        not(feature = "tower-common"),
        ignore = "flatten common(√2,∛2) is slow; default enables tower-common (T4a)"
    )]
    fn common_cache_hits_same_pair() {
        use super::super::field_session::FieldSession;

        let session = FieldSession::new(ExtensionField::rational());
        let a = k1_adjoin_sqrt2();
        let b = k1_adjoin_cbrt2();
        let c1 = session.common_over_q(&a, &b).unwrap();
        assert_eq!(session.common_cache_len(), 1);
        let c2 = session.common_over_q(&b, &a).unwrap();
        assert_eq!(session.common_cache_len(), 1);
        assert_eq!(*c1.field, *c2.field);
    }

    #[test]
    fn adjoin_cache_dedup_same_minpoly() {
        use super::super::field_session::FieldSession;

        let q = ExtensionField::rational();
        let session = FieldSession::new(Arc::clone(&q));
        let k1 = session
            .adjoin_irreducible(&q, minpoly_u2_minus(-2))
            .unwrap();
        let k2 = session
            .adjoin_irreducible(&q, minpoly_u2_minus(-2))
            .unwrap();
        assert_eq!(*k1, *k2);
        assert!(Arc::ptr_eq(&k1, &k2));
    }

    #[test]
    fn r1_duplicate_field_arc_subfield_embedding_parent_semantic_match() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let dup_k1 = duplicate_field_arc_for_test(&k1);
        assert!(!Arc::ptr_eq(&k1, &dup_k1));
        assert_eq!(*k1, *dup_k1);
        let emb = ExtensionField::try_subfield_embedding(&dup_k1, &k2)
            .unwrap()
            .expect("semantic-equal parent embeds in child");
        let orig = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .expect("original parent embeds");
        let alpha = k1.generator_coords();
        assert!(k2
            .element_eq_mod(&emb.apply(&alpha), &orig.apply(&alpha))
            .unwrap());
    }

    #[test]
    #[cfg_attr(
        not(feature = "tower-common"),
        ignore = "tower-common required for semantic common cache hit"
    )]
    fn r1_duplicate_field_arc_common_cache_semantic_key() {
        use super::super::field_session::FieldSession;

        let session = FieldSession::new(ExtensionField::rational());
        let sqrt2 = k1_adjoin_sqrt2();
        let dup_sqrt2 = duplicate_field_arc_for_test(&sqrt2);
        let cbrt2 = k1_adjoin_cbrt2();
        let _ = session.common_over_q(&sqrt2, &cbrt2).unwrap();
        let before = session.common_cache_len();
        let c2 = session.common_over_q(&dup_sqrt2, &cbrt2).unwrap();
        assert_eq!(session.common_cache_len(), before);
        assert_eq!(sqrt2.semantic_key(), dup_sqrt2.semantic_key());
        assert_eq!(c2.field.dimension(), 6);
    }

    #[test]
    fn embedding_for_matches_source_not_operand_order() {
        let sqrt2 = k1_adjoin_sqrt2();
        let cbrt2 = k1_adjoin_cbrt2();
        let (low, high) = if sqrt2.semantic_key() <= cbrt2.semantic_key() {
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
        assert!(Arc::ptr_eq(
            &ExtensionField::embedding_for(low, &pair).unwrap().source,
            low
        ));
        assert!(Arc::ptr_eq(
            &ExtensionField::embedding_for(high, &pair).unwrap().source,
            high
        ));
    }

    #[test]
    fn align_elements_reverse_order_after_cache_warm() {
        use super::super::field_session::FieldSession;

        let session = FieldSession::new(ExtensionField::rational());
        let q = ExtensionField::rational();
        let sqrt2 = k1_adjoin_sqrt2();
        let _ = session.common_over_q(&q, &sqrt2).unwrap();

        let one_q = q.one_coords();
        let alpha = sqrt2.generator_coords();

        let qr = session
            .align_elements(&q, &one_q, &sqrt2, &alpha)
            .unwrap();
        let rq = session
            .align_elements(&sqrt2, &alpha, &q, &one_q)
            .unwrap();

        assert_eq!(*qr.field, *rq.field);
        assert_eq!(qr.field, sqrt2);
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
        assert_eq!(k2.parent_field(), Some(&k1));
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
        assert!(Arc::ptr_eq(&emb.source, &k1));
        assert!(Arc::ptr_eq(&emb.target, &k2));
        assert_eq!(emb.matrix.len(), k2.dimension());
        assert_eq!(emb.matrix[0].len(), k1.dimension());
        let img = emb.apply(&k1.generator_coords());
        assert_eq!(img.len(), k2.dimension());
    }

    #[test]
    fn t4b_common_subfield_is_superfield_not_compositum() {
        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let pair = ExtensionField::common_over_q(&k1, &k2).unwrap();
        assert_eq!(*pair.field, *k2);
        assert_eq!(pair.field.parent_field(), Some(&k1));
        let alpha = k1.generator_coords();
        let emb = ExtensionField::embedding_for(&k1, &pair).unwrap();
        let img = emb.apply(&alpha);
        let alpha_k2 = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .unwrap()
            .apply(&alpha);
        assert!(pair.field.element_eq_mod(&img, &alpha_k2).unwrap());
    }

    #[test]
    fn t2_align_subfield_does_not_grow_common_cache() {
        use super::super::field_session::FieldSession;

        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let session = FieldSession::new(Arc::clone(&k1));
        let before = session.common_cache_len();
        let alpha = k1.generator_coords();
        let alpha_in_k2 = ExtensionField::try_subfield_embedding(&k1, &k2)
            .unwrap()
            .unwrap()
            .apply(&alpha);
        let _ = session
            .align_elements(&k1, &alpha, &k2, &alpha_in_k2)
            .unwrap();
        assert_eq!(session.common_cache_len(), before);
    }

    #[test]
    fn t2_algext_add_same_generator_in_superfield_no_common_cache() {
        use super::super::field_session::FieldSession;
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
        let session = FieldSession::new(Arc::clone(&k1));
        let before = session.common_cache_len();
        let sum = a.add(&b).unwrap();
        assert_eq!(session.common_cache_len(), before);
        assert_eq!(sum.field, a.field);
        let two_b = b.add(&b).unwrap();
        assert!(sum.eq_mod(&two_b).unwrap());
    }

    #[test]
    fn k2_embedded_sqrt2_squared_is_two() {
        use crate::algebra::test_fixtures::t1b_k2_adjoin_sqrt3_over_k1;
        let fix = t1b_k2_adjoin_sqrt3_over_k1();
        let sq = fix
            .k2
            .element_mul(&fix.sqrt2_in_k2, &fix.sqrt2_in_k2)
            .unwrap();
        let two = fix
            .k2
            .embed_rational(&Ratio::from_integer(2.into()));
        assert!(fix.k2.element_eq_mod(&sq, &two).unwrap());
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

    #[cfg(not(feature = "tower-common"))]
    mod flatten_bisect {
        use super::*;
        use crate::algebra::test_fixtures::{k1_adjoin_sqrt2, k1_adjoin_sqrt3};

        /// Fast smoke: `--no-default-features` still builds and common(Q,√2) works.
        #[test]
        fn common_q_sqrt2_without_tower_common() {
            let q = ExtensionField::rational();
            let sqrt2 = k1_adjoin_sqrt2();
            let pair = ExtensionField::common_over_q(&q, &sqrt2).unwrap();
            assert_eq!(pair.field, sqrt2);
        }

        /// DoD: `--no-default-features` restores Phase 0 flatten compositum shape (slow).
        #[test]
        #[ignore = "flatten common(√2,√3) slow; cargo test -p giac-core --no-default-features flatten_bisect -- --ignored"]
        fn common_sqrt2_sqrt3_has_flatten_parent_none() {
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let pair = ExtensionField::common_over_q(&k1, &k3).unwrap();
            assert_eq!(pair.field.dimension(), 4);
            assert!(pair.field.parent_field().is_none());
        }
    }

    #[cfg(feature = "tower-common")]
    mod t4a {
        //! Legacy T4a dispatch — superseded by [GIAC-ext-common-unified-path.md](../../../../.doc/issues/GIAC-ext-common-unified-path.md) U3.
        use super::*;
        use crate::algebra::test_fixtures::{k1_adjoin_cbrt2, k1_adjoin_sqrt2, k1_adjoin_sqrt3};

    #[test]
        fn common_sqrt2_sqrt3_tower_adjoin_has_parent() {
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            // Legacy T4a path (dispatch now prefers U1a flat compositum).
            let pair = compute_common_tower(&k1, &k3).unwrap();
            assert_eq!(pair.field.dimension(), 4);
            let expected_parent = tower_adjoin_parent_for_test(&k1, &k3);
            assert_eq!(pair.field.parent_field(), Some(expected_parent));
        }

    #[test]
        fn align_sqrt2_plus_sqrt3_reverse_order() {
            use crate::algebra::field_session::FieldSession;

            let session = FieldSession::new(ExtensionField::rational());
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let _ = session.common_over_q(&k1, &k3).unwrap();
            let alpha = k1.generator_coords();
            let beta = k3.generator_coords();
            let ab = session
                .align_elements(&k1, &alpha, &k3, &beta)
                .unwrap();
            let ba = session
                .align_elements(&k3, &beta, &k1, &alpha)
                .unwrap();
            assert_eq!(*ab.field, *ba.field);
            let sum_ab = ab.field.element_add(&ab.left, &ab.right).unwrap();
            let sum_ba = ba.field.element_add(&ba.left, &ba.right).unwrap();
            assert!(ab.field.element_eq_mod(&sum_ab, &sum_ba).unwrap());
            assert!(!ab.field.element_is_zero(&sum_ab));
        }

    #[test]
        fn common_cache_hits_after_tower_common() {
            use crate::algebra::field_session::FieldSession;

            let session = FieldSession::new(ExtensionField::rational());
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let c1 = session.common_over_q(&k1, &k3).unwrap();
            assert_eq!(session.common_cache_len(), 1);
            let c2 = session.common_over_q(&k3, &k1).unwrap();
            assert_eq!(session.common_cache_len(), 1);
            assert_eq!(*c1.field, *c2.field);
        }

    #[test]
        fn align_sqrt2_cbrt2_dim_six() {
            let sqrt2 = k1_adjoin_sqrt2();
            let cbrt2 = k1_adjoin_cbrt2();
            let alpha = sqrt2.generator_coords();
            let beta = cbrt2.generator_coords();
            let aligned =
                ExtensionField::align_elements(&sqrt2, &alpha, &cbrt2, &beta).unwrap();
            assert_eq!(aligned.field.dimension(), 6);
            let sum = aligned
                .field
                .element_add(&aligned.left, &aligned.right)
                .unwrap();
            assert!(!aligned.field.element_is_zero(&sum));
        }

    #[test]
        fn tower_common_invariants_sqrt2_sqrt3() {
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let pair = ExtensionField::common_over_q(&k1, &k3).unwrap();
            assert_eq!(pair.field.dimension(), k1.dimension() * k3.dimension());
            let alpha = k1.generator_coords();
            let beta = k3.generator_coords();
            let a = ExtensionField::embedding_for(&k1, &pair).unwrap().apply(&alpha);
            let b = ExtensionField::embedding_for(&k3, &pair).unwrap().apply(&beta);
            let two = pair.field.embed_rational(&Ratio::from_integer(2.into()));
            assert!(pair
                .field
                .element_eq_mod(&pair.field.element_mul(&a, &a).unwrap(), &two)
                .unwrap());
            let three = pair.field.embed_rational(&Ratio::from_integer(3.into()));
            assert!(pair
                .field
                .element_eq_mod(&pair.field.element_mul(&b, &b).unwrap(), &three)
                .unwrap());
        }

    #[test]
        #[ignore = "flatten char-poly cross-check is slow; cargo test -p giac-core -- --ignored"]
        fn tower_common_matches_flatten_minpoly_on_sqrt2_sqrt3() {
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let tower = ExtensionField::common_over_q(&k1, &k3).unwrap();
            let flatten = compute_common_flatten_for_test(&k1, &k3).unwrap();
            assert_eq!(tower.field.dimension(), flatten.field.dimension());
            assert_eq!(
                flatten_min_poly_over_q(&tower.field, None).unwrap(),
                flatten_min_poly_over_q(&flatten.field, None).unwrap()
            );
            let alpha = k1.generator_coords();
            let beta = k3.generator_coords();
            let t_left = ExtensionField::embedding_for(&k1, &tower)
                .unwrap()
                .apply(&alpha);
            let t_right = ExtensionField::embedding_for(&k3, &tower)
                .unwrap()
                .apply(&beta);
            let f_left = ExtensionField::embedding_for(&k1, &flatten)
                .unwrap()
                .apply(&alpha);
            let f_right = ExtensionField::embedding_for(&k3, &flatten)
                .unwrap()
                .apply(&beta);
            let t_sum = tower.field.element_add(&t_left, &t_right).unwrap();
            let f_sum = flatten.field.element_add(&f_left, &f_right).unwrap();
            assert!(!tower.field.element_is_zero(&t_sum));
            assert!(!flatten.field.element_is_zero(&f_sum));
            let t_alpha_sq = tower.field.element_mul(&t_left, &t_left).unwrap();
            let two = tower.field.embed_rational(&Ratio::from_integer(2.into()));
            assert!(tower.field.element_eq_mod(&t_alpha_sq, &two).unwrap());
            let f_alpha_sq = flatten.field.element_mul(&f_left, &f_left).unwrap();
            assert!(flatten.field.element_eq_mod(&f_alpha_sq, &two).unwrap());
        }
    }
}
