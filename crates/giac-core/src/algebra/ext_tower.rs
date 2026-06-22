//! Tower of algebraic extensions K₀=ℚ → K₁ → … → Kₙ and cross-field `common`.
//!
//! Normative model: [GIAC-algext-adoption.md](../../../../.doc/issues/GIAC-algext-adoption.md) §8.2.
//!
//! Phase 0 scope: coefficients in ℚ only; each [`ExtensionField`] is described by a
//! primitive/minimal polynomial over ℚ (single adjoin or `common` composite).
//!
//! T1+: nontrivial adjoins record a true [`ExtensionTower::Adj`] parent chain and
//! `(parent_id, min_poly)` registry keys; see [GIAC-lazy-common-tower-plan.md] T1.

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

use super::field_arith::{
    apply_linear_map, char_poly_matrix, coords_all_zero, embed_in_square_extension,
    generator_coords, kron_left, mat_mul, mult_matrix_of_adjoin_generator,
    mult_matrix_of_element, pad_to_len, poly_add, poly_degree, poly_inv_mod, poly_mul, poly_neg,
    poly_reduce, poly_sub, poly_add_with_coeffs_in_field, poly_inv_mod_with_coeffs_in_field,
    poly_mul_with_coeffs_in_field, poly_neg_with_coeffs_in_field,
    poly_reduce_with_coeffs_in_field, ParentCoeffRing, trim_leading_zero, CoordsQ,
};

static FIELD_ID: AtomicU64 = AtomicU64::new(1);

// **Pipeline private** — `next_field_id`
fn next_field_id() -> u64 {
    FIELD_ID.fetch_add(1, Ordering::Relaxed)
}

// **Pipeline private** — one layer-minpoly coefficient (parent-field element) as `Expr`
fn layer_minpoly_block_to_expr(block: &CoordsQ) -> Result<ExprArc, EvalError> {
    use super::field_arith::{canonical_poly1_expr, coords_to_expr};
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
        self.parent_field
            .as_ref()
            .is_some_and(|parent| !parent.is_base())
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
    fn layer_ext_degree(&self) -> Result<usize, EvalError> {
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
            ExtensionTower::Base => {
                return Err(EvalError::TypeError("expected extension field"));
            }
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

/// **Pipeline private** — T1 layer minpoly coefficient check.
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
enum LayerMinPolyForAdjoin<'a> {
    Rational(&'a CoordsQ),
    ParentBlocks(&'a [CoordsQ]),
}

// **Pipeline private** — `layer_minpoly_coords_for_adjoin`
fn layer_minpoly_coords_for_adjoin(
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

// **Pipeline private** — `compose_min_poly_over_q`
fn compose_min_poly_over_q(
    parent: &ExtensionField,
    layer_min_poly: &CoordsQ,
    session: Option<&super::field_session::FieldSession>,
) -> Result<(CoordsQ, Option<i64>), EvalError> {
    let na = parent.dimension();
    let nb = poly_degree(layer_min_poly);
    let expected = na * nb;
    let parent_mp = flatten_min_poly_over_q(parent, session)?;
    for k in 1i64..=12 {
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
    Err(EvalError::NotImplemented("ExtensionField::adjoin compose minpoly"))
}

// **Pipeline private** — `rational_subfield_embedding`
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

// **Pipeline private** — `compute_common_dispatch`
fn compute_common_dispatch(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    session: Option<&super::field_session::FieldSession>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    if ExtensionField::is_subfield_of(a, b) || ExtensionField::is_subfield_of(b, a) {
        return subfield_common_pair(a, b);
    }
    #[cfg(feature = "tower-common")]
    if tower_common_eligible(a, b) {
        return compute_common_tower(a, b);
    }
    compute_common_flatten(a, b, session)
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
#[cfg(feature = "tower-common")]
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
#[cfg(feature = "tower-common")]
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

/// Phase 0 flatten compositum: primitive element θ = α + k·β, search `k = 1..12`.
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
    for k in 1i64..=12 {
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
        k1_adjoin_cbrt2, k1_adjoin_sqrt2, minpoly_u2_minus,
    };

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

        let mut session = FieldSession::new(ExtensionField::rational());
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
        let mut session = FieldSession::new(Arc::clone(&k1));
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
        use super::*;
        use crate::algebra::test_fixtures::{k1_adjoin_cbrt2, k1_adjoin_sqrt2, k1_adjoin_sqrt3};

    #[test]
        fn common_sqrt2_sqrt3_has_tower_parent() {
            let k1 = k1_adjoin_sqrt2();
            let k3 = k1_adjoin_sqrt3();
            let pair = ExtensionField::common_over_q(&k1, &k3).unwrap();
            assert_eq!(pair.field.dimension(), 4);
            let expected_parent = tower_adjoin_parent_for_test(&k1, &k3);
            assert_eq!(pair.field.parent_field(), Some(expected_parent));
            assert!(pair.field.parent_field().is_some());
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
