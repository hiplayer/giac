//! Ambient **K** + working **L** for `poly_algext_roots` (G5 / PR-B′).
//!
//! Coefficients are lifted and aligned on **L**; `ambient` stays fixed for normalized
//! polynomial input. See [FieldSession plan](../../../../.doc/issues/GIAC-poly-roots-field-session-plan.md).
//!
//! ## Session 契约（F4 / M4）
//!
//! | 不变量 | 规则 |
//! |--------|------|
//! | **ambient K** | 创建后不变；`normalize_coeffs` 将输入系数 embed 到 K |
//! | **working L** | 单调扩域：`bump_to` / `lift` / `align` / `adjoin_*` 仅当新域 ⊇ 当前 L 时增大 L |
//! | **checkpoint** | [`Self::checkpoint`] 快照 L；[`Self::restore`] 回退到快照（放弃更大分支） |
//! | **lift / align** | 两操作数对齐到同一 **L** 后再算术；`align` 后 `bump_to` 到公共超域 |
//! | **开方 / adjoin** | 生产路径经 [`Self::sqrt_in_field`]（先 [`Self::try_sqrt_in_field`] → [`ExtensionField::try_square_root_in_field`]）；[`Self::adjoin_sqrt_new`] 仅 blind fallback |
//!
//! 四次 `t⁴+t+1` 维数门禁（CI 默认）：resolvent 阶段 `dim(L)≤6`；全流程分裂域硬顶 `dim(L)≤24`
//!（S₄）；见 `poly_roots::tests::field_session_dimension_bound_quartic`。
//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::cell::RefCell;

/// Splitting-field hard cap for deg≤4 root pipelines (F1/F4).
pub(crate) const POLY_ROOTS_DIM_HARD: usize = 24;
/// Resolvent cubic stage bound at \(d_K=1\) (F2 A′ + deflate for \(p\neq0\)).
pub(crate) const POLY_ROOTS_DIM_RESOLVENT: usize = 6;
/// Reference degree when Gal(f) ≅ A₄ ([K_f:Q]=12); **not** a CI gate for all quartics (see algorithm spec §5.4).
pub(crate) const POLY_ROOTS_DIM_QUARTIC_OUT: usize = 12;
use std::collections::HashMap;
use std::rc::Rc;
use std::sync::Arc;

use giac_poly::PolyCoeff;
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;

use super::alg_ext::{algext_cube_root, AlgExtData};
use super::alg_ext_c::AlgExtCData;
use super::ext_tower::{self, AlignedElements, CommonCache, ExtensionField};
use super::field_arith::{coords_to_expr, poly_degree, rationalize_poly1, trim_leading_zero, CoordsQ, pad_to_len};
use super::poly_alg_coeff::AlgExtCPolyCoeff;

type AdjoinCache = HashMap<Vec<u8>, Arc<ExtensionField>>;
type BaseByMinPoly = HashMap<Vec<u8>, Arc<ExtensionField>>;
pub(crate) type FlattenCache = HashMap<Vec<u8>, CoordsQ>;

/// Session carrying ambient coefficient field **K** and monotonically growing working field **L**.
/// **Stable (bounded)** — roots pipeline field session
pub struct FieldSession {
    ambient: Arc<ExtensionField>,
    working: RefCell<Arc<ExtensionField>>,
    /// R2: session-local `common_over_q` cache (replaces global registry cache).
    common_cache: Rc<RefCell<CommonCache>>,
    /// R4: optional adjoin dedup (replaces `FieldRegistry::by_adjoin` / `by_min_poly`).
    adjoin_cache: Rc<RefCell<AdjoinCache>>,
    base_by_min_poly: Rc<RefCell<BaseByMinPoly>>,
    /// R6: memo for explicit `flatten_min_poly_over_q`.
    pub(crate) flatten_cache: Rc<RefCell<FlattenCache>>,
}

impl FieldSession {
    /// **Stable** — `FieldSession::new`
    pub fn new(ambient: Arc<ExtensionField>) -> Self {
        let working = Arc::clone(&ambient);
        Self {
            ambient,
            working: RefCell::new(working),
            common_cache: Rc::new(RefCell::new(HashMap::new())),
            adjoin_cache: Rc::new(RefCell::new(HashMap::new())),
            base_by_min_poly: Rc::new(RefCell::new(HashMap::new())),
            flatten_cache: Rc::new(RefCell::new(HashMap::new())),
        }
    }

    /// Fresh K/L while sharing extension caches with `parent` (R5 eval / roots).
    /// **Stable (bounded)** — fork ambient with shared extension cache
    pub fn fork_ambient(&self, ambient: Arc<ExtensionField>) -> Self {
        Self {
            ambient: Arc::clone(&ambient),
            working: RefCell::new(Arc::clone(&ambient)),
            common_cache: Rc::clone(&self.common_cache),
            adjoin_cache: Rc::clone(&self.adjoin_cache),
            base_by_min_poly: Rc::clone(&self.base_by_min_poly),
            flatten_cache: Rc::clone(&self.flatten_cache),
        }
    }

    /// **Pipeline private** — run closure on common cache (reentrant with adjoin caches).
    pub(crate) fn with_common_cache<R>(&self, f: impl FnOnce(&mut CommonCache) -> R) -> R {
        f(&mut self.common_cache.borrow_mut())
    }

    /// **Pipeline private** — clear extension caches (tests)
    #[cfg(test)]
    pub(crate) fn clear_caches_for_test(&self) {
        self.common_cache.borrow_mut().clear();
        self.adjoin_cache.borrow_mut().clear();
        self.base_by_min_poly.borrow_mut().clear();
        self.flatten_cache.borrow_mut().clear();
    }

    /// **Pipeline private** — explicit ℚ-flatten with session memo (R6).
    pub(crate) fn flatten_min_poly_over_q(
        &self,
        field: &ExtensionField,
    ) -> Result<CoordsQ, EvalError> {
        ext_tower::flatten_min_poly_over_q(field, Some(self))
    }

    pub(crate) fn flatten_cache_len(&self) -> usize {
        self.flatten_cache.borrow().len()
    }

    /// Adjoin with session-local dedup (R4).
    /// **Stable (bounded)** — adjoin irreducible
    pub fn adjoin_irreducible(
        &self,
        parent: &Arc<ExtensionField>,
        min_poly_q: CoordsQ,
    ) -> Result<Arc<ExtensionField>, EvalError> {
        let min_poly_q = trim_leading_zero(min_poly_q);
        if poly_degree(&min_poly_q) < 1 {
            return Err(EvalError::TypeError("extension degree >= 1"));
        }
        if !parent.is_base() && !ext_tower::layer_minpoly_rational_constants(&min_poly_q) {
            return Err(EvalError::TypeError(
                "T1 adjoin: layer minpoly must have rational coefficients",
            ));
        }
        let key = ext_tower::adjoin_cache_key_rational(parent, &min_poly_q);
        if let Some(hit) = self.adjoin_cache.borrow().get(&key) {
            return Ok(Arc::clone(hit));
        }
        let mp_key = ext_tower::min_poly_key_bytes(&min_poly_q);
        if parent.is_base() {
            if let Some(hit) = self.base_by_min_poly.borrow().get(&mp_key) {
                self.adjoin_cache
                    .borrow_mut()
                    .insert(key, Arc::clone(hit));
                return Ok(Arc::clone(hit));
            }
        }
        let field = ext_tower::build_adjoin_irreducible(parent, min_poly_q)?;
        self.adjoin_cache
            .borrow_mut()
            .insert(key, Arc::clone(&field));
        if parent.is_base() {
            self.base_by_min_poly
                .borrow_mut()
                .insert(mp_key, Arc::clone(&field));
        }
        Ok(field)
    }

    /// **Stable (bounded)** — adjoin with parent-coeff layer minpoly
    pub fn adjoin_irreducible_parent_coeffs(
        &self,
        parent: &Arc<ExtensionField>,
        layer_blocks: Vec<CoordsQ>,
    ) -> Result<Arc<ExtensionField>, EvalError> {
        if !parent.parent_coeff_ring_capable() {
            return Err(EvalError::TypeError(
                "parent not suitable for parent-coeff adjoin",
            ));
        }
        if parent.is_base() {
            return Err(EvalError::TypeError(
                "parent-coeff adjoin requires nontrivial parent",
            ));
        }
        let nb = layer_blocks
            .len()
            .checked_sub(1)
            .ok_or(EvalError::TypeError("layer minpoly"))?;
        if nb < 1 {
            return Err(EvalError::TypeError("extension degree >= 1"));
        }
        let key = ext_tower::adjoin_cache_key_parent_blocks(parent, &layer_blocks);
        if let Some(hit) = self.adjoin_cache.borrow().get(&key) {
            return Ok(Arc::clone(hit));
        }
        let field = ext_tower::build_adjoin_parent_coeffs(parent, layer_blocks)?;
        self.adjoin_cache
            .borrow_mut()
            .insert(key, Arc::clone(&field));
        Ok(field)
    }

    /// Adjoin a layer over current **L** with parent-field minpoly coefficients.
    /// **Stable (bounded)** — parent-coeff adjoin on working field
    pub fn adjoin_parent_coeff_layer(
        &self,
        layer_blocks: Vec<CoordsQ>,
    ) -> Result<Arc<ExtensionField>, EvalError> {
        let parent = self.working();
        self.adjoin_irreducible_parent_coeffs(&parent, layer_blocks)
    }

    pub(crate) fn adjoin_cache_len(&self) -> usize {
        self.adjoin_cache.borrow().len()
    }

    /// Dedup ℚ(θ) by minimal polynomial bytes (flatten compositum path).
    pub(crate) fn get_or_create_base_by_min_poly(
        &self,
        min_poly_q: CoordsQ,
    ) -> Arc<ExtensionField> {
        let key = ext_tower::min_poly_key_bytes(&min_poly_q);
        if let Some(hit) = self.base_by_min_poly.borrow().get(&key) {
            return Arc::clone(hit);
        }
        let field = ext_tower::build_base_extension_uncached(min_poly_q);
        self.base_by_min_poly
            .borrow_mut()
            .insert(key, Arc::clone(&field));
        field
    }

    /// **Stable** — ambient K (normalized poly coefficients)
    pub fn ambient(&self) -> &Arc<ExtensionField> {
        &self.ambient
    }

    /// **Stable** — current working L
    pub fn working(&self) -> Arc<ExtensionField> {
        Arc::clone(&self.working.borrow())
    }

    /// Entries in the session `common_over_q` cache (tests / diagnostics).
    /// **Stable (bounded)** — common cache length
    pub fn common_cache_len(&self) -> usize {
        self.common_cache.borrow().len()
    }

    /// **Stable (bounded)** — common extension with session cache
    pub fn common_over_q(
        &self,
        a: &Arc<ExtensionField>,
        b: &Arc<ExtensionField>,
    ) -> Result<Arc<ext_tower::CommonFieldPair>, EvalError> {
        ext_tower::common_over_q_in_cache(a, b, &mut self.common_cache.borrow_mut(), Some(self))
    }

    /// Align two field elements; uses this session's common cache (R2).
    /// **Stable (bounded)** — align coords on session cache
    pub fn align_elements(
        &self,
        a_field: &Arc<ExtensionField>,
        a_coords: &super::field_arith::CoordsQ,
        b_field: &Arc<ExtensionField>,
        b_coords: &super::field_arith::CoordsQ,
    ) -> Result<AlignedElements, EvalError> {
        ext_tower::align_elements_in_cache(
            a_field,
            a_coords,
            b_field,
            b_coords,
            &mut self.common_cache.borrow_mut(),
            Some(self),
        )
    }

    /// **Stable** — zero in L
    pub fn zero(&self) -> AlgExtCPolyCoeff {
        let w = self.working();
        coeff_in_field(&w, |f| AlgExtCData::zero(Arc::clone(f)).expect("zero"))
    }

    /// **Stable** — one in L
    pub fn one(&self) -> AlgExtCPolyCoeff {
        let w = self.working();
        coeff_in_field(&w, |f| AlgExtCData::one(Arc::clone(f)).expect("one"))
    }

    /// **Stable** — integer constant in L
    pub fn int(&self, n: i64) -> Result<AlgExtCPolyCoeff, EvalError> {
        rat(
            &self.working(),
            Ratio::from_integer(BigInt::from(n)),
        )
    }

    /// **Stable** — 1/2 in L
    pub fn half(&self) -> Result<AlgExtCPolyCoeff, EvalError> {
        rat(&self.working(), Ratio::new(1.into(), 2.into()))
    }

    /// Embed `c` into current **L** (may bump working).
    /// **Stable** — lift coeff to working field
    pub fn lift(&self, c: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let working = self.working();
        if Arc::ptr_eq(&c.as_inner().field, &working) {
            return Ok(c.clone());
        }
        if ExtensionField::is_subfield_of(&c.as_inner().field, &working) {
            let z = self.zero();
            let (lifted, _) = self.align(c, &z)?;
            return Ok(lifted);
        }
        if ExtensionField::is_subfield_of(&working, &c.as_inner().field) {
            self.bump_to(&c.as_inner().field);
            return Ok(c.clone());
        }
        let z = self.zero();
        let (lifted, _) = self.align(c, &z)?;
        Ok(lifted)
    }

    /// Align two coefficients to a common field; bumps **L** to the align target.
    /// **Stable** — align pair on working field
    pub fn align(
        &self,
        a: &AlgExtCPolyCoeff,
        b: &AlgExtCPolyCoeff,
    ) -> Result<(AlgExtCPolyCoeff, AlgExtCPolyCoeff), EvalError> {
        let (aa, bb) = AlgExtCData::align_pair_in_cache(
            a.as_inner(),
            b.as_inner(),
            &mut self.common_cache.borrow_mut(),
        )?;
        self.bump_to(&aa.field);
        Ok((AlgExtCPolyCoeff::from(aa), AlgExtCPolyCoeff::from(bb)))
    }

    /// Formal i·z with z ∈ K embedded in current **L** (i² = −1, not adjoined to the tower).
    /// **Stable (bounded)** — multiply by formal i
    pub fn mul_formal_i(&self, z: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let z = self.lift(z)?;
        let inner = z.as_inner();
        Ok(AlgExtCPolyCoeff::from(AlgExtCData {
            field: Arc::clone(&inner.field),
            re: coords_to_expr(&inner.field.zero_coords())?,
            im: inner.re.clone(),
            root_index: None,
        }))
    }

    /// Canonical principal square root on **L** (may extend the tower).
    ///
    /// ℚ negative → i√|u|; real u → [`sqrt_in_field`].
    /// **Stable (bounded)** — principal sqrt on working field
    pub fn sqrt_principal(&self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let u = self.lift(u)?;
        if u.coeff_is_zero() {
            return Ok(self.zero());
        }
        if is_negative_rational(&u) {
            let abs = self.neg(&u)?;
            let beta = self.sqrt_in_field(&abs)?;
            return self.mul_formal_i(&beta);
        }
        self.sqrt_in_field(&u)
    }

    /// Return ε ∈ **L** with ε² ≡ u, adjoining only when u is not already a square in L.
    /// **Stable (bounded)** — sqrt in field or adjoin
    pub fn sqrt_in_field(&self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let u = self.lift(u)?;
        if let Some(beta) = self.try_sqrt_in_field(&u)? {
            let sq = self.mul(&beta, &beta)?;
            let (sq_a, u_a) = self.align(&sq, &u)?;
            debug_assert!(
                sq_a.coeff_sub(&u_a).unwrap().coeff_is_zero(),
                "sqrt_in_field: beta^2 must eq_mod u"
            );
            return Ok(beta);
        }
        if self.working().dimension() >= POLY_ROOTS_DIM_HARD {
            return Err(EvalError::NotImplemented("poly roots dim bound"));
        }
        self.adjoin_sqrt_new(&u)
    }

    /// If u is a square in **L**, return some square root (no tower extension).
    /// Delegates to [`ExtensionField::try_square_root_in_field`] on **L**.
    /// **Stable (bounded)** — try existing square root
    pub fn try_sqrt_in_field(
        &self,
        u: &AlgExtCPolyCoeff,
    ) -> Result<Option<AlgExtCPolyCoeff>, EvalError> {
        let u = self.lift(u)?;
        let inner = u.as_inner();
        if !inner.im.iter().all(|e| e.is_zero()) {
            return Ok(None);
        }
        if u.coeff_is_zero() {
            return Ok(Some(self.zero()));
        }
        let field = self.working();
        let re = rationalize_poly1(&inner.re)?;
        let coords = pad_to_len(&re, field.dimension());
        if let Some(root_coords) = ExtensionField::try_square_root_in_field(&field, &coords)? {
            let beta = AlgExtData::from_field_coords(
                Arc::clone(&field),
                coords_to_expr(&root_coords)?,
            )?;
            return Ok(Some(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&beta)?)));
        }
        Ok(None)
    }

    /// Shallow [`Self::try_sqrt_in_field`]: deg-2 layer scan + basis squares only.
    /// **Pipeline private** — Euler second-sqrt fast path before blind adjoin
    pub(crate) fn try_sqrt_in_field_shallow(
        &self,
        u: &AlgExtCPolyCoeff,
    ) -> Result<Option<AlgExtCPolyCoeff>, EvalError> {
        let u = self.lift(u)?;
        let inner = u.as_inner();
        if !inner.im.iter().all(|e| e.is_zero()) {
            return Ok(None);
        }
        if u.coeff_is_zero() {
            return Ok(Some(self.zero()));
        }
        let field = self.working();
        let re = rationalize_poly1(&inner.re)?;
        let coords = pad_to_len(&re, field.dimension());
        if let Some(root_coords) = ExtensionField::try_square_root_in_field_shallow(&field, &coords)? {
            let beta = AlgExtData::from_field_coords(
                Arc::clone(&field),
                coords_to_expr(&root_coords)?,
            )?;
            return Ok(Some(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&beta)?)));
        }
        Ok(None)
    }

    /// Adjoin √u (real part) and return a square root generator in **L**.
    /// **Stable (bounded)** — adjoin sqrt primitive
    pub fn adjoin_sqrt(&self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        self.sqrt_in_field(u)
    }

    /// Always adjoin a new √u layer (skips [`Self::try_sqrt_in_field`]).
    /// **Pipeline private** — blind adjoin sqrt
    pub(crate) fn adjoin_sqrt_new(
        &self,
        u: &AlgExtCPolyCoeff,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let u = self.lift(u)?;
        let inner = u.as_inner();
        if !inner.im.iter().all(|e| e.is_zero()) {
            return Err(EvalError::NotImplemented("adjoin_sqrt complex"));
        }
        let parent = self.working();
        let sqrt_field = if parent.is_base() {
            let re = rationalize_poly1(&inner.re)?;
            let u_val = pad_to_len(&re, 1)[0].clone();
            self.adjoin_irreducible(
                &parent,
                vec![Ratio::one(), Ratio::from_integer(0.into()), -u_val],
            )?
        } else {
            let one = parent.one_coords();
            let zero = parent.zero_coords();
            let u_coords = coords_in_field(self, &u, &parent)?;
            let neg = parent.element_neg(&u_coords)?;
            self.adjoin_parent_coeff_layer(vec![one, zero, neg])?
        };
        let beta = coeff_from_coords(&sqrt_field, &sqrt_field.generator_coords())?;
        self.bump_to(&sqrt_field);
        Ok(beta)
    }

    /// Adjoin ∛u (real part) and return a cube root generator in **L**.
    /// **Stable (bounded)** — adjoin cbrt primitive
    pub fn adjoin_cbrt(&self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let u = self.lift(u)?;
        let inner = u.as_inner();
        if !inner.im.iter().all(|e| e.is_zero()) {
            return Err(EvalError::NotImplemented("adjoin_cbrt complex"));
        }
        let re = AlgExtData::from_field_coords(Arc::clone(&inner.field), inner.re.clone())?;
        let beta = algext_cube_root(&re)?;
        let out = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&beta)?);
        self.bump_to(&out.as_inner().field);
        Ok(out)
    }

    /// Adjoin a primitive cube root of unity ω (ω²+ω+1=0) over **L**.
    /// **Stable (bounded)** — adjoin ω for pure cubic roots
    pub fn adjoin_primitive_cube_root_of_unity(
        &self,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let parent = self.working();
        let one = parent.one_coords();
        let layer = vec![one.clone(), one.clone(), one];
        let omega_field = if parent.dimension() == 1 {
            let one_r = Ratio::one();
            self.adjoin_irreducible(&parent, vec![one_r.clone(), one_r.clone(), one_r])?
        } else {
            self.adjoin_irreducible_parent_coeffs(&parent, layer)?
        };
        let coords = omega_field.generator_coords();
        let omega = AlgExtData::from_field_coords(
            Arc::clone(&omega_field),
            coords_to_expr(&coords)?,
        )?;
        let out = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&omega)?);
        self.bump_to(&out.as_inner().field);
        Ok(out)
    }

    /// Snapshot current **L** for later [`Self::restore`].
    /// **Stable (bounded)** — working-field checkpoint
    pub fn checkpoint(&self) -> Arc<ExtensionField> {
        self.working()
    }

    /// Reset **L** to a prior [`Self::checkpoint`] (abandons extensions after that point).
    /// **Pipeline private** — restore working field
    pub(crate) fn restore(&self, field: &Arc<ExtensionField>) {
        self.set_working(field);
    }

    /// Reset **L** to a prior checkpoint (Euler branch search).
    /// **Pipeline private** — `set_working`
    pub(crate) fn set_working(&self, field: &Arc<ExtensionField>) {
        *self.working.borrow_mut() = Arc::clone(field);
    }

    /// **Stable** — add with auto-align
    pub fn add(
        &self,
        a: &AlgExtCPolyCoeff,
        b: &AlgExtCPolyCoeff,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let (a, b) = self.align(a, b)?;
        a.coeff_add(&b)
    }

    /// **Stable** — multiply with auto-align
    pub fn mul(
        &self,
        a: &AlgExtCPolyCoeff,
        b: &AlgExtCPolyCoeff,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let (a, b) = self.align(a, b)?;
        a.coeff_mul(&b)
    }

    /// **Stable** — divide with auto-align
    pub fn div(
        &self,
        a: &AlgExtCPolyCoeff,
        b: &AlgExtCPolyCoeff,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let (a, b) = self.align(a, b)?;
        a.coeff_div(&b)
    }

    /// **Stable** — negate (no align)
    pub fn neg(&self, a: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        a.coeff_neg()
    }

    /// Grow **L** toward `new` when `new` is a superfield or align common field.
    /// **Pipeline private** — `bump_to`
    pub(crate) fn bump_to(&self, new: &Arc<ExtensionField>) {
        let mut working = self.working.borrow_mut();
        if Arc::ptr_eq(&*working, new) {
            return;
        }
        if ExtensionField::is_subfield_of(&working, new) {
            *working = Arc::clone(new);
            return;
        }
        if ExtensionField::is_subfield_of(new, &working) {
            return;
        }
        if new.dimension() >= working.dimension() {
            *working = Arc::clone(new);
        }
    }
}

// **Pipeline private** — rational constant in field
fn rat(field: &Arc<ExtensionField>, r: Ratio<BigInt>) -> Result<AlgExtCPolyCoeff, EvalError> {
    let coords = field.embed_rational(&r);
    let a = AlgExtData::from_field_coords(Arc::clone(field), coords_to_expr(&coords)?)?;
    Ok(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&a)?))
}

fn coeff_in_field(
    field: &Arc<ExtensionField>,
    f: impl FnOnce(&Arc<ExtensionField>) -> AlgExtCData,
) -> AlgExtCPolyCoeff {
    AlgExtCPolyCoeff::from(f(field))
}

// **Pipeline private** — negative constant in ℚ ⊂ K (for Δ<0 guard)
pub(crate) fn is_negative_rational(c: &AlgExtCPolyCoeff) -> bool {
    let inner = c.as_inner();
    if !inner.im.iter().all(|e| e.is_zero()) {
        return false;
    }
    let Ok(re) = rationalize_poly1(&inner.re) else {
        return false;
    };
    if inner.field.dimension() != 1 {
        return false;
    }
    pad_to_len(&re, 1)[0] < Ratio::zero()
}

// **Pipeline private** — embed real coeff coords into `target`.
pub(crate) fn coords_in_field(
    session: &FieldSession,
    c: &AlgExtCPolyCoeff,
    target: &Arc<ExtensionField>,
) -> Result<CoordsQ, EvalError> {
    let inner = c.as_inner();
    if !inner.im.iter().all(|e| e.is_zero()) {
        return Err(EvalError::TypeError("expected real coefficient"));
    }
    let re = rationalize_poly1(&inner.re)?;
    let aligned = session.align_elements(&inner.field, &re, target, &target.zero_coords())?;
    Ok(pad_to_len(&aligned.left, target.dimension()))
}

// **Pipeline private** — embed coords as coeff in `field`.
pub(crate) fn coeff_from_coords(
    field: &Arc<ExtensionField>,
    coords: &CoordsQ,
) -> Result<AlgExtCPolyCoeff, EvalError> {
    let data = AlgExtData::from_field_coords(
        Arc::clone(field),
        coords_to_expr(&pad_to_len(coords, field.dimension()))?,
    )?;
    Ok(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&data)?))
}

#[cfg(test)]
mod tests {
    use giac_poly::PolyCoeff;

    use crate::algebra::test_fixtures::{k1_adjoin_sqrt2, k1_adjoin_sqrt3, sqrt2_algext};

    use super::*;

    #[test]
    fn r2_session_common_cache_hit_on_second_common() {
        let mut session = FieldSession::new(ExtensionField::rational());
        let sqrt2 = k1_adjoin_sqrt2();
        let sqrt3 = k1_adjoin_sqrt3();
        let _ = session.common_over_q(&sqrt2, &sqrt3).unwrap();
        let before = session.common_cache_len();
        assert!(before >= 1);
        let _ = session.common_over_q(&sqrt3, &sqrt2).unwrap();
        assert_eq!(session.common_cache_len(), before);
    }

    #[test]
    fn restore_checkpoint_discards_later_adjoin() {
        let session = FieldSession::new(ExtensionField::rational());
        let cp = session.checkpoint();
        let two = session.int(2).unwrap();
        let _ = session.adjoin_sqrt(&two).unwrap();
        assert!(session.working().dimension() > 1);
        session.restore(&cp);
        assert_eq!(session.working().dimension(), 1);
        let beta = session.sqrt_in_field(&two).unwrap();
        let sq = session.mul(&beta, &beta).unwrap();
        let (sq_a, two_a) = session.align(&sq, &two).unwrap();
        assert!(sq_a.coeff_sub(&two_a).unwrap().coeff_is_zero());
    }

    #[test]
    fn int_and_one_on_rational() {
        let k = ExtensionField::rational();
        let mut session = FieldSession::new(k);
        let two = session.int(2).unwrap();
        let one = session.one();
        assert!(!two.coeff_is_zero());
        assert!(!one.coeff_is_zero());
        assert!(Arc::ptr_eq(&session.working(), session.ambient()));
    }

    #[test]
    fn lift_and_align_bump_working_on_k1() {
        let k1 = k1_adjoin_sqrt2();
        let mut session = FieldSession::new(Arc::clone(&k1));
        let sqrt2 = AlgExtCPolyCoeff::from(
            AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap(),
        );
        let four = session.int(4).unwrap();
        let disc = session.mul(&sqrt2, &four).unwrap();
        let dim_before = session.working().dimension();
        let beta = session.adjoin_sqrt(&disc).unwrap();
        assert!(
            session.working().dimension() > dim_before,
            "adjoin_sqrt should grow working field"
        );
        let sq = session.mul(&beta, &beta).unwrap();
        let (sq_a, disc_a) = session.align(&sq, &disc).unwrap();
        assert!(sq_a.coeff_sub(&disc_a).unwrap().coeff_is_zero());
    }

    #[test]
    fn adjoin_sqrt_matches_sqrt_layer() {
        let k1 = k1_adjoin_sqrt2();
        let mut session = FieldSession::new(Arc::clone(&k1));
        let sqrt2 = AlgExtCPolyCoeff::from(
            AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap(),
        );
        let four = session.int(4).unwrap();
        let u = session.mul(&sqrt2, &four).unwrap();
        let beta = session.adjoin_sqrt(&u).unwrap();
        let sq = session.mul(&beta, &beta).unwrap();
        let (sq_a, u_a) = session.align(&sq, &u).unwrap();
        assert!(sq_a.coeff_sub(&u_a).unwrap().coeff_is_zero());
    }

    #[test]
    fn adjoin_parent_coeff_quadratic_with_linear_term() {
        let k1 = k1_adjoin_sqrt2();
        let session = FieldSession::new(Arc::clone(&k1));
        let one = k1.one_coords();
        let sqrt2 = k1.generator_coords();
        let three = k1.embed_rational(&Ratio::from_integer(3.into()));
        let field = session
            .adjoin_parent_coeff_layer(vec![one, sqrt2, three])
            .unwrap();
        session.bump_to(&field);
        let beta = coeff_from_coords(&field, &field.generator_coords()).unwrap();
        let beta2 = session.mul(&beta, &beta).unwrap();
        let sqrt2 = AlgExtCPolyCoeff::from(
            AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap(),
        );
        let linear = session.mul(&sqrt2, &beta).unwrap();
        let three = session.int(3).unwrap();
        let sum = session
            .add(&session.add(&beta2, &linear).unwrap(), &three)
            .unwrap();
        assert!(
            sum.as_inner()
                .eq_mod(session.zero().as_inner())
                .unwrap_or(false),
            "generator must satisfy beta^2 + sqrt2*beta + 3"
        );
    }

    #[test]
    fn adjoin_quadratic_over_parent_coeff_cubic_parent() {
        let k1 = k1_adjoin_sqrt2();
        let session = FieldSession::new(Arc::clone(&k1));
        let one = k1.one_coords();
        let zero = k1.zero_coords();
        let sqrt2 = k1.generator_coords();
        let cubic = session
            .adjoin_parent_coeff_layer(vec![one.clone(), zero, sqrt2, one])
            .unwrap();
        session.bump_to(&cubic);
        let gamma = coeff_from_coords(&cubic, &cubic.generator_coords()).unwrap();
        let gamma_coords = coords_in_field(&session, &gamma, &cubic).unwrap();
        let three = cubic.embed_rational(&Ratio::from_integer(3.into()));
        let quad = session
            .adjoin_parent_coeff_layer(vec![cubic.one_coords(), gamma_coords, three])
            .unwrap();
        session.bump_to(&quad);
        let beta = coeff_from_coords(&quad, &quad.generator_coords()).unwrap();
        let beta2 = session.mul(&beta, &beta).unwrap();
        let linear = session.mul(&gamma, &beta).unwrap();
        let sum = session
            .add(&session.add(&beta2, &linear).unwrap(), &session.int(3).unwrap())
            .unwrap();
        assert!(
            sum.as_inner()
                .eq_mod(session.zero().as_inner())
                .unwrap_or(false),
            "generator must satisfy beta^2 + gamma*beta + 3"
        );
    }

    #[test]
    fn r6_flatten_cache_hit_same_semantic_key() {
        use crate::algebra::test_fixtures::minpoly_u2_minus;

        let k1 = k1_adjoin_sqrt2();
        let k2 = ExtensionField::adjoin_irreducible(&k1, minpoly_u2_minus(-3)).unwrap();
        let dup = ext_tower::duplicate_field_arc_for_test(&k2);
        let session = FieldSession::new(ExtensionField::rational());
        assert_eq!(session.flatten_cache_len(), 0);
        let _ = session.flatten_min_poly_over_q(&k2).unwrap();
        assert_eq!(session.flatten_cache_len(), 1);
        let _ = session.flatten_min_poly_over_q(&dup).unwrap();
        assert_eq!(session.flatten_cache_len(), 1);
    }

    #[test]
    fn r6_simple_over_q_skips_flatten_cache() {
        let session = FieldSession::new(ExtensionField::rational());
        let sqrt2 = k1_adjoin_sqrt2();
        let _ = session.flatten_min_poly_over_q(&sqrt2).unwrap();
        assert_eq!(session.flatten_cache_len(), 0);
    }

    // **B** — F4: restore after adjoin abandons the expanded branch.
    #[test]
    fn set_working_restores_adjoin() {
        let k1 = k1_adjoin_sqrt2();
        let session = FieldSession::new(Arc::clone(&k1));
        let cp = session.checkpoint();
        let dim_cp = cp.dimension();
        assert_eq!(dim_cp, 2);
        let three = session.int(3).unwrap();
        let _ = session.sqrt_in_field(&three).unwrap();
        assert!(
            session.working().dimension() > dim_cp,
            "sqrt(3) over Q(sqrt2) should extend L"
        );
        session.restore(&cp);
        assert_eq!(session.working().dimension(), dim_cp);
        assert!(Arc::ptr_eq(&session.working(), &cp));
    }
}
