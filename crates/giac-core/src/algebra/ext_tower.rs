//! Tower of algebraic extensions K₀=ℚ → K₁ → … → Kₙ and cross-field `common`.
//!
//! Normative model: [GIAC-algext-adoption.md](../../../../.doc/issues/GIAC-algext-adoption.md) §8.2.
//!
//! Phase 0 scope: coefficients in ℚ only; each [`ExtensionField`] is described by a
//! primitive/minimal polynomial over ℚ (single adjoin or `common` composite).

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
    poly_degree, poly_inv_mod, poly_mul, poly_neg, poly_reduce, poly_sub, trim_leading_zero,
    CoordsQ,
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
    /// K = parent(γ); `min_poly_q` is the minimal polynomial of γ over ℚ (Phase 0).
    ///
    /// For towers built via [`ExtensionField::common_over_q`], `min_poly_q` is the
    /// primitive-element polynomial over ℚ describing the entire field.
    Adj {
        parent: Arc<ExtensionTower>,
        min_poly_q: CoordsQ,
        ext_degree: usize,
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
}

impl PartialEq for ExtensionField {
    fn eq(&self, other: &Self) -> bool {
        self.tower == other.tower && self.min_poly_over_q == other.min_poly_over_q
    }
}

impl Eq for ExtensionField {}

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

    /// Build ℚ(α) from an irreducible (or at least degree ≥ 1) polynomial over ℚ.
    pub fn adjoin_irreducible_over_q(min_poly_q: CoordsQ) -> Result<Arc<Self>, EvalError> {
        let min_poly_q = trim_leading_zero(min_poly_q);
        let d = poly_degree(&min_poly_q);
        if d < 1 {
            return Err(EvalError::TypeError("extension degree >= 1"));
        }
        field_registry().register_simple(&min_poly_q)
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
        let mut v = self.zero_coords();
        if let Some(c) = v.last_mut() {
            *c = Ratio::one();
        }
        v
    }

    pub(crate) fn generator_coords(&self) -> CoordsQ {
        let mut v = self.zero_coords();
        if !v.is_empty() {
            v[0] = Ratio::one();
        }
        v
    }

    pub fn element_add(&self, a: &CoordsQ, b: &CoordsQ) -> Result<CoordsQ, EvalError> {
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

struct FieldRegistry {
    by_min_poly: Mutex<HashMap<Vec<u8>, Arc<ExtensionField>>>,
    common_cache: Mutex<HashMap<(u64, u64), Arc<CommonFieldPair>>>,
}

impl FieldRegistry {
    fn new() -> Self {
        Self {
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

    fn register_simple(&self, min_poly_q: &CoordsQ) -> Result<Arc<ExtensionField>, EvalError> {
        let key = Self::min_poly_key(min_poly_q);
        if let Some(hit) = self.by_min_poly.lock().unwrap().get(&key) {
            return Ok(Arc::clone(hit));
        }
        let d = poly_degree(min_poly_q);
        let field = Arc::new(ExtensionField {
            id: next_field_id(),
            tower: Arc::new(ExtensionTower::Adj {
                parent: Arc::new(ExtensionTower::Base),
                min_poly_q: min_poly_q.clone(),
                ext_degree: d,
            }),
            min_poly_over_q: min_poly_q.clone(),
        });
        self.by_min_poly
            .lock()
            .unwrap()
            .insert(key, Arc::clone(&field));
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
                }),
                min_poly_over_q: min_poly_q,
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

    fn sqrt2_field() -> Arc<ExtensionField> {
        ExtensionField::adjoin_irreducible_over_q(vec![
            Ratio::one(),
            Ratio::zero(),
            Ratio::from_integer((-2).into()),
        ])
        .unwrap()
    }

    fn cbrt2_field() -> Arc<ExtensionField> {
        ExtensionField::adjoin_irreducible_over_q(vec![
            Ratio::one(),
            Ratio::zero(),
            Ratio::zero(),
            Ratio::from_integer((-2).into()),
        ])
        .unwrap()
    }

    #[test]
    fn rational_field_dimension_one() {
        let q = ExtensionField::rational();
        assert_eq!(q.dimension(), 1);
        assert!(q.element_is_one(&q.one_coords()));
    }

    #[test]
    fn adjoin_sqrt2_dimension_two() {
        let k = sqrt2_field();
        assert_eq!(k.dimension(), 2);
        let alpha = k.generator_coords();
        let sq = k.element_mul(&alpha, &alpha).unwrap();
        let two = pad_to_len(&[Ratio::from_integer(2.into())], 2);
        assert!(k.element_eq_mod(&sq, &two).unwrap());
    }

    #[test]
    #[ignore = "primitive-element common(√2,∛2) char poly is slow; see ext_tower perf follow-up"]
    fn common_sqrt2_cbrt2_has_degree_six() {
        let a = sqrt2_field();
        let b = cbrt2_field();
        let common = ExtensionField::common_over_q(&a, &b).unwrap();
        assert_eq!(common.field.dimension(), 6);
        let ea = common.embed_a.apply(&a.generator_coords());
        let eb = common.embed_b.apply(&b.generator_coords());
        let sum = common.field.element_add(&ea, &eb).unwrap();
        assert!(!common.field.element_is_zero(&sum));
    }

    #[test]
    fn common_cache_identity_is_fast() {
        let a = sqrt2_field();
        let c1 = ExtensionField::common_over_q(&a, &a).unwrap();
        let c2 = ExtensionField::common_over_q(&a, &a).unwrap();
        assert_eq!(c1.field.id(), c2.field.id());
    }

    #[test]
    #[ignore = "primitive-element common(√2,∛2) char poly is slow; see ext_tower perf follow-up"]
    fn common_cache_hits_same_pair() {
        let a = sqrt2_field();
        let b = cbrt2_field();
        let c1 = ExtensionField::common_over_q(&a, &b).unwrap();
        let c2 = ExtensionField::common_over_q(&b, &a).unwrap();
        assert_eq!(c1.field.id(), c2.field.id());
    }

    #[test]
    fn field_registry_dedup_same_minpoly() {
        let k1 = sqrt2_field();
        let k2 = sqrt2_field();
        assert_eq!(k1.id(), k2.id());
    }
}
