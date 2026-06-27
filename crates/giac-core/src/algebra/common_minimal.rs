//! Upstream `common_minimal_POLY` (U1a rational + U1b parent-coeff `mb`).
//!
//! See `.doc/issues/GIAC-ext-common-unified-path.md` and `giac-2.0.0/src/alg_ext.cc`.

use std::cell::Cell;
use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;

use super::{
    eval_rational_poly1_at_field, flatten_min_poly_over_q, gauss_elim_rref,
    get_or_create_base_by_min_poly, mult_matrix_layer_gen, verify_common_pair, CommonFieldPair,
    ExtensionField, FieldEmbedding, LayerMinPolyForAdjoin, PRIMITIVE_K_SEARCH_MAX,
};
use crate::algebra::field_arith::{
    apply_linear_map, generator_coords, identity_matrix, kron_left, mat_mul, mrref_q,
    mult_matrix_of_element, pad_to_len, poly_degree, rationalize_poly1, trim_leading_zero, CoordsQ,
};
use crate::algebra::field_session::{coords_in_field, FieldSession};
use crate::algebra::poly::PolyAlgExt;
use crate::algebra::poly_alg_coeff::AlgExtCPolyCoeff;
use crate::algebra::poly_alg_factor::{
    common_minpoly_var, factor_minpoly_over_field, minpoly_q_to_poly_over_field,
    select_factor_for_common,
};
use crate::algebra::poly_alg_ops::normalize_algext_poly;

thread_local! {
    static COMPOSITUM_SESSION_DEPTH: Cell<u32> = const { Cell::new(0) };
}

/// Active compositum pipeline token (T3). Construct only via [`with_compositum_session`].
pub(crate) struct CompositumSession(());

struct CompositumSessionGuard;

impl CompositumSessionGuard {
    fn enter() -> Self {
        COMPOSITUM_SESSION_DEPTH.with(|d| d.set(d.get() + 1));
        Self
    }
}

impl Drop for CompositumSessionGuard {
    fn drop(&mut self) {
        COMPOSITUM_SESSION_DEPTH.with(|d| d.set(d.get().saturating_sub(1)));
    }
}

/// Run compositum pipeline with an active session (align must not re-enter minimal).
pub(crate) fn with_compositum_session<R>(
    f: impl FnOnce(&CompositumSession) -> Result<R, EvalError>,
) -> Result<R, EvalError> {
    let _guard = CompositumSessionGuard::enter();
    f(&CompositumSession(()))
}

/// True while a [`CompositumSession`] is active (see [`with_compositum_session`]).
pub(crate) fn compositum_session_active() -> bool {
    COMPOSITUM_SESSION_DEPTH.with(|d| d.get() > 0)
}

/// True while [`compute_compositum_upstream`] is on the stack (alias for dispatch).
pub(crate) fn compositum_upstream_active() -> bool {
    compositum_session_active()
}

/// mrref compositum data produced atomically by upstream (T1); do not split fields across calls.
pub(crate) struct MrrefCompositumSpec {
    min_g: CoordsQ,
    k: i64,
    mat_theta: Vec<Vec<Ratio<BigInt>>>,
    w_a: CoordsQ,
    w_b: CoordsQ,
    na: usize,
    nb: usize,
}

impl MrrefCompositumSpec {
    fn new(
        min_g: CoordsQ,
        k: i64,
        mat_theta: Vec<Vec<Ratio<BigInt>>>,
        w_a: CoordsQ,
        w_b: CoordsQ,
        na: usize,
        nb: usize,
    ) -> Result<Self, EvalError> {
        let dim = na.checked_mul(nb).ok_or(EvalError::TypeError(
            "common_minimal: compositum dim overflow",
        ))?;
        if mat_theta.len() != dim || mat_theta.iter().any(|row| row.len() != dim) {
            return Err(EvalError::TypeError(
                "common_minimal: mat_theta dim mismatch",
            ));
        }
        if poly_degree(&min_g) != dim {
            return Err(EvalError::TypeError(
                "common_minimal: min_g degree mismatch",
            ));
        }
        Ok(Self {
            min_g,
            k,
            mat_theta,
            w_a,
            w_b,
            na,
            nb,
        })
    }

    pub(crate) fn min_g(&self) -> &CoordsQ {
        &self.min_g
    }

    pub(crate) fn k(&self) -> i64 {
        self.k
    }

    pub(crate) fn mat_theta(&self) -> &[Vec<Ratio<BigInt>>] {
        &self.mat_theta
    }

    pub(crate) fn w_a(&self) -> &CoordsQ {
        &self.w_a
    }

    pub(crate) fn w_b(&self) -> &CoordsQ {
        &self.w_b
    }

    pub(crate) fn na(&self) -> usize {
        self.na
    }

    pub(crate) fn nb(&self) -> usize {
        self.nb
    }

    pub(crate) fn dim(&self) -> usize {
        self.na * self.nb
    }
}

/// Tensor basis index with `(ja,jb)` = exponents of adjoining gens `a`, `b`.
fn tensor_index_exp(ja: usize, jb: usize, na: usize, nb: usize) -> usize {
    (na - 1 - ja) + (nb - 1 - jb) * na
}

/// Tensor basis index (matches [`field_arith::tower_flat_index`]).
fn tensor_index(i: usize, u_power: usize, na: usize, nb: usize) -> usize {
    (nb - 1 - u_power) * na + i
}

// **Pipeline private** — bivariate reduction mod rational `ma`, `mb` (U1a).
fn add_tensor_term_rational(
    v: &mut [Ratio<BigInt>],
    ja: usize,
    jb: usize,
    c: Ratio<BigInt>,
    ma: &[Ratio<BigInt>],
    mb: &[Ratio<BigInt>],
    na: usize,
    nb: usize,
) {
    if c.is_zero() {
        return;
    }
    if ja >= na {
        for t in 1..=na {
            add_tensor_term_rational(v, ja - t, jb, -c.clone() * ma[t].clone(), ma, mb, na, nb);
        }
        return;
    }
    if jb >= nb {
        for t in 1..=nb {
            add_tensor_term_rational(v, ja, jb - t, -c.clone() * mb[t].clone(), ma, mb, na, nb);
        }
        return;
    }
    v[tensor_index_exp(ja, jb, na, nb)] += c;
}

fn mul_by_gen_a_rational(
    v: &[Ratio<BigInt>],
    ma: &[Ratio<BigInt>],
    mb: &[Ratio<BigInt>],
    na: usize,
    nb: usize,
) -> Vec<Ratio<BigInt>> {
    let dim = na * nb;
    let mut out = vec![Ratio::zero(); dim];
    for jb in 0..nb {
        for ja in 0..na {
            let c = v[tensor_index_exp(ja, jb, na, nb)].clone();
            if !c.is_zero() {
                add_tensor_term_rational(&mut out, ja + 1, jb, c, ma, mb, na, nb);
            }
        }
    }
    out
}

fn mul_by_gen_b_rational(
    v: &[Ratio<BigInt>],
    ma: &[Ratio<BigInt>],
    mb: &[Ratio<BigInt>],
    na: usize,
    nb: usize,
) -> Vec<Ratio<BigInt>> {
    let dim = na * nb;
    let mut out = vec![Ratio::zero(); dim];
    for jb in 0..nb {
        for ja in 0..na {
            let c = v[tensor_index_exp(ja, jb, na, nb)].clone();
            if !c.is_zero() {
                add_tensor_term_rational(&mut out, ja, jb + 1, c, ma, mb, na, nb);
            }
        }
    }
    out
}

fn mul_by_theta_rational(
    v: &[Ratio<BigInt>],
    k: i64,
    ma: &[Ratio<BigInt>],
    mb: &[Ratio<BigInt>],
    na: usize,
    nb: usize,
) -> Vec<Ratio<BigInt>> {
    let av = mul_by_gen_a_rational(v, ma, mb, na, nb);
    let bv = mul_by_gen_b_rational(v, ma, mb, na, nb);
    if k == 0 {
        return bv;
    }
    let k_rat = Ratio::from_integer(BigInt::from(k));
    av.into_iter()
        .zip(bv)
        .map(|(a, b)| a + k_rat.clone() * b)
        .collect()
}

fn build_mrref_matrix_rational(
    ma: &[Ratio<BigInt>],
    mb: &[Ratio<BigInt>],
    na: usize,
    nb: usize,
    k: i64,
) -> Vec<Vec<Ratio<BigInt>>> {
    let dim = na * nb;
    let mut rows = Vec::with_capacity(dim + 3);
    let mut p = vec![Ratio::zero(); dim];
    p[tensor_index_exp(0, 0, na, nb)] = Ratio::one();
    rows.push(p.clone());
    for _ in 1..=dim {
        p = mul_by_theta_rational(&p, k, ma, mb, na, nb);
        rows.push(p.clone());
    }
    let mut gen_b = vec![Ratio::zero(); dim];
    gen_b[tensor_index_exp(0, 1, na, nb)] = Ratio::one();
    rows.push(gen_b);
    let mut gen_a = vec![Ratio::zero(); dim];
    gen_a[tensor_index_exp(1, 0, na, nb)] = Ratio::one();
    rows.push(gen_a);
    let nrows = dim;
    let ncols = dim + 3;
    let mut m = vec![vec![Ratio::zero(); ncols]; nrows];
    for (i, row) in rows.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if j < nrows {
                m[j][i] = c.clone();
            }
        }
    }
    m
}

enum MbSpec<'a> {
    Rational(&'a [Ratio<BigInt>]),
    ParentBlocks {
        parent: &'a ExtensionField,
        blocks: &'a [CoordsQ],
    },
}

fn mat_b_from_spec(spec: &MbSpec<'_>, na: usize) -> Result<(usize, Vec<Vec<Ratio<BigInt>>>), EvalError> {
    match spec {
        MbSpec::Rational(mb) => {
            let nb = poly_degree(mb);
            Ok((
                nb,
                kron_left(
                    &mult_matrix_of_element(&generator_coords(nb), mb),
                    na,
                ),
            ))
        }
        MbSpec::ParentBlocks { parent, blocks } => {
            let nb = blocks.len().checked_sub(1).ok_or(EvalError::TypeError(
                "common_minimal: layer minpoly",
            ))?;
            Ok((nb, mult_matrix_layer_gen(parent, blocks)?))
        }
    }
}

fn build_mat_theta(
    ma: &[Ratio<BigInt>],
    na: usize,
    mb: &MbSpec<'_>,
    k: i64,
) -> Result<(usize, Vec<Vec<Ratio<BigInt>>>), EvalError> {
    let (nb, mat_b) = mat_b_from_spec(mb, na)?;
    let mat_a = kron_left(&mult_matrix_of_element(&generator_coords(na), ma), nb);
    let k_rat = Ratio::from_integer(BigInt::from(k));
    let dim = na * nb;
    let mut mat_theta = vec![vec![Ratio::zero(); dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            mat_theta[i][j] = mat_a[i][j].clone() + k_rat.clone() * mat_b[i][j].clone();
        }
    }
    Ok((nb, mat_theta))
}

fn tensor_unit_1(na: usize, nb: usize) -> Vec<Ratio<BigInt>> {
    let dim = na * nb;
    let mut v = vec![Ratio::zero(); dim];
    v[tensor_index(na - 1, 0, na, nb)] = Ratio::one();
    v
}

fn tensor_gen_a(na: usize, nb: usize) -> Vec<Ratio<BigInt>> {
    let dim = na * nb;
    let mut v = vec![Ratio::zero(); dim];
    if na >= 2 {
        v[tensor_index(na - 2, 0, na, nb)] = Ratio::one();
    }
    v
}

fn tensor_gen_b(na: usize, nb: usize) -> Vec<Ratio<BigInt>> {
    let dim = na * nb;
    let mut v = vec![Ratio::zero(); dim];
    if nb >= 2 {
        v[tensor_index(na - 1, 1, na, nb)] = Ratio::one();
    }
    v
}

fn build_mrref_matrix(
    mat_theta: &[Vec<Ratio<BigInt>>],
    na: usize,
    nb: usize,
) -> Vec<Vec<Ratio<BigInt>>> {
    let dim = na * nb;
    let mut rows = Vec::with_capacity(dim + 3);
    let mut p = tensor_unit_1(na, nb);
    rows.push(p.clone());
    for _ in 1..=dim {
        p = apply_linear_map(mat_theta, &p);
        rows.push(p.clone());
    }
    rows.push(tensor_gen_b(na, nb));
    rows.push(tensor_gen_a(na, nb));
    let nrows = dim;
    let ncols = dim + 3;
    let mut m = vec![vec![Ratio::zero(); ncols]; nrows];
    for (i, row) in rows.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if j < nrows {
                m[j][i] = c.clone();
            }
        }
    }
    m
}

fn mrref_candidate_ok(m: &[Vec<Ratio<BigInt>>], dim: usize) -> bool {
    let row = dim - 1;
    m[row][0..dim].iter().any(|c| !c.is_zero())
}

fn minpoly_from_rref(m: &[Vec<Ratio<BigInt>>], dim: usize) -> CoordsQ {
    let mut v: Vec<Ratio<BigInt>> = (0..dim).map(|i| -m[i][dim].clone()).collect();
    v.push(Ratio::one());
    v.reverse();
    trim_leading_zero(v)
}

fn embed_poly_from_rref_col(m: &[Vec<Ratio<BigInt>>], dim: usize, col: usize) -> CoordsQ {
    let mut w: Vec<Ratio<BigInt>> = (0..dim).map(|i| m[i][col].clone()).collect();
    w.reverse();
    trim_leading_zero(w)
}

fn common_minimal_poly_core(
    ma: &[Ratio<BigInt>],
    na: usize,
    mb: MbSpec<'_>,
    k_start: i64,
) -> Result<(CoordsQ, i64, CoordsQ, CoordsQ), EvalError> {
    let nb = match &mb {
        MbSpec::Rational(m) => poly_degree(m),
        MbSpec::ParentBlocks { blocks, .. } => blocks.len().checked_sub(1).ok_or(
            EvalError::TypeError("common_minimal: layer minpoly"),
        )?,
    };
    if na == 0 || nb == 0 {
        return Err(EvalError::TypeError("common_minimal: degree zero"));
    }
    if nb == 1 {
        if let MbSpec::Rational(mb) = mb {
            let k = 0;
            let mut w_a = vec![Ratio::zero(); na];
            if na > 0 {
                w_a[0] = Ratio::one();
            }
            let c = -mb[nb].clone();
            let w_b = if c.is_zero() { Vec::new() } else { vec![c] };
            return Ok((ma.to_vec(), k, w_a, trim_leading_zero(w_b)));
        }
    }
    let dim = na * nb;
    for k in k_start..=PRIMITIVE_K_SEARCH_MAX {
        let mut m = match &mb {
            MbSpec::Rational(mb) => build_mrref_matrix_rational(ma, mb, na, nb, k),
            MbSpec::ParentBlocks { .. } => {
                let (_, mat_theta) = build_mat_theta(ma, na, &mb, k)?;
                build_mrref_matrix(&mat_theta, na, nb)
            }
        };
        let rank = mrref_q(&mut m);
        if rank < dim || !mrref_candidate_ok(&m, dim) {
            continue;
        }
        let min_g = minpoly_from_rref(&m, dim);
        if poly_degree(&min_g) != dim {
            continue;
        }
        let w_b = embed_poly_from_rref_col(&m, dim, dim + 1);
        let w_a = embed_poly_from_rref_col(&m, dim, dim + 2);
        return Ok((min_g, k, w_a, w_b));
    }
    Err(EvalError::NotImplemented("common_minimal_poly: no k found"))
}

/// U1a — both minpolys over ℚ.
// **Pipeline private** — `common_minimal_poly_over_q`
pub(crate) fn common_minimal_poly_over_q(
    ma: &[Ratio<BigInt>],
    mb: &[Ratio<BigInt>],
) -> Result<(CoordsQ, i64, CoordsQ, CoordsQ), EvalError> {
    let na = poly_degree(ma);
    common_minimal_poly_core(ma, na, MbSpec::Rational(mb), 1)
}

/// U1b — `mb` with parent-field coefficients over simple `parent` (minpoly `ma`).
// **Pipeline private** — `common_minimal_poly_over_parent`
pub(crate) fn common_minimal_poly_over_parent(
    parent: &ExtensionField,
    ma: &[Ratio<BigInt>],
    blocks: &[CoordsQ],
) -> Result<(CoordsQ, i64, CoordsQ, CoordsQ), EvalError> {
    if !parent.tower().is_simple_over_q() {
        return Err(EvalError::TypeError(
            "common_minimal_over_parent: parent not simple over Q",
        ));
    }
    let na = parent.dimension();
    if poly_degree(ma) != na {
        return Err(EvalError::TypeError(
            "common_minimal_over_parent: ma degree mismatch",
        ));
    }
    // ponytail: upstream k_init=1 when mb coeffs reference parent ext
    common_minimal_poly_core(
        ma,
        na,
        MbSpec::ParentBlocks { parent, blocks },
        1,
    )
}

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

// **Pipeline private** — invert square matrix over ℚ (dim ≤ compositum bound).
fn invert_rational_matrix(m: &[Vec<Ratio<BigInt>>]) -> Result<Vec<Vec<Ratio<BigInt>>>, EvalError> {
    let n = m.len();
    if n == 0 || m.iter().any(|row| row.len() != n) {
        return Err(EvalError::TypeError("invert: not square"));
    }
    let mut aug = vec![vec![Ratio::zero(); 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = m[i][j].clone();
        }
        aug[i][n + i] = Ratio::one();
    }
    let (_, inconsistent) = gauss_elim_rref(&mut aug);
    if inconsistent {
        return Err(EvalError::TypeError("invert: singular"));
    }
    Ok((0..n)
        .map(|i| (0..n).map(|j| aug[i][n + j].clone()).collect())
        .collect())
}

/// True when `{1, prim, …, prim^{d-1}}` spans **F** operationally (ponytail: dim≤8).
fn power_basis_invertible(field: &ExtensionField, prim: &CoordsQ) -> Result<bool, EvalError> {
    let dim = field.dimension();
    let mut power_to_ops = vec![vec![Ratio::zero(); dim]; dim];
    let mut p = field.one_coords();
    for i in 0..dim {
        for j in 0..dim {
            power_to_ops[j][i] = p[j].clone();
        }
        if i + 1 < dim {
            p = field.element_mul(&p, prim)?;
        }
    }
    invert_rational_matrix(&power_to_ops).map(|_| true).or(Ok(false))
}

/// Operational coords for flat compositum primitive `parent_gen ± k·layer_gen`.
fn nested_primitive_from_k(field: &ExtensionField, k: i64) -> Result<CoordsQ, EvalError> {
    let parent = field
        .parent_field()
        .filter(|p| !p.is_base())
        .ok_or(EvalError::TypeError("common_minimal: nested primitive needs parent"))?;
    let pd = parent.dimension();
    let mut pg = field.zero_coords();
    let pgen = parent.generator_coords();
    for i in 0..pd {
        if i < pgen.len() {
            pg[i] = pgen[i].clone();
        }
    }
    let g = field.generator_coords();
    let krat = field.embed_rational(&Ratio::from_integer(BigInt::from(k)));
    let kg = field.element_mul(&krat, &g)?;
    field.element_add(&pg, &kg)
}

/// Columns of **P** are operational coords of `prim^i`; returns **P**⁻¹.
fn ops_power_change_inverse(
    field: &ExtensionField,
    prim: &CoordsQ,
) -> Result<Vec<Vec<Ratio<BigInt>>>, EvalError> {
    let dim = field.dimension();
    let mut power_to_ops = vec![vec![Ratio::zero(); dim]; dim];
    let mut p = field.one_coords();
    for i in 0..dim {
        for j in 0..dim {
            power_to_ops[j][i] = p[j].clone();
        }
        if i + 1 < dim {
            p = field.element_mul(&p, prim)?;
        }
    }
    invert_rational_matrix(&power_to_ops)
}

/// Embed tensor raw vector into compositum θ-basis (matches [`ext_tower::embed_in_gamma_vector`]).
fn embed_in_gamma_vector(
    v: &[Ratio<BigInt>],
    mat_theta: &[Vec<Ratio<BigInt>>],
    dim: usize,
) -> Result<CoordsQ, EvalError> {
    let mut pow = identity_matrix(dim);
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

/// Nested **F** ↪ common via mrref block model (β^j at `j·na` in compositum tensor).
fn nested_embedding_from_compositum(
    source: &Arc<ExtensionField>,
    target: &Arc<ExtensionField>,
    spec: &MrrefCompositumSpec,
) -> Result<FieldEmbedding, EvalError> {
    let nb = source.dimension();
    let dim = target.dimension();
    let partner_na = spec.na();
    let mat_theta = spec.mat_theta();
    let compositum_k = Some(spec.k());
    let mut beta_images = Vec::with_capacity(nb);
    for j in 0..nb {
        let mut raw = vec![Ratio::zero(); dim];
        if j * partner_na < raw.len() {
            raw[j * partner_na] = Ratio::one();
        }
        beta_images.push(embed_in_gamma_vector(&raw, mat_theta, dim)?);
    }
    let candidates = flat_primitive_candidates(source.as_ref(), compositum_k)?;
    for prim in candidates {
        let inv = ops_power_change_inverse(source.as_ref(), &prim)?;
        let mut matrix = vec![vec![Ratio::zero(); nb]; dim];
        for k in 0..nb {
            for i in 0..nb {
                let coeff = inv[i][k].clone();
                if coeff.is_zero() {
                    continue;
                }
                for tgt_i in 0..dim {
                    matrix[tgt_i][k] += coeff.clone() * beta_images[i][tgt_i].clone();
                }
            }
        }
        let emb = FieldEmbedding {
            source: Arc::clone(source),
            target: Arc::clone(target),
            matrix,
        };
        let one = emb.apply(&source.one_coords());
        if target
            .element_eq_mod(&one, &target.one_coords())
            .unwrap_or(false)
        {
            return Ok(emb);
        }
    }
    Err(EvalError::NotImplemented("common_minimal: nested embed"))
}

/// Candidate primitives for nested ops↔power change (ponytail: dim≤8 scan).
fn flat_primitive_candidates(
    field: &ExtensionField,
    k_prefer: Option<i64>,
) -> Result<Vec<CoordsQ>, EvalError> {
    let dim = field.dimension();
    let mut out = Vec::new();
    let mut push_if = |v: &CoordsQ| -> Result<(), EvalError> {
        if field.element_is_zero(v) || !power_basis_invertible(field, v)? {
            return Ok(());
        }
        if !out.iter().any(|u| field.element_eq_mod(u, v).unwrap_or(false)) {
            out.push(v.clone());
        }
        Ok(())
    };
    push_if(&field.generator_coords())?;
    if let Some(k) = k_prefer {
        for kk in [k, -k] {
            if kk != 0 {
                if let Ok(c) = nested_primitive_from_k(field, kk) {
                    push_if(&c)?;
                }
            }
        }
    }
    let mut vecs = vec![field.generator_coords()];
    for i in 0..dim {
        let mut e = field.zero_coords();
        if i < e.len() {
            e[i] = Ratio::one();
        }
        vecs.push(e);
    }
    if field.parent_field().is_some_and(|p| !p.is_base()) {
        for k in 1i64..=PRIMITIVE_K_SEARCH_MAX {
            for kk in [k, -k] {
                if let Ok(c) = nested_primitive_from_k(field, kk) {
                    push_if(&c)?;
                }
            }
        }
    }
    for v in &vecs {
        push_if(v)?;
    }
    for &ci in &[-2i32, -1, 1, 2] {
        for &cj in &[-2i32, -1, 1, 2] {
            if ci == 0 && cj == 0 {
                continue;
            }
            for a in 0..vecs.len() {
                for b in 0..vecs.len() {
                    let ca = field.embed_rational(&Ratio::from_integer(BigInt::from(ci)));
                    let cb = field.embed_rational(&Ratio::from_integer(BigInt::from(cj)));
                    let mut v = field.zero_coords();
                    if ci != 0 {
                        v = field.element_add(&v, &field.element_mul(&ca, &vecs[a])?)?;
                    }
                    if cj != 0 {
                        v = field.element_add(&v, &field.element_mul(&cb, &vecs[b])?)?;
                    }
                    push_if(&v)?;
                }
            }
        }
    }
    if out.is_empty() {
        Err(EvalError::NotImplemented("common_minimal: flat primitive"))
    } else {
        Ok(out)
    }
}

/// Some θ ∈ **F** spanning `{θ^i}` operationally (ponytail: dim≤8 combo scan).
fn find_flat_primitive_coords(
    field: &ExtensionField,
    k_prefer: Option<i64>,
) -> Result<CoordsQ, EvalError> {
    flat_primitive_candidates(field, k_prefer)?
        .into_iter()
        .next()
        .ok_or(EvalError::NotImplemented("common_minimal: flat primitive"))
}

fn embedding_from_theta_poly(
    source: &Arc<ExtensionField>,
    target: &Arc<ExtensionField>,
    poly: &CoordsQ,
    compositum_k: Option<i64>,
    mrref_spec: Option<&MrrefCompositumSpec>,
) -> Result<FieldEmbedding, EvalError> {
    if let Some(spec) = mrref_spec {
        return nested_embedding_from_compositum(source, target, spec);
    }
    let src_dim = source.dimension();
    let tgt_dim = target.dimension();
    let gen = target.generator_coords();
    let w = eval_rational_poly1_at_field(target.as_ref(), poly, &gen)?;
    let mut matrix = vec![vec![Ratio::zero(); src_dim]; tgt_dim];
    if source.tower().is_simple_over_q() {
        for j in 0..src_dim {
            let exp = (src_dim - 1 - j) as u64;
            let image = element_pow_coords(target.as_ref(), &w, exp)?;
            for i in 0..tgt_dim {
                matrix[i][j] = image[i].clone();
            }
        }
    } else {
        let prim = find_flat_primitive_coords(source.as_ref(), compositum_k)?;
        let inv = ops_power_change_inverse(source.as_ref(), &prim)?;
        for k in 0..src_dim {
            let mut w_pow = target.one_coords();
            for i in 0..src_dim {
                let coeff = inv[i][k].clone();
                if !coeff.is_zero() {
                    for tgt_i in 0..tgt_dim {
                        matrix[tgt_i][k] += coeff.clone() * w_pow[tgt_i].clone();
                    }
                }
                if i + 1 < src_dim {
                    w_pow = target.element_mul(&w_pow, &w)?;
                }
            }
        }
    }
    Ok(FieldEmbedding {
        source: Arc::clone(source),
        target: Arc::clone(target),
        matrix,
    })
}

fn simple_layer_minpoly(field: &ExtensionField) -> Result<CoordsQ, EvalError> {
    match super::layer_minpoly_coords_for_adjoin(field)? {
        LayerMinPolyForAdjoin::Rational(m) => Ok(m.to_vec()),
        LayerMinPolyForAdjoin::ParentBlocks(_) => Err(EvalError::TypeError(
            "common_minimal: expected rational layer",
        )),
    }
}

/// Layer minpoly for simple fields; explicit ℚ-flatten otherwise (U2 input).
fn minpoly_for_common(
    field: &ExtensionField,
    session: Option<&FieldSession>,
) -> Result<CoordsQ, EvalError> {
    if field.tower().is_simple_over_q() {
        simple_layer_minpoly(field)
    } else {
        flatten_min_poly_over_q(field, session)
    }
}

/// Degree proxy for operand ordering (avoid flatten in `order_operands_for_common`).
fn common_operand_degree(field: &ExtensionField) -> Result<usize, EvalError> {
    if field.tower().is_simple_over_q() {
        Ok(poly_degree(&simple_layer_minpoly(field)?))
    } else {
        Ok(field.dimension())
    }
}

fn order_operands_for_common(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    _session: Option<&FieldSession>,
) -> Result<(Arc<ExtensionField>, Arc<ExtensionField>, bool), EvalError> {
    let da = common_operand_degree(a)?;
    let db = common_operand_degree(b)?;
    let swap = db < da || (db == da && b.semantic_key() < a.semantic_key());
    Ok(if swap {
        (Arc::clone(b), Arc::clone(a), true)
    } else {
        (Arc::clone(a), Arc::clone(b), false)
    })
}

fn is_embedded_rational(base: &ExtensionField, coords: &CoordsQ) -> bool {
    let coords = pad_to_len(coords, base.dimension());
    if base.dimension() == 1 {
        return true;
    }
    coords
        .iter()
        .take(base.dimension().saturating_sub(1))
        .all(|c| c.is_zero())
}

fn coeff_as_embedded_rational(
    base: &ExtensionField,
    c: &AlgExtCPolyCoeff,
) -> Result<Ratio<BigInt>, EvalError> {
    let inner = c.as_inner();
    if !inner.im.iter().all(|e| e.is_zero()) {
        return Err(EvalError::TypeError("expected real coefficient"));
    }
    let re = rationalize_poly1(&inner.re)?;
    if inner.field.is_base() {
        return Ok(pad_to_len(&re, 1)[0].clone());
    }
    let coords = pad_to_len(&re, base.dimension());
    if !is_embedded_rational(base, &coords) {
        return Err(EvalError::TypeError(
            "common_minimal: factor coeff not rational over Q",
        ));
    }
    Ok(coords[base.dimension().saturating_sub(1)].clone())
}

fn factor_to_minpoly_q(
    f: &PolyAlgExt,
    base: &Arc<ExtensionField>,
    _session: &FieldSession,
) -> Result<CoordsQ, EvalError> {
    let var = common_minpoly_var();
    let deg = f.degree_wrt(&var) as usize;
    let mut out = vec![Ratio::zero(); deg + 1];
    for (slot, exp) in (0..=deg).rev().enumerate() {
        let c = giac_poly::scalar_coeff_wrt(f, &var, exp as u64);
        out[slot] = coeff_as_embedded_rational(base.as_ref(), &c)?;
    }
    Ok(trim_leading_zero(out))
}

fn factor_to_parent_blocks(
    f: &PolyAlgExt,
    base: &Arc<ExtensionField>,
    session: &FieldSession,
) -> Result<Vec<CoordsQ>, EvalError> {
    let var = common_minpoly_var();
    let deg = f.degree_wrt(&var) as usize;
    let mut blocks = Vec::with_capacity(deg + 1);
    for exp in (0..=deg).rev() {
        let c = giac_poly::scalar_coeff_wrt(f, &var, exp as u64);
        blocks.push(coords_in_field(session, &c, base)?);
    }
    Ok(blocks)
}

/// Upstream `common_EXT` compositum: factor `b` over **K**\_`a`, then U1 minimal + mrref (U2).
fn compute_compositum_upstream(
    _session: &CompositumSession,
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    field_session: Option<&FieldSession>,
) -> Result<MrrefCompositumSpec, EvalError> {
    compute_compositum_upstream_inner(a, b, field_session)
}

fn compute_compositum_upstream_inner(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    session: Option<&FieldSession>,
) -> Result<MrrefCompositumSpec, EvalError> {
    let (fa, fb, swapped) = order_operands_for_common(a, b, session)?;
    let ma = minpoly_for_common(&fa, session)?;
    let mb_full = minpoly_for_common(&fb, session)?;
    let var = common_minpoly_var();
    let coeff_session = FieldSession::new(Arc::clone(&fa));
    let selected = if fb.tower().is_simple_over_q() {
        let factors = factor_minpoly_over_field(&mb_full, &fa)?;
        select_factor_for_common(&factors, &fb, &fa)?
    } else {
        // ponytail: nested b — Q-factor split without vanishing probe is unreliable; use full mb
        let p_alg = minpoly_q_to_poly_over_field(&mb_full, &fa, &var)?;
        normalize_algext_poly(&p_alg, &coeff_session)?
    };
    let na = poly_degree(&ma);
    let (min_g, k, w_fa, w_fb, mat_theta, nb) = if fa.tower().is_simple_over_q() {
        if let Ok(mb_q) = factor_to_minpoly_q(&selected, &fa, &coeff_session) {
            let (min_g, k, w_a, w_b) = common_minimal_poly_over_q(&ma, &mb_q)?;
            let (nb, mat) = build_mat_theta(&ma, na, &MbSpec::Rational(&mb_q), k)?;
            (min_g, k, w_a, w_b, mat, nb)
        } else {
            let blocks = factor_to_parent_blocks(&selected, &fa, &coeff_session)?;
            let (min_g, k, w_a, w_b) =
                common_minimal_poly_over_parent(fa.as_ref(), &ma, &blocks)?;
            let (nb, mat) = build_mat_theta(
                &ma,
                na,
                &MbSpec::ParentBlocks {
                    parent: fa.as_ref(),
                    blocks: &blocks,
                },
                k,
            )?;
            (min_g, k, w_a, w_b, mat, nb)
        }
    } else {
        let mb_q = factor_to_minpoly_q(&selected, &fa, &coeff_session)?;
        let (min_g, k, w_a, w_b) = common_minimal_poly_over_q(&ma, &mb_q)?;
        let (nb, mat) = build_mat_theta(&ma, na, &MbSpec::Rational(&mb_q), k)?;
        (min_g, k, w_a, w_b, mat, nb)
    };
    let (w_a, w_b) = if swapped { (w_fb, w_fa) } else { (w_fa, w_fb) };
    MrrefCompositumSpec::new(min_g, k, mat_theta, w_a, w_b, na, nb)
}

#[cfg(test)]
pub(crate) fn embedding_from_theta_poly_for_verify_test(
    source: &Arc<ExtensionField>,
    target: &Arc<ExtensionField>,
    poly: &CoordsQ,
    compositum_k: Option<i64>,
    mrref_spec: Option<&MrrefCompositumSpec>,
) -> Result<FieldEmbedding, EvalError> {
    embedding_from_theta_poly(source, target, poly, compositum_k, mrref_spec)
}

/// Compositum via upstream `common_EXT` + minimal-polynomial + mrref embeddings.
// **Pipeline private** — `compute_common_minimal_pair`
pub(crate) fn compute_common_minimal_pair(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    session: Option<&FieldSession>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    with_compositum_session(|compositum| compute_common_minimal_pair_inner(compositum, a, b, session))
}

fn compute_common_minimal_pair_inner(
    compositum: &CompositumSession,
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
    session: Option<&FieldSession>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    if a.is_base() || b.is_base() {
        return Err(EvalError::TypeError("common_minimal: need two extensions"));
    }
    if ExtensionField::is_subfield_of(a, b) || ExtensionField::is_subfield_of(b, a) {
        return Err(EvalError::TypeError("common_minimal: subfield — not compositum"));
    }

    let spec = compute_compositum_upstream(compositum, a, b, session)?;

    let common = match session {
        Some(s) => s.get_or_create_base_by_min_poly(spec.min_g().clone()),
        None => get_or_create_base_by_min_poly(spec.min_g().clone()),
    };
    let k_hint = Some(spec.k());
    let embed_a = embedding_from_theta_poly(
        a,
        &common,
        spec.w_a(),
        if a.tower().is_simple_over_q() {
            None
        } else {
            k_hint
        },
        if a.tower().is_simple_over_q() {
            None
        } else {
            Some(&spec)
        },
    );
    let embed_b = embedding_from_theta_poly(
        b,
        &common,
        spec.w_b(),
        if b.tower().is_simple_over_q() {
            None
        } else {
            k_hint
        },
        if b.tower().is_simple_over_q() {
            None
        } else {
            Some(&spec)
        },
    );
    if let (Ok(embed_a), Ok(embed_b)) = (embed_a, embed_b) {
        let pair = Arc::new(CommonFieldPair {
            field: Arc::clone(&common),
            embed_a,
            embed_b,
        });
        if verify_common_pair(&pair, a, b).is_ok() {
            return Ok(pair);
        }
    }
    // ponytail: tower operational basis ≠ flat mrref tensor; adjoin when one operand is simple-over-ℚ
    if let Ok(pair) = super::try_common_adjoin_one_simple(a, b) {
        verify_common_pair(&pair, a, b)?;
        return Ok(pair);
    }
    Err(EvalError::NotImplemented("common_minimal: compositum embed"))
}

#[cfg(test)]
pub(crate) fn compute_compositum_upstream_for_test(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
) -> Result<MrrefCompositumSpec, EvalError> {
    with_compositum_session(|session| compute_compositum_upstream(session, a, b, None))
}

#[cfg(test)]
pub(crate) fn compute_common_minimal_pair_for_test(
    a: &Arc<ExtensionField>,
    b: &Arc<ExtensionField>,
) -> Result<Arc<CommonFieldPair>, EvalError> {
    compute_common_minimal_pair(a, b, None)
}

#[cfg(test)]
pub(crate) fn common_minimal_poly_over_parent_for_test(
    parent: &ExtensionField,
    ma: &[Ratio<BigInt>],
    blocks: &[CoordsQ],
) -> Result<(CoordsQ, i64, CoordsQ, CoordsQ), EvalError> {
    common_minimal_poly_over_parent(parent, ma, blocks)
}
