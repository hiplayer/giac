//! Univariate polynomial arithmetic over ℚ in giac `poly1` convention (high degree first).
//!
//! Shared by [`super::ext_tower::ExtensionField`] and [`super::alg_ext::AlgExtData`].
//! Dense ℚ ring ops delegate to [`giac_poly::dense::poly1`] (see GIAC-dense-poly1-refactor D2).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::sync::Arc;

use giac_poly::dense::{self, Poly1Order, RatioRingCtx};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc, FuncKind};

/// Coefficient vector in giac `poly1` order (leading term first).
pub type CoordsQ = Vec<Ratio<BigInt>>;
pub type CoordsQSlice<'a> = &'a [Ratio<BigInt>];

const POLY1_Q: Poly1Order = Poly1Order::HighFirst;
const RATIO_CTX: RatioRingCtx = RatioRingCtx;

// **Pipeline private** — `dense_q`
fn dense_q<T>(r: giac_poly::PolyResult<T>) -> T {
    r.expect("dense poly1 over Q is infallible except inv_mod")
}

/// **Stable** — `trim_leading_zero`
pub fn trim_leading_zero(mut v: CoordsQ) -> CoordsQ {
    dense::trim(&RATIO_CTX, &mut v, POLY1_Q);
    v
}

/// **Stable** — `poly_degree`
pub fn poly_degree(p: &[Ratio<BigInt>]) -> usize {
    dense::poly_degree(p)
}

/// **Stable** — `poly_add`
pub fn poly_add(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::add(&RATIO_CTX, a, b, POLY1_Q))
}

/// **Stable** — `poly_sub`
pub fn poly_sub(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::sub(&RATIO_CTX, a, b, POLY1_Q))
}

/// **Stable** — `poly_mul`
pub fn poly_mul(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::mul(&RATIO_CTX, a, b, POLY1_Q))
}

/// **Stable** — `poly_neg`
pub fn poly_neg(a: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::neg(&RATIO_CTX, a, POLY1_Q))
}

/// **Stable** — `poly_scale`
pub fn poly_scale(a: &[Ratio<BigInt>], s: &Ratio<BigInt>) -> CoordsQ {
    dense_q(dense::scale(&RATIO_CTX, a, s, POLY1_Q))
}

/// Reduce `p` modulo monic `m` (high-degree-first, leading coeff of `m` is ±1).
/// **Stable** — `poly_reduce`
pub fn poly_reduce(p: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::reduce_mod_monic(&RATIO_CTX, p, m, POLY1_Q))
}

/// **Stable** — `poly_inv_mod`
pub fn poly_inv_mod(a: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> Result<CoordsQ, EvalError> {
    dense::inv_mod(&RATIO_CTX, a, m, POLY1_Q)
}

/// Extended GCD: returns `(g, s)` with `s*a + t*b = g` (only `s` needed for inverse).
/// **Stable** — `poly_ext_gcd`
pub fn poly_ext_gcd(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (CoordsQ, CoordsQ) {
    dense_q(dense::ext_gcd(&RATIO_CTX, a, b, POLY1_Q))
}

/// **Stable** — `poly_divrem`
pub fn poly_divrem(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (CoordsQ, CoordsQ) {
    dense_q(dense::div_rem(&RATIO_CTX, a, b, POLY1_Q))
}

/// **Stable** — `pad_to_len`
pub fn pad_to_len(v: &[Ratio<BigInt>], n: usize) -> CoordsQ {
    if v.len() >= n {
        return v[v.len() - n..].to_vec();
    }
    let mut out = vec![Ratio::zero(); n - v.len()];
    out.extend(v.iter().cloned());
    out
}

/// **Stable** — `generator_coords`
pub fn generator_coords(n: usize) -> CoordsQ {
    let mut v = vec![Ratio::zero(); n];
    if n > 0 {
        v[0] = Ratio::one();
    }
    v
}

/// **Stable** — `expr_to_ratio`
pub fn expr_to_ratio(e: &Expr) -> Result<Ratio<BigInt>, EvalError> {
    match e {
        Expr::Int(n) => Ok(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Ok(r.clone()),
        _ => Err(EvalError::TypeError("rational coeff expected")),
    }
}

/// **Stable** — `rationalize_poly1`
pub fn rationalize_poly1(items: &[ExprArc]) -> Result<CoordsQ, EvalError> {
    items.iter().map(|e| expr_to_ratio(e.as_ref())).collect()
}

/// **Stable** — `coords_to_expr`
pub fn coords_to_expr(coords: &[Ratio<BigInt>]) -> Result<Vec<ExprArc>, EvalError> {
    Ok(coords
        .iter()
        .map(|r| {
            if r.denom() == &BigInt::one() {
                Expr::int(r.numer().to_string().parse().unwrap_or(0))
            } else {
                Expr::rat(
                    r.numer().to_string().parse().unwrap_or(0),
                    r.denom().to_string().parse().unwrap_or(1),
                )
            }
        })
        .collect())
}

/// **Stable** — `ratio_to_expr_arc`
pub fn ratio_to_expr_arc(r: &Ratio<BigInt>) -> ExprArc {
    crate::ratio_to_expr(r)
}

/// **Stable** — `min_poly_exprs_to_q`
pub fn min_poly_exprs_to_q(min_poly: &[ExprArc]) -> Result<CoordsQ, EvalError> {
    rationalize_poly1(min_poly)
}

/// **Stable** — `poly1_coeffs`
pub fn poly1_coeffs(e: &ExprArc) -> Result<Vec<ExprArc>, EvalError> {
    match e.as_ref() {
        Expr::Func(FuncKind::Poly1, args) => match args.first().map(|a| a.as_ref()) {
            Some(Expr::Seq(items)) => Ok(items.clone()),
            Some(other) => Ok(vec![Arc::new(other.clone())]),
            None => Err(EvalError::TypeError("poly1")),
        },
        Expr::Seq(items) => Ok(items.clone()),
        _ => Err(EvalError::TypeError("poly1 expected")),
    }
}

/// Trim leading zero coefficients while keeping at least one term.
/// **Stable** — `canonical_poly1_expr`
pub fn canonical_poly1_expr(coeffs: &[ExprArc]) -> Vec<ExprArc> {
    if coeffs.is_empty() {
        return vec![Expr::int(0)];
    }
    let start = coeffs
        .iter()
        .position(|c| !c.is_zero())
        .unwrap_or(coeffs.len() - 1);
    coeffs[start..].to_vec()
}

/// **Stable** — `minpoly_at_square`
pub fn minpoly_at_square(m: &[Ratio<BigInt>]) -> CoordsQ {
    let n = poly_degree(m);
    let mut out = vec![Ratio::zero(); 2 * n + 1];
    for i in 0..=n {
        out[2 * i] += m[i].clone();
    }
    trim_leading_zero(out)
}

/// **Stable** — `minpoly_at_cube`
pub fn minpoly_at_cube(m: &[Ratio<BigInt>]) -> CoordsQ {
    let n = poly_degree(m);
    let mut out = vec![Ratio::zero(); 3 * n + 1];
    for i in 0..=n {
        out[3 * i] += m[i].clone();
    }
    trim_leading_zero(out)
}

/// **Stable** — `embed_in_cube_extension`
pub fn embed_in_cube_extension(
    f: &[Ratio<BigInt>],
    n_old: usize,
    n_new: usize,
) -> CoordsQ {
    let mut out = vec![Ratio::zero(); n_new];
    for (i, c) in f.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let exp = 3 * (n_old - 1 - i);
        if exp < n_new {
            out[n_new - 1 - exp] += c.clone();
        }
    }
    out
}

/// **Stable** — `embed_in_square_extension`
pub fn embed_in_square_extension(
    f: &[Ratio<BigInt>],
    n_old: usize,
    n_new: usize,
) -> CoordsQ {
    let mut out = vec![Ratio::zero(); n_new];
    for (i, c) in f.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let exp = 2 * (n_old - 1 - i);
        if exp < n_new {
            out[n_new - 1 - exp] += c.clone();
        }
    }
    out
}

/// **Stable** — `identity_matrix`
pub fn identity_matrix(n: usize) -> Vec<Vec<Ratio<BigInt>>> {
    let mut out = vec![vec![Ratio::zero(); n]; n];
    for i in 0..n {
        out[i][i] = Ratio::one();
    }
    out
}

/// **Stable** — `mat_mul`
pub fn mat_mul(a: &[Vec<Ratio<BigInt>>], b: &[Vec<Ratio<BigInt>>]) -> Vec<Vec<Ratio<BigInt>>> {
    let n = a.len();
    let mut out = vec![vec![Ratio::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                out[i][j] += a[i][k].clone() * b[k][j].clone();
            }
        }
    }
    out
}

/// **Stable** — `kron_left`
pub fn kron_left(a: &[Vec<Ratio<BigInt>>], nb: usize) -> Vec<Vec<Ratio<BigInt>>> {
    let na = a.len();
    let dim = na * nb;
    let mut out = vec![vec![Ratio::zero(); dim]; dim];
    for jb in 0..nb {
        for ja in 0..na {
            for ia in 0..na {
                out[ja + jb * na][ia + jb * na] += a[ja][ia].clone();
            }
        }
    }
    out
}

/// **Stable** — `mult_matrix_of_element`
pub fn mult_matrix_of_element(u: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> Vec<Vec<Ratio<BigInt>>> {
    let n = poly_degree(m);
    let mut mat = vec![vec![Ratio::zero(); n]; n];
    for j in 0..n {
        let mut basis = vec![Ratio::zero(); n];
        basis[n - 1 - j] = Ratio::one();
        let prod = poly_reduce(&poly_mul(u, &basis), m);
        let prod = pad_to_len(&prod, n);
        for i in 0..n {
            mat[i][j] = prod[i].clone();
        }
    }
    mat
}

/// Parent operational basis vector `e_i` (length = `parent.dimension()`).
/// **Stable** — `parent_basis_vector`
pub fn parent_basis_vector(parent: &super::ext_tower::ExtensionField, i: usize) -> CoordsQ {
    let mut v = parent.zero_coords();
    if i < v.len() {
        v[i] = Ratio::one();
    }
    v
}

// **Pipeline private** — `tower_flat_index`
fn tower_flat_index(i: usize, u_power: usize, na: usize, nb: usize) -> usize {
    (nb - 1 - u_power) * na + i
}

// **Pipeline private** — `tower_blocks_to_flat`
fn tower_blocks_to_flat(blocks: &[CoordsQ], na: usize, nb: usize) -> CoordsQ {
    let mut out = vec![Ratio::zero(); na * nb];
    for j in 0..nb {
        let block = pad_to_len(&blocks[j], na);
        for i in 0..na {
            out[tower_flat_index(i, j, na, nb)] = block[i].clone();
        }
    }
    out
}

/// Multiply tower basis `e_i · u^j` by the layer generator `u` (T3+ compose helper).
// **Pipeline private** — `mul_tower_basis_by_layer_generator`
fn mul_tower_basis_by_layer_generator(
    parent: &super::ext_tower::ExtensionField,
    layer_blocks: &[CoordsQ],
    i: usize,
    j: usize,
    ring: &ParentCoeffRing<'_>,
) -> Result<CoordsQ, EvalError> {
    let na = parent.dimension();
    let nb = layer_blocks.len().checked_sub(1).ok_or(EvalError::TypeError("layer minpoly"))?;
    let ei = parent_basis_vector(parent, i);
    if j + 1 < nb {
        let mut blocks = vec![ring.zero.clone(); nb];
        blocks[j + 1] = ei;
        return Ok(tower_blocks_to_flat(&blocks, na, nb));
    }
    let mut blocks = vec![ring.zero.clone(); nb];
    for k in 1..=nb {
        let term = (ring.mul)(&ei, &layer_blocks[k])?;
        let neg = (ring.neg)(&term)?;
        let u_power = nb - k;
        blocks[u_power] = (ring.add)(&blocks[u_power], &neg)?;
    }
    Ok(tower_blocks_to_flat(&blocks, na, nb))
}

/// Multiplication-by-`u` matrix for `parent(u) / parent` with parent-coefficient layer minpoly.
// **Stable** — `mult_matrix_of_adjoin_generator`
pub(crate) fn mult_matrix_of_adjoin_generator(
    parent: &super::ext_tower::ExtensionField,
    layer_blocks: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<Vec<Ratio<BigInt>>>, EvalError> {
    let na = parent.dimension();
    let nb = layer_blocks
        .len()
        .checked_sub(1)
        .ok_or(EvalError::TypeError("layer minpoly"))?;
    if nb < 1 {
        return Err(EvalError::TypeError("extension degree >= 1"));
    }
    if !parent.element_eq_mod(&layer_blocks[0], &ring.one)? {
        return Err(EvalError::TypeError("layer minpoly must be monic"));
    }
    let dim = na * nb;
    let mut mat = vec![vec![Ratio::zero(); dim]; dim];
    for j in 0..nb {
        for i in 0..na {
            let col = tower_flat_index(i, j, na, nb);
            let image = mul_tower_basis_by_layer_generator(parent, layer_blocks, i, j, ring)?;
            for row in 0..dim {
                mat[row][col] = image[row].clone();
            }
        }
    }
    Ok(mat)
}

/// **Stable** — `char_poly_matrix`
pub fn char_poly_matrix(mat: &[Vec<Ratio<BigInt>>]) -> CoordsQ {
    let n = mat.len();
    let mut pow = identity_matrix(n);
    let mut traces = Vec::with_capacity(n);
    for _ in 0..n {
        pow = mat_mul(&pow, mat);
        let tr = (0..n).map(|i| pow[i][i].clone()).fold(Ratio::zero(), |a, b| a + b);
        traces.push(tr);
    }
    newton_char_poly(&traces, n)
}

// **Pipeline private** — `newton_char_poly`
fn newton_char_poly(traces: &[Ratio<BigInt>], n: usize) -> CoordsQ {
    let mut e = vec![Ratio::zero(); n + 1];
    e[0] = Ratio::one();
    for k in 1..=n {
        let mut ek = Ratio::zero();
        for j in 1..k {
            ek -= e[j].clone() * traces[k - j - 1].clone();
        }
        ek -= traces[k - 1].clone();
        ek /= Ratio::from_integer(BigInt::from(k as i64));
        e[k] = ek;
    }
    let mut out = vec![Ratio::one(); n + 1];
    for k in 1..=n {
        out[k] = if k % 2 == 1 {
            -e[k].clone()
        } else {
            e[k].clone()
        };
    }
    trim_leading_zero(out)
}

/// **Stable** — `apply_linear_map`
pub fn apply_linear_map(matrix: &[Vec<Ratio<BigInt>>], v: &[Ratio<BigInt>]) -> CoordsQ {
    let dim = matrix.len();
    let mut out = vec![Ratio::zero(); dim];
    for i in 0..dim {
        for j in 0..dim {
            if j < v.len() {
                out[i] += matrix[i][j].clone() * v[j].clone();
            }
        }
    }
    out
}

/// **Stable** — `coords_all_zero`
pub fn coords_all_zero(v: &[Ratio<BigInt>]) -> bool {
    v.iter().all(|c| c.is_zero())
}

/// **Stable** — `coords_is_one`
pub fn coords_is_one(v: &[Ratio<BigInt>]) -> bool {
    v.len() == 1 && v[0].is_one()
}

/// Coefficient ring for univariate polys over an extension (T3 parent-coeff arithmetic).
pub struct ParentCoeffRing<'a> {
    pub zero: CoordsQ,
    pub one: CoordsQ,
    pub add: &'a (dyn Fn(&CoordsQ, &CoordsQ) -> Result<CoordsQ, EvalError> + 'a),
    pub sub: &'a (dyn Fn(&CoordsQ, &CoordsQ) -> Result<CoordsQ, EvalError> + 'a),
    pub mul: &'a (dyn Fn(&CoordsQ, &CoordsQ) -> Result<CoordsQ, EvalError> + 'a),
    pub neg: &'a (dyn Fn(&CoordsQ) -> Result<CoordsQ, EvalError> + 'a),
    pub inv: &'a (dyn Fn(&CoordsQ) -> Result<CoordsQ, EvalError> + 'a),
    pub is_zero: &'a (dyn Fn(&CoordsQ) -> bool + 'a),
}

/// Adapter: parent-field coordinates as dense poly1 coefficients (T3).
pub struct ParentBlockRing<'a> {
    inner: &'a ParentCoeffRing<'a>,
}

impl<'a> ParentBlockRing<'a> {
    /// **Stable** — `new`
    pub fn new(inner: &'a ParentCoeffRing<'a>) -> Self {
        Self { inner }
    }
}

impl<'a> giac_poly::dense::Poly1RingCtx for ParentBlockRing<'a> {
    type Coeff = CoordsQ;

    // **Stable** — Poly zero
    fn zero(&self) -> CoordsQ {
        self.inner.zero.clone()
    }

    // **Stable** — Poly one
    fn one(&self) -> CoordsQ {
        self.inner.one.clone()
    }

    // **Stable** — Poly is zero
    fn is_zero(&self, c: &CoordsQ) -> bool {
        (self.inner.is_zero)(c)
    }

    // **Stable** — Poly addition
    fn add(&self, a: &CoordsQ, b: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.add)(a, b)
    }

    // **Stable** — Poly subtraction
    fn sub(&self, a: &CoordsQ, b: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.sub)(a, b)
    }

    // **Stable** — Poly negation
    fn neg(&self, c: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.neg)(c)
    }

    // **Stable** — Poly multiplication
    fn mul(&self, a: &CoordsQ, b: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.mul)(a, b)
    }

    // **Pipeline private** — `inv`
    fn inv(&self, c: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.inv)(c)
    }
}

/// Add polynomials with coefficients in `ring`.
/// **Stable** — `poly_add_with_coeffs_in_field`
pub fn poly_add_with_coeffs_in_field(
    a: &[CoordsQ],
    b: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::add(&ParentBlockRing::new(ring), a, b, POLY1_Q)
}

/// Multiply polynomials with coefficients in `ring`.
/// **Stable** — `poly_mul_with_coeffs_in_field`
pub fn poly_mul_with_coeffs_in_field(
    a: &[CoordsQ],
    b: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::mul(&ParentBlockRing::new(ring), a, b, POLY1_Q)
}

/// Reduce `p` modulo monic `m` (leading parent-coeff is `one`).
///
/// T3 layer minpoly from [`ExtensionField::layer_minpoly_parent_coeffs`](super::ext_tower::ExtensionField)
/// embeds ℚ minpoly into parent operational coords; **leading block is always `ring.one`**.
/// Generic [`giac_poly::dense::reduce_mod_monic`] also accepts leading **−1** (ℚ path); that branch
/// is unreachable for T3 modulus in practice — see
/// [GIAC-dense-poly1-refactor](.doc/issues/GIAC-dense-poly1-refactor.md) §4.4.
/// **Stable** — `poly_reduce_with_coeffs_in_field`
pub fn poly_reduce_with_coeffs_in_field(
    p: &[CoordsQ],
    m: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::reduce_mod_monic(&ParentBlockRing::new(ring), p, m, POLY1_Q)
}

/// **Stable** — `poly_sub_with_coeffs_in_field`
pub fn poly_sub_with_coeffs_in_field(
    a: &[CoordsQ],
    b: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::sub(&ParentBlockRing::new(ring), a, b, POLY1_Q)
}

/// **Stable** — `poly_neg_with_coeffs_in_field`
pub fn poly_neg_with_coeffs_in_field(
    a: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::neg(&ParentBlockRing::new(ring), a, POLY1_Q)
}

/// **Stable** — `poly_inv_mod_with_coeffs_in_field`
pub fn poly_inv_mod_with_coeffs_in_field(
    a: &[CoordsQ],
    m: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::inv_mod(&ParentBlockRing::new(ring), a, m, POLY1_Q)
}
