//! Univariate polynomial arithmetic over ℚ in giac `poly1` convention (high degree first).
//!
//! Shared by [`super::ext_tower::ExtensionField`] and [`super::alg_ext::AlgExtData`].
//! Dense ℚ ring ops delegate to [`giac_poly::dense::poly1`] (see GIAC-dense-poly1-refactor D2).

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

fn dense_q<T>(r: giac_poly::PolyResult<T>) -> T {
    r.expect("dense poly1 over Q is infallible except inv_mod")
}

pub fn trim_leading_zero(mut v: CoordsQ) -> CoordsQ {
    dense::trim(&RATIO_CTX, &mut v, POLY1_Q);
    v
}

pub fn poly_degree(p: &[Ratio<BigInt>]) -> usize {
    dense::poly_degree(p)
}

pub fn poly_add(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::add(&RATIO_CTX, a, b, POLY1_Q))
}

pub fn poly_sub(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::sub(&RATIO_CTX, a, b, POLY1_Q))
}

pub fn poly_mul(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::mul(&RATIO_CTX, a, b, POLY1_Q))
}

pub fn poly_neg(a: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::neg(&RATIO_CTX, a, POLY1_Q))
}

pub fn poly_scale(a: &[Ratio<BigInt>], s: &Ratio<BigInt>) -> CoordsQ {
    dense_q(dense::scale(&RATIO_CTX, a, s, POLY1_Q))
}

/// Reduce `p` modulo monic `m` (high-degree-first, leading coeff of `m` is ±1).
pub fn poly_reduce(p: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> CoordsQ {
    dense_q(dense::reduce_mod_monic(&RATIO_CTX, p, m, POLY1_Q))
}

pub fn poly_inv_mod(a: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> Result<CoordsQ, EvalError> {
    dense::inv_mod(&RATIO_CTX, a, m, POLY1_Q)
}

/// Extended GCD: returns `(g, s)` with `s*a + t*b = g` (only `s` needed for inverse).
pub fn poly_ext_gcd(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (CoordsQ, CoordsQ) {
    dense_q(dense::ext_gcd(&RATIO_CTX, a, b, POLY1_Q))
}

pub fn poly_divrem(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (CoordsQ, CoordsQ) {
    dense_q(dense::div_rem(&RATIO_CTX, a, b, POLY1_Q))
}

pub fn pad_to_len(v: &[Ratio<BigInt>], n: usize) -> CoordsQ {
    if v.len() >= n {
        return v[v.len() - n..].to_vec();
    }
    let mut out = vec![Ratio::zero(); n - v.len()];
    out.extend(v.iter().cloned());
    out
}

pub fn generator_coords(n: usize) -> CoordsQ {
    let mut v = vec![Ratio::zero(); n];
    if n > 0 {
        v[0] = Ratio::one();
    }
    v
}

pub fn expr_to_ratio(e: &Expr) -> Result<Ratio<BigInt>, EvalError> {
    match e {
        Expr::Int(n) => Ok(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Ok(r.clone()),
        _ => Err(EvalError::TypeError("rational coeff expected")),
    }
}

pub fn rationalize_poly1(items: &[ExprArc]) -> Result<CoordsQ, EvalError> {
    items.iter().map(|e| expr_to_ratio(e.as_ref())).collect()
}

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

pub fn ratio_to_expr_arc(r: &Ratio<BigInt>) -> ExprArc {
    crate::ratio_to_expr(r)
}

pub fn min_poly_exprs_to_q(min_poly: &[ExprArc]) -> Result<CoordsQ, EvalError> {
    rationalize_poly1(min_poly)
}

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

pub fn minpoly_at_square(m: &[Ratio<BigInt>]) -> CoordsQ {
    let n = poly_degree(m);
    let mut out = vec![Ratio::zero(); 2 * n + 1];
    for i in 0..=n {
        out[2 * i] += m[i].clone();
    }
    trim_leading_zero(out)
}

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

pub fn identity_matrix(n: usize) -> Vec<Vec<Ratio<BigInt>>> {
    let mut out = vec![vec![Ratio::zero(); n]; n];
    for i in 0..n {
        out[i][i] = Ratio::one();
    }
    out
}

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

pub fn coords_all_zero(v: &[Ratio<BigInt>]) -> bool {
    v.iter().all(|c| c.is_zero())
}

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
    pub fn new(inner: &'a ParentCoeffRing<'a>) -> Self {
        Self { inner }
    }
}

impl<'a> giac_poly::dense::Poly1RingCtx for ParentBlockRing<'a> {
    type Coeff = CoordsQ;

    fn zero(&self) -> CoordsQ {
        self.inner.zero.clone()
    }

    fn one(&self) -> CoordsQ {
        self.inner.one.clone()
    }

    fn is_zero(&self, c: &CoordsQ) -> bool {
        (self.inner.is_zero)(c)
    }

    fn add(&self, a: &CoordsQ, b: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.add)(a, b)
    }

    fn sub(&self, a: &CoordsQ, b: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.sub)(a, b)
    }

    fn neg(&self, c: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.neg)(c)
    }

    fn mul(&self, a: &CoordsQ, b: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.mul)(a, b)
    }

    fn inv(&self, c: &CoordsQ) -> giac_poly::PolyResult<CoordsQ> {
        (self.inner.inv)(c)
    }
}

/// Add polynomials with coefficients in `ring`.
pub fn poly_add_with_coeffs_in_field(
    a: &[CoordsQ],
    b: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::add(&ParentBlockRing::new(ring), a, b, POLY1_Q)
}

/// Multiply polynomials with coefficients in `ring`.
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
pub fn poly_reduce_with_coeffs_in_field(
    p: &[CoordsQ],
    m: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::reduce_mod_monic(&ParentBlockRing::new(ring), p, m, POLY1_Q)
}

pub fn poly_sub_with_coeffs_in_field(
    a: &[CoordsQ],
    b: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::sub(&ParentBlockRing::new(ring), a, b, POLY1_Q)
}

pub fn poly_neg_with_coeffs_in_field(
    a: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::neg(&ParentBlockRing::new(ring), a, POLY1_Q)
}

pub fn poly_inv_mod_with_coeffs_in_field(
    a: &[CoordsQ],
    m: &[CoordsQ],
    ring: &ParentCoeffRing<'_>,
) -> Result<Vec<CoordsQ>, EvalError> {
    dense::inv_mod(&ParentBlockRing::new(ring), a, m, POLY1_Q)
}
