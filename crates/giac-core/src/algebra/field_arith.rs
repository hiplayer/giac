//! Univariate polynomial arithmetic over ℚ in giac `poly1` convention (high degree first).
//!
//! Shared by [`super::ext_tower::ExtensionField`] and [`super::alg_ext::AlgExtData`].

use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc, FuncKind};

/// Coefficient vector in giac `poly1` order (leading term first).
pub type CoordsQ = Vec<Ratio<BigInt>>;
pub type CoordsQSlice<'a> = &'a [Ratio<BigInt>];

pub fn trim_leading_zero(mut v: CoordsQ) -> CoordsQ {
    while v.len() > 1 && matches!(v.first(), Some(c) if c.is_zero()) {
        v.remove(0);
    }
    if v.is_empty() {
        vec![Ratio::zero()]
    } else {
        v
    }
}

pub fn poly_degree(p: &[Ratio<BigInt>]) -> usize {
    p.len().saturating_sub(1)
}

pub fn poly_add(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    let da = poly_degree(a);
    let db = poly_degree(b);
    let d = da.max(db);
    let mut out = vec![Ratio::zero(); d + 1];
    for (i, c) in a.iter().enumerate() {
        out[d - da + i] += c;
    }
    for (i, c) in b.iter().enumerate() {
        out[d - db + i] += c;
    }
    trim_leading_zero(out)
}

pub fn poly_sub(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    let da = poly_degree(a);
    let db = poly_degree(b);
    let d = da.max(db);
    let mut out = vec![Ratio::zero(); d + 1];
    for (i, c) in a.iter().enumerate() {
        out[d - da + i] += c;
    }
    for (i, c) in b.iter().enumerate() {
        out[d - db + i] -= c;
    }
    trim_leading_zero(out)
}

pub fn poly_mul(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> CoordsQ {
    if a.is_empty() || b.is_empty() {
        return vec![Ratio::zero()];
    }
    let da = poly_degree(a);
    let db = poly_degree(b);
    let mut out = vec![Ratio::zero(); da + db + 1];
    for (i, ca) in a.iter().enumerate() {
        for (j, cb) in b.iter().enumerate() {
            out[i + j] += ca * cb;
        }
    }
    trim_leading_zero(out)
}

pub fn poly_neg(a: &[Ratio<BigInt>]) -> CoordsQ {
    trim_leading_zero(a.iter().map(|c| -c).collect())
}

pub fn poly_scale(a: &[Ratio<BigInt>], s: &Ratio<BigInt>) -> CoordsQ {
    if s.is_zero() {
        return vec![Ratio::zero()];
    }
    trim_leading_zero(a.iter().map(|c| c * s).collect())
}

/// Reduce `p` modulo monic `m` (high-degree-first, leading coeff of `m` is ±1).
pub fn poly_reduce(p: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> CoordsQ {
    let mut r = p.to_vec();
    let dm = poly_degree(m);
    if dm == 0 {
        return trim_leading_zero(r);
    }
    if !matches!(m.first(), Some(c) if *c == Ratio::one() || *c == Ratio::from_integer((-1).into())) {
        return trim_leading_zero(r);
    }
    loop {
        r = trim_leading_zero(r);
        let dr = poly_degree(&r);
        if dr < dm {
            break;
        }
        let q = r.first().cloned().unwrap_or_else(Ratio::zero)
            / m.first().cloned().unwrap_or_else(Ratio::one);
        for i in 0..=dm {
            if i < r.len() {
                r[i] -= &q * &m[i];
            }
        }
    }
    trim_leading_zero(r)
}

pub fn poly_inv_mod(a: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> Result<CoordsQ, EvalError> {
    let (_, bezout) = poly_ext_gcd(a, m);
    let inv = poly_reduce(&bezout, m);
    let check = poly_reduce(&poly_mul(a, &inv), m);
    if check.len() == 1 && check[0].is_one() {
        Ok(inv)
    } else {
        Err(EvalError::NotImplemented("field inverse"))
    }
}

/// Extended GCD: returns `(g, s)` with `s*a + t*b = g` (only `s` needed for inverse).
pub fn poly_ext_gcd(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (CoordsQ, CoordsQ) {
    let mut r_prev = trim_leading_zero(b.to_vec());
    let mut r = trim_leading_zero(a.to_vec());
    let mut s_prev = vec![Ratio::zero()];
    let mut s = vec![Ratio::one()];
    while !r.iter().all(|c| c.is_zero()) {
        let (q, _) = poly_divrem(&r_prev, &r);
        let qr = poly_mul(&q, &r);
        let r_next = poly_sub(&r_prev, &qr);
        let sr = poly_mul(&q, &s);
        let s_next = poly_sub(&s_prev, &sr);
        r_prev = r;
        r = trim_leading_zero(r_next);
        s_prev = s;
        s = s_next;
    }
    let g = r_prev;
    if let Some(lc) = g.first().cloned() {
        if !lc.is_zero() && lc != Ratio::one() {
            let inv_lc = Ratio::one() / lc;
            let scale = |p: &[Ratio<BigInt>]| p.iter().map(|c| c * &inv_lc).collect::<Vec<_>>();
            return (scale(&g), scale(&s_prev));
        }
    }
    (g, s_prev)
}

pub fn poly_divrem(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (CoordsQ, CoordsQ) {
    let mut rem = trim_leading_zero(a.to_vec());
    let b = trim_leading_zero(b.to_vec());
    if b.iter().all(|c| c.is_zero()) {
        return (vec![Ratio::zero()], rem);
    }
    let db = poly_degree(&b);
    let da = poly_degree(&rem);
    if da < db {
        return (vec![Ratio::zero()], rem);
    }
    let lc_b = b.first().cloned().unwrap_or_else(Ratio::one);
    let orig_da = da;
    let mut quot = vec![Ratio::zero(); da - db + 1];
    while poly_degree(&rem) >= db && !rem.iter().all(|c| c.is_zero()) {
        let dr = poly_degree(&rem);
        let lc_r = rem.first().cloned().unwrap_or_else(Ratio::zero);
        if lc_r.is_zero() {
            rem = trim_leading_zero(rem);
            continue;
        }
        let q = lc_r / lc_b.clone();
        let qi = orig_da - dr;
        if qi < quot.len() {
            quot[qi] = q.clone();
        }
        for i in 0..=db {
            if i < rem.len() {
                rem[i] -= &q * &b[i];
            }
        }
        rem = trim_leading_zero(rem);
    }
    (trim_leading_zero(quot), rem)
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
    if r.denom() == &BigInt::one() {
        Expr::int(r.numer().to_string().parse().unwrap_or(0))
    } else {
        Expr::rat(
            r.numer().to_string().parse().unwrap_or(0),
            r.denom().to_string().parse().unwrap_or(1),
        )
    }
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
