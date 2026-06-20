//! Dense univariate polynomial arithmetic over a generic coefficient ring.
//!
//! Coefficient vectors use giac **`poly1`** layout when [`Poly1Order::HighFirst`]:
//! index `0` is the leading (highest-degree) term; the constant term is last.
//! See [GIAC-algext-adoption](../../.doc/issues/GIAC-algext-adoption.md) and
//! [GIAC-dense-poly1-refactor](../../.doc/issues/GIAC-dense-poly1-refactor.md).
//!
//! **Tier:** Stable (crate-internal).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use crate::error::{EvalError, PolyResult};

/// Coefficient order for dense univariate polynomials.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Poly1Order {
    /// giac `poly1` / `ext_tower` [`CoordsQ`](../../giac-core/src/algebra/field_arith.rs): high degree first.
    HighFirst,
    /// giac-poly univariate ascending layout: constant term first.
    Ascending,
}

/// Ring context for dense poly1 coefficients (instance-based: supports parent-block rings).
///
/// **Stable (crate-internal)** — distinct from sparse [`PolyCoeff`](crate::PolyCoeff).
pub trait Poly1RingCtx {
    type Coeff: Clone;
    // **Stable** — Poly zero
    fn zero(&self) -> Self::Coeff;
    // **Stable** — Poly one
    fn one(&self) -> Self::Coeff;
    // **Stable** — Poly is zero
    fn is_zero(&self, c: &Self::Coeff) -> bool;
    // **Stable** — Poly addition
    fn add(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff>;
    // **Stable** — Poly subtraction
    fn sub(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff>;
    // **Stable** — Poly negation
    fn neg(&self, c: &Self::Coeff) -> PolyResult<Self::Coeff>;
    // **Stable** — Poly multiplication
    fn mul(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff>;
    // **Pipeline private** — `inv`
    fn inv(&self, c: &Self::Coeff) -> PolyResult<Self::Coeff>;
}

// **Stable** — Poly is one
fn is_one<R: Poly1RingCtx>(ctx: &R, c: &R::Coeff) -> bool {
    ctx.sub(c, &ctx.one())
        .map(|d| ctx.is_zero(&d))
        .unwrap_or(false)
}

// **Pipeline private** — `div_coeff`
fn div_coeff<R: Poly1RingCtx>(ctx: &R, a: &R::Coeff, b: &R::Coeff) -> PolyResult<R::Coeff> {
    ctx.mul(a, &ctx.inv(b)?)
}

// **Pipeline private** — `to_high_first`
fn to_high_first<C: Clone>(p: &[C], order: Poly1Order) -> Vec<C> {
    match order {
        Poly1Order::HighFirst => p.to_vec(),
        Poly1Order::Ascending => p.iter().rev().cloned().collect(),
    }
}

// **Pipeline private** — `from_high_first`
fn from_high_first<C: Clone>(p: Vec<C>, order: Poly1Order) -> Vec<C> {
    match order {
        Poly1Order::HighFirst => p,
        Poly1Order::Ascending => p.into_iter().rev().collect(),
    }
}

/// **Stable (crate-internal)** — degree in the given order (zero poly has degree 0).
pub fn poly_degree<C: Clone>(p: &[C]) -> usize {
    p.len().saturating_sub(1)
}

// **Pipeline private** — `trim_high_first`
fn trim_high_first<R: Poly1RingCtx>(ctx: &R, v: &mut Vec<R::Coeff>) {
    while v.len() > 1 && v.first().is_some_and(|c| ctx.is_zero(c)) {
        v.remove(0);
    }
    if v.is_empty() {
        v.push(ctx.zero());
    }
}

/// **Stable (crate-internal)** — drop redundant leading zeros; keep at least one coefficient.
pub fn trim<R: Poly1RingCtx>(ctx: &R, p: &mut Vec<R::Coeff>, order: Poly1Order) {
    let mut work = to_high_first(p, order);
    trim_high_first(ctx, &mut work);
    *p = from_high_first(work, order);
}

// **Pipeline private** — `trim_high_first_collect`
fn trim_high_first_collect<R: Poly1RingCtx>(ctx: &R, mut v: Vec<R::Coeff>) -> Vec<R::Coeff> {
    trim_high_first(ctx, &mut v);
    v
}

/// **Stable (crate-internal)**
pub fn add<R: Poly1RingCtx>(
    ctx: &R,
    a: &[R::Coeff],
    b: &[R::Coeff],
    order: Poly1Order,
) -> PolyResult<Vec<R::Coeff>> {
    let a = to_high_first(a, order);
    let b = to_high_first(b, order);
    let da = poly_degree(&a);
    let db = poly_degree(&b);
    let d = da.max(db);
    let mut out = vec![ctx.zero(); d + 1];
    for (i, c) in a.iter().enumerate() {
        out[d - da + i] = ctx.add(&out[d - da + i], c)?;
    }
    for (i, c) in b.iter().enumerate() {
        out[d - db + i] = ctx.add(&out[d - db + i], c)?;
    }
    trim_high_first(ctx, &mut out);
    Ok(from_high_first(out, order))
}

/// **Stable (crate-internal)**
pub fn sub<R: Poly1RingCtx>(
    ctx: &R,
    a: &[R::Coeff],
    b: &[R::Coeff],
    order: Poly1Order,
) -> PolyResult<Vec<R::Coeff>> {
    let a = to_high_first(a, order);
    let b = to_high_first(b, order);
    let da = poly_degree(&a);
    let db = poly_degree(&b);
    let d = da.max(db);
    let mut out = vec![ctx.zero(); d + 1];
    for (i, c) in a.iter().enumerate() {
        out[d - da + i] = ctx.add(&out[d - da + i], c)?;
    }
    for (i, c) in b.iter().enumerate() {
        out[d - db + i] = ctx.sub(&out[d - db + i], c)?;
    }
    trim_high_first(ctx, &mut out);
    Ok(from_high_first(out, order))
}

/// **Stable (crate-internal)**
pub fn mul<R: Poly1RingCtx>(
    ctx: &R,
    a: &[R::Coeff],
    b: &[R::Coeff],
    order: Poly1Order,
) -> PolyResult<Vec<R::Coeff>> {
    let a = to_high_first(a, order);
    let b = to_high_first(b, order);
    if a.is_empty() || b.is_empty() {
        return Ok(from_high_first(vec![ctx.zero()], order));
    }
    let da = poly_degree(&a);
    let db = poly_degree(&b);
    let mut out = vec![ctx.zero(); da + db + 1];
    for (i, ca) in a.iter().enumerate() {
        for (j, cb) in b.iter().enumerate() {
            let prod = ctx.mul(ca, cb)?;
            out[i + j] = ctx.add(&out[i + j], &prod)?;
        }
    }
    trim_high_first(ctx, &mut out);
    Ok(from_high_first(out, order))
}

/// **Stable (crate-internal)**
pub fn neg<R: Poly1RingCtx>(ctx: &R, a: &[R::Coeff], order: Poly1Order) -> PolyResult<Vec<R::Coeff>> {
    let mut out: Vec<R::Coeff> = a.iter().map(|c| ctx.neg(c)).collect::<PolyResult<_>>()?;
    trim(ctx, &mut out, order);
    Ok(out)
}

/// **Stable (crate-internal)**
pub fn scale<R: Poly1RingCtx>(
    ctx: &R,
    a: &[R::Coeff],
    s: &R::Coeff,
    order: Poly1Order,
) -> PolyResult<Vec<R::Coeff>> {
    if ctx.is_zero(s) {
        return Ok(vec![ctx.zero()]);
    }
    let mut out: Vec<R::Coeff> = a.iter().map(|c| ctx.mul(c, s)).collect::<PolyResult<_>>()?;
    trim(ctx, &mut out, order);
    Ok(out)
}

/// **Stable (crate-internal)** — polynomial division; quotient and remainder in `order`.
pub fn div_rem<R: Poly1RingCtx>(
    ctx: &R,
    a: &[R::Coeff],
    b: &[R::Coeff],
    order: Poly1Order,
) -> PolyResult<(Vec<R::Coeff>, Vec<R::Coeff>)> {
    let mut rem = to_high_first(a, order);
    let b = to_high_first(b, order);
    trim_high_first(ctx, &mut rem);
    if b.iter().all(|c| ctx.is_zero(c)) {
        return Ok((
            from_high_first(vec![ctx.zero()], order),
            from_high_first(rem, order),
        ));
    }
    let db = poly_degree(&b);
    let da = poly_degree(&rem);
    if da < db {
        return Ok((
            from_high_first(vec![ctx.zero()], order),
            from_high_first(rem, order),
        ));
    }
    let lc_b = b.first().cloned().unwrap_or_else(|| ctx.zero());
    let orig_da = da;
    let mut quot = vec![ctx.zero(); da - db + 1];
    while poly_degree(&rem) >= db && !rem.iter().all(|c| ctx.is_zero(c)) {
        let dr = poly_degree(&rem);
        let lc_r = rem.first().cloned().unwrap_or_else(|| ctx.zero());
        if ctx.is_zero(&lc_r) {
            trim_high_first(ctx, &mut rem);
            continue;
        }
        let q = div_coeff(ctx, &lc_r, &lc_b)?;
        let qi = orig_da - dr;
        if qi < quot.len() {
            quot[qi] = q.clone();
        }
        for i in 0..=db {
            if i < rem.len() {
                let sub = ctx.mul(&q, &b[i])?;
                rem[i] = ctx.sub(&rem[i], &sub)?;
            }
        }
        trim_high_first(ctx, &mut rem);
    }
    Ok((
        from_high_first(trim_high_first_collect(ctx, quot), order),
        from_high_first(rem, order),
    ))
}

/// **Stable (crate-internal)** — reduce `p` modulo monic `m` (leading coeff of `m` is ±1).
pub fn reduce_mod_monic<R: Poly1RingCtx>(
    ctx: &R,
    p: &[R::Coeff],
    m: &[R::Coeff],
    order: Poly1Order,
) -> PolyResult<Vec<R::Coeff>> {
    let mut r = to_high_first(p, order);
    let m = to_high_first(m, order);
    let dm = poly_degree(&m);
    if dm == 0 {
        trim_high_first(ctx, &mut r);
        return Ok(from_high_first(r, order));
    }
    let lead = m.first().cloned().unwrap_or_else(|| ctx.zero());
    if !is_one(ctx, &lead) {
        let neg_one = ctx.neg(&ctx.one())?;
        if !ctx.is_zero(&ctx.sub(&lead, &neg_one)?) {
            trim_high_first(ctx, &mut r);
            return Ok(from_high_first(r, order));
        }
    }
    loop {
        trim_high_first(ctx, &mut r);
        let dr = poly_degree(&r);
        if dr < dm {
            break;
        }
        let lc_r = r.first().cloned().unwrap_or_else(|| ctx.zero());
        let q = div_coeff(ctx, &lc_r, &lead)?;
        for i in 0..=dm {
            if i < r.len() {
                let sub = ctx.mul(&q, &m[i])?;
                r[i] = ctx.sub(&r[i], &sub)?;
            }
        }
    }
    trim_high_first(ctx, &mut r);
    Ok(from_high_first(r, order))
}

/// **Stable (crate-internal)** — extended GCD: `(g, s)` with `s*a + t*b = g` (`t` omitted).
pub fn ext_gcd<R: Poly1RingCtx>(
    ctx: &R,
    a: &[R::Coeff],
    b: &[R::Coeff],
    order: Poly1Order,
) -> PolyResult<(Vec<R::Coeff>, Vec<R::Coeff>)> {
    let mut r_prev = to_high_first(b, order);
    let mut r = to_high_first(a, order);
    trim_high_first(ctx, &mut r_prev);
    trim_high_first(ctx, &mut r);
    let mut s_prev = vec![ctx.zero()];
    let mut s = vec![ctx.one()];
    while !r.iter().all(|c| ctx.is_zero(c)) {
        let (q, _) = div_rem(ctx, &r_prev, &r, Poly1Order::HighFirst)?;
        let qr = mul(ctx, &q, &r, Poly1Order::HighFirst)?;
        let r_next = sub(ctx, &r_prev, &qr, Poly1Order::HighFirst)?;
        let sr = mul(ctx, &q, &s, Poly1Order::HighFirst)?;
        let s_next = sub(ctx, &s_prev, &sr, Poly1Order::HighFirst)?;
        r_prev = r;
        r = trim_high_first_collect(ctx, r_next);
        s_prev = s;
        s = s_next;
    }
    let g = r_prev;
    if let Some(lc) = g.first() {
        if !ctx.is_zero(lc) && !is_one(ctx, lc) {
            let inv_lc = ctx.inv(lc)?;
            let scale_vec = |p: &[R::Coeff]| -> PolyResult<Vec<R::Coeff>> {
                p.iter()
                    .map(|c| ctx.mul(c, &inv_lc))
                    .collect::<PolyResult<_>>()
            };
            return Ok((
                from_high_first(scale_vec(&g)?, order),
                from_high_first(scale_vec(&s_prev)?, order),
            ));
        }
    }
    Ok((
        from_high_first(g, order),
        from_high_first(s_prev, order),
    ))
}

/// **Stable (crate-internal)** — multiplicative inverse of `a` modulo monic `m`.
pub fn inv_mod<R: Poly1RingCtx>(
    ctx: &R,
    a: &[R::Coeff],
    m: &[R::Coeff],
    order: Poly1Order,
) -> PolyResult<Vec<R::Coeff>> {
    let (_, bezout) = ext_gcd(ctx, a, m, order)?;
    let inv = reduce_mod_monic(ctx, &bezout, m, order)?;
    let prod = mul(ctx, a, &inv, order)?;
    let check = reduce_mod_monic(ctx, &prod, m, order)?;
    if check.len() == 1 && is_one(ctx, &check[0]) {
        Ok(inv)
    } else {
        Err(EvalError::NotImplemented("field inverse"))
    }
}
