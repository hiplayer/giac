//! Algebraic extension elements (upstream `_EXT` / `ref_algext`).
//!
//! An `AlgExtData` is an element of ℚ(α) where α is a root of `min_poly`.
//! Coefficients are stored as a univariate polynomial in α (high degree first,
//! giac `poly1` convention).

use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc, FuncKind};

/// Element of an algebraic extension ℚ(α) ≅ ℚ[x]/(min_poly).
#[derive(Clone, Debug, PartialEq)]
pub struct AlgExtData {
    /// Minimum polynomial coefficients, high degree first (`poly1` / giac style).
    pub min_poly: Vec<ExprArc>,
    /// Element as polynomial in α, same coefficient order.
    pub coords: Vec<ExprArc>,
    /// Optional `rootof` branch index (upstream `select_root` hint).
    pub root_index: Option<u32>,
}

impl AlgExtData {
    /// Constant `0` in ℚ(α).
    pub fn zero(min_poly: Vec<ExprArc>) -> Self {
        Self {
            min_poly,
            coords: vec![Expr::int(0)],
            root_index: None,
        }
    }

    /// Constant `1` in ℚ(α).
    pub fn one(min_poly: Vec<ExprArc>) -> Self {
        Self {
            min_poly,
            coords: vec![Expr::int(1)],
            root_index: None,
        }
    }

    pub fn into_expr(self) -> ExprArc {
        Arc::new(Expr::AlgExt(Arc::new(self)))
    }

    /// Build from evaluated `rootof(num, minpoly)` arguments.
    pub fn from_rootof(num: &ExprArc, min_poly: &ExprArc) -> Result<Self, EvalError> {
        let min = poly1_coeffs(min_poly)?;
        if min.len() < 2 {
            return Err(EvalError::TypeError("alg ext min poly degree >= 1"));
        }
        let coords = match num.as_ref() {
            Expr::Seq(items) => rationalize_poly1(items)?,
            Expr::Int(_) | Expr::Rat(_) => rationalize_poly1(std::slice::from_ref(num))?,
            _ => return Err(EvalError::TypeError("rootof numerator")),
        };
        Ok(Self {
            min_poly: min,
            coords: coords_to_expr(&coords)?,
            root_index: None,
        })
    }

    pub fn add(&self, other: &Self) -> Result<Self, EvalError> {
        ensure_same_field(self, other)?;
        let a = rationalize_poly1_expr(&self.coords)?;
        let b = rationalize_poly1_expr(&other.coords)?;
        let m = rationalize_poly1_expr(&self.min_poly)?;
        let sum = poly_add(&a, &b);
        Ok(Self {
            min_poly: self.min_poly.clone(),
            coords: coords_to_expr(&poly_reduce(&sum, &m))?,
            root_index: None,
        })
    }

    pub fn sub(&self, other: &Self) -> Result<Self, EvalError> {
        ensure_same_field(self, other)?;
        let a = rationalize_poly1_expr(&self.coords)?;
        let b = rationalize_poly1_expr(&other.coords)?;
        let m = rationalize_poly1_expr(&self.min_poly)?;
        let diff = poly_sub(&a, &b);
        Ok(Self {
            min_poly: self.min_poly.clone(),
            coords: coords_to_expr(&poly_reduce(&diff, &m))?,
            root_index: None,
        })
    }

    pub fn mul(&self, other: &Self) -> Result<Self, EvalError> {
        ensure_same_field(self, other)?;
        let a = rationalize_poly1_expr(&self.coords)?;
        let b = rationalize_poly1_expr(&other.coords)?;
        let m = rationalize_poly1_expr(&self.min_poly)?;
        let prod = poly_mul(&a, &b);
        Ok(Self {
            min_poly: self.min_poly.clone(),
            coords: coords_to_expr(&poly_reduce(&prod, &m))?,
            root_index: None,
        })
    }

    /// Multiplicative inverse in ℚ(α) (`inv_EXT`).
    pub fn inv(&self) -> Result<Self, EvalError> {
        if self.is_zero() {
            return Err(EvalError::DivisionByZero);
        }
        let a = rationalize_poly1_expr(&self.coords)?;
        let m = rationalize_poly1_expr(&self.min_poly)?;
        let inv = poly_inv_mod(&a, &m)?;
        Ok(Self {
            min_poly: self.min_poly.clone(),
            coords: coords_to_expr(&inv)?,
            root_index: None,
        })
    }

    /// Same extension field and equal coordinates mod `min_poly`.
    pub fn eq_mod(&self, other: &Self) -> Result<bool, EvalError> {
        if self.min_poly != other.min_poly {
            return Ok(false);
        }
        let a = rationalize_poly1_expr(&self.coords)?;
        let b = rationalize_poly1_expr(&other.coords)?;
        let m = rationalize_poly1_expr(&self.min_poly)?;
        Ok(poly_reduce(&poly_sub(&a, &b), &m)
            .iter()
            .all(|c| c.is_zero()))
    }

    /// Scale by a rational scalar (same `min_poly`).
    pub fn mul_rational(&self, r: &Ratio<BigInt>) -> Result<Self, EvalError> {
        if r.is_zero() {
            return Ok(Self::zero(self.min_poly.clone()));
        }
        if r.is_one() {
            return Ok(self.clone());
        }
        let m = rationalize_poly1_expr(&self.min_poly)?;
        let mut coords = rationalize_poly1_expr(&self.coords)?;
        for c in &mut coords {
            *c *= r;
        }
        Ok(Self {
            min_poly: self.min_poly.clone(),
            coords: coords_to_expr(&poly_reduce(&coords, &m))?,
            root_index: None,
        })
    }

    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|c| c.is_zero())
    }

    pub fn is_one(&self) -> bool {
        matches!(self.coords.as_slice(), [c] if c.is_one())
    }

    /// Giac-compatible `rootof([...], poly1[...])` when representable.
    pub fn to_rootof_expr(&self) -> ExprArc {
        let coords = canonical_poly1_expr(&self.coords);
        let min_poly = canonical_poly1_expr(&self.min_poly);
        let num = Arc::new(Expr::Seq(coords));
        let min = if min_poly.len() == 1
            && matches!(min_poly[0].as_ref(), Expr::Func(FuncKind::Poly1, _))
        {
            Arc::clone(&min_poly[0])
        } else {
            Arc::new(Expr::Func(
                FuncKind::Poly1,
                vec![Arc::new(Expr::Seq(min_poly))],
            ))
        };
        Expr::func(FuncKind::RootOf, vec![num, min])
    }
}

/// Trim leading zero coefficients while keeping at least one term.
fn canonical_poly1_expr(coeffs: &[ExprArc]) -> Vec<ExprArc> {
    if coeffs.is_empty() {
        return vec![Expr::int(0)];
    }
    let start = coeffs
        .iter()
        .position(|c| !c.is_zero())
        .unwrap_or(coeffs.len() - 1);
    coeffs[start..].to_vec()
}

fn ensure_same_field(a: &AlgExtData, b: &AlgExtData) -> Result<(), EvalError> {
    if a.min_poly != b.min_poly {
        return Err(EvalError::TypeError("algebraic extension mismatch"));
    }
    Ok(())
}

fn poly1_coeffs(e: &ExprArc) -> Result<Vec<ExprArc>, EvalError> {
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

fn expr_to_ratio(e: &Expr) -> Result<Ratio<BigInt>, EvalError> {
    match e {
        Expr::Int(n) => Ok(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Ok(r.clone()),
        _ => Err(EvalError::TypeError("rational coeff expected")),
    }
}

fn rationalize_poly1(items: &[ExprArc]) -> Result<Vec<Ratio<BigInt>>, EvalError> {
    items.iter().map(|e| expr_to_ratio(e.as_ref())).collect()
}

fn rationalize_poly1_expr(items: &[ExprArc]) -> Result<Vec<Ratio<BigInt>>, EvalError> {
    rationalize_poly1(items)
}

fn coords_to_expr(coords: &[Ratio<BigInt>]) -> Result<Vec<ExprArc>, EvalError> {
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

fn trim_leading_zero(mut v: Vec<Ratio<BigInt>>) -> Vec<Ratio<BigInt>> {
    while v.len() > 1 && matches!(v.first(), Some(c) if c.is_zero()) {
        v.remove(0);
    }
    if v.is_empty() {
        vec![Ratio::zero()]
    } else {
        v
    }
}

fn poly_degree(p: &[Ratio<BigInt>]) -> usize {
    p.len().saturating_sub(1)
}

fn poly_add(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>> {
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

fn poly_sub(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>> {
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

fn poly_mul(a: &[Ratio<BigInt>], b: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>> {
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

fn poly_inv_mod(a: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> Result<Vec<Ratio<BigInt>>, EvalError> {
    let (_, bezout) = poly_ext_gcd(a, m);
    let inv = poly_reduce(&bezout, m);
    let check = poly_reduce(&poly_mul(a, &inv), m);
    if check.len() == 1 && check[0].is_one() {
        Ok(inv)
    } else {
        Err(EvalError::NotImplemented("alg ext inverse"))
    }
}

/// Extended GCD for `a`, `b`: returns `(g, s)` with `s*a + t*b = g` (only `s` needed for inv).
fn poly_ext_gcd(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (Vec<Ratio<BigInt>>, Vec<Ratio<BigInt>>) {
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

fn poly_divrem(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (Vec<Ratio<BigInt>>, Vec<Ratio<BigInt>>) {
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
        // Subtract q * x^(dr-db) * b: term b[i]*x^(db-i) shifts to degree dr-i → index i.
        for i in 0..=db {
            if i < rem.len() {
                rem[i] -= &q * &b[i];
            }
        }
        rem = trim_leading_zero(rem);
    }
    (trim_leading_zero(quot), rem)
}

/// Reduce `p` modulo monic `m` (high-degree-first, leading coeff of `m` is 1).
fn poly_reduce(p: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>> {
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

/// Fold a sum that may contain `AlgExt` leaves (same `min_poly` only).
pub fn fold_algext_sum(terms: &[ExprArc]) -> Result<ExprArc, EvalError> {
    let mut groups: Vec<(Vec<ExprArc>, AlgExtData)> = Vec::new();
    let mut rat_sum = Ratio::<BigInt>::zero();
    let mut rest = Vec::new();
    for t in terms {
        match t.as_ref() {
            Expr::AlgExt(a) => {
                if let Some((_, acc)) = groups
                    .iter_mut()
                    .find(|(k, _)| k == &a.min_poly)
                {
                    *acc = acc.add(a)?;
                } else {
                    groups.push((a.min_poly.clone(), (**a).clone()));
                }
            }
            Expr::Int(n) => rat_sum += Ratio::from_integer(n.clone()),
            Expr::Rat(r) => rat_sum += r.clone(),
            _ => rest.push(Arc::clone(t)),
        }
    }
    if !rat_sum.is_zero() {
        if groups.len() == 1 {
            let (key, acc) = &mut groups[0];
            let m = rationalize_poly1_expr(&acc.min_poly)?;
            let mut coords = rationalize_poly1_expr(&acc.coords)?;
            if let Some(c0) = coords.last_mut() {
                *c0 += rat_sum;
            } else {
                coords.push(rat_sum);
            }
            *acc = AlgExtData {
                min_poly: key.clone(),
                coords: coords_to_expr(&poly_reduce(&coords, &m))?,
                root_index: acc.root_index,
            };
        } else {
            rest.push(ratio_to_expr_arc(&rat_sum));
        }
    }
    for (_, a) in groups {
        if !a.is_zero() {
            rest.push(a.into_expr());
        }
    }
    match rest.len() {
        0 => Ok(Expr::int(0)),
        1 => Ok(Arc::clone(&rest[0])),
        _ => Ok(Expr::add(rest)),
    }
}

/// Fold a product that may contain `AlgExt` factors (same `min_poly` only).
pub fn fold_algext_product(factors: &[ExprArc]) -> Result<ExprArc, EvalError> {
    let mut acc_ext: Option<AlgExtData> = None;
    let mut rat_prod = Ratio::<BigInt>::one();
    let mut rest = Vec::new();
    for f in factors {
        match f.as_ref() {
            Expr::AlgExt(a) => {
                acc_ext = Some(match acc_ext {
                    None => (**a).clone(),
                    Some(e) => e.mul(a)?,
                });
            }
            Expr::Int(n) => {
                if n.is_zero() {
                    return Ok(Expr::int(0));
                }
                rat_prod *= Ratio::from_integer(n.clone());
            }
            Expr::Rat(r) => {
                if r.is_zero() {
                    return Ok(Expr::int(0));
                }
                rat_prod *= r.clone();
            }
            _ => rest.push(Arc::clone(f)),
        }
    }
    if let Some(mut e) = acc_ext {
        if !rat_prod.is_one() {
            let m = rationalize_poly1_expr(&e.min_poly)?;
            let mut coords = rationalize_poly1_expr(&e.coords)?;
            for c in &mut coords {
                *c *= &rat_prod;
            }
            e.coords = coords_to_expr(&poly_reduce(&coords, &m))?;
        }
        if !e.is_one() {
            rest.push(e.into_expr());
        }
    } else if !rat_prod.is_one() {
        rest.push(ratio_to_expr_arc(&rat_prod));
    }
    match rest.len() {
        0 => Ok(Expr::int(1)),
        1 => Ok(Arc::clone(&rest[0])),
        _ => Ok(Expr::mul(rest)),
    }
}

fn ratio_to_expr_arc(r: &Ratio<BigInt>) -> ExprArc {
    if r.denom() == &BigInt::one() {
        Expr::int(r.numer().to_string().parse().unwrap_or(0))
    } else {
        Expr::rat(
            r.numer().to_string().parse().unwrap_or(0),
            r.denom().to_string().parse().unwrap_or(1),
        )
    }
}

/// Try to construct `AlgExt` from a `rootof(...)` call.
pub fn try_rootof_to_algext(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("rootof"));
    }
    Ok(AlgExtData::from_rootof(&args[0], &args[1])?.into_expr())
}

/// True when `e` contains an `Expr::AlgExt` leaf.
pub fn contains_algext(e: &Expr) -> bool {
    match e {
        Expr::AlgExt(_) => true,
        Expr::Add(terms) | Expr::Mul(terms) | Expr::Seq(terms) | Expr::List(terms) => {
            terms.iter().any(|t| contains_algext(t.as_ref()))
        }
        Expr::Frac(n, d) | Expr::Pow(n, d) | Expr::Complex(n, d) | Expr::Mod(n, d) => {
            contains_algext(n.as_ref()) || contains_algext(d.as_ref())
        }
        Expr::Func(_, args) => args.iter().any(|a| contains_algext(a.as_ref())),
        _ => false,
    }
}

/// Parse a concrete `AlgExt` value from `AlgExt` or rational `rootof`.
pub fn try_as_algext_data(e: &Expr) -> Option<AlgExtData> {
    match e {
        Expr::AlgExt(a) => Some((**a).clone()),
        Expr::Func(FuncKind::RootOf, args) if args.len() == 2 => {
            AlgExtData::from_rootof(&args[0], &args[1]).ok()
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q_minpoly() -> ExprArc {
        Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ))
    }

    #[test]
    fn algext_mul_squares_to_two() {
        let min = q_minpoly();
        let alpha = AlgExtData::from_rootof(&Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])), &min)
            .unwrap();
        let sq = alpha.mul(&alpha).unwrap();
        let two = AlgExtData {
            min_poly: alpha.min_poly.clone(),
            coords: vec![Expr::int(2)],
            root_index: None,
        };
        assert_eq!(sq.coords, two.coords);
    }

    #[test]
    fn algext_add_neg_cancels() {
        let min = q_minpoly();
        let pos = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap();
        let neg = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(-1), Expr::int(0)])),
            &min,
        )
        .unwrap();
        let z = pos.add(&neg).unwrap();
        assert!(z.is_zero());
    }

    #[test]
    fn algext_to_rootof_roundtrip_display() {
        let min = q_minpoly();
        let e = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap()
        .to_rootof_expr();
        let s = crate::format_expr(e.as_ref());
        assert!(s.contains("rootof"), "{s}");
        assert_eq!(s, "rootof([1,0],poly1[1,0,-2])");
    }

    #[test]
    fn algext_inv_divides_to_one() {
        let min = q_minpoly();
        let alpha = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap();
        let inv = alpha.inv().unwrap();
        let one = alpha.mul(&inv).unwrap();
        assert!(one.is_one());
    }

    #[test]
    fn algext_frac_via_eval() {
        use crate::{eval, Context};
        let ctx = Context::default();
        let min = q_minpoly();
        let alpha_data = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap();
        let alpha = alpha_data.clone().into_expr();
        let frac = Arc::new(Expr::Frac(Expr::int(2), alpha));
        let r = eval(frac.as_ref(), &ctx).unwrap();
        let quotient = match r.as_ref() {
            Expr::AlgExt(a) => (**a).clone(),
            other => panic!("expected AlgExt, got {other:?}"),
        };
        let prod = quotient.mul(&alpha_data).unwrap();
        let two = AlgExtData {
            min_poly: alpha_data.min_poly.clone(),
            coords: vec![Expr::int(2)],
            root_index: None,
        };
        assert!(prod.eq_mod(&two).unwrap());
    }
}
