//! Parametric resultant Res_x(P(x,t), Q(x)) eliminating `x`, yielding a polynomial in `t`.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, resultant, univariate_degree};
use crate::univariate::univariate_derivative;

/// Build `num(x) - t * den'(x)` as a polynomial in `(x, t)`.
/// **Stable** — RT numerator derivative
pub fn num_minus_t_derivative(num: &Poly, den: &Poly, x: &Var, t: &Var) -> Poly {
    let dp = univariate_derivative(den, x);
    let mut p1 = embed_univariate_x(num, x);
    for (m, c) in &dp.terms {
        let xe = m.exp_of(x);
        let mut term = Poly::var(t.clone()).mul_scalar(c);
        if xe > 0 {
            term = term.mul(&Poly::var(x.clone()).pow(xe));
        }
        p1 = p1.sub(&term);
    }
    p1
}

// **Pipeline private** — `embed_univariate_x`
fn embed_univariate_x(p: &Poly, x: &Var) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        if m.iter().all(|(v, _)| v == x) || m.is_const() {
            out.terms.insert(m.clone(), c.clone());
        }
    }
    out
}

/// Resultant eliminating `x`; result is univariate in `t` (via interpolation).
/// **Stable** — eliminate x via t-resultant
pub fn tresultant_eliminate_x(a: &Poly, b: &Poly, x: &Var, t: &Var) -> PolyResult<Poly> {
    let da = univariate_degree(a, x);
    let db = univariate_degree(b, x);
    if da == 0 || db == 0 {
        return Err(PolyError::TypeError("degenerate tresultant"));
    }
    let max_deg = da.max(db) as usize + 1;
    let mut samples = Vec::new();
    for k in 1..=(max_deg + 1) {
        let alpha = Ratio::from_integer(BigInt::from(k as i64));
        let ae = eval_param_poly(a, t, &alpha, x);
        let val = resultant(&ae, b, x)?;
        let c = val
            .terms
            .get(&crate::monomial::Monomial::one())
            .cloned()
            .unwrap_or_else(|| coeff_at(&val, x, 0));
        samples.push((alpha, c));
    }
    Ok(lagrange_poly(&samples, t))
}

// **Pipeline private** — `lagrange_poly`
fn lagrange_poly(samples: &[(Ratio<BigInt>, Ratio<BigInt>)], t: &Var) -> Poly {
    let mut out = Poly::zero();
    let n = samples.len();
    for i in 0..n {
        let (xi, yi) = &samples[i];
        let mut basis = Poly::constant(yi.clone());
        for j in 0..n {
            if i == j {
                continue;
            }
            let (xj, _) = &samples[j];
            let linear = Poly::var(t.clone()).sub(&Poly::constant(xj.clone()));
            let scale = Ratio::one() / (xi.clone() - xj.clone());
            basis = basis.mul(&linear.mul_scalar(&scale));
        }
        out = out.add(&basis);
    }
    out
}

/// Evaluate `p(x,t)` at `t = alpha`, returning a univariate polynomial in `x`.
/// **Stable** — substitute parameter in Poly
pub fn eval_param_poly(p: &Poly, t: &Var, alpha: &Ratio<BigInt>, x: &Var) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        let te = m.exp_of(t);
        let mut scale = c.clone();
        for _ in 0..te {
            scale *= alpha.clone();
        }
        let xe = m.exp_of(x);
        let xp = if xe == 0 {
            Poly::one()
        } else {
            Poly::var(x.clone()).pow(xe)
        };
        out = out.add(&xp.mul_scalar(&scale));
    }
    out
}

/// Rational roots of a univariate polynomial in `t` (degree ≤ 4).
/// **Stable** — rational roots in parameter t
pub fn rational_roots_in_t(p: &Poly, t: &Var) -> PolyResult<Vec<Ratio<BigInt>>> {
    let d = univariate_degree(p, t);
    if d == 0 {
        return Ok(vec![]);
    }
    if d == 1 {
        let a = coeff_at(p, t, 1);
        let b = coeff_at(p, t, 0);
        if a.is_zero() {
            return Ok(vec![]);
        }
        return Ok(vec![-b / a]);
    }
    let candidates = rational_root_candidates(p, t);
    let mut roots = Vec::new();
    for r in candidates {
        if eval_univariate_at_t(p, t, &r).is_zero() && !roots.contains(&r) {
            roots.push(r);
        }
    }
    Ok(roots)
}

// **Pipeline private** — `eval_univariate_at_t`
fn eval_univariate_at_t(p: &Poly, t: &Var, val: &Ratio<BigInt>) -> Ratio<BigInt> {
    let mut sum = Ratio::zero();
    for (m, c) in &p.terms {
        if !m.iter().all(|(v, _)| v == t) && !m.is_const() {
            continue;
        }
        let e = m.exp_of(t);
        let mut pow = Ratio::one();
        for _ in 0..e {
            pow *= val.clone();
        }
        sum += c.clone() * pow;
    }
    sum
}

// **Pipeline private** — `rational_root_candidates`
fn rational_root_candidates(p: &Poly, t: &Var) -> Vec<Ratio<BigInt>> {
    let d = univariate_degree(p, t);
    if d == 0 {
        return vec![];
    }
    let lc = coeff_at(p, t, d);
    let ac = coeff_at(p, t, 0);
    let mut numer_divs = vec![BigInt::from(1)];
    let mut denom_divs = vec![BigInt::from(1)];
    push_divisors(ac.numer(), &mut numer_divs);
    push_divisors(ac.denom(), &mut denom_divs);
    push_divisors(lc.numer(), &mut denom_divs);
    push_divisors(lc.denom(), &mut numer_divs);
    let mut out = Vec::new();
    for p in &numer_divs {
        for q in &denom_divs {
            if q.is_zero() {
                continue;
            }
            for sign in [1i64, -1] {
                let r = Ratio::new(p.clone() * sign, q.clone());
                if !out.contains(&r) {
                    out.push(r);
                }
            }
        }
    }
    out.sort();
    out
}

// **Pipeline private** — `push_divisors`
fn push_divisors(n: &BigInt, out: &mut Vec<BigInt>) {
    if n.is_zero() {
        return;
    }
    let a = n.abs();
    let mut i = BigInt::one();
    while &i * &i <= a {
        if (&a % &i).is_zero() {
            if !out.contains(&i) {
                out.push(i.clone());
            }
            let q = &a / &i;
            if !out.contains(&q) {
                out.push(q);
            }
        }
        i += BigInt::one();
    }
}

/// Element `(re + re_b·√rad) + (im_a + im_b·√rad)·i` with rational coefficients.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlgebraicRt {
    pub re: Ratio<BigInt>,
    pub im_a: Ratio<BigInt>,
    pub re_b: Ratio<BigInt>,
    pub im_b: Ratio<BigInt>,
    pub rad: u64,
}

/// Conjugate pair with `im(α) > 0` (in the `im_a + im_b·√rad` sense).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConjugatePair {
    pub alpha: AlgebraicRt,
}

/// Conjugate-root pairs for even `Res_t` of degree 4 with no rational roots.
/// **Partial** — RT biquadratic resolvent
pub fn biquadratic_res_conjugate_pairs(p: &Poly, t: &Var) -> PolyResult<Vec<ConjugatePair>> {
    if univariate_degree(p, t) != 4 {
        return Ok(vec![]);
    }
    if !coeff_at(p, t, 3).is_zero() || !coeff_at(p, t, 1).is_zero() {
        return Ok(vec![]);
    }
    if !rational_roots_in_t(p, t)?.is_empty() {
        return Ok(vec![]);
    }
    let a4 = coeff_at(p, t, 4);
    let a2 = coeff_at(p, t, 2);
    let a0 = coeff_at(p, t, 0);
    if a4.is_zero() || a0.is_zero() {
        return Ok(vec![]);
    }
    if a0 == Ratio::one() {
        if let Some(k) = ratio_perfect_sqrt(&a4) {
            let m_sq = Ratio::from_integer(BigInt::from(2)) * k.clone() - a2.clone();
            if let Some(m) = ratio_perfect_sqrt(&m_sq) {
                return Ok(pairs_from_symmetric_res_factors(&k, &m));
            }
        }
    }
    if a2.is_zero() {
        return pure_biquartic_res_pairs(&a4, &a0);
    }
    Ok(vec![])
}

/// Legacy alias for pure `a·t⁴ + c` resultants.
/// **Partial** — RT biquartic resolvent
pub fn biquartic_conjugate_pairs(p: &Poly, t: &Var) -> PolyResult<Vec<ConjugatePair>> {
    biquadratic_res_conjugate_pairs(p, t)
}

// **Pipeline private** — `pairs_from_symmetric_res_factors`
fn pairs_from_symmetric_res_factors(k: &Ratio<BigInt>, m: &Ratio<BigInt>) -> Vec<ConjugatePair> {
    let disc = m.clone() * m.clone() - Ratio::from_integer(BigInt::from(4)) * k.clone();
    if disc >= Ratio::zero() {
        return vec![];
    }
    let neg_disc = -disc;
    let two_k = Ratio::from_integer(BigInt::from(2)) * k.clone();
    let re_pos = m.clone() / two_k.clone();
    let (sqrt_c, rad) = sqrt_rational_coeff_radicand(&neg_disc);
    let imag = sqrt_c / two_k;
    let mk = |re: Ratio<BigInt>| {
        let (im_a, im_b, rad) = if rad == 1 {
            (imag.clone(), Ratio::zero(), 1)
        } else {
            (Ratio::zero(), imag.clone(), rad)
        };
        ConjugatePair {
            alpha: AlgebraicRt {
                re,
                im_a,
                re_b: Ratio::zero(),
                im_b,
                rad,
            },
        }
    };
    vec![mk(re_pos.clone()), mk(-re_pos)]
}

/// Write `√r = coeff · √rad` with squarefree `rad`.
// **Pipeline private** — `sqrt_rational_coeff_radicand`
fn sqrt_rational_coeff_radicand(r: &Ratio<BigInt>) -> (Ratio<BigInt>, u64) {
    if let Some(s) = ratio_perfect_sqrt(r) {
        return (s, 1);
    }
    let num = r.numer().abs();
    let den = r.denom().abs();
    let (out_num, sf_num) = extract_sqrt_factor(&num);
    let (out_den, sf_den) = extract_sqrt_factor(&den);
    let coeff = Ratio::new(out_num, out_den);
    let rad_int = sf_num * sf_den;
    let rad = rad_int.to_string().parse().unwrap_or(1);
    (coeff, rad)
}

// **Pipeline private** — `extract_sqrt_factor`
fn extract_sqrt_factor(n: &BigInt) -> (BigInt, BigInt) {
    let mut outer = BigInt::one();
    let mut inner = n.abs();
    let mut p = BigInt::from(2);
    while &p * &p <= inner {
        if (&inner % &p).is_zero() {
            let mut count = 0u32;
            while (&inner % &p).is_zero() {
                inner /= &p;
                count += 1;
            }
            if count % 2 == 1 {
                inner *= &p;
            }
            if count >= 2 {
                outer *= p.pow(count / 2);
            }
        }
        p += BigInt::one();
    }
    (outer, inner)
}

// **Pipeline private** — `pure_biquartic_res_pairs`
fn pure_biquartic_res_pairs(a: &Ratio<BigInt>, c: &Ratio<BigInt>) -> PolyResult<Vec<ConjugatePair>> {
    let t4 = -c.clone() / a.clone();
    if t4 >= Ratio::zero() {
        return Ok(vec![]);
    }
    let rad = t4.abs();
    let mag = ratio_perfect_sqrt(&rad).ok_or(PolyError::NotImplemented(
        "biquartic algebraic roots",
    ))?;
    let quarter_mag = mag.clone() / Ratio::from_integer(BigInt::from(4));
    let re_im = ratio_perfect_sqrt(&quarter_mag).ok_or(PolyError::NotImplemented(
        "biquartic algebraic roots",
    ))?;
    let mut pairs = Vec::new();
    for re_sign in [1i64, -1i64] {
        let re_b = re_im.clone() * Ratio::from_integer(BigInt::from(re_sign));
        if re_b.is_zero() {
            continue;
        }
        pairs.push(ConjugatePair {
            alpha: AlgebraicRt {
                re: Ratio::zero(),
                im_a: Ratio::zero(),
                re_b: re_b.clone(),
                im_b: re_im.clone(),
                rad: 2,
            },
        });
    }
    Ok(pairs)
}

// **Stable** — detect perfect square Ratio
fn ratio_perfect_sqrt(r: &Ratio<BigInt>) -> Option<Ratio<BigInt>> {
    if r.is_zero() {
        return Some(Ratio::zero());
    }
    let sn = integer_perfect_sqrt(r.numer())?;
    let sd = integer_perfect_sqrt(r.denom())?;
    Some(Ratio::new(sn, sd))
}

// **Pipeline private** — `integer_perfect_sqrt`
fn integer_perfect_sqrt(n: &BigInt) -> Option<BigInt> {
    if n.is_negative() {
        return None;
    }
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    let mut lo = BigInt::zero();
    let mut hi = n.clone() + BigInt::one();
    while lo < hi {
        let mid = (&lo + &hi) / BigInt::from(2);
        let sq = &mid * &mid;
        match sq.cmp(n) {
            std::cmp::Ordering::Equal => return Some(mid),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn x() -> Var {
        Var::from("x")
    }

    fn t() -> Var {
        Var::from("__rt")
    }

    #[test]
    fn biquartic_pairs_one_over_x_fourth_plus_one() {
        let tv = t();
        let r = Poly::var("__rt")
            .pow(4)
            .mul_scalar(&Ratio::from_integer(256.into()))
            .add(&Poly::one());
        let pairs = biquartic_conjugate_pairs(&r, &tv).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].alpha.re_b, Ratio::new(1.into(), 8.into()));
        assert_eq!(pairs[0].alpha.im_b, Ratio::new(1.into(), 8.into()));
        assert_eq!(pairs[0].alpha.rad, 2);
        assert_eq!(pairs[1].alpha.re_b, Ratio::new((-1).into(), 8.into()));
    }

    #[test]
    fn biquadratic_pairs_x_fourth_plus_four_res() {
        let tv = t();
        let r = Poly::var("__rt")
            .pow(4)
            .mul_scalar(&Ratio::from_integer(16384.into()))
            .add(&Poly::one());
        let pairs = biquadratic_res_conjugate_pairs(&r, &tv).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].alpha.re, Ratio::new(1.into(), 16.into()));
        assert_eq!(pairs[0].alpha.im_a, Ratio::new(1.into(), 16.into()));
    }

    #[test]
    fn biquadratic_pairs_x_fourth_plus_x_squared_plus_one_res() {
        let tv = t();
        let r = Poly::var("__rt")
            .pow(4)
            .mul_scalar(&Ratio::from_integer(144.into()))
            .add(&Poly::var("__rt").pow(2).mul_scalar(&Ratio::from_integer((-12).into())))
            .add(&Poly::one());
        let pairs = biquadratic_res_conjugate_pairs(&r, &tv).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].alpha.re, Ratio::new(1.into(), 4.into()));
        assert_eq!(pairs[0].alpha.im_b, Ratio::new(1.into(), 12.into()));
        assert_eq!(pairs[0].alpha.rad, 3);
    }

    #[test]
    fn tresultant_one_over_x_squared_plus_one() {
        let xv = x();
        let tv = t();
        let den = Poly::var("x").pow(2).add(&Poly::one());
        let num = Poly::one();
        let p1 = num_minus_t_derivative(&num, &den, &xv, &tv);
        let r = tresultant_eliminate_x(&p1, &den, &xv, &tv).unwrap();
        let roots = rational_roots_in_t(&r, &tv).unwrap();
        assert!(!roots.is_empty());
    }

    #[test]
    fn tresultant_one_over_x_fourth_plus_one() {
        let xv = x();
        let tv = t();
        let den = Poly::var("x").pow(4).add(&Poly::one());
        let num = Poly::one();
        let p1 = num_minus_t_derivative(&num, &den, &xv, &tv);
        let r = tresultant_eliminate_x(&p1, &den, &xv, &tv).unwrap();
        assert_eq!(univariate_degree(&r, &tv), 4);
        assert!(rational_roots_in_t(&r, &tv).unwrap().is_empty());
    }
}
