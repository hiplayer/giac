//! GIAC-228a/228b: algebraic Rothstein–Trager log part via conjugate root pairing.
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md) §5.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Partial** | `integrate_monic_x4_plus_one`, `try_algebraic_rt_log_part`, `is_monic_x4_plus_one`, `is_monic_even_quartic` |
//! | **Pipeline private** | factorization, conjugate pairing, sqrt/ratio helpers |

use std::sync::Arc;

use giac_core::{bigint_to_i64, Expr, ExprArc, FuncKind, Ident};
use giac_poly::{
    biquadratic_res_conjugate_pairs, coeff_at, num_minus_t_derivative, tresultant_eliminate_x,
    univariate_degree, AlgebraicRt, Poly, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::integrate::{ln_abs_expr, var_to_expr};

struct QuadraticFactor {
    g: ExprArc,
    b_lin: ExprArc,
    c: Ratio<BigInt>,
    /// `4c - b²` for atan denominator (same for ±u factors).
    atan_disc: Ratio<BigInt>,
}

/// **Partial** — integrate `k/(x^4+1)` via RT conjugate pairing. **退役：** `try_algebraic_rt_log_part` general path.
pub fn integrate_monic_x4_plus_one(k: &Ratio<BigInt>, x: &Ident) -> ExprArc {
    let q = Poly::var("x").pow(4).add(&Poly::one());
    let var = Var::from("x");
    let factors = even_quartic_quadratic_factors(&q, &var, x).expect("x^4+1 factors");
    let pairs = default_x4_plus_one_pairs();
    integrate_from_pairs(k, &pairs, &factors, x)
}

/// **Partial** — RT log part when `Res_t` has algebraic conjugate pairs on an even monic quartic. **退役：** unified RT resultant handler.
pub fn try_algebraic_rt_log_part(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
    res_t: &Poly,
    t_var: &Var,
) -> Option<ExprArc> {
    if !is_monic_even_quartic(factor, var) {
        return None;
    }
    if univariate_degree(numer, var) > 0 {
        return None;
    }
    let k = coeff_at(numer, var, 0);
    if k.is_zero() {
        return Some(Expr::int(0));
    }
    // Conjugate data comes from `Res(1 - t·Q', Q)`; scale the primitive by `k` afterward.
    let unit_res = if k == Ratio::one() {
        res_t.clone()
    } else {
        let p1 = num_minus_t_derivative(&Poly::one(), factor, var, t_var);
        tresultant_eliminate_x(&p1, factor, var, t_var).ok()?
    };
    let mut pairs = biquadratic_res_conjugate_pairs(&unit_res, t_var).ok()?;
    if pairs.is_empty() {
        return None;
    }
    let mut factors = even_quartic_quadratic_factors(factor, var, x)?;
    if factors.len() != pairs.len() {
        return None;
    }
    sort_pairs_and_factors(&mut pairs, &mut factors);
    Some(integrate_from_pairs(&k, &pairs, &factors, x))
}

/// **Partial** — monic quartic `x^4+1` shape predicate.
pub fn is_monic_x4_plus_one(p: &Poly, var: &Var) -> bool {
    is_monic_even_quartic(p, var)
        && coeff_at(p, var, 2).is_zero()
        && coeff_at(p, var, 0) == Ratio::one()
}

/// **Partial** — monic even quartic shape predicate (odd coefficients zero).
pub fn is_monic_even_quartic(p: &Poly, var: &Var) -> bool {
    if univariate_degree(p, var) != 4 {
        return false;
    }
    coeff_at(p, var, 4) == Ratio::one()
        && coeff_at(p, var, 3).is_zero()
        && coeff_at(p, var, 1).is_zero()
}

// **Pipeline private** — factor monic even quartic into two quadratics.
fn even_quartic_quadratic_factors(q: &Poly, var: &Var, x: &Ident) -> Option<Vec<QuadraticFactor>> {
    if !is_monic_even_quartic(q, var) {
        return None;
    }
    let a2 = coeff_at(q, var, 2);
    let c0 = coeff_at(q, var, 0);
    let v = ratio_perfect_sqrt(&c0)?;
    let u_sq = Ratio::from_integer(BigInt::from(2)) * v.clone() - a2.clone();
    let u_lin = linear_sqrt_coeff_expr(&u_sq)?;
    let neg_u = negate_linear_expr(&u_lin);
    let atan_disc = Ratio::from_integer(BigInt::from(4)) * v.clone() - u_sq;
    let xv = var_to_expr(x);
    let q_plus = quadratic_expr(&xv, &u_lin, &v);
    let q_minus = quadratic_expr(&xv, &neg_u, &v);
    Some(vec![
        QuadraticFactor {
            g: q_plus,
            b_lin: u_lin,
            c: v.clone(),
            atan_disc: atan_disc.clone(),
        },
        QuadraticFactor {
            g: q_minus,
            b_lin: neg_u,
            c: v,
            atan_disc,
        },
    ])
}

// **Pipeline private** — linear coefficient from `√(u²)` rational or surd.
fn linear_sqrt_coeff_expr(u_sq: &Ratio<BigInt>) -> Option<ExprArc> {
    if let Some(r) = ratio_perfect_sqrt(u_sq) {
        return Some(ratio_to_expr(&r));
    }
    let (coeff, rad) = sqrt_rational_coeff_radicand(u_sq);
    if rad == 1 {
        return Some(ratio_to_expr(&coeff));
    }
    Some(Expr::mul(vec![
        ratio_to_expr(&coeff),
        Expr::func(FuncKind::Sqrt, vec![Expr::int(rad as i64)]),
    ]))
}

// **Pipeline private** — split rational into outer coeff and inner radicand.
fn sqrt_rational_coeff_radicand(r: &Ratio<BigInt>) -> (Ratio<BigInt>, u64) {
    if let Some(s) = ratio_perfect_sqrt(r) {
        return (s, 1);
    }
    let num = r.numer().abs();
    let den = r.denom().abs();
    let (out_num, sf_num) = extract_sqrt_factor_bigint(&num);
    let (out_den, sf_den) = extract_sqrt_factor_bigint(&den);
    let coeff = Ratio::new(out_num, out_den);
    let rad_int = sf_num * sf_den;
    let rad = rad_int.to_string().parse().unwrap_or(1);
    (coeff, rad)
}

// **Pipeline private** — extract perfect-square factor from integer.
fn extract_sqrt_factor_bigint(n: &BigInt) -> (BigInt, BigInt) {
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

// **Pipeline private** — negate linear expression coefficient.
fn negate_linear_expr(e: &ExprArc) -> ExprArc {
    match e.as_ref() {
        Expr::Rat(r) => ratio_to_expr(&(-r.clone())),
        Expr::Int(n) => bigint_to_i64(n)
            .map(Expr::int)
            .map(|e| Expr::mul(vec![Expr::int(-1), e]))
            .unwrap_or_else(|_| Expr::mul(vec![Expr::int(-1), e.clone()])),
        Expr::Mul(factors) => {
            if let Some(_first) = factors.first() {
                return Expr::mul(std::iter::once(Expr::int(-1))
                    .chain(factors.iter().cloned())
                    .collect());
            }
            Expr::int(-1)
        }
        _ => Expr::mul(vec![Expr::int(-1), e.clone()]),
    }
}

// **Pipeline private** — build `x² + b·x + c` expression.
fn quadratic_expr(x: &ExprArc, b_lin: &ExprArc, c: &Ratio<BigInt>) -> ExprArc {
    Expr::add(vec![
        Expr::pow(x.clone(), Expr::int(2)),
        Expr::mul(vec![b_lin.clone(), x.clone()]),
        ratio_to_expr(c),
    ])
}

// **Pipeline private** — align conjugate pairs with quadratic factors.
fn sort_pairs_and_factors(pairs: &mut [giac_poly::ConjugatePair], factors: &mut [QuadraticFactor]) {
    pairs.sort_by(|a, b| alpha_re_key(&a.alpha).cmp(&alpha_re_key(&b.alpha)));
    factors.sort_by(|a, b| b_lin_sign_key(&a.b_lin).cmp(&b_lin_sign_key(&b.b_lin)));
}

// **Pipeline private** — sort key from linear `b` coefficient sign.
fn b_lin_sign_key(b_lin: &ExprArc) -> (i8, String) {
    let sign = match b_lin.as_ref() {
        Expr::Rat(r) if r.is_negative() => 1,
        Expr::Int(n) if n.is_negative() => 1,
        Expr::Mul(fs) if fs.first().map(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative())).unwrap_or(false) => {
            1
        }
        _ => 0,
    };
    (sign, format!("{:?}", b_lin))
}

// **Pipeline private** — sort key from conjugate root real part.
fn alpha_re_key(a: &AlgebraicRt) -> (i8, Ratio<BigInt>, Ratio<BigInt>) {
    let sign = if a.re > Ratio::zero() || a.re_b > Ratio::zero() {
        0
    } else {
        1
    };
    (sign, a.re.clone(), a.re_b.clone())
}

// **Pipeline private** — sum log/atan contributions from conjugate pairs.
fn integrate_from_pairs(
    k: &Ratio<BigInt>,
    pairs: &[giac_poly::ConjugatePair],
    factors: &[QuadraticFactor],
    x: &Ident,
) -> ExprArc {
    let mut parts = Vec::new();
    for (pair, qf) in pairs.iter().zip(factors.iter()) {
        parts.push(conjugate_pair_log_contribution(k, &pair.alpha, qf, x));
    }
    Expr::add(parts)
}

// **Pipeline private** — hard-coded conjugate pairs for `x^4+1`.
fn default_x4_plus_one_pairs() -> Vec<giac_poly::ConjugatePair> {
    vec![
        giac_poly::ConjugatePair {
            alpha: AlgebraicRt {
                re: Ratio::zero(),
                im_a: Ratio::zero(),
                re_b: Ratio::new(1.into(), 8.into()),
                im_b: Ratio::new(1.into(), 8.into()),
                rad: 2,
            },
        },
        giac_poly::ConjugatePair {
            alpha: AlgebraicRt {
                re: Ratio::zero(),
                im_a: Ratio::zero(),
                re_b: Ratio::new((-1).into(), 8.into()),
                im_b: Ratio::new(1.into(), 8.into()),
                rad: 2,
            },
        },
    ]
}

// **Pipeline private** — `Re(α)·ln|G| + 2·Im(α)·atan((2x+b)/√(4c-b²))`.
fn conjugate_pair_log_contribution(
    k: &Ratio<BigInt>,
    alpha: &AlgebraicRt,
    qf: &QuadraticFactor,
    x: &Ident,
) -> ExprArc {
    let k_expr = ratio_to_expr(k);
    let re_part = alpha_re_expr(alpha);
    let im_part = alpha_im_expr(alpha);
    let ln_term = Expr::mul(vec![k_expr.clone(), re_part, ln_abs_expr(qf.g.clone())]);
    let disc = qf.atan_disc.clone();
    let sqrt_disc = sqrt_ratio_expr(&disc);
    let two_x = Expr::mul(vec![Expr::int(2), var_to_expr(x)]);
    let numer = Expr::add(vec![two_x, qf.b_lin.clone()]);
    let atan_arg = Arc::new(Expr::Frac(numer, sqrt_disc));
    let atan_term = Expr::mul(vec![
        k_expr,
        Expr::int(2),
        im_part,
        Expr::func(FuncKind::Atan, vec![atan_arg]),
    ]);
    Expr::add(vec![ln_term, atan_term])
}

// **Pipeline private** — real part of algebraic root as expression.
fn alpha_re_expr(a: &AlgebraicRt) -> ExprArc {
    let mut parts = Vec::new();
    if !a.re.is_zero() {
        parts.push(ratio_to_expr(&a.re));
    }
    if !a.re_b.is_zero() {
        parts.push(Expr::mul(vec![
            ratio_to_expr(&a.re_b),
            Expr::func(FuncKind::Sqrt, vec![Expr::int(a.rad as i64)]),
        ]));
    }
    if parts.is_empty() {
        Expr::int(0)
    } else {
        Expr::add(parts)
    }
}

// **Pipeline private** — imaginary part of algebraic root as expression.
fn alpha_im_expr(a: &AlgebraicRt) -> ExprArc {
    let mut parts = Vec::new();
    if !a.im_a.is_zero() {
        parts.push(ratio_to_expr(&a.im_a));
    }
    if !a.im_b.is_zero() {
        parts.push(Expr::mul(vec![
            ratio_to_expr(&a.im_b),
            Expr::func(FuncKind::Sqrt, vec![Expr::int(a.rad as i64)]),
        ]));
    }
    if parts.is_empty() {
        Expr::int(0)
    } else {
        Expr::add(parts)
    }
}

// **Pipeline private** — `√r` as expression (perfect square or surd ratio).
fn sqrt_ratio_expr(r: &Ratio<BigInt>) -> ExprArc {
    if r.is_zero() {
        return Expr::int(0);
    }
    if let Some(s) = ratio_perfect_sqrt(r) {
        return ratio_to_expr(&s);
    }
    let num = r.numer().abs();
    let den = r.denom().abs();
    let sqrt_num = Expr::func(
        FuncKind::Sqrt,
        vec![ratio_to_expr(&Ratio::new(num, BigInt::one()))],
    );
    if den == BigInt::one() {
        return sqrt_num;
    }
    Arc::new(Expr::Frac(
        sqrt_num,
        Expr::func(
            FuncKind::Sqrt,
            vec![ratio_to_expr(&Ratio::new(den, BigInt::one()))],
        ),
    ))
}

// **Pipeline private** — perfect rational square root if exists.
fn ratio_perfect_sqrt(r: &Ratio<BigInt>) -> Option<Ratio<BigInt>> {
    if r.is_zero() {
        return Some(Ratio::zero());
    }
    let sn = integer_perfect_sqrt(r.numer())?;
    let sd = integer_perfect_sqrt(r.denom())?;
    Some(Ratio::new(sn, sd))
}

// **Pipeline private** — integer perfect square root via binary search.
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

// **Pipeline private** — `Ratio<BigInt>` to `ExprArc`.
fn ratio_to_expr(r: &Ratio<BigInt>) -> ExprArc {
    if *r.denom() == BigInt::one() {
        bigint_to_i64(r.numer())
            .map(Expr::int)
            .unwrap_or_else(|_| Arc::new(Expr::Rat(r.clone())))
    } else {
        Arc::new(Expr::Rat(r.clone()))
    }
}

#[cfg(test)]
mod tests {
    use giac_core::format_expr;
    use giac_poly::{num_minus_t_derivative, tresultant_eliminate_x, Poly, Var};

    use super::*;

    #[test]
    fn algebraic_rt_one_over_x4_plus_one_shape() {
        let x = Ident::new("x");
        let r = integrate_monic_x4_plus_one(&Ratio::one(), &x);
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln"), "got {s}");
        assert!(s.contains("atan"), "got {s}");
    }

    #[test]
    fn algebraic_rt_one_over_x4_plus_four_via_res() {
        let x = Ident::new("x");
        let var = Var::from("x");
        let tv = Var::from("__rt");
        let den = Poly::var("x").pow(4).add(&Poly::constant(Ratio::from_integer(4.into())));
        let p1 = num_minus_t_derivative(&Poly::one(), &den, &var, &tv);
        let res_t = tresultant_eliminate_x(&p1, &den, &var, &tv).unwrap();
        let r = try_algebraic_rt_log_part(&Poly::one(), &den, &var, &x, &res_t, &tv).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln"), "got {s}");
        assert!(s.contains("atan"), "got {s}");
    }

    #[test]
    fn algebraic_rt_one_over_x4_plus_x2_plus_one_via_res() {
        let x = Ident::new("x");
        let var = Var::from("x");
        let tv = Var::from("__rt");
        let den = Poly::var("x")
            .pow(4)
            .add(&Poly::var("x").pow(2))
            .add(&Poly::one());
        let p1 = num_minus_t_derivative(&Poly::one(), &den, &var, &tv);
        let res_t = tresultant_eliminate_x(&p1, &den, &var, &tv).unwrap();
        let r = try_algebraic_rt_log_part(&Poly::one(), &den, &var, &x, &res_t, &tv).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln"), "got {s}");
        assert!(s.contains("atan"), "got {s}");
    }
}
