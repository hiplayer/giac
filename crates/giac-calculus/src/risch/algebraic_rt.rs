//! GIAC-228a: algebraic Rothstein–Trager log part via conjugate root pairing.

use std::sync::Arc;

use giac_core::{bigint_to_i64, Expr, ExprArc, FuncKind, Ident};
use giac_poly::{
    biquartic_conjugate_pairs, coeff_at, univariate_degree, AlgebraicRt, ConjugatePair, Poly, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::integrate::{ln_abs_expr, var_to_expr};

struct LogFactor {
    pair: ConjugatePair,
    g: ExprArc,
    atan_offset: i64,
}

/// Integrate `k/(x^4+1)` via RT conjugate pairing.
pub fn integrate_monic_x4_plus_one(k: &Ratio<BigInt>, x: &Ident) -> ExprArc {
    let mut parts = Vec::new();
    for lf in x4_plus_one_log_factors(x) {
        parts.push(conjugate_pair_log_contribution(k, &lf.pair.alpha, &lf.g, lf.atan_offset, x));
    }
    Expr::add(parts)
}

/// RT log part when `Res_t` has algebraic conjugate pairs and `Q` is monic `x^4+1`.
pub fn try_algebraic_rt_log_part(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
    res_t: &Poly,
    t_var: &Var,
) -> Option<ExprArc> {
    if !is_monic_x4_plus_one(factor, var) {
        return None;
    }
    if univariate_degree(numer, var) > 0 {
        return None;
    }
    let k = coeff_at(numer, var, 0);
    if k.is_zero() {
        return Some(Expr::int(0));
    }
    let pairs = biquartic_conjugate_pairs(res_t, t_var).ok()?;
    if pairs.len() != 2 {
        return None;
    }
    let mut log_factors = x4_plus_one_log_factors(x);
    log_factors.sort_by_key(|lf| {
        if lf.pair.alpha.re > Ratio::zero() {
            0
        } else {
            1
        }
    });
    let mut sorted_pairs = pairs;
    sorted_pairs.sort_by_key(|p| {
        if p.alpha.re > Ratio::zero() {
            0
        } else {
            1
        }
    });
    let mut parts = Vec::new();
    for (pair, lf) in sorted_pairs.iter().zip(log_factors.iter()) {
        parts.push(conjugate_pair_log_contribution(
            &k,
            &pair.alpha,
            &lf.g,
            lf.atan_offset,
            x,
        ));
    }
    Some(Expr::add(parts))
}

pub fn is_monic_x4_plus_one(p: &Poly, var: &Var) -> bool {
    if univariate_degree(p, var) != 4 {
        return false;
    }
    coeff_at(p, var, 4) == Ratio::one()
        && coeff_at(p, var, 3).is_zero()
        && coeff_at(p, var, 2).is_zero()
        && coeff_at(p, var, 1).is_zero()
        && coeff_at(p, var, 0) == Ratio::one()
}

fn x4_plus_one_log_factors(x: &Ident) -> Vec<LogFactor> {
    let rt2 = Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]);
    let xv = var_to_expr(x);
    let q_plus = Expr::add(vec![
        Expr::pow(xv.clone(), Expr::int(2)),
        Expr::mul(vec![rt2.clone(), xv.clone()]),
        Expr::int(1),
    ]);
    let q_minus = Expr::add(vec![
        Expr::pow(xv.clone(), Expr::int(2)),
        Expr::mul(vec![Expr::int(-1), rt2, xv]),
        Expr::int(1),
    ]);
    vec![
        LogFactor {
            pair: ConjugatePair {
                alpha: AlgebraicRt {
                    re: Ratio::new(1.into(), 8.into()),
                    im: Ratio::new(1.into(), 8.into()),
                    ext: 2,
                },
            },
            g: q_plus,
            atan_offset: 1,
        },
        LogFactor {
            pair: ConjugatePair {
                alpha: AlgebraicRt {
                    re: Ratio::new((-1).into(), 8.into()),
                    im: Ratio::new(1.into(), 8.into()),
                    ext: 2,
                },
            },
            g: q_minus,
            atan_offset: -1,
        },
    ]
}

/// `Re(α)·ln|G| + 2·Im(α)·atan(√2·x + offset)` for `α = (re + im·i)·√ext`.
fn conjugate_pair_log_contribution(
    k: &Ratio<BigInt>,
    alpha: &AlgebraicRt,
    g: &ExprArc,
    atan_offset: i64,
    x: &Ident,
) -> ExprArc {
    let rt_ext = Expr::func(FuncKind::Sqrt, vec![Expr::int(alpha.ext as i64)]);
    let re_part = Expr::mul(vec![ratio_to_expr(&alpha.re), rt_ext.clone()]);
    let im_part = Expr::mul(vec![ratio_to_expr(&alpha.im), rt_ext]);
    let k_expr = ratio_to_expr(k);
    let ln_term = Expr::mul(vec![k_expr.clone(), re_part, ln_abs_expr(g.clone())]);
    let rt2 = Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]);
    let xv = var_to_expr(x);
    let atan_arg = Expr::add(vec![
        Expr::mul(vec![rt2, xv]),
        Expr::int(atan_offset),
    ]);
    let atan_term = Expr::mul(vec![
        k_expr,
        Expr::int(2),
        im_part,
        Expr::func(FuncKind::Atan, vec![atan_arg]),
    ]);
    Expr::add(vec![ln_term, atan_term])
}

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

    use super::*;

    #[test]
    fn algebraic_rt_one_over_x4_plus_one_shape() {
        let x = Ident::new("x");
        let r = integrate_monic_x4_plus_one(&Ratio::one(), &x);
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln"), "got {s}");
        assert!(s.contains("atan"), "got {s}");
    }
}
