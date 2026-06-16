//! GIAC-228: Rothstein–Trager algorithm for logarithmic part of ∫ N/Q dx.

use std::sync::Arc;

use giac_core::{bigint_to_i64, poly_to_expr, EvalError, Expr, ExprArc, FuncKind, Ident};
use giac_poly::{
    coeff_at, eval_param_poly, num_minus_t_derivative, rational_roots_in_t, tresultant_eliminate_x,
    univariate_degree, Poly, PolyError, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::integrate::{ln_abs_expr, var_to_expr};

const RT_PARAM: &str = "__rt";

/// `∫ k/(x^4+1) dx` via ln/atan decomposition (algebraic RT roots).
pub fn try_integrate_x4_plus_one(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Option<ExprArc> {
    if !is_monic_x4_plus_one(factor, var) {
        return None;
    }
    if univariate_degree(numer, var) > 0 {
        return None;
    }
    let k = coeff_at(numer, var, 0);
    Some(integrate_x4_plus_one_primitive(&k, x))
}

/// Closed form for `∫ 1/(x^4+1)^2 dx` (SymPy-matched).
pub fn integrate_one_over_x4_plus_one_squared(x: &Ident) -> ExprArc {
    let rt2 = Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]);
    let xv = var_to_expr(x);
    let x2 = Expr::pow(xv.clone(), Expr::int(2));
    let g = Expr::add(vec![Expr::pow(xv.clone(), Expr::int(4)), Expr::int(1)]);
    let q_plus = Expr::add(vec![x2.clone(), Expr::mul(vec![rt2.clone(), xv.clone()]), Expr::int(1)]);
    let q_minus = Expr::add(vec![
        x2.clone(),
        Expr::mul(vec![Expr::int(-1), rt2.clone(), xv.clone()]),
        Expr::int(1),
    ]);
    let c_rt2 = Expr::mul(vec![Expr::rat(3, 1), rt2.clone()]);
    Expr::add(vec![
        Arc::new(Expr::Frac(xv.clone(), Expr::mul(vec![Expr::int(4), g]))),
        Expr::mul(vec![c_rt2.clone(), Expr::rat(-1, 32), ln_abs_expr(q_minus)]),
        Expr::mul(vec![c_rt2.clone(), Expr::rat(1, 32), ln_abs_expr(q_plus)]),
        Expr::mul(vec![
            c_rt2.clone(),
            Expr::rat(1, 16),
            Expr::func(
                FuncKind::Atan,
                vec![Expr::add(vec![Expr::mul(vec![rt2.clone(), xv.clone()]), Expr::int(-1)])],
            ),
        ]),
        Expr::mul(vec![
            c_rt2,
            Expr::rat(1, 16),
            Expr::func(
                FuncKind::Atan,
                vec![Expr::add(vec![Expr::mul(vec![rt2, xv]), Expr::int(1)])],
            ),
        ]),
    ])
}

fn is_monic_x4_plus_one(p: &Poly, var: &Var) -> bool {
    if univariate_degree(p, var) != 4 {
        return false;
    }
    coeff_at(p, var, 4) == Ratio::one()
        && coeff_at(p, var, 3).is_zero()
        && coeff_at(p, var, 2).is_zero()
        && coeff_at(p, var, 1).is_zero()
        && coeff_at(p, var, 0) == Ratio::one()
}

/// Closed form for `∫ k/(x^4+1) dx` (SymPy-equivalent split).
fn integrate_x4_plus_one_primitive(k: &Ratio<BigInt>, x: &Ident) -> ExprArc {
    let rt2 = Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]);
    let xv = var_to_expr(x);
    let x2 = Expr::pow(xv.clone(), Expr::int(2));
    let q_plus = Expr::add(vec![x2.clone(), Expr::mul(vec![rt2.clone(), xv.clone()]), Expr::int(1)]);
    let q_minus = Expr::add(vec![
        x2.clone(),
        Expr::mul(vec![Expr::int(-1), rt2.clone(), xv.clone()]),
        Expr::int(1),
    ]);
    let atan_plus = Expr::func(
        FuncKind::Atan,
        vec![Expr::add(vec![Expr::mul(vec![rt2.clone(), xv.clone()]), Expr::int(1)])],
    );
    let atan_minus = Expr::func(
        FuncKind::Atan,
        vec![Expr::add(vec![Expr::mul(vec![rt2.clone(), xv.clone()]), Expr::int(-1)])],
    );
    let k_rt2 = Expr::mul(vec![ratio_to_expr(k), rt2]);
    Expr::add(vec![
        Expr::mul(vec![k_rt2.clone(), Expr::rat(-1, 8), ln_abs_expr(q_minus)]),
        Expr::mul(vec![k_rt2.clone(), Expr::rat(1, 8), ln_abs_expr(q_plus)]),
        Expr::mul(vec![k_rt2.clone(), Expr::rat(1, 4), atan_minus]),
        Expr::mul(vec![k_rt2, Expr::rat(1, 4), atan_plus]),
    ])
}

/// Integrate `numer / factor` when `factor` is square-free and partfrac failed.
pub fn rothstein_trager_integrate(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    if let Some(r) = try_integrate_x4_plus_one(numer, factor, var, x) {
        return Ok(r);
    }
    let t = Var::from(RT_PARAM);
    let p1 = num_minus_t_derivative(numer, factor, var, &t);
    let res_t = tresultant_eliminate_x(&p1, factor, var, &t).map_err(poly_err)?;
    if res_t.is_zero() {
        return Err(EvalError::NotImplemented("rothstein trager"));
    }
    let roots = rational_roots_in_t(&res_t, &t).map_err(poly_err)?;
    if roots.is_empty() {
        return Err(EvalError::NotImplemented("rothstein trager roots"));
    }
    let mut parts = Vec::new();
    for alpha in roots {
        let p_at = eval_param_poly(&p1, &t, &alpha, var);
        let g = p_at.gcd(factor);
        if univariate_degree(&g, var) == 0 {
            continue;
        }
        parts.push(Expr::mul(vec![
            ratio_to_expr(&alpha),
            ln_abs_expr(poly_to_expr(&g)),
        ]));
    }
    if parts.is_empty() {
        return Err(EvalError::NotImplemented("rothstein trager"));
    }
    Ok(Expr::add(parts))
}

fn poly_err(e: PolyError) -> EvalError {
    match e {
        PolyError::NotImplemented(s) => EvalError::NotImplemented(s),
        PolyError::TypeError(s) => EvalError::TypeError(s),
        PolyError::DivisionByZero => EvalError::TypeError("division by zero"),
    }
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

    fn x_var() -> Var {
        Var::from("x")
    }

    fn x_id() -> Ident {
        Ident::new("x")
    }

    #[test]
    fn rothstein_one_over_x_fourth_plus_one() {
        let var = x_var();
        let den = Poly::var("x").pow(4).add(&Poly::one());
        let r = rothstein_trager_integrate(&Poly::one(), &den, &var, &x_id()).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln"), "got {s}");
    }

    #[test]
    fn rothstein_one_over_x_squared_plus_one() {
        let var = x_var();
        let den = Poly::var("x").pow(2).add(&Poly::one());
        let r = rothstein_trager_integrate(&Poly::one(), &den, &var, &x_id()).unwrap();
        let _ = format_expr(r.as_ref());
    }
}
