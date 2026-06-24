//! GIAC-228: Rothstein–Trager algorithm for logarithmic part of ∫ N/Q dx.
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md) §5.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Partial** | `rothstein_trager_integrate`, `try_algebraic_rt_even_quartic`, `try_integrate_x4_plus_one` |
//! | **Pipeline private** | `ratio_to_expr` |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{bigint_to_i64, poly_to_expr, ratio_to_expr, EvalError, Expr, ExprArc, Ident};
use giac_poly::{
    eval_param_poly, num_minus_t_derivative, rational_roots_in_t, tresultant_eliminate_x,
    univariate_degree, Poly, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::integrate::ln_abs_expr;

use super::algebraic_rt::{is_monic_even_quartic, try_algebraic_rt_log_part};

const RT_PARAM: &str = "__rt";

/// **Partial** — algebraic RT for monic even quartics with constant numerator. **退役：** merge into full RT pipeline.
pub fn try_algebraic_rt_even_quartic(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Option<ExprArc> {
    if !is_monic_even_quartic(factor, var) || univariate_degree(numer, var) > 0 {
        return None;
    }
    let t = Var::from(RT_PARAM);
    let p1 = num_minus_t_derivative(numer, factor, var, &t);
    let res_t = tresultant_eliminate_x(&p1, factor, var, &t).ok()?;
    try_algebraic_rt_log_part(numer, factor, var, x, &res_t, &t)
}

/// **Partial** — `∫ k/(x^4+1) dx` via algebraic RT conjugate pairing. **退役：** general even-quartic RT.
pub fn try_integrate_x4_plus_one(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Option<ExprArc> {
    try_algebraic_rt_even_quartic(numer, factor, var, x)
}

/// **Partial** — integrate `numer / factor` when `factor` is square-free and partfrac failed. **退役：** full Rothstein–Trager with algebraic extensions.
pub fn rothstein_trager_integrate(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    if let Some(r) = try_algebraic_rt_even_quartic(numer, factor, var, x) {
        return Ok(r);
    }
    let t = Var::from(RT_PARAM);
    let p1 = num_minus_t_derivative(numer, factor, var, &t);
    let res_t = tresultant_eliminate_x(&p1, factor, var, &t)?;
    if res_t.is_zero() {
        return Err(EvalError::NotImplemented("rothstein trager"));
    }
    if let Some(r) = try_algebraic_rt_log_part(numer, factor, var, x, &res_t, &t) {
        return Ok(r);
    }
    let roots = rational_roots_in_t(&res_t, &t)?;
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


#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{format_expr, Expr, FuncKind};

    use crate::plugin::xcas_default;
    use crate::test_verify::assert_deriv_equals_integrand;

    use super::*;

    fn x_var() -> Var {
        Var::from("x")
    }

    fn x_id() -> Ident {
        Ident::new("x")
    }

    fn integrand_one_over_x4_plus_one() -> ExprArc {
        Arc::new(Expr::Frac(
            Expr::int(1),
            Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(4)),
                Expr::int(1),
            ]),
        ))
    }

    #[test]
    // smoke-until B-T3: delete when `rothstein_deriv_equals_integrand` green
    fn rothstein_one_over_x_fourth_plus_one_smoke() {
        let var = x_var();
        let den = Poly::var("x").pow(4).add(&Poly::one());
        let r = rothstein_trager_integrate(&Poly::one(), &den, &var, &x_id()).unwrap();
        assert!(!format_expr(r.as_ref()).is_empty());
    }

    #[test]
    #[ignore = "GIAC-expr-api T3: diff(antiderivative) not implemented for ln/atan shape"]
    fn rothstein_deriv_equals_integrand() {
        let ctx = xcas_default();
        let var = x_var();
        let den = Poly::var("x").pow(4).add(&Poly::one());
        let r = rothstein_trager_integrate(&Poly::one(), &den, &var, &x_id()).unwrap();
        assert_deriv_equals_integrand(
            &integrand_one_over_x4_plus_one(),
            &x_id(),
            &r,
            &ctx,
        );
    }

    #[test]
    fn rothstein_one_over_x_squared_plus_one() {
        let var = x_var();
        let den = Poly::var("x").pow(2).add(&Poly::one());
        let r = rothstein_trager_integrate(&Poly::one(), &den, &var, &x_id());
        assert!(matches!(r, Err(EvalError::NotImplemented("rothstein trager"))));
    }
}
