//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-simplify-api-stability.md`.
//!
//!
use std::sync::Arc;

use num_traits::Signed;

use giac_core::{
    algext_poly_to_expr, expr_contains_alg_coeff, expr_to_poly, expr_to_rational_polys,
    factor_into_algext, factor_into_via_algext, poly_alg_from_expr, poly_to_expr,
    Context, EvalError, Expr, ExprArc,
};
use giac_poly::{factor_into, factor_poly};

use crate::expand::normal;
use crate::ifactor::ifactor;

/// **Stable (bounded)** — structural (`Mul`/`Pow`/`Frac`) then polynomial factorization.
///
/// Multivariate Hensel / sparse fallback gaps are in `giac-poly` (FAC-G1–G3).
pub fn factor(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    factor_expr(expr, ctx)
}

// **Pipeline private** — recursive factor on Mul/Pow/Frac
fn factor_expr(e: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    match e {
        Expr::Mul(factors) => {
            let mut out = Vec::new();
            for f in factors {
                let ff = factor_expr(f.as_ref(), ctx)?;
                out.extend(flatten_mul(ff.as_ref()));
            }
            Ok(Expr::mul(out))
        }
        Expr::Pow(base, exp) => {
            if let Expr::Int(n) = exp.as_ref() {
                if n.is_negative() {
                    let b = factor_expr(base.as_ref(), ctx)?;
                    return Ok(Expr::pow(b, Arc::clone(exp)));
                }
                if let Ok(e_u) = giac_core::bigint_to_nonneg_u32(n) {
                    if e_u == 0 {
                        return Ok(Expr::int(1));
                    }
                    let fb = factor_expr(base.as_ref(), ctx)?;
                    let mut out = Vec::new();
                    for _ in 0..e_u {
                        out.extend(flatten_mul(fb.as_ref()));
                    }
                    return Ok(Expr::mul(out));
                }
            }
            factor_poly_form(e, ctx)
        }
        Expr::Frac(num, den) => {
            let fn_ = factor_expr(num.as_ref(), ctx)?;
            let fd = factor_expr(den.as_ref(), ctx)?;
            Ok(Expr::mul(vec![
                fn_,
                Expr::pow(fd, Expr::int(-1)),
            ]))
        }
        Expr::Int(n) => Ok(ifactor(n)),
        _ => factor_poly_form(e, ctx),
    }
}

// **Pipeline private** — flatten Mul to factor vec
fn flatten_mul(e: &Expr) -> Vec<ExprArc> {
    match e {
        Expr::Mul(fs) => fs.clone(),
        other => vec![Arc::new(other.clone())],
    }
}

// **Pipeline private** — normal→poly→factor_into chain
fn factor_poly_form(e: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    if let Ok((num, den)) = expr_to_rational_polys(e) {
        if !den.is_one() {
            let fn_ = factor_poly_form(&poly_to_expr(&num), ctx)?;
            let fd = factor_poly_form(&poly_to_expr(&den), ctx)?;
            return Ok(Expr::mul(vec![
                fn_,
                Expr::pow(fd, Expr::int(-1)),
            ]));
        }
    }
    let n = normal(e, ctx)?;
    if expr_contains_alg_coeff(n.as_ref()) {
        return factor_algext_form(n.as_ref());
    }
    let p = expr_to_poly(n.as_ref())?;
    if let Some(factors) = factor_into(&p) {
        if factors.len() > 1 {
            return Ok(Expr::mul(
                factors.into_iter().map(|f| poly_to_expr(&f)).collect(),
            ));
        }
    }
    if let Some(factors) = factor_into(&p) {
        if factors.len() == 1 {
            return Ok(poly_to_expr(&factors[0]));
        }
    }
    if ctx.with_sqrt {
        if let Some(factors) = factor_into_via_algext(&p)? {
            if factors.len() > 1 {
                return Ok(Expr::mul(
                    factors
                        .iter()
                        .map(algext_poly_to_expr)
                        .collect::<Result<Vec<_>, _>>()?,
                ));
            }
        }
    }
    Ok(poly_to_expr(&factor_poly(&p)))
}

// **Pipeline private** — path B: factor in K[var] when expr has algebraic coefficients.
fn factor_algext_form(e: &Expr) -> Result<ExprArc, EvalError> {
    let p = poly_alg_from_expr(e)?;
    if let Some(factors) = factor_into_algext(&p)? {
        if factors.len() > 1 {
            return Ok(Expr::mul(
                factors
                    .iter()
                    .map(algext_poly_to_expr)
                    .collect::<Result<Vec<_>, _>>()?,
            ));
        }
    }
    algext_poly_to_expr(&p)
}
