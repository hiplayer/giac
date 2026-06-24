//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-simplify-api-stability.md`.
//!
//!
use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use giac_core::{
    algext_poly_to_expr, expr_contains_alg_coeff, expr_to_poly, expr_to_rational_polys,
    factor_into_algext, factor_into_via_algext, poly_alg_from_expr, poly_to_expr, ratio_to_expr,
    univariate_poly_to_poly1_expr, AlgExtData, Context, EvalError, Expr, ExprArc, FuncKind,
};
use giac_poly::{factor_into, factor_poly, quadratic_abc, ratio_perfect_sqrt, vars_in, Poly};

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
    if let Some(factors) = try_factor_via_algext(&p)? {
        return Ok(Expr::mul(factors));
    }
    if let Some(factors) = try_factor_quadratic_rootof(&p) {
        return Ok(Expr::mul(factors));
    }
    if let Some(factors) = factor_into(&p) {
        if factors.len() == 1 {
            return Ok(poly_to_expr(&factors[0]));
        }
    }
    if ctx.with_sqrt {
        if let Some(factors) = try_factor_quadratic_sqrt(&p) {
            return Ok(Expr::mul(factors));
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

// **Pipeline private** — lift ℚ[x] to K[x] and factor when split exists (T2-2).
fn try_factor_via_algext(p: &Poly) -> Result<Option<Vec<ExprArc>>, EvalError> {
    let Some(factors) = factor_into_via_algext(p)? else {
        return Ok(None);
    };
    if factors.len() <= 1 {
        return Ok(None);
    }
    factors
        .iter()
        .map(algext_poly_to_expr)
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

// **Temporary** — Partial internal: quadratic → rootof when discriminant non-square.
fn try_factor_quadratic_rootof(p: &Poly) -> Option<Vec<ExprArc>> {
    let vars = vars_in(p);
    if vars.len() != 1 {
        return None;
    }
    let var = &vars[0];
    let (a, b, c) = quadratic_abc(p, var)?;
    let disc = &b * &b - Ratio::from_integer(BigInt::from(4)) * &a * &c;
    if disc.is_zero() || ratio_perfect_sqrt(&disc).is_some() {
        return None;
    }
    let minpoly = univariate_poly_to_poly1_expr(p, var);
    let pos = giac_core::AlgExtData::from_rootof(
        &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
        &minpoly,
    )
    .ok()?
    .into_expr();
    let neg = giac_core::AlgExtData::from_rootof(
        &Arc::new(Expr::Seq(vec![Expr::int(-1), Expr::int(0)])),
        &minpoly,
    )
    .ok()?
    .into_expr();
    let x = poly_to_expr(&Poly::var(var.clone()));
    let mut factors = vec![
        Expr::add(vec![x.clone(), Expr::mul(vec![Expr::int(-1), pos])]),
        Expr::add(vec![x, Expr::mul(vec![Expr::int(-1), neg])]),
    ];
    if !a.is_one() {
        factors.insert(0, poly_to_expr(&Poly::constant(a)));
    }
    Some(factors)
}

// **Temporary** — Partial internal: `ctx.with_sqrt` quadratic sqrt factors.
fn try_factor_quadratic_sqrt(p: &Poly) -> Option<Vec<ExprArc>> {
    let vars = vars_in(p);
    if vars.len() != 1 {
        return None;
    }
    let var = &vars[0];
    let (a, b, c) = quadratic_abc(p, var)?;
    let disc = &b * &b - Ratio::from_integer(BigInt::from(4)) * &a * &c;
    if disc.is_zero() {
        return None;
    }
    if ratio_perfect_sqrt(&disc).is_some() {
        return None;
    }
    let sqrt_d = Expr::func(FuncKind::Sqrt, vec![ratio_to_expr(&disc)]);
    let two_a = ratio_to_expr(&(Ratio::from_integer(BigInt::from(2)) * &a));
    let neg_b = ratio_to_expr(&(-&b));
    let r1 = Arc::new(Expr::Frac(
        Expr::add(vec![neg_b.clone(), sqrt_d.clone()]),
        two_a.clone(),
    ));
    let r2 = Arc::new(Expr::Frac(
        Expr::add(vec![neg_b, Expr::mul(vec![Expr::int(-1), sqrt_d])]),
        two_a,
    ));
    let x = poly_to_expr(&Poly::var(var.clone()));
    Some(vec![
        Expr::add(vec![x.clone(), Expr::mul(vec![Expr::int(-1), r1])]),
        Expr::add(vec![x, Expr::mul(vec![Expr::int(-1), r2])]),
    ])
}
