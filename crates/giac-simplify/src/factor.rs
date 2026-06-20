//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-simplify-api-stability.md`.
//!
//!
use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use giac_core::{Context, EvalError, Expr, ExprArc, FuncKind};
use giac_poly::{factor_into, factor_poly, ratio_perfect_sqrt, univariate_degree, vars_in, coeff_at};
use giac_poly::Poly;

use crate::ifactor::ifactor;
use crate::expand::normal;
use giac_core::{expr_to_poly, poly_to_expr, univariate_poly_to_poly1_expr, ratio_to_expr};

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
    if let Ok((num, den)) = rational_num_den(e) {
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
    let p = expr_to_poly(n.as_ref())?;
    if let Some(factors) = factor_into(&p) {
        if factors.len() > 1 {
            return Ok(Expr::mul(
                factors.into_iter().map(|f| poly_to_expr(&f)).collect(),
            ));
        }
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

// **Temporary** — Partial internal: quadratic → rootof when discriminant non-square.
fn try_factor_quadratic_rootof(p: &Poly) -> Option<Vec<ExprArc>> {
    let vars = vars_in(p);
    if vars.len() != 1 {
        return None;
    }
    let var = &vars[0];
    if univariate_degree(p, var) != 2 {
        return None;
    }
    let mut a = Ratio::zero();
    let mut b = Ratio::zero();
    let mut c = Ratio::zero();
    for (m, coeff) in &p.terms {
        match m.exp_of(var) {
            2 => a += coeff,
            1 => b += coeff,
            0 => c += coeff,
            _ => return None,
        }
    }
    if a.is_zero() {
        return None;
    }
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
    if univariate_degree(p, var) != 2 {
        return None;
    }
    let mut a = Ratio::zero();
    let mut b = Ratio::zero();
    let mut c = Ratio::zero();
    for (m, coeff) in &p.terms {
        match m.exp_of(var) {
            2 => a += coeff,
            1 => b += coeff,
            0 => c += coeff,
            _ => return None,
        }
    }
    if a.is_zero() {
        return None;
    }
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

// **Pipeline private** — Expr leaf to (num,den) Poly
fn rational_num_den(e: &Expr) -> Result<(Poly, Poly), EvalError> {
    match e {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) => {
            Ok((Poly::one(), expr_to_poly(base.as_ref())?))
        }
        Expr::Frac(n, d) => Ok((
            expr_to_poly(n.as_ref())?,
            expr_to_poly(d.as_ref())?,
        )),
        other => Ok((expr_to_poly(other)?, Poly::one())),
    }
}
