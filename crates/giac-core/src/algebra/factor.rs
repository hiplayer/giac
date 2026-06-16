use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::Signed;

use crate::{Context, EvalError, Expr, ExprArc};
use giac_poly::{factor_into, factor_poly};
use giac_poly::Poly;

use super::normal::normal;
use super::poly::{expr_to_poly, poly_to_expr};

/// Full factorization: structural (`Mul`/`Pow`/`Frac`) then polynomial engine.
pub fn factor(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    factor_expr(expr, ctx)
}

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
                if let Ok(e_u) = crate::num_util::bigint_to_nonneg_u32(n) {
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
        _ => factor_poly_form(e, ctx),
    }
}

fn flatten_mul(e: &Expr) -> Vec<ExprArc> {
    match e {
        Expr::Mul(fs) => fs.clone(),
        other => vec![Arc::new(other.clone())],
    }
}

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
        if factors.len() == 1 {
            return Ok(poly_to_expr(&factors[0]));
        }
    }
    Ok(poly_to_expr(&factor_poly(&p)))
}

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
