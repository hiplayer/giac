use crate::{Context, EvalError, Expr, ExprArc};
use giac_poly::{as_perfect_power, factor_into, factor_poly};

use super::normal::normal;
use super::poly::{expr_to_poly, poly_to_expr};

/// Basic polynomial factorization (perfect powers, x^n-1 for small n).
pub fn factor(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    let n = normal(expr, ctx)?;
    let p = expr_to_poly(n.as_ref())?;
    if let Some((base, exp)) = as_perfect_power(&p) {
        return Ok(Expr::pow(poly_to_expr(&base), Expr::int(exp as i64)));
    }
    if let Some(factors) = factor_into(&p) {
        if factors.len() > 1 {
            return Ok(Expr::mul(
                factors.into_iter().map(|f| poly_to_expr(&f)).collect(),
            ));
        }
    }
    Ok(poly_to_expr(&factor_poly(&p)))
}
