//! Expr → (numerator, denominator) as `Poly` over ℚ.
//!
//! **Tier:** **Stable (bounded)** — structural decomposition; eval variant for builtin adapters.

use num_bigint::BigInt;
use num_traits::Signed;

use giac_poly::Poly;

use crate::context::Context;
use crate::error::EvalError;
use crate::eval::eval;
use crate::expr::Expr;

use crate::algebra::poly::expr_to_poly;

/// Decompose an expression into `(num, den)` polynomials (no subexpression `eval`).
/// **Stable (bounded)** — `Pow`⁻¹ / `Frac` / leaf → poly
pub fn expr_to_rational_polys(e: &Expr) -> Result<(Poly, Poly), EvalError> {
    match e {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) => {
            Ok((Poly::one(), expr_to_poly(base)?))
        }
        Expr::Frac(n, d) => Ok((expr_to_poly(n.as_ref())?, expr_to_poly(d.as_ref())?)),
        other => Ok((expr_to_poly(other)?, Poly::one())),
    }
}

/// Like [`expr_to_rational_polys`] but evaluates leaves through `ctx` first.
/// **Stable (bounded)** — for `froot` / partially evaluated rationals
pub fn expr_to_rational_polys_eval(e: &Expr, ctx: &Context) -> Result<(Poly, Poly), EvalError> {
    match e {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) => Ok((
            Poly::one(),
            expr_to_poly(eval(base.as_ref(), ctx)?.as_ref())?,
        )),
        Expr::Frac(num, den) => Ok((
            expr_to_poly(eval(num.as_ref(), ctx)?.as_ref())?,
            expr_to_poly(eval(den.as_ref(), ctx)?.as_ref())?,
        )),
        Expr::Mul(factors) => {
            if factors.len() == 2 {
                if let Expr::Pow(base, exp) = factors[1].as_ref() {
                    if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) {
                        return Ok((
                            expr_to_poly(eval(factors[0].as_ref(), ctx)?.as_ref())?,
                            expr_to_poly(eval(base.as_ref(), ctx)?.as_ref())?,
                        ));
                    }
                }
            }
            Ok((expr_to_poly(eval(e, ctx)?.as_ref())?, Poly::one()))
        }
        other => Ok((expr_to_poly(eval(other, ctx)?.as_ref())?, Poly::one())),
    }
}
