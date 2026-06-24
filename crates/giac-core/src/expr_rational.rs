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
/// **Stable (bounded)** — `Pow`⁻¹ / `Frac` / `num·base⁻¹` / leaf → poly
pub fn expr_to_rational_polys(e: &Expr) -> Result<(Poly, Poly), EvalError> {
    match e {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) => {
            Ok((Poly::one(), expr_to_poly(base)?))
        }
        Expr::Frac(n, d) => Ok((expr_to_poly(n.as_ref())?, expr_to_poly(d.as_ref())?)),
        Expr::Mul(factors) => {
            if factors.len() == 2 {
                if let Some((num, den)) = mul_inverse_rational_polys(factors[0].as_ref(), factors[1].as_ref())? {
                    return Ok((num, den));
                }
            }
            Ok((expr_to_poly(e)?, Poly::one()))
        }
        other => Ok((expr_to_poly(other)?, Poly::one())),
    }
}

// **Pipeline private** — `a * b^(-1)` or `b^(-1) * a` with poly-shaped factors.
fn mul_inverse_rational_polys(a: &Expr, b: &Expr) -> Result<Option<(Poly, Poly)>, EvalError> {
    if let Expr::Pow(base, exp) = b {
        if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
            return Ok(Some((expr_to_poly(a)?, expr_to_poly(base)?)));
        }
    }
    if let Expr::Pow(base, exp) = a {
        if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
            return Ok(Some((expr_to_poly(b)?, expr_to_poly(base)?)));
        }
    }
    Ok(None)
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

#[cfg(test)]
mod tests {
    use super::*;
    use giac_poly::Poly;
    use num_bigint::BigInt;
    use num_rational::Ratio;

    // **B** — structural `num * den^(-1)` decomposition.
    #[test]
    fn rational_polys_x_over_x_squared_minus_two() {
        let e = Expr::mul(vec![
            Expr::sym("x"),
            Expr::pow(
                Expr::add(vec![
                    Expr::pow(Expr::sym("x"), Expr::int(2)),
                    Expr::int(-2),
                ]),
                Expr::int(-1),
            ),
        ]);
        let (num, den) = expr_to_rational_polys(&e).unwrap();
        assert_eq!(giac_poly::univariate_degree(&num, &giac_poly::Var::from("x")), 1);
        assert_eq!(giac_poly::univariate_degree(&den, &giac_poly::Var::from("x")), 2);
        assert_eq!(
            coeff_at_wrt(&den, "x", 0),
            Ratio::from_integer(BigInt::from(-2))
        );
    }

    fn coeff_at_wrt(p: &Poly, var: &str, exp: u64) -> Ratio<BigInt> {
        giac_poly::coeff_at(p, &giac_poly::Var::from(var), exp)
    }
}
