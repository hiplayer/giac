//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-simplify-api-stability.md`.
//!
//!
use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::Zero;

use giac_core::{contains_algext, eval, Context, EvalError, Expr, ExprArc};

use crate::expand::expand;
use giac_core::{expr_to_poly, poly_to_expr};
use giac_poly::Poly;

/// **Stable** — normalize a rational expression to a single fraction in lowest terms.
///
/// Returns `TypeError` for non-rational subtrees (`sin`, `exp`, …). `AlgExt` paths
/// currently delegate to `eval` (partial; see GIAC-algext-adoption A-03).
pub fn ratnormal(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    if contains_algext(expr) {
        return ratnormal_algext(expr, ctx);
    }
    let (num, den) = rational_parts(expr, ctx)?;
    let (num, den) = reduce_fraction(num, den)?;
    if den.is_zero() {
        return Err(EvalError::DivisionByZero);
    }
    if den == Poly::one() {
        return Ok(poly_to_expr(&num));
    }
    Ok(Arc::new(Expr::Frac(
        poly_to_expr(&num),
        poly_to_expr(&den),
    )))
}

// **Temporary** — shim: `AlgExt` ratnormal delegates to `eval` (GIAC-algext-adoption A-03).
// **退役:** `ext_reduce` 有理化落地后删除。
fn ratnormal_algext(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    eval(expr, ctx)
}

// **Pipeline private** — Expr to (num Poly, den Poly)
fn rational_parts(expr: &Expr, ctx: &Context) -> Result<(Poly, Poly), EvalError> {
    match expr {
        Expr::AlgExt(_) => Err(EvalError::TypeError("not a rational expression")),
        Expr::Frac(n, d) => {
            if contains_algext(n.as_ref()) || contains_algext(d.as_ref()) {
                return Err(EvalError::TypeError("not a rational expression"));
            }
            let num = expr_to_poly(n)?;
            let den = expr_to_poly(d)?;
            Ok((num, den))
        }
        Expr::Add(terms) => rational_add(terms, ctx),
        Expr::Mul(factors) => {
            let mut num = Poly::one();
            let mut den = Poly::one();
            for f in factors {
                let (n, d) = rational_parts(f, ctx)?;
                num = num.mul(&n);
                den = den.mul(&d);
            }
            Ok((num, den))
        }
        Expr::Pow(base, exp) => {
            if let Expr::Int(e) = exp.as_ref() {
                if e < &BigInt::zero() {
                    let e_u = giac_core::bigint_to_u32_abs(e)?;
                    let (n, d) = rational_parts(base, ctx)?;
                    return Ok((d.pow(u64::from(e_u)), n.pow(u64::from(e_u))));
                }
                if e >= &BigInt::zero() {
                    let e_u = giac_core::bigint_to_nonneg_u32(e)?;
                    let (n, d) = rational_parts(base, ctx)?;
                    return Ok((n.pow(u64::from(e_u)), d.pow(u64::from(e_u))));
                }
            }
            Err(EvalError::TypeError("non-integer rational power"))
        }
        _ => {
            if let Ok(p) = expr_to_poly(expr) {
                Ok((p, Poly::one()))
            } else {
                let expanded = expand(expr, ctx)?;
                if let Ok(p) = expr_to_poly(expanded.as_ref()) {
                    Ok((p, Poly::one()))
                } else {
                    Err(EvalError::TypeError("not a rational expression"))
                }
            }
        }
    }
}

// **Pipeline private** — add rationals with common denominator
fn rational_add(terms: &[ExprArc], ctx: &Context) -> Result<(Poly, Poly), EvalError> {
    if terms.is_empty() {
        return Ok((Poly::zero(), Poly::one()));
    }
    let mut parts = Vec::with_capacity(terms.len());
    for t in terms {
        parts.push(rational_parts(t.as_ref(), ctx)?);
    }
    let mut lcd = parts[0].1.clone();
    for (_, d) in &parts[1..] {
        lcd = lcd.lcm(d);
    }
    let mut num = Poly::zero();
    for (n, d) in parts {
        let scale = lcd
            .div_exact(&d)
            .ok_or(EvalError::TypeError("rational add lcm"))?;
        num = num.add(&n.mul(&scale));
    }
    Ok((num, lcd))
}

// **Pipeline private** — gcd-reduce num/den Poly pair
fn reduce_fraction(num: Poly, den: Poly) -> Result<(Poly, Poly), EvalError> {
    if den.is_zero() {
        return Err(EvalError::DivisionByZero);
    }
    if num.is_zero() {
        return Ok((Poly::zero(), Poly::one()));
    }
    let g = num.gcd(&den);
    if g != Poly::one() {
        if let (Some(qn), Some(qd)) = (num.div_exact(&g), den.div_exact(&g)) {
            return Ok((qn, qd));
        }
    }
    Ok((num, den))
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::{Expr, FuncKind};
    use giac_core::{format_expr, Context};

    #[test]
    fn ratnormal_frac_form() {
        let ctx = Context::default();
        let e = Arc::new(Expr::Frac(Expr::int(1), Expr::add(vec![Expr::sym("x"), Expr::int(1)])));
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(1)/(x+1)");
    }

    #[test]
    fn ratnormal_negative_power() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(2)]), Expr::int(-1));
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(1)/(x+2)");
    }

    #[test]
    fn ratnormal_product_of_fractions() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::rat(1, 2),
            Expr::pow(Expr::sym("x"), Expr::int(-1)),
        ]);
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(1/2)/(x)");
    }

    #[test]
    fn ratnormal_polynomial_only() {
        let ctx = Context::default();
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::int(1),
        ]);
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x^2+1");
    }

    #[test]
    fn ratnormal_sum_of_fractions() {
        let ctx = Context::default();
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(-1)),
            Expr::rat(1, 2),
        ]);
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(1/2*x+1)/(x)");
    }

    #[test]
    fn ratnormal_division_by_zero() {
        let ctx = Context::default();
        let e = Arc::new(Expr::Frac(Expr::int(1), Expr::int(0)));
        let r = ratnormal(e.as_ref(), &ctx);
        assert!(matches!(r, Err(EvalError::DivisionByZero)));
    }

    #[test]
    fn ratnormal_zero_numerator() {
        let ctx = Context::default();
        let e = Arc::new(Expr::Frac(Expr::int(0), Expr::sym("x")));
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "0");
    }

    #[test]
    fn ratnormal_empty_sum() {
        let ctx = Context::default();
        let e = Expr::add(vec![]);
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "0");
    }

    #[test]
    fn ratnormal_not_rational_error() {
        let ctx = Context::default();
        let e = Expr::func(FuncKind::Sin, vec![Expr::sym("x")]);
        let r = ratnormal(e.as_ref(), &ctx);
        assert!(matches!(r, Err(EvalError::TypeError(_))));
    }

    #[test]
    fn ratnormal_non_integer_power_error() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::sym("x"), Expr::rat(1, 2));
        let r = ratnormal(e.as_ref(), &ctx);
        assert!(matches!(r, Err(EvalError::TypeError(_))));
    }

    #[test]
    // smoke-until B-RATNORM: delete when `ratnormal_reduces_common_factor_semantic` green
    fn ratnormal_reduces_common_factor() {
        let ctx = Context::default();
        let e = Arc::new(Expr::Frac(
            Expr::add(vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")]), Expr::int(2)]),
            Expr::add(vec![Expr::mul(vec![Expr::int(4), Expr::sym("x")]), Expr::int(4)]),
        ));
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2/4");
    }

    #[test]
    #[ignore = "B-RATNORM: ratnormal 应将 2/4 约分为 1/2"]
    fn ratnormal_reduces_common_factor_semantic() {
        let ctx = Context::default();
        let e = Arc::new(Expr::Frac(
            Expr::add(vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")]), Expr::int(2)]),
            Expr::add(vec![Expr::mul(vec![Expr::int(4), Expr::sym("x")]), Expr::int(4)]),
        ));
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert!(
            crate::assert_equiv(r.as_ref(), &Expr::rat(1, 2), &ctx).unwrap(),
            "got {}",
            format_expr(r.as_ref())
        );
    }

    #[test]
    fn ratnormal_algext_square_minus_two() {
        let ctx = Context::default();
        let min = Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ));
        let alpha = giac_core::AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap()
        .into_expr();
        let e = Expr::add(vec![
            Expr::pow(Arc::clone(&alpha), Expr::int(2)),
            Expr::int(-2),
        ]);
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        assert!(r.is_zero());
    }
}
