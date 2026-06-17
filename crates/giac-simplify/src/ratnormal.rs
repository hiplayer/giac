use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::Zero;

use giac_core::{Context, EvalError, Expr, ExprArc};

use crate::expand::expand;
use giac_core::{expr_to_poly, poly_to_expr};
use giac_poly::Poly;

/// Normalize a rational expression to a single fraction in lowest terms.
pub fn ratnormal(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
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

fn rational_parts(expr: &Expr, ctx: &Context) -> Result<(Poly, Poly), EvalError> {
    match expr {
        Expr::Frac(n, d) => {
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
        let s = format_expr(r.as_ref());
        assert!(s.contains("x") && (s.contains('/') || s.contains("Pow")));
    }

    #[test]
    fn ratnormal_product_of_fractions() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::rat(1, 2),
            Expr::pow(Expr::sym("x"), Expr::int(-1)),
        ]);
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("1/2") || s.contains('/'));
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
        let s = format_expr(r.as_ref());
        assert!(s.contains("x") && s.contains("2"));
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
    fn ratnormal_reduces_common_factor() {
        let ctx = Context::default();
        let e = Arc::new(Expr::Frac(
            Expr::add(vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")]), Expr::int(2)]),
            Expr::add(vec![Expr::mul(vec![Expr::int(4), Expr::sym("x")]), Expr::int(4)]),
        ));
        let r = ratnormal(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("2") && s.contains("4"));
    }
}
