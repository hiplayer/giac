use std::sync::Arc;

use crate::{Context, EvalError, Expr, ExprArc};
use num_bigint::BigInt;
use num_traits::Zero;

use super::poly::{expr_to_poly, poly_to_expr};

/// Distribute products over sums and expand powers of sums.
pub fn expand(expr: &Expr, _ctx: &Context) -> Result<ExprArc, EvalError> {
    match expr {
        Expr::Add(terms) => {
            let expanded: Result<Vec<_>, _> =
                terms.iter().map(|t| expand(t, _ctx)).collect();
            Ok(Expr::add(expanded?))
        }
        Expr::Mul(factors) => {
            let mut acc = expand(&factors[0], _ctx)?;
            for f in &factors[1..] {
                acc = expand_mul_pair(acc.as_ref(), expand(f, _ctx)?.as_ref(), _ctx)?;
            }
            Ok(acc)
        }
        Expr::Pow(base, exp) => expand_pow(base, exp, _ctx),
        Expr::Func(_, _) | Expr::Symbol(_) | Expr::Int(_) | Expr::Rat(_) => {
            Ok(Arc::new(expr.clone()))
        }
        other => Ok(Arc::new(other.clone())),
    }
}

fn expand_mul_pair(lhs: &Expr, rhs: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    match (lhs, rhs) {
        (Expr::Add(terms), other) | (other, Expr::Add(terms)) => {
            let (terms, other) = if matches!(lhs, Expr::Add(_)) {
                (terms, other)
            } else {
                // rhs is Add
                match rhs {
                    Expr::Add(t) => (t, lhs),
                    _ => unreachable!(),
                }
            };
            let parts: Result<Vec<_>, _> = terms
                .iter()
                .map(|t| expand_mul_pair(t.as_ref(), other, ctx))
                .collect();
            Ok(Expr::add(parts?))
        }
        _ => Ok(Expr::mul(vec![Arc::new(lhs.clone()), Arc::new(rhs.clone())])),
    }
}

fn expand_pow(base: &ExprArc, exp: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let base_e = expand(base, ctx)?;
    if let Expr::Int(n) = exp.as_ref() {
        if n >= &BigInt::zero() && n <= &BigInt::from(20) {
            let e = crate::num_util::bigint_to_nonneg_u32(n)?;
            if e == 0 {
                return Ok(Expr::int(1));
            }
            if let Expr::Add(terms) = base_e.as_ref() {
                return Ok(expand_binomial(terms, e));
            }
            if e == 1 {
                return Ok(base_e);
            }
            let mut result = base_e.clone();
            for _ in 1..e {
                result = expand_mul_pair(result.as_ref(), base_e.as_ref(), ctx)?;
            }
            return Ok(result);
        }
    }
    Ok(Expr::pow(base_e, Arc::clone(exp)))
}

fn expand_binomial(terms: &[ExprArc], n: u32) -> ExprArc {
    // (a+b+...)^n — use polynomial algebra to avoid deep recursive expand
    let sum = Expr::add(terms.to_vec());
    if let Ok(p) = expr_to_poly(sum.as_ref()) {
        return poly_to_expr(&p.pow(n));
    }
    Expr::pow(sum, Expr::int(n as i64))
}

/// Expand then collect into polynomial form.
pub fn normal(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    let expanded = expand(expr, ctx)?;
    if let Ok(p) = expr_to_poly(expanded.as_ref()) {
        return Ok(poly_to_expr(&p));
    }
    Ok(expanded)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{format_expr, Context, Expr};

    #[test]
    fn expand_square_of_sum() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(3)]), Expr::int(4));
        let r = normal(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("x^4"));
        assert!(s.contains("108"));
    }

    #[test]
    fn expand_distribute_mul_over_add() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::sym("y"),
        ]);
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x*y+1*y");
    }

    #[test]
    fn expand_pow_single_base() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::sym("x"), Expr::int(1));
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(r, Expr::sym("x"));
    }

    #[test]
    fn expand_pow_cube_of_symbol() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::sym("x"), Expr::int(3));
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(x*x)*x");
    }

    #[test]
    fn expand_complex_and_non_poly_normal() {
        let ctx = Context::default();
        let c = Expr::Complex(Expr::sym("a"), Expr::sym("b"));
        let r = expand(&c, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "a+b*i");

        let trig = Expr::func(crate::expr::FuncKind::Sin, vec![Expr::sym("x")]);
        let n = normal(trig.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(n.as_ref()), "sin(x)");
    }

    #[test]
    fn expand_rhs_add_and_zero_power() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::sym("y"),
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
        ]);
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x*y+1*y");

        let z = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(0));
        assert_eq!(expand(z.as_ref(), &ctx).unwrap(), Expr::int(1));
    }

    #[test]
    fn expand_binomial_fallback_for_non_poly() {
        let ctx = Context::default();
        let e = Expr::pow(
            Expr::add(vec![Expr::func(crate::expr::FuncKind::Sin, vec![Expr::sym("x")]), Expr::int(1)]),
            Expr::int(2),
        );
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert!(format_expr(r.as_ref()).contains("sin(x)"));
    }
}
