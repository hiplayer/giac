//! Mathematical equivalence helpers (`assert_equiv` per conformance-testing.md §3).

use std::sync::Arc;

use num_rational::Ratio;
use num_traits::{One, Zero};

use giac_core::Context;
use giac_core::EvalError;
use giac_core::{Expr, ExprArc, FuncKind};

use crate::expand::normal;

/// Construct `a - b` as an expression tree.
pub fn sub(a: &Expr, b: &Expr) -> Result<ExprArc, EvalError> {
    Ok(Expr::add(vec![
        Arc::new(a.clone()),
        Expr::mul(vec![Expr::int(-1), Arc::new(b.clone())]),
    ]))
}

/// True when `normal(e)` simplifies to zero.
pub fn is_zero(e: &Expr, ctx: &Context) -> Result<bool, EvalError> {
    match e {
        Expr::Rat(r) => return Ok(r.is_zero()),
        Expr::Complex(re, im) => {
            return Ok(is_zero(re.as_ref(), ctx)? && is_zero(im.as_ref(), ctx)?);
        }
        _ => {}
    }
    let n = normal(e, ctx)?;
    Ok(n.is_zero())
}

/// True when `a` and `b` are mathematically equivalent under `normal`.
pub fn assert_equiv(a: &Expr, b: &Expr, ctx: &Context) -> Result<bool, EvalError> {
    if a == b {
        return Ok(true);
    }
    if let (Expr::Rat(ra), Expr::Rat(rb)) = (a, b) {
        return Ok(ra == rb);
    }
    if let (Expr::Complex(ar, ai), Expr::Complex(br, bi)) = (a, b) {
        return Ok(
            assert_equiv(ar.as_ref(), br.as_ref(), ctx)?
                && assert_equiv(ai.as_ref(), bi.as_ref(), ctx)?,
        );
    }
    let ca = canonical_radical(a);
    let cb = canonical_radical(b);
    if ca == cb {
        return Ok(true);
    }
    let d = sub(ca.as_ref(), cb.as_ref())?;
    is_zero(d.as_ref(), ctx)
}

/// `1/sqrt(n)` → `sqrt(n)/n` so `sqrt(n)/2` and `1/sqrt(n)` share a shape when equal.
fn canonical_radical(e: &Expr) -> ExprArc {
    match e {
        Expr::Pow(base, exp) => {
            if let Expr::Int(e) = exp.as_ref() {
                if e == &num_bigint::BigInt::from(-1) {
                    if let Some(m) = inv_sqrt_to_mul(base.as_ref()) {
                        return m;
                    }
                }
            }
            Arc::new(e.clone())
        }
        Expr::Mul(factors) => {
            let mut scalar = Ratio::one();
            let mut rest = Vec::new();
            for f in factors {
                match canonical_radical(f.as_ref()).as_ref() {
                    Expr::Int(n) => {
                        scalar *= Ratio::from(n.clone());
                    }
                    Expr::Rat(r) => scalar *= r.clone(),
                    other => rest.push(other.clone()),
                }
            }
            if scalar.is_zero() {
                return Expr::int(0);
            }
            if scalar.is_one() && rest.len() == 1 {
                return Arc::new(rest.remove(0));
            }
            if scalar.is_one() {
                return Expr::mul(rest.into_iter().map(Arc::new).collect());
            }
            let mut out: Vec<ExprArc> = vec![ratio_to_expr(&scalar)];
            out.extend(rest.into_iter().map(Arc::new));
            Expr::mul(out)
        }
        _ => Arc::new(e.clone()),
    }
}

fn ratio_to_expr(r: &Ratio<num_bigint::BigInt>) -> ExprArc {
    let n: i64 = r.numer().to_string().parse().unwrap_or(0);
    let d: i64 = r.denom().to_string().parse().unwrap_or(1);
    if d == 1 {
        Expr::int(n)
    } else {
        Expr::rat(n, d)
    }
}

fn inv_sqrt_to_mul(base: &Expr) -> Option<ExprArc> {
    let args = match base {
        Expr::Func(FuncKind::Sqrt, args) => args,
        _ => return None,
    };
    let n = match args.first()?.as_ref() {
        Expr::Int(v) => v.to_string().parse::<i64>().ok()?,
        _ => return None,
    };
    if n <= 0 {
        return None;
    }
    Some(Expr::mul(vec![
        Expr::rat(1, n),
        Arc::new(base.clone()),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::format_expr;

    #[test]
    fn equiv_commutative_add() {
        let ctx = Context::default();
        let a = Expr::add(vec![Expr::sym("x"), Expr::int(1)]);
        let b = Expr::add(vec![Expr::int(1), Expr::sym("x")]);
        assert!(assert_equiv(a.as_ref(), b.as_ref(), &ctx).unwrap());
    }

    #[test]
    fn equiv_expanded_square() {
        let ctx = Context::default();
        let a = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2));
        let b = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
            Expr::int(1),
        ]);
        assert!(assert_equiv(a.as_ref(), b.as_ref(), &ctx).unwrap());
    }

    #[test]
    fn equiv_sqrt_half_forms() {
        let ctx = Context::default();
        let a = Expr::mul(vec![
            Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]),
            Expr::rat(1, 2),
        ]);
        let b = Expr::pow(
            Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]),
            Expr::int(-1),
        );
        assert!(assert_equiv(a.as_ref(), b.as_ref(), &ctx).unwrap());
    }

    #[test]
    fn equiv_rational_reduced() {
        let ctx = Context::default();
        assert!(assert_equiv(
            Expr::rat(1, 2).as_ref(),
            Expr::rat(2, 4).as_ref(),
            &ctx
        )
        .unwrap());
    }

    #[test]
    fn is_zero_complex() {
        let ctx = Context::default();
        let z = Expr::Complex(Expr::int(0), Expr::int(0));
        assert!(is_zero(&z, &ctx).unwrap());
        let c = Expr::Complex(Expr::int(1), Expr::int(0));
        assert!(!is_zero(&c, &ctx).unwrap());
    }

    #[test]
    fn sub_fraction_difference() {
        let ctx = Context::default();
        let d = sub(
            Expr::rat(1, 2).as_ref(),
            Expr::rat(2, 4).as_ref(),
        )
        .unwrap();
        assert!(is_zero(d.as_ref(), &ctx).unwrap());
    }

    #[test]
    fn equiv_mod_difference() {
        let ctx = Context::default();
        let a = Expr::Mod(Expr::int(7), Expr::int(3));
        let b = Expr::int(1);
        assert!(!assert_equiv(&a, b.as_ref(), &ctx).unwrap());
        let diff = sub(&a, b.as_ref()).unwrap();
        assert!(!is_zero(diff.as_ref(), &ctx).unwrap());
        let _ = format_expr(diff.as_ref());
    }

    #[test]
    fn not_equiv_different() {
        let ctx = Context::default();
        let a = Expr::sym("x");
        let b = Expr::add(vec![Expr::sym("x"), Expr::int(1)]);
        assert!(!assert_equiv(a.as_ref(), b.as_ref(), &ctx).unwrap());
    }
}
