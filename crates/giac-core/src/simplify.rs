use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::context::Context;
use crate::error::EvalError;
use crate::expr::{Expr, ExprArc};

/// Flatten nested Add/Mul and combine numeric coefficients.
pub fn simplify(expr: &Expr, _ctx: &Context) -> Result<ExprArc, EvalError> {
    match expr {
        Expr::Add(terms) => simplify_add(terms),
        Expr::Mul(factors) => simplify_mul(factors),
        Expr::Pow(base, exp) => simplify_pow(base, exp),
        Expr::Func(kind, args) => {
            let args: Result<Vec<_>, _> = args.iter().map(|a| simplify(a, _ctx)).collect();
            Ok(Expr::func(*kind, args?))
        }
        Expr::Complex(re, im) => {
            let re = simplify(re, _ctx)?;
            let im = simplify(im, _ctx)?;
            Ok(Arc::new(Expr::Complex(re, im)))
        }
        other => Ok(Arc::new(other.clone())),
    }
}

fn simplify_add(terms: &[ExprArc]) -> Result<ExprArc, EvalError> {
    let mut flat = Vec::new();
    for t in terms {
        match t.as_ref() {
            Expr::Add(inner) => flat.extend(inner.iter().cloned()),
            other => flat.push(Arc::new(other.clone())),
        }
    }

    let mut int_sum = BigInt::zero();
    let mut rat_sum = Ratio::<BigInt>::zero();
    let mut symbolic = Vec::new();

    for t in flat {
        match t.as_ref() {
            Expr::Int(n) => int_sum += n,
            Expr::Rat(r) => rat_sum += r,
            _ => symbolic.push(t),
        }
    }

    if !int_sum.is_zero() {
        symbolic.insert(0, Expr::int(int_to_i64(&int_sum)?));
    }
    if !rat_sum.is_zero() {
        if rat_sum.denom() == &BigInt::one() {
            if int_sum.is_zero() {
                symbolic.insert(0, Expr::int(int_to_i64(rat_sum.numer())?));
            }
        } else {
            symbolic.push(rat_to_expr(&rat_sum)?);
        }
    }

    match symbolic.len() {
        0 => Ok(Expr::int(0)),
        1 => Ok(Arc::clone(&symbolic[0])),
        _ => Ok(Expr::add(symbolic)),
    }
}

fn simplify_mul(factors: &[ExprArc]) -> Result<ExprArc, EvalError> {
    let mut flat = Vec::new();
    for f in factors {
        match f.as_ref() {
            Expr::Mul(inner) => flat.extend(inner.iter().cloned()),
            other => flat.push(Arc::new(other.clone())),
        }
    }

    let mut int_prod = BigInt::one();
    let mut rat_prod = Ratio::<BigInt>::one();
    let mut symbolic = Vec::new();

    for f in flat {
        match f.as_ref() {
            Expr::Int(n) => int_prod *= n,
            Expr::Rat(r) => rat_prod *= r,
            _ => symbolic.push(f),
        }
    }

    if int_prod != BigInt::one() {
        symbolic.insert(0, Expr::int(int_to_i64(&int_prod)?));
    }
    if rat_prod != Ratio::one() {
        if rat_prod.denom() == &BigInt::one() && int_prod == BigInt::one() {
            symbolic.insert(0, Expr::int(int_to_i64(rat_prod.numer())?));
        } else if rat_prod != Ratio::one() {
            symbolic.push(rat_to_expr(&rat_prod)?);
        }
    }

    match symbolic.len() {
        0 => Ok(Expr::int(1)),
        1 => Ok(Arc::clone(&symbolic[0])),
        _ => Ok(Expr::mul(symbolic)),
    }
}

fn simplify_pow(base: &ExprArc, exp: &ExprArc) -> Result<ExprArc, EvalError> {
    match (base.as_ref(), exp.as_ref()) {
        (Expr::Int(b), Expr::Int(e)) if e >= &BigInt::from(0) && e <= &BigInt::from(20) => {
            let e = e.to_string().parse::<u32>().unwrap();
            Ok(Expr::int(int_to_i64(&b.pow(e))?))
        }
        (Expr::Int(b), Expr::Int(e)) if e == &BigInt::from(2) => {
            let sq = b * b;
            Ok(Expr::int(int_to_i64(&sq)?))
        }
        _ => Ok(Expr::pow(Arc::clone(base), Arc::clone(exp))),
    }
}

fn int_to_i64(n: &BigInt) -> Result<i64, EvalError> {
    n.to_string()
        .parse()
        .map_err(|_| EvalError::TypeError("integer out of i64 range"))
}

fn rat_to_expr(r: &Ratio<BigInt>) -> Result<ExprArc, EvalError> {
    if r.denom() == &BigInt::one() {
        Ok(Expr::int(int_to_i64(r.numer())?))
    } else {
        Ok(Arc::new(Expr::Rat(r.clone())))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Expr;

    #[test]
    fn flatten_add() {
        let ctx = Context::new();
        let e = Expr::add(vec![
            Expr::add(vec![Expr::int(1), Expr::int(2)]),
            Expr::int(3),
        ]);
        let s = simplify(e.as_ref(), &ctx).unwrap();
        assert_eq!(s, Expr::int(6));
    }

    #[test]
    fn flatten_mul() {
        let ctx = Context::new();
        let e = Expr::mul(vec![Expr::int(2), Expr::int(3), Expr::int(4)]);
        let s = simplify(e.as_ref(), &ctx).unwrap();
        assert_eq!(s, Expr::int(24));
    }
}
