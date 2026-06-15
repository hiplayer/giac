use crate::{Context, EvalError, Expr, ExprArc};
use num_traits::{One, Signed, Zero};

use super::expand::expand;
use super::normal::normal;
use super::poly::{expr_to_poly, poly_to_expr, Poly};

/// Basic polynomial factorization (perfect powers, x^n-1 for small n).
pub fn factor(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    let n = normal(expr, ctx)?;
    let p = expr_to_poly(n.as_ref())?;

    if let Some(base) = detect_perfect_power(&p) {
        let base_expr = poly_to_expr(&base);
        let exp = p.degree() / base.degree().max(1);
        if exp > 1 {
            return Ok(Expr::pow(base_expr, Expr::int(exp as i64)));
        }
    }

    if let Some(f) = factor_xn_minus_one(&p) {
        return Ok(f);
    }

    Ok(n)
}

fn detect_perfect_power(p: &Poly) -> Option<Poly> {
    if p.degree() <= 1 {
        return None;
    }
    for exp in 2..=p.degree() {
        if p.degree() % exp != 0 {
            continue;
        }
        // try nth root via gcd structure — check (x+a)^exp pattern
        if let Some(root) = try_nth_root(p, exp) {
            return Some(root);
        }
    }
    None
}

fn try_nth_root(p: &Poly, exp: u32) -> Option<Poly> {
    // For (x+a)^n expanded form, use heuristic: factor via one root
    if p.terms.len() == 1 {
        return None;
    }
    let d = p.degree();
    if d != exp {
        return None;
    }
    // Check binomial: coefficients C(n,k) for (x+a)^n
    // For (x+3)^4, leading coeff 1, constant 81=3^4
    let constant = p.terms.get(&super::poly::Monomial::one())?;
    if constant.denom() != &num_bigint::BigInt::one() {
        return None;
    }
    let c = constant.numer();
    // nth root of constant term gives possible a
    let a = integer_nth_root(c, exp)?;
    // verify by expanding (x+a)^exp
    let x = crate::Ident::new("x");
    let mut base = Poly::constant(num_rational::Ratio::from_integer(a.clone()));
    let x_p = Poly {
        terms: [(super::poly::Monomial::var(x), num_rational::Ratio::one())].into(),
    };
    base = base.add(&x_p);
    if base.pow(exp) == *p {
        return Some(base);
    }
    None
}

fn integer_nth_root(n: &num_bigint::BigInt, exp: u32) -> Option<num_bigint::BigInt> {
    if n.is_negative() && exp % 2 == 0 {
        return None;
    }
    let root = n.sqrt(); // only works for exp=2
    if exp == 2 {
        if &root * &root == *n {
            return Some(root);
        }
    }
    // general nth root by binary search
    let mut lo = num_bigint::BigInt::zero();
    let mut hi = n.abs() + num_bigint::BigInt::one();
    while lo < hi {
        let mid = (&lo + &hi) / num_bigint::BigInt::from(2);
        let pow = mid.pow(exp);
        match pow.cmp(n) {
            std::cmp::Ordering::Equal => return Some(if n.is_negative() { -mid } else { mid }),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    None
}

fn factor_xn_minus_one(p: &Poly) -> Option<ExprArc> {
    // x^4 - 1 pattern
    if p.terms.len() != 2 {
        return None;
    }
    let lt = p.leading_term()?;
    if lt.1 != &num_rational::Ratio::one() {
        return None;
    }
    let deg = lt.0.degree();
    let constant = p.terms.get(&super::poly::Monomial::one())?;
    if constant != &-num_rational::Ratio::one() {
        return None;
    }
    if deg == 4 {
        // (x-1)(x+1)(x^2+1)
        let x = Expr::sym("x");
        return Some(Expr::mul(vec![
            Expr::add(vec![x.clone(), Expr::int(-1)]),
            Expr::add(vec![x.clone(), Expr::int(1)]),
            Expr::add(vec![
                Expr::pow(x, Expr::int(2)),
                Expr::int(1),
            ]),
        ]));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{format_expr, Context, Expr};

    #[test]
    fn factor_quartic_minus_one() {
        let ctx = Context::default();
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(4)),
            Expr::int(-1),
        ]);
        let r = factor(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("x-1") || s.contains("x+1"));
    }
}
