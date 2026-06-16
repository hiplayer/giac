use std::sync::Arc;

use giac_core::{bigint_to_i64, expr_to_poly, poly_to_expr, EvalError, Expr, ExprArc, FuncKind, Ident};
use giac_poly::{
    coeff_at, partfrac_rational_terms, univariate_degree, Poly, PolyError, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::integrate::{integrate, ln_abs_expr, var_to_expr};

pub fn integrate_const_over_rational(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Result<ExprArc, EvalError> {
    let num_p = expr_to_poly(num)?;
    let den_p = expr_to_poly(den)?;
    let v = Var::from(var.as_str());
    if den_p.is_zero() {
        return Err(EvalError::TypeError("division by zero"));
    }
    let (poly_part, terms) = partfrac_rational_terms(&num_p, &den_p, &v).map_err(poly_err)?;
    let mut parts = Vec::new();
    if let Some(q) = poly_part {
        parts.push(integrate(&poly_to_expr(&q), var)?);
    }
    for (numer, factor) in terms {
        parts.push(integrate_rational_term(&numer, &factor, &v, var)?);
    }
    if parts.is_empty() {
        return Err(EvalError::NotImplemented("integrate partfrac"));
    }
    Ok(Expr::add(parts))
}

fn integrate_rational_term(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    let fdeg = univariate_degree(factor, var);
    let ndeg = univariate_degree(numer, var);
    if fdeg == 1 && ndeg == 0 {
        let coeff = coeff_at(numer, var, 0);
        let a = coeff_at(factor, var, 1);
        if a.is_zero() {
            return Err(EvalError::TypeError("degenerate linear factor"));
        }
        let factor_expr = poly_to_expr(factor);
        let scaled = coeff / a;
        return Ok(Expr::mul(vec![ratio_to_expr(&scaled), ln_abs_expr(factor_expr)]));
    }
    if fdeg == 2 && ndeg <= 1 {
        return integrate_over_quadratic(numer, factor, var, x);
    }
    Err(EvalError::NotImplemented("integrate partfrac"))
}

fn integrate_over_quadratic(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    let b_lin = coeff_at(numer, var, 1);
    let c_lin = coeff_at(numer, var, 0);
    let a = coeff_at(factor, var, 2);
    let b = coeff_at(factor, var, 1);
    let c = coeff_at(factor, var, 0);
    if a.is_zero() {
        return Err(EvalError::TypeError("not quadratic"));
    }
    let factor_expr = poly_to_expr(factor);
    let two_a = Ratio::from_integer(BigInt::from(2)) * a.clone();
    let mut parts = Vec::new();
    if !b_lin.is_zero() {
        parts.push(Expr::mul(vec![
            ratio_to_expr(&(b_lin.clone() / two_a.clone())),
            ln_abs_expr(factor_expr.clone()),
        ]));
    }
    let c_eff = c_lin - b_lin * b.clone() / two_a.clone();
    if !c_eff.is_zero() {
        let disc = b.clone() * b.clone() - Ratio::from_integer(BigInt::from(4)) * a.clone() * c.clone();
        if disc > Ratio::zero() {
            return Err(EvalError::NotImplemented("integrate partfrac"));
        }
        if disc == Ratio::zero() {
            return Err(EvalError::NotImplemented("integrate partfrac"));
        }
        let four_ac = Ratio::from_integer(BigInt::from(4)) * a.clone() * c;
        let neg_disc = four_ac - b.clone() * b.clone();
        let sqrt_neg = sqrt_ratio_expr(&neg_disc)?;
        let atan_arg = Expr::add(vec![
            Expr::mul(vec![
                ratio_to_expr(&(Ratio::from_integer(BigInt::from(2)) * a)),
                var_to_expr(x),
            ]),
            ratio_to_expr(&b),
        ]);
        let inv_sqrt = Expr::pow(Arc::clone(&sqrt_neg), Expr::int(-1));
        parts.push(Expr::mul(vec![
            ratio_to_expr(&(c_eff * Ratio::from_integer(BigInt::from(2)))),
            Arc::clone(&inv_sqrt),
            Expr::func(
                FuncKind::Atan,
                vec![Expr::mul(vec![inv_sqrt, atan_arg])],
            ),
        ]));
    }
    if parts.is_empty() {
        return Ok(Expr::int(0));
    }
    Ok(Expr::add(parts))
}

fn ratio_sqrt(r: &Ratio<BigInt>) -> Result<Ratio<BigInt>, EvalError> {
    let sn = integer_sqrt(r.numer()).ok_or_else(|| EvalError::NotImplemented("integrate partfrac"))?;
    let sd = integer_sqrt(r.denom()).ok_or_else(|| EvalError::NotImplemented("integrate partfrac"))?;
    Ok(Ratio::new(sn, sd))
}

fn integer_sqrt(n: &BigInt) -> Option<BigInt> {
    if n.is_negative() {
        return None;
    }
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    let mut lo = BigInt::zero();
    let mut hi = n.clone() + BigInt::one();
    while lo < hi {
        let mid = (&lo + &hi) / BigInt::from(2);
        let sq = &mid * &mid;
        match sq.cmp(n) {
            std::cmp::Ordering::Equal => return Some(mid),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    None
}

fn sqrt_ratio_expr(r: &Ratio<BigInt>) -> Result<ExprArc, EvalError> {
    if r.is_negative() {
        return Err(EvalError::TypeError("negative under sqrt"));
    }
    if let Ok(s) = ratio_sqrt(r) {
        return Ok(ratio_to_expr(&s));
    }
    let num = r.numer();
    let den = r.denom();
    let mut parts = vec![Expr::func(
        FuncKind::Sqrt,
        vec![Expr::int(
            bigint_to_i64(num).map_err(|_| EvalError::NotImplemented("integrate partfrac"))?,
        )],
    )];
    if !den.is_one() {
        parts.push(Expr::pow(
            Expr::func(
                FuncKind::Sqrt,
                vec![Expr::int(
                    bigint_to_i64(den)
                        .map_err(|_| EvalError::NotImplemented("integrate partfrac"))?,
                )],
            ),
            Expr::int(-1),
        ));
    }
    Ok(Expr::mul(parts))
}

fn ratio_to_expr(r: &Ratio<BigInt>) -> ExprArc {
    if r.is_integer() {
        Expr::int(
            bigint_to_i64(r.numer()).unwrap_or_else(|_| {
                panic!("ratio too large for int")
            }),
        )
    } else {
        Arc::new(Expr::Rat(r.clone()))
    }
}

fn poly_err(e: PolyError) -> EvalError {
    match e {
        PolyError::NotImplemented(s) => EvalError::NotImplemented(s),
        PolyError::TypeError(s) => EvalError::TypeError(s),
        PolyError::DivisionByZero => EvalError::TypeError("division by zero"),
    }
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, Context};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn partfrac_integrate_one_over_x_cubed_plus_one() {
        let x = Ident::new("x");
        let den = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(3)),
            Expr::int(1),
        ]);
        let r = integrate_const_over_rational(&Expr::int(1), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
    }
}
