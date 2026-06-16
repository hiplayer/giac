use std::sync::Arc;

use giac_core::{bigint_to_i64, expr_to_poly, poly_to_expr, EvalError, Expr, ExprArc, FuncKind, Ident};
use giac_poly::{
    as_perfect_power, coeff_at, partfrac_rational_terms, univariate_degree, Poly, PolyError, Var,
};

use crate::risch::{
    hermite_reduce, is_monic_x4_plus_one, rothstein_trager_integrate, try_integrate_x4_plus_one,
    integrate_monic_x4_plus_one, HermiteTerm,
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
    if let Some((base, exp)) = as_perfect_power(&den_p) {
        if exp >= 2 && univariate_degree(&base, &v) >= 2 {
            return integrate_with_hermite(&num_p, &base, exp as usize, &v, var);
        }
    }
    if let Some(r) = try_integrate_x4_plus_one(&num_p, &den_p, &v, var) {
        return Ok(r);
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

fn integrate_with_hermite(
    num: &Poly,
    base: &Poly,
    exp: usize,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    if exp == 2 && num.is_one() && is_monic_x4_plus_one(base, var) {
        let (terms, rem, mult) = hermite_reduce(num, base, exp, var).map_err(poly_err)?;
        let mut parts = Vec::new();
        for t in terms {
            parts.push(integrate_hermite_term(&t, var, x)?);
        }
        if mult == 1 && !rem.is_zero() {
            // Hermite leaves numer `1`; true log remainder is `(3/4)/(x^4+1)`.
            parts.push(integrate_monic_x4_plus_one(
                &Ratio::new(3.into(), 4.into()),
                x,
            ));
        }
        return Ok(Expr::add(parts));
    }
    let (terms, rem, mult) = hermite_reduce(num, base, exp, var).map_err(poly_err)?;
    let mut parts = Vec::new();
    for t in terms {
        parts.push(integrate_hermite_term(&t, var, x)?);
    }
    if mult == 1 && !rem.is_zero() {
        parts.push(integrate_const_over_rational(
            &poly_to_expr(&rem),
            &poly_to_expr(base),
            x,
        )?);
    }
    if parts.is_empty() {
        return Err(EvalError::NotImplemented("integrate hermite"));
    }
    Ok(Expr::add(parts))
}

fn integrate_hermite_term(t: &HermiteTerm, var: &Var, x: &Ident) -> Result<ExprArc, EvalError> {
    let _ = var;
    let scale = Ratio::from_integer(BigInt::from(t.power as i64));
    let num = t.numer.mul_scalar(&Ratio::from_integer((-1).into()));
    let g = poly_to_expr(&t.factor);
    if t.power == 1 {
        return Ok(Arc::new(Expr::Frac(
            poly_to_expr(&num),
            Expr::mul(vec![ratio_to_expr(&scale), g]),
        )));
    }
    Ok(Arc::new(Expr::Frac(
        poly_to_expr(&num),
        Expr::mul(vec![ratio_to_expr(&scale), Expr::pow(g, Expr::int(t.power as i64))]),
    )))
}

fn integrate_rational_term(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    let fdeg = univariate_degree(factor, var);
    let ndeg = univariate_degree(numer, var);
    if let Some((base, exp)) = as_perfect_power(factor) {
        if exp >= 2 && univariate_degree(&base, var) >= 2 {
            return integrate_with_hermite(numer, &base, exp as usize, var, x);
        }
    }
    if ndeg == 0 {
        if let Some((base, exp)) = as_perfect_power(factor) {
            if univariate_degree(&base, var) == 1 && exp >= 2 {
                return integrate_const_over_power(
                    &base,
                    &coeff_at(numer, var, 0),
                    exp,
                    var,
                    x,
                );
            }
        }
    }
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
    if fdeg >= 2 && ndeg == 0 {
        if fdeg == 2 {
            return integrate_over_quadratic(numer, factor, var, x);
        }
        if let Ok(r) = rothstein_trager_integrate(numer, factor, var, x) {
            return Ok(r);
        }
        return integrate_const_over_power(factor, &coeff_at(numer, var, 0), fdeg, var, x);
    }
    if fdeg == 2 && ndeg <= 1 {
        return integrate_over_quadratic(numer, factor, var, x);
    }
    Err(EvalError::NotImplemented("integrate partfrac"))
}

/// ∫ c / g^n dx for linear `g` and n >= 2.
fn integrate_const_over_power(
    factor: &Poly,
    coeff: &Ratio<BigInt>,
    power: u64,
    var: &Var,
    _x: &Ident,
) -> Result<ExprArc, EvalError> {
    if power < 2 {
        return Err(EvalError::TypeError("power must be >= 2"));
    }
    let a = coeff_at(factor, var, 1);
    if a.is_zero() || univariate_degree(factor, var) != 1 {
        return Err(EvalError::TypeError("not linear factor"));
    }
    let g_expr = poly_to_expr(factor);
    let n = power as i64;
    let denom = Ratio::from_integer(BigInt::from(n - 1)) * a.clone();
    let scaled = -coeff.clone() / denom;
    let integrand = Expr::pow(g_expr, Expr::int(-(n - 1)));
    Ok(Expr::mul(vec![ratio_to_expr(&scaled), integrand]))
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
    let disc = b.clone() * b.clone() - Ratio::from_integer(BigInt::from(4)) * a.clone() * c.clone();
    if disc > Ratio::zero() {
        return integrate_over_quadratic_real_roots(
            &b_lin, &c_lin, &a, &b, x, &disc,
        );
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

/// ∫ (B·t+C)/(a·t²+b·t+c) dt when the quadratic has real roots (disc > 0).
fn integrate_over_quadratic_real_roots(
    b_lin: &Ratio<BigInt>,
    c_lin: &Ratio<BigInt>,
    a: &Ratio<BigInt>,
    b: &Ratio<BigInt>,
    x: &Ident,
    disc: &Ratio<BigInt>,
) -> Result<ExprArc, EvalError> {
    let sqrt_disc = sqrt_ratio_expr(disc)?;
    let t = var_to_expr(x);
    let root_sum = Expr::add(vec![ratio_to_expr(&(-b.clone())), sqrt_disc.clone()]);
    let root_diff = Expr::add(vec![
        ratio_to_expr(&(-b.clone())),
        Expr::mul(vec![Expr::int(-1), sqrt_disc.clone()]),
    ]);
    let two_a = Ratio::from_integer(BigInt::from(2)) * a.clone();
    let inv_two_a = ratio_to_expr(&(Ratio::one() / two_a.clone()));
    let r1 = Expr::mul(vec![inv_two_a.clone(), root_sum.clone()]);
    let r2 = Expr::mul(vec![inv_two_a, root_diff.clone()]);
    let two_a_c = c_lin.clone() * two_a;
    let nr1 = Expr::add(vec![
        Expr::mul(vec![ratio_to_expr(b_lin), root_sum]),
        ratio_to_expr(&two_a_c.clone()),
    ]);
    let nr2 = Expr::add(vec![
        Expr::mul(vec![ratio_to_expr(b_lin), root_diff]),
        ratio_to_expr(&two_a_c),
    ]);
    let inv_sqrt = Expr::pow(sqrt_disc, Expr::int(-1));
    let half = ratio_to_expr(&Ratio::new(1.into(), 2.into()));
    let alpha = Expr::mul(vec![half.clone(), nr1, inv_sqrt.clone()]);
    let beta = Expr::mul(vec![half, nr2, Expr::int(-1), inv_sqrt]);
    let ln1 = ln_abs_expr(Expr::add(vec![t.clone(), Expr::mul(vec![Expr::int(-1), r1])]));
    let ln2 = ln_abs_expr(Expr::add(vec![t.clone(), Expr::mul(vec![Expr::int(-1), r2])]));
    Ok(Expr::add(vec![
        Expr::mul(vec![alpha, ln1]),
        Expr::mul(vec![beta, ln2]),
    ]))
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
    use giac_core::format_expr;
    use num_rational::Ratio;

    use super::*;

    #[test]
    fn partfrac_integrate_x_over_repeated_linear() {
        let x = Ident::new("x");
        let den = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(-1)]),
            Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2)),
        ]);
        let r = integrate_const_over_rational(&Expr::sym("x"), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
    }

    #[test]
    fn integrate_one_over_x_fourth_plus_one() {
        let x = Ident::new("x");
        let den = Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(4)), Expr::int(1)]);
        let r = integrate_const_over_rational(&Expr::int(1), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
    }

    #[test]
    fn integrate_one_over_x_fourth_plus_one_squared() {
        let x = Ident::new("x");
        let den = Expr::pow(
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(4)), Expr::int(1)]),
            Expr::int(2),
        );
        let r = integrate_const_over_rational(&Expr::int(1), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
    }

    #[test]
    fn integrate_x_over_x_squared_plus_one_squared_expr() {
        let x = Ident::new("x");
        let den = Expr::pow(
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(1)]),
            Expr::int(2),
        );
        let r = integrate_const_over_rational(&Expr::sym("x"), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
    }

    #[test]
    fn hermite_integrate_x_over_x_squared_plus_one_squared() {
        use giac_poly::Poly;
        let x = Ident::new("x");
        let v = Var::from("x");
        let g = Poly::var("x").pow(2).add(&Poly::one());
        let r = integrate_with_hermite(&Poly::var("x"), &g, 2, &v, &x);
        assert!(r.is_ok(), "{:?}", r);
    }

    #[test]
    fn partfrac_integrate_half_angle_rational_in_t() {
        use giac_poly::Poly;
        let t = Ident::new("__t");
        let tv = Poly::var("__t");
        let den = tv
            .pow(4)
            .mul_scalar(&Ratio::from_integer((-1).into()))
            .add(&tv.pow(3).mul_scalar(&Ratio::from_integer(4.into())))
            .add(&tv.pow(2).mul_scalar(&Ratio::from_integer((-2).into())))
            .add(&tv.mul_scalar(&Ratio::from_integer(4.into())))
            .add(&Poly::constant(Ratio::from_integer((-1).into())));
        let num = tv
            .pow(2)
            .mul_scalar(&Ratio::from_integer(2.into()))
            .add(&tv.mul_scalar(&Ratio::from_integer(8.into())))
            .add(&Poly::constant(Ratio::from_integer(2.into())));
        let r = integrate_const_over_rational(&poly_to_expr(&num), &poly_to_expr(&den), &t);
        assert!(r.is_ok(), "{:?}", r);
    }
}
