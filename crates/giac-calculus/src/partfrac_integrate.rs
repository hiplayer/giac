//! Partial-fraction and Hermite integration for rational denominators.
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../.doc/giac-calculus-api-stability.md) §7.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `integrate_one_over_quadratic`, `integrate_const_over_rational` |
//! | **Pipeline private** | `integrate_*`, `den_*`, `hermite_*`, `ratio_*`, `sqrt_*`, `poly_err` |

use std::sync::Arc;

use giac_core::{bigint_to_i64, expr_to_poly, poly_to_expr, EvalError, Expr, ExprArc, FuncKind, Ident};
use giac_poly::{
    as_perfect_power, coeff_at, partfrac_rational_terms, substitute_univariate, try_linear_power,
    univariate_degree, Poly, PolyError, Var,
};

use crate::risch::{
    hermite_reduce, rothstein_trager_integrate, try_algebraic_rt_even_quartic, HermiteTerm,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::integrate::{integrate, ln_abs_expr, var_to_expr};

/// **Stable** — ∫ 1/(ax²+bx+c) dx for constant-coefficient denominator (degree 1 or 2).
pub fn integrate_one_over_quadratic(den: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    let v = Var::from(var.as_str());
    let den_p = expr_to_poly(den)?;
    if den_p.is_zero() {
        return Err(EvalError::TypeError("division by zero"));
    }
    let deg = univariate_degree(&den_p, &v);
    match deg {
        0 => Err(EvalError::TypeError("constant denominator")),
        1 => {
            let a = coeff_at(&den_p, &v, 1);
            if a.is_zero() {
                return Err(EvalError::TypeError("degenerate linear denominator"));
            }
            Ok(Expr::mul(vec![
                ratio_to_expr(&(Ratio::one() / a)),
                ln_abs_expr(poly_to_expr(&den_p)),
            ]))
        }
        2 => {
            let a = coeff_at(&den_p, &v, 2);
            let b = coeff_at(&den_p, &v, 1);
            let c = coeff_at(&den_p, &v, 0);
            let disc = b.clone() * b - Ratio::from_integer(BigInt::from(4)) * a * c;
            if disc > Ratio::zero() {
                return integrate_rational_partfrac(&Poly::one(), &den_p, var);
            }
            integrate_over_quadratic(&Poly::one(), &den_p, &v, var)
        }
        _ => Err(EvalError::NotImplemented("integrate quadratic")),
    }
}

// **Pipeline private** — integrate rational via partial fractions.
fn integrate_rational_partfrac(
    num: &Poly,
    den: &Poly,
    var: &Ident,
) -> Result<ExprArc, EvalError> {
    let v = Var::from(var.as_str());
    let (poly_part, terms) = partfrac_rational_terms(num, den, &v).map_err(poly_err)?;
    let mut parts = Vec::new();
    if let Some(q) = poly_part {
        parts.push(integrate(&poly_to_expr(&q), var)?);
    }
    for (numer, factor) in terms {
        if numer.is_zero() {
            continue;
        }
        parts.push(integrate_rational_term(&numer, &factor, &v, var)?);
    }
    if parts.is_empty() {
        return Err(EvalError::NotImplemented("integrate partfrac"));
    }
    Ok(Expr::add(parts))
}

/// **Stable** — ∫ num/den dx for rational expressions (Hermite, Rothstein–Trager, partfrac).
pub fn integrate_const_over_rational(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Result<ExprArc, EvalError> {
    let v = Var::from(var.as_str());
    let num_p = expr_to_poly(num)?;
    if let Some((base, exp)) = den_perfect_power_expr(den) {
        let base_p = expr_to_poly(&base)?;
        if exp >= 2 && univariate_degree(&base_p, &v) >= 2 {
            return integrate_with_hermite(&num_p, &base_p, exp, &v, var);
        }
    }
    let den_p = expr_to_poly(den)?;
    if den_p.is_zero() {
        return Err(EvalError::TypeError("division by zero"));
    }
    if let Some((base, exp)) = as_perfect_power(&den_p) {
        if exp >= 2 && univariate_degree(&base, &v) >= 2 {
            return integrate_with_hermite(&num_p, &base, exp as usize, &v, var);
        }
    }
    if let Some(r) = try_algebraic_rt_even_quartic(&num_p, &den_p, &v, var) {
        return Ok(r);
    }
    integrate_rational_partfrac(&num_p, &den_p, var)
}

// **Pipeline private** — detect `base^exp` denominator with `exp >= 2`.
fn den_perfect_power_expr(den: &ExprArc) -> Option<(ExprArc, usize)> {
    match den.as_ref() {
        Expr::Pow(base, exp) => {
            let n = match exp.as_ref() {
                Expr::Int(i) => bigint_to_i64(i).ok()?,
                _ => return None,
            };
            if n >= 2 {
                Some((Arc::clone(base), n as usize))
            } else {
                None
            }
        }
        _ => None,
    }
}

// **Pipeline private** — sign correction for Hermite quadratic factors.
fn hermite_factor_sign(factor: &Poly, var: &Var) -> Ratio<BigInt> {
    if coeff_at(factor, var, 0) < Ratio::zero() {
        Ratio::from_integer((-1).into())
    } else {
        Ratio::one()
    }
}

// **Pipeline private** — integrate via Hermite reduction on repeated quadratics.
fn integrate_with_hermite(
    num: &Poly,
    base: &Poly,
    exp: usize,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    let (terms, rem, mult) = hermite_reduce(num, base, exp, var).map_err(poly_err)?;
    let mut parts = Vec::new();
    for t in terms {
        parts.push(integrate_hermite_term(&t, var, x)?);
    }
    if mult == 1 && !rem.is_zero() {
        let mut rem_num = rem.clone();
        if hermite_factor_sign(base, var) < Ratio::zero() {
            rem_num = rem_num.mul_scalar(&Ratio::from_integer((-1).into()));
        }
        parts.push(integrate_const_over_rational(
            &poly_to_expr(&rem_num),
            &poly_to_expr(base),
            x,
        )?);
    }
    if parts.is_empty() {
        return Err(EvalError::NotImplemented("integrate hermite"));
    }
    Ok(Expr::add(parts))
}

// **Pipeline private** — integrate one Hermite reduction term.
fn integrate_hermite_term(t: &HermiteTerm, var: &Var, x: &Ident) -> Result<ExprArc, EvalError> {
    let scale = Ratio::from_integer(BigInt::from(t.power as i64));
    let sign = hermite_factor_sign(&t.factor, var);
    let num = t.numer.mul_scalar(&(-sign));
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

// **Pipeline private** — integrate one partial-fraction term.
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
    if try_linear_power(factor, var)
        .is_some_and(|(base, exp)| univariate_degree(&base, var) == 1 && exp >= 2)
    {
        return integrate_polynomial_over_linear_power(numer, factor, var, x);
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
    if fdeg == 1 && ndeg >= 1 {
        return integrate_poly_over_linear(numer, factor, var, x);
    }
    if fdeg >= 2 && ndeg == 0 {
        if fdeg == 2 {
            return integrate_over_quadratic(numer, factor, var, x);
        }
        if let Ok(r) = rothstein_trager_integrate(numer, factor, var, x) {
            return Ok(r);
        }
        return Err(EvalError::NotImplemented("integrate partfrac"));
    }
    if fdeg == 2 && ndeg <= 1 {
        return integrate_over_quadratic(numer, factor, var, x);
    }
    Err(EvalError::NotImplemented("integrate partfrac"))
}

// **Pipeline private** — ∫ P(x)/(cx+d) dx with linear denominator.
fn integrate_poly_over_linear(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    let a = coeff_at(factor, var, 1);
    if a.is_zero() || univariate_degree(factor, var) != 1 {
        return Err(EvalError::TypeError("not linear factor"));
    }
    let (q, r) = numer.div_rem(factor);
    if univariate_degree(&r, var) > 0 {
        return Err(EvalError::TypeError("non-constant remainder"));
    }
    let mut parts = Vec::new();
    if !q.is_zero() {
        parts.push(integrate(&poly_to_expr(&q), x)?);
    }
    if !r.is_zero() {
        let coeff = coeff_at(&r, var, 0);
        let scaled = coeff / a;
        parts.push(Expr::mul(vec![
            ratio_to_expr(&scaled),
            ln_abs_expr(poly_to_expr(factor)),
        ]));
    }
    if parts.is_empty() {
        Ok(Expr::int(0))
    } else {
        Ok(Expr::add(parts))
    }
}

// **Pipeline private** — ∫ c / g^n dx for linear `g` and n >= 2.
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

// **Pipeline private** — ∫ P(x)/g^n dx for `g = c·(ax+b)^n`.
fn integrate_polynomial_over_linear_power(
    numer: &Poly,
    factor: &Poly,
    var: &Var,
    x: &Ident,
) -> Result<ExprArc, EvalError> {
    let (lin, power) = try_linear_power(factor, var).ok_or(EvalError::TypeError(
        "not linear power",
    ))?;
    if univariate_degree(&lin, var) != 1 || power < 2 {
        return Err(EvalError::TypeError("not linear power"));
    }
    let slope = coeff_at(&lin, var, 1);
    if slope.is_zero() {
        return Err(EvalError::TypeError("not linear factor"));
    }
    let lin_pow = lin.pow(power);
    let (quotient, rem) = factor.div_rem(&lin_pow);
    if !rem.is_zero() || univariate_degree(&quotient, var) != 0 {
        return Err(EvalError::TypeError("not linear power"));
    }
    let extra = coeff_at(&quotient, var, 0);
    let root = -coeff_at(&lin, var, 0) / slope.clone();
    let shifted = substitute_univariate(
        numer,
        var,
        &Poly::var(var.clone()).add(&Poly::constant(root.clone())),
    );
    let denom_scale = ratio_pow(&slope, power) * extra;
    let u = Expr::add(vec![
        var_to_expr(x),
        Expr::mul(vec![Expr::rat(-1, 1), ratio_to_expr(&root)]),
    ]);
    let mut parts = Vec::new();
    for j in 0..=univariate_degree(&shifted, var) {
        let c = coeff_at(&shifted, var, j);
        if c.is_zero() {
            continue;
        }
        let coeff = c / denom_scale.clone();
        let exp = j as i64 - power as i64;
        if exp == -1 {
            parts.push(Expr::mul(vec![ratio_to_expr(&coeff), ln_abs_expr(u.clone())]));
        } else if exp < -1 {
            let new_exp = exp + 1;
            let scaled = coeff / Ratio::from_integer(BigInt::from(new_exp));
            parts.push(Expr::mul(vec![
                ratio_to_expr(&scaled),
                Expr::pow(u.clone(), Expr::int(new_exp)),
            ]));
        }
    }
    if parts.is_empty() {
        return Err(EvalError::NotImplemented("integrate partfrac"));
    }
    Ok(Expr::add(parts))
}

// **Pipeline private** — raise a rational to an integer power.
fn ratio_pow(r: &Ratio<BigInt>, n: u64) -> Ratio<BigInt> {
    let mut out = Ratio::one();
    for _ in 0..n {
        out *= r.clone();
    }
    out
}

// **Pipeline private** — ∫ rational over irreducible or repeated quadratic.
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
    let c_eff = c_lin.clone() - b_lin.clone() * b.clone() / two_a.clone();
    if !c_eff.is_zero() {
        if disc == Ratio::zero() {
            let root = -b.clone() / two_a;
            if b_lin.is_zero() {
                let u = if root.is_zero() {
                    var_to_expr(x)
                } else {
                    Expr::add(vec![
                        var_to_expr(x),
                        Expr::mul(vec![Expr::int(-1), ratio_to_expr(&root)]),
                    ])
                };
                let scaled = -c_lin.clone() / a.clone();
                return Ok(Expr::mul(vec![ratio_to_expr(&scaled), Expr::pow(u, Expr::int(-1))]));
            }
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

// **Pipeline private** — ∫ (B·t+C)/(a·t²+b·t+c) dt when the quadratic has real roots (disc > 0).
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

// **Pipeline private** — exact square root of a perfect-square rational.
fn ratio_sqrt(r: &Ratio<BigInt>) -> Result<Ratio<BigInt>, EvalError> {
    let sn = integer_sqrt(r.numer()).ok_or_else(|| EvalError::NotImplemented("integrate partfrac"))?;
    let sd = integer_sqrt(r.denom()).ok_or_else(|| EvalError::NotImplemented("integrate partfrac"))?;
    Ok(Ratio::new(sn, sd))
}

// **Pipeline private** — integer square root by binary search.
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

// **Pipeline private** — build `sqrt(r)` as `Expr` (exact or nested `sqrt`).
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

// **Pipeline private** — convert `Ratio<BigInt>` to `ExprArc`.
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

// **Pipeline private** — map `PolyError` to `EvalError`.
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
    fn integrate_poly_over_linear_term() {
        let x = Ident::new("x");
        let v = Var::from("x");
        let numer = Poly::var("x").add(&Poly::constant(Ratio::from_integer(2.into())));
        let factor = Poly::var("x").sub(&Poly::constant(Ratio::from_integer(3.into())));
        let r = integrate_rational_term(&numer, &factor, &v, &x);
        assert!(r.is_ok(), "{:?}", r);
    }

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
    fn integrate_one_over_x_fourth_minus_one_squared() {
        let x = Ident::new("x");
        let den = Expr::pow(
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(4)), Expr::int(-1)]),
            Expr::int(2),
        );
        let r = integrate_const_over_rational(&Expr::int(1), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
        eprintln!("{}", format_expr(r.unwrap().as_ref()));
    }

    #[test]
    fn integrate_x_over_x_plus_one_times_x_fourth_minus_one() {
        let x = Ident::new("x");
        let den = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(4)), Expr::int(-1)]),
        ]);
        let r = integrate_const_over_rational(&Expr::sym("x"), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
        eprintln!("{}", format_expr(r.unwrap().as_ref()));
    }

    #[test]
    fn integrate_one_over_x_fourth_plus_one_fourth_power() {
        let x = Ident::new("x");
        let den = Expr::pow(
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(4)), Expr::int(1)]),
            Expr::int(4),
        );
        let r = integrate_const_over_rational(&Expr::int(1), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
        eprintln!("{}", format_expr(r.unwrap().as_ref()));
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

    #[test]
    #[test]
    fn integrate_one_over_quadratic_one_minus_x_squared() {
        let x = Ident::new("x");
        let den = Expr::add(vec![
            Expr::int(1),
            Expr::mul(vec![Expr::int(-1), Expr::pow(Expr::sym("x"), Expr::int(2))]),
        ]);
        let r = integrate_one_over_quadratic(&den, &x);
        assert!(r.is_ok(), "{:?}", r);
    }

    #[test]
    fn integrate_ck_int_05_reciprocal() {
        let x = Ident::new("x");
        let den = Expr::mul(vec![
            Expr::int(3),
            Expr::sym("x"),
            Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::sym("x"),
                Expr::int(1),
            ]),
            Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(-1)]), Expr::int(3)),
        ]);
        let r = integrate_const_over_rational(&Expr::int(1), &den, &x);
        assert!(r.is_ok(), "{:?}", r);
    }
}
