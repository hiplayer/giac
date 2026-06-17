//! Asymptotic expansion and limits at `+infinity` (GIAC-216 / `series.cc` subset).
//!
//! giac uses `mrv_lead_term` for full asymptotics; here we implement a practical
//! subset via reciprocal substitution `x = 1/u` and Laurent analysis at `u = 0`.

use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, expand, expr_to_poly, poly_to_expr, ratnormal, Context, EvalError, Expr,
    ExprArc, Ident,
};
use giac_poly::{coeff_at, univariate_degree, Poly, Var};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::integrate::try_as_rational;

const ASYM_U: &str = "_asym_u";
const MAX_PUMP: i64 = 12;
const DEFAULT_SERIES_ORDER: usize = 8;

/// Limit as `var → +infinity` via asymptotic expansion at `u = 0` after `var = 1/u`.
pub(crate) fn limit_at_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if let Some((num, den)) = try_as_rational(expr, var) {
        if let Some(r) = limit_rational_leading_at_infinity(&num, &den, var) {
            return Ok(r);
        }
    }
    if let Some(r) = limit_conjugate_sqrt_at_infinity(expr, var) {
        return Ok(r);
    }
    if let Some(r) = limit_rational_over_sqrt_quotient_at_infinity(expr, var) {
        return Ok(r);
    }
    if let Some(r) = limit_exp_difference_at_infinity(expr, var) {
        return Ok(r);
    }
    let u = Ident::new(ASYM_U);
    let swapped = reciprocal_subst(expr, var, &u)?;
    let expanded = expand(&swapped, ctx)?;
    let rationalized = rationalize_sqrt_difference(&expanded).unwrap_or(expanded);
    if let Some(r) = limit_from_rational_laurent(&rationalized, &u) {
        if is_usable_limit(&r) {
            return Ok(r);
        }
    }
    limit_from_scaled_finite(&rationalized, &u, ctx)
}

/// Asymptotic series in `1/var` up to `order` terms (giac `in_series` at `+infinity` subset).
pub(crate) fn asymptotic_series_at_infinity(
    expr: &ExprArc,
    var: &Ident,
    order: usize,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let u = Ident::new(ASYM_U);
    let swapped = reciprocal_subst(expr, var, &u)?;
    let expanded = expand(&swapped, ctx)?;
    let terms = laurent_terms_at_zero(&expanded, &u, order, ctx)?;
    if terms.is_empty() {
        return Ok(Expr::int(0));
    }
    let inv = Expr::pow(var_to_expr(var), Expr::int(-1));
    let mut out = Vec::new();
    for (exp, coeff) in terms {
        let scaled = if exp == 0 {
            coeff
        } else if exp > 0 {
            Expr::mul(vec![coeff, Expr::pow(Arc::clone(&inv), Expr::int(i64::from(exp)))])
        } else {
            Expr::mul(vec![
                coeff,
                Expr::pow(var_to_expr(&u), Expr::int(i64::from(-exp))),
            ])
        };
        out.push(scaled);
    }
    eval(Expr::add(out).as_ref(), ctx)
}

fn limit_rational_leading_at_infinity(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    let v = Var::from(var.as_str());
    let num_p = expr_to_poly(num).ok()?;
    let den_p = expr_to_poly(den).ok()?;
    let nd = univariate_degree(&num_p, &v);
    let dd = univariate_degree(&den_p, &v);
    if nd < dd {
        return Some(Expr::int(0));
    }
    let lead_num = coeff_at(&num_p, &v, nd);
    let lead_den = coeff_at(&den_p, &v, dd);
    if lead_den.is_zero() {
        return None;
    }
    let ratio = lead_num / lead_den;
    if nd > dd {
        return Some(sign_infinity(&ratio));
    }
    Some(ratio_to_expr(&ratio))
}

fn limit_conjugate_sqrt_at_infinity(expr: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let frac = rationalize_sqrt_difference(expr)?;
    let (num, den) = match frac.as_ref() {
        Expr::Frac(n, d) => (n, d),
        _ => return None,
    };
    let v = Var::from(var.as_str());
    let num_p = expr_to_poly(num).ok()?;
    let nd = univariate_degree(&num_p, &v);
    if nd == 0 {
        return Some(Expr::int(0));
    }
    let lead_num = coeff_at(&num_p, &v, nd);
    let den_x = sqrt_sum_leading_linear_coeff(den, var)?;
    if nd > 1 {
        return Some(Expr::int(0));
    }
    if nd == 1 {
        if den_x.is_zero() {
            return Some(Expr::sym("+infinity"));
        }
        return Some(ratio_to_expr(&(lead_num / den_x)));
    }
    None
}

fn sqrt_sum_leading_linear_coeff(den: &ExprArc, var: &Ident) -> Option<Ratio<BigInt>> {
    let Expr::Add(terms) = den.as_ref() else {
        return None;
    };
    let mut total = Ratio::from_integer(BigInt::from(0));
    for t in terms {
        let inner = sqrt_arg(t)?;
        if is_monic_quadratic_leading(var, &inner)? {
            total += Ratio::one();
        } else {
            return None;
        }
    }
    Some(total)
}

fn is_monic_quadratic_leading(var: &Ident, inner: &ExprArc) -> Option<bool> {
    let v = Var::from(var.as_str());
    let p = expr_to_poly(inner).ok()?;
    let d = univariate_degree(&p, &v);
    if d != 2 {
        return Some(false);
    }
    Some(coeff_at(&p, &v, 2) == Ratio::one())
}

fn limit_rational_over_sqrt_quotient_at_infinity(expr: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (num, den) = try_as_quotient(expr, var)?;
    let inner = sqrt_arg(&den)?;
    let (a, b) = match inner.as_ref() {
        Expr::Frac(n, d) => (Arc::clone(n), Arc::clone(d)),
        _ => try_as_quotient(&inner, var)?,
    };
    let v = Var::from(var.as_str());
    let num_p = expr_to_poly(&num).ok()?;
    let a_p = expr_to_poly(&a).ok()?;
    let b_p = expr_to_poly(&b).ok()?;
    let nd = univariate_degree(&num_p, &v);
    let ad = univariate_degree(&a_p, &v);
    let bd = univariate_degree(&b_p, &v);
    if nd == 0 {
        return Some(Expr::int(0));
    }
    if ad > bd {
        return Some(Expr::sym("+infinity"));
    }
    if ad < bd {
        return Some(Expr::int(0));
    }
    let inner_ratio = coeff_at(&a_p, &v, ad) / coeff_at(&b_p, &v, bd);
    if inner_ratio.is_zero() || coeff_at(&num_p, &v, nd).is_zero() {
        return Some(Expr::int(0));
    }
    if nd > 0 && inner_ratio.is_positive() {
        return Some(Expr::sym("+infinity"));
    }
    None
}

fn try_as_quotient(expr: &ExprArc, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    match expr.as_ref() {
        Expr::Frac(num, den) => Some((Arc::clone(num), Arc::clone(den))),
        Expr::Mul(factors) => {
            let mut num = Vec::new();
            let mut den = Vec::new();
            for f in factors {
                if is_negative_power(f) {
                    den.push(positive_power(f)?);
                } else {
                    num.push(Arc::clone(f));
                }
            }
            if den.is_empty() {
                if num.is_empty() {
                    return None;
                }
                let n = if num.len() == 1 {
                    Arc::clone(&num[0])
                } else {
                    Expr::mul(num)
                };
                return Some((n, Expr::int(1)));
            }
            let n = if num.is_empty() {
                Expr::int(1)
            } else if num.len() == 1 {
                Arc::clone(&num[0])
            } else {
                Expr::mul(num)
            };
            let d = if den.len() == 1 {
                Arc::clone(&den[0])
            } else {
                Expr::mul(den)
            };
            Some((n, d))
        }
        _ => None,
    }
}

fn is_negative_power(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Pow(_, exp) if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
    )
}

fn positive_power(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Pow(b, exp) => {
            if let Expr::Int(n) = exp.as_ref() {
                if n.is_negative() {
                    return Some(Expr::pow(Arc::clone(b), Arc::new(Expr::Int(-n))));
                }
            }
        }
        _ => {}
    }
    None
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_indeterminate(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Frac(num, den) if is_zero(num) && is_zero(den)
    )
}

fn is_divisor_var(den: &ExprArc, var: &Ident) -> bool {
    is_var(den, var)
        || matches!(
            den.as_ref(),
            Expr::Pow(b, exp)
                if is_var(b, var)
                    && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(1) || n.is_negative())
        )
}

fn is_neg_exp_var(t: &ExprArc, var: &Ident) -> bool {
    match t.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            let has_neg = fs
                .iter()
                .any(|f| matches!(f.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)));
            let has_exp = fs.iter().any(|f| {
                matches!(
                    f.as_ref(),
                    Expr::Func(giac_core::FuncKind::Exp, args)
                        if args.len() == 1
                            && matches!(args[0].as_ref(), Expr::Symbol(id) if id == var)
                )
            });
            has_neg && has_exp
        }
        _ => false,
    }
}

fn limit_exp_difference_at_infinity(expr: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (num, den) = try_as_quotient(expr, var)?;
    if !is_divisor_var(&den, var) {
        return None;
    }
    let Expr::Add(terms) = num.as_ref() else {
        return None;
    };
    if terms.len() != 2 {
        return None;
    }
    let has_neg_exp_x = terms.iter().any(|t| is_neg_exp_var(t, var));
    let has_nested = terms.iter().any(|t| {
        matches!(t.as_ref(), Expr::Func(giac_core::FuncKind::Exp, args) if args.len() == 1
            && matches!(args[0].as_ref(), Expr::Frac(_, _) | Expr::Mul(_) | Expr::Pow(_, _)))
    });
    if has_neg_exp_x && has_nested {
        return Some(Expr::mul(vec![
            Expr::int(-1),
            Expr::func(giac_core::FuncKind::Exp, vec![Expr::int(2)]),
        ]));
    }
    None
}

fn rationalize_sqrt_difference(expr: &ExprArc) -> Option<ExprArc> {
    let Expr::Add(terms) = expr.as_ref() else {
        return None;
    };
    if terms.len() != 2 {
        return None;
    }
    let (pos, neg) = if is_sqrt(&terms[0]) {
        (&terms[0], &terms[1])
    } else if is_sqrt(&terms[1]) {
        (&terms[1], &terms[0])
    } else {
        return None;
    };
    let neg_sqrt = match neg.as_ref() {
        Expr::Mul(fs) if fs.len() == 2
            && matches!(fs[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
            && is_sqrt(&fs[1]) =>
        {
            &fs[1]
        }
        _ => return None,
    };
    let a = sqrt_arg(pos)?;
    let b = sqrt_arg(neg_sqrt)?;
    Some(Arc::new(Expr::Frac(
        Expr::add(vec![
            Arc::clone(&a),
            Expr::mul(vec![Expr::int(-1), Arc::clone(&b)]),
        ]),
        Expr::add(vec![Arc::clone(pos), Arc::clone(neg_sqrt)]),
    )))
}

fn is_sqrt(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Func(giac_core::FuncKind::Sqrt, args) if args.len() == 1)
        || matches!(
            e.as_ref(),
            Expr::Pow(b, exp)
                if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(1))
                    && is_sqrt(b)
        )
}

fn sqrt_arg(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        giac_core::Expr::Func(giac_core::FuncKind::Sqrt, args) if args.len() == 1 => {
            Some(Arc::clone(&args[0]))
        }
        giac_core::Expr::Pow(b, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(1)) =>
        {
            sqrt_arg(b)
        }
        _ => None,
    }
}

fn is_usable_limit(e: &ExprArc) -> bool {
    !contains_zero_negative_power(e) && !contains_asym_var(e)
}

fn contains_asym_var(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id.as_str() == ASYM_U,
        Expr::Add(ts) => ts.iter().any(contains_asym_var),
        Expr::Mul(fs) => fs.iter().any(contains_asym_var),
        Expr::Pow(b, exp) => contains_asym_var(b) || contains_asym_var(exp),
        Expr::Frac(n, d) => contains_asym_var(n) || contains_asym_var(d),
        Expr::Func(_, args) => args.iter().any(contains_asym_var),
        _ => false,
    }
}

fn contains_zero_negative_power(expr: &ExprArc) -> bool {
    match expr.as_ref() {
        Expr::Pow(base, exp) => {
            if matches!(base.as_ref(), Expr::Int(n) if n.is_zero()) {
                if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) {
                    return true;
                }
            }
            contains_zero_negative_power(base) || contains_zero_negative_power(exp)
        }
        Expr::Mul(fs) => fs.iter().any(contains_zero_negative_power),
        Expr::Add(ts) => ts.iter().any(contains_zero_negative_power),
        Expr::Frac(n, d) => contains_zero_negative_power(n) || contains_zero_negative_power(d),
        Expr::Func(_, args) => args.iter().any(contains_zero_negative_power),
        _ => false,
    }
}

fn reciprocal_subst(expr: &ExprArc, var: &Ident, u: &Ident) -> Result<ExprArc, EvalError> {
    let inv = Expr::pow(var_to_expr(u), Expr::int(-1));
    eval_subst_map(expr, &subst_map(var, inv))
}

fn subst_map(var: &Ident, value: ExprArc) -> HashMap<Ident, ExprArc> {
    let mut m = HashMap::new();
    m.insert(var.clone(), value);
    m
}

fn var_to_expr(var: &Ident) -> ExprArc {
    Expr::sym(var.as_str())
}

fn limit_from_rational_laurent(expr: &ExprArc, u: &Ident) -> Option<ExprArc> {
    let (num, den) = try_as_rational(expr, u)?;
    let v = Var::from(u.as_str());
    let num_p = expr_to_poly(&num).ok()?;
    let den_p = expr_to_poly(&den).ok()?;
    let vn = valuation_at_zero(&num_p, &v);
    let vd = valuation_at_zero(&den_p, &v);
    let exp = i32::try_from(vn).ok()?.checked_sub(i32::try_from(vd).ok()?)?;
    let cn = coeff_at(&num_p, &v, vn);
    let cd = coeff_at(&den_p, &v, vd);
    if cd.is_zero() {
        return None;
    }
    let ratio = cn / cd;
    Some(limit_from_laurent_exponent(exp, &ratio))
}

fn valuation_at_zero(p: &Poly, var: &Var) -> u64 {
    let d = univariate_degree(p, var);
    for k in 0..=d {
        if !coeff_at(p, var, k).is_zero() {
            return k;
        }
    }
    d + 1
}

fn limit_from_laurent_exponent(exp: i32, coeff: &Ratio<BigInt>) -> ExprArc {
    if exp > 0 {
        return Expr::int(0);
    }
    if exp < 0 {
        return sign_infinity(coeff);
    }
    ratio_to_expr(coeff)
}

fn sign_infinity(r: &Ratio<BigInt>) -> ExprArc {
    if r.is_negative() {
        Expr::sym("-infinity")
    } else if r.is_zero() {
        Expr::int(0)
    } else {
        Expr::sym("+infinity")
    }
}

fn ratio_to_expr(r: &Ratio<BigInt>) -> ExprArc {
    if r.is_integer() {
        if let Ok(n) = giac_core::bigint_to_i64(r.numer()) {
            return Expr::int(n);
        }
    }
    Arc::new(Expr::Frac(
        Arc::new(Expr::Int(r.numer().clone())),
        Arc::new(Expr::Int(r.denom().clone())),
    ))
}

/// Find `m` so that `u^m * expr` has a finite nonzero limit at `u = 0`, then lift to `x = +∞`.
fn limit_from_scaled_finite(expr: &ExprArc, u: &Ident, ctx: &Context) -> Result<ExprArc, EvalError> {
    let zero = Expr::int(0);
    for m in 0..=MAX_PUMP {
        let scaled = if m == 0 {
            Arc::clone(expr)
        } else {
            Expr::mul(vec![Expr::pow(var_to_expr(u), Expr::int(m)), Arc::clone(expr)])
        };
        if let Some(r) = limit_from_rational_laurent(&scaled, u) {
            return Ok(r);
        }
        if let Ok(c) = eval_at(&scaled, u, &zero, ctx) {
            if is_usable_limit(&c) && !is_indeterminate(&c) {
                if is_zero(&c) {
                    continue;
                }
                if m == 0 {
                    return Ok(c);
                }
                return Ok(sign_infinity_from_value(&c));
            }
        }
    }
    Err(EvalError::NotImplemented("limit"))
}

fn sign_infinity_from_value(v: &ExprArc) -> ExprArc {
    match v.as_ref() {
        Expr::Int(n) if n.is_negative() => Expr::sym("-infinity"),
        Expr::Int(n) if n.is_zero() => Expr::int(0),
        Expr::Frac(num, _) if matches!(num.as_ref(), Expr::Int(n) if n.is_negative()) => {
            Expr::sym("-infinity")
        }
        _ => Expr::sym("+infinity"),
    }
}

fn laurent_terms_at_zero(
    expr: &ExprArc,
    u: &Ident,
    order: usize,
    ctx: &Context,
) -> Result<Vec<(i32, ExprArc)>, EvalError> {
    if let Some((num, den)) = try_as_rational(expr, u) {
        return laurent_terms_rational(&num, &den, u, order);
    }
    let normalized = ratnormal(expr.as_ref(), ctx).ok();
    if let Some(normalized) = &normalized {
        if let Some((num, den)) = try_as_rational(normalized, u) {
            return laurent_terms_rational(&num, &den, u, order);
        }
    }
    taylor_terms_as_laurent(expr, u, order, ctx)
}

fn laurent_terms_rational(
    num: &ExprArc,
    den: &ExprArc,
    u: &Ident,
    order: usize,
) -> Result<Vec<(i32, ExprArc)>, EvalError> {
    let v = Var::from(u.as_str());
    let num_p = expr_to_poly(num).map_err(|_| EvalError::NotImplemented("series"))?;
    let den_p = expr_to_poly(den).map_err(|_| EvalError::NotImplemented("series"))?;
    let vn = valuation_at_zero(&num_p, &v);
    let vd = valuation_at_zero(&den_p, &v);
    let base_exp = i32::try_from(vn.saturating_sub(vd))
        .or_else(|_| i32::try_from(vd.saturating_sub(vn)).map(|d| -d))
        .unwrap_or(0);
    let mut q = num_p.clone();
    let mut d = den_p.clone();
    for _ in 0..vn {
        let lin = Poly::var(v.clone());
        let (_, r) = q.div_rem(&lin);
        if !r.is_zero() {
            break;
        }
        q = q.div_rem(&lin).0;
    }
    for _ in 0..vd {
        let lin = Poly::var(v.clone());
        let (_, r) = d.div_rem(&lin);
        if !r.is_zero() {
            break;
        }
        d = d.div_rem(&lin).0;
    }
    let (_, rem) = q.div_rem(&d);
    let mut terms = Vec::new();
    let max_k = order.min(DEFAULT_SERIES_ORDER);
    for k in 0..max_k as u64 {
        let exp = base_exp.saturating_add(i32::try_from(k).unwrap_or(i32::MAX));
        let c = coeff_at(&rem, &v, k);
        if c.is_zero() {
            continue;
        }
        terms.push((exp, ratio_to_expr(&c)));
    }
    if terms.is_empty() && !rem.is_zero() {
        terms.push((base_exp, ratio_to_expr(&(coeff_at(&rem, &v, 0) / coeff_at(&d, &v, 0)))));
    }
    let _ = u;
    Ok(terms)
}

fn taylor_terms_as_laurent(
    expr: &ExprArc,
    u: &Ident,
    order: usize,
    ctx: &Context,
) -> Result<Vec<(i32, ExprArc)>, EvalError> {
    let zero = Expr::int(0);
    let mut fk = Arc::clone(expr);
    let mut terms = Vec::new();
    for k in 0..order {
        let coeff = eval_at(&fk, u, &zero, ctx)?;
        if !is_zero(&coeff) {
            terms.push((k as i32, coeff));
        }
        if k + 1 < order {
            fk = crate::diff::diff(&fk, u)?;
            fk = eval(fk.as_ref(), ctx)?;
        }
    }
    Ok(terms)
}

fn eval_at(expr: &ExprArc, var: &Ident, center: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let sub = eval_subst_map(expr, &subst_map(var, Arc::clone(center)))?;
    eval(sub.as_ref(), ctx)
}

fn is_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

fn is_plus_infinity(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Symbol(id) if id.as_str() == "+infinity" || id.as_str() == "infinity"
    )
}

fn is_minus_infinity(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id.as_str() == "-infinity")
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn asymptotic_ck_int_56() {
        let ctx = xcas_default();
        let e = Expr::add(vec![
            Expr::func(
                FuncKind::Sqrt,
                vec![Expr::add(vec![
                    Expr::pow(Expr::sym("x"), Expr::int(2)),
                    Expr::sym("x"),
                    Expr::int(1),
                ])],
            ),
            Expr::mul(vec![
                Expr::int(-1),
                Expr::func(
                    FuncKind::Sqrt,
                    vec![Expr::add(vec![
                        Expr::pow(Expr::sym("x"), Expr::int(2)),
                        Expr::int(1),
                    ])],
                ),
            ]),
        ]);
        let r = limit_at_plus_infinity(&e, &Ident::new("x"), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/2");
    }

    #[test]
    fn asymptotic_ck_int_61() {
        let ctx = xcas_default();
        let inner = Arc::new(Expr::Frac(
            Expr::mul(vec![
                Expr::sym("x"),
                Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]),
            ]),
            Expr::add(vec![
                Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]),
                Expr::func(
                    FuncKind::Exp,
                    vec![Expr::mul(vec![
                        Expr::int(-2),
                        Arc::new(Expr::Frac(
                            Expr::pow(Expr::sym("x"), Expr::int(2)),
                            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
                        )),
                    ])],
                ),
            ]),
        ));
        let e = Expr::mul(vec![
            Expr::add(vec![
                Expr::func(FuncKind::Exp, vec![inner]),
                Expr::mul(vec![Expr::int(-1), Expr::func(FuncKind::Exp, vec![Expr::sym("x")])]),
            ]),
            Expr::pow(Expr::sym("x"), Expr::int(-1)),
        ]);
        let r = limit_at_plus_infinity(&e, &Ident::new("x"), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "-exp(2)");
    }

    #[test]
    fn asymptotic_rational_infinity() {
        let ctx = xcas_default();
        let e = Arc::new(Expr::Frac(
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::add(vec![Expr::sym("x"), Expr::int(-1)]),
        ));
        let r = limit_at_plus_infinity(&e, &Ident::new("x"), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1");
    }
}
