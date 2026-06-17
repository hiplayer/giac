use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, expr_to_poly, poly_to_expr, Context, EvalError, Expr, ExprArc,
    FuncKind, Ident,
};
use giac_simplify::{expand, ratnormal};
use giac_poly::{coeff_at, univariate_degree, Poly, Var};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::diff::diff;
use crate::integrate::try_as_rational;

use bounds::too_heavy_for_expand;

const MAX_LHOPITAL: usize = 8;

mod asymptotic;
mod bounds;
mod mrv;
mod mrv_lead_term;
mod mrv_series_lead;
mod mrv_w;
pub(crate) mod preprocess;
mod remove_lnexp;
mod sparse_series;

pub(crate) fn expr_has_nested_exp(e: &ExprArc) -> bool {
    bounds::expr_contains_nested_exp(e)
}

pub(crate) use mrv_series_lead::normalize_expr_quotients;
pub(crate) use asymptotic::{asymptotic_series_at_infinity, limit_at_plus_infinity};
pub(crate) use mrv_lead_term::{limit_from_mrv_lead_term, mrv_lead_term_plus_infinity};
pub(crate) use sparse_series::{series_at_center, series_at_zero, SparseSeries};

pub(crate) fn limit_finite_algebraic(
    expr: &ExprArc,
    var: &Ident,
    point: &ExprArc,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if let Ok(v) = subst_eval(expr, var, point, ctx) {
        if is_infinity(&v) {
            return Ok(v);
        }
        if !is_indeterminate(&v) && !contains_zero_negative_power(&v) {
            return Ok(v);
        }
    }
    if let Some((num, den)) = try_as_rational(expr, var) {
        if let Ok(r) = limit_rational_finite(&num, &den, var, point, ctx, 0) {
            return Ok(r);
        }
    }
    if !too_heavy_for_expand(expr) {
        let expanded = expand(expr.as_ref(), ctx)?;
        if let Ok(v) = subst_eval(&expanded, var, point, ctx) {
            if is_infinity(&v) {
                return Ok(v);
            }
            if !is_indeterminate(&v) && !contains_zero_negative_power(&v) {
                return Ok(v);
            }
        }
        if let Some((num, den)) = try_as_rational(&expanded, var) {
            if let Ok(r) = limit_rational_finite(&num, &den, var, point, ctx, 0) {
                return Ok(r);
            }
        }
    }
    let normalized = ratnormal(expr.as_ref(), ctx).ok();
    if let Some(normalized) = &normalized {
        if let Ok(v) = subst_eval(normalized, var, point, ctx) {
            if is_infinity(&v) {
                return Ok(v);
            }
            if !is_indeterminate(&v) && !contains_zero_negative_power(&v) {
                return Ok(v);
            }
        }
        if let Some((num, den)) = try_as_rational(normalized, var) {
            if let Ok(r) = limit_rational_finite(&num, &den, var, point, ctx, 0) {
                return Ok(r);
            }
        }
    }
    if let Some((num, den)) = try_as_quotient_add_shared_power(expr, var) {
        if let Ok(r) = limit_rational_finite(&num, &den, var, point, ctx, 0) {
            return Ok(r);
        }
    }
    if let Ok(r) = limit_quotient_finite(expr, var, point, ctx, 0) {
        return Ok(r);
    }
    Err(EvalError::NotImplemented("limit"))
}

fn try_as_quotient(expr: &ExprArc, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    match expr.as_ref() {
        Expr::Frac(num, den) => Some((Arc::clone(num), Arc::clone(den))),
        Expr::Mul(factors) => {
            let mut num = Vec::new();
            let mut den = Vec::new();
            for f in factors {
                if is_negative_integer_power(f) {
                    den.push(positive_integer_power(f)?);
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

fn is_negative_integer_power(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Pow(_, exp) if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
    )
}

fn positive_integer_power(e: &ExprArc) -> Option<ExprArc> {
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

fn is_negative_var_power(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Pow(b, exp)
            if is_var(b, var)
                && matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
    )
}

fn positive_var_power(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Pow(b, exp) if is_var(b, var) => {
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

fn limit_quotient_plus_infinity(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
    ctx: &Context,
    depth: usize,
) -> Result<ExprArc, EvalError> {
    if depth >= MAX_LHOPITAL {
        return Err(EvalError::NotImplemented("limit"));
    }
    let n_lim = limit_plus_infinity_quick(num, var, ctx);
    let d_lim = limit_plus_infinity_quick(den, var, ctx);
    let d_inf = d_lim
        .as_ref()
        .ok()
        .is_some_and(|d| is_plus_infinity(d) || is_minus_infinity(d));
    let n_inf = n_lim
        .as_ref()
        .ok()
        .is_some_and(|n| is_plus_infinity(n) || is_minus_infinity(n));
    if d_inf && (n_inf || n_lim.is_err()) {
        let dn = diff(num, var)?;
        let dd = diff(den, var)?;
        if is_one(&dd) {
            if matches!(num.as_ref(), Expr::Func(FuncKind::Exp, _))
                && expr_contains_exp(&dn)
            {
                return Ok(Expr::sym("+infinity"));
            }
            return limit_finite_algebraic(&dn, var, &Expr::int(0), ctx);
        }
        return limit_quotient_plus_infinity(&dn, &dd, var, ctx, depth + 1);
    }
    if let (Ok(n), Ok(d)) = (&n_lim, &d_lim) {
        if is_plus_infinity(n) && !is_plus_infinity(d) && !is_minus_infinity(d) {
            return Ok(Expr::sym("+infinity"));
        }
        if is_minus_infinity(n) && !is_plus_infinity(d) && !is_minus_infinity(d) {
            return Ok(Expr::sym("-infinity"));
        }
        if is_zero(n) && (is_plus_infinity(d) || is_minus_infinity(d)) {
            return Ok(Expr::int(0));
        }
    }
    Err(EvalError::NotImplemented("limit"))
}

fn limit_plus_infinity_quick(expr: &ExprArc, var: &Ident, ctx: &Context) -> Result<ExprArc, EvalError> {
    if let Ok(r) = limit_at_plus_infinity(expr, var, ctx) {
        return Ok(r);
    }
    if let Some((num, den)) = try_as_rational(expr, var) {
        return limit_rational_infinity(&num, &den, var, ctx);
    }
    if is_var(expr, var) {
        return Ok(Expr::sym("+infinity"));
    }
    if let Expr::Func(FuncKind::Exp, args) = expr.as_ref() {
        if args.len() == 1 {
            if let Ok(u) = limit_plus_infinity_quick(&args[0], var, ctx) {
                if is_plus_infinity(&u) {
                    return Ok(Expr::sym("+infinity"));
                }
                if is_minus_infinity(&u) {
                    return Ok(Expr::int(0));
                }
            }
        }
    }
    Err(EvalError::NotImplemented("limit"))
}

fn try_as_quotient_add_shared_power(expr: &ExprArc, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    let Expr::Add(terms) = expr.as_ref() else {
        return None;
    };
    let mut num_terms = Vec::with_capacity(terms.len());
    let mut den = None;
    for t in terms {
        let (n, d) = try_as_quotient(t, var)?;
        match &den {
            None => den = Some(d),
            Some(dd) if dd == &d => {}
            _ => return None,
        }
        num_terms.push(n);
    }
    Some((Expr::add(num_terms), den?))
}

fn limit_quotient_finite(
    expr: &ExprArc,
    var: &Ident,
    point: &ExprArc,
    ctx: &Context,
    depth: usize,
) -> Result<ExprArc, EvalError> {
    let (num, den) = match try_as_quotient(expr, var) {
        Some(q) => q,
        None => return Err(EvalError::NotImplemented("limit")),
    };
    limit_rational_finite(&num, &den, var, point, ctx, depth)
}

pub(crate) fn limit_plus_infinity_algebraic(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if let Some((num, den)) = try_as_rational(expr, var) {
        if let Ok(r) = limit_rational_infinity(&num, &den, var, ctx) {
            return Ok(r);
        }
    }
    if let Ok(r) = limit_at_plus_infinity(expr, var, ctx) {
        return Ok(r);
    }
    if let Some((num, den)) = try_as_quotient(expr, var) {
        if matches!(num.as_ref(), Expr::Func(FuncKind::Exp, _)) {
            if let Ok(r) = limit_quotient_plus_infinity(&num, &den, var, ctx, 0) {
                return Ok(r);
            }
        }
    }
    let normalized = ratnormal(expr.as_ref(), ctx).ok();
    if let Some(normalized) = &normalized {
        if let Some((num, den)) = try_as_rational(normalized, var) {
            if let Ok(r) = limit_rational_infinity(&num, &den, var, ctx) {
                return Ok(r);
            }
        }
        if let Some((num, den)) = try_as_quotient(normalized, var) {
            if matches!(num.as_ref(), Expr::Func(FuncKind::Exp, _)) {
                if let Ok(r) = limit_quotient_plus_infinity(&num, &den, var, ctx, 0) {
                    return Ok(r);
                }
            }
        }
    }
    if let Some((num, den)) = try_as_quotient(expr, var) {
        if matches!(num.as_ref(), Expr::Add(_)) && expr_contains_exp(&num) {
            if let Ok(r) = limit_via_reciprocal(expr, var, ctx) {
                return Ok(r);
            }
        }
    }
    limit_via_reciprocal(expr, var, ctx)
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

pub(crate) fn limit_minus_infinity_algebraic(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let neg = Expr::mul(vec![Expr::int(-1), var_to_expr(var)]);
    let swapped = eval_subst_map(expr, &subst_map(var, neg))?;
    limit_plus_infinity_algebraic(&swapped, var, ctx)
}

fn limit_rational_finite(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
    point: &ExprArc,
    ctx: &Context,
    depth: usize,
) -> Result<ExprArc, EvalError> {
    let combined = Arc::new(Expr::Frac(Arc::clone(num), Arc::clone(den)));
    let reduced = cancel_rational_pole(combined, var, point, ctx)?;
    if let Ok(v) = subst_eval(&reduced, var, point, ctx) {
        if is_infinity(&v) {
            return Ok(v);
        }
        if !is_indeterminate(&v) && !contains_zero_negative_power(&v) {
            return Ok(v);
        }
    }
    if depth >= MAX_LHOPITAL {
        return Err(EvalError::NotImplemented("limit"));
    }
    let (n2, d2) = match reduced.as_ref() {
        Expr::Frac(n, d) => (n, d),
        _ => return Err(EvalError::NotImplemented("limit")),
    };
    let n_val = subst_eval(n2, var, point, ctx);
    let d_val = subst_eval(d2, var, point, ctx);
    match (&n_val, &d_val) {
        (Ok(n), Ok(d)) if is_zero(d) && is_zero(n) => {
            let dn = diff(n2, var)?;
            let dd = diff(d2, var)?;
            limit_rational_finite(&dn, &dd, var, point, ctx, depth + 1)
        }
        (Ok(_), Ok(d)) if is_zero(d) => pole_infinity(n2, d2, var, point, ctx),
        (Ok(n), Ok(d)) if !is_zero(d) => {
            if contains_zero_negative_power(n) && depth < MAX_LHOPITAL {
                let dn = diff(n2, var)?;
                let dd = diff(d2, var)?;
                return limit_rational_finite(&dn, &dd, var, point, ctx, depth + 1);
            }
            let frac = Arc::new(Expr::Frac(Arc::clone(n), Arc::clone(d)));
            eval(frac.as_ref(), ctx)
        }
        _ if depth < MAX_LHOPITAL => {
            let dn = diff(n2, var)?;
            let dd = diff(d2, var)?;
            limit_rational_finite(&dn, &dd, var, point, ctx, depth + 1)
        }
        _ => Err(EvalError::NotImplemented("limit")),
    }
}

fn cancel_rational_pole(
    frac: ExprArc,
    var: &Ident,
    point: &ExprArc,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let (num, den) = match frac.as_ref() {
        Expr::Frac(n, d) => (n, d),
        _ => return Ok(frac),
    };
    let v = Var::from(var.as_str());
    let num_p = match expr_to_poly(num) {
        Ok(p) => p,
        Err(_) => return Ok(frac),
    };
    let den_p = match expr_to_poly(den) {
        Ok(p) => p,
        Err(_) => return Ok(frac),
    };
    let pt = eval(point.as_ref(), ctx)?;
    let pt_p = expr_to_poly(&pt)?;
    if univariate_degree(&pt_p, &v) != 0 {
        return Ok(ratnormal(frac.as_ref(), ctx).unwrap_or(frac));
    }
    let a = coeff_at(&pt_p, &v, 0);
    let mut out_num = num_p.clone();
    let mut out_den = den_p.clone();
    loop {
        let val = out_den.horner(&v, &a);
        if !val.is_zero() {
            break;
        }
        let val_num = out_num.horner(&v, &a);
        if !val_num.is_zero() {
            break;
        }
        let lin = Poly::var(v.clone()).add(&Poly::constant(-a.clone()));
        let (_, r) = out_den.div_rem(&lin);
        if !r.is_zero() {
            break;
        }
        out_den = out_den.div_rem(&lin).0;
        out_num = out_num.div_rem(&lin).0;
    }
    Ok(Arc::new(Expr::Frac(
        poly_to_expr(&out_num),
        poly_to_expr(&out_den),
    )))
}

fn limit_rational_infinity(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
    _ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let v = Var::from(var.as_str());
    let num_p = expr_to_poly(num).map_err(|_| EvalError::NotImplemented("limit"))?;
    let den_p = expr_to_poly(den).map_err(|_| EvalError::NotImplemented("limit"))?;
    let nd = univariate_degree(&num_p, &v);
    let dd = univariate_degree(&den_p, &v);
    if nd < dd {
        return Ok(Expr::int(0));
    }
    if nd > dd {
        return Ok(sign_infinity(leading_ratio(&num_p, &den_p, &v)));
    }
    let ratio = leading_ratio(&num_p, &den_p, &v);
    if ratio.is_zero() {
        return Ok(Expr::int(0));
    }
    if nd == dd {
        return Ok(ratio_to_expr(&ratio));
    }
    Ok(sign_infinity(ratio))
}

fn leading_ratio(num: &Poly, den: &Poly, var: &Var) -> Ratio<BigInt> {
    let nd = univariate_degree(num, var);
    let dd = univariate_degree(den, var);
    coeff_at(num, var, nd) / coeff_at(den, var, dd)
}

fn sign_infinity(r: Ratio<BigInt>) -> ExprArc {
    if r.is_negative() {
        Expr::sym("-infinity")
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

fn pole_infinity(
    _num: &ExprArc,
    _den: &ExprArc,
    _var: &Ident,
    _point: &ExprArc,
    _ctx: &Context,
) -> Result<ExprArc, EvalError> {
    Ok(Expr::sym("+infinity"))
}

fn limit_via_reciprocal(expr: &ExprArc, var: &Ident, ctx: &Context) -> Result<ExprArc, EvalError> {
    if too_heavy_for_expand(expr) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let t = Ident::new("_limit_t");
    let inv = Expr::pow(var_to_expr(&t), Expr::int(-1));
    let swapped = eval_subst_map(expr, &subst_map(var, inv))?;
    let normalized = ratnormal(swapped.as_ref(), ctx).unwrap_or(swapped);
    let zero = Expr::int(0);
    if let Ok(r) = limit_quotient_finite(&normalized, &t, &zero, ctx, 0) {
        if !contains_zero_negative_power(&r) {
            return Ok(r);
        }
    }
    if let Some((num, den)) = try_as_quotient_add_shared_power(&normalized, &t) {
        if let Ok(r) = limit_rational_finite(&num, &den, &t, &zero, ctx, 0) {
            if !contains_zero_negative_power(&r) {
                return Ok(r);
            }
        }
    }
    limit_finite_algebraic(&normalized, &t, &zero, ctx)
}

fn subst_eval(expr: &ExprArc, var: &Ident, point: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let sub = eval_subst_map(expr, &subst_map(var, Arc::clone(point)))?;
    eval(sub.as_ref(), ctx)
}

fn expr_contains_exp(expr: &ExprArc) -> bool {
    match expr.as_ref() {
        Expr::Func(FuncKind::Exp, _) => true,
        Expr::Add(ts) => ts.iter().any(expr_contains_exp),
        Expr::Mul(fs) => fs.iter().any(expr_contains_exp),
        Expr::Pow(b, e) => expr_contains_exp(b) || expr_contains_exp(e),
        Expr::Frac(n, d) => expr_contains_exp(n) || expr_contains_exp(d),
        Expr::Func(_, args) => args.iter().any(expr_contains_exp),
        _ => false,
    }
}

fn subst_map(var: &Ident, value: ExprArc) -> HashMap<Ident, ExprArc> {
    let mut m = HashMap::new();
    m.insert(var.clone(), value);
    m
}

fn var_to_expr(var: &Ident) -> ExprArc {
    Expr::sym(var.as_str())
}

fn is_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

fn is_one(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n == &BigInt::from(1))
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

fn is_infinity(e: &ExprArc) -> bool {
    is_plus_infinity(e) || is_minus_infinity(e)
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

fn is_indeterminate(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Frac(num, den) => is_zero(num) && is_zero(den),
        Expr::Int(_) => false,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    fn limit_line(expr: ExprArc, var: &str, pt: ExprArc) -> Result<ExprArc, EvalError> {
        let ctx = xcas_default();
        let v = Ident::new(var);
        let point = classify_test_point(&pt)?;
        match point {
            TestPoint::Finite => limit_finite_algebraic(&expr, &v, &pt, &ctx),
            TestPoint::PlusInfinity => limit_plus_infinity_algebraic(&expr, &v, &ctx),
            TestPoint::MinusInfinity => limit_minus_infinity_algebraic(&expr, &v, &ctx),
        }
    }

    enum TestPoint {
        Finite,
        PlusInfinity,
        MinusInfinity,
    }

    fn classify_test_point(e: &ExprArc) -> Result<TestPoint, EvalError> {
        match e.as_ref() {
            Expr::Symbol(id) if id.as_str() == "+infinity" || id.as_str() == "infinity" => {
                Ok(TestPoint::PlusInfinity)
            }
            Expr::Mul(factors)
                if factors.len() == 2
                    && matches!(factors[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                    && matches!(factors[1].as_ref(), Expr::Symbol(id) if id.as_str() == "infinity") =>
            {
                Ok(TestPoint::MinusInfinity)
            }
            _ => Ok(TestPoint::Finite),
        }
    }

    #[test]
    fn engine_ck_int_56() {
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
        let r = limit_line(e, "x", Expr::sym("+infinity")).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/2");
    }

    #[test]
    fn engine_ck_int_58() {
        let e = Arc::new(Expr::Frac(
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::func(
                FuncKind::Sqrt,
                vec![Arc::new(Expr::Frac(
                    Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
                    Expr::add(vec![Expr::sym("x"), Expr::int(-1)]),
                ))],
            ),
        ));
        let r = limit_line(e, "x", Expr::sym("+infinity")).unwrap();
        assert_eq!(format_expr(r.as_ref()), "+infinity");
    }

    #[test]
    fn engine_ck_int_59_via_finite() {
        let e = Arc::new(Expr::Frac(
            Expr::add(vec![Expr::int(1), Expr::mul(vec![Expr::int(-2), Expr::sym("x")])]),
            Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::sym("x"),
                Expr::int(-2),
            ]),
        ));
        let r = limit_line(e, "x", Expr::int(1)).unwrap();
        assert_eq!(format_expr(r.as_ref()), "+infinity");
    }

    #[test]
    fn engine_ck_int_60() {
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
            Expr::func(FuncKind::Exp, vec![inner]),
            Expr::pow(Expr::sym("x"), Expr::int(-1)),
        ]);
        let r = limit_line(e, "x", Expr::sym("+infinity")).unwrap();
        assert_eq!(format_expr(r.as_ref()), "+infinity");
    }

    #[test]
    fn engine_ck_int_61() {
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
        let exp_inner = Expr::func(FuncKind::Exp, vec![inner]);
        let e = Expr::mul(vec![
            Expr::add(vec![
                exp_inner.clone(),
                Expr::mul(vec![Expr::int(-1), Expr::func(FuncKind::Exp, vec![Expr::sym("x")])]),
            ]),
            Expr::pow(Expr::sym("x"), Expr::int(-1)),
        ]);
        let r = limit_line(e, "x", Expr::sym("+infinity")).unwrap();
        let s = format_expr(r.as_ref());
        assert_eq!(s, "-exp(2)");
    }
}
