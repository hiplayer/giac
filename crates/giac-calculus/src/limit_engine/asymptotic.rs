//! Asymptotic expansion and limits at `+infinity` (GIAC-216 / `series.cc` subset).
//!
//! giac uses `mrv_lead_term` for full asymptotics; here we implement a practical
//! subset via reciprocal substitution `x = 1/u` and Laurent analysis at `u = 0`.

use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{bigint_to_i64, eval, eval_subst_map, expr_to_poly, Context, EvalError, Expr,
    ExprArc, FuncKind, Ident,
};
use giac_simplify::{normal, ratnormal};
use giac_poly::{coeff_at, univariate_degree, Poly, Var};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use super::bounds::{mrv_limit_eligible, too_heavy_for_expand, MAX_SERIES_EXPANSION_ORDER};
use super::exp_diff::{
    classify_exp_at_plus_infinity, classify_signed_exp_at_plus_infinity, simplify_add_sum,
    unwrap_signed_frac,
};
use super::mrv::{try_const_f64, vanishes_faster_than_at_plus_infinity};
use super::mrv_lead_term::limit_unidirectional_plus_infinity;
use super::mrv_series_lead::try_as_quotient;
use super::preprocess::{limit_preprocess_plus_infinity, limit_preprocess_struct};
use super::sparse_series::series_at_zero_order;

use crate::integrate::try_as_rational;
use crate::expr_util::depends_on_var;

const ASYM_U: &str = "_asym_u";
const MAX_PUMP: i64 = 12;
const DEFAULT_SERIES_ORDER: usize = 8;

/// Limit as `var → +infinity`: preprocess → algebraic lead → MRV → fallback.
pub(crate) fn limit_at_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let pre = limit_preprocess_struct(expr, var);
    if let Some(r) = limit_preprocessed_at_plus_infinity(&pre, var, ctx) {
        return Ok(normalize_limit_result(&r, ctx));
    }
    if let Some(r) = limit_var_over_x_pow_ln(expr, var)
        .or_else(|| limit_exp_sum_nth_root(expr, var))
        .or_else(|| limit_poly_over_sqrt_at_infinity(expr, var, ctx))
    {
        return Ok(normalize_limit_result(&r, ctx));
    }
    limit_at_plus_infinity_fallback(expr, var, ctx)
        .map(|r| normalize_limit_result(&r, ctx))
}

/// After `limit_preprocess_struct`: classify reduced forms before MRV.
fn limit_preprocessed_at_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    limit_add_at_plus_infinity(expr, var, ctx)
        .or_else(|| classify_signed_exp_at_plus_infinity(expr, var, ctx))
        .or_else(|| {
            if let Expr::Func(FuncKind::Exp, args) = expr.as_ref() {
                if args.len() == 1 {
                    return classify_exp_at_plus_infinity(&args[0], var, ctx);
                }
            }
            None
        })
        .or_else(|| limit_exp_of_vanishing_frac_argument(expr, var, ctx))
        .or_else(|| limit_exp_times_frac_quotient_at_plus_infinity(expr, var, ctx))
        .or_else(|| {
            if !mrv_limit_eligible(expr) {
                return None;
            }
            limit_unidirectional_plus_infinity(expr, var, ctx)
                .ok()
                .filter(|r| is_usable_limit(r))
        })
}

/// `exp(N/D)` with `N/D → 0` at `+∞` after balance → `1`.
fn limit_exp_of_vanishing_frac_argument(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let arg = match expr.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => &args[0],
        _ => return None,
    };
    if let Some((n, d)) = unwrap_signed_frac(arg) {
        if vanishes_faster_than_at_plus_infinity(&n, &d, var, ctx) {
            return Some(Expr::int(1));
        }
    }
    let lim = limit_at_plus_infinity_fallback(arg, var, ctx).ok()?;
    if matches!(lim.as_ref(), Expr::Int(n) if n.is_zero()) {
        return Some(Expr::int(1));
    }
    None
}

/// `±exp(L)·n/d` with `exp(L)·n/d → exp(c)` at `+∞` (CK-INT-61 after preprocess).
fn limit_exp_times_frac_quotient_at_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let (neg, scale, n, d) = parse_signed_exp_frac_product(expr, var)?;
    let num = Expr::mul(vec![Expr::func(FuncKind::Exp, vec![scale]), n]);
    let num_ln = dominant_ln_exponent_at_plus_infinity(&num, var, ctx)?;
    let den_ln = dominant_ln_exponent_at_plus_infinity(&d, var, ctx)?;
    let diff = simplify_add_sum(&Expr::add(vec![
        num_ln,
        Expr::mul(vec![Expr::int(-1), den_ln]),
    ]));
    let c = limit_const_rational_at_plus_infinity(&diff, var, ctx)?;
    let mut out = Expr::func(FuncKind::Exp, vec![float_to_expr(c)?]);
    if neg {
        out = Expr::mul(vec![Expr::int(-1), out]);
    }
    eval(out.as_ref(), ctx).ok()
}

pub(crate) fn parse_signed_exp_frac_product(
    expr: &ExprArc,
    var: &Ident,
) -> Option<(bool, ExprArc, ExprArc, ExprArc)> {
    let flat = flatten_top_mul(expr);
    let mut neg = false;
    let mut exp_logs = Vec::new();
    let mut frac = None;
    let mut var_pow = 0isize;
    for f in &flat {
        match f.as_ref() {
            Expr::Int(n) if n.is_negative() => neg = !neg,
            Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
                exp_logs.push(Arc::clone(&args[0]));
            }
            Expr::Frac(n, d) => {
                if frac.is_some() {
                    return None;
                }
                frac = Some((Arc::clone(n), Arc::clone(d)));
            }
            Expr::Pow(b, e) if is_limit_var(b, var) => {
                let k = match e.as_ref() {
                    Expr::Int(n) => bigint_to_i64(n).ok()? as isize,
                    _ => return None,
                };
                var_pow += k;
            }
            Expr::Symbol(id) if id == var => var_pow += 1,
            _ => return None,
        }
    }
    let (mut n, d) = frac?;
    let (n_neg, n_body) = extract_mul_sign(&n);
    neg ^= n_neg;
    n = n_body;
    if var_pow != 0 {
        n = Expr::mul(vec![
            n,
            Expr::pow(Expr::sym(var.as_str()), Expr::int(var_pow as i64)),
        ]);
        n = cancel_var_power_in_mul(&n, var);
    }
    let scale = match exp_logs.len() {
        0 => return None,
        1 => exp_logs.into_iter().next().unwrap(),
        _ => Expr::add(exp_logs),
    };
    Some((neg, scale, n, d))
}

fn cancel_var_power_in_mul(e: &ExprArc, var: &Ident) -> ExprArc {
    let flat = flatten_top_mul(e);
    let mut net: isize = 0;
    let mut rest = Vec::new();
    for f in &flat {
        if is_limit_var(f, var) {
            net += 1;
        } else if let Expr::Pow(b, exp) = f.as_ref() {
            if is_limit_var(b, var) {
                if let Expr::Int(n) = exp.as_ref() {
                    if let Ok(k) = bigint_to_i64(n) {
                        net += k as isize;
                        continue;
                    }
                }
            }
            rest.push(Arc::clone(f));
        } else {
            rest.push(Arc::clone(f));
        }
    }
    if net != 0 {
        rest.push(Expr::pow(Expr::sym(var.as_str()), Expr::int(net as i64)));
    }
    match rest.len() {
        0 => Expr::int(1),
        1 => rest.into_iter().next().unwrap(),
        _ => Expr::mul(rest),
    }
}

fn extract_mul_sign(e: &ExprArc) -> (bool, ExprArc) {
    match e.as_ref() {
        Expr::Mul(fs) => {
            let mut neg = false;
            let mut rest = Vec::new();
            for f in fs {
                if matches!(f.as_ref(), Expr::Int(n) if n.is_negative()) {
                    neg = !neg;
                } else {
                    rest.push(Arc::clone(f));
                }
            }
            let body = match rest.len() {
                0 => Expr::int(1),
                1 => rest.into_iter().next().unwrap(),
                _ => Expr::mul(rest),
            };
            (neg, body)
        }
        Expr::Int(n) if n.is_negative() => (true, Expr::int(1)),
        _ => (false, Arc::clone(e)),
    }
}

fn flatten_top_mul(expr: &ExprArc) -> Vec<ExprArc> {
    match expr.as_ref() {
        Expr::Mul(fs) => fs.iter().flat_map(flatten_top_mul).collect(),
        _ => vec![Arc::clone(expr)],
    }
}

fn is_limit_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

pub(crate) fn dominant_ln_exponent_at_plus_infinity(
    e: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => Some(Arc::clone(&args[0])),
        Expr::Mul(fs) => {
            let mut parts = Vec::new();
            for f in fs {
                if let Expr::Func(FuncKind::Exp, args) = f.as_ref() {
                    if args.len() == 1 {
                        parts.push(Arc::clone(&args[0]));
                        continue;
                    }
                }
                if let Some(p) = dominant_ln_exponent_at_plus_infinity(f, var, ctx) {
                    parts.push(p);
                } else if depends_on_var(f, var) {
                    return None;
                }
            }
            if parts.is_empty() {
                None
            } else if parts.len() == 1 {
                Some(parts.into_iter().next().unwrap())
            } else {
                Some(Expr::add(parts))
            }
        }
        Expr::Add(ts) => ts
            .iter()
            .filter_map(|t| dominant_ln_exponent_at_plus_infinity(t, var, ctx))
            .max_by(|a, b| {
                super::mrv::mrv_compare(a, b, var, ctx)
            }),
        _ => None,
    }
}

fn limit_const_rational_at_plus_infinity(
    e: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<f64> {
    if let Some((n, d)) = try_add_rational_to_frac(e, var) {
        let n = normal(n.as_ref(), ctx).ok()?;
        return super::mrv::limit_rational_const_at_plus_infinity(&n, &d, var);
    }
    match e.as_ref() {
        Expr::Frac(n, d) => super::mrv::limit_rational_const_at_plus_infinity(n, d, var),
        Expr::Add(ts) => {
            let mut sum = 0.0;
            let mut any = false;
            for t in ts {
                if let Some(v) = limit_const_rational_at_plus_infinity(t, var, ctx) {
                    sum += v;
                    any = true;
                } else if depends_on_var(t, var) {
                    return None;
                }
            }
            any.then_some(sum)
        }
        _ => super::mrv::try_const_f64(e),
    }
}

fn try_add_rational_to_frac(e: &ExprArc, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    let Expr::Add(ts) = e.as_ref() else {
        return None;
    };
    let mut common_den = None;
    for t in ts {
        if let Some((_, d)) = as_frac_form(t, var) {
            if common_den.is_none() {
                common_den = Some(d);
            } else if common_den.as_ref() != Some(&d) {
                return None;
            }
        }
    }
    let d = common_den?;
    let mut lin_coeff = Expr::int(0);
    let mut num = Expr::int(0);
    for t in ts {
        if let Some((n, df)) = as_frac_form(t, var) {
            if df.as_ref() != d.as_ref() {
                return None;
            }
            num = Expr::add(vec![num, n]);
        } else if let Some(c) = linear_var_term(t, var) {
            lin_coeff = Expr::add(vec![lin_coeff, c]);
        } else {
            return None;
        }
    }
    if !matches!(lin_coeff.as_ref(), Expr::Int(n) if n.is_zero()) {
        num = Expr::add(vec![
            num,
            Expr::mul(vec![
                lin_coeff,
                Expr::sym(var.as_str()),
                Arc::clone(&d),
            ]),
        ]);
    }
    Some((num, d))
}

fn split_mul_leading_coeff(e: &ExprArc) -> (ExprArc, ExprArc) {
    let flat = flatten_top_mul(e);
    let mut coeff = Expr::int(1);
    let mut rest = Vec::new();
    for f in &flat {
        match f.as_ref() {
            Expr::Int(_) | Expr::Rat(_) => coeff = Expr::mul(vec![coeff, Arc::clone(f)]),
            _ => rest.push(Arc::clone(f)),
        }
    }
    let body = match rest.len() {
        0 => Expr::int(1),
        1 => rest.into_iter().next().unwrap(),
        _ => Expr::mul(rest),
    };
    (coeff, body)
}

fn as_frac_form(e: &ExprArc, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    let _ = var;
    if let Expr::Frac(n, d) = e.as_ref() {
        return Some((Arc::clone(n), Arc::clone(d)));
    }
    let (coeff, body) = split_mul_leading_coeff(e);
    if let Expr::Frac(n, d) = body.as_ref() {
        return Some((Expr::mul(vec![coeff, Arc::clone(n)]), Arc::clone(d)));
    }
    let flat = flatten_top_mul(&body);
    let mut inv_den = None;
    let mut num_parts = Vec::new();
    for f in &flat {
        if let Expr::Pow(b, exp) = f.as_ref() {
            if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) {
                inv_den = Some(Arc::clone(b));
                continue;
            }
        }
        num_parts.push(Arc::clone(f));
    }
    let d = inv_den?;
    let n_body = match num_parts.len() {
        0 => Expr::int(1),
        1 => num_parts.into_iter().next().unwrap(),
        _ => Expr::mul(num_parts),
    };
    Some((Expr::mul(vec![coeff, n_body]), d))
}

fn linear_var_term(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let flat = flatten_top_mul(e);
    if !flat.iter().any(|f| is_limit_var(f, var)) {
        return None;
    }
    let mut coeff = Expr::int(1);
    for f in &flat {
        if is_limit_var(f, var) {
            continue;
        }
        match f.as_ref() {
            Expr::Int(_) | Expr::Rat(_) => coeff = Expr::mul(vec![coeff, Arc::clone(f)]),
            _ => return None,
        }
    }
    Some(coeff)
}

fn float_to_expr(c: f64) -> Option<ExprArc> {
    if (c - c.round()).abs() < f64::EPSILON {
        Some(Expr::int(c as i64))
    } else {
        None
    }
}

fn limit_add_at_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let Expr::Add(ts) = expr.as_ref() else {
        return None;
    };
    let mut sum = Expr::int(0);
    for t in ts {
        let lim = limit_term_at_plus_infinity(t, var, ctx)?;
        sum = Expr::add(vec![sum, lim]);
    }
    eval(sum.as_ref(), ctx).ok()
}

fn limit_term_at_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    classify_signed_exp_at_plus_infinity(expr, var, ctx).or_else(|| {
        if let Expr::Func(FuncKind::Exp, args) = expr.as_ref() {
            if args.len() == 1 {
                return classify_exp_at_plus_infinity(&args[0], var, ctx);
            }
        }
        None
    }).or_else(|| {
        if !mrv_limit_eligible(expr) {
            return None;
        }
        limit_unidirectional_plus_infinity(expr, var, ctx)
            .ok()
            .filter(|r| is_usable_limit(r))
    })
}

/// `x = 1/u` then limit at `u = 0` (upstream finite-point substitution for `+infinity`).
pub(crate) fn limit_at_plus_infinity_fallback(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if too_heavy_for_expand(expr) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let pre = limit_preprocess_plus_infinity(expr, var, ctx)?;
    if let Some((num, den)) = try_as_rational(&pre, var) {
        if let Some(r) = limit_rational_leading_at_infinity(&num, &den, var) {
            return Ok(normalize_limit_result(&r, ctx));
        }
    }
    if let Some(r) = limit_rational_over_sqrt_quotient_at_infinity(&pre, var)
        .or_else(|| limit_poly_over_sqrt_at_infinity(&pre, var, ctx))
        .or_else(|| limit_sqrt_sum_quotient_at_infinity(&pre, var))
    {
        return Ok(normalize_limit_result(&r, ctx));
    }
    let u = Ident::new(ASYM_U);
    let swapped = reciprocal_subst(&pre, var, &u)?;
    let swapped = peel_shared_u_inv_in_frac(&swapped, &u);
    limit_at_zero_fallback(&swapped, &u, ctx)
}

pub(crate) fn limit_at_zero_fallback(
    expr: &ExprArc,
    u: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let rationalized = ratnormal(expr.as_ref(), ctx).unwrap_or_else(|_| Arc::clone(expr));
    if let Some(r) = limit_from_rational_laurent(&rationalized, u) {
        if is_usable_limit(&r) {
            return Ok(normalize_limit_result(&r, ctx));
        }
    }
    if let Some(r) = limit_at_zero_from_series_escalating(&rationalized, u, ctx) {
        if is_usable_limit(&r) {
            return Ok(normalize_limit_result(&r, ctx));
        }
    }
    if let Some(r) = limit_at_zero_rational_lead_escalating(&rationalized, u, ctx) {
        if is_usable_limit(&r) {
            return Ok(normalize_limit_result(&r, ctx));
        }
    }
    if let Some(r) = limit_from_fractional_u_valuation(&rationalized, u, ctx) {
        if is_usable_limit(&r) {
            return Ok(normalize_limit_result(&r, ctx));
        }
    }
    limit_from_scaled_finite(&rationalized, u, ctx)
        .map(|r| normalize_limit_result(&r, ctx))
}

/// `num(u)/den(u)` at `u=0` when `den(0) != 0` and `num` has a series lead term.
fn limit_at_zero_rational_lead(expr: &ExprArc, u: &Ident, order: usize, ctx: &Context) -> Option<ExprArc> {
    let (num, den) = try_as_rational(expr, u)?;
    let zero = Expr::int(0);
    let den0 = eval_at(&den, u, &zero, ctx).ok()?;
    let den0 = collapse_unit_powers(&den0);
    let den0 = eval(den0.as_ref(), ctx).ok()?;
    if is_zero(&den0) || is_indeterminate(&den0) {
        return None;
    }
    let s = series_at_zero_order(&num, u, order, MAX_SERIES_EXPANSION_ORDER, ctx).ok()?;
    let (exp, coeff) = s.lead()?;
    if exp != 0 {
        return None;
    }
    let quot = Arc::new(Expr::Frac(coeff, den0));
    eval(collapse_unit_powers(&quot).as_ref(), ctx).ok()
}

fn collapse_unit_powers(e: &ExprArc) -> ExprArc {
    match e.as_ref() {
        Expr::Frac(n, d) if matches!(d.as_ref(), Expr::Int(n) if n.is_one()) => {
            collapse_unit_powers(n)
        }
        Expr::Pow(b, exp) if is_half_exponent(exp) && matches!(b.as_ref(), Expr::Int(n) if n.is_one()) => {
            Expr::int(1)
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(collapse_unit_powers).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(collapse_unit_powers).collect()),
        Expr::Pow(b, exp) => Expr::pow(collapse_unit_powers(b), collapse_unit_powers(exp)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(collapse_unit_powers(n), collapse_unit_powers(d))),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(collapse_unit_powers).collect()),
        _ => Arc::clone(e),
    }
}

fn limit_at_zero_rational_lead_escalating(expr: &ExprArc, u: &Ident, ctx: &Context) -> Option<ExprArc> {
    series_ordre_escalation(|order| limit_at_zero_rational_lead(expr, u, order, ctx))
}

/// Leading term of `series_at_zero` as a limit at `u = 0`.
fn limit_at_zero_from_series(expr: &ExprArc, u: &Ident, order: usize, ctx: &Context) -> Option<ExprArc> {
    let s = series_at_zero_order(expr, u, order, MAX_SERIES_EXPANSION_ORDER, ctx).ok()?;
    let (exp, coeff) = s.lead()?;
    if exp > 0 {
        return Some(Expr::int(0));
    }
    if exp < 0 {
        return Some(sign_infinity_from_value(&eval(coeff.as_ref(), ctx).ok()?));
    }
    eval(coeff.as_ref(), ctx).ok()
}

fn limit_at_zero_from_series_escalating(expr: &ExprArc, u: &Ident, ctx: &Context) -> Option<ExprArc> {
    series_ordre_escalation(|order| limit_at_zero_from_series(expr, u, order, ctx))
}

fn series_ordre_escalation(mut try_order: impl FnMut(usize) -> Option<ExprArc>) -> Option<ExprArc> {
    let mut ordre = DEFAULT_SERIES_ORDER as f64;
    let cap = MAX_SERIES_EXPANSION_ORDER as f64;
    while ordre < cap {
        let try_ord = (ordre as usize).min(MAX_SERIES_EXPANSION_ORDER);
        if let Some(r) = try_order(try_ord) {
            return Some(r);
        }
        ordre = ordre * 1.5 + 1.0;
    }
    None
}

/// Asymptotic series in `1/var` up to `order` terms (GIAC-216d).
pub(crate) fn asymptotic_series_at_infinity(
    expr: &ExprArc,
    var: &Ident,
    order: usize,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if too_heavy_for_expand(expr) {
        return Err(EvalError::NotImplemented("series"));
    }
    let u = Ident::new(ASYM_U);
    let swapped = reciprocal_subst(expr, var, &u)?;
    let normalized = ratnormal(swapped.as_ref(), ctx).unwrap_or(swapped);
    let series = series_at_zero_order(&normalized, &u, order, MAX_SERIES_EXPANSION_ORDER, ctx)?;
    if series.is_zero() {
        return Ok(Expr::int(0));
    }
    let inv = Expr::pow(var_to_expr(var), Expr::int(-1));
    let mut out = Vec::new();
    for (exp, coeff) in series.iter_terms() {
        let scaled = if exp == 0 {
            Arc::clone(coeff)
        } else if exp > 0 {
            Expr::mul(vec![
                Arc::clone(coeff),
                Expr::pow(Arc::clone(&inv), Expr::int(i64::from(exp))),
            ])
        } else {
            Expr::mul(vec![
                Arc::clone(coeff),
                Expr::pow(var_to_expr(&u), Expr::int(i64::from(-exp))),
            ])
        };
        out.push(scaled);
    }
    eval(Expr::add(out).as_ref(), ctx)
}

/// Polynomial degree ratio at `+∞` (equivalent to Laurent after reciprocal; not a shape table).
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

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_sqrt(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Func(giac_core::FuncKind::Sqrt, args) if args.len() == 1)
        || matches!(e.as_ref(), Expr::Pow(b, exp) if is_half_exponent(exp) && !is_sqrt(b))
        || matches!(
            e.as_ref(),
            Expr::Pow(b, exp)
                if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(1)) && is_sqrt(b)
        )
}

fn sqrt_arg(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Func(giac_core::FuncKind::Sqrt, args) if args.len() == 1 => {
            Some(Arc::clone(&args[0]))
        }
        Expr::Pow(b, exp) if is_half_exponent(exp) => Some(Arc::clone(b)),
        Expr::Pow(b, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(1)) =>
        {
            sqrt_arg(b)
        }
        _ => None,
    }
}

fn limit_rational_over_sqrt_quotient_at_infinity(expr: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (num, den) = try_as_quotient(expr)?;
    let inner = sqrt_arg(&den)?;
    let (a, b) = match inner.as_ref() {
        Expr::Frac(n, d) => (Arc::clone(n), Arc::clone(d)),
        _ => try_as_quotient(&inner)?,
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

/// After `normalize_sqrt_conjugates`, `poly/(sqrt+…)` leading term at `+∞`.
fn limit_sqrt_sum_quotient_at_infinity(expr: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (num, den) = try_as_quotient(expr)?;
    let v = Var::from(var.as_str());
    let num_p = expr_to_poly(&num).ok()?;
    let nd = univariate_degree(&num_p, &v);
    if nd == 0 {
        return Some(Expr::int(0));
    }
    let den_x = sqrt_sum_leading_linear_coeff(&den, var)?;
    if nd > 1 {
        return Some(Expr::int(0));
    }
    if nd == 1 {
        if den_x.is_zero() {
            return Some(Expr::sym("+infinity"));
        }
        return Some(ratio_to_expr(&(coeff_at(&num_p, &v, nd) / den_x)));
    }
    None
}

fn sqrt_sum_leading_linear_coeff(den: &ExprArc, var: &Ident) -> Option<Ratio<BigInt>> {
    let Expr::Add(terms) = den.as_ref() else {
        let inner = sqrt_arg(den)?;
        return match is_monic_quadratic_leading(var, &inner) {
            Some(true) => Some(Ratio::one()),
            _ => None,
        };
    };
    let mut total = Ratio::from_integer(BigInt::from(0));
    for t in terms {
        if is_var(t, var) {
            total += Ratio::one();
            continue;
        }
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

fn is_indeterminate(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Frac(num, den) if is_zero(num) && is_zero(den)
    )
}

/// `(a/u)/(b/u) → a/b` after `x=1/u` (cancels common `u^-1` factor).
pub(crate) fn peel_shared_u_inv_in_frac(expr: &ExprArc, u: &Ident) -> ExprArc {
    let expr = simplify_reciprocal_sqrt(&rewrite_u_inv_sums(expr, u), u);
    let Some((num, den)) = try_as_rational(&expr, u) else {
        return expr;
    };
    let den = cancel_u_factors(&den, u);
    let num = peel_u_inv_factor(&num, u).unwrap_or(num);
    let den = peel_u_inv_factor(&den, u).unwrap_or(den);
    Arc::new(Expr::Frac(num, den))
}

fn cancel_u_factors(expr: &ExprArc, u: &Ident) -> ExprArc {
    let factors = flatten_mul(expr);
    let mut net = 0i64;
    let mut rest = Vec::new();
    for f in factors {
        if is_u_var(&f, u) {
            net += 1;
        } else if is_u_inv(&f, u) {
            net -= 1;
        } else {
            rest.push(Arc::clone(&f));
        }
    }
    for _ in 0..net {
        rest.push(var_to_expr(u));
    }
    for _ in 0..(-net).max(0) {
        rest.push(Expr::pow(var_to_expr(u), Expr::int(-1)));
    }
    if rest.is_empty() {
        Expr::int(1)
    } else if rest.len() == 1 {
        Arc::clone(&rest[0])
    } else {
        Expr::mul(rest)
    }
}

fn flatten_mul(expr: &ExprArc) -> Vec<ExprArc> {
    match expr.as_ref() {
        Expr::Mul(fs) => fs.iter().flat_map(flatten_mul).collect(),
        _ => vec![Arc::clone(expr)],
    }
}

fn normalize_limit_result(expr: &ExprArc, ctx: &Context) -> ExprArc {
    let mut out = collapse_unit_powers(expr);
    for _ in 0..4 {
        let evaluated = eval(out.as_ref(), ctx).unwrap_or_else(|_| Arc::clone(&out));
        let normalized = ratnormal(evaluated.as_ref(), ctx).unwrap_or(evaluated);
        let collapsed = collapse_unit_powers(&normalized);
        if collapsed == out {
            break;
        }
        out = collapsed;
    }
    out
}

fn expr_to_ratio(e: &ExprArc) -> Option<Ratio<BigInt>> {
    match e.as_ref() {
        Expr::Int(n) => Some(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Some(r.clone()),
        Expr::Frac(num, den) => {
            let nn = expr_to_ratio(num)?;
            let dd = expr_to_ratio(den)?;
            if dd.is_zero() {
                None
            } else {
                Some(nn / dd)
            }
        }
        _ => None,
    }
}

/// Net power of `u` in a multiplicative form (`u^{-1/2}` etc.).
fn u_exponent_bound(expr: &ExprArc, u: &Ident) -> Option<Ratio<BigInt>> {
    match expr.as_ref() {
        Expr::Pow(b, e) if is_u_var(b, u) => expr_to_ratio(e),
        Expr::Pow(b, e) if is_half_exponent(e) => {
            let inner = u_exponent_bound(b, u)?;
            Some(inner / Ratio::from_integer(2.into()))
        }
        Expr::Pow(b, e) if !depends_on_var(b, u) => {
            let _ = e;
            Some(Ratio::zero())
        }
        Expr::Mul(fs) => fs.iter().try_fold(Ratio::zero(), |acc, f| {
            Some(acc + u_exponent_bound(f, u)?)
        }),
        Expr::Frac(n, d) => Some(u_exponent_bound(n, u)? - u_exponent_bound(d, u)?),
        _ if !depends_on_var(expr, u) => Some(Ratio::zero()),
        _ => None,
    }
}

fn limit_from_fractional_u_valuation(expr: &ExprArc, u: &Ident, ctx: &Context) -> Option<ExprArc> {
    let rat = ratnormal(expr.as_ref(), ctx).unwrap_or_else(|_| Arc::clone(expr));
    let exp = u_exponent_bound(&rat, u)?;
    if exp > Ratio::zero() {
        return Some(Expr::int(0));
    }
    if exp < Ratio::zero() {
        return Some(Expr::sym("+infinity"));
    }
    None
}

fn is_inv_var(exp: &ExprArc, var: &Ident) -> bool {
    matches!(
        exp.as_ref(),
        Expr::Pow(b, e)
            if is_var(b, var) && matches!(e.as_ref(), Expr::Int(n) if n.is_negative())
    ) || matches!(
        exp.as_ref(),
        Expr::Frac(n, d) if matches!(n.as_ref(), Expr::Int(nn) if nn.is_one()) && is_var(d, var)
    )
}

/// `x/(x^ln(x)) → 0` at `+∞` (`x^ln(x)` grows faster than any polynomial).
fn limit_var_over_x_pow_ln(expr: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (num, den) = try_as_rational(expr, var)?;
    if !is_var(&num, var) {
        return None;
    }
    let Expr::Pow(base, exp) = den.as_ref() else {
        return None;
    };
    if !is_var(base, var) {
        return None;
    }
    match exp.as_ref() {
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && is_var(&args[0], var) => {
            Some(Expr::int(0))
        }
        _ => None,
    }
}

/// `poly/sqrt(poly)` at `+∞` when numerator degree ≥ 1 (e.g. `(1+x)/(sqrt(x+1)+1)`).
fn limit_poly_over_sqrt_at_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let (num, den) = try_as_quotient(expr)?;
    let num = ratnormal(num.as_ref(), ctx).unwrap_or(num);
    let v = Var::from(var.as_str());
    let num_p = expr_to_poly(&num).ok()?;
    if univariate_degree(&num_p, &v) < 1 {
        return None;
    }
    let inner = sqrt_term_in_expr(&den)?;
    let inner_p = expr_to_poly(&inner).ok()?;
    let nd = univariate_degree(&num_p, &v);
    let id = univariate_degree(&inner_p, &v);
    if nd <= id / 2 {
        return None;
    }
    if id >= 1 {
        return Some(Expr::sym("+infinity"));
    }
    None
}

fn sqrt_term_in_expr(e: &ExprArc) -> Option<ExprArc> {
    if let Some(inner) = sqrt_arg(e) {
        return Some(inner);
    }
    match e.as_ref() {
        Expr::Add(ts) | Expr::Mul(ts) => {
            for t in ts {
                if let Some(inner) = sqrt_term_in_expr(t) {
                    return Some(inner);
                }
            }
        }
        _ => {}
    }
    None
}

/// `(a^x + b^x + …)^(1/x) → max(a,b,…)` for positive constants.
fn limit_exp_sum_nth_root(expr: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let Expr::Pow(base, exp) = expr.as_ref() else {
        return None;
    };
    if !is_inv_var(exp, var) {
        return None;
    }
    let Expr::Add(terms) = base.as_ref() else {
        return None;
    };
    let mut max_base = None::<f64>;
    for t in terms {
        let Expr::Pow(b, e) = t.as_ref() else {
            return None;
        };
        if !is_var(e, var) {
            return None;
        }
        let c = try_const_f64(b)?;
        if c <= 0.0 {
            return None;
        }
        max_base = Some(max_base.map_or(c, |m| m.max(c)));
    }
    let m = max_base?;
    if (m - m.round()).abs() < 1e-12 {
        Some(Expr::int(m.round() as i64))
    } else {
        None
    }
}

/// `a + 1/u → (1+u)/u` so [`peel_u_inv_factor`] can cancel `u^-1` in numerators/denominators.
fn rewrite_u_inv_sums(expr: &ExprArc, u: &Ident) -> ExprArc {
    match expr.as_ref() {
        Expr::Add(ts) if ts.iter().any(|t| is_u_inv(t, u)) => {
            let mut sum = Vec::new();
            for t in ts {
                if is_u_inv(t, u) {
                    sum.push(Expr::int(1));
                } else {
                    sum.push(Expr::mul(vec![Arc::clone(t), var_to_expr(u)]));
                }
            }
            Expr::mul(vec![Expr::pow(var_to_expr(u), Expr::int(-1)), Expr::add(sum)])
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| rewrite_u_inv_sums(t, u)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| rewrite_u_inv_sums(t, u)).collect()),
        Expr::Pow(b, e) => Expr::pow(rewrite_u_inv_sums(b, u), rewrite_u_inv_sums(e, u)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            rewrite_u_inv_sums(n, u),
            rewrite_u_inv_sums(d, u),
        )),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(|t| rewrite_u_inv_sums(t, u)).collect()),
        _ => Arc::clone(expr),
    }
}

fn peel_u_inv_factor(e: &ExprArc, u: &Ident) -> Option<ExprArc> {
    if is_u_inv(e, u) {
        return Some(Expr::int(1));
    }
    match e.as_ref() {
        Expr::Mul(fs) => {
            let mut rest = Vec::new();
            let mut found = false;
            for f in fs {
                if is_u_inv(f, u) {
                    found = true;
                } else {
                    rest.push(Arc::clone(f));
                }
            }
            if !found {
                return None;
            }
            Some(if rest.is_empty() {
                Expr::int(1)
            } else if rest.len() == 1 {
                Arc::clone(&rest[0])
            } else {
                Expr::mul(rest)
            })
        }
        _ => None,
    }
}

fn is_u_inv(e: &ExprArc, u: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Pow(b, exp)
            if is_u_var(b, u) && matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
    ) || matches!(
        e.as_ref(),
        Expr::Frac(n, d) if matches!(n.as_ref(), Expr::Int(nn) if nn.is_one()) && is_u_var(d, u)
    )
}

fn is_u_var(e: &ExprArc, u: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == u)
}

fn is_half_exponent(exp: &ExprArc) -> bool {
    matches!(exp.as_ref(), Expr::Rat(r) if *r == Ratio::new(1.into(), 2.into()))
        || matches!(
            exp.as_ref(),
            Expr::Frac(n, d)
                if matches!(n.as_ref(), Expr::Int(nn) if nn.is_one())
                    && matches!(d.as_ref(), Expr::Int(dd) if dd == &BigInt::from(2))
        )
}

fn is_one_plus_u_inv_sq(b: &ExprArc, u: &Ident) -> bool {
    let Expr::Add(ts) = b.as_ref() else {
        return false;
    };
    ts.iter().any(|t| matches!(t.as_ref(), Expr::Int(n) if n.is_one()))
        && ts.iter().any(|t| u_negative_power_degree(t, u) == Some(-2))
}

fn u_negative_power_degree(e: &ExprArc, u: &Ident) -> Option<i64> {
    match e.as_ref() {
        Expr::Pow(b, exp) if is_u_var(b, u) => match exp.as_ref() {
            Expr::Int(n) => bigint_to_i64(n).ok(),
            _ => None,
        },
        Expr::Pow(inner, exp) => {
            let outer = match exp.as_ref() {
                Expr::Int(n) => bigint_to_i64(n).ok()?,
                _ => return None,
            };
            let inner_deg = u_negative_power_degree(inner, u)?;
            inner_deg.checked_mul(outer)
        }
        _ => None,
    }
}

/// `sqrt(1+1/u^2) → sqrt(1+u^2)/u` after reciprocal substitution.
fn simplify_reciprocal_sqrt(expr: &ExprArc, u: &Ident) -> ExprArc {
    match expr.as_ref() {
        Expr::Pow(b, exp) if is_half_exponent(exp) && is_one_plus_u_inv_sq(b, u) => Expr::mul(vec![
            Expr::pow(
                Expr::add(vec![Expr::int(1), Expr::pow(var_to_expr(u), Expr::int(2))]),
                Arc::clone(exp),
            ),
            Expr::pow(var_to_expr(u), Expr::int(-1)),
        ]),
        Expr::Mul(fs) if fs.len() == 2 => {
            for (i, j) in [(0, 1), (1, 0)] {
                if is_u_var(&fs[i], u) {
                    if let Expr::Pow(b, exp) = fs[j].as_ref() {
                        if is_half_exponent(exp) && is_one_plus_u_inv_sq(b, u) {
                            return Expr::pow(
                                Expr::add(vec![
                                    Expr::int(1),
                                    Expr::pow(var_to_expr(u), Expr::int(2)),
                                ]),
                                Arc::clone(exp),
                            );
                        }
                    }
                }
            }
            Expr::mul(fs.iter().map(|f| simplify_reciprocal_sqrt(f, u)).collect())
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| simplify_reciprocal_sqrt(t, u)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| simplify_reciprocal_sqrt(t, u)).collect()),
        Expr::Pow(b, e) => Expr::pow(simplify_reciprocal_sqrt(b, u), simplify_reciprocal_sqrt(e, u)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            simplify_reciprocal_sqrt(n, u),
            simplify_reciprocal_sqrt(d, u),
        )),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(|t| simplify_reciprocal_sqrt(t, u)).collect()),
        _ => Arc::clone(expr),
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
    fn asymptotic_series_exp_at_infinity() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let e = Expr::rat(1, 2);
        let r = asymptotic_series_at_infinity(&e, &var, 4, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/2");
    }

    #[test]
    fn asymptotic_maxima_sqrt_conjugate_shape() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let stmts =
            giac_parse::parse_program("x*(sqrt(1+x^2)-x);", &ctx).expect("parse");
        let giac_core::Stmt::ExprStmt(e) = stmts.first().expect("stmt") else {
            panic!("expected expr");
        };
        let r = limit_at_plus_infinity(e, &var, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/2");
    }

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
        use crate::limit_engine::ck_int_gruntz_fixture::ck_int_61;
        let ctx = xcas_default();
        let e = ck_int_61();
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

    /// Phase 0: reciprocal + series path only (skip MRV).
    mod fallback_only {
        use super::*;

        fn assert_fallback(input: &str, expected: &str) {
            let ctx = xcas_default();
            let var = Ident::new("x");
            let stmts = giac_parse::parse_program(&format!("{input};"), &ctx).expect("parse");
            let giac_core::Stmt::ExprStmt(e) = stmts.first().expect("stmt") else {
                panic!("expected expr");
            };
            let r = limit_at_plus_infinity_fallback(e, &var, &ctx).unwrap();
            assert_eq!(format_expr(r.as_ref()), expected, "input: {input}");
        }

        #[test]
        fn rational_one() {
            assert_fallback("(x+1)/(x-1)", "1");
        }

        #[test]
        fn sqrt_conjugate_x() {
            assert_fallback("x*(sqrt(1+x^2)-x)", "1/2");
        }

        #[test]
        fn ck_int_58() {
            assert_fallback("(x+1)/sqrt((x+1)/(x-1))", "+infinity");
        }

        #[test]
        fn atan_over_x_plus_one() {
            assert_fallback("x*atan(x)/(x+1)", "pi/2");
        }

        #[test]
        fn one_plus_one_over_x_sqrt() {
            assert_fallback("(1+1/x)*(sqrt(x+1)+1)", "+infinity");
        }
    }
}
