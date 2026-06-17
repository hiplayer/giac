//! Sparse series in one variable (GIAC-216a / `sparse_poly1` subset).
//!
//! Bounded: no `expand`, capped term count, no Taylor on heavy `exp` forms.

use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, expr_to_poly, Context, EvalError, Expr, ExprArc, FuncKind,
    Ident,
};
use giac_simplify::ratnormal;
use giac_poly::{coeff_at, univariate_degree, Poly, Var};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use super::bounds::{MAX_SERIES_DEPTH, MAX_SERIES_EXPANSION_ORDER, MAX_SERIES_ORDER, MAX_SERIES_TERMS};
use super::mrv_w::{
    decompose_ln_w_coeff, expr_contains_ln_w, is_expr_one, is_expr_zero as mrv_is_zero, is_mrv_w_var,
    is_neg_ln_w_expr, mrv_ln_w_expr,
};
use super::remove_lnexp::{expr_contains_exp_or_ln, remove_lnexp};
use crate::integrate::try_as_rational;
use crate::risch::depends_on_var;

/// `sum coeff * var^exponent` with integer exponents, sorted ascending.
#[derive(Clone, Debug, Default)]
pub(crate) struct SparseSeries {
    terms: Vec<(i32, ExprArc)>,
}

impl SparseSeries {
    pub(crate) fn constant(c: ExprArc) -> Self {
        if is_expr_zero(&c) {
            return Self::default();
        }
        Self {
            terms: vec![(0, c)],
        }
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub(crate) fn lead(&self) -> Option<(i32, ExprArc)> {
        self.terms.first().map(|(e, c)| (*e, Arc::clone(c)))
    }

    pub(crate) fn term_count(&self) -> usize {
        self.terms.len()
    }

    pub(crate) fn iter_terms(&self) -> impl Iterator<Item = (i32, &ExprArc)> + '_ {
        self.terms.iter().map(|(e, c)| (*e, c))
    }

    pub(crate) fn to_expr(&self, var: &Ident) -> ExprArc {
        if self.terms.is_empty() {
            return Expr::int(0);
        }
        let v = var_to_expr(var);
        let parts: Vec<ExprArc> = self
            .terms
            .iter()
            .map(|(exp, coeff)| {
                if *exp == 0 {
                    Arc::clone(coeff)
                } else {
                    Expr::mul(vec![
                        Arc::clone(coeff),
                        Expr::pow(Arc::clone(&v), Expr::int(i64::from(*exp))),
                    ])
                }
            })
            .collect();
        Expr::add(parts)
    }

    pub(crate) fn add(&self, other: &Self, ctx: &Context) -> Result<Self, EvalError> {
        let mut map: HashMap<i32, ExprArc> = HashMap::new();
        for (e, c) in &self.terms {
            merge_term(&mut map, *e, Arc::clone(c));
        }
        for (e, c) in &other.terms {
            merge_term(&mut map, *e, Arc::clone(c));
        }
        normalize_map(map, ctx)
    }

    pub(crate) fn mul_with_cap(
        &self,
        other: &Self,
        max_order: usize,
        order_cap: usize,
        ctx: &Context,
    ) -> Result<Self, EvalError> {
        if self.is_zero() || other.is_zero() {
            return Ok(Self::default());
        }
        if self.term_count() * other.term_count() > MAX_SERIES_TERMS * 4 {
            return Err(EvalError::NotImplemented("series"));
        }
        let cap = max_order.min(order_cap);
        let mut map: HashMap<i32, ExprArc> = HashMap::new();
        for (e1, c1) in &self.terms {
            for (e2, c2) in &other.terms {
                let exp = e1.saturating_add(*e2);
                if (exp as usize) > cap.saturating_add(2) {
                    continue;
                }
                merge_term(&mut map, exp, Expr::mul(vec![Arc::clone(c1), Arc::clone(c2)]));
                if map.len() > MAX_SERIES_TERMS * 2 {
                    return Err(EvalError::NotImplemented("series"));
                }
            }
        }
        let mut s = normalize_map(map, ctx)?;
        s.truncate_to_order_cap(cap);
        Ok(s)
    }

    pub(crate) fn mul(&self, other: &Self, max_order: usize, ctx: &Context) -> Result<Self, EvalError> {
        self.mul_with_cap(other, max_order, MAX_SERIES_ORDER, ctx)
    }

    fn truncate_to_order(&mut self, max_order: usize) {
        self.truncate_to_order_cap(max_order.min(MAX_SERIES_ORDER));
    }

    fn truncate_to_order_cap(&mut self, max_order: usize) {
        self.terms.sort_by(|a, b| a.0.cmp(&b.0));
        self.terms.retain(|(e, _)| {
            if *e < 0 {
                true
            } else {
                (*e as usize) <= max_order.saturating_add(2)
            }
        });
        self.terms.truncate(MAX_SERIES_TERMS);
    }

    pub(crate) fn map_coeffs<F>(&self, f: F) -> Self
    where
        F: Fn(&ExprArc) -> ExprArc,
    {
        Self {
            terms: self
                .terms
                .iter()
                .map(|(e, c)| (*e, f(c)))
                .collect(),
        }
    }

    pub(crate) fn from_rational_laurent(
        num: &ExprArc,
        den: &ExprArc,
        var: &Ident,
        order: usize,
    ) -> Result<Self, EvalError> {
        let v = Var::from(var.as_str());
        let num_p = expr_to_poly(num).map_err(|_| EvalError::NotImplemented("series"))?;
        let den_p = expr_to_poly(den).map_err(|_| EvalError::NotImplemented("series"))?;
        let vn = valuation_at_zero(&num_p, &v);
        let vd = valuation_at_zero(&den_p, &v);
        let base_exp = i32::try_from(vn.saturating_sub(vd)).unwrap_or(0);
        let mut q = num_p.clone();
        let mut d = den_p.clone();
        for _ in 0..vn {
            let lin = Poly::var(v.clone());
            if !q.div_rem(&lin).1.is_zero() {
                break;
            }
            q = q.div_rem(&lin).0;
        }
        for _ in 0..vd {
            let lin = Poly::var(v.clone());
            if !d.div_rem(&lin).1.is_zero() {
                break;
            }
            d = d.div_rem(&lin).0;
        }
        let (_, rem) = q.div_rem(&d);
        let mut terms = Vec::new();
        let lim = order.min(MAX_SERIES_ORDER);
        for k in 0..lim as u64 {
            let exp = base_exp.saturating_add(i32::try_from(k).unwrap_or(i32::MAX));
            let c = coeff_at(&rem, &v, k);
            if c.is_zero() {
                continue;
            }
            terms.push((exp, ratio_to_expr(&c)));
        }
        if terms.is_empty() && !rem.is_zero() {
            let c = coeff_at(&rem, &v, 0) / coeff_at(&d, &v, 0);
            terms.push((base_exp, ratio_to_expr(&c)));
        }
        Ok(Self { terms })
    }
}

/// Taylor series of `expr` in `var` about `center` (GIAC-216d).
pub(crate) fn series_at_center(
    expr: &ExprArc,
    var: &Ident,
    center: &ExprArc,
    order: usize,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if is_expr_zero(center) {
        let s = series_at_zero(expr, var, order, ctx)?;
        return series_sparse_to_expr(&s, var, center);
    }
    let t = Ident::new("_series_t");
    let x_sub = Expr::add(vec![var_to_expr(&t), Arc::clone(center)]);
    let mut subs = HashMap::new();
    subs.insert(var.clone(), x_sub);
    let shifted = eval_subst_map(expr, &subs)?;
    let s = series_at_zero(&shifted, &t, order, ctx)?;
    let delta = Expr::add(vec![
        var_to_expr(var),
        Expr::mul(vec![Expr::int(-1), Arc::clone(center)]),
    ]);
    let mut out = Vec::new();
    for (exp, coeff) in s.iter_terms() {
        let term = if exp == 0 {
            Arc::clone(coeff)
        } else {
            Expr::mul(vec![
                Arc::clone(coeff),
                Expr::pow(Arc::clone(&delta), Expr::int(i64::from(exp))),
            ])
        };
        out.push(term);
    }
    if out.is_empty() {
        return Ok(Expr::int(0));
    }
    eval(Expr::add(out).as_ref(), ctx)
}

fn series_sparse_to_expr(
    s: &SparseSeries,
    var: &Ident,
    center: &ExprArc,
) -> Result<ExprArc, EvalError> {
    let delta = if is_expr_zero(center) {
        var_to_expr(var)
    } else {
        Expr::add(vec![
            var_to_expr(var),
            Expr::mul(vec![Expr::int(-1), Arc::clone(center)]),
        ])
    };
    let mut out = Vec::new();
    for (exp, coeff) in s.iter_terms() {
        let term = if exp == 0 {
            Arc::clone(coeff)
        } else {
            Expr::mul(vec![
                Arc::clone(coeff),
                Expr::pow(Arc::clone(&delta), Expr::int(i64::from(exp))),
            ])
        };
        out.push(term);
    }
    if out.is_empty() {
        Ok(Expr::int(0))
    } else {
        Ok(Expr::add(out))
    }
}

/// Taylor/Laurent series of `expr` in `var` at `var = 0` (bounded).
pub(crate) fn series_at_zero(
    expr: &ExprArc,
    var: &Ident,
    order: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    series_at_zero_order(expr, var, order, MAX_SERIES_ORDER, ctx)
}

/// MRV/limit path: allow `order_cap` up to [`MAX_SERIES_EXPANSION_ORDER`].
pub(crate) fn series_at_zero_order(
    expr: &ExprArc,
    var: &Ident,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    let order_cap = order_cap.min(MAX_SERIES_EXPANSION_ORDER);
    let order = order.min(order_cap);
    series_at_zero_depth(expr, var, order, order_cap, 0, ctx)
}

/// giac `spdiv(sparse_poly1(1), p)` for MRV lead extraction.
pub(crate) fn series_spdiv_one(
    den: &SparseSeries,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    series_div(
        &SparseSeries::constant(Expr::int(1)),
        den,
        order,
        order_cap,
        ctx,
    )
}

fn series_at_zero_depth(
    expr: &ExprArc,
    var: &Ident,
    order: usize,
    order_cap: usize,
    depth: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    if depth > MAX_SERIES_DEPTH {
        return Err(EvalError::NotImplemented("series"));
    }
    if order == 0 {
        return Ok(SparseSeries::default());
    }
    if !depends_on_var(expr, var) {
        return Ok(SparseSeries::constant(eval(expr.as_ref(), ctx)?));
    }
    if !is_mrv_w_var(var) {
        if let Some((num, den)) = try_as_rational(expr, var) {
            return SparseSeries::from_rational_laurent(&num, &den, var, order);
        }
        if let Some(normalized) = ratnormal(expr.as_ref(), ctx).ok() {
            if let Some((num, den)) = try_as_rational(&normalized, var) {
                return SparseSeries::from_rational_laurent(&num, &den, var, order);
            }
        }
    }
    match expr.as_ref() {
        Expr::Symbol(id) if id == var => Ok(SparseSeries {
            terms: vec![(1, Expr::int(1))],
        }),
        Expr::Add(ts) => {
            if is_mrv_w_var(var) {
                let combined = remove_lnexp(
                    &Expr::add(ts.iter().map(|t| remove_lnexp(t, ctx)).collect()),
                    ctx,
                );
                let combined = ratnormal(combined.as_ref(), ctx).unwrap_or(combined);
                return match combined.as_ref() {
                    Expr::Add(ts2) => {
                        let mut acc = SparseSeries::default();
                        for t in ts2 {
                            acc = acc.add(
                                &series_at_zero_depth(t, var, order, order_cap, depth + 1, ctx)?,
                                ctx,
                            )?;
                        }
                        acc.truncate_to_order_cap(order.min(order_cap));
                        Ok(acc)
                    }
                    _ => series_at_zero_depth(&combined, var, order, order_cap, depth + 1, ctx),
                };
            }
            let mut acc = SparseSeries::default();
            for t in ts {
                acc = acc.add(
                    &series_at_zero_depth(t, var, order, order_cap, depth + 1, ctx)?,
                    ctx,
                )?;
            }
            acc.truncate_to_order_cap(order.min(order_cap));
            Ok(acc)
        }
        Expr::Mul(fs) => {
            let mut acc = SparseSeries::constant(Expr::int(1));
            for f in fs {
                acc = acc.mul_with_cap(
                    &series_at_zero_depth(f, var, order, order_cap, depth + 1, ctx)?,
                    order,
                    order_cap,
                    ctx,
                )?;
            }
            Ok(acc)
        }
        Expr::Frac(n, d) => {
            let num_s = series_at_zero_depth(n, var, order, order_cap, depth + 1, ctx)?;
            let den_s = series_at_zero_depth(d, var, order, order_cap, depth + 1, ctx)?;
            series_div(&num_s, &den_s, order, order_cap, ctx)
        }
        Expr::Pow(b, exp) => {
            if is_mrv_w_var(var) && is_mrv_symbolic_pow(b, exp) {
                return Ok(SparseSeries::constant(Expr::pow(
                    Arc::clone(b),
                    Arc::clone(exp),
                )));
            }
            if is_half_exponent(exp) {
                let base = series_at_zero_depth(b, var, order, order_cap, depth + 1, ctx)?;
                return series_sqrt(&base, order, order_cap, ctx);
            }
            if let Expr::Int(n) = exp.as_ref() {
                let base = series_at_zero_depth(b, var, order, order_cap, depth + 1, ctx)?;
                return series_pow_int(&base, bigint_to_i64(n)?, order, order_cap, ctx);
            }
            Err(EvalError::NotImplemented("series"))
        }
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            let arg = series_at_zero_depth(&args[0], var, order, order_cap, depth + 1, ctx)?;
            if is_mrv_w_var(var) {
                if arg_has_symbolic_ln_w(&arg) {
                    return Ok(SparseSeries::constant(Expr::func(
                        FuncKind::Exp,
                        vec![arg.to_expr(var)],
                    )));
                }
                return series_exp_mrv(&arg, order, order_cap, ctx).or_else(|_| {
                    Ok(SparseSeries::constant(Expr::func(
                        FuncKind::Exp,
                        vec![arg.to_expr(var)],
                    )))
                });
            }
            if arg.lead().is_some_and(|(e, _)| e < 0) {
                return Err(EvalError::NotImplemented("series"));
            }
            series_exp(&arg, order, order_cap, ctx)
        }
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => {
            if is_mrv_w_var(var)
                && matches!(args[0].as_ref(), Expr::Symbol(id) if is_mrv_w_var(id))
            {
                return Ok(SparseSeries::constant(mrv_ln_w_expr()));
            }
            if is_series_var(&args[0], var) && !is_mrv_w_var(var) {
                return Ok(SparseSeries::constant(Expr::func(
                    FuncKind::Ln,
                    vec![var_to_expr(var)],
                )));
            }
            let arg = series_at_zero_depth(&args[0], var, order, order_cap, depth + 1, ctx)?;
            if is_mrv_w_var(var) {
                series_ln_mrv(&arg, order, ctx)
            } else {
                Err(EvalError::NotImplemented("series"))
            }
        }
        Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_series_var(&args[0], var) => {
            series_sin(order)
        }
        Expr::Func(FuncKind::Cos, args) if args.len() == 1 && is_series_var(&args[0], var) => {
            series_cos(order)
        }
        Expr::Func(FuncKind::Sqrt, args) if args.len() == 1 => {
            let arg = series_at_zero_depth(&args[0], var, order, order_cap, depth + 1, ctx)?;
            series_sqrt(&arg, order, order_cap, ctx)
        }
        Expr::Func(FuncKind::Atan, args) if args.len() == 1 => {
            if is_inv_series_var(&args[0], var) {
                return series_atan_of_inv(order);
            }
            if is_series_var(&args[0], var) {
                return series_atan(order);
            }
            Err(EvalError::NotImplemented("series"))
        }
        _ => Err(EvalError::NotImplemented("series")),
    }
}

fn series_sin(order: usize) -> Result<SparseSeries, EvalError> {
    let mut terms = Vec::new();
    let lim = order.min(MAX_SERIES_ORDER);
    let mut k = 1usize;
    let mut sign = 1i64;
    let mut fact = 1i64;
    while k < lim {
        terms.push((k as i32, Expr::rat(sign, fact)));
        sign = -sign;
        k += 2;
        if k < lim {
            fact = fact
                .checked_mul((k - 1) as i64)
                .and_then(|f| f.checked_mul(k as i64))
                .ok_or(EvalError::NotImplemented("series"))?;
        }
    }
    Ok(SparseSeries { terms })
}

fn series_cos(order: usize) -> Result<SparseSeries, EvalError> {
    let mut terms = vec![(0, Expr::int(1))];
    let lim = order.min(MAX_SERIES_ORDER);
    let mut k = 2usize;
    let mut sign = -1i64;
    let mut fact = 1i64;
    while k < lim {
        fact = fact
            .checked_mul((k - 1) as i64)
            .and_then(|f| f.checked_mul(k as i64))
            .ok_or(EvalError::NotImplemented("series"))?;
        terms.push((k as i32, Expr::rat(sign, fact)));
        sign = -sign;
        k += 2;
    }
    Ok(SparseSeries { terms })
}

fn series_sqrt(
    arg: &SparseSeries,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    let one = SparseSeries::constant(Expr::int(1));
    let neg_one = SparseSeries::constant(Expr::int(-1));
    let delta = arg.add(&neg_one, ctx)?;
    let mut acc = one;
    let mut term = delta.clone();
    let mut binom = Ratio::new(BigInt::from(1), BigInt::from(2));
    let lim = order.min(order_cap).min(MAX_SERIES_EXPANSION_ORDER);
    for k in 1..lim {
        let scale = ratio_to_expr(&binom);
        let scaled = term.map_coeffs(|c| Expr::mul(vec![Arc::clone(c), Arc::clone(&scale)]));
        acc = acc.add(&scaled, ctx)?;
        term = term.mul_with_cap(&delta, order, order_cap, ctx)?;
        binom = binom * Ratio::new(BigInt::from(3 - 2 * (k as i64)), BigInt::from(2 * (k as i64)));
    }
    Ok(acc)
}

/// `atan(u) = u - u^3/3 + u^5/5 - ...`
fn series_atan(order: usize) -> Result<SparseSeries, EvalError> {
    let mut terms = Vec::new();
    let lim = order.min(MAX_SERIES_ORDER);
    let mut k = 1usize;
    let mut sign = 1i64;
    while k < lim {
        terms.push((k as i32, Expr::rat(sign, k as i64)));
        sign = -sign;
        k += 2;
    }
    Ok(SparseSeries { terms })
}

/// `atan(1/u) = pi/2 - u + u^3/3 - u^5/5 + ...` for `u → 0+`.
fn series_atan_of_inv(order: usize) -> Result<SparseSeries, EvalError> {
    let mut terms = vec![(0, Arc::new(Expr::Frac(Expr::sym("pi"), Expr::int(2))))];
    let lim = order.min(MAX_SERIES_ORDER);
    let mut k = 1usize;
    let mut sign = -1i64;
    while k < lim {
        terms.push((k as i32, Expr::rat(sign, k as i64)));
        sign = -sign;
        k += 2;
    }
    Ok(SparseSeries { terms })
}

fn is_inv_series_var(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Pow(b, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) && is_series_var(b, var) =>
        {
            true
        }
        Expr::Frac(n, d) if matches!(n.as_ref(), Expr::Int(n) if n.is_one()) && is_series_var(d, var) => {
            true
        }
        _ => false,
    }
}

/// `exp(arg)` when `arg` may contain `k*ln(w)` in the constant term (MRV series).
fn series_exp_mrv(
    arg: &SparseSeries,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    let min_e = arg.lead().map(|(e, _)| e).unwrap_or(0);
    if min_e < 0 {
        let shifted = shift_series_exponents(arg, -min_e);
        let exp_s = series_exp_mrv_positive(&shifted, order, order_cap, ctx)?;
        let w_part = SparseSeries {
            terms: vec![(min_e, Expr::int(1))],
        };
        return w_part.mul_with_cap(&exp_s, order, order_cap, ctx);
    }
    series_exp_mrv_positive(arg, order, order_cap, ctx)
}

fn shift_series_exponents(s: &SparseSeries, delta: i32) -> SparseSeries {
    SparseSeries {
        terms: s
            .terms
            .iter()
            .map(|(e, c)| (e.saturating_add(delta), Arc::clone(c)))
            .collect(),
    }
}

fn series_exp_mrv_positive(
    arg: &SparseSeries,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    let mut const_term = Expr::int(0);
    let mut hi_terms = Vec::new();
    for (e, c) in &arg.terms {
        if *e == 0 {
            const_term = Expr::add(vec![const_term, Arc::clone(c)]);
        } else if *e < 0 {
            return Err(EvalError::NotImplemented("series"));
        } else {
            hi_terms.push((*e, Arc::clone(c)));
        }
    }
    let (k, a_rest) = decompose_ln_w_coeff(&const_term);
    let w_part = if k == 0 {
        SparseSeries::constant(Expr::int(1))
    } else {
        SparseSeries {
            terms: vec![(k, Expr::int(1))],
        }
    };
    let exp_a = if mrv_is_zero(&a_rest) || is_expr_one(&a_rest) {
        SparseSeries::constant(Expr::int(1))
    } else {
        SparseSeries::constant(Expr::func(FuncKind::Exp, vec![a_rest]))
    };
    let hi = SparseSeries { terms: hi_terms };
    let taylor = if hi.is_zero() {
        SparseSeries::constant(Expr::int(1))
    } else {
        series_exp(&hi, order, order_cap, ctx)?
    };
    Ok(w_part
        .mul_with_cap(&exp_a, order, order_cap, ctx)?
        .mul_with_cap(&taylor, order, order_cap, ctx)?)
}

fn series_ln_mrv(arg: &SparseSeries, order: usize, ctx: &Context) -> Result<SparseSeries, EvalError> {
    let (min_e, lead_c) = arg.lead().ok_or(EvalError::NotImplemented("series"))?;
    if min_e > 0 {
        let mut shifted = Vec::new();
        for (e, c) in &arg.terms {
            shifted.push((e - min_e, Arc::clone(c)));
        }
        let shifted_s = SparseSeries { terms: shifted };
        let (k2, c_rest) = shifted_s
            .lead()
            .map(|(_, c)| decompose_ln_w_coeff(&c))
            .unwrap_or((0, Expr::int(1)));
        let total_k = i32::try_from(min_e).unwrap_or(i32::MAX).saturating_add(k2);
        let mut parts = Vec::new();
        if total_k != 0 {
            parts.push(Expr::mul(vec![Expr::int(i64::from(total_k)), mrv_ln_w_expr()]));
        }
        if !is_expr_one(&c_rest) && !mrv_is_zero(&c_rest) {
            parts.push(Expr::func(FuncKind::Ln, vec![c_rest]));
        }
        if parts.is_empty() {
            return Ok(SparseSeries::constant(Expr::int(0)));
        }
        return Ok(SparseSeries::constant(Expr::add(parts)));
    }
    if min_e == 0 {
        let (k, rest) = decompose_ln_w_coeff(&lead_c);
        let mut parts = Vec::new();
        if k != 0 {
            parts.push(Expr::mul(vec![Expr::int(i64::from(k)), mrv_ln_w_expr()]));
        }
        if !is_expr_one(&rest) && !mrv_is_zero(&rest) {
            parts.push(Expr::func(FuncKind::Ln, vec![rest]));
        }
        if parts.is_empty() {
            return Ok(SparseSeries::constant(Expr::int(0)));
        }
        return Ok(SparseSeries::constant(Expr::add(parts)));
    }
    let _ = order;
    let _ = ctx;
    Err(EvalError::NotImplemented("series"))
}

fn series_exp(
    arg: &SparseSeries,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    let mut acc = SparseSeries::constant(Expr::int(1));
    let mut term = SparseSeries::constant(Expr::int(1));
    for k in 1..order.min(8) {
        term = term.mul_with_cap(arg, order, order_cap, ctx)?;
        let coeff = Expr::rat(1, factorial(k)?);
        let scaled = SparseSeries {
            terms: term
                .terms
                .iter()
                .map(|(e, c)| (*e, Expr::mul(vec![Arc::clone(c), coeff.clone()])))
                .collect(),
        };
        acc = acc.add(&scaled, ctx)?;
        if acc.term_count() > MAX_SERIES_TERMS {
            break;
        }
    }
    acc.truncate_to_order_cap(order.min(order_cap));
    Ok(acc)
}

fn series_pow_int(
    base: &SparseSeries,
    n: i64,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    if n == 0 {
        return Ok(SparseSeries::constant(Expr::int(1)));
    }
    if n < 0 {
        let pos = series_pow_int(base, -n, order, order_cap, ctx)?;
        return series_div(
            &SparseSeries::constant(Expr::int(1)),
            &pos,
            order,
            order_cap,
            ctx,
        );
    }
    if n > 8 {
        return Err(EvalError::NotImplemented("series"));
    }
    let mut acc = SparseSeries::constant(Expr::int(1));
    let mut b = base.clone();
    let mut exp = n;
    while exp > 0 {
        if exp % 2 == 1 {
            acc = acc.mul_with_cap(&b, order, order_cap, ctx)?;
        }
        b = b.mul_with_cap(&b, order, order_cap, ctx)?;
        exp /= 2;
    }
    Ok(acc)
}

fn series_div(
    num: &SparseSeries,
    den: &SparseSeries,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    num.mul_with_cap(&series_inv(den, order, order_cap, ctx)?, order, order_cap, ctx)
}

fn series_inv(
    den: &SparseSeries,
    order: usize,
    order_cap: usize,
    ctx: &Context,
) -> Result<SparseSeries, EvalError> {
    let (_, lead_coeff) = den.lead().ok_or(EvalError::NotImplemented("series"))?;
    let inv_lead = eval(
        Expr::pow(Arc::clone(&lead_coeff), Expr::int(-1)).as_ref(),
        ctx,
    )
    .unwrap_or_else(|_| Expr::pow(Arc::clone(&lead_coeff), Expr::int(-1)));
    let mut acc = SparseSeries::constant(inv_lead.clone());
    let mut rest = den.clone();
    if let Some((_, c)) = rest.terms.first_mut() {
        *c = Expr::mul(vec![Arc::clone(c), inv_lead]);
    }
    if rest.terms.first().is_some_and(|(e, _)| *e == 0) {
        rest.terms.remove(0);
    }
    for _ in 0..order.min(8) {
        let prod = acc.mul_with_cap(&rest, order, order_cap, ctx)?;
        let neg = SparseSeries {
            terms: prod
                .terms
                .iter()
                .map(|(e, c)| (*e, Expr::mul(vec![Expr::int(-1), Arc::clone(c)])))
                .collect(),
        };
        acc = acc.add(&neg, ctx)?;
        if acc.term_count() > MAX_SERIES_TERMS {
            break;
        }
    }
    Ok(acc)
}

fn merge_term(map: &mut HashMap<i32, ExprArc>, exp: i32, coeff: ExprArc) {
    map.entry(exp)
        .and_modify(|c| *c = Expr::add(vec![Arc::clone(c), Arc::clone(&coeff)]))
        .or_insert(coeff);
}

fn normalize_map(map: HashMap<i32, ExprArc>, ctx: &Context) -> Result<SparseSeries, EvalError> {
    let mut terms: Vec<(i32, ExprArc)> = map
        .into_iter()
        .filter_map(|(e, c)| {
            let ev = if expr_contains_exp_or_ln(&c) {
                let r = remove_lnexp(&c, ctx);
                ratnormal(r.as_ref(), ctx)
                    .ok()
                    .or_else(|| eval(r.as_ref(), ctx).ok())
            } else {
                ratnormal(c.as_ref(), ctx)
                    .ok()
                    .or_else(|| eval(c.as_ref(), ctx).ok())
            }?;
            if is_expr_zero(&ev) {
                None
            } else {
                Some((e, ev))
            }
        })
        .collect();
    if terms.len() > MAX_SERIES_TERMS {
        return Err(EvalError::NotImplemented("series"));
    }
    terms.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(SparseSeries { terms })
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

fn bigint_to_i64(n: &BigInt) -> Result<i64, EvalError> {
    giac_core::bigint_to_i64(n).map_err(|_| EvalError::TypeError("int"))
}

fn factorial(n: usize) -> Result<i64, EvalError> {
    let mut acc = 1_i64;
    for i in 2..=n {
        acc = acc
            .checked_mul(i as i64)
            .ok_or(EvalError::NotImplemented("series"))?;
    }
    Ok(acc)
}

fn var_to_expr(var: &Ident) -> ExprArc {
    Expr::sym(var.as_str())
}

fn is_expr_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

fn is_series_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn arg_has_symbolic_ln_w(arg: &SparseSeries) -> bool {
    arg.iter_terms()
        .any(|(_, c)| expr_contains_ln_w(c) || is_neg_ln_w_expr(c))
}

/// `(-ln(w))^k` and similar non-Taylor powers stay symbolic in MRV series coeffs.
fn is_half_exponent(exp: &ExprArc) -> bool {
    use num_bigint::BigInt;
    use num_rational::Ratio;
    matches!(exp.as_ref(), Expr::Rat(r) if *r == Ratio::new(1.into(), 2.into()))
        || matches!(
            exp.as_ref(),
            Expr::Frac(n, d)
                if matches!(n.as_ref(), Expr::Int(nn) if nn.is_one())
                    && matches!(d.as_ref(), Expr::Int(dd) if dd == &BigInt::from(2))
        )
}

fn is_mrv_symbolic_pow(base: &ExprArc, exp: &ExprArc) -> bool {
    if !matches!(exp.as_ref(), Expr::Int(_)) {
        return false;
    }
    is_neg_ln_w_expr(base)
        || matches!(
            base.as_ref(),
            Expr::Func(FuncKind::Ln, args)
                if args.len() == 1 && matches!(args[0].as_ref(), Expr::Symbol(id) if is_mrv_w_var(id))
        )
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn sparse_series_rational_at_zero() {
        let ctx = xcas_default();
        let var = Ident::new("w");
        let e = Expr::add(vec![var_to_expr(&var), Expr::int(1)]);
        let s = series_at_zero(&e, &var, 4, &ctx).unwrap();
        let (exp, coeff) = s.lead().unwrap();
        assert_eq!(exp, 0);
        assert_eq!(format_expr(coeff.as_ref()), "1");
    }

    #[test]
    fn sparse_series_sin_at_zero() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let e = Expr::func(FuncKind::Sin, vec![var_to_expr(&var)]);
        let s = series_at_zero(&e, &var, 5, &ctx).unwrap();
        let out = s.to_expr(&var);
        let text = format_expr(out.as_ref());
        assert!(text.contains("x") && !text.contains("x^5"), "got {text}");
    }

    #[test]
    fn sparse_series_at_center_shift() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let e = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2));
        let r = series_at_center(&e, &var, &Expr::int(0), 3, &ctx).unwrap();
        let text = format_expr(r.as_ref());
        assert!(text.contains("x^2") || text.contains("2*x"), "got {text}");
    }

    #[test]
    fn sparse_series_atan_inv_over_one_plus_u() {
        let ctx = xcas_default();
        let var = Ident::new("u");
        let atan = Expr::func(FuncKind::Atan, vec![Expr::pow(var_to_expr(&var), Expr::int(-1))]);
        let s = series_at_zero(&atan, &var, 12, &ctx).unwrap();
        let (exp, coeff) = s.lead().unwrap();
        assert_eq!(exp, 0);
        assert_eq!(format_expr(coeff.as_ref()), "pi/2");
    }
}
