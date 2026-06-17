//! `mrv_lead_term` and limit from MRV series (GIAC-216c / GIAC-216e).

use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, expr_to_poly, poly_to_expr, Context, EvalError, Expr, ExprArc,
    FuncKind, Ident,
};
use giac_simplify::ratnormal;
use num_traits::{Signed, Zero};

use super::bounds::{
    mrv_limit_eligible, mrv_rewrite_bounded, MAX_SERIES_EXPANSION_ORDER, MAX_SERIES_ORDER,
};
use super::mrv::{choose_mrv_w, linear_coeff_in_var, mrv_at_plus_infinity, is_negative_const_expr};
use super::mrv_w::{
    decompose_ln_w_coeff, expr_contains_ln_w, is_expr_zero as mrv_is_zero, is_mrv_w_var, MRV_W,
};
use super::mrv_series_lead::normalize_expr_quotients;
use super::preprocess::limit_preprocess_plus_infinity;
use super::remove_lnexp::{divide_lead_coeffs, remove_lnexp};
use super::sparse_series::{series_at_zero_order, series_spdiv_one, SparseSeries};

#[derive(Clone, Debug)]
pub(crate) struct MrvLeadTerm {
    pub coeff: ExprArc,
    pub exponent: i32,
}

/// giac `mrv_lead_term` at `+infinity` (limit mode, bounded).
pub(crate) fn mrv_lead_term_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    if !mrv_limit_eligible(expr) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let pre = limit_preprocess_plus_infinity(expr, var, ctx)?;
    let pre = upscale_while_var_in_mrv(&pre, var, ctx);
    let pre = normalize_expr_quotients(&pre);
    if !mrv_limit_eligible(&pre) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let set = mrv_at_plus_infinity(&pre, var, ctx);
    if set.is_empty() {
        if crate::risch::depends_on_var(&pre, var) {
            return Err(EvalError::NotImplemented("limit"));
        }
        return Ok(MrvLeadTerm {
            coeff: ratnormal(pre.as_ref(), ctx).unwrap_or(pre),
            exponent: 0,
        });
    }
    let (omega, _) = choose_mrv_w(&set, var, ctx).ok_or(EvalError::NotImplemented("limit"))?;
    let w = Ident::new(MRV_W);
    let swapped = rewrite_in_mrv_w(&pre, var, &omega, &w);
    if !mrv_rewrite_bounded(&swapped) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let swapped = ratnormal(swapped.as_ref(), ctx).unwrap_or(swapped);
    let g = g_from_omega(&omega, var, &w);
    let dont_invert = omega_tends_to_zero_at_plus_inf(&omega, var, ctx);
    let omega_is_exp = matches!(omega.as_ref(), Expr::Func(FuncKind::Exp, _));
    mrv_series_lead_loop(
        &swapped,
        &w,
        &g,
        dont_invert,
        omega_is_exp,
        MAX_SERIES_ORDER,
        ctx,
    )
}

pub(crate) fn limit_from_mrv_lead_term(
    lead: &MrvLeadTerm,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if lead.exponent > 0 {
        return Ok(Expr::int(0));
    }
    if lead.exponent < 0 {
        return Ok(sign_infinity(&eval(lead.coeff.as_ref(), ctx)?));
    }
    let mut coeff = Arc::clone(&lead.coeff);
    if contains_ln_w(&coeff) {
        coeff = rewrite_ln_w(&coeff, var);
        if crate::risch::depends_on_var(&coeff, var) {
            return Err(EvalError::NotImplemented("limit"));
        }
    }
    let coeff = eval(coeff.as_ref(), ctx)?;
    if contains_w(&coeff) {
        return Err(EvalError::NotImplemented("limit"));
    }
    Ok(coeff)
}

/// upstream `unidirectional_limit` at `+infinity`: `mrv_lead_term` then recurse on coeff.
pub(crate) fn limit_unidirectional_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let mut e_copy = Arc::clone(expr);
    for _ in 0..8 {
        let lead = mrv_lead_term_plus_infinity(&e_copy, var, ctx)?;
        if lead.exponent > 0 {
            return Ok(Expr::int(0));
        }
        if lead.exponent < 0 {
            return Ok(sign_infinity(&eval(lead.coeff.as_ref(), ctx)?));
        }
        let coeff = limit_from_mrv_lead_term(&lead, var, ctx)?;
        if !crate::risch::depends_on_var(&coeff, var) {
            return Ok(coeff);
        }
        e_copy = coeff;
    }
    Err(EvalError::NotImplemented("limit"))
}

/// Replace `exp(±x)` with `w`/`w^-1` and `x` with `-ln(w)` when `omega = exp(-x)`.
fn subst_symbol(expr: &ExprArc, from: &Ident, to: &ExprArc) -> ExprArc {
    match expr.as_ref() {
        Expr::Symbol(id) if id == from => Arc::clone(to),
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| subst_symbol(t, from, to)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| subst_symbol(t, from, to)).collect()),
        Expr::Pow(b, e) => Expr::pow(subst_symbol(b, from, to), subst_symbol(e, from, to)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(subst_symbol(n, from, to), subst_symbol(d, from, to))),
        Expr::Func(kind, args) => Expr::func(
            *kind,
            args.iter().map(|a| subst_symbol(a, from, to)).collect(),
        ),
        _ => Arc::clone(expr),
    }
}

/// giac `upscale`: `ln(x)→x`, `x→exp(x)` while `x` remains in the MRV set.
fn upscale_while_var_in_mrv(expr: &ExprArc, var: &Ident, ctx: &Context) -> ExprArc {
    let ln_x = Expr::func(FuncKind::Ln, vec![Expr::sym(var.as_str())]);
    let exp_x = Expr::func(FuncKind::Exp, vec![Expr::sym(var.as_str())]);
    let mut out = Arc::clone(expr);
    for _ in 0..4 {
        let set = mrv_at_plus_infinity(&out, var, ctx);
        if !set.faster.iter().any(|f| is_var(f, var)) {
            break;
        }
        out = subst_symbol(&out, var, &exp_x);
        out = subst_expr_once(&out, &ln_x, &Expr::sym(var.as_str()));
    }
    out
}

fn rewrite_in_mrv_w(expr: &ExprArc, var: &Ident, omega: &ExprArc, w: &Ident) -> ExprArc {
    let w_expr = Expr::sym(w.as_str());
    if expr_eq(expr, omega) {
        return w_expr;
    }
    if is_var(expr, var) && is_exp_neg_var(omega, var) {
        return Expr::mul(vec![Expr::int(-1), Expr::func(FuncKind::Ln, vec![w_expr])]);
    }
    if is_exp_pos_var(expr, var) && is_exp_neg_var(omega, var) {
        return Expr::pow(w_expr, Expr::int(-1));
    }
    match expr.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| rewrite_in_mrv_w(t, var, omega, w)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| rewrite_in_mrv_w(t, var, omega, w)).collect()),
        Expr::Pow(b, e) => Expr::pow(rewrite_in_mrv_w(b, var, omega, w), rewrite_in_mrv_w(e, var, omega, w)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            rewrite_in_mrv_w(n, var, omega, w),
            rewrite_in_mrv_w(d, var, omega, w),
        )),
        Expr::Func(kind, args) => Expr::func(
            *kind,
            args.iter()
                .map(|a| rewrite_in_mrv_w(a, var, omega, w))
                .collect(),
        ),
        _ => Arc::clone(expr),
    }
}

fn expr_eq(a: &ExprArc, b: &ExprArc) -> bool {
    a == b
}

fn is_exp_neg_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 && is_neg_var(&args[0], var)
    )
}

fn is_exp_pos_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 && matches!(args[0].as_ref(), Expr::Symbol(id) if id == var)
    )
}

fn is_neg_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Mul(fs) if fs.len() == 2
            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Symbol(id) if id == var))
    )
}

fn rewrite_ln_w(expr: &ExprArc, var: &Ident) -> ExprArc {
    if expr_contains_ln_w(expr) {
        let (k, rest) = decompose_ln_w_coeff(expr);
        if k != 0 {
            let mut parts = vec![Expr::mul(vec![
                Expr::int(i64::from(k)),
                Expr::mul(vec![Expr::int(-1), Expr::sym(var.as_str())]),
            ])];
            if !super::mrv_w::is_expr_one(&rest) && !super::mrv_w::is_expr_zero(&rest) {
                parts.push(rewrite_ln_w(&rest, var));
            }
            return Expr::add(parts);
        }
    }
    match expr.as_ref() {
        Expr::Func(FuncKind::Ln, args)
            if args.len() == 1
                && matches!(args[0].as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)) =>
        {
            Expr::mul(vec![Expr::int(-1), Expr::sym(var.as_str())])
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| rewrite_ln_w(t, var)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| rewrite_ln_w(t, var)).collect()),
        Expr::Pow(b, e) => Expr::pow(rewrite_ln_w(b, var), rewrite_ln_w(e, var)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(rewrite_ln_w(n, var), rewrite_ln_w(d, var))),
        Expr::Func(kind, args) => Expr::func(
            *kind,
            args.iter().map(|a| rewrite_ln_w(a, var)).collect(),
        ),
        _ => Arc::clone(expr),
    }
}

fn g_from_omega(omega: &ExprArc, var: &Ident, w: &Ident) -> ExprArc {
    match omega.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            rewrite_in_mrv_w(&args[0], var, omega, w)
        }
        _ => Expr::int(0),
    }
}

fn omega_tends_to_zero_at_plus_inf(omega: &ExprArc, var: &Ident, ctx: &Context) -> bool {
    match omega.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            if is_neg_var_linear(&args[0], var) {
                return true;
            }
            linear_coeff_in_var(&args[0], var)
                .is_some_and(|c| is_negative_const_expr(&c, ctx))
        }
        _ => false,
    }
}

fn is_neg_var_linear(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Mul(fs) if fs.len() == 2
            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
            && fs.iter().any(|f| is_var(f, var))
    )
}

fn subst_w_inv(expr: &ExprArc, w: &Ident) -> Result<ExprArc, EvalError> {
    let mut m = HashMap::new();
    m.insert(w.clone(), Expr::pow(Expr::sym(w.as_str()), Expr::int(-1)));
    eval_subst_map(expr, &m)
}

fn is_ln_w(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Ln, args)
            if args.len() == 1 && matches!(args[0].as_ref(), Expr::Symbol(id) if is_mrv_w_var(id))
    )
}

fn subst_ln_w_expr(expr: &ExprArc, replacement: &ExprArc) -> ExprArc {
    if is_ln_w(expr) {
        return Arc::clone(replacement);
    }
    match expr.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| subst_ln_w_expr(t, replacement)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| subst_ln_w_expr(t, replacement)).collect()),
        Expr::Pow(b, e) => Expr::pow(subst_ln_w_expr(b, replacement), subst_ln_w_expr(e, replacement)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            subst_ln_w_expr(n, replacement),
            subst_ln_w_expr(d, replacement),
        )),
        Expr::Func(kind, args) => Expr::func(
            *kind,
            args.iter()
                .map(|a| subst_ln_w_expr(a, replacement))
                .collect(),
        ),
        _ => Arc::clone(expr),
    }
}

fn collect_ln_exprs(expr: &ExprArc, out: &mut Vec<ExprArc>) {
    if is_ln_w(expr) {
        out.push(Arc::clone(expr));
        return;
    }
    match expr.as_ref() {
        Expr::Add(ts) | Expr::Mul(ts) => {
            for t in ts {
                collect_ln_exprs(t, out);
            }
        }
        Expr::Pow(b, e) => {
            collect_ln_exprs(b, out);
            collect_ln_exprs(e, out);
        }
        Expr::Frac(n, d) => {
            collect_ln_exprs(n, out);
            collect_ln_exprs(d, out);
        }
        Expr::Func(_, args) => {
            for a in args {
                collect_ln_exprs(a, out);
            }
        }
        _ => {}
    }
}

/// upstream `ln(exp(g)^k*...) -> k*g + ln(...)` when MRV element is `exp`.
fn rewrite_ln_exp_in_f(
    f: &ExprArc,
    w: &Ident,
    g: &ExprArc,
    dont_invert: bool,
    begin_ordre: usize,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let order_cap = MAX_SERIES_EXPANSION_ORDER;
    let mut ln_nodes = Vec::new();
    collect_ln_exprs(f, &mut ln_nodes);
    if ln_nodes.is_empty() {
        return Ok(Arc::clone(f));
    }
    let g_use = if dont_invert {
        Arc::clone(g)
    } else {
        Expr::mul(vec![Expr::int(-1), Arc::clone(g)])
    };
    let mut out = Arc::clone(f);
    for ln_e in ln_nodes {
        let Expr::Func(FuncKind::Ln, args) = ln_e.as_ref() else {
            continue;
        };
        let argln = &args[0];
        let Ok(s) = series_at_zero_order(argln, w, begin_ordre, order_cap, ctx) else {
            continue;
        };
        let Some((lead_exp, lead_c)) = s.lead() else {
            continue;
        };
        if is_series_coeff_undef(&lead_c) {
            continue;
        }
        let mut arg = Arc::clone(argln);
        if lead_exp != 0 {
            arg = Expr::mul(vec![
                Arc::clone(argln),
                Expr::pow(Expr::sym(w.as_str()), Expr::int(-i64::from(lead_exp))),
            ]);
        }
        let new_ln = Expr::add(vec![
            Expr::mul(vec![Expr::int(i64::from(lead_exp)), Arc::clone(&g_use)]),
            Expr::func(FuncKind::Ln, vec![arg]),
        ]);
        out = subst_expr_once(&out, &ln_e, &new_ln);
    }
    Ok(out)
}

fn subst_expr_once(expr: &ExprArc, from: &ExprArc, to: &ExprArc) -> ExprArc {
    if expr == from {
        return Arc::clone(to);
    }
    match expr.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| subst_expr_once(t, from, to)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| subst_expr_once(t, from, to)).collect()),
        Expr::Pow(b, e) => Expr::pow(subst_expr_once(b, from, to), subst_expr_once(e, from, to)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            subst_expr_once(n, from, to),
            subst_expr_once(d, from, to),
        )),
        Expr::Func(kind, args) => Expr::func(
            *kind,
            args.iter().map(|a| subst_expr_once(a, from, to)).collect(),
        ),
        _ => Arc::clone(expr),
    }
}

fn peel_neg_ln_w_inv(expr: &ExprArc, w: &Ident) -> Option<(ExprArc, ExprArc)> {
    if let Expr::Frac(n, d) = expr.as_ref() {
        if super::mrv_w::is_neg_ln_w_expr(d) {
            let ln_inv = Expr::pow(Arc::clone(d), Expr::int(-1));
            return Some((Arc::clone(n), ln_inv));
        }
    }
    let Expr::Mul(fs) = expr.as_ref() else {
        return None;
    };
    if fs.len() != 2 {
        return None;
    }
    let ln_inv = |e: &ExprArc| {
        matches!(
            e.as_ref(),
            Expr::Pow(base, exp)
                if matches!(exp.as_ref(), Expr::Int(n) if n == &-num_bigint::BigInt::from(1))
                    && super::mrv_w::is_neg_ln_w_expr(base)
        )
    };
    if ln_inv(&fs[0]) {
        return Some((Arc::clone(&fs[1]), Arc::clone(&fs[0])));
    }
    if ln_inv(&fs[1]) {
        return Some((Arc::clone(&fs[0]), Arc::clone(&fs[1])));
    }
    None
}

fn combine_lead_with_ln_inv(
    lead: &MrvLeadTerm,
    ln_inv: &ExprArc,
    ctx: &Context,
) -> MrvLeadTerm {
    let den = match ln_inv.as_ref() {
        Expr::Pow(base, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) =>
        {
            Arc::clone(base)
        }
        _ => Arc::clone(ln_inv),
    };
    MrvLeadTerm {
        exponent: lead.exponent,
        coeff: divide_lead_coeffs(&lead.coeff, &den, ctx),
    }
}

fn is_series_coeff_undef(c: &ExprArc) -> bool {
    matches!(c.as_ref(), Expr::Undefined) || expr_contains_ln_w(c)
}

fn normalize_series_coeff(c: &ExprArc, ctx: &Context) -> ExprArc {
    let c = remove_lnexp(c, ctx);
    ratnormal(c.as_ref(), ctx).unwrap_or(c)
}

fn pnormal_series(p: &SparseSeries, ctx: &Context) -> SparseSeries {
    p.map_coeffs(|c| ratnormal(c.as_ref(), ctx).unwrap_or_else(|_| Arc::clone(c)))
}

/// upstream `mrv_lead_term` ordre loop: `series__SPOL1`, `ln(w)→±g`, `spdiv`.
fn mrv_series_lead_loop(
    swapped: &ExprArc,
    w: &Ident,
    g: &ExprArc,
    dont_invert: bool,
    omega_is_exp: bool,
    begin_ordre: usize,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    if let Some((core, ln_inv)) = peel_neg_ln_w_inv(swapped, w) {
        return lead_from_peeled_core(&core, &ln_inv, w, g, dont_invert, begin_ordre, ctx);
    }
    let mut f = Arc::clone(swapped);
    if !dont_invert {
        f = subst_w_inv(&f, w)?;
    }
    if omega_is_exp {
        f = rewrite_ln_exp_in_f(&f, w, g, dont_invert, begin_ordre, ctx)?;
    }
    if let Some((core, ln_inv)) = peel_neg_ln_w_inv(&f, w) {
        return lead_from_peeled_core(&core, &ln_inv, w, g, dont_invert, begin_ordre, ctx);
    }
    mrv_series_lead_loop_inner(&f, w, g, dont_invert, begin_ordre, ctx)
}

fn lead_from_peeled_core(
    core: &ExprArc,
    ln_inv: &ExprArc,
    w: &Ident,
    g: &ExprArc,
    dont_invert: bool,
    begin_ordre: usize,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    let core = ratnormal(remove_lnexp(core, ctx).as_ref(), ctx).unwrap_or_else(|_| remove_lnexp(core, ctx));
    if let Ok(s) = series_at_zero_order(&core, w, begin_ordre, MAX_SERIES_EXPANSION_ORDER, ctx) {
        if let Some((exp, coeff)) = s.lead() {
            if !mrv_is_zero(&coeff) {
                let lead = MrvLeadTerm {
                    exponent: exp,
                    coeff: normalize_series_coeff(&coeff, ctx),
                };
                return Ok(combine_lead_with_ln_inv(&lead, ln_inv, ctx));
            }
        }
    }
    let lead = mrv_series_lead_loop_inner(&core, w, g, dont_invert, begin_ordre, ctx)?;
    Ok(combine_lead_with_ln_inv(&lead, ln_inv, ctx))
}

fn mrv_series_lead_loop_inner(
    f: &ExprArc,
    w: &Ident,
    g: &ExprArc,
    dont_invert: bool,
    begin_ordre: usize,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    let order_cap = MAX_SERIES_EXPANSION_ORDER;
    let g_subst = if dont_invert {
        Arc::clone(g)
    } else {
        Expr::mul(vec![Expr::int(-1), Arc::clone(g)])
    };
    let mut ordre = begin_ordre as f64;

    while ordre < order_cap as f64 {
        let try_ord = (ordre as usize).min(order_cap);
        let mut inv = false;

        let mut p = match series_at_zero_order(f, w, try_ord, order_cap, ctx) {
            Ok(s) => s,
            Err(EvalError::NotImplemented(_)) => {
                ordre = ordre * 1.5 + 1.0;
                continue;
            }
            Err(e) => return Err(e),
        };

        if p.is_zero() {
            return Ok(MrvLeadTerm {
                coeff: Expr::int(0),
                exponent: 0,
            });
        }

        if let Some((lead_exp, coeff)) = p.lead() {
            let coeff = normalize_series_coeff(&coeff, ctx);
            let needs_inv = is_series_coeff_undef(&coeff) || !lead_coeff_ready(&coeff, w);
            if needs_inv {
                let substituted = subst_ln_w_expr(&coeff, &g_subst);
                let tmp = ratnormal(substituted.as_ref(), ctx).unwrap_or(substituted);
                if is_series_coeff_undef(&tmp) || !lead_coeff_ready(&tmp, w) {
                    inv = true;
                    p = series_spdiv_one(&p, try_ord, order_cap, ctx)?;
                    p = pnormal_series(&p, ctx);
                }
            }
        }

        p = p.map_coeffs(|c| subst_ln_w_expr(c, &g_subst));

        if inv {
            p = series_spdiv_one(&p, try_ord, order_cap, ctx)?;
            p = pnormal_series(&p, ctx);
        }

        if let Some((exp, coeff)) = p.lead() {
            let coeff = normalize_series_coeff(&coeff, ctx);
            if !is_series_coeff_undef(&coeff) && lead_coeff_ready(&coeff, w) && !mrv_is_zero(&coeff) {
                return Ok(MrvLeadTerm {
                    coeff,
                    exponent: exp,
                });
            }
        }

        ordre = ordre * 1.5 + 1.0;
    }
    Err(EvalError::NotImplemented("series"))
}

fn lead_coeff_ready(coeff: &ExprArc, w: &Ident) -> bool {
    if expr_contains_ln_w(coeff) {
        return false;
    }
    !contains_w(coeff) && !depends_on_w(coeff, w)
}

fn depends_on_w(e: &ExprArc, w: &Ident) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id == w,
        Expr::Add(ts) => ts.iter().any(|t| depends_on_w(t, w)),
        Expr::Mul(fs) => fs.iter().any(|f| depends_on_w(f, w)),
        Expr::Pow(b, exp) => depends_on_w(b, w) || depends_on_w(exp, w),
        Expr::Frac(n, d) => depends_on_w(n, w) || depends_on_w(d, w),
        Expr::Func(_, args) => args.iter().any(|a| depends_on_w(a, w)),
        _ => false,
    }
}

fn contains_ln_w(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && contains_w(&args[0]) => true,
        Expr::Add(ts) => ts.iter().any(contains_ln_w),
        Expr::Mul(fs) => fs.iter().any(contains_ln_w),
        Expr::Pow(b, exp) => contains_ln_w(b) || contains_ln_w(exp),
        Expr::Frac(n, d) => contains_ln_w(n) || contains_ln_w(d),
        Expr::Func(_, args) => args.iter().any(contains_ln_w),
        _ => false,
    }
}

fn contains_w(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id.as_str() == MRV_W,
        Expr::Add(ts) => ts.iter().any(contains_w),
        Expr::Mul(fs) => fs.iter().any(contains_w),
        Expr::Pow(b, exp) => contains_w(b) || contains_w(exp),
        Expr::Frac(n, d) => contains_w(n) || contains_w(d),
        Expr::Func(_, args) => args.iter().any(contains_w),
        _ => false,
    }
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn sign_infinity(coeff: &ExprArc) -> ExprArc {
    match coeff.as_ref() {
        Expr::Int(n) if n.is_negative() => Expr::sym("-infinity"),
        Expr::Int(n) if n.is_zero() => Expr::int(0),
        Expr::Frac(num, _) if matches!(num.as_ref(), Expr::Int(n) if n.is_negative()) => {
            Expr::sym("-infinity")
        }
        _ => Expr::sym("+infinity"),
    }
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn mrv_series_expansion_order_cap() {
        assert_eq!(MAX_SERIES_ORDER, 10);
        assert!(
            MAX_SERIES_EXPANSION_ORDER > MAX_SERIES_ORDER,
            "MRV ordre loop must escalate beyond default Taylor cap"
        );
    }

    #[test]
    fn limit_seven_pow_n_over_eight() {
        let ctx = xcas_default();
        let var = Ident::new("n");
        let e = Arc::new(Expr::Frac(
            Expr::pow(Expr::int(7), Expr::sym("n")),
            Expr::pow(Expr::int(8), Expr::sym("n")),
        ));
        let r = limit_unidirectional_plus_infinity(&e, &var, &ctx).unwrap();
        let got = format_expr(r.as_ref());
        assert_eq!(got, "0", "got {got}");
    }

    #[test]
    fn mrv_lead_ck_int_61() {
        let ctx = xcas_default();
        let var = Ident::new("x");
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
        let lead = mrv_lead_term_plus_infinity(&e, &var, &ctx).expect("mrv lead");
        let lim = limit_from_mrv_lead_term(&lead, &var, &ctx).expect("limit from mrv");
        assert_eq!(format_expr(lim.as_ref()), "-exp(2)");
    }
}
