//! `mrv_lead_term` and limit from MRV series (GIAC-216c / GIAC-216e).
//!
//! **API 分层：** [`giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md)
//! **专项契约：** [`limit-engine-expr-api.md`](../../../../../.doc/limit-engine-expr-api.md)
//!
//! | 层级 | 内容 |
//! |------|------|
//! | **Stable** | `mrv_lead_term_plus_infinity`、`limit_from_mrv_lead_term`、`limit_unidirectional_plus_infinity` |
//! | **Pipeline private** | `peel_neg_ln_w_inv`、`rewrite_in_mrv_w`、`mrv_series_lead_loop*` |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, Context, EvalError, Expr, ExprArc,
    FuncKind, Ident,
};
use giac_simplify::ratnormal;
use num_traits::{Signed, Zero};

use crate::expr_util::{depends_on_var, is_var};
use super::bounds::{
    mrv_limit_eligible, mrv_rewrite_bounded, MAX_SERIES_EXPANSION_ORDER, MAX_SERIES_ORDER,
};
use super::mrv::{choose_mrv_w, linear_coeff_in_var, mrv_at_plus_infinity, vanishes_faster_than_at_plus_infinity};
use super::exp_diff::{simplify_add_sum, unwrap_signed_frac, vanishes_at_plus_infinity};
use super::simplify_util::{is_negative_const_expr, simplify_limit_expr};
use super::mrv_w::{
    decompose_ln_w_coeff, decompose_mrv_coeff, expr_contains_ln_w, is_expr_zero as mrv_is_zero,
    expr_contains_w_var, is_mrv_w_var, is_neg_ln_first_power, neg_ln_w_expr, neg_ln_w_inv_expr,
    MRV_W,
};
use super::mrv_series_lead::normalize_expr_quotients;
use super::preprocess::limit_preprocess_mrv;
use super::remove_lnexp::{divide_lead_coeffs, remove_lnexp};
use super::sparse_series::{series_at_zero_order, series_spdiv_one, SparseSeries};

#[derive(Clone, Debug)]
pub(crate) struct MrvLeadTerm {
    pub coeff: ExprArc,
    pub exponent: i32,
}

// **Pipeline private** — mrv lead exp vanishing plus infinity
fn mrv_lead_exp_vanishing_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<MrvLeadTerm> {
    let Expr::Func(FuncKind::Exp, args) = expr.as_ref() else {
        return None;
    };
    if args.len() != 1 {
        return None;
    }
    let f = simplify_add_sum(&args[0]);
    if let Some((n, d)) = unwrap_signed_frac(&f) {
        if vanishes_faster_than_at_plus_infinity(&n, &d, var, ctx) {
            return Some(MrvLeadTerm {
                coeff: Expr::int(1),
                exponent: 0,
            });
        }
    }
    if vanishes_at_plus_infinity(&f, var) {
        return Some(MrvLeadTerm {
            coeff: Expr::int(1),
            exponent: 0,
        });
    }
    None
}

/// **Stable** — MRV lead 管线入口（`+∞`）
pub(crate) fn mrv_lead_term_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    if !mrv_limit_eligible(expr) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let pre = limit_preprocess_mrv(expr, var);
    let pre = ratnormal(pre.as_ref(), ctx).unwrap_or(pre);
    if let Some(lead) = mrv_lead_exp_vanishing_plus_infinity(&pre, var, ctx) {
        return Ok(lead);
    }
    let pre = upscale_while_var_in_mrv(&pre, var, ctx);
    let pre = normalize_expr_quotients(&pre);
    if !mrv_limit_eligible(&pre) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let set = mrv_at_plus_infinity(&pre, var, ctx);
    if set.is_empty() {
        if depends_on_var(&pre, var) {
            return Err(EvalError::NotImplemented("limit"));
        }
        return Ok(MrvLeadTerm {
            coeff: ratnormal(pre.as_ref(), ctx).unwrap_or(pre),
            exponent: 0,
        });
    }
    let (omega, _) = choose_mrv_w(&set, var, ctx).ok_or(EvalError::NotImplemented("limit"))?;
    let w = Ident::new(MRV_W);
    let swapped = rewrite_in_mrv_w(&pre, var, &omega, &w, ctx);
    if !mrv_rewrite_bounded(&swapped) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let g = g_from_omega(&omega, var, &w, ctx);
    let dont_invert = omega_tends_to_zero_at_plus_inf(&omega, var, ctx);
    let omega_is_exp = matches!(omega.as_ref(), Expr::Func(FuncKind::Exp, _));
    let swapped_canon = super::mrv_w::canonical_mrv_coeff(&swapped);
    if let Some((core, ln_inv)) = peel_neg_ln_w_inv(&swapped_canon, &w) {
        return lead_from_peeled_core(
            &core,
            &ln_inv,
            &w,
            &g,
            dont_invert,
            MAX_SERIES_ORDER,
            omega_is_exp,
            var,
            ctx,
        );
    }
    let swapped = simplify_limit_expr(&swapped, ctx);
    let swapped = ratnormal(swapped.as_ref(), ctx).unwrap_or(swapped);
    mrv_series_lead_loop(
        &swapped,
        &w,
        &g,
        dont_invert,
        omega_is_exp,
        MAX_SERIES_ORDER,
        var,
        ctx,
    )
}

/// **Stable** — 由 `MrvLeadTerm` 求极限
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
    if expr_contains_ln_w(&coeff) {
        coeff = rewrite_ln_w(&coeff, var);
        if depends_on_var(&coeff, var) {
            return Err(EvalError::NotImplemented("limit"));
        }
    }
    let coeff = eval(coeff.as_ref(), ctx)?;
    if expr_contains_w_var(&coeff) {
        return Err(EvalError::NotImplemented("limit"));
    }
    Ok(coeff)
}

/// **Stable** — 单向 `+∞` 极限（递归 lead）
pub(crate) fn limit_unidirectional_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if let Some(r) = super::mrv::limit_const_pow_quotient_at_plus_infinity(expr, var) {
        return Ok(r);
    }
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
        if !depends_on_var(&coeff, var) {
            return Ok(coeff);
        }
        e_copy = coeff;
    }
    Err(EvalError::NotImplemented("limit"))
}

// **Pipeline private** — subst symbol
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

// **Pipeline private** — upscale while var in mrv
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

// **Pipeline private** — omega neg linear coeff
fn omega_neg_linear_coeff(omega: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let Expr::Func(FuncKind::Exp, args) = omega.as_ref() else {
        return None;
    };
    if args.len() != 1 {
        return None;
    }
    let c = linear_coeff_in_var(&args[0], var)?;
    if is_negative_const_expr(&c) {
        Some(c)
    } else {
        None
    }
}

// **Pipeline private** — neg ln w from omega
fn neg_ln_w_from_omega(omega_lin: ExprArc, w_expr: ExprArc) -> ExprArc {
    if matches!(omega_lin.as_ref(), Expr::Int(n) if n == &-num_bigint::BigInt::from(1)) {
        return super::mrv_w::neg_ln_w_expr();
    }
    Expr::mul(vec![
        Expr::pow(omega_lin, Expr::int(-1)),
        Expr::func(FuncKind::Ln, vec![w_expr]),
    ])
}

// **Pipeline private** — linear inner coeff only
fn linear_inner_coeff_only(inner: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let c = linear_coeff_in_var(inner, var)?;
    let mut map = HashMap::new();
    map.insert(var.clone(), Expr::int(0));
    let at0 = eval_subst_map(inner, &map).ok()?;
    if mrv_is_zero(&at0) {
        Some(c)
    } else {
        None
    }
}

// **Pipeline private** — rewrite in mrv w
fn rewrite_in_mrv_w(
    expr: &ExprArc,
    var: &Ident,
    omega: &ExprArc,
    w: &Ident,
    ctx: &Context,
) -> ExprArc {
    let w_expr = Expr::sym(w.as_str());
    if expr_eq(expr, omega) {
        return w_expr;
    }
    if let Some(omega_lin) = omega_neg_linear_coeff(omega, var) {
        if is_var(expr, var) {
            return neg_ln_w_from_omega(omega_lin, w_expr);
        }
        if let Expr::Pow(base, exp) = expr.as_ref() {
            if is_var(base, var) {
                if let Expr::Int(_) = exp.as_ref() {
                    let x_as_ln = neg_ln_w_from_omega(omega_lin, w_expr);
                    return Expr::pow(x_as_ln, Arc::clone(exp));
                }
            }
        }
        if let Expr::Func(FuncKind::Exp, args) = expr.as_ref() {
            if args.len() == 1 {
                if let Some(e_lin) = linear_inner_coeff_only(&args[0], var) {
                    let exp = Arc::new(Expr::Frac(e_lin, omega_lin));
                    return Expr::pow(w_expr, exp);
                }
            }
        }
    }
    if is_var(expr, var) && is_exp_neg_var(omega, var) {
        return Expr::mul(vec![Expr::int(-1), Expr::func(FuncKind::Ln, vec![w_expr])]);
    }
    if is_exp_pos_var(expr, var) && is_exp_neg_var(omega, var) {
        return Expr::pow(w_expr, Expr::int(-1));
    }
    match expr.as_ref() {
        Expr::Add(ts) => Expr::add(
            ts.iter()
                .map(|t| rewrite_in_mrv_w(t, var, omega, w, ctx))
                .collect(),
        ),
        Expr::Mul(fs) => Expr::mul(
            fs.iter()
                .map(|t| rewrite_in_mrv_w(t, var, omega, w, ctx))
                .collect(),
        ),
        Expr::Pow(b, e) => Expr::pow(
            rewrite_in_mrv_w(b, var, omega, w, ctx),
            rewrite_in_mrv_w(e, var, omega, w, ctx),
        ),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            rewrite_in_mrv_w(n, var, omega, w, ctx),
            rewrite_in_mrv_w(d, var, omega, w, ctx),
        )),
        Expr::Func(kind, args) => Expr::func(
            *kind,
            args.iter()
                .map(|a| rewrite_in_mrv_w(a, var, omega, w, ctx))
                .collect(),
        ),
        _ => Arc::clone(expr),
    }
}

// **Pipeline private** — expr eq
fn expr_eq(a: &ExprArc, b: &ExprArc) -> bool {
    a == b
}

// **Pipeline private** — is exp neg var
fn is_exp_neg_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 && is_neg_var(&args[0], var)
    )
}

// **Pipeline private** — is exp pos var
fn is_exp_pos_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 && matches!(args[0].as_ref(), Expr::Symbol(id) if id == var)
    )
}

// **Pipeline private** — is neg var
fn is_neg_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Mul(fs) if fs.len() == 2
            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Symbol(id) if id == var))
    )
}

// **Pipeline private** — rewrite ln w
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

// **Pipeline private** — g from omega
fn g_from_omega(omega: &ExprArc, var: &Ident, w: &Ident, ctx: &Context) -> ExprArc {
    match omega.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            rewrite_in_mrv_w(&args[0], var, omega, w, ctx)
        }
        _ => Expr::int(0),
    }
}

// **Pipeline private** — omega tends to zero at plus inf
fn omega_tends_to_zero_at_plus_inf(omega: &ExprArc, var: &Ident, _ctx: &Context) -> bool {
    if omega_neg_linear_coeff(omega, var).is_some() {
        return true;
    }
    match omega.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            if is_neg_var_linear(&args[0], var) {
                return true;
            }
            linear_coeff_in_var(&args[0], var)
                .is_some_and(|c| is_negative_const_expr(&c))
        }
        _ => false,
    }
}

// **Pipeline private** — is neg var linear
fn is_neg_var_linear(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Mul(fs) if fs.len() == 2
            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
            && fs.iter().any(|f| is_var(f, var))
    )
}

// **Pipeline private** — subst w inv
fn subst_w_inv(expr: &ExprArc, w: &Ident) -> Result<ExprArc, EvalError> {
    let mut m = HashMap::new();
    m.insert(w.clone(), Expr::pow(Expr::sym(w.as_str()), Expr::int(-1)));
    eval_subst_map(expr, &m)
}

// **Pipeline private** — is ln w
fn is_ln_w(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Ln, args)
            if args.len() == 1 && matches!(args[0].as_ref(), Expr::Symbol(id) if is_mrv_w_var(id))
    )
}

// **Pipeline private** — subst ln w expr
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

// **Pipeline private** — collect ln exprs
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

// **Pipeline private** — rewrite ln exp in f
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

// **Pipeline private** — subst expr once
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

// **Pipeline private** — peel neg ln w inv
fn peel_neg_ln_w_inv(expr: &ExprArc, _w: &Ident) -> Option<(ExprArc, ExprArc)> {
    let inv = neg_ln_w_inv_expr();
    let neg_ln = neg_ln_w_expr();
    let e = super::mrv_w::canonical_mrv_coeff(expr);
    if let Expr::Frac(n, d) = e.as_ref() {
        let d = super::mrv_w::canonical_mrv_coeff(d);
        if is_neg_ln_first_power(&d) {
            return Some((Arc::clone(n), inv));
        }
        if d == inv {
            return Some((Expr::mul(vec![Arc::clone(n), neg_ln]), inv));
        }
    }
    if let Expr::Mul(fs) = e.as_ref() {
        if let Some(i) = fs.iter().position(|f| f == &inv) {
            let mut rest: Vec<_> = fs.to_vec();
            rest.remove(i);
            let core = if rest.is_empty() {
                Expr::int(1)
            } else if rest.len() == 1 {
                rest.pop().unwrap()
            } else {
                Expr::mul(rest)
            };
            return Some((core, inv));
        }
    }
    let parts = super::mrv_w::decompose_mrv_coeff(&e);
    if !parts.has_neg_ln_inv() {
        return None;
    }
    let core = if parts.neg_ln_pow == -1 {
        parts.rest
    } else {
        Expr::mul(vec![
            parts.rest,
            Expr::pow(
                neg_ln,
                Expr::int(i64::from(parts.neg_ln_pow + 1)),
            ),
        ])
    };
    Some((core, inv))
}

// **Pipeline private** — combine lead with ln inv
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

// **Pipeline private** — is series coeff undef
fn is_series_coeff_undef(c: &ExprArc) -> bool {
    matches!(c.as_ref(), Expr::Undefined) || decompose_mrv_coeff(c).ln_w_pow != 0
}

// **Pipeline private** — normalize series coeff
fn normalize_series_coeff(c: &ExprArc, ctx: &Context) -> ExprArc {
    simplify_limit_expr(&remove_lnexp(c, ctx), ctx)
}

// **Pipeline private** — pnormal series
fn pnormal_series(p: &SparseSeries, ctx: &Context) -> SparseSeries {
    p.map_coeffs(|c| ratnormal(c.as_ref(), ctx).unwrap_or_else(|_| Arc::clone(c)))
}

/// **Stable** — `w=0` 级数 lead
pub(crate) fn series_lead_at_zero(
    f: &ExprArc,
    w: &Ident,
    g: &ExprArc,
    dont_invert: bool,
    begin_ordre: usize,
    var: &Ident,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    mrv_series_lead_loop_inner(f, w, g, dont_invert, begin_ordre, var, ctx)
}

// **Pipeline private** — mrv series lead loop
fn mrv_series_lead_loop(
    swapped: &ExprArc,
    w: &Ident,
    g: &ExprArc,
    dont_invert: bool,
    omega_is_exp: bool,
    begin_ordre: usize,
    var: &Ident,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    let swapped = super::mrv_w::canonical_mrv_coeff(swapped);
    if let Some((core, ln_inv)) = peel_neg_ln_w_inv(&swapped, w) {
        return lead_from_peeled_core(
            &core,
            &ln_inv,
            w,
            g,
            dont_invert,
            begin_ordre,
            omega_is_exp,
            var,
            ctx,
        );
    }
    let mut f = remove_lnexp(&swapped, ctx);
    f = ratnormal(f.as_ref(), ctx).unwrap_or(f);
    if let Some((core, ln_inv)) = peel_neg_ln_w_inv(&f, w) {
        return lead_from_peeled_core(
            &core,
            &ln_inv,
            w,
            g,
            dont_invert,
            begin_ordre,
            omega_is_exp,
            var,
            ctx,
        );
    }
    if !dont_invert {
        f = subst_w_inv(&f, w)?;
    }
    if omega_is_exp {
        f = rewrite_ln_exp_in_f(&f, w, g, dont_invert, begin_ordre, ctx)?;
    }
    if let Some((core, ln_inv)) = peel_neg_ln_w_inv(&f, w) {
        return lead_from_peeled_core(
            &core,
            &ln_inv,
            w,
            g,
            dont_invert,
            begin_ordre,
            omega_is_exp,
            var,
            ctx,
        );
    }
    mrv_series_lead_loop_inner(&f, w, g, dont_invert, begin_ordre, var, ctx)
}

// **Pipeline private** — lead from peeled core
#[allow(clippy::too_many_arguments)] // ponytail: mirrors upstream peel/lead arity
fn lead_from_peeled_core(
    core: &ExprArc,
    ln_inv: &ExprArc,
    w: &Ident,
    g: &ExprArc,
    dont_invert: bool,
    begin_ordre: usize,
    omega_is_exp: bool,
    var: &Ident,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    let mut core = ratnormal(remove_lnexp(core, ctx).as_ref(), ctx).unwrap_or_else(|_| remove_lnexp(core, ctx));
    core = super::mrv_w::canonical_mrv_coeff(&core);
    let (k, rest) = decompose_ln_w_coeff(&core);
    if k != 0 && lead_coeff_ready(&rest, var, w) && !mrv_is_zero(&rest) {
        let lead = MrvLeadTerm {
            exponent: 0,
            coeff: normalize_series_coeff(&core, ctx),
        };
        return Ok(combine_lead_with_ln_inv(&lead, ln_inv, ctx));
    }
    if omega_is_exp {
        core = rewrite_ln_exp_in_f(&core, w, g, dont_invert, begin_ordre, ctx)?;
    }
    core = super::mrv_w::canonical_mrv_coeff(&core);
    let core = ratnormal(remove_lnexp(&core, ctx).as_ref(), ctx).unwrap_or_else(|_| remove_lnexp(&core, ctx));
    let (k, rest) = decompose_ln_w_coeff(&core);
    if k != 0 && lead_coeff_ready(&rest, var, w) && !mrv_is_zero(&rest) {
        let lead = MrvLeadTerm {
            exponent: 0,
            coeff: normalize_series_coeff(&core, ctx),
        };
        return Ok(combine_lead_with_ln_inv(&lead, ln_inv, ctx));
    }
    if let Ok(s) = series_at_zero_order(&core, w, begin_ordre, MAX_SERIES_EXPANSION_ORDER, ctx) {
        if let Some((exp, coeff)) = s.lead() {
            let coeff = normalize_series_coeff(&coeff, ctx);
            if !mrv_is_zero(&coeff) && lead_coeff_ready(&coeff, var, w) {
                let lead = MrvLeadTerm { exponent: exp, coeff };
                return Ok(combine_lead_with_ln_inv(&lead, ln_inv, ctx));
            }
        }
    }
    let lead = mrv_series_lead_loop_inner(&core, w, g, dont_invert, begin_ordre, var, ctx)?;
    Ok(combine_lead_with_ln_inv(&lead, ln_inv, ctx))
}

// **Pipeline private** — mrv series lead loop inner
fn mrv_series_lead_loop_inner(
    f: &ExprArc,
    w: &Ident,
    g: &ExprArc,
    dont_invert: bool,
    begin_ordre: usize,
    var: &Ident,
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

        if let Some((_lead_exp, coeff)) = p.lead() {
            let coeff = normalize_series_coeff(&coeff, ctx);
            let needs_inv = is_series_coeff_undef(&coeff) || !lead_coeff_ready(&coeff, var, w);
            if needs_inv {
                let substituted = subst_ln_w_expr(&coeff, &g_subst);
                let normalized = ratnormal(substituted.as_ref(), ctx).unwrap_or(substituted);
                if is_series_coeff_undef(&normalized) || !lead_coeff_ready(&normalized, var, w) {
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
            if !is_series_coeff_undef(&coeff) && lead_coeff_ready(&coeff, var, w) && !mrv_is_zero(&coeff) {
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

// **Pipeline private** — lead coeff ready
fn lead_coeff_ready(coeff: &ExprArc, var: &Ident, w: &Ident) -> bool {
    let parts = decompose_mrv_coeff(coeff);
    if parts.pending_for_series() {
        return false;
    }
    if depends_on_var(&parts.rest, var) {
        return false;
    }
    !expr_contains_w_var(&parts.rest) && !depends_on_w(&parts.rest, w)
}

// **Pipeline private** — depends on w
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

// **Pipeline private** — sign infinity
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
    use std::sync::Arc;

    use giac_core::{format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    #[allow(clippy::assertions_on_constants)]
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
    fn rewrite_exp_x_to_w_inv() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let w = Ident::new(MRV_W);
        let omega = Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]);
        let expx = Expr::func(FuncKind::Exp, vec![Expr::sym("x")]);
        let r = rewrite_in_mrv_w(&expx, &var, &omega, &w, &ctx);
        assert!(crate::limit_engine::mrv_w::is_neg_w_inv(&r));
    }

    #[test]
    fn rewrite_x_powers_in_mrv_w() {
        
        let ctx = xcas_default();
        let var = Ident::new("x");
        let w = Ident::new(MRV_W);
        let omega = Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]);
        let x1 = Expr::pow(Expr::sym("x"), Expr::int(1));
        let xm1 = Expr::pow(Expr::sym("x"), Expr::int(-1));
        let r1 = rewrite_in_mrv_w(&x1, &var, &omega, &w, &ctx);
        let rm1 = rewrite_in_mrv_w(&xm1, &var, &omega, &w, &ctx);
        use crate::limit_engine::mrv_w::{canonical_mrv_coeff, neg_ln_w_expr, neg_ln_w_inv_expr};
        assert_eq!(
            canonical_mrv_coeff(&r1),
            neg_ln_w_expr(),
            "x^1 -> -ln(w)"
        );
        assert_eq!(
            canonical_mrv_coeff(&rm1),
            neg_ln_w_inv_expr(),
            "x^-1 -> (-ln(w))^-1"
        );
    }

    #[test]
    fn rewrite_ck61_normed_in_mrv_w() {
        use crate::limit_engine::ck_int_gruntz_fixture::ck_int_61;
        use crate::limit_engine::preprocess::{limit_preprocess_mrv, limit_preprocess_struct};
        use crate::limit_engine::mrv::{choose_mrv_w, mrv_at_plus_infinity};
        let ctx = xcas_default();
        let var = Ident::new("x");
        let pre = limit_preprocess_struct(&ck_int_61(), &var);
        let mrv_pre = limit_preprocess_mrv(&pre, &var);
        let upscaled = upscale_while_var_in_mrv(&mrv_pre, &var, &ctx);
        let normed = normalize_expr_quotients(&upscaled);
        let set = mrv_at_plus_infinity(&normed, &var, &ctx);
        let (omega, _) = choose_mrv_w(&set, &var, &ctx).expect("omega");
        let w = Ident::new(MRV_W);
        let swapped = rewrite_in_mrv_w(&normed, &var, &omega, &w, &ctx);
        assert!(
            super::super::mrv_w::expr_contains_ln_w(&swapped),
            "got {}",
            format_expr(swapped.as_ref())
        );
        assert!(!depends_on_var(&swapped, &var), "x should be eliminated");
    }

    #[test]
    fn ratio_mrv_lead_exp_vanishing_after_preprocess() {
        use crate::limit_engine::ck_int_gruntz_fixture::ratio;
        use crate::limit_engine::preprocess::limit_preprocess_mrv;
        let ctx = xcas_default();
        let var = Ident::new("x");
        let pre = limit_preprocess_mrv(&ratio(), &var);
        let lead = mrv_lead_exp_vanishing_plus_infinity(&pre, &var, &ctx)
            .unwrap_or_else(|| panic!("expected vanishing exp lead, pre={}", format_expr(pre.as_ref())));
        assert_eq!(lead.exponent, 0);
        assert_eq!(format_expr(lead.coeff.as_ref()), "1");
    }

    #[test]
    fn mrv_lead_ck_int_60_ratio() {
        use crate::limit_engine::ck_int_gruntz_fixture::ratio;
        let ctx = xcas_default();
        let var = Ident::new("x");
        let r = limit_unidirectional_plus_infinity(&ratio(), &var, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1");
    }

    #[test]
    fn mrv_lead_ck_int_61() {
        use crate::limit_engine::ck_int_gruntz_fixture::ck_int_61;
        let ctx = xcas_default();
        let var = Ident::new("x");
        let lead = mrv_lead_term_plus_infinity(&ck_int_61(), &var, &ctx).expect("mrv lead");
        let lim = limit_from_mrv_lead_term(&lead, &var, &ctx).expect("limit from mrv");
        assert_eq!(format_expr(lim.as_ref()), "-exp(2)");
    }
}
