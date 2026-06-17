//! `mrv_lead_term` and limit from MRV series (GIAC-216c).

use std::sync::Arc;

use giac_core::{
    eval, ratnormal, Context, EvalError, Expr, ExprArc, FuncKind, Ident,
};
use num_traits::{Signed, Zero};

use super::bounds::{mrv_series_eligible, MAX_SERIES_ORDER};
use super::mrv::{choose_mrv_w, mrv_at_plus_infinity};
use super::preprocess::limit_preprocess_plus_infinity;
use super::sparse_series::series_at_zero;

const MRV_W: &str = "_mrv_w";

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
    if !mrv_series_eligible(expr) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let pre = limit_preprocess_plus_infinity(expr, var, ctx)?;
    if !mrv_series_eligible(&pre) {
        return Err(EvalError::NotImplemented("limit"));
    }
    let set = mrv_at_plus_infinity(&pre, var);
    if set.is_empty() {
        return Ok(MrvLeadTerm {
            coeff: ratnormal(pre.as_ref(), ctx).unwrap_or(pre),
            exponent: 0,
        });
    }
    let (omega, _) = choose_mrv_w(&set, var).ok_or(EvalError::NotImplemented("limit"))?;
    let w = Ident::new(MRV_W);
    let swapped = rewrite_in_mrv_w(&pre, var, &omega, &w);
    let swapped = ratnormal(swapped.as_ref(), ctx).unwrap_or(swapped);
    series_lead_at_zero(&swapped, &w, MAX_SERIES_ORDER, ctx)
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

/// Replace `exp(±x)` with `w` / `w^-1` only — do not rewrite bare `x` to `ln(w)`.
fn rewrite_in_mrv_w(expr: &ExprArc, var: &Ident, omega: &ExprArc, w: &Ident) -> ExprArc {
    let w_expr = Expr::sym(w.as_str());
    if expr_eq(expr, omega) {
        return w_expr;
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
    match expr.as_ref() {
        Expr::Func(FuncKind::Ln, args)
            if args.len() == 1
                && matches!(args[0].as_ref(), Expr::Symbol(id) if id.as_str() == MRV_W) =>
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

fn series_lead_at_zero(
    expr: &ExprArc,
    w: &Ident,
    order: usize,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    let series = series_at_zero(expr, w, order, ctx)?;
    let (exp, coeff) = series.lead().ok_or(EvalError::NotImplemented("limit"))?;
    Ok(MrvLeadTerm {
        coeff: ratnormal(coeff.as_ref(), ctx).unwrap_or(coeff),
        exponent: exp,
    })
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
