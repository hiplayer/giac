//! `mrv_lead_term` and limit from MRV series (GIAC-216c).

use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, expr_to_poly, poly_to_expr, Context, EvalError, Expr, ExprArc,
    FuncKind, Ident,
};
use giac_simplify::ratnormal;
use num_traits::{Signed, Zero};

use super::bounds::{mrv_rewrite_bounded, mrv_series_eligible, MAX_SERIES_ORDER};
use super::mrv::{choose_mrv_w, mrv_at_plus_infinity};
use super::mrv_w::{
    decompose_ln_w_coeff, expr_contains_ln_w, is_expr_zero as mrv_is_zero, is_mrv_w_var, MRV_W,
};
use super::mrv_series_lead::mrv_lead_term_at_zero;
use super::mrv_series_lead::normalize_expr_quotients;
use super::preprocess::limit_preprocess_plus_infinity;
use super::sparse_series::series_at_zero;

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
    if !mrv_rewrite_bounded(&swapped) {
        return Err(EvalError::NotImplemented("limit"));
    }
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

/// Replace `exp(±x)` with `w`/`w^-1` and `x` with `-ln(w)` when `omega = exp(-x)`.
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

fn series_lead_at_zero(
    expr: &ExprArc,
    w: &Ident,
    order: usize,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    if is_mrv_w_var(w) && needs_exp_difference_cancellation(expr) {
        return mrv_lead_fallback(expr, w, ctx);
    }
    let series = series_at_zero(expr, w, order, ctx)?;
    let (exp, coeff) = series.lead().ok_or(EvalError::NotImplemented("limit"))?;
    Ok(MrvLeadTerm {
        coeff: ratnormal(coeff.as_ref(), ctx).unwrap_or(coeff),
        exponent: exp,
    })
}

fn mrv_lead_fallback(
    expr: &ExprArc,
    w: &Ident,
    ctx: &Context,
) -> Result<MrvLeadTerm, EvalError> {
    let (exp, coeff) = mrv_lead_term_at_zero(expr, w, ctx)?;
    let coeff = ratnormal(coeff.as_ref(), ctx).unwrap_or(coeff);
    let coeff = eval(coeff.as_ref(), ctx).unwrap_or(coeff);
    Ok(MrvLeadTerm {
        coeff: ratnormal(coeff.as_ref(), ctx).unwrap_or(coeff),
        exponent: exp,
    })
}

/// `exp(f) - w^{-1}` style cancellation (CK-INT-61); full `SparseSeries` on `Add` blows up.
fn needs_exp_difference_cancellation(expr: &ExprArc) -> bool {
    if add_has_exp_w_inv_difference(expr) {
        return true;
    }
    match expr.as_ref() {
        Expr::Frac(n, _) => add_has_exp_w_inv_difference(n),
        Expr::Mul(fs) => fs.iter().any(needs_exp_difference_cancellation),
        _ => false,
    }
}

fn add_has_exp_w_inv_difference(expr: &ExprArc) -> bool {
    let Expr::Add(terms) = expr.as_ref() else {
        return false;
    };
    if terms.len() != 2 {
        return false;
    }
    let has_exp = terms.iter().any(|t| {
        matches!(t.as_ref(), Expr::Func(FuncKind::Exp, args) if args.len() == 1)
    });
    let has_w_inv = terms.iter().any(is_neg_w_inv);
    has_exp && has_w_inv
}

fn is_neg_w_inv(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Pow(base, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
                && matches!(base.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)) =>
        {
            true
        }
        Expr::Mul(fs)
            if fs.len() == 2
                && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
                && fs.iter().any(|f| {
                    matches!(f.as_ref(), Expr::Pow(b, e)
                        if matches!(e.as_ref(), Expr::Int(n) if n.is_negative())
                            && matches!(b.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)))
                }) =>
        {
            true
        }
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
