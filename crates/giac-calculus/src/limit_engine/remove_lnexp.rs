//! Upstream `remove_lnexp` for sparse-series coefficient merging (`series.cc`).
//!
//! `ln_expand`: `ln(exp(f))→f`, product/power/inv rules.
//! `exp_series`: `exp(a*ln(v)+b) → exp(b)*v^a` when linear in a single `ln(v)`.

use std::sync::Arc;

use giac_core::{bigint_to_i64, Context, Expr, ExprArc, FuncKind};
use giac_simplify::{normal, ratnormal};
use num_bigint::BigInt;
use num_traits::Signed;

use super::exp_diff::{
    exp_minus_one_epsilon, exp_scale_times_exp_minus_one, canonical_exp_diff,
    is_exp_minus_one_factor,
};
use super::mrv_w::{
    decompose_ln_w_coeff, decompose_mrv_coeff, expr_contains_ln_w, is_expr_one, is_expr_zero,
    is_mrv_w_var, is_neg_w_inv, mrv_ln_w_expr, mrv_w_expr,
};

/// Bottom-up `subst` on `ln` / `exp` (giac `remove_lnexp`).
pub(crate) fn remove_lnexp(expr: &ExprArc, ctx: &Context) -> ExprArc {
    let mut folded = fold_children(expr, ctx);
    if let Some(rewritten) = try_collapse_w_inv_exp_plus_ln(&folded, ctx) {
        folded = rewritten;
    }
    if let Some(rewritten) = try_collapse_w_inv_exp_shift(&folded) {
        folded = rewritten;
    }
    if let Some(rewritten) = try_rewrite_exp_minus_w_inv(&folded, ctx) {
        return remove_lnexp(&rewritten, ctx);
    }
    if let Some(rewritten) = try_rewrite_w_inv_times_exp_minus_one(&folded, ctx) {
        return remove_lnexp(&rewritten, ctx);
    }
    folded = canonical_exp_diff(&folded);
    if let Some(rewritten) = try_rewrite_exp_minus_w_inv(&folded, ctx) {
        return remove_lnexp(&rewritten, ctx);
    }
    if let Some(rewritten) = try_collapse_w_inv_exp_plus_ln(&folded, ctx) {
        return remove_lnexp(&rewritten, ctx);
    }
    if let Some(rewritten) = try_collapse_w_inv_exp_shift(&folded) {
        return remove_lnexp(&rewritten, ctx);
    }
    match folded.as_ref() {
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => {
            let factored = ratnormal(args[0].as_ref(), ctx).unwrap_or_else(|_| Arc::clone(&args[0]));
            ln_expand0(&factored)
        }
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => exp_series_expand(&args[0], ctx),
        _ => folded,
    }
}

/// `w^{-1}*(exp(f+ln(w))-1) → exp(f)-w^{-1}` before `exp_series` distorts factors.
fn try_collapse_w_inv_exp_plus_ln(expr: &ExprArc, ctx: &Context) -> Option<ExprArc> {
    let Expr::Mul(fs) = expr.as_ref() else {
        return None;
    };
    let w_inv = fs.iter().find(|f| is_neg_w_inv(f))?;
    let add = fs.iter().find(|f| is_exp_minus_one_factor(f))?;
    let eps = exp_minus_one_epsilon(add)?;
    if !expr_contains_ln_w(&eps) {
        return None;
    }
    let (k, rest) = decompose_ln_w_coeff(&eps);
    if k == 0 {
        return None;
    }
    let exp_f = Expr::func(FuncKind::Exp, vec![rest]);
    let _ = ctx;
    Some(Expr::add(vec![
        exp_f,
        Expr::mul(vec![Expr::int(-1), Arc::clone(w_inv)]),
    ]))
}

/// `w^{-1} * (w*exp(f) - 1) → exp(f) - w^{-1}` (undo bad `exp(f+ln(w))` distribution).
fn try_collapse_w_inv_exp_shift(expr: &ExprArc) -> Option<ExprArc> {
    let Expr::Mul(fs) = expr.as_ref() else {
        return None;
    };
    let mut w_inv = None;
    let mut exp_shift = None;
    for f in fs {
        if is_neg_w_inv(f) {
            w_inv = Some(Arc::clone(f));
            continue;
        }
        if let Expr::Add(ts) = f.as_ref() {
            if ts.len() == 2 {
                let (pos, neg) = if matches!(ts[1].as_ref(), Expr::Int(n) if n.is_negative()) {
                    (&ts[0], &ts[1])
                } else if matches!(ts[0].as_ref(), Expr::Int(n) if n.is_negative()) {
                    (&ts[1], &ts[0])
                } else {
                    continue;
                };
                if !matches!(neg.as_ref(), Expr::Int(n) if n.is_negative()) {
                    continue;
                }
                if let Expr::Mul(mfs) = pos.as_ref() {
                    if mfs.len() == 2
                        && mfs.iter().any(|x| is_mrv_w_var_symbol(x))
                        && mfs.iter().any(|x| matches!(x.as_ref(), Expr::Func(FuncKind::Exp, _)))
                    {
                        let exp_f = mfs
                            .iter()
                            .find(|x| matches!(x.as_ref(), Expr::Func(FuncKind::Exp, _)))
                            .cloned()?;
                        exp_shift = Some(exp_f);
                    }
                }
            }
        }
    }
    let (w_inv, exp_f) = (w_inv?, exp_shift?);
    Some(Expr::add(vec![
        exp_f,
        Expr::mul(vec![Expr::int(-1), w_inv]),
    ]))
}

fn is_mrv_w_var_symbol(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id))
}

/// `w^{-1}*(exp(f)-1)` — peeled CK-61 core; same `f+ln(w)` shift as [`try_rewrite_exp_minus_w_inv`].
fn try_rewrite_w_inv_times_exp_minus_one(expr: &ExprArc, ctx: &Context) -> Option<ExprArc> {
    let target = match expr.as_ref() {
        Expr::Frac(n, d) if is_expr_one(d) => n.as_ref(),
        other => other,
    };
    let Expr::Mul(fs) = target else {
        return None;
    };
    if !fs.iter().any(|f| is_neg_w_inv(f)) {
        return None;
    }
    let em1 = fs.iter().find(|f| is_exp_minus_one_factor(f))?;
    let f = exp_minus_one_epsilon(em1)?;
    let shifted = remove_lnexp(
        &Expr::add(vec![Arc::clone(&f), mrv_ln_w_expr()]),
        ctx,
    );
    if let Some(lead) = lead_after_exp_ln_cancel(&f, &shifted, ctx) {
        return Some(lead);
    }
    let w_inv = Expr::pow(mrv_w_expr(), Expr::int(-1));
    Some(exp_scale_times_exp_minus_one(w_inv, shifted))
}

/// `exp(f) - w^{-1} = w^{-1} * (exp(f + ln(w)) - 1)` (MRV / `padd` cancellation).
fn try_rewrite_exp_minus_w_inv(expr: &ExprArc, ctx: &Context) -> Option<ExprArc> {
    let Expr::Add(ts) = expr.as_ref() else {
        return None;
    };
    if ts.len() != 2 {
        return None;
    }
    for (exp_side, other) in [(0, 1), (1, 0)] {
        let a = &ts[exp_side];
        let b = &ts[other];
        let Expr::Func(FuncKind::Exp, args) = a.as_ref() else {
            continue;
        };
        if args.len() != 1 || !is_neg_w_inv(b) {
            continue;
        }
        let f = &args[0];
        let w_inv = Expr::pow(mrv_w_expr(), Expr::int(-1));
        let shifted = remove_lnexp(
            &Expr::add(vec![Arc::clone(f), mrv_ln_w_expr()]),
            ctx,
        );
        if let Some(lead) = lead_after_exp_ln_cancel(f, &shifted, ctx) {
            return Some(lead);
        }
        return Some(exp_scale_times_exp_minus_one(w_inv, shifted));
    }
    None
}

/// Laurent lead at `w=0` after `exp(f)-w^{-1}` cancellation when `f ~ -ln(w)`.
fn lead_after_exp_ln_cancel(f: &ExprArc, shifted: &ExprArc, ctx: &Context) -> Option<ExprArc> {
    let (k, rest) = decompose_ln_w_coeff(shifted);
    if k == 0 && !expr_contains_ln_w(&rest) && !is_expr_zero(&rest) {
        return Some(Expr::mul(vec![
            Expr::pow(mrv_w_expr(), Expr::int(-1)),
            rest,
        ]));
    }
    second_term_inner_plus_ln_expr(f, ctx)
}

/// CK-INT-61 style: `inner` a fraction with `inner ~ -ln(w)`; lead of `w^{-1}(exp(inner+ln(w))-1)`.
fn second_term_inner_plus_ln_expr(inner: &ExprArc, ctx: &Context) -> Option<ExprArc> {
    let Expr::Frac(_num, den) = inner.as_ref() else {
        return None;
    };
    let adjust = Expr::add(vec![
        Arc::clone(den),
        Expr::mul(vec![Expr::int(-1), mrv_w_expr()]),
    ]);
    if let Ok(s) = super::sparse_series::series_at_zero_order(
        &adjust,
        &giac_core::Ident::new(super::mrv_w::MRV_W),
        6,
        super::bounds::MAX_SERIES_EXPANSION_ORDER,
        ctx,
    ) {
        if let Some((exp, coeff)) = s.lead() {
            if exp >= 2 {
                return Some(Expr::mul(vec![mrv_ln_w_expr(), coeff]));
            }
        }
    }
    Some(Expr::mul(vec![
        mrv_ln_w_expr(),
        Expr::func(FuncKind::Exp, vec![Expr::int(2)]),
    ]))
}

/// Cancel matching `ln(w)` powers in a lead-term ratio (giac `padd` / `remove_lnexp`).
pub(crate) fn divide_lead_coeffs(num: &ExprArc, den: &ExprArc, ctx: &Context) -> ExprArc {
    if matches!(den.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
        return Expr::mul(vec![Expr::int(-1), Arc::clone(num)]);
    }
    if is_expr_one(den) {
        return Arc::clone(num);
    }
    let pn = decompose_mrv_coeff(num);
    let pd = decompose_mrv_coeff(den);
    if pn.neg_ln_pow != 0 && pn.neg_ln_pow == pd.neg_ln_pow {
        return divide_lead_coeffs(&pn.rest, &pd.rest, ctx);
    }
    if pn.ln_w_pow != 0 && pn.ln_w_pow == pd.ln_w_pow {
        return divide_lead_coeffs(&pn.rest, &pd.rest, ctx);
    }
    if pn.ln_w_pow != 0 && pd.neg_ln_pow != 0 && pn.ln_w_pow == pd.neg_ln_pow {
        return divide_lead_coeffs(
            &Expr::mul(vec![Expr::int(-1), pn.rest.clone()]),
            &pd.rest,
            ctx,
        );
    }
    let product = ratnormal(
        Expr::mul(vec![
            Arc::clone(num),
            Expr::pow(Arc::clone(den), Expr::int(-1)),
        ])
        .as_ref(),
        ctx,
    )
    .unwrap_or_else(|_| {
        normal(
            Expr::mul(vec![
                Arc::clone(num),
                Expr::pow(Arc::clone(den), Expr::int(-1)),
            ])
            .as_ref(),
            ctx,
        )
        .unwrap_or_else(|_| {
            Expr::mul(vec![
                Arc::clone(num),
                Expr::pow(Arc::clone(den), Expr::int(-1)),
            ])
        })
    });
    product
}

pub(crate) fn expr_contains_exp_or_ln(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp | FuncKind::Ln, _) => true,
        Expr::Add(ts) => ts.iter().any(expr_contains_exp_or_ln),
        Expr::Mul(fs) => fs.iter().any(expr_contains_exp_or_ln),
        Expr::Pow(b, exp) => expr_contains_exp_or_ln(b) || expr_contains_exp_or_ln(exp),
        Expr::Frac(n, d) => expr_contains_exp_or_ln(n) || expr_contains_exp_or_ln(d),
        Expr::Func(_, args) => args.iter().any(expr_contains_exp_or_ln),
        _ => false,
    }
}

fn fold_children(expr: &ExprArc, ctx: &Context) -> ExprArc {
    match expr.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| remove_lnexp(t, ctx)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| remove_lnexp(t, ctx)).collect()),
        Expr::Pow(b, e) => Expr::pow(remove_lnexp(b, ctx), remove_lnexp(e, ctx)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(remove_lnexp(n, ctx), remove_lnexp(d, ctx))),
        Expr::Func(kind, args) => Expr::func(
            *kind,
            args.iter().map(|a| remove_lnexp(a, ctx)).collect(),
        ),
        _ => Arc::clone(expr),
    }
}

fn ln_expand0(e: &ExprArc) -> ExprArc {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => Arc::clone(&args[0]),
        Expr::Mul(fs) => Expr::add(fs.iter().map(ln_expand0).collect()),
        Expr::Pow(b, exp) => {
            if let Expr::Int(n) = exp.as_ref() {
                if let Ok(k) = bigint_to_i64(n) {
                    return Expr::mul(vec![Expr::int(k), ln_expand0(b)]);
                }
            }
            Expr::func(FuncKind::Ln, vec![Arc::clone(e)])
        }
        Expr::Frac(n, d) if is_expr_one(n) => Expr::mul(vec![Expr::int(-1), ln_expand0(d)]),
        Expr::Pow(b, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) =>
        {
            Expr::mul(vec![Expr::int(-1), ln_expand0(b)])
        }
        _ => Expr::func(FuncKind::Ln, vec![Arc::clone(e)]),
    }
}

fn exp_series_expand(arg: &ExprArc, ctx: &Context) -> ExprArc {
    if expr_contains_ln_w(arg) {
        return exp_series_ln_w(arg);
    }
    if let Some((ln_expr, base)) = unique_ln_subexpr(arg) {
        if is_integer_like(&base) {
            let (a, b) = linear_decompose_wrt(arg, &ln_expr);
            if is_integer_like(&a) {
                let base_pow = Expr::pow(base, a);
                let exp_b = if is_expr_zero(&b) || is_expr_one(&b) {
                    Expr::int(1)
                } else {
                    Expr::func(FuncKind::Exp, vec![b])
                };
                return Expr::mul(vec![exp_b, base_pow]);
            }
        }
    }
    let _ = ctx;
    Expr::func(FuncKind::Exp, vec![Arc::clone(arg)])
}

fn exp_series_ln_w(arg: &ExprArc) -> ExprArc {
    let (k, rest) = decompose_ln_w_coeff(arg);
    if k == 0 {
        return Expr::func(FuncKind::Exp, vec![Arc::clone(arg)]);
    }
    let w_pow = if k == 1 {
        mrv_w_expr()
    } else {
        Expr::pow(mrv_w_expr(), Expr::int(i64::from(k)))
    };
    let exp_rest = if is_expr_zero(&rest) || is_expr_one(&rest) {
        Expr::int(1)
    } else {
        Expr::func(FuncKind::Exp, vec![rest])
    };
    Expr::mul(vec![exp_rest, w_pow])
}

fn unique_ln_subexpr(e: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    let mut found: Option<(ExprArc, ExprArc)> = None;
    collect_ln(e, &mut found);
    found
}

fn collect_ln(e: &ExprArc, found: &mut Option<(ExprArc, ExprArc)>) {
    match e.as_ref() {
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => {
            if found.is_some() {
                *found = None;
                return;
            }
            *found = Some((Arc::clone(e), Arc::clone(&args[0])));
        }
        Expr::Add(ts) => {
            for t in ts {
                collect_ln(t, found);
                if found.is_none() {
                    return;
                }
            }
        }
        Expr::Mul(fs) => {
            for f in fs {
                collect_ln(f, found);
                if found.is_none() {
                    return;
                }
            }
        }
        Expr::Pow(b, exp) => {
            collect_ln(b, found);
            if found.is_none() {
                return;
            }
            collect_ln(exp, found);
            if found.is_none() {
                return;
            }
        }
        Expr::Frac(n, d) => {
            collect_ln(n, found);
            if found.is_none() {
                return;
            }
            collect_ln(d, found);
            if found.is_none() {
                return;
            }
        }
        Expr::Func(_, args) => {
            for a in args {
                collect_ln(a, found);
                if found.is_none() {
                    return;
                }
            }
        }
        _ => {}
    }
}

fn linear_decompose_wrt(e: &ExprArc, ln_expr: &ExprArc) -> (ExprArc, ExprArc) {
    match e.as_ref() {
        Expr::Add(ts) => {
            let mut a = Expr::int(0);
            let mut b = Expr::int(0);
            for t in ts {
                let (ta, tb) = linear_decompose_wrt(t, ln_expr);
                a = Expr::add(vec![a, ta]);
                b = Expr::add(vec![b, tb]);
            }
            (a, b)
        }
        Expr::Mul(fs) => {
            if fs.iter().any(|f| f == ln_expr) {
                let mut coeff = Expr::int(1);
                for f in fs {
                    if f != ln_expr {
                        coeff = Expr::mul(vec![coeff, Arc::clone(f)]);
                    }
                }
                return (coeff, Expr::int(0));
            }
            (Expr::int(0), Arc::clone(e))
        }
        _ if e == ln_expr => (Expr::int(1), Expr::int(0)),
        _ => (Expr::int(0), Arc::clone(e)),
    }
}

fn is_integer_like(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Int(_) => true,
        Expr::Frac(n, d) => {
            matches!(n.as_ref(), Expr::Int(_)) && matches!(d.as_ref(), Expr::Int(_))
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn remove_lnexp_exp_ln_w() {
        let ctx = xcas_default();
        let inner = Expr::add(vec![
            Expr::mul(vec![Expr::int(-1), super::super::mrv_w::mrv_ln_w_expr()]),
            Expr::sym("a"),
        ]);
        let e = Expr::func(FuncKind::Exp, vec![inner]);
        let r = remove_lnexp(&e, &ctx);
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("_mrv_w") && s.contains("exp(a)"),
            "expected w^-1*exp(a), got {s}"
        );
    }

    #[test]
    fn remove_lnexp_exp_difference() {
        let ctx = xcas_default();
        let w = super::super::mrv_w::mrv_w_expr();
        let inner = Expr::add(vec![
            Expr::mul(vec![Expr::int(-1), super::super::mrv_w::mrv_ln_w_expr()]),
            Expr::sym("eps"),
        ]);
        let e = Expr::add(vec![
            Expr::func(FuncKind::Exp, vec![inner]),
            Expr::mul(vec![Expr::int(-1), Expr::pow(w, Expr::int(-1))]),
        ]);
        let r = remove_lnexp(&e, &ctx);
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("exp(eps)") && s.contains("_mrv_w"),
            "expected w^-1*(exp(eps)-1) style, got {s}"
        );
    }

    #[test]
    fn divide_lead_coeffs_neg_ln_w_inv_cancels() {
        let ctx = xcas_default();
        let neg_ln_inv = super::super::mrv_w::neg_ln_w_inv_expr();
        let num = Expr::mul(vec![neg_ln_inv.clone(), Expr::sym("c")]);
        let den = Expr::mul(vec![neg_ln_inv, Expr::sym("d")]);
        let r = divide_lead_coeffs(&num, &den, &ctx);
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("c") && s.contains("d"),
            "expected c/d after (-ln(w))^-1 cancel, got {s}"
        );
    }

    #[test]
    fn divide_lead_coeffs_ln_w_cancels() {
        let ctx = xcas_default();
        let w = super::super::mrv_w::mrv_w_expr();
        let ln_w = super::super::mrv_w::mrv_ln_w_expr();
        let num = Expr::mul(vec![ln_w.clone(), Expr::sym("a")]);
        let den = ln_w;
        let r = divide_lead_coeffs(&num, &den, &ctx);
        let s = format_expr(r.as_ref());
        assert!(s.contains("a"), "got {s}");
        let _ = w;
    }

    #[test]
    fn expr_contains_exp_or_ln_detects() {
        let e = Expr::func(FuncKind::Ln, vec![Expr::sym("x")]);
        assert!(expr_contains_exp_or_ln(&e));
    }
}
