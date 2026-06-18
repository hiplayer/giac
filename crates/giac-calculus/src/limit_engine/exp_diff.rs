//! Shared `exp` difference normalization (GIAC-limit-exp-diff P1).
//!
//! Upstream `remove_lnexp` (`series.cc`) applies `ln_expand` / `exp_series` on coefficients;
//! the identity `exp(f) - w^{-1} = w^{-1}(exp(f + ln(w)) - 1)` is the MRV (`w`) layer.
//!
//! `fold_exp_shifted_difference` is the same ε→0 pattern on the original variable:
//! all shapes reduce to `exp(scale_log) * (exp(ε) - 1)` via [`detect_exp_difference`] + [`emit_exp_difference`].

use std::cmp::Ordering;
use std::sync::Arc;

use giac_core::{Context, Expr, ExprArc, FuncKind, Ident};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use crate::expr_util::depends_on_var;

/// Normal form: `exp(scale_log) * (exp(epsilon) - 1)` [× optional complement].
#[derive(Clone, Debug)]
struct ExpDiffForm {
    scale_log: ExprArc,
    epsilon: ExprArc,
    complement: Option<ExprArc>,
}

/// Bottom-up `exp` difference factorization (x-layer; used by preprocess).
pub(crate) fn fold_exp_shifted_difference(expr: &ExprArc) -> ExprArc {
    let rebuilt = match expr.as_ref() {
        Expr::Mul(fs) => Expr::mul(fs.iter().map(fold_exp_shifted_difference).collect()),
        Expr::Add(ts) => Expr::add(ts.iter().map(fold_exp_shifted_difference).collect()),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            fold_exp_shifted_difference(n),
            fold_exp_shifted_difference(d),
        )),
        Expr::Pow(b, e) => Expr::pow(
            fold_exp_shifted_difference(b),
            fold_exp_shifted_difference(e),
        ),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter().map(fold_exp_shifted_difference).collect(),
        ),
        _ => Arc::clone(expr),
    };
    detect_exp_difference(&rebuilt)
        .map(emit_exp_difference)
        .unwrap_or(rebuilt)
}

/// `scale * (exp(ε) - 1)` — common target form (P1).
pub(crate) fn exp_scale_times_exp_minus_one(scale: ExprArc, epsilon: ExprArc) -> ExprArc {
    Expr::mul(vec![
        scale,
        Expr::add(vec![Expr::func(FuncKind::Exp, vec![epsilon]), Expr::int(-1)]),
    ])
}

fn emit_exp_difference(form: ExpDiffForm) -> ExprArc {
    let core = exp_scale_times_exp_minus_one(
        Expr::func(FuncKind::Exp, vec![form.scale_log]),
        form.epsilon,
    );
    match form.complement {
        None => core,
        Some(c) => Expr::mul(vec![core, c]),
    }
}

/// Detect `exp(scale_log) * (exp(ε) - 1)` on an already-folded expression tree.
fn detect_exp_difference(expr: &ExprArc) -> Option<ExpDiffForm> {
    match expr.as_ref() {
        Expr::Mul(_) => detect_exp_difference_mul(&flatten_mul_factors(expr)),
        Expr::Add(ts) if ts.len() == 2 => detect_exp_difference_add(&ts[0], &ts[1]),
        _ => None,
    }
}

fn flatten_mul_factors(expr: &ExprArc) -> Vec<ExprArc> {
    match expr.as_ref() {
        Expr::Mul(fs) => {
            let mut out = Vec::new();
            for f in fs {
                out.extend(flatten_mul_factors(f));
            }
            out
        }
        _ => vec![Arc::clone(expr)],
    }
}

fn detect_exp_difference_mul(factors: &[ExprArc]) -> Option<ExpDiffForm> {
    let mut expanded = Vec::new();
    for f in factors {
        match f.as_ref() {
            Expr::Mul(inner) => {
                if let Some(form) = detect_exp_difference_mul(&flatten_mul_from_vec(inner.clone())) {
                    expanded.push(emit_exp_difference(form));
                } else {
                    expanded.push(Arc::clone(f));
                }
            }
            Expr::Add(ts) if ts.len() == 2 => {
                if let Some(form) = detect_exp_difference_add(&ts[0], &ts[1]) {
                    expanded.push(emit_exp_difference(form));
                } else {
                    expanded.push(Arc::clone(f));
                }
            }
            _ => expanded.push(Arc::clone(f)),
        }
    }
    let expanded = flatten_mul_from_vec(expanded);

    let mut minus_one_idx = None;
    for (i, f) in expanded.iter().enumerate() {
        if exp_minus_one_inner_arg(f).is_some() {
            minus_one_idx = Some(i);
            break;
        }
    }
    let mi = minus_one_idx?;
    let epsilon = exp_minus_one_inner_arg(&expanded[mi])?;
    let mut scale_logs = Vec::new();
    let mut complement = Vec::new();
    for (i, f) in expanded.iter().enumerate() {
        if i == mi {
            continue;
        }
        if let Some(log) = exp_func_arg(f) {
            scale_logs.push(log);
        } else {
            complement.push(Arc::clone(f));
        }
    }
    if scale_logs.is_empty() {
        return None;
    }
    Some(ExpDiffForm {
        scale_log: sum_exp_logs(scale_logs),
        epsilon,
        complement: match complement.len() {
            0 => None,
            1 => Some(complement.into_iter().next().unwrap()),
            _ => Some(Expr::mul(complement)),
        },
    })
}

fn flatten_mul_from_vec(factors: Vec<ExprArc>) -> Vec<ExprArc> {
    let mut out = Vec::new();
    for f in factors {
        if let Expr::Mul(inner) = f.as_ref() {
            out.extend(flatten_mul_from_vec(inner.clone()));
        } else {
            out.push(f);
        }
    }
    out
}

fn sum_exp_logs(logs: Vec<ExprArc>) -> ExprArc {
    match logs.len() {
        0 => Expr::int(0),
        1 => logs.into_iter().next().unwrap(),
        _ => Expr::add(logs),
    }
}

/// `exp(A) - exp(B) → exp(B) * (exp(ε) - 1)` with `ε = A - B`.
/// When `B` is an add subterm of `A`, `remove_add_term` keeps `ε` compact (L1 / gruntz).
fn detect_exp_difference_add(a: &ExprArc, b: &ExprArc) -> Option<ExpDiffForm> {
    let (pos_arg, neg_arg) = exp_signed_pair(a, b)?;
    let epsilon = remove_add_term(&pos_arg, &neg_arg)
        .unwrap_or_else(|| difference_add_exprs(&pos_arg, &neg_arg));
    Some(ExpDiffForm {
        scale_log: neg_arg,
        epsilon,
        complement: None,
    })
}

fn difference_add_exprs(a: &ExprArc, b: &ExprArc) -> ExprArc {
    let mut terms: Vec<(bool, ExprArc)> = signed_add_terms(a)
        .into_iter()
        .map(|(p, t)| (p, canonical_add_term(&t)))
        .collect();
    for (pos, t) in signed_add_terms(b) {
        let t = canonical_add_term(&t);
        let pos = !pos;
        if let Some(i) = terms.iter().position(|(p, tt)| *p == pos && tt == &t) {
            terms.remove(i);
        } else if let Some(i) = terms.iter().position(|(p, tt)| *p != pos && tt == &t) {
            terms.remove(i);
        } else {
            terms.push((pos, t));
        }
    }
    rebuild_signed_add(terms)
}

fn canonical_add_term(t: &ExprArc) -> ExprArc {
    match t.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 && matches!(fs[0].as_ref(), Expr::Int(n) if n.is_one()) => {
            canonical_add_term(&fs[1])
        }
        Expr::Mul(fs) if fs.len() == 2 && matches!(fs[1].as_ref(), Expr::Int(n) if n.is_one()) => {
            canonical_add_term(&fs[0])
        }
        _ => Arc::clone(t),
    }
}

fn signed_add_terms(e: &ExprArc) -> Vec<(bool, ExprArc)> {
    match e.as_ref() {
        Expr::Add(ts) => ts.iter().flat_map(signed_add_terms).collect(),
        Expr::Mul(fs) => {
            if let Some(rest) = peel_unit_negative_factor(fs) {
                vec![(false, canonical_add_term(&rest))]
            } else if fs.len() == 2 {
                if matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative() && n.abs() == BigInt::one()) {
                    vec![(false, canonical_add_term(&fs[1]))]
                } else if matches!(fs[1].as_ref(), Expr::Int(n) if n.is_negative() && n.abs() == BigInt::one()) {
                    vec![(false, canonical_add_term(&fs[0]))]
                } else {
                    vec![(true, canonical_add_term(e))]
                }
            } else {
                vec![(true, canonical_add_term(e))]
            }
        }
        _ => vec![(true, canonical_add_term(e))],
    }
}

fn peel_unit_negative_factor(fs: &[ExprArc]) -> Option<ExprArc> {
    let pos = fs
        .iter()
        .position(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative() && n.abs() == BigInt::one()))?;
    let rest: Vec<ExprArc> = fs
        .iter()
        .enumerate()
        .filter_map(|(i, f)| if i == pos { None } else { Some(Arc::clone(f)) })
        .collect();
    Some(match rest.len() {
        0 => Expr::int(1),
        1 => rest.into_iter().next().unwrap(),
        _ => Expr::mul(rest),
    })
}

fn rebuild_mul(fs: Vec<ExprArc>) -> ExprArc {
    match fs.len() {
        0 => Expr::int(1),
        1 => fs.into_iter().next().unwrap(),
        _ => Expr::mul(fs),
    }
}

fn distribute_mul_over_add(e: &ExprArc) -> ExprArc {
    match e.as_ref() {
        Expr::Mul(fs) => {
            let expanded: Vec<ExprArc> = fs.iter().map(distribute_mul_over_add).collect();
            if let Some(i) = expanded
                .iter()
                .position(|f| matches!(f.as_ref(), Expr::Add(_)))
            {
                let Expr::Add(ts) = expanded[i].as_ref() else {
                    unreachable!();
                };
                let ts: Vec<ExprArc> = ts.iter().cloned().collect();
                let mut factors = expanded;
                factors.remove(i);
                let flat: Vec<ExprArc> = factors
                    .into_iter()
                    .flat_map(|f| flatten_mul_factors(&f))
                    .collect();
                return Expr::add(
                    ts.into_iter()
                        .map(|t| {
                            let mut parts = flat.clone();
                            parts.push(t);
                            rebuild_mul(parts)
                        })
                        .collect(),
                );
            }
            Expr::mul(expanded)
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(distribute_mul_over_add).collect()),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            distribute_mul_over_add(n),
            distribute_mul_over_add(d),
        )),
        Expr::Pow(b, exp) => Expr::pow(distribute_mul_over_add(b), distribute_mul_over_add(exp)),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter().map(distribute_mul_over_add).collect(),
        ),
        _ => Arc::clone(e),
    }
}

fn simplify_balanced_frac(f: &ExprArc) -> ExprArc {
    match f.as_ref() {
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            simplify_add_sum(&distribute_mul_over_add(n)),
            Arc::clone(d),
        )),
        _ => Arc::clone(f),
    }
}

fn rebuild_signed_add(mut terms: Vec<(bool, ExprArc)>) -> ExprArc {
    terms.retain(|(_, t)| !matches!(t.as_ref(), Expr::Int(n) if n.is_zero()));
    if terms.is_empty() {
        return Expr::int(0);
    }
    let parts: Vec<ExprArc> = terms
        .into_iter()
        .map(|(pos, t)| {
            if pos {
                t
            } else {
                Expr::mul(vec![Expr::int(-1), t])
            }
        })
        .collect();
    match parts.len() {
        0 => Expr::int(0),
        1 => parts.into_iter().next().unwrap(),
        _ => Expr::add(parts),
    }
}

fn exp_signed_pair(a: &ExprArc, b: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    let (n0, a0) = unwrap_exp_arg(a)?;
    let (n1, a1) = unwrap_exp_arg(b)?;
    if !n0 && n1 {
        Some((a0, a1))
    } else if n0 && !n1 {
        Some((a1, a0))
    } else {
        None
    }
}

/// `exp(f) - scale^{-1} → scale^{-1} * (exp(f + ln(scale)) - 1)` (w-layer / MRV).
pub(crate) fn rewrite_exp_minus_scale_inv(
    f: &ExprArc,
    scale_inv: ExprArc,
    ln_scale: ExprArc,
) -> ExprArc {
    let shifted = Expr::add(vec![Arc::clone(f), ln_scale]);
    exp_scale_times_exp_minus_one(scale_inv, shifted)
}

pub(crate) fn simplify_add_sum(e: &ExprArc) -> ExprArc {
    let mut terms: Vec<(bool, ExprArc)> = signed_add_terms(e)
        .into_iter()
        .map(|(p, t)| (p, canonical_add_term(&t)))
        .collect();
    let mut i = 0;
    while i < terms.len() {
        let (pos_i, ref t_i) = terms[i];
        let mut cancelled = false;
        let mut j = i + 1;
        while j < terms.len() {
            if &terms[j].1 == t_i && terms[j].0 != pos_i {
                terms.remove(j);
                terms.remove(i);
                cancelled = true;
                break;
            }
            j += 1;
        }
        if !cancelled {
            i += 1;
        }
    }
    rebuild_signed_add(terms)
}

/// Cancel opposing addends inside `exp` arguments (e.g. `x+1/x-x` → `1/x`).
pub(crate) fn simplify_exp_argument_adds(expr: &ExprArc) -> ExprArc {
    match expr.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => Expr::func(
            FuncKind::Exp,
            vec![simplify_add_sum(&simplify_exp_argument_adds(&args[0]))],
        ),
        Expr::Add(ts) => Expr::add(
            ts.iter()
                .map(simplify_exp_argument_adds)
                .collect(),
        ),
        Expr::Mul(fs) => Expr::mul(
            fs.iter()
                .map(simplify_exp_argument_adds)
                .collect(),
        ),
        Expr::Pow(b, e) => Expr::pow(
            simplify_exp_argument_adds(b),
            simplify_exp_argument_adds(e),
        ),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            simplify_exp_argument_adds(n),
            simplify_exp_argument_adds(d),
        )),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter()
                .map(simplify_exp_argument_adds)
                .collect(),
        ),
        _ => Arc::clone(expr),
    }
}

/// `exp(N/D ± x)` → `exp((N ∓ x·D)/D)` inside `exp` args (CK-INT-60 ratio balance).
pub(crate) fn balance_exp_arguments_frac_var(expr: &ExprArc, var: &Ident) -> ExprArc {
    match expr.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            let arg = balance_frac_minus_var(&args[0], var);
            Expr::func(FuncKind::Exp, vec![arg])
        }
        Expr::Add(ts) => Expr::add(
            ts.iter()
                .map(|t| balance_exp_arguments_frac_var(t, var))
                .collect(),
        ),
        Expr::Mul(fs) => Expr::mul(
            fs.iter()
                .map(|t| balance_exp_arguments_frac_var(t, var))
                .collect(),
        ),
        Expr::Pow(b, e) => Expr::pow(
            balance_exp_arguments_frac_var(b, var),
            balance_exp_arguments_frac_var(e, var),
        ),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            balance_exp_arguments_frac_var(n, var),
            balance_exp_arguments_frac_var(d, var),
        )),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter()
                .map(|a| balance_exp_arguments_frac_var(a, var))
                .collect(),
        ),
        _ => Arc::clone(expr),
    }
}

fn balance_frac_minus_var(f: &ExprArc, var: &Ident) -> ExprArc {
    if let Some(balanced) = try_balance_frac_minus_var(f, var) {
        return balanced;
    }
    match f.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| balance_frac_minus_var(t, var)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|t| balance_frac_minus_var(t, var)).collect()),
        Expr::Pow(b, e) => Expr::pow(balance_frac_minus_var(b, var), balance_frac_minus_var(e, var)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            balance_frac_minus_var(n, var),
            balance_frac_minus_var(d, var),
        )),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter().map(|a| balance_frac_minus_var(a, var)).collect(),
        ),
        _ => Arc::clone(f),
    }
}

fn try_balance_frac_minus_var(f: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let Expr::Add(ts) = f.as_ref() else {
        return None;
    };
    let mut frac = None;
    let mut var_term = None;
    for t in ts {
        if matches!(t.as_ref(), Expr::Frac(_, _)) {
            if frac.is_some() {
                return None;
            }
            frac = Some(Arc::clone(t));
        } else if linear_coeff_of_var(t, var).is_some() {
            if var_term.is_some() {
                return None;
            }
            var_term = Some(Arc::clone(t));
        } else if depends_on_var(t, var) {
            return None;
        }
    }
    let frac = frac?;
    let var_term = var_term?;
    let Expr::Frac(n, d) = frac.as_ref() else {
        return None;
    };
    // `N/D + t` with linear `t` → `(N + t·D)/D` (e.g. `inner - x` for CK-INT-60 ratio).
    Some(simplify_balanced_frac(&Arc::new(Expr::Frac(
        Expr::add(vec![
            Arc::clone(n),
            Expr::mul(vec![var_term, Arc::clone(d)]),
        ]),
        Arc::clone(d),
    ))))
}

fn linear_coeff_of_var(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
    if is_var(e, var) {
        return Some(Expr::int(1));
    }
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            if is_var(&fs[0], var) && matches!(fs[1].as_ref(), Expr::Int(n) if n.is_one()) {
                return Some(Expr::int(1));
            }
            if is_var(&fs[1], var) && matches!(fs[0].as_ref(), Expr::Int(n) if n.is_one()) {
                return Some(Expr::int(1));
            }
            if is_var(&fs[0], var) {
                return Some(Arc::clone(&fs[1]));
            }
            if is_var(&fs[1], var) {
                return Some(Arc::clone(&fs[0]));
            }
        }
        _ => {}
    }
    None
}

/// Classify `exp(f)` at `+∞`: `0`, `1`, or `None`.
pub(crate) fn classify_exp_at_plus_infinity(
    f: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let f = simplify_add_sum(f);
    if let Some((n, d)) = unwrap_signed_frac(&f) {
        if super::mrv::vanishes_faster_than_at_plus_infinity(&n, &d, var, ctx) {
            return Some(Expr::int(1));
        }
    }
    if vanishes_at_plus_infinity(&f, var) {
        return Some(Expr::int(1));
    }
    if exp_argument_tends_to_negative_infinity(&f, var) {
        return Some(Expr::int(0));
    }
    None
}

/// `±num/den` or `num/den` for vanishing-ratio checks after preprocess.
pub(crate) fn unwrap_signed_frac(e: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    if let Expr::Frac(n, d) = e.as_ref() {
        return Some((Arc::clone(n), Arc::clone(d)));
    }
    if let Expr::Mul(fs) = e.as_ref() {
        if let Some(i) = fs.iter().position(|f| matches!(f.as_ref(), Expr::Frac(_, _))) {
            let Expr::Frac(n, d) = fs[i].as_ref() else {
                unreachable!();
            };
            let mut num = Arc::clone(n);
            let mut scalar = Expr::int(1);
            for (j, f) in fs.iter().enumerate() {
                if j == i {
                    continue;
                }
                scalar = Expr::mul(vec![scalar, Arc::clone(f)]);
            }
            if !matches!(scalar.as_ref(), Expr::Int(n) if n.is_one()) {
                num = Expr::mul(vec![scalar, num]);
            }
            return Some((num, Arc::clone(d)));
        }
    }
    None
}

/// Signed `exp(f)` at `+∞`.
pub(crate) fn classify_signed_exp_at_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let (neg, inner) = unwrap_exp_arg(expr)?;
    classify_exp_at_plus_infinity(&inner, var, ctx)
        .map(|v| if neg { Expr::mul(vec![Expr::int(-1), v]) } else { v })
}

fn exp_argument_tends_to_negative_infinity(f: &ExprArc, var: &Ident) -> bool {
    let f = simplify_add_sum(f);
    if let Some(c) = super::mrv::linear_coeff_in_var(&f, var) {
        if super::mrv::try_const_f64(&c).is_some_and(|x| x < 0.0) {
            return true;
        }
        if super::mrv::try_const_f64(&c).is_some_and(|x| x > 0.0) {
            return false;
        }
    }
    signed_add_terms(&f).iter().any(|(pos, t)| {
        !pos && (is_positive_var_power_ge2(t, var) || is_neg_exp_of_positive_growth(t, var))
    })
}

fn is_positive_var_power_ge2(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Pow(b, exp) if is_var(b, var) => {
            matches!(exp.as_ref(), Expr::Int(n) if *n >= BigInt::from(2))
        }
        Expr::Mul(fs) if fs.len() == 2 => {
            let neg = fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()));
            neg && fs.iter().any(|f| {
                matches!(
                    f.as_ref(),
                    Expr::Pow(b, exp)
                        if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if *n >= BigInt::from(2))
                )
            })
        }
        _ => false,
    }
}

fn is_neg_exp_of_positive_growth(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            is_var(&args[0], var)
                || super::mrv::linear_coeff_in_var(&args[0], var)
                    .and_then(|c| super::mrv::try_const_f64(&c))
                    .is_some_and(|x| x > 0.0)
        }
        Expr::Mul(fs) if fs.len() == 2 => {
            let neg = fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()));
            let growth = fs.iter().any(|f| exp_arg_grows_at_plus_infinity(f, var));
            neg && growth
        }
        _ => false,
    }
}

fn exp_arg_grows_at_plus_infinity(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            super::mrv::linear_coeff_in_var(&args[0], var)
                .and_then(|c| super::mrv::try_const_f64(&c))
                .is_some_and(|x| x > 0.0)
                || is_var(&args[0], var)
        }
        _ => is_var(e, var)
            || super::mrv::linear_coeff_in_var(e, var)
                .and_then(|c| super::mrv::try_const_f64(&c))
                .is_some_and(|x| x > 0.0),
    }
}

/// `exp(L) * ε` with `ε = c*exp(g)` → `c*exp(L+g)`; `exp(L)*Add(…)` → `Add(exp(L)*…)`.
pub(crate) fn algebraize_exp_vanishing_products(expr: &ExprArc, var: &Ident) -> ExprArc {
    let rebuilt = match expr.as_ref() {
        Expr::Add(ts) => Expr::add(
            ts.iter()
                .map(|t| algebraize_exp_vanishing_products(t, var))
                .collect(),
        ),
        Expr::Mul(fs) => {
            let folded: Vec<ExprArc> = fs
                .iter()
                .map(|t| algebraize_exp_vanishing_products(t, var))
                .collect();
            algebraize_exp_mul_factors(&folded)
        }
        Expr::Pow(b, e) => Expr::pow(
            algebraize_exp_vanishing_products(b, var),
            algebraize_exp_vanishing_products(e, var),
        ),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            algebraize_exp_vanishing_products(n, var),
            algebraize_exp_vanishing_products(d, var),
        )),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter()
                .map(|a| algebraize_exp_vanishing_products(a, var))
                .collect(),
        ),
        _ => Arc::clone(expr),
    };
    rebuilt
}

fn algebraize_exp_mul_factors(factors: &[ExprArc]) -> ExprArc {
    let flat = flatten_mul_factors_slice(factors);
    let product = rebuild_mul(flat.clone());
    if detect_exp_difference(&product).is_some() {
        return product;
    }
    let mut logs = Vec::new();
    let mut coeff = BigInt::one();
    let mut add_factor = None;
    let mut other = Vec::new();

    for f in &flat {
        if let Some(log) = exp_func_arg(f) {
            logs.push(log);
        } else if let Expr::Int(n) = f.as_ref() {
            coeff *= n;
        } else if let Some((neg, g)) = unwrap_exp_arg(f) {
            if neg {
                coeff = -coeff;
            }
            logs.push(g);
        } else if let Expr::Add(_) = f.as_ref() {
            if add_factor.is_none() && logs.len() == 1 {
                add_factor = Some(Arc::clone(f));
            } else {
                other.push(Arc::clone(f));
            }
        } else {
            other.push(Arc::clone(f));
        }
    }

    if logs.is_empty() || !other.is_empty() {
        return Expr::mul(flat);
    }

    let log_sum = logs
        .into_iter()
        .reduce(|a, b| simplify_add_sum(&Expr::add(vec![a, b])))
        .unwrap();

    if let Some(add) = add_factor {
        let Expr::Add(ts) = add.as_ref() else {
            return Expr::mul(flat);
        };
        return Expr::add(ts.iter().map(|t| algebraize_exp_log_pair(&log_sum, t)).collect());
    }

    let body = Expr::func(FuncKind::Exp, vec![log_sum]);
    if coeff.is_zero() {
        Expr::int(0)
    } else if coeff == BigInt::one() {
        body
    } else if coeff == -BigInt::one() {
        Expr::mul(vec![Expr::int(-1), body])
    } else {
        Expr::mul(vec![Arc::new(Expr::Int(coeff)), body])
    }
}

fn algebraize_exp_log_pair(log: &ExprArc, factor: &ExprArc) -> ExprArc {
    if let Some((neg, g)) = unwrap_exp_arg(factor) {
        let merged = simplify_add_sum(&Expr::add(vec![Arc::clone(log), g]));
        let body = Expr::func(FuncKind::Exp, vec![merged]);
        if neg {
            Expr::mul(vec![Expr::int(-1), body])
        } else {
            body
        }
    } else {
        Expr::mul(vec![
            Expr::func(FuncKind::Exp, vec![Arc::clone(log)]),
            Arc::clone(factor),
        ])
    }
}

fn flatten_mul_factors_slice(factors: &[ExprArc]) -> Vec<ExprArc> {
    let mut out = Vec::new();
    for f in factors {
        match f.as_ref() {
            Expr::Mul(inner) => out.extend(flatten_mul_factors_slice(inner)),
            _ => out.push(Arc::clone(f)),
        }
    }
    out
}

/// Gruntz first-order (P1): `exp(ε)-1 ~ ε` when `ε→0` at `+∞`.
pub(crate) fn first_order_exp_vanishing_epsilon(expr: &ExprArc, var: &Ident) -> ExprArc {
    if let Some(form) = detect_exp_difference(expr) {
        let epsilon = balance_epsilon_expr(&form.epsilon, var);
        if epsilon_vanishes_at_plus_infinity(&epsilon, var) {
            let mut factors = vec![
                Expr::func(FuncKind::Exp, vec![form.scale_log]),
                epsilon,
            ];
            if let Some(c) = form.complement {
                factors.push(c);
            }
            return Expr::mul(factors);
        }
    }
    match expr.as_ref() {
        Expr::Add(ts) => Expr::add(
            ts.iter()
                .map(|t| first_order_exp_vanishing_epsilon(t, var))
                .collect(),
        ),
        Expr::Mul(fs) => Expr::mul(
            fs.iter()
                .map(|t| first_order_exp_vanishing_epsilon(t, var))
                .collect(),
        ),
        Expr::Pow(b, e) => Expr::pow(
            first_order_exp_vanishing_epsilon(b, var),
            first_order_exp_vanishing_epsilon(e, var),
        ),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            first_order_exp_vanishing_epsilon(n, var),
            first_order_exp_vanishing_epsilon(d, var),
        )),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter()
                .map(|a| first_order_exp_vanishing_epsilon(a, var))
                .collect(),
        ),
        _ => Arc::clone(expr),
    }
}

pub(crate) fn match_exp_times_exp_minus_one(expr: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    detect_exp_difference(expr).map(|form| (form.scale_log, form.epsilon))
}

pub(crate) fn unwrap_exp_arg(e: &ExprArc) -> Option<(bool, ExprArc)> {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => Some((false, Arc::clone(&args[0]))),
        Expr::Mul(fs) if fs.len() == 2 => {
            let (neg, inner) = if matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative()) {
                (true, &fs[1])
            } else if matches!(fs[1].as_ref(), Expr::Int(n) if n.is_negative()) {
                (true, &fs[0])
            } else {
                return None;
            };
            let Expr::Func(FuncKind::Exp, args) = inner.as_ref() else {
                return None;
            };
            if args.len() == 1 {
                Some((neg, Arc::clone(&args[0])))
            } else {
                None
            }
        }
        _ => None,
    }
}

fn exp_func_arg(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => Some(Arc::clone(&args[0])),
        _ => None,
    }
}

pub(crate) fn exp_minus_one_epsilon(e: &ExprArc) -> Option<ExprArc> {
    exp_minus_one_inner_arg(e)
}

/// `exp(ε) - 1` factor (not mergeable with sibling `exp` terms).
pub(crate) fn is_exp_minus_one_factor(e: &ExprArc) -> bool {
    exp_minus_one_inner_arg(e).is_some()
}

fn exp_minus_one_inner_arg(e: &ExprArc) -> Option<ExprArc> {
    let Expr::Add(ts) = e.as_ref() else {
        return None;
    };
    if ts.len() != 2 {
        return None;
    }
    let (pos, neg) = if matches!(ts[1].as_ref(), Expr::Int(n) if n.is_negative()) {
        (&ts[0], &ts[1])
    } else if matches!(ts[0].as_ref(), Expr::Int(n) if n.is_negative()) {
        (&ts[1], &ts[0])
    } else {
        return None;
    };
    if !matches!(neg.as_ref(), Expr::Int(n) if n.is_negative()) {
        return None;
    }
    exp_func_arg(pos)
}

fn remove_add_term(sum: &ExprArc, term: &ExprArc) -> Option<ExprArc> {
    if sum == term {
        return Some(Expr::int(0));
    }
    let Expr::Add(ts) = sum.as_ref() else {
        return None;
    };
    if !ts.iter().any(|t| t == term) {
        return None;
    }
    let rest: Vec<ExprArc> = ts.iter().filter(|t| *t != term).cloned().collect();
    Some(match rest.len() {
        0 => Expr::int(0),
        1 => Arc::clone(&rest[0]),
        _ => Expr::add(rest),
    })
}

pub(crate) fn vanishes_at_plus_infinity(e: &ExprArc, var: &Ident) -> bool {
    epsilon_vanishes_at_plus_infinity(e, var)
}

/// True when `e → 0` at `+∞` (for Gruntz ε in `exp(ε)-1 ~ ε`).
pub(crate) fn epsilon_vanishes_at_plus_infinity(e: &ExprArc, var: &Ident) -> bool {
    let balanced = balance_epsilon_expr(e, var);
    match balanced.as_ref() {
        Expr::Pow(b, exp) if is_var(b, var) => {
            matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
        }
        Expr::Frac(n, d) => {
            (matches!(n.as_ref(), Expr::Int(nn) if nn.is_one()) && is_var(d, var))
                || super::mrv::ratio_tends_to_zero_at_plus_infinity(n, d, var)
        }
        Expr::Add(ts) => ts
            .iter()
            .all(|t| epsilon_vanishes_at_plus_infinity(t, var) || !depends_on_var(t, var)),
        Expr::Mul(fs) => fs
            .iter()
            .all(|f| epsilon_vanishes_at_plus_infinity(f, var) || !depends_on_var(f, var)),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            exp_inner_vanishes_at_plus_infinity(&args[0], var)
                || exp_of_neg_exp_growth(&balanced, var)
        }
        _ => {
            if let Some((n, d)) = unwrap_signed_frac(&balanced) {
                return super::mrv::ratio_tends_to_zero_at_plus_infinity(&n, &d, var);
            }
            false
        }
    }
}

fn balance_epsilon_expr(e: &ExprArc, var: &Ident) -> ExprArc {
    try_balance_frac_minus_var(e, var).unwrap_or_else(|| Arc::clone(e))
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn exp_inner_vanishes_at_plus_infinity(arg: &ExprArc, var: &Ident) -> bool {
    if super::mrv::linear_coeff_in_var(arg, var)
        .and_then(|c| super::mrv::try_const_f64(&c))
        .is_some_and(|x| x < 0.0)
    {
        return true;
    }
    if is_neg_var_exp(arg, var) {
        return true;
    }
    if is_negated_var_power_mul(arg, var) {
        return true;
    }
    match arg.as_ref() {
        Expr::Mul(fs) if fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative())) => {
            fs.iter().any(|f| is_negated_var_power_mul(f, var) || {
                matches!(
                    f.as_ref(),
                    Expr::Func(FuncKind::Exp, a)
                        if a.len() == 1
                            && super::mrv::linear_coeff_in_var(&a[0], var)
                                .and_then(|c| super::mrv::try_const_f64(&c))
                                .is_some_and(|x| x > 0.0)
                )
            })
        }
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            is_neg_var_exp(&args[0], var)
                || super::mrv::linear_coeff_in_var(&args[0], var)
                    .and_then(|c| super::mrv::try_const_f64(&c))
                    .is_some_and(|x| x < 0.0)
                || exp_inner_vanishes_at_plus_infinity(&args[0], var)
        }
        _ => false,
    }
}

fn exp_of_neg_exp_growth(e: &ExprArc, var: &Ident) -> bool {
    fn arg_grows(f: &ExprArc, var: &Ident) -> bool {
        match f.as_ref() {
            Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
                super::mrv::linear_coeff_in_var(&args[0], var)
                    .and_then(|c| super::mrv::try_const_f64(&c))
                    .is_some_and(|x| x > 0.0)
                    || is_var(&args[0], var)
            }
            _ => is_var(f, var)
                || super::mrv::linear_coeff_in_var(f, var)
                    .and_then(|c| super::mrv::try_const_f64(&c))
                    .is_some_and(|x| x > 0.0),
        }
    }
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => match args[0].as_ref() {
            Expr::Mul(fs) if fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative())) => {
                fs.iter().any(|f| arg_grows(f, var))
            }
            _ => false,
        },
        _ => false,
    }
}

/// `-x^n` as `Mul(-1, Pow(x,n))` inside an `exp` argument.
fn is_negated_var_power_mul(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            let neg = fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()));
            neg && fs.iter().any(|f| {
                matches!(
                    f.as_ref(),
                    Expr::Pow(b, exp)
                        if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n.is_positive())
                )
            })
        }
        _ => false,
    }
}

fn is_neg_var_exp(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Mul(fs) if fs.len() == 2
            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
            && fs.iter().any(|f| is_var(f, var))
    ) || matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp, args)
            if args.len() == 1
                && matches!(
                    args[0].as_ref(),
                    Expr::Mul(fs) if fs.len() == 2
                        && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
                        && fs.iter().any(|f| is_var(f, var))
                )
    )
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, FuncKind};
    use std::sync::Arc;

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn fold_mul_shifted_difference_gruntz() {
        let ctx = xcas_default();
        let stmts = giac_parse::parse_program(
            "exp(x)*(exp(1/x-exp(-x))-exp(1/x));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let r = fold_exp_shifted_difference(e);
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("exp(-exp(-x))") || s.contains("exp(-exp(-1*x))"),
            "got {s}"
        );
    }

    #[test]
    fn rewrite_exp_minus_w_inv_matches_remove_lnexp() {
        use super::super::mrv_w::{mrv_ln_w_expr, mrv_w_expr};
        let w = mrv_w_expr();
        let inner = Expr::add(vec![
            Expr::mul(vec![Expr::int(-1), mrv_ln_w_expr()]),
            Expr::sym("eps"),
        ]);
        let r = exp_scale_times_exp_minus_one(Expr::pow(w, Expr::int(-1)), inner);
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("exp(eps)") || s.contains("exp(-ln(_mrv_w)+eps)"),
            "expected w^-1*(exp(eps)-1) style, got {s}"
        );
        assert!(s.contains("_mrv_w"), "got {s}");
    }

    #[test]
    fn fold_add_exp_difference() {
        let a = Expr::func(FuncKind::Exp, vec![Expr::add(vec![Expr::sym("x"), Expr::int(1)])]);
        let b = Expr::mul(vec![
            Expr::int(-1),
            Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
        ]);
        let e = Expr::add(vec![a, b]);
        let r = fold_exp_shifted_difference(&e);
        let s = format_expr(r.as_ref());
        assert!(s.contains("exp(x)") && s.contains("-1"), "got {s}");
    }

    #[test]
    fn balance_frac_minus_var_rewrites_inner_minus_x() {
        let var = Ident::new("x");
        let inner = Arc::new(Expr::Frac(
            Expr::sym("a"),
            Expr::sym("b"),
        ));
        let e = Expr::add(vec![
            Arc::clone(&inner),
            Expr::mul(vec![Expr::int(-1), Expr::sym("x")]),
        ]);
        let r = balance_exp_arguments_frac_var(
            &Expr::func(FuncKind::Exp, vec![e]),
            &var,
        );
        let s = format_expr(r.as_ref());
        assert!(
            s.contains("a") && s.contains("b") && s.contains("x"),
            "expected (a-x*b)/b inside exp, got {s}"
        );
        assert!(!s.contains("a+-1*b") && !s.contains("a+(-1)*b"), "wrong coeff-only balance: {s}");
    }

    #[test]
    fn first_order_vanishing_epsilon_rewrite() {
        let var = Ident::new("x");
        let e = Expr::mul(vec![
            Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
            Expr::add(vec![
                Expr::func(FuncKind::Exp, vec![Expr::pow(Expr::sym("x"), Expr::int(-1))]),
                Expr::int(-1),
            ]),
        ]);
        let r = first_order_exp_vanishing_epsilon(&e, &var);
        let s = format_expr(r.as_ref());
        assert!(
            !match_exp_times_exp_minus_one(&r).is_some(),
            "should not keep exp()-1: {s}"
        );
        assert!(
            s.contains("x^-1") || s.contains("1*x^-1") || s.contains("1/x"),
            "expected exp(x)*x^-1 style, got {s}"
        );
    }

    #[test]
    fn limit_preprocessed_gruntz_minus_one() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let stmts = giac_parse::parse_program(
            "exp(x)*(exp(1/x-exp(-x))-exp(1/x));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let pre = crate::limit_engine::preprocess::limit_preprocess_struct(e, &var);
        let r = super::super::asymptotic::limit_at_plus_infinity(&pre, &var, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "-1");
    }

    #[test]
    fn limit_exp_over_exp_ck_int_ratio() {
        use super::super::ck_int_gruntz_fixture::ratio;
        let ctx = xcas_default();
        let var = Ident::new("x");
        let r = super::super::asymptotic::limit_at_plus_infinity(&ratio(), &var, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1");
    }

    #[test]
    fn debug_ck61_const() {
        use crate::limit_engine::ck_int_gruntz_fixture::ck_int_61;
        use crate::limit_engine::preprocess::limit_preprocess_struct;
        use crate::plugin::xcas_default;
        let ctx = xcas_default();
        let var = Ident::new("x");
        let pre = limit_preprocess_struct(&ck_int_61(), &var);
        let r = super::super::asymptotic::limit_at_plus_infinity(&pre, &var, &ctx);
        eprintln!("limit: {r:?} -> {}", r.as_ref().map(|x| format_expr(x.as_ref())).unwrap_or_default());
    }

    #[test]
    fn fold_ck_int_61_shape() {
        use crate::limit_engine::ck_int_gruntz_fixture::{ck_int_61, exp_inner};
        let var = Ident::new("x");
        let e = ck_int_61();
        let add = Expr::add(vec![
            exp_inner(),
            Expr::mul(vec![Expr::int(-1), Expr::func(FuncKind::Exp, vec![Expr::sym("x")])]),
        ]);
        let add_fold = fold_exp_shifted_difference(&add);
        assert!(
            match_exp_times_exp_minus_one(&add_fold).is_some(),
            "add fold: {}",
            format_expr(add_fold.as_ref())
        );
        let pre = crate::limit_engine::preprocess::limit_preprocess_struct(&e, &var);
        assert!(
            format_expr(pre.as_ref()).contains("exp(x)"),
            "struct preprocess: {}",
            format_expr(pre.as_ref())
        );
    }

    #[test]
    fn fold_nested_gruntz_add_difference() {
        let ctx = xcas_default();
        let stmts = giac_parse::parse_program(
            "exp(1/x+exp(-x)+exp(-(x^2)))-exp(1/x-exp(-exp(x)));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let giac_core::Expr::Add(ts) = e.as_ref() else {
            panic!("expected add, got {}", format_expr(e.as_ref()));
        };
        assert_eq!(ts.len(), 2);
        let form = detect_exp_difference_add(&ts[0], &ts[1]).expect("nested gruntz exp diff");
        let s = format_expr(&emit_exp_difference(form));
        assert!(s.contains("-1") && s.contains("exp("), "got {s}");
    }

    #[test]
    fn fold_nested_gruntz_full_expr() {
        let ctx = xcas_default();
        let stmts = giac_parse::parse_program(
            "exp(x)*(exp(1/x+exp(-x)+exp(-(x^2)))-exp(1/x-exp(-exp(x))));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let folded = fold_exp_shifted_difference(e);
        let s = format_expr(folded.as_ref());
        assert!(
            match_exp_times_exp_minus_one(&folded).is_some(),
            "expected exp(L)*(exp(S)-1) form, got {s}"
        );
        let _ = ctx;
    }

    #[test]
    fn parse_neg_x_squared_in_exp() {
        use giac_core::format_expr;
        let ctx = xcas_default();
        let stmts = giac_parse::parse_program("exp(-(x^2));", &ctx).unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let s = format_expr(e.as_ref());
        assert!(
            s.contains("-x^2") || s.contains("-1*x^2"),
            "expected -(x^2) in exp arg, got {s}"
        );
    }

    #[test]
    fn limit_nested_gruntz_preprocessed() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let stmts = giac_parse::parse_program(
            "exp(x)*(exp(1/x+exp(-x)+exp(-(x^2)))-exp(1/x-exp(-exp(x))));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let pre = crate::limit_engine::preprocess::limit_preprocess_struct(e, &var);
        let r = super::super::asymptotic::limit_at_plus_infinity(&pre, &var, &ctx)
            .expect("limit_at_plus_infinity");
        assert_eq!(format_expr(r.as_ref()), "1");
    }

    #[test]
    fn exp_scale_times_exp_minus_one_shape() {
        let r = exp_scale_times_exp_minus_one(Expr::sym("s"), Expr::sym("eps"));
        let s = format_expr(r.as_ref());
        assert!(s.contains("exp(eps)") && s.contains("-1"), "got {s}");
    }
}

