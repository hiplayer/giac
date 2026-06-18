//! Shared `exp` difference normalization (GIAC-limit-exp-diff P1).
//!
//! Upstream `remove_lnexp` (`series.cc`) applies `ln_expand` / `exp_series` on coefficients;
//! the identity `exp(f) - w^{-1} = w^{-1}(exp(f + ln(w)) - 1)` is the MRV (`w`) layer.
//!
//! `fold_exp_shifted_difference` is the same ε→0 pattern on the original variable:
//! all shapes reduce to `exp(scale_log) * (exp(ε) - 1)` via [`detect_exp_difference`] + [`emit_exp_difference`].

use std::sync::Arc;

use giac_core::{eval, Context, Expr, ExprArc, FuncKind, Ident};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use crate::expr_util::depends_on_var;

/// Normal form: `exp(scale_log) * (exp(epsilon) - 1)`.
#[derive(Clone, Debug)]
struct ExpDiffForm {
    scale_log: ExprArc,
    epsilon: ExprArc,
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
    exp_scale_times_exp_minus_one(
        Expr::func(FuncKind::Exp, vec![form.scale_log]),
        form.epsilon,
    )
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
    for (i, f) in expanded.iter().enumerate() {
        if i == mi {
            continue;
        }
        scale_logs.push(exp_func_arg(f)?);
    }
    Some(ExpDiffForm {
        scale_log: sum_exp_logs(scale_logs),
        epsilon,
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

/// `exp(A) - exp(B) → exp(B) * (exp(ε) - 1)` when `ε = A-B` cancels a shared add term.
/// General `exp(A)-exp(B)` without shared subterms stays for MRV / `remove_lnexp` (CK-INT-61).
fn detect_exp_difference_add(a: &ExprArc, b: &ExprArc) -> Option<ExpDiffForm> {
    let (pos_arg, neg_arg) = exp_signed_pair(a, b)?;
    let epsilon = remove_add_term(&pos_arg, &neg_arg)?;
    Some(ExpDiffForm {
        scale_log: neg_arg,
        epsilon,
    })
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

/// `exp(a)/exp(b)` @ `+∞`: `(exp(a)-exp(b))/v * v/exp(b)` when both sub-limits exist.
pub(crate) fn try_limit_exp_over_exp_via_quotient(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let (a, b) = match_exp_over_exp(expr)?;
    let v = Expr::sym(var.as_str());
    let diff_q = Expr::mul(vec![
        Expr::add(vec![
            Expr::func(FuncKind::Exp, vec![Arc::clone(&a)]),
            Expr::mul(vec![
                Expr::int(-1),
                Expr::func(FuncKind::Exp, vec![Arc::clone(&b)]),
            ]),
        ]),
        Expr::pow(v.clone(), Expr::int(-1)),
    ]);
    let scale = Expr::mul(vec![
        v,
        Expr::pow(Expr::func(FuncKind::Exp, vec![b.clone()]), Expr::int(-1)),
    ]);
    if try_limit_var_over_exp(&b, var, ctx).is_some_and(|s| is_zero_expr(&s)) {
        return Some(Expr::int(1));
    }
    let l = super::asymptotic::limit_at_plus_infinity(&diff_q, var, ctx).ok()?;
    let s = try_limit_var_over_exp(&b, var, ctx)
        .or_else(|| super::asymptotic::limit_at_plus_infinity(&scale, var, ctx).ok())?;
    if depends_on_var(&l, var) || depends_on_var(&s, var) {
        return None;
    }
    if is_zero_expr(&s) {
        return Some(Expr::int(1));
    }
    eval(
        Expr::add(vec![Expr::int(1), Expr::mul(vec![l, s])]).as_ref(),
        ctx,
    )
    .ok()
}

/// `x/exp(b) → 0` at `+∞` when `b` is `x` or grows faster than `x` (e.g. `exp(x)`).
fn try_limit_var_over_exp(b: &ExprArc, var: &Ident, ctx: &Context) -> Option<ExprArc> {
    if is_var(b, var) {
        return Some(Expr::int(0));
    }
    if super::mrv::linear_coeff_in_var(b, var)
        .and_then(|c| super::mrv::try_const_f64(&c))
        .is_some_and(|x| x > 0.0)
    {
        return Some(Expr::int(0));
    }
    let v = Expr::sym(var.as_str());
    let scale = Expr::mul(vec![
        v,
        Expr::pow(Expr::func(FuncKind::Exp, vec![Arc::clone(b)]), Expr::int(-1)),
    ]);
    if vanishes_at_plus_infinity(&scale, var) {
        return Some(Expr::int(0));
    }
    let _ = ctx;
    None
}

fn is_zero_expr(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

fn match_exp_over_exp(expr: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    if let Expr::Frac(num, den) = expr.as_ref() {
        if let (Expr::Func(FuncKind::Exp, na), Expr::Func(FuncKind::Exp, da)) =
            (num.as_ref(), den.as_ref())
        {
            if na.len() == 1 && da.len() == 1 {
                return Some((Arc::clone(&na[0]), Arc::clone(&da[0])));
            }
        }
    }
    let Expr::Mul(fs) = expr.as_ref() else {
        return None;
    };
    if fs.len() != 2 {
        return None;
    }
    for (i, j) in [(0, 1), (1, 0)] {
        if let Some(a) = exp_func_arg(&fs[i]) {
            if let Some(b) = exp_func_arg_from_inv(&fs[j]) {
                return Some((a, b));
            }
        }
    }
    None
}

fn exp_func_arg_from_inv(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Pow(base, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) =>
        {
            exp_func_arg(base)
        }
        _ => None,
    }
}

/// After [`super::preprocess::limit_preprocess_struct`]: `exp(L)*(exp(S)-1)` @ `+∞` (P2).
pub(crate) fn try_limit_exp_times_exp_minus_one_preprocessed(
    pre: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let (outer_arg, small_shift) = match_exp_times_exp_minus_one(pre)?;
    if vanishes_at_plus_infinity(&small_shift, var)
        || super::mrv_lead_term::limit_unidirectional_plus_infinity(&small_shift, var, ctx)
            .ok()
            .is_some_and(|r| is_zero_expr(&r))
    {
        let approx = Expr::mul(vec![
            Expr::func(FuncKind::Exp, vec![Arc::clone(&outer_arg)]),
            Arc::clone(&small_shift),
        ]);
        if let Ok(r) = super::mrv_lead_term::limit_unidirectional_plus_infinity(&approx, var, ctx) {
            if !depends_on_var(&r, var) {
                return eval(r.as_ref(), ctx).ok();
            }
        }
    }
    let rest = exp_rest_after_unit_var(&outer_arg, var)?;
    if !is_neg_exp_of_neg_var(&small_shift, var) {
        return None;
    }
    limit_neg_exp_rest_at_plus_infinity(&rest, var, ctx)
}

/// `limit(exp(f), x, +∞)` when `limit(f, x, +∞) = L` is finite: result `exp(L)`.
pub(crate) fn try_limit_exp_finite_exponent_at_plus_infinity(
    pre: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Option<ExprArc> {
    let arg = match pre.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => &args[0],
        _ => return None,
    };
    let lim_arg = super::mrv_lead_term::limit_unidirectional_plus_infinity(arg, var, ctx).ok()?;
    if depends_on_var(&lim_arg, var) || is_infinity_symbol(&lim_arg) {
        return None;
    }
    eval(
        Expr::func(FuncKind::Exp, vec![lim_arg]).as_ref(),
        ctx,
    )
    .ok()
}

fn is_infinity_symbol(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Symbol(id) if id.as_str() == "+infinity" || id.as_str() == "infinity" || id.as_str() == "-infinity"
    )
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

fn exp_rest_after_unit_var(arg: &ExprArc, var: &Ident) -> Option<ExprArc> {
    if is_unit_var_term(arg, var) {
        return Some(Expr::int(0));
    }
    let Expr::Add(ts) = arg.as_ref() else {
        return None;
    };
    let mut rest = Vec::new();
    let mut saw_unit_var = false;
    for t in ts {
        if is_unit_var_term(t, var) {
            if saw_unit_var {
                return None;
            }
            saw_unit_var = true;
        } else {
            rest.push(Arc::clone(t));
        }
    }
    if !saw_unit_var {
        return None;
    }
    Some(match rest.len() {
        0 => Expr::int(0),
        1 => Arc::clone(&rest[0]),
        _ => Expr::add(rest),
    })
}

fn is_unit_var_term(e: &ExprArc, var: &Ident) -> bool {
    if is_var(e, var) {
        return true;
    }
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            (is_var(&fs[0], var) && is_one(&fs[1])) || (is_var(&fs[1], var) && is_one(&fs[0]))
        }
        _ => false,
    }
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_one(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n == &BigInt::from(1))
}

fn limit_neg_exp_rest_at_plus_infinity(rest: &ExprArc, var: &Ident, ctx: &Context) -> Option<ExprArc> {
    if !depends_on_var(rest, var) {
        return eval(
            Expr::mul(vec![
                Expr::int(-1),
                Expr::func(FuncKind::Exp, vec![Arc::clone(rest)]),
            ])
            .as_ref(),
            ctx,
        )
        .ok();
    }
    if vanishes_at_plus_infinity(rest, var) {
        return Some(Expr::int(-1));
    }
    None
}

fn vanishes_at_plus_infinity(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Pow(b, exp) if is_var(b, var) => {
            matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
        }
        Expr::Frac(n, d) => {
            matches!(n.as_ref(), Expr::Int(nn) if nn.is_one()) && is_var(d, var)
        }
        Expr::Mul(fs) => fs
            .iter()
            .all(|f| vanishes_at_plus_infinity(f, var) || !depends_on_var(f, var)),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => super::mrv::linear_coeff_in_var(
            &args[0],
            var,
        )
        .and_then(|c| super::mrv::try_const_f64(&c))
        .is_some_and(|x| x < 0.0),
        _ => false,
    }
}

fn is_neg_exp_of_neg_var(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            let neg = fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()));
            let exp_inner = fs.iter().find_map(|f| {
                if let Expr::Func(FuncKind::Exp, args) = f.as_ref() {
                    if args.len() == 1 {
                        return Some(&args[0]);
                    }
                }
                None
            });
            neg && exp_inner.is_some_and(|a| is_neg_var_exp(a, var))
        }
        _ => is_neg_var_exp(e, var),
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
        let r = try_limit_exp_times_exp_minus_one_preprocessed(&pre, &var, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "-1");
    }

    #[test]
    fn limit_exp_over_exp_ck_int_ratio() {
        use super::super::ck_int_gruntz_fixture::ratio;
        let ctx = xcas_default();
        let var = Ident::new("x");
        assert!(
            match_exp_over_exp(&ratio()).is_some(),
            "ratio should match exp/exp"
        );
        let r = try_limit_exp_over_exp_via_quotient(&ratio(), &var, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1");
    }

    #[test]
    fn fold_nested_gruntz_add_difference() {
        let ctx = xcas_default();
        let stmts = giac_parse::parse_program(
            "exp(1/x+exp(-x)+exp(-x^2))-exp(1/x-exp(-exp(x)));",
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
        let form = detect_exp_difference_add(&ts[0], &ts[1]);
        assert!(
            form.is_none(),
            "x-layer fold skips general exp diff without shared subterm: {}",
            format_expr(e.as_ref())
        );
    }

    #[test]
    fn fold_nested_gruntz_full_expr() {
        let ctx = xcas_default();
        let stmts = giac_parse::parse_program(
            "exp(x)*(exp(1/x+exp(-x)+exp(-x^2))-exp(1/x-exp(-exp(x))));",
            &ctx,
        )
        .unwrap();
        let giac_core::Stmt::ExprStmt(e) = stmts.first().unwrap() else {
            panic!();
        };
        let folded = fold_exp_shifted_difference(e);
        let s = format_expr(folded.as_ref());
        assert!(
            s.contains("-1") && s.contains("exp("),
            "expected exp(L)*(exp(S)-1) form, got {s}"
        );
        let _ = ctx;
    }

    #[test]
    fn exp_scale_times_exp_minus_one_shape() {
        let r = exp_scale_times_exp_minus_one(Expr::sym("s"), Expr::sym("eps"));
        let s = format_expr(r.as_ref());
        assert!(s.contains("exp(eps)") && s.contains("-1"), "got {s}");
    }
}
