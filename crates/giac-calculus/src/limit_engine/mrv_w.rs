//! Helpers for MRV auxiliary variable `w` (`_mrv_w`) in asymptotic series.

use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind, Ident};
use num_bigint::BigInt;
use num_traits::{One, Zero};

pub(crate) const MRV_W: &str = "_mrv_w";

pub(crate) fn is_mrv_w_var(id: &Ident) -> bool {
    id.as_str() == MRV_W
}

pub(crate) fn mrv_w_expr() -> ExprArc {
    Expr::sym(MRV_W)
}

pub(crate) fn mrv_ln_w_expr() -> ExprArc {
    Expr::func(FuncKind::Ln, vec![mrv_w_expr()])
}

pub(crate) fn expr_contains_w_var(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => is_mrv_w_var(id),
        Expr::Add(ts) => ts.iter().any(expr_contains_w_var),
        Expr::Mul(fs) => fs.iter().any(expr_contains_w_var),
        Expr::Pow(b, exp) => expr_contains_w_var(b) || expr_contains_w_var(exp),
        Expr::Frac(n, d) => expr_contains_w_var(n) || expr_contains_w_var(d),
        Expr::Func(_, args) => args.iter().any(expr_contains_w_var),
        _ => false,
    }
}

pub(crate) fn expr_contains_ln_w(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && expr_contains_w_var(&args[0]) => true,
        Expr::Add(ts) => ts.iter().any(expr_contains_ln_w),
        Expr::Mul(fs) => fs.iter().any(expr_contains_ln_w),
        Expr::Pow(b, exp) => expr_contains_ln_w(b) || expr_contains_ln_w(exp),
        Expr::Frac(n, d) => expr_contains_ln_w(n) || expr_contains_ln_w(d),
        Expr::Func(_, args) => args.iter().any(expr_contains_ln_w),
        _ => false,
    }
}

fn is_ln_w(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && expr_contains_w_var(&args[0])
    )
}

pub(crate) fn is_expr_one(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_one())
}

pub(crate) fn is_expr_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

/// Write `e` as `k * ln(w) + rest` with `rest` free of `ln(w)`.
pub(crate) fn decompose_ln_w_coeff(e: &ExprArc) -> (i32, ExprArc) {
    match e.as_ref() {
        Expr::Add(ts) => {
            let mut k = 0i32;
            let mut rest = Vec::new();
            for t in ts {
                let (tk, tr) = decompose_ln_w_coeff(t);
                k = k.saturating_add(tk);
                if !is_expr_one(&tr) && !is_expr_zero(&tr) {
                    rest.push(tr);
                }
            }
            (k, rest_expr(rest))
        }
        Expr::Mul(fs) => {
            let mut k = 0i32;
            let mut rest = Vec::new();
            for f in fs {
                if is_ln_w(f) {
                    k += 1;
                    continue;
                }
                let (fk, fr) = decompose_ln_w_coeff(f);
                k = k.saturating_add(fk);
                if !is_expr_one(&fr) && !is_expr_zero(&fr) {
                    rest.push(fr);
                }
            }
            (k, rest_expr(rest))
        }
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && expr_contains_w_var(&args[0]) => {
            (1, Expr::int(1))
        }
        Expr::Int(n) if n.is_zero() => (0, Expr::int(0)),
        Expr::Int(n) if n.is_one() => (0, Expr::int(1)),
        Expr::Int(n) if n == &-BigInt::from(1) => (0, Expr::int(-1)),
        other => (0, Arc::new(other.clone())),
    }
}

fn rest_expr(mut parts: Vec<ExprArc>) -> ExprArc {
    if parts.is_empty() {
        Expr::int(1)
    } else if parts.len() == 1 {
        parts.pop().unwrap()
    } else {
        Expr::mul(parts)
    }
}
