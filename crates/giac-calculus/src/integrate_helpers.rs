//! Shared integration helpers (breaks `integrate` ↔ `integrate_heuristics` cycle).
use std::sync::Arc;

use num_bigint::BigInt;

use giac_core::{Expr, ExprArc, FuncKind, Ident};

use crate::expr_util::is_const_wrt;

pub(crate) use crate::expr_util::{is_var, var_to_expr};

// **Pipeline private** — linear coefficient.
pub(crate) fn linear_coefficient(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
    if is_var(e, var) {
        return Some(Expr::int(1));
    }
    if let Expr::Mul(factors) = e.as_ref() {
        if factors.len() == 2 {
            if is_var(&factors[0], var) && is_const_wrt(&factors[1], var) {
                return Some(Arc::clone(&factors[1]));
            }
            if is_var(&factors[1], var) && is_const_wrt(&factors[0], var) {
                return Some(Arc::clone(&factors[0]));
            }
            if matches!(factors[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
                if let Some(c) = linear_coefficient(&factors[1], var) {
                    return Some(Expr::mul(vec![Expr::int(-1), c]));
                }
            }
            if matches!(factors[1].as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
                if let Some(c) = linear_coefficient(&factors[0], var) {
                    return Some(Expr::mul(vec![Expr::int(-1), c]));
                }
            }
        }
    }
    None
}

// **Pipeline private** — affine var coeff.
pub(crate) fn affine_var_coeff(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Symbol(id) if id == var => Some(Expr::int(1)),
        Expr::Mul(factors) if factors.len() == 2 => {
            if is_var(&factors[0], var) && is_const_wrt(&factors[1], var) {
                return Some(Arc::clone(&factors[1]));
            }
            if is_var(&factors[1], var) && is_const_wrt(&factors[0], var) {
                return Some(Arc::clone(&factors[0]));
            }
            None
        }
        Expr::Add(terms) => {
            let mut k = Expr::int(0);
            for t in terms {
                if is_var(t, var) {
                    k = Expr::add(vec![k, Expr::int(1)]);
                } else if let Some(c) = linear_coefficient(t, var) {
                    k = Expr::add(vec![k, c]);
                } else if !is_const_wrt(t, var) {
                    return None;
                }
            }
            if k.as_ref().is_zero() {
                None
            } else {
                Some(k)
            }
        }
        _ if is_const_wrt(e, var) => None,
        _ => None,
    }
}
/// **Stable** — detect `exp(var)` form.
pub(crate) fn is_exp_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Exp, args) if args.len() == 1 && is_var(&args[0], var))
}
// **Pipeline private** — is x squared.
pub(crate) fn is_x_squared(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2)))
}

// **Pipeline private** — is x squared plus const.
pub(crate) fn is_x_squared_plus_const(den: &ExprArc, var: &Ident) -> bool {
    match den.as_ref() {
        Expr::Add(terms) if terms.len() == 2 => {
            terms.iter().any(|t| is_x_squared(t, var))
                && terms.iter().all(|t| is_x_squared(t, var) || is_const_wrt(t, var))
        }
        _ => false,
    }
}

// **Pipeline private** — var coefficient.
pub(crate) fn var_coefficient(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Symbol(id) if id == var => Some(Expr::int(1)),
        Expr::Mul(factors) if factors.len() == 2 => {
            if is_var(&factors[0], var) && is_const_wrt(&factors[1], var) {
                Some(Arc::clone(&factors[1]))
            } else if is_var(&factors[1], var) && is_const_wrt(&factors[0], var) {
                Some(Arc::clone(&factors[0]))
            } else {
                None
            }
        }
        _ => None,
    }
}
/// **Stable** — build `ln(abs(arg))` expression.
pub(crate) fn ln_abs_expr(arg: ExprArc) -> ExprArc {
    Expr::func(FuncKind::Ln, vec![Expr::func(FuncKind::Abs, vec![arg])])
}

// **Pipeline private** — is one.
pub(crate) fn is_one(e: &ExprArc) -> bool {
    e.is_one()
}
