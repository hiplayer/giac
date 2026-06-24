//! `try_integrate_*` heuristic rules (moved from `integrate.rs`, D2-1).
use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;

use giac_core::{Expr, ExprArc, FuncKind, Ident};

use crate::expr_util::{is_const_wrt, is_cos_of_var, is_sin_of_var};
use crate::integrate_helpers::{
    affine_var_coeff, is_exp_of_var, is_one, is_x_squared, is_x_squared_plus_const,
    linear_coefficient, ln_abs_expr, var_coefficient, is_var, var_to_expr,
};

// **Pipeline private** — is tan of var.
fn is_tan_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Tan, args) if args.len() == 1 && affine_var_coeff(&args[0], var).is_some()
    )
}

/// **Partial** — heuristic `try_integrate_tan_plus_tan_cubed`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_tan_plus_tan_cubed(terms: &[ExprArc], var: &Ident) -> Option<ExprArc> {
    if terms.len() != 2 {
        return None;
    }
    let mut tan_lin = None;
    let mut tan_cubed = None;
    for t in terms {
        if is_tan_of_var(t, var) {
            tan_lin = Some(t);
        } else if let Expr::Pow(base, exp) = t.as_ref() {
            if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(3)) && is_tan_of_var(base, var) {
                tan_cubed = Some(t);
            }
        }
    }
    if tan_lin.is_none() || tan_cubed.is_none() {
        return None;
    }
    let arg = match tan_lin?.as_ref() {
        Expr::Func(FuncKind::Tan, args) => Arc::clone(&args[0]),
        _ => return None,
    };
    Some(Expr::mul(vec![
        Expr::rat(1, 2),
        Expr::pow(Expr::func(FuncKind::Tan, vec![arg]), Expr::int(2)),
    ]))
}

/// **Partial** — heuristic `try_integrate_exp_over_linear_exp`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_exp_over_linear_exp(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (exp_x, den) = if is_exp_of_var(a, var) {
        (a, b)
    } else if is_exp_of_var(b, var) {
        (b, a)
    } else {
        return None;
    };
    let (c, d) = const_plus_const_times_exp(den, exp_x, var)?;
    Some(Expr::mul(vec![
        Expr::pow(Arc::clone(&d), Expr::int(-1)),
        ln_abs_expr(Expr::add(vec![c, Expr::mul(vec![d, Arc::clone(exp_x)])])),
    ]))
}

/// **Partial** — heuristic `try_integrate_exp_over_one_plus_exp2`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_exp_over_one_plus_exp2(
    exp_x: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    if !is_exp_of_var(exp_x, var) {
        return None;
    }
    let Expr::Add(ts) = den.as_ref() else {
        return None;
    };
    let mut has_one = false;
    let mut has_exp2x = false;
    for t in ts {
        if is_one(t) {
            has_one = true;
        } else if is_exp_of_double_x(t, var) || is_exp_x_squared(exp_x, t) {
            has_exp2x = true;
        } else {
            return None;
        }
    }
    if has_one && has_exp2x {
        Some(Expr::func(FuncKind::Atan, vec![Arc::clone(exp_x)]))
    } else {
        None
    }
}

// **Pipeline private** — is exp of double x.
fn is_exp_of_double_x(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1
            && matches!(
                linear_coefficient(&args[0], var).as_ref().map(|k| k.as_ref()),
                Some(Expr::Int(n)) if n == &BigInt::from(2)
            )
    )
}

// **Pipeline private** — is exp x squared.
fn is_exp_x_squared(exp_x: &ExprArc, e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Pow(base, exp) if base == exp_x && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2))
    )
}

/// **Partial** — heuristic `try_integrate_one_over_cos_squared`; **退役：** Risch / partfrac.
/// **Partial** — ∫ 1/cos(x)² dx = tan(x). **退役：** Risch / partfrac.
pub(crate) fn try_integrate_one_over_cos_squared(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    if !is_one(num) {
        return None;
    }
    let cosx = match den.as_ref() {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2)) => {
            if is_cos_of_var(base, var) {
                base
            } else {
                return None;
            }
        }
        _ => return None,
    };
    let arg = match cosx.as_ref() {
        Expr::Func(FuncKind::Cos, args) if args.len() == 1 => Arc::clone(&args[0]),
        _ => return None,
    };
    Some(Expr::func(FuncKind::Tan, vec![arg]))
}

/// **Partial** — heuristic `try_integrate_sin_over_cos_sq_frac`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_sin_over_cos_sq_frac(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    if !is_sin_of_var(num, var) {
        return None;
    }
    let cosx = match den.as_ref() {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2)) => {
            if is_cos_of_var(base, var) {
                base
            } else {
                return None;
            }
        }
        _ => return None,
    };
    Some(Expr::pow(Arc::clone(cosx), Expr::int(-1)))
}

// **Pipeline private** — const plus const times exp.
fn const_plus_const_times_exp(
    den: &ExprArc,
    exp_x: &ExprArc,
    var: &Ident,
) -> Option<(ExprArc, ExprArc)> {
    let Expr::Add(terms) = den.as_ref() else {
        return None;
    };
    let mut c = Expr::int(0);
    let mut d = None;
    for t in terms {
        if is_exp_of_var(t, var) {
            if d.is_some() {
                return None;
            }
            d = Some(Expr::int(1));
        } else if let Expr::Mul(fs) = t.as_ref() {
            if fs.len() == 2 && is_exp_of_var(&fs[0], var) && is_const_wrt(&fs[1], var) {
                if d.is_some() {
                    return None;
                }
                d = Some(Arc::clone(&fs[1]));
            } else if fs.len() == 2 && is_exp_of_var(&fs[1], var) && is_const_wrt(&fs[0], var) {
                if d.is_some() {
                    return None;
                }
                d = Some(Arc::clone(&fs[0]));
            } else {
                return None;
            }
        } else if is_const_wrt(t, var) {
            c = Expr::add(vec![c, Arc::clone(t)]);
        } else {
            return None;
        }
    }
    let d = d?;
    let _ = exp_x;
    Some((c, d))
}

// **Pipeline private** — is sin double angle.
fn is_sin_double_angle(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Sin, args) if args.len() == 1
            && linear_coefficient(&args[0], var).as_ref().map(|k| matches!(k.as_ref(), Expr::Int(n) if n == &BigInt::from(2))).unwrap_or(false)
    )
}

/// **Partial** — heuristic `try_integrate_sin2x_cos`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_sin2x_cos(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (sin2x, _cosx) = if is_sin_double_angle(a, var) && is_cos_of_var(b, var) {
        (a, b)
    } else if is_sin_double_angle(b, var) && is_cos_of_var(a, var) {
        (b, a)
    } else {
        return None;
    };
    let _ = sin2x;
    let x = var_to_expr(var);
    Some(Expr::mul(vec![
        Expr::rat(-2, 3),
        Expr::pow(Expr::func(FuncKind::Cos, vec![x]), Expr::int(3)),
    ]))
}

/// **Partial** — heuristic `try_integrate_var_shifted_sqrt`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_var_shifted_sqrt(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (k, inner) = var_times_x_squared_plus_const(a, var)?;
    let sqrt_base = shifted_sqrt_base(b, var)?;
    if inner != sqrt_base {
        return None;
    }
    Some(Expr::mul(vec![
        k,
        Expr::rat(2, 1),
        Expr::func(FuncKind::Sqrt, vec![inner]),
    ]))
}

// **Pipeline private** — var times x squared plus const.
fn var_times_x_squared_plus_const(e: &ExprArc, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    if is_var(e, var) {
        return Some((Expr::int(1), Expr::add(vec![
            Expr::pow(var_to_expr(var), Expr::int(2)),
            Expr::int(-1),
        ])));
    }
    if let Expr::Mul(fs) = e.as_ref() {
        if fs.len() == 2 && is_var(&fs[0], var) && is_const_wrt(&fs[1], var) {
            return Some((
                Arc::clone(&fs[1]),
                Expr::add(vec![
                    Expr::pow(var_to_expr(var), Expr::int(2)),
                    Expr::int(-1),
                ]),
            ));
        }
    }
    None
}

// **Pipeline private** — shifted sqrt base.
fn shifted_sqrt_base(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let base = match e.as_ref() {
        Expr::Pow(b, exp) => {
            if matches!(exp.as_ref(), Expr::Rat(r) if *r == Ratio::new((-1).into(), 2.into()))
                || matches!(exp.as_ref(), Expr::Frac(n, d) if n.is_one() && matches!(d.as_ref(), Expr::Int(i) if i == &BigInt::from(2)))
            {
                Arc::clone(b)
            } else {
                return None;
            }
        }
        Expr::Func(FuncKind::Sqrt, args) if args.len() == 1 => {
            return shifted_sqrt_base(
                &Expr::pow(Arc::clone(&args[0]), Expr::rat(1, 2)),
                var,
            );
        }
        _ => return None,
    };
    match base.as_ref() {
        Expr::Add(ts) if ts.len() == 2 => {
            if ts.iter().any(|t| is_x_squared(t, var)) && ts.iter().all(|t| is_x_squared(t, var) || is_const_wrt(t, var)) {
                Some(base)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// **Partial** — heuristic `try_integrate_sin_over_cos_squared`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_sin_over_cos_squared(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (sinx, cos_inv_sq) = if is_sin_of_var(a, var) {
        (a, b)
    } else if is_sin_of_var(b, var) {
        (b, a)
    } else {
        return None;
    };
    let cos_sq = match cos_inv_sq.as_ref() {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(2)) => base,
        _ => return None,
    };
    if !is_cos_of_var(cos_sq, var) {
        return None;
    }
    let _ = sinx;
    let x = var_to_expr(var);
    Some(Expr::mul(vec![
        Expr::int(-1),
        Expr::pow(Expr::func(FuncKind::Cos, vec![x]), Expr::int(-1)),
    ]))
}
/// **Partial** — heuristic `try_integrate_tanh_exp_form`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_tanh_exp_form(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
    if !is_exp_minus_exp_neg(a, var) || !is_exp_plus_exp_neg(b, var) {
        return None;
    }
    let _ = (a, b);
    let x = var_to_expr(var);
    let sum = Expr::add(vec![
        Expr::func(FuncKind::Exp, vec![x.clone()]),
        Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), x])]),
    ]);
    Some(ln_abs_expr(sum))
}

// **Pipeline private** — is exp minus exp neg.
fn is_exp_minus_exp_neg(e: &ExprArc, var: &Ident) -> bool {
    let Expr::Add(ts) = e.as_ref() else {
        return false;
    };
    let mut pos = false;
    let mut neg = false;
    for t in ts {
        match exp_term_sign(t, var) {
            Some(true) => pos = true,
            Some(false) => neg = true,
            None => return false,
        }
    }
    pos && neg
}

// **Pipeline private** — exp term sign.
// **Pipeline private** — `Some(true)` for `+exp(x)`, `Some(false)` for `-exp(x)` / `exp(-x)` terms.
fn exp_term_sign(e: &ExprArc, var: &Ident) -> Option<bool> {
    if is_exp_of_var(e, var) {
        return Some(true);
    }
    if is_exp_neg_of_var(e, var) {
        return Some(false);
    }
    if let Expr::Mul(fs) = e.as_ref() {
        if fs.len() == 2 {
            if matches!(fs[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) && (is_exp_of_var(&fs[1], var) || is_exp_neg_of_var(&fs[1], var)) {
                return Some(false);
            }
            if matches!(fs[1].as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) && (is_exp_of_var(&fs[0], var) || is_exp_neg_of_var(&fs[0], var)) {
                return Some(false);
            }
        }
    }
    None
}

// **Pipeline private** — is exp neg of var.
fn is_exp_neg_of_var(e: &ExprArc, var: &Ident) -> bool {
    if let Expr::Pow(base, exp) = e.as_ref() {
        if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
            return is_exp_of_var(base, var);
        }
    }
    if let Expr::Func(FuncKind::Exp, args) = e.as_ref() {
        if args.len() != 1 {
            return false;
        }
        if let Expr::Mul(fs) = args[0].as_ref() {
            if fs.len() == 2 {
                return (matches!(fs[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                    && is_var(&fs[1], var))
                    || (is_var(&fs[0], var)
                        && matches!(fs[1].as_ref(), Expr::Int(n) if n == &-BigInt::from(1)));
            }
        }
    }
    if let Expr::Mul(fs) = e.as_ref() {
        if fs.len() == 2
            && matches!(fs[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
            && is_exp_of_var(&fs[1], var)
        {
            return true;
        }
    }
    false
}

// **Pipeline private** — is exp plus exp neg.
fn is_exp_plus_exp_neg(e: &ExprArc, var: &Ident) -> bool {
    is_exp_minus_exp_neg(e, var)
}
/// **Partial** — heuristic `try_integrate_exp_trig`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_exp_trig(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
    if is_exp_of_var(a, var) && is_sin_of_var(b, var) {
        return Some(integrate_exp_sin(var));
    }
    if is_exp_of_var(b, var) && is_sin_of_var(a, var) {
        return Some(integrate_exp_sin(var));
    }
    if is_exp_of_var(a, var) && is_cos_of_var(b, var) {
        return Some(integrate_exp_cos(var));
    }
    if is_exp_of_var(b, var) && is_cos_of_var(a, var) {
        return Some(integrate_exp_cos(var));
    }
    None
}

// **Pipeline private** — integrate exp sin.
fn integrate_exp_sin(var: &Ident) -> ExprArc {
    let x = var_to_expr(var);
    let exp_x = Expr::func(FuncKind::Exp, vec![x.clone()]);
    Expr::mul(vec![
        Expr::rat(1, 2),
        exp_x,
        Expr::add(vec![
            Expr::func(FuncKind::Sin, vec![x.clone()]),
            Expr::mul(vec![
                Expr::int(-1),
                Expr::func(FuncKind::Cos, vec![x]),
            ]),
        ]),
    ])
}

// **Pipeline private** — integrate exp cos.
fn integrate_exp_cos(var: &Ident) -> ExprArc {
    let x = var_to_expr(var);
    let exp_x = Expr::func(FuncKind::Exp, vec![x.clone()]);
    Expr::mul(vec![
        Expr::rat(1, 2),
        exp_x,
        Expr::add(vec![
            Expr::func(FuncKind::Sin, vec![x.clone()]),
            Expr::func(FuncKind::Cos, vec![x]),
        ]),
    ])
}
/// **Partial** — heuristic `try_integrate_var_over_quadratic_squared`; **退役：** Risch / partfrac.
pub(crate) fn try_integrate_var_over_quadratic_squared(
    factors: &[ExprArc],
    var: &Ident,
) -> Option<ExprArc> {
    if factors.len() != 2 {
        return None;
    }
    for (num_idx, den_idx) in [(0, 1), (1, 0)] {
        let num = &factors[num_idx];
        let den_factor = &factors[den_idx];
        let k = var_coefficient(num, var)?;
        let Expr::Pow(base, exp) = den_factor.as_ref() else {
            continue;
        };
        if !matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
            continue;
        }
        if is_x_squared_plus_const(base, var) {
            return Some(Expr::mul(vec![
                Expr::rat(1, 2),
                k,
                ln_abs_expr(Arc::clone(base)),
            ]));
        }
    }
    None
}
