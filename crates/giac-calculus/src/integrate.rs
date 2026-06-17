use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;

use giac_core::{
    bigint_to_i64, Context, EvalError, Expr, ExprArc, FuncKind, Ident,
};
use giac_simplify::expand;

/// Basic integration rules (Phase 1 / GIAC-110 subset).
///
/// ## Supported
///
/// - Constants, `x`, and `x^n` for integer `n ≠ -1`
/// - Sums and constant multiples (`integrate(c*f) = c*integrate(f)`)
/// - `1/x`, `1/(ax²+bx+c)` for constant `a,b,c` (degree ≤ 2), `x/(x²+1)`
/// - Definite bounds via `eval_integrate` (4-arg `integrate(f,x,a,b)`)
///
/// ## Still `NotImplemented`
///
/// | Message | Trigger |
/// |---------|---------|
/// | `"integrate"` | Unknown top-level forms (e.g. `sin(x)`) |
/// | `"integrate frac"` | Non-constant numerator in `num/den` |
/// | `"integrate reciprocal"` | Unsupported denominator shape |
/// | `"integrate product"` | Product with multiple non-constant factors after expand |
/// | `"integrate pow"` | General power bases |
/// | `"integrate quadratic"` | Unsupported quadratic denominators |
pub fn integrate(expr: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    if let Some(r) = crate::integrate_heuristics::try_integrate_heuristic(expr, var) {
        return r;
    }
    if let Some((num, den)) = try_as_rational(expr, var) {
        if let Ok(r) = integrate_frac(&num, &den, var) {
            return Ok(r);
        }
    }
    match expr.as_ref() {
        Expr::Symbol(id) if id == var => integrate_pow(expr, &Expr::int(1), var),
        Expr::Pow(base, exp) => integrate_pow(base, exp, var),
        Expr::Frac(num, den) => integrate_frac(num, den, var),
        Expr::Mul(factors) => {
            if let Some(r) = try_integrate_var_over_quadratic_squared(factors, var) {
                return Ok(r);
            }
            integrate_mul(factors, var)
        }
        Expr::Add(terms) => {
            if let Some(r) = try_integrate_tan_plus_tan_cubed(terms, var) {
                return Ok(r);
            }
            let parts: Result<Vec<_>, _> = terms.iter().map(|t| integrate(t, var)).collect();
            Ok(Expr::add(parts?))
        }
        Expr::Func(kind, args) => integrate_func(*kind, args, var),
        _ if is_const_wrt(expr, var) => Ok(Expr::mul(vec![Arc::clone(expr), var_to_expr(var)])),
        _ => Err(EvalError::NotImplemented("integrate")),
    }
}

pub(crate) fn try_as_rational(expr: &ExprArc, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    let _ = var;
    match expr.as_ref() {
        Expr::Frac(num, den) => Some((Arc::clone(num), Arc::clone(den))),
        Expr::Pow(base, exp) => {
            let n = match exp.as_ref() {
                Expr::Int(i) => bigint_to_i64(i).ok()?,
                _ => return None,
            };
            if n >= 0 {
                return None;
            }
            let den = if n == -1 {
                Arc::clone(base)
            } else {
                Expr::pow(Arc::clone(base), Expr::int(-n))
            };
            Some((Expr::int(1), den))
        }
        Expr::Mul(factors) => {
            let mut num_parts = Vec::new();
            let mut den_parts = Vec::new();
            for f in factors {
                if let Expr::Pow(base, exp) = f.as_ref() {
                    if let Expr::Int(n) = exp.as_ref() {
                        if let Ok(ni) = bigint_to_i64(n) {
                            if ni < 0 {
                                den_parts.push(if ni == -1 {
                                    Arc::clone(base)
                                } else {
                                    Expr::pow(Arc::clone(base), Expr::int(-ni))
                                });
                                continue;
                            }
                        }
                    }
                }
                num_parts.push(Arc::clone(f));
            }
            if den_parts.is_empty() {
                return None;
            }
            let num = if num_parts.is_empty() {
                Expr::int(1)
            } else if num_parts.len() == 1 {
                num_parts.remove(0)
            } else {
                Expr::mul(num_parts)
            };
            let den = if den_parts.len() == 1 {
                den_parts.remove(0)
            } else {
                Expr::mul(den_parts)
            };
            Some((num, den))
        }
        _ => None,
    }
}

pub(crate) fn integrate_frac(num: &ExprArc, den: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    if is_exp_of_var(num, var) {
        if let Some(r) = try_integrate_exp_over_linear_exp(num, den, var) {
            return Ok(r);
        }
        if let Some(r) = try_integrate_exp_over_one_plus_exp2(num, den, var) {
            return Ok(r);
        }
    }
    if let Some(r) = try_integrate_tanh_exp_form(num, den, var) {
        return Ok(r);
    }
    if let Some(r) = try_integrate_sin_over_cos_sq_frac(num, den, var) {
        return Ok(r);
    }
    if is_const_wrt(num, var) {
        if let Ok(inner) = integrate_reciprocal(den, var) {
            if num.is_one() {
                return Ok(inner);
            }
            return Ok(Expr::mul(vec![Arc::clone(num), inner]));
        }
    }
    if let Ok(r) = crate::partfrac_integrate::integrate_const_over_rational(num, den, var) {
        return Ok(r);
    }
    if let Some(k) = var_coefficient(num, var) {
        if is_x_squared_plus_const(den, var) {
            return Ok(Expr::mul(vec![
                Expr::rat(1, 2),
                k,
                ln_abs_expr(Arc::clone(den)),
            ]));
        }
    }
    Err(EvalError::NotImplemented("integrate frac"))
}

fn integrate_reciprocal(den: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    if is_var(den, var) {
        return Ok(ln_abs(var));
    }
    if let Expr::Pow(base, exp) = den.as_ref() {
        if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
            return integrate_pow(base, exp, var);
        }
        if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2)) {
            if let Expr::Func(FuncKind::Cos, args) = base.as_ref() {
                if args.len() == 1 && is_var(&args[0], var) {
                    return Ok(Expr::func(FuncKind::Tan, vec![Arc::clone(&args[0])]));
                }
            }
        }
    }
    if let Ok(r) = crate::partfrac_integrate::integrate_one_over_quadratic(den, var) {
        return Ok(r);
    }
    if let Ok(r) = crate::partfrac_integrate::integrate_const_over_rational(&Expr::int(1), den, var)
    {
        return Ok(r);
    }
    Err(EvalError::NotImplemented("integrate reciprocal"))
}

fn integrate_func(
    kind: FuncKind,
    args: &[ExprArc],
    var: &Ident,
) -> Result<ExprArc, EvalError> {
    match (kind, args) {
        (FuncKind::Tan, [arg]) if is_var(arg, var) => Ok(Expr::mul(vec![
            Expr::int(-1),
            ln_abs_expr(Expr::func(FuncKind::Cos, vec![Arc::clone(arg)])),
        ])),
        (FuncKind::Ln, [arg]) if is_var(arg, var) => {
            let x = var_to_expr(var);
            Ok(Expr::add(vec![
                Expr::mul(vec![x.clone(), ln_abs_expr(x.clone())]),
                Expr::mul(vec![Expr::int(-1), x]),
            ]))
        }
        (FuncKind::Sin, [arg]) => {
            if let Some(k) = linear_coefficient(arg, var).or_else(|| affine_var_coeff(arg, var)) {
                return Ok(Expr::mul(vec![
                    Expr::int(-1),
                    Expr::func(FuncKind::Cos, vec![Arc::clone(arg)]),
                    Expr::pow(k, Expr::int(-1)),
                ]));
            }
            Err(EvalError::NotImplemented("integrate"))
        }
        (FuncKind::Cos, [arg]) => {
            if let Some(k) = linear_coefficient(arg, var).or_else(|| affine_var_coeff(arg, var)) {
                return Ok(Expr::mul(vec![
                    Expr::func(FuncKind::Sin, vec![Arc::clone(arg)]),
                    Expr::pow(k, Expr::int(-1)),
                ]));
            }
            Err(EvalError::NotImplemented("integrate"))
        }
        _ => Err(EvalError::NotImplemented("integrate")),
    }
}

fn linear_coefficient(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
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

fn affine_var_coeff(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
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

fn is_tan_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Tan, args) if args.len() == 1 && affine_var_coeff(&args[0], var).is_some()
    )
}

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

fn is_exp_x_squared(exp_x: &ExprArc, e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Pow(base, exp) if base == exp_x && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2))
    )
}

/// ∫ 1/cos(x)² dx = tan(x).
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

fn is_sin_double_angle(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Sin, args) if args.len() == 1
            && linear_coefficient(&args[0], var).as_ref().map(|k| matches!(k.as_ref(), Expr::Int(n) if n == &BigInt::from(2))).unwrap_or(false)
    )
}

fn try_integrate_sin2x_cos(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let (sin2x, cosx) = if is_sin_double_angle(a, var) && is_cos_of_var(b, var) {
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

fn try_integrate_var_shifted_sqrt(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
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

fn try_integrate_sin_over_cos_squared(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
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

/// `Some(true)` for `+exp(x)`, `Some(false)` for `-exp(x)` / `exp(-x)` terms.
fn exp_term_sign(e: &ExprArc, var: &Ident) -> Option<bool> {
    if is_exp_of_var(e, var) {
        return Some(true);
    }
    if is_exp_neg_of_var(e, var) {
        return Some(false);
    }
    if let Expr::Mul(fs) = e.as_ref() {
        if fs.len() == 2 {
            if matches!(fs[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
                if is_exp_of_var(&fs[1], var) || is_exp_neg_of_var(&fs[1], var) {
                    return Some(false);
                }
            }
            if matches!(fs[1].as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
                if is_exp_of_var(&fs[0], var) || is_exp_neg_of_var(&fs[0], var) {
                    return Some(false);
                }
            }
        }
    }
    None
}

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

fn is_exp_plus_exp_neg(e: &ExprArc, var: &Ident) -> bool {
    is_exp_minus_exp_neg(e, var)
}

fn integrate_sin_squared(var: &Ident) -> ExprArc {
    let x = var_to_expr(var);
    Expr::add(vec![
        Expr::mul(vec![Expr::rat(1, 2), x.clone()]),
        Expr::mul(vec![
            Expr::rat(-1, 4),
            Expr::func(FuncKind::Sin, vec![Expr::mul(vec![Expr::int(2), x])]),
        ]),
    ])
}

fn integrate_cos_squared(var: &Ident) -> ExprArc {
    let x = var_to_expr(var);
    Expr::add(vec![
        Expr::mul(vec![Expr::rat(1, 2), x.clone()]),
        Expr::mul(vec![
            Expr::rat(1, 4),
            Expr::func(FuncKind::Sin, vec![Expr::mul(vec![Expr::int(2), x])]),
        ]),
    ])
}

fn integrate_sin_cos_product(var: &Ident) -> ExprArc {
    let x = var_to_expr(var);
    Expr::mul(vec![
        Expr::rat(1, 2),
        Expr::pow(Expr::func(FuncKind::Sin, vec![x]), Expr::int(2)),
    ])
}

fn integrate_x_ln(var: &Ident) -> ExprArc {
    let x = var_to_expr(var);
    Expr::add(vec![
        Expr::mul(vec![
            Expr::rat(1, 2),
            Expr::pow(x.clone(), Expr::int(2)),
            ln_abs_expr(x.clone()),
        ]),
        Expr::mul(vec![Expr::rat(-1, 4), Expr::pow(x, Expr::int(2))]),
    ])
}

fn try_integrate_exp_trig(a: &ExprArc, b: &ExprArc, var: &Ident) -> Option<ExprArc> {
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

pub(crate) fn is_exp_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Exp, args) if args.len() == 1 && is_var(&args[0], var))
}

fn integrate_mul(factors: &[ExprArc], var: &Ident) -> Result<ExprArc, EvalError> {
    if factors.len() == 1 {
        return integrate(&factors[0], var);
    }
    if let Some(r) = try_integrate_var_over_quadratic_squared(factors, var) {
        return Ok(r);
    }
    if factors.len() == 2 {
        if (is_sin_of_var(&factors[0], var) && is_cos_of_var(&factors[1], var))
            || (is_sin_of_var(&factors[1], var) && is_cos_of_var(&factors[0], var))
        {
            return Ok(integrate_sin_cos_product(var));
        }
        if (is_var(&factors[0], var) && is_ln_of_var(&factors[1], var))
            || (is_var(&factors[1], var) && is_ln_of_var(&factors[0], var))
        {
            return Ok(integrate_x_ln(var));
        }
        if let Some(r) = try_integrate_exp_trig(&factors[0], &factors[1], var) {
            return Ok(r);
        }
        if let Some(r) = try_integrate_exp_over_linear_exp(&factors[0], &factors[1], var) {
            return Ok(r);
        }
        if let Some(r) = try_integrate_sin2x_cos(&factors[0], &factors[1], var) {
            return Ok(r);
        }
        if let Some(r) = try_integrate_var_shifted_sqrt(&factors[0], &factors[1], var) {
            return Ok(r);
        }
        if let Some(r) = try_integrate_sin_over_cos_squared(&factors[0], &factors[1], var) {
            return Ok(r);
        }
        if let Some(r) = try_integrate_tanh_exp_form(&factors[0], &factors[1], var) {
            return Ok(r);
        }
    }
    let const_part: Vec<ExprArc> = factors
        .iter()
        .filter(|f| is_const_wrt(f, var))
        .cloned()
        .collect();
    let var_part: Vec<ExprArc> = factors
        .iter()
        .filter(|f| !is_const_wrt(f, var))
        .cloned()
        .collect();
    match var_part.len() {
        0 => Ok(Expr::mul(vec![Expr::mul(factors.to_vec()), var_to_expr(var)])),
        1 => {
            let inner = integrate(&var_part[0], var)?;
            if const_part.is_empty() {
                Ok(inner)
            } else {
                Ok(Expr::mul(std::iter::once(inner).chain(const_part).collect()))
            }
        }
        _ => {
            for (i, f) in factors.iter().enumerate() {
                if let Expr::Pow(base, exp) = f.as_ref() {
                    if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
                        let others: Vec<ExprArc> = factors
                            .iter()
                            .enumerate()
                            .filter(|(j, _)| *j != i)
                            .map(|(_, g)| Arc::clone(g))
                            .collect();
                        if others.iter().all(|g| is_const_wrt(g, var)) {
                            let num = if others.is_empty() {
                                Expr::int(1)
                            } else if others.len() == 1 {
                                Arc::clone(&others[0])
                            } else {
                                Expr::mul(others)
                            };
                            return integrate_frac(&num, base, var);
                        }
                    }
                }
            }
            if var_part.len() == 2 {
                if let (Expr::Add(as_), Expr::Add(bs)) =
                    (var_part[0].as_ref(), var_part[1].as_ref())
                {
                    let mut terms = Vec::new();
                    for a in as_ {
                        for b in bs {
                            terms.push(Expr::mul(vec![Arc::clone(a), Arc::clone(b)]));
                        }
                    }
                    return integrate(&Expr::add(terms), var);
                }
                if var_part[0] == var_part[1] {
                    return integrate_pow(&var_part[0], &Expr::int(2), var);
                }
            }
            let product = Expr::mul(var_part);
            let ctx = Context::default();
            let expanded = expand(product.as_ref(), &ctx)?;
            let var_factors = count_var_factors(expanded.as_ref(), var);
            if var_factors > 1 {
                return Err(EvalError::NotImplemented("integrate product"));
            }
            integrate(&expanded, var)
        }
    }
}

fn count_var_factors(e: &Expr, var: &Ident) -> usize {
    match e {
        Expr::Mul(fs) => fs.iter().filter(|f| !is_const_wrt(f, var)).count(),
        _ => if is_const_wrt(&Arc::new(e.clone()), var) {
            0
        } else {
            1
        },
    }
}

fn integrate_pow(base: &ExprArc, exp: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
        return integrate_reciprocal(base, var);
    }
    if let Expr::Func(FuncKind::Sin, args) = base.as_ref() {
        if args.len() == 1
            && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2))
            && is_var(&args[0], var)
        {
            return Ok(integrate_sin_squared(var));
        }
    }
    if let Expr::Func(FuncKind::Cos, args) = base.as_ref() {
        if args.len() == 1 && is_var(&args[0], var) {
            if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2)) {
                return Ok(integrate_cos_squared(var));
            }
            if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(2)) {
                return Ok(Expr::func(FuncKind::Tan, vec![Arc::clone(&args[0])]));
            }
        }
    }
    if is_var(base, var) {
        if let Expr::Int(n) = exp.as_ref() {
            let n = bigint_to_i64(n)?;
            if n == -1 {
                return Ok(ln_abs(var));
            }
            if n >= 0 {
                return Ok(Expr::mul(vec![
                    Expr::rat(1, n + 1),
                    Expr::pow(Arc::clone(base), Expr::int(n + 1)),
                ]));
            }
        }
    }
    if let Expr::Add(_) = base.as_ref() {
        if let Expr::Int(n) = exp.as_ref() {
            if bigint_to_i64(n)? >= 0 {
                let ctx = Context::default();
                let powered = Expr::pow(Arc::clone(base), Arc::clone(exp));
                let expanded = expand(powered.as_ref(), &ctx)?;
                return integrate(&expanded, var);
            }
        }
    }
    Err(EvalError::NotImplemented("integrate pow"))
}

fn is_x_squared(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2)))
}

fn is_x_squared_plus_const(den: &ExprArc, var: &Ident) -> bool {
    match den.as_ref() {
        Expr::Add(terms) if terms.len() == 2 => {
            terms.iter().any(|t| is_x_squared(t, var))
                && terms.iter().all(|t| is_x_squared(t, var) || is_const_wrt(t, var))
        }
        _ => false,
    }
}

fn var_coefficient(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
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

fn try_integrate_var_over_quadratic_squared(
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

fn ln_abs(var: &Ident) -> ExprArc {
    ln_abs_expr(var_to_expr(var))
}

pub(crate) fn ln_abs_expr(arg: ExprArc) -> ExprArc {
    Expr::func(FuncKind::Ln, vec![Expr::func(FuncKind::Abs, vec![arg])])
}

pub(crate) fn var_to_expr(var: &Ident) -> ExprArc {
    Expr::sym(var.as_str())
}

fn is_sin_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_var(&args[0], var))
}

fn is_cos_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Cos, args) if args.len() == 1 && is_var(&args[0], var))
}

fn is_ln_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Ln, args) if args.len() == 1 && is_var(&args[0], var))
}

pub(crate) fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_one(e: &ExprArc) -> bool {
    e.is_one()
}

pub(crate) fn is_const_wrt(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id != var,
        Expr::Int(_) | Expr::Rat(_) => true,
        Expr::Frac(n, d) => is_const_wrt(n, var) && is_const_wrt(d, var),
        Expr::Add(ts) => ts.iter().all(|t| is_const_wrt(t, var)),
        Expr::Mul(fs) => fs.iter().all(|f| is_const_wrt(f, var)),
        Expr::Pow(b, _) => is_const_wrt(b, var),
        Expr::Func(_, args) => args.iter().all(|a| is_const_wrt(a, var)),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::{eval, format_expr, Context, FuncKind};

    #[test]
    fn giac223_tanh_exp_frac() {
        let x = Ident::new("x");
        let ex = Expr::func(FuncKind::Exp, vec![Expr::sym("x")]);
        let ex_neg = Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]);
        let num = Expr::add(vec![
            Arc::clone(&ex),
            Expr::mul(vec![Expr::int(-1), ex_neg.clone()]),
        ]);
        let den = Expr::add(vec![ex, ex_neg]);
        let e = Arc::new(Expr::Frac(num, den));
        assert!(integrate(&e, &x).is_ok());
    }

    #[test]
    fn giac223_exp_over_linear() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![
            Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
            Expr::pow(
                Expr::add(vec![Expr::int(3), Expr::mul(vec![Expr::int(2), Expr::func(FuncKind::Exp, vec![Expr::sym("x")])])]),
                Expr::int(-1),
            ),
        ]);
        assert!(integrate(&e, &x).is_ok());
    }

    #[test]
    fn integrate_reciprocal() {
        let x = Ident::new("x");
        let e = Expr::pow(Expr::sym("x"), Expr::int(-1));
        let r = integrate(&e, &x).unwrap();
        let ev = eval(r.as_ref(), &Context::default()).unwrap();
        assert_eq!(format_expr(ev.as_ref()), "ln(abs(x))");
    }

    #[test]
    fn integrate_one_minus_x_fourth() {
        let x = Ident::new("x");
        let den = Expr::add(vec![
            Expr::int(1),
            Expr::mul(vec![Expr::int(-1), Expr::pow(Expr::sym("x"), Expr::int(4))]),
        ]);
        let e = Expr::pow(den, Expr::int(-1));
        let r = integrate(&e, &x).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln(abs(x-1))"));
        assert!(s.contains("atan"));
    }

    #[test]
    fn integrate_constant_and_sum() {
        let x = Ident::new("x");
        let c = Expr::int(5);
        let r = integrate(&c, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "5*x");

        let sum = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(-1)),
            Expr::int(1),
        ]);
        let r = integrate(&sum, &x).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln(abs(x))"));
        assert!(s.contains("x"));
    }

    #[test]
    fn integrate_const_times_x() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![Expr::int(3), Expr::sym("x")]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(1/2*x^2)*3");
    }

    #[test]
    fn integrate_x_and_x_squared() {
        let x = Ident::new("x");
        let r = integrate(&Expr::sym("x"), &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/2*x^2");
        let r = integrate(&Expr::pow(Expr::sym("x"), Expr::int(2)), &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/3*x^3");
    }

    #[test]
    fn integrate_one() {
        let x = Ident::new("x");
        let r = integrate(&Expr::int(1), &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1*x");
    }

    #[test]
    fn integrate_one_over_one_plus_x_squared() {
        let x = Ident::new("x");
        let den = Expr::add(vec![Expr::int(1), Expr::pow(Expr::sym("x"), Expr::int(2))]);
        let e = Expr::pow(den, Expr::int(-1));
        let r = integrate(&e, &x).unwrap();
        assert!(format_expr(r.as_ref()).contains("atan"));
    }

    #[test]
    fn integrate_frac_with_constant_numerator() {
        let x = Ident::new("x");
        let e = Arc::new(Expr::Frac(Expr::int(2), Expr::sym("x")));
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2*ln(abs(x))");
    }

    #[test]
    fn integrate_one_over_x_squared_plus_one_term_order() {
        let x = Ident::new("x");
        let den = Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(1)]);
        let e = Expr::pow(den, Expr::int(-1));
        let r = integrate(&e, &x).unwrap();
        assert!(format_expr(r.as_ref()).contains("atan"));
    }

    #[test]
    fn integrate_x_over_x_squared_plus_one() {
        let x = Ident::new("x");
        let den = Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(1)]);
        let e = Arc::new(Expr::Frac(Expr::sym("x"), den));
        let r = integrate(&e, &x).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln(abs("));
        assert!(s.contains("x^2+1"));
    }

    #[test]
    fn integrate_unsupported_returns_not_implemented() {
        let x = Ident::new("x");
        let e = Expr::func(FuncKind::Atan, vec![Expr::sym("x")]);
        assert!(matches!(
            integrate(&e, &x),
            Err(EvalError::NotImplemented(_))
        ));
    }

    #[test]
    fn integrate_const_over_quadratic() {
        let x = Ident::new("x");
        let den = Expr::add(vec![Expr::int(4), Expr::pow(Expr::sym("x"), Expr::int(2))]);
        let e = Expr::pow(den, Expr::int(-1));
        let r = integrate(&e, &x).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("atan"));
    }

    #[test]
    fn integrate_const_times_reciprocal() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![Expr::int(3), Expr::pow(Expr::sym("x"), Expr::int(-1))]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "3*ln(abs(x))");
    }

    #[test]
    fn integrate_product_of_consts_only() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![Expr::int(2), Expr::int(3)]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(2*3)*x");
    }

    #[test]
    fn integrate_definite_bounds() {
        let ctx = crate::plugin::xcas_default();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![
                Expr::pow(Expr::add(vec![Expr::int(1), Expr::sym("x")]), Expr::int(2)),
                Expr::sym("x"),
                Expr::int(-1),
                Expr::int(1),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "8/3");
    }

    #[test]
    fn integrate_product_one_plus_x_squared() {
        let x = Ident::new("x");
        let factor = Expr::add(vec![Expr::int(1), Expr::sym("x")]);
        let e = Expr::mul(vec![factor.clone(), factor]);
        let r = integrate(&e, &x).unwrap();
        assert!(format_expr(r.as_ref()).contains("x^3"));
    }

    #[test]
    fn integrate_x_over_x_squared_plus_one_squared() {
        let x = Ident::new("x");
        let e = Arc::new(Expr::Frac(
            Expr::sym("x"),
            Expr::pow(
                Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(1)]),
                Expr::int(2),
            ),
        ));
        let r = integrate(&e, &x).expect("integrate");
        let s = format_expr(r.as_ref());
        assert!(s.contains("x^2+1") || s.contains("x^2 + 1"), "got {s}");
    }

    #[test]
    fn integrate_x_squared_reciprocal() {
        let x = Ident::new("x");
        let den = Expr::pow(Expr::sym("x"), Expr::int(2));
        let e = Expr::pow(den, Expr::int(-1));
        match integrate(&e, &x) {
            Ok(r) => {
                let s = format_expr(r.as_ref());
                assert!(
                    s.contains("x^(-1)") || s.contains("x^-1") || s.contains("1/x") || s.contains("^-1"),
                    "got {s}"
                );
            }
            Err(e) => panic!("integrate 1/x^2 failed: {e:?}"),
        }
    }

    #[test]
    fn integrate_sin_times_cos() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![
            Expr::func(FuncKind::Sin, vec![Expr::sym("x")]),
            Expr::func(FuncKind::Cos, vec![Expr::sym("x")]),
        ]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/2*sin(x)^2");
    }

    #[test]
    fn integrate_one_over_cos_squared() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![
            Expr::int(1),
            Expr::pow(
                Expr::pow(Expr::func(FuncKind::Cos, vec![Expr::sym("x")]), Expr::int(2)),
                Expr::int(-1),
            ),
        ]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "tan(x)");
    }

    #[test]
    fn integrate_not_implemented_messages() {
        let x = Ident::new("x");
        let frac = Arc::new(Expr::Frac(
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::sym("x"),
        ));
        assert!(integrate(&frac, &x).is_ok());
        let e = Expr::pow(
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(2)]),
            Expr::int(-1),
        );
        assert!(integrate(&e, &x).is_ok());
        let e = Expr::pow(
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(3)), Expr::int(2)]),
            Expr::int(-1),
        );
        assert!(matches!(
            integrate(&e, &x),
            Err(EvalError::NotImplemented(_))
        ));
    }
}
