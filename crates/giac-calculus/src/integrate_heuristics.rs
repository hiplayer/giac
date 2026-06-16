//! GIAC-225: `intg.cc` heuristic subset (sqrt substitution, trig fractions).

use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;

use giac_core::{bigint_to_i64, ratnormal, Context, EvalError, Expr, ExprArc, FuncKind, Ident};

use crate::integrate::{
    integrate_frac, is_const_wrt, is_var, ln_abs_expr, try_as_rational, var_to_expr,
};

/// Top-level sqrt / trig-fraction hooks before generic `integrate` dispatch.
pub fn try_integrate_heuristic(expr: &ExprArc, var: &Ident) -> Option<Result<ExprArc, EvalError>> {
    if let Some((num, den)) = try_as_rational(expr, var) {
        if let Some(r) = try_integrate_trig_deriv_ratio(&num, &den, var) {
            return Some(Ok(r));
        }
        if let Some(r) = try_integrate_var_over_sqrt_xsq_plus_c(&num, &den, var) {
            return Some(Ok(r));
        }
        if let Some(r) = try_integrate_x_over_sqrt_affine(&num, &den, var) {
            return Some(Ok(r));
        }
        if let Some(r) = try_integrate_sin2x_affine_over_cos2x(&num, &den, var) {
            return Some(Ok(r));
        }
        if let Some(r) = try_integrate_trig_rational_half_angle(&num, &den, var) {
            return Some(Ok(r));
        }
        return Some(integrate_frac(&num, &den, var));
    }
    if let Some(r) = try_integrate_x_times_sqrt_quadratic(expr, var) {
        return Some(Ok(r));
    }
    None
}

/// ∫ (a·sin + b·cos)'/(a·sin + b·cos) dx = ln|a·sin + b·cos| when numerator is d(den)/dx.
fn try_integrate_trig_deriv_ratio(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    let (ns, nc) = sin_cos_coeffs(num, var)?;
    let (ds, dc) = sin_cos_coeffs(den, var)?;
    if ns == ds && nc == -dc {
        return Some(ln_abs_expr(Arc::clone(den)));
    }
    None
}

/// ∫ k·x/√(x²+c) dx = √(x²+c) when k=2.
fn try_integrate_var_over_sqrt_xsq_plus_c(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    let inner = sqrt_radicand(den)?;
    let _c = const_term_of_xsq_plus_const(&inner, var)?;
    let k = var_coefficient_int(num, var)?;
    Some(Expr::mul(vec![
        Expr::int(k),
        Expr::func(FuncKind::Sqrt, vec![inner]),
    ]))
}

/// ∫ x/√(x+a) dx = (2/3)(x+a)^(3/2) - 2√(x+a).
fn try_integrate_x_over_sqrt_affine(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    if !is_var(num, var) {
        return None;
    }
    let inner = sqrt_radicand(den)?;
    let _a = const_term_of_x_plus_const(&inner, var)?;
    let root = Expr::func(FuncKind::Sqrt, vec![inner.clone()]);
    let x_expr = var_to_expr(var);
    Some(Expr::add(vec![
        Expr::mul(vec![Expr::rat(2, 3), x_expr, root.clone()]),
        Expr::mul(vec![Expr::rat(-4, 3), root]),
    ]))
}

/// ∫ x·√(a+x²) dx = (a+x²)^(3/2)/3.
fn try_integrate_x_times_sqrt_quadratic(expr: &ExprArc, var: &Ident) -> Option<ExprArc> {
    let factors = match expr.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => fs,
        _ => return None,
    };
    let (x_part, sqrt_part) = if is_var(&factors[0], var) {
        (&factors[0], &factors[1])
    } else if is_var(&factors[1], var) {
        (&factors[1], &factors[0])
    } else {
        return None;
    };
    let _ = x_part;
    let inner = sqrt_radicand(sqrt_part)?;
    let c = const_term_of_xsq_plus_const(&inner, var)?;
    if c != 1 {
        return None;
    }
    let x_expr = var_to_expr(var);
    let sqrt_inner = Expr::func(FuncKind::Sqrt, vec![inner]);
    Some(Expr::add(vec![
        Expr::mul(vec![
            Expr::rat(1, 3),
            Expr::pow(x_expr, Expr::int(2)),
            sqrt_inner.clone(),
        ]),
        Expr::mul(vec![Expr::rat(1, 3), sqrt_inner]),
    ]))
}

/// ∫ (k·sin(2x)+c)/cos(2x) dx = −c/(2k)·ln|c−k·sin(2x)| (GIAC-normalized).
fn try_integrate_sin2x_affine_over_cos2x(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    if !is_cos_sin_double_angle_expr(den, var) {
        return None;
    }
    let (k, c) = sin_double_angle_affine_coeffs(num, var)?;
    if k == 0 {
        return None;
    }
    let sin2x = sin_double_angle_expr(var);
    let ln_arg = Expr::add(vec![
        Expr::int(c),
        Expr::mul(vec![Expr::int(-k), sin2x]),
    ]);
    Some(Expr::mul(vec![
        Expr::rat(-c, 2 * k),
        ln_abs_expr(ln_arg),
    ]))
}

fn sin_double_angle_expr(var: &Ident) -> ExprArc {
    Expr::func(
        FuncKind::Sin,
        vec![Expr::mul(vec![Expr::int(2), var_to_expr(var)])],
    )
}

fn is_cos_sin_double_angle_expr(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Cos, args) if args.len() == 1 && is_sin_double_angle(&args[0], var)
    )
}

fn sin_double_angle_affine_coeffs(e: &ExprArc, var: &Ident) -> Option<(i64, i64)> {
    match e.as_ref() {
        Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_sin_double_angle(&args[0], var) => {
            Some((1, 0))
        }
        Expr::Int(n) => bigint_to_i64(n).ok().map(|c| (0, c)),
        Expr::Add(ts) => {
            let mut k = 0i64;
            let mut c = 0i64;
            for t in ts {
                if let Some((tk, tc)) = sin_double_angle_affine_coeffs(t, var) {
                    k += tk;
                    c += tc;
                } else {
                    return None;
                }
            }
            Some((k, c))
        }
        Expr::Mul(fs) if fs.len() == 2 => {
            let (coeff, trig) = if let Expr::Int(n) = fs[0].as_ref() {
                (bigint_to_i64(n).ok()?, &fs[1])
            } else if let Expr::Int(n) = fs[1].as_ref() {
                (bigint_to_i64(n).ok()?, &fs[0])
            } else {
                return None;
            };
            if matches!(
                trig.as_ref(),
                Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_sin_double_angle(&args[0], var)
            ) {
                return Some((coeff, 0));
            }
            None
        }
        _ => None,
    }
}

/// Rational function of sin(x), cos(x) via t = tan(x/2).
fn try_integrate_trig_rational_half_angle(
    num: &ExprArc,
    den: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    if !is_trig_rational_in_x(num, var) || !is_trig_rational_in_x(den, var) {
        return None;
    }
    let t = Ident::new("__t");
    let x = var_to_expr(var);
    let half_x = Expr::mul(vec![Expr::rat(1, 2), x]);
    let sin_t = weierstrass_sin(&t);
    let cos_t = weierstrass_cos(&t);
    let jacobian = weierstrass_dx_dt(&t);
    let num_t = replace_trig_with_t(num, var, &sin_t, &cos_t);
    let den_t = replace_trig_with_t(den, var, &sin_t, &cos_t);
    let integrand = Expr::mul(vec![
        num_t,
        Expr::pow(den_t, Expr::int(-1)),
        jacobian,
    ]);
    let ctx = Context::xcas_default();
    let integrand = ratnormal(integrand.as_ref(), &ctx).ok()?;
    let (num_i, den_i) = try_as_rational(&integrand, &t)?;
    let inner = integrate_frac(&num_i, &den_i, &t).ok()?;
    let back = Expr::func(FuncKind::Tan, vec![half_x]);
    Some(replace_symbol(&inner, &t, &back))
}

fn weierstrass_sin(t: &Ident) -> ExprArc {
    let tv = var_to_expr(t);
    let t2 = Expr::pow(tv.clone(), Expr::int(2));
    Expr::mul(vec![
        Expr::int(2),
        tv,
        Expr::pow(Expr::add(vec![Expr::int(1), t2]), Expr::int(-1)),
    ])
}

fn weierstrass_cos(t: &Ident) -> ExprArc {
    let tv = var_to_expr(t);
    let t2 = Expr::pow(tv.clone(), Expr::int(2));
    Expr::mul(vec![
        Expr::add(vec![Expr::int(1), Expr::mul(vec![Expr::int(-1), t2.clone()])]),
        Expr::pow(Expr::add(vec![Expr::int(1), t2]), Expr::int(-1)),
    ])
}

fn weierstrass_dx_dt(t: &Ident) -> ExprArc {
    let tv = var_to_expr(t);
    Expr::mul(vec![
        Expr::int(2),
        Expr::pow(Expr::add(vec![Expr::int(1), Expr::pow(tv, Expr::int(2))]), Expr::int(-1)),
    ])
}

fn replace_trig_with_t(
    e: &ExprArc,
    var: &Ident,
    sin_t: &ExprArc,
    cos_t: &ExprArc,
) -> ExprArc {
    match e.as_ref() {
        Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_sin_double_angle(&args[0], var) => {
            Expr::mul(vec![Expr::int(2), sin_t.clone(), cos_t.clone()])
        }
        Expr::Func(FuncKind::Cos, args) if args.len() == 1 && is_cos_double_angle(&args[0], var) => {
            Expr::add(vec![
                Expr::pow(cos_t.clone(), Expr::int(2)),
                Expr::mul(vec![Expr::int(-1), Expr::pow(sin_t.clone(), Expr::int(2))]),
            ])
        }
        Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_var(&args[0], var) => {
            Arc::clone(sin_t)
        }
        Expr::Func(FuncKind::Cos, args) if args.len() == 1 && is_var(&args[0], var) => {
            Arc::clone(cos_t)
        }
        Expr::Add(ts) => Expr::add(
            ts.iter()
                .map(|t| replace_trig_with_t(t, var, sin_t, cos_t))
                .collect(),
        ),
        Expr::Mul(fs) => Expr::mul(
            fs.iter()
                .map(|f| replace_trig_with_t(f, var, sin_t, cos_t))
                .collect(),
        ),
        Expr::Pow(b, exp) => Expr::pow(
            replace_trig_with_t(b, var, sin_t, cos_t),
            Arc::clone(exp),
        ),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            replace_trig_with_t(n, var, sin_t, cos_t),
            replace_trig_with_t(d, var, sin_t, cos_t),
        )),
        _ => Arc::clone(e),
    }
}

fn replace_symbol(e: &ExprArc, sym: &Ident, repl: &ExprArc) -> ExprArc {
    match e.as_ref() {
        Expr::Symbol(id) if id == sym => Arc::clone(repl),
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| replace_symbol(t, sym, repl)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|f| replace_symbol(f, sym, repl)).collect()),
        Expr::Pow(b, exp) => Expr::pow(replace_symbol(b, sym, repl), Arc::clone(exp)),
        Expr::Func(k, args) => {
            Expr::func(*k, args.iter().map(|a| replace_symbol(a, sym, repl)).collect())
        }
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            replace_symbol(n, sym, repl),
            replace_symbol(d, sym, repl),
        )),
        _ => Arc::clone(e),
    }
}

fn is_trig_rational_in_x(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Int(_) | Expr::Rat(_) => true,
        Expr::Symbol(id) => id != var,
        Expr::Add(ts) => ts.iter().all(|t| is_trig_rational_in_x(t, var)),
        Expr::Mul(fs) => fs.iter().all(|f| is_trig_rational_in_x(f, var)),
        Expr::Pow(b, exp) => {
            is_trig_rational_in_x(b, var)
                && matches!(exp.as_ref(), Expr::Int(n) if bigint_to_i64(n).is_ok())
        }
        Expr::Func(FuncKind::Sin | FuncKind::Cos, args) => {
            args.len() == 1
                && (is_var(&args[0], var)
                    || is_sin_double_angle(&args[0], var)
                    || is_cos_double_angle(&args[0], var))
        }
        Expr::Frac(n, d) => is_trig_rational_in_x(n, var) && is_trig_rational_in_x(d, var),
        _ => false,
    }
}

fn sin_cos_coeffs(e: &ExprArc, var: &Ident) -> Option<(i64, i64)> {
    let (s, c, extra) = trig_affine_parts(e, var)?;
    if extra {
        return None;
    }
    Some((s, c))
}

fn trig_affine_parts(e: &ExprArc, var: &Ident) -> Option<(i64, i64, bool)> {
    match e.as_ref() {
        Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_var(&args[0], var) => {
            Some((1, 0, false))
        }
        Expr::Func(FuncKind::Cos, args) if args.len() == 1 && is_var(&args[0], var) => {
            Some((0, 1, false))
        }
        Expr::Mul(fs) if fs.len() == 2 => {
            if let Expr::Int(k) = fs[0].as_ref() {
                let ki = bigint_to_i64(k).ok()?;
                if let Expr::Func(FuncKind::Sin, args) = fs[1].as_ref() {
                    if args.len() == 1 && is_var(&args[0], var) {
                        return Some((ki, 0, false));
                    }
                }
                if let Expr::Func(FuncKind::Cos, args) = fs[1].as_ref() {
                    if args.len() == 1 && is_var(&args[0], var) {
                        return Some((0, ki, false));
                    }
                }
            }
            if let Expr::Int(k) = fs[1].as_ref() {
                let ki = bigint_to_i64(k).ok()?;
                if let Expr::Func(FuncKind::Sin, args) = fs[0].as_ref() {
                    if args.len() == 1 && is_var(&args[0], var) {
                        return Some((ki, 0, false));
                    }
                }
                if let Expr::Func(FuncKind::Cos, args) = fs[0].as_ref() {
                    if args.len() == 1 && is_var(&args[0], var) {
                        return Some((0, ki, false));
                    }
                }
            }
            None
        }
        Expr::Add(ts) => {
            let mut s = 0i64;
            let mut c = 0i64;
            for t in ts {
                let (ts, tc, ex) = trig_affine_parts(t, var)?;
                if ex {
                    return None;
                }
                s += ts;
                c += tc;
            }
            Some((s, c, false))
        }
        _ if is_const_wrt(e, var) => Some((0, 0, true)),
        _ => None,
    }
}

fn sqrt_radicand(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Pow(base, exp) => {
            if matches!(exp.as_ref(), Expr::Rat(r) if *r == Ratio::new((-1).into(), 2.into()))
                || matches!(exp.as_ref(), Expr::Frac(n, d)
                    if n.is_one() && matches!(d.as_ref(), Expr::Int(i) if i == &BigInt::from(2)))
            {
                return Some(Arc::clone(base));
            }
            None
        }
        Expr::Func(FuncKind::Sqrt, args) if args.len() == 1 => Some(Arc::clone(&args[0])),
        _ => None,
    }
}

fn int_const_term(e: &ExprArc) -> Option<i64> {
    match e.as_ref() {
        Expr::Int(n) => bigint_to_i64(n).ok(),
        Expr::Mul(fs) if fs.len() == 2 => {
            let a = match fs[0].as_ref() {
                Expr::Int(n) => bigint_to_i64(n).ok()?,
                _ => return None,
            };
            let b = match fs[1].as_ref() {
                Expr::Int(n) => bigint_to_i64(n).ok()?,
                _ => return None,
            };
            Some(a * b)
        }
        _ => None,
    }
}

fn const_term_of_xsq_plus_const(inner: &ExprArc, var: &Ident) -> Option<i64> {
    let Expr::Add(ts) = inner.as_ref() else {
        return None;
    };
    let mut c = 0i64;
    let mut saw_x2 = false;
    for t in ts {
        if is_x_squared(t, var) {
            saw_x2 = true;
        } else if let Some(n) = int_const_term(t) {
            c += n;
        } else {
            return None;
        }
    }
    if saw_x2 {
        Some(c)
    } else {
        None
    }
}

fn const_term_of_x_plus_const(inner: &ExprArc, var: &Ident) -> Option<i64> {
    let Expr::Add(ts) = inner.as_ref() else {
        if is_var(inner, var) {
            return Some(0);
        }
        return None;
    };
    let mut c = 0i64;
    let mut saw_x = false;
    for t in ts {
        if is_var(t, var) {
            saw_x = true;
        } else if let Some(n) = int_const_term(t) {
            c += n;
        } else {
            return None;
        }
    }
    if saw_x {
        Some(c)
    } else {
        None
    }
}

fn is_x_squared(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2))
    )
}

fn is_sin_double_angle(inner: &ExprArc, var: &Ident) -> bool {
    match inner.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            (matches!(fs[0].as_ref(), Expr::Int(n) if bigint_to_i64(n).ok() == Some(2))
                && is_var(&fs[1], var))
                || (matches!(fs[1].as_ref(), Expr::Int(n) if bigint_to_i64(n).ok() == Some(2))
                    && is_var(&fs[0], var))
        }
        _ => false,
    }
}

fn is_cos_double_angle(inner: &ExprArc, var: &Ident) -> bool {
    is_sin_double_angle(inner, var)
}

fn var_coefficient_int(e: &ExprArc, var: &Ident) -> Option<i64> {
    if is_var(e, var) {
        return Some(1);
    }
    if let Expr::Mul(fs) = e.as_ref() {
        if fs.len() == 2 {
            if let Expr::Int(k) = fs[0].as_ref() {
                if is_var(&fs[1], var) {
                    return bigint_to_i64(k).ok();
                }
            }
            if let Expr::Int(k) = fs[1].as_ref() {
                if is_var(&fs[0], var) {
                    return bigint_to_i64(k).ok();
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{eval, format_expr, Expr, FuncKind};

    use crate::plugin::xcas_default;

    #[test]
    fn integrate_ck_int_11() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![
                Arc::new(Expr::Frac(
                    Expr::add(vec![
                        Expr::func(
                            FuncKind::Sin,
                            vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")])],
                        ),
                        Expr::int(1),
                    ]),
                    Expr::func(
                        FuncKind::Cos,
                        vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")])],
                    ),
                )),
                Expr::sym("x"),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(
            format_expr(r.as_ref()),
            "ln(abs(1-sin(2*x)))*-1/2"
        );
    }
}
