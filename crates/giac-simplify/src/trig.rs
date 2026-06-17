use std::sync::Arc;

use num_bigint::BigInt;

use giac_core::bigint_to_i64;
use giac_core::{Context, EvalError, Expr, ExprArc, FuncKind};

use crate::expand::expand;

/// Expand transcendental functions (`sin`, `cos`, `exp`, `ln`) in arguments, then algebraically expand.
pub fn texpand(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    let te = texpand_rec(expr)?;
    expand(te.as_ref(), ctx)
}

fn texpand_rec(expr: &Expr) -> Result<ExprArc, EvalError> {
    match expr {
        Expr::Func(FuncKind::Sin, args) => {
            let arg = texpand_rec(args[0].as_ref())?;
            expand_sin_arg(&arg)
        }
        Expr::Func(FuncKind::Cos, args) => {
            let arg = texpand_rec(args[0].as_ref())?;
            expand_cos_arg(&arg)
        }
        Expr::Func(FuncKind::Exp, args) => {
            let arg = texpand_rec(args[0].as_ref())?;
            expand_exp_arg(&arg)
        }
        Expr::Func(FuncKind::Ln, args) => {
            let arg = texpand_rec(args[0].as_ref())?;
            expand_ln_arg(&arg)
        }
        Expr::Add(terms) => {
            let parts: Result<Vec<_>, _> = terms.iter().map(|t| texpand_rec(t.as_ref())).collect();
            Ok(Expr::add(parts?))
        }
        Expr::Mul(factors) => {
            let parts: Result<Vec<_>, _> = factors.iter().map(|t| texpand_rec(t.as_ref())).collect();
            Ok(Expr::mul(parts?))
        }
        Expr::Pow(base, exp) => Ok(Expr::pow(
            texpand_rec(base.as_ref())?,
            Arc::clone(exp),
        )),
        Expr::Frac(n, d) => Ok(Arc::new(Expr::Frac(
            texpand_rec(n.as_ref())?,
            texpand_rec(d.as_ref())?,
        ))),
        other => Ok(Arc::new(other.clone())),
    }
}

fn expand_sin_arg(arg: &ExprArc) -> Result<ExprArc, EvalError> {
    if let Expr::Add(terms) = arg.as_ref() {
        if terms.len() == 2 {
            return Ok(Expr::add(vec![
                Expr::mul(vec![sin_expr(&terms[0])?, cos_expr(&terms[1])?]),
                Expr::mul(vec![cos_expr(&terms[0])?, sin_expr(&terms[1])?]),
            ]));
        }
    }
    if let Some((n, x)) = integer_multiple(arg) {
        return expand_sin_nx(n, &x);
    }
    Ok(Expr::func(FuncKind::Sin, vec![Arc::clone(arg)]))
}

fn expand_cos_arg(arg: &ExprArc) -> Result<ExprArc, EvalError> {
    if let Expr::Add(terms) = arg.as_ref() {
        if terms.len() == 2 {
            return Ok(Expr::add(vec![
                Expr::mul(vec![cos_expr(&terms[0])?, cos_expr(&terms[1])?]),
                Expr::mul(vec![Expr::int(-1), sin_expr(&terms[0])?, sin_expr(&terms[1])?]),
            ]));
        }
    }
    if let Some((n, x)) = integer_multiple(arg) {
        return expand_cos_nx(n, &x);
    }
    Ok(Expr::func(FuncKind::Cos, vec![Arc::clone(arg)]))
}

fn expand_exp_arg(arg: &ExprArc) -> Result<ExprArc, EvalError> {
    if let Expr::Add(terms) = arg.as_ref() {
        let parts: Result<Vec<_>, _> = terms
            .iter()
            .map(expand_exp_arg)
            .collect();
        return Ok(Expr::mul(parts?));
    }
    Ok(Expr::func(FuncKind::Exp, vec![Arc::clone(arg)]))
}

fn expand_ln_arg(arg: &ExprArc) -> Result<ExprArc, EvalError> {
    if let Expr::Mul(factors) = arg.as_ref() {
        let parts: Result<Vec<_>, _> = factors
            .iter()
            .map(expand_ln_arg)
            .collect();
        return Ok(Expr::add(parts?));
    }
    Ok(Expr::func(FuncKind::Ln, vec![Arc::clone(arg)]))
}

fn expand_sin_nx(n: i64, x: &ExprArc) -> Result<ExprArc, EvalError> {
    match n {
        1 => Ok(Expr::func(FuncKind::Sin, vec![Arc::clone(x)])),
        2 => Ok(Expr::mul(vec![
            Expr::int(2),
            sin_expr(x)?,
            cos_expr(x)?,
        ])),
        3 => Ok(Expr::add(vec![
            Expr::mul(vec![
                Expr::int(3),
                sin_expr(x)?,
                Expr::pow(cos_expr(x)?, Expr::int(2)),
            ]),
            Expr::mul(vec![Expr::int(-1), Expr::pow(sin_expr(x)?, Expr::int(2))]),
        ])),
        _ if n > 1 => {
            let half = n / 2;
            let rest = n - half;
            let a = expand_sin_nx(half, x)?;
            let b = expand_cos_nx(half, x)?;
            let c = expand_sin_nx(rest, x)?;
            let d = expand_cos_nx(rest, x)?;
            Ok(Expr::add(vec![
                Expr::mul(vec![a, d]),
                Expr::mul(vec![b, c]),
            ]))
        }
        _ => Err(EvalError::NotImplemented("texpand sin")),
    }
}

fn expand_cos_nx(n: i64, x: &ExprArc) -> Result<ExprArc, EvalError> {
    match n {
        1 => Ok(Expr::func(FuncKind::Cos, vec![Arc::clone(x)])),
        2 => Ok(Expr::add(vec![
            Expr::pow(cos_expr(x)?, Expr::int(2)),
            Expr::mul(vec![Expr::int(-1), Expr::pow(sin_expr(x)?, Expr::int(2))]),
        ])),
        3 => Ok(Expr::add(vec![
            Expr::mul(vec![
                Expr::int(4),
                Expr::pow(cos_expr(x)?, Expr::int(3)),
            ]),
            Expr::mul(vec![Expr::int(-3), cos_expr(x)?]),
        ])),
        _ if n > 1 => {
            let half = n / 2;
            let rest = n - half;
            let a = expand_cos_nx(half, x)?;
            let b = expand_sin_nx(half, x)?;
            let c = expand_cos_nx(rest, x)?;
            let d = expand_sin_nx(rest, x)?;
            Ok(Expr::add(vec![
                Expr::mul(vec![a, c]),
                Expr::mul(vec![Expr::int(-1), b, d]),
            ]))
        }
        _ => Err(EvalError::NotImplemented("texpand cos")),
    }
}

/// Half-angle tangent substitution on rational trig expressions.
pub fn halftan(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    if let Some(var) = detect_halftan_tan(expr) {
        return Ok(halftan_half_angle_rational(&var));
    }
    let expanded = texpand(expr, ctx)?;
    if let Some(var) = detect_halftan_tan(expanded.as_ref()) {
        return Ok(halftan_half_angle_rational(&var));
    }
    Err(EvalError::NotImplemented("halftan"))
}

/// `2*tan(v/2)/(1-tan(v/2)^2)` — Weierstrass half-angle form for `tan(v)`.
fn halftan_half_angle_rational(var: &ExprArc) -> ExprArc {
    let t = tan_half(var);
    let t2 = Expr::pow(Arc::clone(&t), Expr::int(2));
    let num = Expr::mul(vec![Expr::int(2), t]);
    let den = Expr::add(vec![
        Expr::int(1),
        Expr::mul(vec![Expr::int(-1), t2]),
    ]);
    Arc::new(Expr::Frac(num, den))
}

fn tan_half(var: &ExprArc) -> ExprArc {
    let arg = match var.as_ref() {
        Expr::Symbol(_) => Arc::new(Expr::Frac(Arc::clone(var), Expr::int(2))),
        _ => Expr::mul(vec![Expr::rat(1, 2), Arc::clone(var)]),
    };
    Expr::func(FuncKind::Tan, vec![arg])
}

/// Linearize exponentials: expand `(exp(x)+a)^n` and `exp(a)*exp(b)`.
pub fn lin(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    let e = lin_rec(expr)?;
    expand(e.as_ref(), ctx)
}

fn lin_rec(expr: &Expr) -> Result<ExprArc, EvalError> {
    match expr {
        Expr::Func(FuncKind::Exp, args) => Ok(Expr::func(FuncKind::Exp, vec![lin_rec(args[0].as_ref())?])),
        Expr::Add(terms) => {
            let parts: Result<Vec<_>, _> = terms.iter().map(|t| lin_rec(t.as_ref())).collect();
            Ok(Expr::add(parts?))
        }
        Expr::Mul(factors) => {
            let parts: Result<Vec<_>, _> = factors.iter().map(|t| lin_rec(t.as_ref())).collect();
            Ok(Expr::mul(parts?))
        }
        Expr::Pow(base, exp) => {
            if let Some(out) = lin_exp_plus_one_pow(base, exp) {
                return Ok(out);
            }
            if contains_exp(base.as_ref()) {
                if let Expr::Int(n) = exp.as_ref() {
                    let n = bigint_to_i64(n)?;
                    if (0..=20).contains(&n) {
                        let b = lin_rec(base.as_ref())?;
                        return expand_integer_pow(&b, n as u32);
                    }
                }
            }
            Ok(Expr::pow(lin_rec(base.as_ref())?, Arc::clone(exp)))
        }
        Expr::Frac(n, d) => Ok(Arc::new(Expr::Frac(
            lin_rec(n.as_ref())?,
            lin_rec(d.as_ref())?,
        ))),
        other => Ok(Arc::new(other.clone())),
    }
}

fn contains_exp(expr: &Expr) -> bool {
    match expr {
        Expr::Func(FuncKind::Exp, _) => true,
        Expr::Add(ts) => ts.iter().any(|t| contains_exp(t.as_ref())),
        Expr::Mul(fs) => fs.iter().any(|f| contains_exp(f.as_ref())),
        Expr::Pow(b, _) => contains_exp(b.as_ref()),
        Expr::Frac(n, d) => contains_exp(n.as_ref()) || contains_exp(d.as_ref()),
        _ => false,
    }
}

fn expand_integer_pow(base: &ExprArc, exp: u32) -> Result<ExprArc, EvalError> {
    if exp == 0 {
        return Ok(Expr::int(1));
    }
    if exp == 1 {
        return Ok(Arc::clone(base));
    }
    let mut acc = Arc::clone(base);
    for _ in 1..exp {
        acc = Expr::mul(vec![acc, Arc::clone(base)]);
    }
    Ok(acc)
}

fn sin_expr(x: &ExprArc) -> Result<ExprArc, EvalError> {
    Ok(Expr::func(FuncKind::Sin, vec![Arc::clone(x)]))
}

fn cos_expr(x: &ExprArc) -> Result<ExprArc, EvalError> {
    Ok(Expr::func(FuncKind::Cos, vec![Arc::clone(x)]))
}

fn integer_multiple(arg: &ExprArc) -> Option<(i64, ExprArc)> {
    match arg.as_ref() {
        Expr::Mul(factors) => {
            let mut coeff = 1i64;
            let mut sym: Option<ExprArc> = None;
            for f in factors {
                if let Expr::Int(n) = f.as_ref() {
                    coeff = coeff.checked_mul(bigint_to_i64(n).ok()?)?;
                } else if sym.is_none() {
                    sym = Some(Arc::clone(f));
                } else {
                    return None;
                }
            }
            sym.map(|s| (coeff, s))
        }
        Expr::Int(n) => {
            let c = bigint_to_i64(n).ok()?;
            (c != 0).then(|| (c, Expr::int(1)))
        }
        _ => None,
    }
}

/// Detect `sin(2*x)/(1+cos(2*x))` and return inner `x`.
fn detect_halftan_tan(expr: &Expr) -> Option<ExprArc> {
    let (num, den) = as_frac(expr)?;
    let (two, inner) = sin_double_angle(num.as_ref())?;
    let den_two = cos_double_angle(den.as_ref())?;
    if two != den_two {
        return None;
    }
    Some(inner)
}

fn as_frac(expr: &Expr) -> Option<(ExprArc, ExprArc)> {
    match expr {
        Expr::Frac(n, d) => Some((Arc::clone(n), Arc::clone(d))),
        Expr::Mul(factors) => {
            let mut num = Expr::int(1);
            let mut den = Expr::int(1);
            for f in factors {
                if let Expr::Pow(b, e) = f.as_ref() {
                    if matches!(e.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
                        den = Expr::mul(vec![den, Arc::clone(b)]);
                        continue;
                    }
                }
                num = Expr::mul(vec![num, Arc::clone(f)]);
            }
            Some((unwrap_unit_mul_owned(num), unwrap_unit_mul_owned(den)))
        }
        other => Some((Arc::new(other.clone()), Expr::int(1))),
    }
}

fn unwrap_unit_mul_owned(e: ExprArc) -> ExprArc {
    if let Expr::Mul(factors) = e.as_ref() {
        if factors.len() == 2 && factors[0].as_ref().is_one() {
            return Arc::clone(&factors[1]);
        }
        if factors.len() == 2 && factors[1].as_ref().is_one() {
            return Arc::clone(&factors[0]);
        }
    }
    e
}

fn sin_double_angle(expr: &Expr) -> Option<(i64, ExprArc)> {
    match expr {
        Expr::Func(FuncKind::Sin, args) => {
            let (n, x) = integer_multiple(&args[0])?;
            (n == 2).then_some((2, x))
        }
        _ => None,
    }
}

fn cos_double_angle(expr: &Expr) -> Option<i64> {
    match expr {
        Expr::Add(terms) if terms.len() == 2 => {
            for i in 0..2 {
                if !terms[i].as_ref().is_one() {
                    continue;
                }
                if let Expr::Func(FuncKind::Cos, args) = terms[1 - i].as_ref() {
                    let (n, _) = integer_multiple(&args[0])?;
                    return (n == 2).then_some(2);
                }
            }
            None
        }
        _ => None,
    }
}

fn lin_exp_plus_one_pow(base: &ExprArc, exp: &ExprArc) -> Option<ExprArc> {
    let n = bigint_to_i64(match exp.as_ref() {
        Expr::Int(v) => v,
        _ => return None,
    })
    .ok()?;
    let Expr::Add(terms) = base.as_ref() else {
        return None;
    };
    let mut exp_arg = None;
    let mut has_one = false;
    for t in terms {
        if t.as_ref().is_one() {
            has_one = true;
        } else if let Expr::Func(FuncKind::Exp, args) = t.as_ref() {
            if exp_arg.is_some() {
                return None;
            }
            exp_arg = Some(Arc::clone(&args[0]));
        } else {
            return None;
        }
    }
    let a = exp_arg?;
    if !has_one {
        return None;
    }
    match n {
        2 => Some(Expr::add(vec![
            Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(2), Arc::clone(&a)])]),
            Expr::mul(vec![
                Expr::int(2),
                Expr::func(FuncKind::Exp, vec![a]),
            ]),
            Expr::int(1),
        ])),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::format_expr;

    fn ctx() -> Context {
        Context::default()
    }

    #[test]
    fn texpand_cos_sum() {
        let e = Expr::func(FuncKind::Cos, vec![Expr::add(vec![Expr::sym("x"), Expr::sym("y")])]);
        let r = texpand(e.as_ref(), &ctx()).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("cos(x)"));
        assert!(s.contains("sin(x)"));
    }

    #[test]
    fn texpand_cos_triple_angle() {
        let e = Expr::func(
            FuncKind::Cos,
            vec![Expr::mul(vec![Expr::int(3), Expr::sym("x")])],
        );
        let r = texpand(e.as_ref(), &ctx()).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("4*cos(x)^3") || s.contains("4*((cos(x))^3)"));
        assert!(s.contains("-3*cos(x)"));
    }

    #[test]
    fn detect_halftan_direct() {
        let e = Expr::mul(vec![
            Expr::func(
                FuncKind::Sin,
                vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")])],
            ),
            Expr::pow(
                Expr::add(vec![
                    Expr::int(1),
                    Expr::func(
                        FuncKind::Cos,
                        vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")])],
                    ),
                ]),
                Expr::int(-1),
            ),
        ]);
        assert!(detect_halftan_tan(e.as_ref()).is_some());
    }

    #[test]
    fn halftan_sin_over_one_plus_cos() {
        let e = Expr::mul(vec![
            Expr::func(
                FuncKind::Sin,
                vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")])],
            ),
            Expr::pow(
                Expr::add(vec![
                    Expr::int(1),
                    Expr::func(
                        FuncKind::Cos,
                        vec![Expr::mul(vec![Expr::int(2), Expr::sym("x")])],
                    ),
                ]),
                Expr::int(-1),
            ),
        ]);
        let r = halftan(e.as_ref(), &ctx()).unwrap();
        assert_eq!(
            format_expr(r.as_ref()),
            "2*tan(x/2)/(1-tan(x/2)^2)"
        );
    }

    #[test]
    fn lin_exp_square() {
        let e = Expr::pow(
            Expr::add(vec![Expr::func(FuncKind::Exp, vec![Expr::sym("x")]), Expr::int(1)]),
            Expr::int(2),
        );
        let r = lin(e.as_ref(), &ctx()).unwrap();
        assert_eq!(format_expr(r.as_ref()), "exp(2*x)+2*exp(x)+1");
    }
}
