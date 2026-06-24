//! Basic symbolic integration (`integrate`).
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../.doc/giac-calculus-api-stability.md) §4.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `integrate`, `integrate_frac`, `try_as_rational`, `ln_abs_expr`, `var_to_expr`, `is_var`, `is_const_wrt`, `is_exp_of_var` |
//! | **Partial** | `try_integrate_*` 启发式规则 |
//! | **Pipeline private** | `integrate_*`, `is_*`, 分部 / 换元辅助 |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::sync::Arc;

use num_bigint::BigInt;


use giac_core::{
    bigint_to_i64, Context, EvalError, Expr, ExprArc, FuncKind, Ident,
};
use giac_simplify::expand;

use crate::expr_util::{is_cos_of_var, is_ln_of_var, is_sin_of_var};
pub(crate) use crate::expr_util::{is_const_wrt, is_var, var_to_expr};
pub(crate) use crate::integrate_helpers::{
    affine_var_coeff, is_exp_of_var, is_x_squared_plus_const, linear_coefficient, ln_abs_expr,
    var_coefficient,
};
use crate::integrate_try_rules::{
    try_integrate_exp_over_linear_exp, try_integrate_exp_over_one_plus_exp2,
    try_integrate_exp_trig, try_integrate_sin2x_cos, try_integrate_sin_over_cos_squared,
    try_integrate_sin_over_cos_sq_frac, try_integrate_tanh_exp_form,
    try_integrate_tan_plus_tan_cubed, try_integrate_var_over_quadratic_squared,
    try_integrate_var_shifted_sqrt,
};

/// **Stable** — basic integration rules (Phase 1 / GIAC-110 subset).
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
/// **Stable** — symbolic integration
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

/// **Stable** — normalize `Expr` to `(num, den)` rational form.
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

/// **Stable** — integrate rational `num/den` w.r.t. `var`.
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

// **Pipeline private** — integrate reciprocal.
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

// **Pipeline private** — integrate func.
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




// **Pipeline private** — integrate sin squared.
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

// **Pipeline private** — integrate cos squared.
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

// **Pipeline private** — integrate sin cos product.
fn integrate_sin_cos_product(var: &Ident) -> ExprArc {
    let x = var_to_expr(var);
    Expr::mul(vec![
        Expr::rat(1, 2),
        Expr::pow(Expr::func(FuncKind::Sin, vec![x]), Expr::int(2)),
    ])
}

// **Pipeline private** — integrate x ln.
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



// **Pipeline private** — integrate mul.
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

// **Pipeline private** — count var factors.
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

// **Pipeline private** — integrate pow.
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


// **Pipeline private** — ln abs.
fn ln_abs(var: &Ident) -> ExprArc {
    ln_abs_expr(var_to_expr(var))
}

#[cfg(test)]
mod tests {
    //! Test tiers: **A** / **B** / **B+C** / **C** — `.doc/test-writing-spec.md`
    //! Audit: `.doc/issues/GIAC-expr-api-test-audit.md` §5
    //! Note: rule-class tests pair **B** `integrate()` with **A** `eval(Integrate)` where `assert_equiv` applies.

    use std::sync::Arc;

    use super::*;
    use giac_core::{eval, format_expr, Context, FuncKind};
    use giac_simplify::assert_equiv;

    fn xcas() -> Context {
        crate::plugin::xcas_default()
    }

    fn assert_integrate_matches(integrand: &Expr, var: &Ident, expected: &Expr) {
        let ctx = Context::default();
        let integrand_arc = Arc::new(integrand.clone());
        let r = integrate(&integrand_arc, var).expect("integrate");
        assert!(
            assert_equiv(r.as_ref(), expected, &ctx).expect("assert_equiv"),
            "integrate mismatch: got {}",
            format_expr(r.as_ref())
        );
    }

    // **A** — eval(Integrate) indefinite; same semantics as **B** `assert_integrate_matches`.
    fn eval_integrate_matches(integrand: &Expr, var: &Ident, expected: &Expr) {
        let ctx = xcas();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![Arc::new(integrand.clone()), Expr::sym(var.as_str())],
        );
        let r = eval(e.as_ref(), &ctx).expect("eval(Integrate)");
        assert!(
            assert_equiv(r.as_ref(), expected, &ctx).expect("assert_equiv"),
            "eval(Integrate) mismatch: got {}",
            format_expr(r.as_ref())
        );
    }

    // **B** — smoke: integrate does not fail on tanh-like integrand.
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

    // **B** — smoke: exp-over-linear integrand returns Ok.
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

    // **B** + **A** — integrate(1/x); assert_equiv vs ln(abs(x)).
    #[test]
    fn integrate_reciprocal() {
        let x = Ident::new("x");
        let e = Expr::pow(Expr::sym("x"), Expr::int(-1));
        let expected_ln = ln_abs_expr(Expr::sym("x"));
        let expected = expected_ln.as_ref();
        assert_integrate_matches(&e, &x, expected);
        eval_integrate_matches(&e, &x, expected);
    }

    // **B** — partfrac smoke; full antiderivative shape not pinned (partfrac/heuristic).
    #[test]
    fn integrate_one_minus_x_fourth() {
        let x = Ident::new("x");
        let den = Expr::add(vec![
            Expr::int(1),
            Expr::mul(vec![Expr::int(-1), Expr::pow(Expr::sym("x"), Expr::int(4))]),
        ]);
        let e = Arc::new(Expr::pow(den, Expr::int(-1)));
        assert!(integrate(&e, &x).is_ok());
    }

    // **B+C** — constant **B** display; sum uses display golden for ln+poly term.
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
        assert_eq!(format_expr(r.as_ref()), "ln(abs(x))+1*x");
    }

    // **B+C** — direct integrate(); display golden term order.
    #[test]
    fn integrate_const_times_x() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![Expr::int(3), Expr::sym("x")]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(1/2*x^2)*3");
    }

    // **B+C** — power rule; display golden + **A** eval(Integrate).
    #[test]
    fn integrate_x_and_x_squared() {
        let x = Ident::new("x");
        let r = integrate(&Expr::sym("x"), &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/2*x^2");
        eval_integrate_matches(&Expr::sym("x"), &x, r.as_ref());

        let e2 = Expr::pow(Expr::sym("x"), Expr::int(2));
        let r2 = integrate(&e2, &x).unwrap();
        assert_eq!(format_expr(r2.as_ref()), "1/3*x^3");
        eval_integrate_matches(&e2, &x, r2.as_ref());
    }

    // **B+C** — ∫1 dx display + **A**.
    #[test]
    fn integrate_one() {
        let x = Ident::new("x");
        let r = integrate(&Expr::int(1), &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1*x");
        eval_integrate_matches(&Expr::int(1), &x, r.as_ref());
    }

    // **B+C** — arctan integral; display golden + **A** assert_equiv.
    #[test]
    fn integrate_one_over_one_plus_x_squared() {
        let x = Ident::new("x");
        let den = Expr::add(vec![Expr::int(1), Expr::pow(Expr::sym("x"), Expr::int(2))]);
        let e = Expr::pow(den, Expr::int(-1));
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2*2^-1*atan(2^-1*(2*x+0))");
        eval_integrate_matches(&e, &x, r.as_ref());
    }

    // **B+C** — ∫(2/x); display golden.
    #[test]
    fn integrate_frac_with_constant_numerator() {
        let x = Ident::new("x");
        let e = Arc::new(Expr::Frac(Expr::int(2), Expr::sym("x")));
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2*ln(abs(x))");
    }

    // **B+C** — arctan; display term-order snapshot.
    #[test]
    fn integrate_one_over_x_squared_plus_one_term_order() {
        let x = Ident::new("x");
        let den = Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(1)]);
        let e = Expr::pow(den, Expr::int(-1));
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2*2^-1*atan(2^-1*(2*x+0))");
    }

    // **B** — ∫x/(x²+1); assert_equiv (heuristic / partfrac path).
    #[test]
    fn integrate_x_over_x_squared_plus_one() {
        let x = Ident::new("x");
        let den = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::int(1),
        ]);
        let e = Arc::new(Expr::Frac(Expr::sym("x"), Arc::clone(&den)));
        let expected = Expr::mul(vec![Expr::rat(1, 2), ln_abs_expr(den)]);
        assert_integrate_matches(e.as_ref(), &x, expected.as_ref());
        eval_integrate_matches(e.as_ref(), &x, expected.as_ref());
    }

    // **B** — unsupported integrand returns NotImplemented.
    #[test]
    fn integrate_unsupported_returns_not_implemented() {
        let x = Ident::new("x");
        let e = Expr::func(FuncKind::Atan, vec![Expr::sym("x")]);
        assert!(matches!(
            integrate(&e, &x),
            Err(EvalError::NotImplemented(_))
        ));
    }

    // **B+C** — ∫1/(4+x²); display golden (atan argument canonical form).
    #[test]
    fn integrate_const_over_quadratic() {
        let x = Ident::new("x");
        let den = Expr::add(vec![Expr::int(4), Expr::pow(Expr::sym("x"), Expr::int(2))]);
        let e = Expr::pow(den, Expr::int(-1));
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2*4^-1*atan(4^-1*(2*x+0))");
    }

    // **B+C** — ∫3/x display golden.
    #[test]
    fn integrate_const_times_reciprocal() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![Expr::int(3), Expr::pow(Expr::sym("x"), Expr::int(-1))]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "3*ln(abs(x))");
    }

    // **B+C** — constant product integrand.
    #[test]
    fn integrate_product_of_consts_only() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![Expr::int(2), Expr::int(3)]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(2*3)*x");
    }

    // **A** — eval(Integrate, bounds); definite integral via full plugin path.
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

    // **B** — ∫(1+x)²; assert_equiv vs expanded antiderivative.
    #[test]
    fn integrate_product_one_plus_x_squared() {
        let x = Ident::new("x");
        let factor = Expr::add(vec![Expr::int(1), Expr::sym("x")]);
        let e = Expr::mul(vec![factor.clone(), factor]);
        let expected = Expr::add(vec![
            Expr::sym("x"),
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::rat(1, 3), Expr::pow(Expr::sym("x"), Expr::int(3))]),
        ]);
        assert_integrate_matches(&e, &x, expected.as_ref());
    }

    // **B** — ∫x/(x²+1)² smoke (partfrac/heuristic); no display shape pin.
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
        assert!(integrate(&e, &x).is_ok());
    }

    // **B** — ∫1/x²; assert_equiv vs -1/x.
    #[test]
    fn integrate_x_squared_reciprocal() {
        let x = Ident::new("x");
        let den = Expr::pow(Expr::sym("x"), Expr::int(2));
        let e = Expr::pow(den, Expr::int(-1));
        let expected = Expr::mul(vec![Expr::int(-1), Expr::pow(Expr::sym("x"), Expr::int(-1))]);
        assert_integrate_matches(&e, &x, expected.as_ref());
        eval_integrate_matches(&e, &x, expected.as_ref());
    }

    // **B+C** — sin·cos → sin²/2 display.
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

    // **B+C** — sec² integrand → tan(x) display.
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

    // **B** — Ok vs NotImplemented paths on selected integrands.
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

    // **B+C** — sin/cos primitives; display golden with 1^-1 factors.
    #[test]
    fn integrate_cos_and_sin() {
        let x = Ident::new("x");
        let cos_r = integrate(
            &Expr::func(FuncKind::Cos, vec![Expr::sym("x")]),
            &x,
        )
        .unwrap();
        assert_eq!(format_expr(cos_r.as_ref()), "sin(x)*1^-1");
        let sin_r = integrate(
            &Expr::func(FuncKind::Sin, vec![Expr::sym("x")]),
            &x,
        )
        .unwrap();
        assert_eq!(format_expr(sin_r.as_ref()), "-1*cos(x)*1^-1");
    }

    // **B+C** — ∫x³ display.
    #[test]
    fn integrate_x_cubed() {
        let x = Ident::new("x");
        let e = Expr::pow(Expr::sym("x"), Expr::int(3));
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/4*x^4");
    }
}
