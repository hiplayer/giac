use std::sync::Arc;

use num_bigint::BigInt;

use giac_core::{
    bigint_to_i64, expand, Context, EvalError, Expr, ExprArc, FuncKind, Ident,
};

/// Basic integration rules (Phase 1 / GIAC-110 subset).
///
/// ## Supported
///
/// - Constants, `x`, and `x^n` for integer `n ≠ -1`
/// - Sums and constant multiples (`integrate(c*f) = c*integrate(f)`)
/// - `1/x`, `1/(x^2+1)`, `x/(x^2+1)`, `1/(1-x^2)`, `1/(4+x^2)` (partial)
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
            let parts: Result<Vec<_>, _> = terms.iter().map(|t| integrate(t, var)).collect();
            Ok(Expr::add(parts?))
        }
        Expr::Func(kind, args) => integrate_func(*kind, args, var),
        _ if is_const_wrt(expr, var) => Ok(Expr::mul(vec![Arc::clone(expr), var_to_expr(var)])),
        _ => Err(EvalError::NotImplemented("integrate")),
    }
}

fn integrate_frac(num: &ExprArc, den: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    if is_const_wrt(num, var) {
        let inner = integrate_reciprocal(den, var)?;
        if num.is_one() {
            return Ok(inner);
        }
        return Ok(Expr::mul(vec![Arc::clone(num), inner]));
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
    if let Expr::Add(terms) = den.as_ref() {
        return integrate_reciprocal_quadratic(terms, var);
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
            if linear_coefficient(arg, var).is_some() {
                let k = linear_coefficient(arg, var).unwrap();
                return Ok(Expr::mul(vec![
                    Expr::int(-1),
                    Expr::func(FuncKind::Cos, vec![Arc::clone(arg)]),
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
        }
    }
    None
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
        if let Expr::Add(terms) = base.as_ref() {
            return integrate_reciprocal_quadratic(terms, var);
        }
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

fn integrate_reciprocal_quadratic(terms: &[ExprArc], var: &Ident) -> Result<ExprArc, EvalError> {
    let x = var_to_expr(var);
    if terms.len() == 2 {
        if terms.iter().any(is_one) && terms.iter().any(|t| is_neg_x_power(t, var, 4)) {
            return Ok(Expr::add(vec![
                Expr::mul(vec![
                    Expr::rat(-1, 4),
                    ln_abs_expr(Expr::add(vec![x.clone(), Expr::int(-1)])),
                ]),
                Expr::mul(vec![
                    Expr::rat(1, 4),
                    ln_abs_expr(Expr::add(vec![x.clone(), Expr::int(1)])),
                ]),
                Expr::mul(vec![Expr::rat(1, 2), Expr::func(FuncKind::Atan, vec![x])]),
            ]));
        }
        if terms.iter().any(is_one) && terms.iter().any(|t| is_x_squared(t, var)) {
            return Ok(Expr::func(FuncKind::Atan, vec![x]));
        }
        if terms.iter().any(|t| matches!(t.as_ref(), Expr::Int(n) if n == &num_bigint::BigInt::from(4)))
            && terms.iter().any(|t| is_x_squared(t, var))
        {
            return Ok(Expr::mul(vec![
                Expr::rat(1, 2),
                Expr::func(FuncKind::Atan, vec![Expr::mul(vec![x, Expr::rat(1, 2)])]),
            ]));
        }
    }
    Err(EvalError::NotImplemented("integrate quadratic"))
}

fn is_neg_x_power(e: &ExprArc, var: &Ident, pow: i64) -> bool {
    match e.as_ref() {
        Expr::Mul(factors) if factors.len() == 2 => {
            matches!(factors[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                && matches!(
                    factors[1].as_ref(),
                    Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(pow))
                )
        }
        Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(pow)) => {
            true
        }
        _ => false,
    }
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

fn ln_abs_expr(arg: ExprArc) -> ExprArc {
    Expr::func(FuncKind::Ln, vec![Expr::func(FuncKind::Abs, vec![arg])])
}

fn var_to_expr(var: &Ident) -> ExprArc {
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

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_one(e: &ExprArc) -> bool {
    e.is_one()
}

fn is_const_wrt(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id != var,
        Expr::Int(_) | Expr::Rat(_) => true,
        Expr::Add(ts) => ts.iter().all(|t| is_const_wrt(t, var)),
        Expr::Mul(fs) => fs.iter().all(|f| is_const_wrt(f, var)),
        Expr::Pow(b, _) => is_const_wrt(b, var),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::{eval, format_expr, Context, FuncKind};

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
        assert!(s.contains("atan(x)"));
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
        assert_eq!(format_expr(r.as_ref()), "atan(x)");
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
        assert_eq!(format_expr(r.as_ref()), "atan(x)");
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
        assert_eq!(format_expr(r.as_ref()), "1/2*atan(x*1/2)");
    }

    #[test]
    fn integrate_const_times_reciprocal() {
        let x = Ident::new("x");
        let e = Expr::mul(vec![Expr::int(3), Expr::pow(Expr::sym("x"), Expr::int(-1))]);
        let r = integrate(&e, &x).unwrap();
        assert_eq!(format_expr(r.as_ref()), "ln(abs(x))*3");
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
    fn integrate_x_squared_reciprocal() {
        let x = Ident::new("x");
        let den = Expr::pow(Expr::sym("x"), Expr::int(2));
        let e = Expr::pow(den, Expr::int(-1));
        assert!(matches!(
            integrate(&e, &x),
            Err(EvalError::NotImplemented(_))
        ));
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
        assert_eq!(format_expr(r.as_ref()), "tan(x)*1");
    }

    #[test]
    fn integrate_not_implemented_messages() {
        let x = Ident::new("x");
        let frac = Arc::new(Expr::Frac(
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::sym("x"),
        ));
        assert!(matches!(
            integrate(&frac, &x),
            Err(EvalError::NotImplemented("integrate frac"))
        ));
        let e = Expr::pow(
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(2)]),
            Expr::int(-1),
        );
        assert!(matches!(
            integrate(&e, &x),
            Err(EvalError::NotImplemented("integrate quadratic"))
        ));
    }
}
