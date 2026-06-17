use std::sync::Arc;

use giac_core::{
    bigint_to_i64, EvalError, Expr, ExprArc, FuncKind, Ident,
};

/// Symbolic differentiation (GIAC-113 / `giac-calculus`).
pub fn diff(expr: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    match expr.as_ref() {
        Expr::Add(terms) => {
            let parts: Result<Vec<_>, _> = terms.iter().map(|t| diff(t, var)).collect();
            Ok(Expr::add(parts?))
        }
        Expr::Mul(factors) => diff_mul(factors, var),
        Expr::Pow(base, exp) => diff_pow(base, exp, var),
        Expr::Frac(num, den) => diff_quotient(num, den, var),
        Expr::Func(FuncKind::Sin, args) => diff_sin(&args[0], var),
        Expr::Func(FuncKind::Cos, args) => diff_cos(&args[0], var),
        Expr::Func(FuncKind::Ln, args) => diff_ln(&args[0], var),
        Expr::Func(FuncKind::Exp, args) => diff_exp(&args[0], var),
        Expr::Func(FuncKind::Tan, args) => diff_tan(&args[0], var),
        Expr::Func(FuncKind::Atan, args) => diff_atan(&args[0], var),
        Expr::Symbol(id) if id == var => Ok(Expr::int(1)),
        _ if is_const_wrt(expr, var) => Ok(Expr::int(0)),
        _ => Err(EvalError::NotImplemented("diff")),
    }
}

fn diff_mul(factors: &[ExprArc], var: &Ident) -> Result<ExprArc, EvalError> {
    if factors.is_empty() {
        return Ok(Expr::int(0));
    }
    if factors.len() == 1 {
        return diff(&factors[0], var);
    }
    let mut terms = Vec::with_capacity(factors.len());
    for (i, f) in factors.iter().enumerate() {
        let df = diff(f, var)?;
        let others: Vec<ExprArc> = factors
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, g)| Arc::clone(g))
            .collect();
        let rest = if others.is_empty() {
            Expr::int(1)
        } else if others.len() == 1 {
            Arc::clone(&others[0])
        } else {
            Expr::mul(others)
        };
        terms.push(Expr::mul(vec![df, rest]));
    }
    Ok(Expr::add(terms))
}

fn diff_pow(base: &ExprArc, exp: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    if is_var(base, var) {
        if let Expr::Int(n) = exp.as_ref() {
            let n = bigint_to_i64(n)?;
            return match n {
                0 => Ok(Expr::int(0)),
                1 => Ok(Expr::int(1)),
                _ => Ok(Expr::mul(vec![
                    Expr::int(n),
                    Expr::pow(Arc::clone(base), Expr::int(n - 1)),
                ])),
            };
        }
    }
    if is_const_wrt(exp, var) {
        if let Expr::Int(n) = exp.as_ref() {
            let n = bigint_to_i64(n)?;
            if n == 0 {
                return Ok(Expr::int(0));
            }
            let fp = diff(base, var)?;
            return Ok(Expr::mul(vec![
                Expr::int(n),
                Expr::pow(Arc::clone(base), Expr::int(n - 1)),
                fp,
            ]));
        }
    }
    Err(EvalError::NotImplemented("diff pow"))
}

fn diff_quotient(num: &ExprArc, den: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    let np = diff(num, var)?;
    let dp = diff(den, var)?;
    let num_expr = Expr::add(vec![
        Expr::mul(vec![np, Arc::clone(den)]),
        Expr::mul(vec![Expr::int(-1), Arc::clone(num), dp]),
    ]);
    Ok(Arc::new(Expr::Frac(
        num_expr,
        Expr::pow(Arc::clone(den), Expr::int(2)),
    )))
}

fn diff_sin(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![
        Expr::func(FuncKind::Cos, vec![Arc::clone(arg)]),
        diff(arg, var)?,
    ]))
}

fn diff_cos(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![
        Expr::int(-1),
        Expr::func(FuncKind::Sin, vec![Arc::clone(arg)]),
        diff(arg, var)?,
    ]))
}

fn diff_ln(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Arc::new(Expr::Frac(diff(arg, var)?, Arc::clone(arg))))
}

fn diff_exp(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![
        Expr::func(FuncKind::Exp, vec![Arc::clone(arg)]),
        diff(arg, var)?,
    ]))
}

fn diff_atan(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![
        diff(arg, var)?,
        Expr::pow(
            Expr::add(vec![Expr::int(1), Expr::pow(Arc::clone(arg), Expr::int(2))]),
            Expr::int(-1),
        ),
    ]))
}

fn diff_tan(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    let cos = Expr::func(FuncKind::Cos, vec![Arc::clone(arg)]);
    Ok(Expr::mul(vec![
        diff(arg, var)?,
        Expr::pow(cos, Expr::int(-2)),
    ]))
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_const_wrt(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id != var,
        Expr::Int(_) | Expr::Rat(_) => true,
        Expr::Add(ts) => ts.iter().all(|t| is_const_wrt(t, var)),
        Expr::Mul(fs) => fs.iter().all(|f| is_const_wrt(f, var)),
        Expr::Pow(b, exp) => is_const_wrt(b, var) && is_const_wrt(exp, var),
        Expr::Func(_, args) => args.iter().all(|a| is_const_wrt(a, var)),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::{format_expr, simplify, Context};
    use giac_simplify::assert_equiv;

    fn x() -> Ident {
        Ident::new("x")
    }

    fn diff_simplified(e: ExprArc) -> ExprArc {
        let ctx = Context::default();
        simplify(diff(&e, &x()).unwrap().as_ref(), &ctx).unwrap()
    }

    #[test]
    fn diff_x_squared() {
        let e = Expr::pow(Expr::sym("x"), Expr::int(2));
        let r = diff_simplified(e);
        let expected = Expr::mul(vec![Expr::int(2), Expr::sym("x")]);
        let ctx = Context::default();
        assert!(assert_equiv(r.as_ref(), expected.as_ref(), &ctx).unwrap());
    }

    #[test]
    fn diff_sin_x_squared() {
        let e = Expr::func(FuncKind::Sin, vec![Expr::pow(Expr::sym("x"), Expr::int(2))]);
        let r = diff_simplified(e);
        let s = format_expr(r.as_ref());
        assert!(s.contains("cos(x^2)"));
        assert!(s.contains("2") && s.contains("x"));
    }

    #[test]
    fn diff_ln_times_x_squared() {
        let e = Expr::mul(vec![
            Expr::func(FuncKind::Ln, vec![Expr::sym("x")]),
            Expr::pow(Expr::sym("x"), Expr::int(2)),
        ]);
        let r = diff_simplified(e);
        let s = format_expr(r.as_ref());
        assert!(s.contains("ln(x)"));
        assert!(s.contains("x^2") || s.contains("x^2"));
    }

    #[test]
    fn diff_constant_is_zero() {
        let r = diff_simplified(Expr::int(5));
        assert_eq!(format_expr(r.as_ref()), "0");
    }

    #[test]
    fn diff_atan_x() {
        let e = Expr::func(FuncKind::Atan, vec![Expr::sym("x")]);
        let r = diff(&e, &x());
        assert!(r.is_ok(), "{:?}", r);
    }
}
