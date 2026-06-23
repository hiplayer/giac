//! Symbolic differentiation (`diff`).
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../.doc/giac-calculus-api-stability.md) §7.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `diff` |
//! | **Pipeline private** | `diff_*`, `is_var`, `is_const_wrt` |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{
    bigint_to_i64, EvalError, Expr, ExprArc, FuncKind, Ident,
};

use crate::expr_util::is_var;

/// **Stable** — symbolic differentiation (GIAC-113 / `giac-calculus`).
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

// **Pipeline private** — product rule for `Mul`.
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

// **Pipeline private** — power rule (integer exponent cases).
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

// **Pipeline private** — quotient rule via `Frac`.
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

// **Pipeline private** — chain rule for `sin`.
fn diff_sin(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![
        Expr::func(FuncKind::Cos, vec![Arc::clone(arg)]),
        diff(arg, var)?,
    ]))
}

// **Pipeline private** — chain rule for `cos`.
fn diff_cos(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![
        Expr::int(-1),
        Expr::func(FuncKind::Sin, vec![Arc::clone(arg)]),
        diff(arg, var)?,
    ]))
}

// **Pipeline private** — chain rule for `ln`.
fn diff_ln(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Arc::new(Expr::Frac(diff(arg, var)?, Arc::clone(arg))))
}

// **Pipeline private** — chain rule for `exp`.
fn diff_exp(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![
        Expr::func(FuncKind::Exp, vec![Arc::clone(arg)]),
        diff(arg, var)?,
    ]))
}

// **Pipeline private** — chain rule for `atan`.
fn diff_atan(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![
        diff(arg, var)?,
        Expr::pow(
            Expr::add(vec![Expr::int(1), Expr::pow(Arc::clone(arg), Expr::int(2))]),
            Expr::int(-1),
        ),
    ]))
}

// **Pipeline private** — chain rule for `tan`.
fn diff_tan(arg: &ExprArc, var: &Ident) -> Result<ExprArc, EvalError> {
    let cos = Expr::func(FuncKind::Cos, vec![Arc::clone(arg)]);
    Ok(Expr::mul(vec![
        diff(arg, var)?,
        Expr::pow(cos, Expr::int(-2)),
    ]))
}

// **Pipeline private** — syntactic constness w.r.t. `var` (local copy).
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
    //! Test tiers — `.doc/test-writing-spec.md` · audit §5b (`diff.rs`)

    use std::sync::Arc;

    use giac_core::{eval, format_expr, simplify, Context, Expr, FuncKind};
    use giac_simplify::assert_equiv;

    use super::*;

    fn x() -> Ident {
        Ident::new("x")
    }

    fn xcas() -> Context {
        crate::plugin::xcas_default()
    }

    fn diff_simplified(e: ExprArc) -> ExprArc {
        let ctx = xcas();
        simplify(diff(&e, &x()).unwrap().as_ref(), &ctx).unwrap()
    }

    fn diff_matches(e: &Expr, expected: &Expr) {
        let ctx = xcas();
        let r = diff_simplified(Arc::new(e.clone()));
        assert!(
            assert_equiv(r.as_ref(), expected, &ctx).expect("assert_equiv"),
            "diff mismatch: got {}",
            format_expr(r.as_ref())
        );
    }

    fn eval_diff_matches(e: &Expr, expected: &Expr) {
        let ctx = xcas();
        let call = Expr::func(
            FuncKind::Diff,
            vec![Arc::new(e.clone()), Expr::sym("x")],
        );
        let r = eval(call.as_ref(), &ctx).unwrap();
        assert!(
            assert_equiv(r.as_ref(), expected, &ctx).expect("assert_equiv"),
            "eval(Diff) mismatch: got {}",
            format_expr(r.as_ref())
        );
    }

    // **B** + **A** — power rule.
    #[test]
    fn diff_x_squared() {
        let e = Expr::pow(Expr::sym("x"), Expr::int(2));
        let expected = Expr::mul(vec![Expr::int(2), Expr::sym("x")]);
        diff_matches(&e, expected.as_ref());
        eval_diff_matches(&e, expected.as_ref());
    }

    // **B** + **A** — chain rule sin(x²); factor order matches `diff` output.
    #[test]
    fn diff_sin_x_squared() {
        let e = Expr::func(FuncKind::Sin, vec![Expr::pow(Expr::sym("x"), Expr::int(2))]);
        let expected = Expr::mul(vec![
            Expr::int(2),
            Expr::func(FuncKind::Cos, vec![Expr::pow(Expr::sym("x"), Expr::int(2))]),
            Expr::sym("x"),
        ]);
        diff_matches(&e, expected.as_ref());
        eval_diff_matches(&e, expected.as_ref());
    }

    // **B** + **A** — product rule ln(x)·x²; C display + eval agrees with `diff`.
    #[test]
    fn diff_ln_times_x_squared() {
        let e = Expr::mul(vec![
            Expr::func(FuncKind::Ln, vec![Expr::sym("x")]),
            Expr::pow(Expr::sym("x"), Expr::int(2)),
        ]);
        let r = diff_simplified(Arc::clone(&e));
        let s = format_expr(r.as_ref());
        assert_eq!(s, "(1)/(x)*x^2+(2*x^1)*ln(x)");
        eval_diff_matches(e.as_ref(), r.as_ref());
    }

    // **B** + **A** — constant → 0.
    #[test]
    fn diff_constant_is_zero() {
        let expected = Expr::int(0);
        diff_matches(&Expr::int(5), expected.as_ref());
        eval_diff_matches(&Expr::int(5), expected.as_ref());
    }

    // **B** — atan(x) supported path smoke.
    #[test]
    fn diff_atan_x() {
        let e = Expr::func(FuncKind::Atan, vec![Expr::sym("x")]);
        assert!(diff(&e, &x()).is_ok());
    }

    // **B** + **A** — d/dx(1/x); display golden on quotient form from `diff`.
    #[test]
    fn diff_frac_one_over_x() {
        let e = Expr::Frac(Expr::int(1), Expr::sym("x"));
        let e_arc = Arc::new(e.clone());
        let r = diff_simplified(e_arc);
        assert_eq!(format_expr(r.as_ref()), "(0*x-1*1*1)/(x^2)");
        eval_diff_matches(&e, r.as_ref());
    }
}
