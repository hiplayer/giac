//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    eval, eval_subst_map, Context, EvalError, Expr, ExprArc, Ident, RelOp,
};
use giac_linalg::f64_to_expr_numeric;

/// **Partial** — Newton numeric solve stub
pub fn eval_fsolve(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("fsolve"));
    }
    let var = ident_from_expr(args[1].as_ref())?;
    let f = equation_to_expr(&args[0])?;
    let root = newton(&f, &var, 1.0, ctx)?;
    Ok(f64_to_expr_numeric(root))
}

// **Pipeline private** — `equation_to_expr`
fn equation_to_expr(e: &ExprArc) -> Result<ExprArc, EvalError> {
    match e.as_ref() {
        Expr::Relation(RelOp::Eq, lhs, rhs) => Ok(Expr::add(vec![
            Arc::clone(lhs),
            Expr::mul(vec![Expr::int(-1), Arc::clone(rhs)]),
        ])),
        _ => Ok(Arc::clone(e)),
    }
}

// **Pipeline private** — `ident_from_expr`
fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("variable name expected")),
    }
}

// **Pipeline private** — `eval_at`
fn eval_at(f: &ExprArc, var: &Ident, x: f64, ctx: &Context) -> Result<f64, EvalError> {
    let mut subs = HashMap::new();
    subs.insert(var.clone(), f64_to_expr_numeric(x));
    let sub = eval_subst_map(f, &subs)?;
    let ev = eval(sub.as_ref(), ctx)?;
    expr_to_f64(&ev)
}

// **Pipeline private** — `expr_to_f64`
fn expr_to_f64(e: &ExprArc) -> Result<f64, EvalError> {
    match e.as_ref() {
        Expr::Int(n) => n
            .to_string()
            .parse()
            .map_err(|_| EvalError::TypeError("float conversion")),
        Expr::Rat(r) => {
            let num: f64 = r
                .numer()
                .to_string()
                .parse()
                .map_err(|_| EvalError::TypeError("float conversion"))?;
            let den: f64 = r
                .denom()
                .to_string()
                .parse()
                .map_err(|_| EvalError::TypeError("float conversion"))?;
            Ok(num / den)
        }
        _ => Err(EvalError::TypeError("numeric expected")),
    }
}

// **Pipeline private** — `newton`
fn newton(f: &ExprArc, var: &Ident, mut x: f64, ctx: &Context) -> Result<f64, EvalError> {
    const TOL: f64 = 1e-12;
    const H: f64 = 1e-8;
    for _ in 0..64 {
        let fx = eval_at(f, var, x, ctx)?;
        if fx.abs() < TOL {
            return Ok(x);
        }
        let fp = (eval_at(f, var, x + H, ctx)? - eval_at(f, var, x - H, ctx)?) / (2.0 * H);
        if fp.abs() < 1e-15 {
            break;
        }
        x -= fx / fp;
    }
    let fx = eval_at(f, var, x, ctx)?;
    if fx.abs() < 1e-8 {
        Ok(x)
    } else {
        Err(EvalError::NotImplemented("fsolve"))
    }
}

#[cfg(test)]
mod tests {
    //! Test tiers — `.doc/test-writing-spec.md` · audit §2

    use giac_core::{eval, Expr, FuncKind, RelOp};

    use super::*;
    use crate::plugin::xcas_default;

    // **A** — eval(Fsolve); numeric root satisfies x²=2.
    #[test]
    fn fsolve_x_squared_minus_two() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Fsolve,
            vec![
                Arc::new(Expr::Relation(
                    RelOp::Eq,
                    Expr::add(vec![
                        Expr::pow(Expr::sym("x"), Expr::int(2)),
                        Expr::int(-2),
                    ]),
                    Expr::int(0),
                )),
                Expr::sym("x"),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let v = expr_to_f64(&r).unwrap();
        assert!((v * v - 2.0).abs() < 1e-6, "got {v}");
    }
}
