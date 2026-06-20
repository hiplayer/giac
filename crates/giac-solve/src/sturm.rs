//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{
    eval, expr_to_poly, poly_to_expr, Context, EvalError, Expr, ExprArc, Ident,
};
use giac_poly::{
    sturm_sequence, sturmab_count, univariate_degree, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;

/// **Stable** — Sturm sequence for univariate poly
pub fn eval_sturm(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    let (poly, var) = poly_and_var(args, ctx)?;
    let mut q = giac_poly::odd_multiplicity_part(&poly, &var)?;
    if univariate_degree(&q, &var) == 0 {
        if univariate_degree(&poly, &var) == 0 {
            return Err(EvalError::TypeError("constant polynomial"));
        }
        q = giac_poly::square_free_part(&poly, &var)?;
        if univariate_degree(&q, &var) == 0 {
            return Err(EvalError::TypeError("constant polynomial"));
        }
    }
    let seq = sturm_sequence(&q, &var)?;
    let items: Vec<ExprArc> = seq.into_iter().map(|p| poly_to_expr(&p)).collect();
    Ok(Arc::new(Expr::List(items)))
}

/// **Stable** — root count in (a,b) via Sturm
pub fn eval_sturmab(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 4 {
        return Err(EvalError::TooFewArgs("sturmab"));
    }
    let poly = expr_to_poly(eval(args[0].as_ref(), ctx)?.as_ref())?;
    let var = ident_from_expr(args[1].as_ref())?;
    let v = Var::from(var.as_str());
    if univariate_degree(&poly, &v) == 0 {
        return Err(EvalError::TypeError("constant polynomial"));
    }
    let a = eval_to_rational(&args[2], ctx)?;
    let b = eval_to_rational(&args[3], ctx)?;
    let count = sturmab_count(&poly, &v, &a, &b)?;
    Ok(Expr::int(i64::try_from(count).unwrap_or(0)))
}

// **Pipeline private** — `poly_and_var`
fn poly_and_var(args: &[ExprArc], ctx: &Context) -> Result<(giac_poly::Poly, Var), EvalError> {
    match args.len() {
        1 => {
            let poly = expr_to_poly(eval(args[0].as_ref(), ctx)?.as_ref())?;
            Ok((poly, Var::from("x")))
        }
        2 => {
            let poly = expr_to_poly(eval(args[0].as_ref(), ctx)?.as_ref())?;
            let var = ident_from_expr(args[1].as_ref())?;
            Ok((poly, Var::from(var.as_str())))
        }
        _ => Err(EvalError::TooFewArgs("sturm")),
    }
}

// **Pipeline private** — `ident_from_expr`
fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("variable name expected")),
    }
}

// **Pipeline private** — `eval_to_rational`
fn eval_to_rational(e: &ExprArc, ctx: &Context) -> Result<Ratio<BigInt>, EvalError> {
    let ev = eval(e.as_ref(), ctx)?;
    match ev.as_ref() {
        Expr::Int(n) => Ok(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Ok(r.clone()),
        _ => Err(EvalError::TypeError("rational bound expected")),
    }
}


#[cfg(test)]
mod tests {
    //! Test tiers — `.doc/test-writing-spec.md` · audit §2

    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    // **A** — eval(Sturm); list length only (TODO: semantic Sturm sequence).
    #[test]
    fn sturm_x_cubed_plus_one_squared() {
        let ctx = xcas_default();
        let p = Expr::pow(
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(3)), Expr::int(1)]),
            Expr::int(2),
        );
        let e = Expr::func(FuncKind::Sturm, vec![p]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert!(matches!(r.as_ref(), Expr::List(items) if items.len() >= 2));
    }

    // **A** — eval(Sturm); sequence length (TODO: coefficient checks).
    #[test]
    fn sturm_x_cubed_plus_one() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Sturm,
            vec![
                Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(3)), Expr::int(1)]),
                Expr::sym("x"),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert!(matches!(r.as_ref(), Expr::List(items) if items.len() == 3));
    }

    // **A** + **C** — eval(Sturmab); display golden root count in interval.
    #[test]
    fn sturmab_counts_root_in_interval() {
        let ctx = xcas_default();
        let p = Expr::mul(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(3)), Expr::int(2)]),
        ]);
        let e = Expr::func(
            FuncKind::Sturmab,
            vec![p, Expr::sym("x"), Expr::int(-2), Expr::int(0)],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1");
    }
}
