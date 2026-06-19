use std::sync::Arc;

use giac_core::{eval, expr_to_poly, poly_to_expr, Context, EvalError, Expr, ExprArc, FuncKind, Ident};
use giac_poly::{
    coeff_at, roots, square_free_factorization, univariate_degree, Poly, PolyError, Var,
};
use num_rational::Ratio;
use num_traits::{Signed, Zero};

/// `froot(p)` or `froot(p,x)` — factor roots with signed multiplicities (upstream `misc.cc` `_froot`).
pub fn eval_froot(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    let (expr, var) = parse_froot_args(args, ctx)?;
    let ev = eval(expr.as_ref(), ctx)?;
    let (num, den) = rational_num_den(ev.as_ref(), ctx)?;
    let v = Var::from(var.as_str());
    let mut out = Vec::new();
    append_factor_roots(&num, &v, 1, &mut out, ctx)?;
    append_factor_roots(&den, &v, -1, &mut out, ctx)?;
    let items: Vec<ExprArc> = out
        .into_iter()
        .flat_map(|(root, mult)| vec![root, Expr::int(i64::from(mult))])
        .collect();
    Ok(Arc::new(Expr::List(items)))
}

fn parse_froot_args(args: &[ExprArc], ctx: &Context) -> Result<(ExprArc, Ident), EvalError> {
    match args.len() {
        1 => Ok((Arc::clone(&args[0]), Ident::new("x"))),
        2 => {
            let var = match args[1].as_ref() {
                Expr::Symbol(id) => id.clone(),
                _ => return Err(EvalError::TypeError("variable name expected")),
            };
            Ok((Arc::clone(&args[0]), var))
        }
        _ => Err(EvalError::TooFewArgs("froot")),
    }
}

fn rational_num_den(expr: &Expr, ctx: &Context) -> Result<(Poly, Poly), EvalError> {
    match expr {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) => Ok((
            Poly::one(),
            expr_to_poly(eval(base.as_ref(), ctx)?.as_ref())?,
        )),
        Expr::Frac(num, den) => Ok((
            expr_to_poly(eval(num.as_ref(), ctx)?.as_ref())?,
            expr_to_poly(eval(den.as_ref(), ctx)?.as_ref())?,
        )),
        Expr::Mul(factors) => {
            if factors.len() == 2 {
                if let Expr::Pow(base, exp) = factors[1].as_ref() {
                    if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) {
                        return Ok((
                            expr_to_poly(eval(factors[0].as_ref(), ctx)?.as_ref())?,
                            expr_to_poly(eval(base.as_ref(), ctx)?.as_ref())?,
                        ));
                    }
                }
            }
            Ok((expr_to_poly(eval(expr, ctx)?.as_ref())?, Poly::one()))
        }
        other => Ok((expr_to_poly(eval(other, ctx)?.as_ref())?, Poly::one())),
    }
}

fn append_factor_roots(
    poly: &Poly,
    var: &Var,
    sign: i32,
    out: &mut Vec<(ExprArc, i32)>,
    _ctx: &Context,
) -> Result<(), EvalError> {
    if poly.is_zero() {
        return Ok(());
    }
    let factors = square_free_factorization(poly, var).map_err(poly_err)?;
    for (factor, mult) in factors {
        let signed = mult as i32 * sign;
        for root in solve_factor_roots(&factor, var)? {
            out.push((root, signed));
        }
    }
    Ok(())
}

fn solve_factor_roots(factor: &Poly, var: &Var) -> Result<Vec<ExprArc>, EvalError> {
    let d = univariate_degree(factor, var);
    match d {
        0 => Ok(vec![]),
        1 => {
            let a = coeff_at(factor, var, 1);
            if a.is_zero() {
                return Err(EvalError::TypeError("degenerate linear factor"));
            }
            let root = -coeff_at(factor, var, 0) / a;
            Ok(vec![poly_to_expr(&Poly::constant(root))])
        }
        2 | 3 => roots(factor, var)
            .map_err(poly_err)?
            .into_iter()
            .map(|p| Ok(poly_to_expr(&p)))
            .collect(),
        _ => Err(EvalError::NotImplemented("froot")),
    }
}

fn poly_err(e: PolyError) -> EvalError {
    match e {
        PolyError::NotImplemented(s) => EvalError::NotImplemented(s),
        PolyError::TypeError(s) => EvalError::TypeError(s),
        PolyError::DivisionByZero => EvalError::TypeError("division by zero"),
    }
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn froots_flanex_line195() {
        let ctx = xcas_default();
        let num = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(5)),
            Expr::mul(vec![Expr::int(-2), Expr::pow(Expr::sym("x"), Expr::int(4))]),
            Expr::pow(Expr::sym("x"), Expr::int(3)),
        ]);
        let den = Expr::add(vec![Expr::sym("x"), Expr::int(-2)]);
        let rat = Expr::mul(vec![num, Expr::pow(den, Expr::int(-1))]);
        let e = Expr::func(FuncKind::Froots, vec![rat]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("0"), "got {s}");
        assert!(s.contains("1"), "got {s}");
        assert!(s.contains("2"), "got {s}");
    }

    #[test]
    fn froot_linear_factor() {
        let ctx = xcas_default();
        let p = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(-3)]),
            Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2)),
        ]);
        let e = Expr::func(FuncKind::Froot, vec![p]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("3"), "got {s}");
        assert!(s.contains("-1"), "got {s}");
    }
}
