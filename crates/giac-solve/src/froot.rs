//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{eval, expr_to_rational_polys_eval, Context, EvalError, Expr, ExprArc, Ident};
use giac_poly::{square_free_factorization, Poly, Var};

/// `froot(p)` or `froot(p,x)` — factor roots with signed multiplicities (upstream `misc.cc` `_froot`).
/// **Stable (bounded)** — rational roots of univariate poly
pub fn eval_froot(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    let (expr, var) = parse_froot_args(args, ctx)?;
    let ev = eval(expr.as_ref(), ctx)?;
    let (num, den) = expr_to_rational_polys_eval(ev.as_ref(), ctx)?;
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

// **Pipeline private** — `parse_froot_args`
fn parse_froot_args(args: &[ExprArc], _ctx: &Context) -> Result<(ExprArc, Ident), EvalError> {
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

// **Pipeline private** — `append_factor_roots`
fn append_factor_roots(
    poly: &Poly,
    var: &Var,
    sign: i32,
    out: &mut Vec<(ExprArc, i32)>,
    ctx: &Context,
) -> Result<(), EvalError> {
    if poly.is_zero() {
        return Ok(());
    }
    let factors = square_free_factorization(poly, var)?;
    for (factor, mult) in factors {
        let signed = mult as i32 * sign;
        for root in solve_factor_roots(&factor, var, ctx)? {
            out.push((root, signed));
        }
    }
    Ok(())
}

// **Pipeline private** — S3: deg≤4 via S0 kernel; deg≥5 rootof branch
fn solve_factor_roots(
    factor: &Poly,
    var: &Var,
    ctx: &Context,
) -> Result<Vec<ExprArc>, EvalError> {
    crate::solve_poly::solve_irreducible_factor(factor, var, ctx)
}

#[cfg(test)]
mod tests {
    //! Test tiers — `.doc/test-writing-spec.md` · audit §2

    use giac_core::{eval, Expr, FuncKind};

    use crate::plugin::xcas_default;
    use crate::test_verify::{assert_froot_has_root, list_items};

    // **A** — eval(Froots); assert_froot_has_root.
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
        let items = list_items(&r);
        for root in [Expr::int(0), Expr::int(1), Expr::int(2)] {
            assert_froot_has_root(items, &root, &ctx);
        }
    }

    // **B** — S3: deg-4 factor via poly_algext_roots (same kernel as solve).
    #[test]
    fn froot_quartic_t4_plus_t_plus_1() {
        let ctx = xcas_default();
        let p = Expr::add(vec![
            Expr::pow(Expr::sym("t"), Expr::int(4)),
            Expr::sym("t"),
            Expr::int(1),
        ]);
        let e = Expr::func(FuncKind::Froot, vec![p, Expr::sym("t")]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let items = list_items(&r);
        assert_eq!(items.len(), 8, "4 roots × (root, mult) pairs");
    }

    // **A** — eval(Froot); rational roots in flat list.
    #[test]
    fn froot_linear_factor() {
        let ctx = xcas_default();
        let p = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(-3)]),
            Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2)),
        ]);
        let e = Expr::func(FuncKind::Froot, vec![p]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let items = list_items(&r);
        assert_froot_has_root(items, &Expr::int(3), &ctx);
        assert_froot_has_root(items, &Expr::int(-1), &ctx);
    }
}
