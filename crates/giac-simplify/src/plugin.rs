//! giac-simplify plugin: wires polynomial normalization into `giac-core::Context`.
//!
//! **API inventory:** inline tier on every function; full table in
//! `.doc/giac-simplify-api-stability.md` (Per-file inventory).
//!
use std::sync::Arc;

use num_bigint::BigInt;

use giac_core::{AlgebraPlugin, Context, EvalError, Expr, ExprArc};

use crate::{
    expand, factor, halftan, ifactor, lin, normal, ratnormal, texpand,
};

/// **Stable** — default [`AlgebraPlugin`] wiring all simplify entry points.
pub struct DefaultAlgebraPlugin;

impl AlgebraPlugin for DefaultAlgebraPlugin {
    // **Stable** — expand then polynomial collect
    fn normal(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
        normal(expr, ctx)
    }

    // **Stable** — single fraction in lowest terms
    fn ratnormal(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
        ratnormal(expr, ctx)
    }

    // **Stable** — expand with Full policy
    fn expand(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
        expand(expr, ctx)
    }

    // **Stable (bounded)** — structural then giac-poly factor
    fn factor(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
        factor(expr, ctx)
    }

    // **Stable** — integer prime factorization display
    fn ifactor(&self, n: &BigInt) -> ExprArc {
        ifactor(n)
    }

    // **Partial** — trig/exp/ln arg expand then algebraic expand
    fn texpand(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
        texpand(expr, ctx)
    }

    // **Partial** — half-angle tan on narrow pattern
    fn halftan(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
        halftan(expr, ctx)
    }

    // **Partial** — linearize exp products and (exp+1)^2
    fn lin(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
        lin(expr, ctx)
    }
}

/// **Stable** — install [`DefaultAlgebraPlugin`] on `ctx`.
pub fn install_simplify(ctx: &mut Context) {
    ctx.set_algebra_plugin(Arc::new(DefaultAlgebraPlugin));
}

/// **Stable** — `giac-core::Context::xcas_default()` with simplification enabled.
pub fn xcas_default() -> Context {
    let mut ctx = giac_core::Context::xcas_default();
    install_simplify(&mut ctx);
    ctx
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::xcas_default;

    #[test]
    fn eval_normal_via_plugin() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Normal,
            vec![Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2))],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x^2+2*x+1");
    }

    #[test]
    fn eval_expand_binomial_via_plugin() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Normal,
            vec![Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(3)]), Expr::int(4))],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(
            format_expr(r.as_ref()),
            "x^4+12*x^3+54*x^2+108*x+81"
        );
    }

    #[test]
    fn eval_factor_via_plugin() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Factor,
            vec![Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(4)),
                Expr::int(-1),
            ])],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("(x-1)") && s.contains("(x+1)"));
        assert_eq!(s, "(x-1)*(x+1)*(x^2+1)");
    }

    #[test]
    fn eval_factor_x_squared_minus_two_rootof() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Factor,
            vec![Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::int(-2),
            ])],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("rootof"), "got {s}");
    }
}
