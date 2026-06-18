//! giac-calculus plugin: wires integration/differentiation into `giac-core::Context`.
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../.doc/giac-calculus-api-stability.md) §2.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `install_calculus`, `xcas_default`, `DefaultCalculusPlugin` trait hooks |

use std::sync::Arc;

use giac_core::{CalculusPlugin, Context, EvalError, ExprArc};

use crate::eval_diff::eval_diff;
use crate::eval_integrate::eval_integrate;
use crate::limit::eval_limit;
use crate::risch::eval_risch;
use crate::series::eval_series;

/// Default implementation of [`CalculusPlugin`].
pub struct DefaultCalculusPlugin;

impl CalculusPlugin for DefaultCalculusPlugin {
    /// **Stable** — delegate to [`eval_integrate`].
    fn eval_integrate(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_integrate(args, ctx)
    }

    /// **Stable** — delegate to [`eval_diff`].
    fn eval_diff(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_diff(args, ctx)
    }

    /// **Stable** — delegate to [`eval_limit`].
    fn eval_limit(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_limit(args, ctx)
    }

    /// **Stable** — delegate to [`eval_series`].
    fn eval_series(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_series(args, ctx)
    }

    /// **Partial** — delegate to [`eval_risch`] (narrow Risch subset).
    fn eval_risch(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_risch(args, ctx)
    }
}

/// **Stable** — install the default calculus plugin on `ctx`.
pub fn install_calculus(ctx: &mut Context) {
    ctx.set_calculus_plugin(Arc::new(DefaultCalculusPlugin));
}

/// **Stable** — full CAS context: linear algebra, solving, simplification, and calculus.
pub fn xcas_default() -> Context {
    let mut ctx = giac_solve::xcas_default();
    giac_simplify::install_simplify(&mut ctx);
    install_calculus(&mut ctx);
    ctx
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::xcas_default;

    #[test]
    fn eval_integrate_via_plugin() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![Expr::sym("x"), Expr::sym("x")],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("1/2") && s.contains("x^2"), "got {s}");
    }

    #[test]
    fn eval_limit_via_plugin() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Limit,
            vec![
                Expr::pow(Expr::sym("x"), Expr::int(-1)),
                Expr::sym("x"),
                Expr::sym("+infinity"),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "0");
    }

    #[test]
    fn eval_series_via_plugin() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Series,
            vec![
                Expr::sym("x"),
                Expr::sym("x"),
                Expr::int(0),
                Expr::int(3),
            ],
        );
        assert!(eval(e.as_ref(), &ctx).is_ok());
    }
}
