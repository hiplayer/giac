//! giac-ode plugin: wires `desolve` into `giac-core::Context`.

use std::sync::Arc;

use giac_core::{Context, EvalError, ExprArc, OdePlugin};

use crate::desolve::eval_desolve;

/// Default implementation of [`OdePlugin`].
pub struct DefaultOdePlugin;

impl OdePlugin for DefaultOdePlugin {
    fn eval_desolve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_desolve(args, ctx)
    }
}

/// Install the default ODE plugin on `ctx`.
pub fn install_ode(ctx: &mut Context) {
    ctx.set_ode_plugin(Arc::new(DefaultOdePlugin));
}

/// Full CAS context: linear algebra, solve, calculus, and ODE.
pub fn xcas_default() -> Context {
    let mut ctx = giac_calculus::xcas_default();
    install_ode(&mut ctx);
    ctx
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{eval, format_expr, Expr, FuncKind, RelOp};

    use super::xcas_default;

    #[test]
    fn desolve_via_plugin() {
        let ctx = xcas_default();
        let eq = Arc::new(Expr::Relation(
            RelOp::Eq,
            Expr::func(FuncKind::Prime, vec![Expr::sym("y"), Expr::int(1)]),
            Expr::mul(vec![Expr::sym("x"), Expr::sym("y")]),
        ));
        let e = Expr::func(
            FuncKind::Desolve,
            vec![
                eq,
                Expr::func(FuncKind::Apply, vec![Expr::sym("y"), Expr::sym("x")]),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert!(format_expr(r.as_ref()).contains("c0"));
    }
}
