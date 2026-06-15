//! giac-ode plugin: wires `desolve` into `giac-core::Context`.

use std::sync::Arc;

use giac_core::{Context, EvalError, ExprArc, OdePlugin};

use crate::stubs::eval_desolve;

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
    use giac_core::{eval, Expr, FuncKind};

    use super::xcas_default;

    #[test]
    fn desolve_stub_returns_not_implemented() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Desolve,
            vec![Expr::sym("y"), Expr::sym("x")],
        );
        assert!(eval(e.as_ref(), &ctx).is_err());
    }
}
