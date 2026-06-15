//! Plugin trait so `giac-solve` can implement equation solving without a
//! circular `giac-core` ↔ `giac-solve` dependency.

use std::sync::Arc;

use crate::{Context, EvalError, ExprArc};

/// Polynomial and system solving (`solve`, `linsolve`, …).
pub trait SolvePlugin: Send + Sync {
    fn eval_solve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_linsolve(
        &self,
        eqs: &ExprArc,
        vars: &ExprArc,
        ctx: &Context,
    ) -> Result<ExprArc, EvalError>;
    fn eval_fsolve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_sturm(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_realroot(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
}

impl Context {
    pub fn set_solve_plugin(&mut self, plugin: Arc<dyn SolvePlugin>) {
        self.solve_plugin = Some(plugin);
    }

    pub(crate) fn solve(&self) -> Result<&Arc<dyn SolvePlugin>, EvalError> {
        self.solve_plugin
            .as_ref()
            .ok_or(EvalError::NotImplemented("solve plugin not installed"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EvalError, Expr, FuncKind};

    #[test]
    fn solve_missing_plugin_returns_not_implemented() {
        let ctx = Context::xcas_default();
        let err = crate::eval(
            Expr::func(
                FuncKind::Solve,
                vec![Expr::sym("x"), Expr::sym("x")],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap_err();
        assert!(matches!(err, EvalError::NotImplemented(_)));
    }
}
