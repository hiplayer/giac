//! Plugin trait so `giac-calculus` can implement calculus ops without a
//! circular `giac-core` ↔ `giac-calculus` dependency.

use std::sync::Arc;

use crate::{Context, EvalError, ExprArc};

/// Symbolic integration, differentiation, limits, and series.
pub trait CalculusPlugin: Send + Sync {
    fn eval_integrate(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_diff(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_limit(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_series(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_risch(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
}

impl Context {
    pub fn set_calculus_plugin(&mut self, plugin: Arc<dyn CalculusPlugin>) {
        self.calculus_plugin = Some(plugin);
    }

    pub(crate) fn calculus(&self) -> Result<&Arc<dyn CalculusPlugin>, EvalError> {
        self.calculus_plugin
            .as_ref()
            .ok_or(EvalError::NotImplemented("calculus plugin not installed"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EvalError, Expr, FuncKind};

    #[test]
    fn calculus_missing_plugin_returns_not_implemented() {
        let ctx = Context::xcas_default();
        let err = crate::eval(
            Expr::func(
                FuncKind::Integrate,
                vec![Expr::sym("x"), Expr::sym("x")],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap_err();
        assert!(matches!(err, EvalError::NotImplemented(_)));
    }
}
