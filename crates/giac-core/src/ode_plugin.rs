//! Plugin trait so `giac-ode` can implement ODE solving without a
//! circular `giac-core` ↔ `giac-ode` dependency.

use std::sync::Arc;

use crate::{Context, EvalError, ExprArc};

/// Ordinary differential equation solving (`desolve`, …).
pub trait OdePlugin: Send + Sync {
    fn eval_desolve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
}

impl Context {
    pub fn set_ode_plugin(&mut self, plugin: Arc<dyn OdePlugin>) {
        self.ode_plugin = Some(plugin);
    }

    pub(crate) fn ode(&self) -> Result<&Arc<dyn OdePlugin>, EvalError> {
        self.ode_plugin
            .as_ref()
            .ok_or(EvalError::NotImplemented("ode plugin not installed"))
    }
}
