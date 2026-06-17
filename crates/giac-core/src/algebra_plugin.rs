//! Plugin trait so `giac-simplify` can implement algebra ops without a
//! circular `giac-core` ↔ `giac-simplify` dependency.

use std::sync::Arc;

use num_bigint::BigInt;

use crate::{Context, EvalError, Expr, ExprArc};

/// Polynomial normalization, expansion, and factorization.
pub trait AlgebraPlugin: Send + Sync {
    fn normal(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn ratnormal(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn expand(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn factor(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn ifactor(&self, n: &BigInt) -> ExprArc;
    fn texpand(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn halftan(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn lin(&self, expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError>;
}

impl Context {
    pub fn set_algebra_plugin(&mut self, plugin: Arc<dyn AlgebraPlugin>) {
        self.algebra_plugin = Some(plugin);
    }

    pub(crate) fn algebra(&self) -> Result<&Arc<dyn AlgebraPlugin>, EvalError> {
        self.algebra_plugin
            .as_ref()
            .ok_or(EvalError::NotImplemented("algebra plugin not installed"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EvalError, Expr, FuncKind};

    #[test]
    fn algebra_missing_plugin_returns_not_implemented() {
        let ctx = Context::xcas_default();
        let err = crate::eval(
            Expr::func(FuncKind::Normal, vec![Expr::sym("x")]).as_ref(),
            &ctx,
        )
        .unwrap_err();
        assert!(matches!(err, EvalError::NotImplemented(_)));
    }
}
