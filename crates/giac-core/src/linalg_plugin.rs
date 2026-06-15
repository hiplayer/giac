//! Plugin trait so `giac-linalg` can implement symbolic/numeric matrix ops without a
//! circular `giac-core` ↔ `giac-linalg` dependency.

use std::sync::Arc;

use crate::{Context, EvalError, ExprArc, Ident};

/// Symbolic and numeric linear algebra operations on `Expr::Matrix`.
pub trait LinalgPlugin: Send + Sync {
    fn eval_matrix_mul(&self, a: &ExprArc, b: &ExprArc) -> Result<ExprArc, EvalError>;
    fn eval_matrix_pow(&self, m: &ExprArc, exp: u32, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_idn(&self, n: usize) -> ExprArc;
    fn eval_rref(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_inv(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_det(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_tran(&self, m: &ExprArc) -> Result<ExprArc, EvalError>;
    fn eval_ker(&self, m: &ExprArc) -> Result<ExprArc, EvalError>;
    fn eval_image(&self, m: &ExprArc) -> Result<ExprArc, EvalError>;
    fn eval_pcar(&self, m: &ExprArc) -> Result<ExprArc, EvalError>;
    fn eval_charpoly(
        &self,
        m: &ExprArc,
        var: &Ident,
        ctx: &Context,
    ) -> Result<ExprArc, EvalError>;
    fn eval_linsolve(
        &self,
        eqs: &ExprArc,
        vars: &ExprArc,
        ctx: &Context,
    ) -> Result<ExprArc, EvalError>;
    fn eval_jordan(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_egv(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_lu(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_qr(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_svd(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_gramschmidt(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError>;
    fn eval_trace(&self, m: &ExprArc) -> Result<ExprArc, EvalError>;
}

impl Context {
    pub fn set_linalg_plugin(&mut self, plugin: Arc<dyn LinalgPlugin>) {
        self.linalg_plugin = Some(plugin);
    }

    pub(crate) fn linalg(&self) -> Result<&Arc<dyn LinalgPlugin>, EvalError> {
        self.linalg_plugin
            .as_ref()
            .ok_or(EvalError::NotImplemented("linalg plugin not installed"))
    }
}
