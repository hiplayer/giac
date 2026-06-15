//! Placeholders for ODE solving (Phase 4).

use giac_core::{Context, EvalError, ExprArc};

pub fn eval_desolve(_args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    Err(EvalError::NotImplemented("desolve"))
}
