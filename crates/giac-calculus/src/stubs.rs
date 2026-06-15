//! Placeholders for `limit` and `series` (Phase 5).

use giac_core::{Context, EvalError, ExprArc};

pub fn eval_risch(_args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    Err(EvalError::NotImplemented("risch"))
}

pub fn eval_limit(_args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    Err(EvalError::NotImplemented("limit"))
}

pub fn eval_series(_args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    Err(EvalError::NotImplemented("series"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::EvalError;

    #[test]
    fn stubs_return_not_implemented() {
        let ctx = giac_core::Context::xcas_default();
        assert!(matches!(
            eval_limit(&[], &ctx),
            Err(EvalError::NotImplemented(_))
        ));
        assert!(matches!(
            eval_series(&[], &ctx),
            Err(EvalError::NotImplemented(_))
        ));
    }
}
