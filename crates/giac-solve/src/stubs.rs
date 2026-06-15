//! Placeholders for numeric / Sturm root finders (Phase 5).

use giac_core::{Context, EvalError, ExprArc};

pub fn eval_fsolve(_args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    Err(EvalError::NotImplemented("fsolve"))
}

pub fn eval_sturm(_args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    Err(EvalError::NotImplemented("sturm"))
}

pub fn eval_realroot(_args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    Err(EvalError::NotImplemented("realroot"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::Expr;

    #[test]
    fn stubs_return_not_implemented() {
        let ctx = giac_core::Context::xcas_default();
        assert!(matches!(
            eval_fsolve(&[], &ctx),
            Err(EvalError::NotImplemented(_))
        ));
        assert!(matches!(
            eval_sturm(&[], &ctx),
            Err(EvalError::NotImplemented(_))
        ));
        assert!(matches!(
            eval_realroot(&[Expr::sym("x")], &ctx),
            Err(EvalError::NotImplemented(_))
        ));
    }
}
