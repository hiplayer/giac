//! Placeholders for `realroot` (GIAC-207).

use giac_core::{Context, EvalError, ExprArc};

pub use crate::fsolve::eval_fsolve;
pub use crate::sturm::{eval_sturm, eval_sturmab};

pub fn eval_realroot(_args: &[ExprArc], _ctx: &Context) -> Result<ExprArc, EvalError> {
    Err(EvalError::NotImplemented("realroot"))
}
