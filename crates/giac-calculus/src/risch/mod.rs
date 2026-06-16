mod tower;
mod pow2expln;
mod hermite;

pub use tower::{depends_on_var, risch_tower, rlvarx, RischTowerError};
pub use pow2expln::pow2expln;
pub use hermite::{hermite_reduce, HermiteTerm};

use giac_core::{Context, EvalError, ExprArc};

use crate::eval_integrate::eval_integrate;

/// `risch(f,x)` — minimal subset: delegate to `integrate` (GIAC-217).
pub fn eval_risch(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    eval_integrate(args, ctx)
}
