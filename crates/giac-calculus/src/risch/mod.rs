mod tower;
mod pow2expln;
mod hermite;
mod rothstein_trager;

pub use tower::{depends_on_var, risch_tower, rlvarx, RischTowerError};
pub use pow2expln::pow2expln;
pub use hermite::{hermite_reduce, HermiteTerm};
pub use rothstein_trager::{
    integrate_one_over_x4_plus_one_squared, rothstein_trager_integrate, try_integrate_x4_plus_one,
};

use giac_core::{Context, EvalError, ExprArc};

use crate::eval_integrate::eval_integrate;

/// `risch(f,x)` — minimal subset: delegate to `integrate` (GIAC-217).
pub fn eval_risch(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    eval_integrate(args, ctx)
}
