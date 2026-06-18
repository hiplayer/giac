mod tower;
mod pow2expln;
mod hermite;
mod algebraic_rt;
mod rothstein_trager;

pub use tower::{risch_tower, rlvarx, RischTowerError};
pub use crate::expr_util::depends_on_var;
pub use pow2expln::pow2expln;
pub use hermite::{hermite_reduce, HermiteTerm};
pub use algebraic_rt::{
    integrate_monic_x4_plus_one, is_monic_even_quartic, is_monic_x4_plus_one,
    try_algebraic_rt_log_part,
};
pub use rothstein_trager::{
    rothstein_trager_integrate, try_algebraic_rt_even_quartic, try_integrate_x4_plus_one,
};

use giac_core::{Context, EvalError, ExprArc};

use crate::eval_integrate::eval_integrate;

/// `risch(f,x)` — minimal subset: delegate to `integrate` (GIAC-217).
pub fn eval_risch(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    eval_integrate(args, ctx)
}
