//! Risch integration subset (tower, Hermite, Rothstein–Trager).
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../.doc/giac-calculus-api-stability.md) §5.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `hermite_reduce`, `pow2expln`, `rlvarx`, `risch_tower` |
//! | **Partial** | `eval_risch`, `rothstein_trager_integrate`, `try_algebraic_rt_*`, `try_integrate_x4_plus_one` |
//! | **Pipeline private** | `algebraic_rt.rs` / `tower.rs` 内部分解 |

mod tower;
mod pow2expln;
mod hermite;
mod algebraic_rt;
mod rothstein_trager;

pub use tower::{risch_tower, rlvarx, RischTowerError};

pub use pow2expln::pow2expln;
pub use hermite::{hermite_reduce, HermiteTerm};

pub use rothstein_trager::{
    rothstein_trager_integrate, try_algebraic_rt_even_quartic,
};

use giac_core::{Context, EvalError, ExprArc};

use crate::eval_integrate::eval_integrate;

/// **Partial** — `risch(f,x)` minimal subset; delegates to `integrate` (GIAC-217). **退役：** full Risch decision procedure.
pub fn eval_risch(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    eval_integrate(args, ctx)
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, Expr, FuncKind};

    use crate::plugin::xcas_default;

    

    #[test]
    fn eval_risch_delegates_to_integrate() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Risch,
            vec![Expr::sym("x"), Expr::sym("x")],
        );
        assert!(eval(e.as_ref(), &ctx).is_ok());
    }
}
