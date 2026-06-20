//! giac-solve plugin: wires equation solving into `giac-core::Context`.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{Context, EvalError, ExprArc, SolvePlugin};

use crate::froot::eval_froot;
use crate::realroot::eval_realroot;
use crate::solve::eval_solve;
use crate::sturm::{eval_sturm, eval_sturmab};
use crate::fsolve::eval_fsolve;

/// Default implementation of [`SolvePlugin`].
pub struct DefaultSolvePlugin;

impl SolvePlugin for DefaultSolvePlugin {
    // **Stable (bounded)** — solve via poly roots, rootof, or linsolve
    fn eval_solve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_solve(args, ctx)
    }

    // **Stable** — `Poly::eval_linsolve`
    fn eval_linsolve(
        &self,
        eqs: &ExprArc,
        vars: &ExprArc,
        ctx: &Context,
    ) -> Result<ExprArc, EvalError> {
        giac_linalg::eval_linsolve(eqs, vars, ctx)
    }

    // **Partial** — Newton numeric solve stub
    fn eval_fsolve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_fsolve(args, ctx)
    }

    // **Stable** — Sturm sequence for univariate poly
    fn eval_sturm(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_sturm(args, ctx)
    }

    // **Stable** — root count in (a,b) via Sturm
    fn eval_sturmab(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_sturmab(args, ctx)
    }

    // **Stable (bounded)** — real roots via Sturm isolation
    fn eval_realroot(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_realroot(args, ctx)
    }

    // **Stable (bounded)** — rational roots of univariate poly
    fn eval_froot(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_froot(args, ctx)
    }
}

/// Install the default solve plugin on `ctx`.
/// **Stable** — register DefaultSolvePlugin
pub fn install_solve(ctx: &mut Context) {
    ctx.set_solve_plugin(Arc::new(DefaultSolvePlugin));
}

/// `giac_linalg::xcas_default()` with equation solving enabled.
/// **Stable** — Context with simplify plugin
pub fn xcas_default() -> Context {
    let mut ctx = giac_linalg::xcas_default();
    install_solve(&mut ctx);
    ctx
}

#[cfg(test)]
mod tests {
    //! Test tiers — `.doc/test-writing-spec.md` · audit §2

    use std::sync::Arc;

    use giac_core::{eval, format_expr, Expr, FuncKind, RelOp};

    use super::xcas_default;

    // **A** — eval(Solve) via plugin; display snapshot of solution list.
    #[test]
    fn eval_solve_via_plugin() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Solve,
            vec![
                Arc::new(Expr::Relation(
                    RelOp::Eq,
                    Expr::add(vec![
                        Expr::pow(Expr::sym("x"), Expr::int(2)),
                        Expr::mul(vec![Expr::int(-2), Expr::sym("x")]),
                        Expr::int(1),
                    ]),
                    Expr::int(0),
                )),
                Expr::sym("x"),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "[1]");
    }
}
