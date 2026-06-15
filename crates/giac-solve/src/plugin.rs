//! giac-solve plugin: wires equation solving into `giac-core::Context`.

use std::sync::Arc;

use giac_core::{Context, EvalError, ExprArc, SolvePlugin};

use crate::solve::eval_solve;
use crate::stubs::{eval_fsolve, eval_realroot, eval_sturm};

/// Default implementation of [`SolvePlugin`].
pub struct DefaultSolvePlugin;

impl SolvePlugin for DefaultSolvePlugin {
    fn eval_solve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_solve(args, ctx)
    }

    fn eval_linsolve(
        &self,
        eqs: &ExprArc,
        vars: &ExprArc,
        ctx: &Context,
    ) -> Result<ExprArc, EvalError> {
        giac_linalg::eval_linsolve(eqs, vars, ctx)
    }

    fn eval_fsolve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_fsolve(args, ctx)
    }

    fn eval_sturm(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_sturm(args, ctx)
    }

    fn eval_realroot(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_realroot(args, ctx)
    }
}

/// Install the default solve plugin on `ctx`.
pub fn install_solve(ctx: &mut Context) {
    ctx.set_solve_plugin(Arc::new(DefaultSolvePlugin));
}

/// `giac_linalg::xcas_default()` with equation solving enabled.
pub fn xcas_default() -> Context {
    let mut ctx = giac_linalg::xcas_default();
    install_solve(&mut ctx);
    ctx
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{eval, format_expr, Expr, FuncKind, RelOp};

    use super::xcas_default;

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
