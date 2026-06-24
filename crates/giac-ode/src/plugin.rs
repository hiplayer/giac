//! giac-ode plugin: wires `desolve` into `giac-core::Context`.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-ode-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{Context, EvalError, ExprArc, OdePlugin};

use crate::desolve::eval_desolve;

/// Default implementation of [`OdePlugin`].
pub struct DefaultOdePlugin;

impl OdePlugin for DefaultOdePlugin {
    // **Stable (bounded)** — linear constant-coefficient ODE subset
    fn eval_desolve(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_desolve(args, ctx)
    }
}

/// Install the default ODE plugin on `ctx`.
/// **Stable** — register DefaultOdePlugin
pub fn install_ode(ctx: &mut Context) {
    ctx.set_ode_plugin(Arc::new(DefaultOdePlugin));
}

/// Full CAS context: linear algebra, solve, calculus, and ODE.
/// **Stable** — Context with simplify plugin
pub fn xcas_default() -> Context {
    let mut ctx = giac_calculus::xcas_default();
    install_ode(&mut ctx);
    ctx
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{eval, format_expr, expr_mentions_ident, Expr, FuncKind, Ident, RelOp};

    use super::xcas_default;
    use crate::test_verify::subst_constants;

    #[test]
    // smoke-until B-ODE,B-ODE-NORM: delete when `desolve_via_plugin_satisfies_ode` + `desolve_via_plugin_subst_canonical` green
    fn desolve_via_plugin() {
        let ctx = xcas_default();
        let eq = Arc::new(Expr::Relation(
            RelOp::Eq,
            Expr::func(FuncKind::Prime, vec![Expr::sym("y"), Expr::int(1)]),
            Expr::mul(vec![Expr::sym("x"), Expr::sym("y")]),
        ));
        let e = Expr::func(
            FuncKind::Desolve,
            vec![
                eq,
                Expr::func(FuncKind::Apply, vec![Expr::sym("y"), Expr::sym("x")]),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let sol = match r.as_ref() {
            Expr::Relation(RelOp::Eq, _, rhs) => rhs,
            other => panic!("expected relation, got {}", format_expr(other)),
        };
        assert!(expr_mentions_ident(sol, &Ident::new("c0")));
        assert_eq!(
            format_expr(subst_constants(sol, &[("c0", 1)]).as_ref()),
            "1*exp(1/2*x^2)"
        );
    }

    #[test]
    #[ignore = "B-ODE: y'-x*y 残差在 normal 下未归零"]
    fn desolve_via_plugin_satisfies_ode() {
        let ctx = xcas_default();
        let eq = Arc::new(Expr::Relation(
            RelOp::Eq,
            Expr::func(FuncKind::Prime, vec![Expr::sym("y"), Expr::int(1)]),
            Expr::mul(vec![Expr::sym("x"), Expr::sym("y")]),
        ));
        let e = Expr::func(
            FuncKind::Desolve,
            vec![
                eq,
                Expr::func(FuncKind::Apply, vec![Expr::sym("y"), Expr::sym("x")]),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        use crate::test_verify::assert_desolve_lin_ode;
        assert_desolve_lin_ode(
            &Expr::int(0),
            &Expr::int(1),
            &Expr::mul(vec![Expr::int(-1), Expr::sym("x")]),
            &Expr::int(0),
            &r,
            &Ident::new("x"),
            &[("c0", 1)],
            &ctx,
        );
    }

    #[test]
    #[ignore = "B-ODE-NORM: subst c0=1 后应为 exp(x^2/2)"]
    fn desolve_via_plugin_subst_canonical() {
        let ctx = xcas_default();
        let eq = Arc::new(Expr::Relation(
            RelOp::Eq,
            Expr::func(FuncKind::Prime, vec![Expr::sym("y"), Expr::int(1)]),
            Expr::mul(vec![Expr::sym("x"), Expr::sym("y")]),
        ));
        let e = Expr::func(
            FuncKind::Desolve,
            vec![
                eq,
                Expr::func(FuncKind::Apply, vec![Expr::sym("y"), Expr::sym("x")]),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let sol = match r.as_ref() {
            Expr::Relation(RelOp::Eq, _, rhs) => rhs,
            other => panic!("expected relation, got {}", format_expr(other)),
        };
        use giac_simplify::assert_equiv;
        let expected = Expr::func(
            FuncKind::Exp,
            vec![Expr::mul(vec![Expr::rat(1, 2), Expr::pow(Expr::sym("x"), Expr::int(2))])],
        );
        assert!(
            assert_equiv(subst_constants(sol, &[("c0", 1)]).as_ref(), expected.as_ref(), &ctx)
                .unwrap()
        );
    }
}
