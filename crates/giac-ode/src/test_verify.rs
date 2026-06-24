//! Semantic verification for giac-ode tests (`.doc/conformance-testing.md` §3).
//! See `.doc/issues/GIAC-expr-api-test-contains-cleanup.md` H3.

use std::collections::HashMap;
use std::sync::Arc;

use giac_calculus::diff;
use giac_core::{eval_subst_map, format_expr, Context, Expr, ExprArc, Ident, RelOp};
use giac_simplify::{assert_equiv, normal};

/// Verify `a2*y'' + a1*y' + a0*y - forcing ≡ 0` for concrete `sol(x)`.
pub fn assert_lin_ode_solution(
    a2: &Expr,
    a1: &Expr,
    a0: &Expr,
    forcing: &Expr,
    sol: &ExprArc,
    indep: &Ident,
    ctx: &Context,
) {
    let d1 = diff(sol, indep).expect("diff");
    let d2 = diff(&d1, indep).expect("diff2");
    let residual = normal(
        Expr::add(vec![
            Expr::mul(vec![Arc::new(a2.clone()), d2]),
            Expr::mul(vec![Arc::new(a1.clone()), d1]),
            Expr::mul(vec![Arc::new(a0.clone()), Arc::clone(sol)]),
            Expr::mul(vec![Expr::int(-1), Arc::new(forcing.clone())]),
        ])
        .as_ref(),
        ctx,
    )
    .expect("normal");
    assert!(
        assert_equiv(residual.as_ref(), &Expr::int(0), ctx).expect("assert_equiv"),
        "ODE residual not zero: {}",
        format_expr(residual.as_ref())
    );
}

fn desolve_sol(result: &ExprArc) -> &ExprArc {
    match result.as_ref() {
        Expr::Relation(RelOp::Eq, _, rhs) => rhs,
        other => panic!("expected desolve relation, got {}", format_expr(other)),
    }
}

pub fn subst_constants(sol: &ExprArc, samples: &[(&str, i64)]) -> ExprArc {
    let mut subs = HashMap::new();
    for (name, val) in samples {
        subs.insert(Ident::new(name), Expr::int(*val));
    }
    eval_subst_map(sol, &subs).expect("subst constants")
}

/// `desolve` returns `y(x)=sol`; check the linear ODE for sampled constant values.
pub fn assert_desolve_lin_ode(
    a2: &Expr,
    a1: &Expr,
    a0: &Expr,
    forcing: &Expr,
    result: &ExprArc,
    indep: &Ident,
    const_samples: &[(&str, i64)],
    ctx: &Context,
) {
    let sol = desolve_sol(result);
    let inst = subst_constants(sol, const_samples);
    assert_lin_ode_solution(a2, a1, a0, forcing, &inst, indep, ctx);
}
