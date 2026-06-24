//! K-route `partfrac` / `eval_partfrac` — no literal golden (rootof display varies).
//!
//! Exit criterion: `expand(partfrac(f,x)) ≡ f` via `assert_equiv` (conformance §3).

use std::sync::Arc;

use giac_core::{eval, Context, Expr, FuncKind};
use giac_simplify::{assert_equiv, install_simplify};

#[test]
fn eval_partfrac_one_over_x_squared_minus_two() {
    let mut ctx = Context::xcas_default();
    install_simplify(&mut ctx);
    let orig = Expr::pow(
        Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::int(-2),
        ]),
        Expr::int(-1),
    );
    let pf = eval(
        Expr::func(
            FuncKind::Partfrac,
            vec![Arc::clone(&orig), Expr::sym("x")],
        )
        .as_ref(),
        &ctx,
    )
    .unwrap();
    let expanded = eval(
        Expr::func(FuncKind::Expand, vec![pf]).as_ref(),
        &ctx,
    )
    .unwrap();
    assert_equiv(expanded.as_ref(), orig.as_ref(), &ctx).unwrap();
}
