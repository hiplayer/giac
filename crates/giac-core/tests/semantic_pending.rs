//! Semantic tests blocked on `assert_equiv` / canonical `normal` (see issue §8 B-*).
//! Active crate unit tests keep golden/smoke; these are the target contracts.

use giac_core::{eval, format_expr, Context, Expr, FuncKind};
use giac_simplify::{assert_equiv, install_simplify};

fn ctx() -> Context {
    let mut ctx = Context::xcas_default();
    install_simplify(&mut ctx);
    ctx
}

#[test]
#[ignore = "B-MUL-FMT: (1+i)*x 应 assert_equiv x+x*i"]
fn eval_mul_mixed_complex_symbolic_equiv() {
    let ctx = ctx();
    let i = Expr::sym("i");
    let r = eval(
        Expr::mul(vec![
            Expr::add(vec![Expr::int(1), i.clone()]),
            Expr::sym("x"),
        ])
        .as_ref(),
        &ctx,
    )
    .unwrap();
    let expected = Expr::add(vec![
        Expr::sym("x"),
        Expr::mul(vec![i, Expr::sym("x")]),
    ]);
    assert!(
        assert_equiv(r.as_ref(), expected.as_ref(), &ctx).unwrap(),
        "got {}",
        format_expr(r.as_ref())
    );
}

#[test]
#[ignore = "B-ARG: arg(-1+i) 应 assert_equiv 3*pi/4"]
fn eval_arg_second_quadrant_canonical() {
    let ctx = ctx();
    let i = Expr::sym("i");
    let r = eval(
        Expr::func(
            FuncKind::Arg,
            vec![Expr::add(vec![Expr::int(-1), i])],
        )
        .as_ref(),
        &ctx,
    )
    .unwrap();
    let expected = Expr::mul(vec![Expr::rat(3, 4), Expr::sym("pi")]);
    assert!(
        assert_equiv(r.as_ref(), expected.as_ref(), &ctx).unwrap(),
        "got {}",
        format_expr(r.as_ref())
    );
}
