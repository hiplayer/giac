//! Semantic verification for giac-solve unit tests (see `.doc/conformance-testing.md` §3).
//! Test tier labels on each `#[test]`: `.doc/test-writing-spec.md` · audit §2.
#![allow(clippy::expect_used)] // test helper: expect carries assert context

use std::sync::Arc;

use giac_core::{
    eval, try_as_algext_data, Context, EvalError, Expr, ExprArc, FuncKind, Ident, RelOp,
};
use giac_simplify::{assert_equiv, is_zero};

pub fn list_items(e: &ExprArc) -> &[ExprArc] {
    match e.as_ref() {
        Expr::List(v) => v,
        other => panic!("expected list, got {other:?}"),
    }
}

pub fn eval_at(
    expr: &ExprArc,
    var: &Ident,
    value: &ExprArc,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let mut c = ctx.clone();
    c.set(var.clone(), Arc::clone(value));
    eval(expr, &c)
}

pub fn assert_roots_zero_poly(poly: &ExprArc, var: &Ident, roots: &[ExprArc], ctx: &Context) {
    for r in roots {
        let v = eval_at(poly, var, r, ctx).expect("eval poly at root");
        assert!(
            is_zero(v.as_ref(), ctx).expect("is_zero"),
            "root does not zero polynomial"
        );
    }
}

pub fn assert_is_algext_or_rootof(e: &ExprArc) {
    let ok = try_as_algext_data(e.as_ref()).is_some()
        || matches!(e.as_ref(), Expr::AlgExt(_))
        || matches!(e.as_ref(), Expr::Func(FuncKind::RootOf, _));
    assert!(ok, "expected AlgExt/rootof, got {e:?}");
}

pub fn assert_equation_solutions(eq: &ExprArc, var: &Ident, solutions: &ExprArc, ctx: &Context) {
    let (lhs, rhs) = match eq.as_ref() {
        Expr::Relation(RelOp::Eq, l, r) => (l, r),
        other => panic!("expected equation, got {other:?}"),
    };
    for sol in list_items(solutions) {
        let lv = eval_at(lhs, var, sol, ctx).expect("eval lhs");
        let rv = eval_at(rhs, var, sol, ctx).expect("eval rhs");
        assert!(
            assert_equiv(lv.as_ref(), rv.as_ref(), ctx).expect("assert_equiv"),
            "solution does not satisfy equation"
        );
    }
}

#[allow(dead_code)] // ponytail: reserved for list golden migration
pub fn assert_list_has_equiv(expected: &Expr, items: &[ExprArc], ctx: &Context) {
    assert!(
        items
            .iter()
            .any(|it| assert_equiv(expected, it.as_ref(), ctx).unwrap()),
        "no entry equivalent to {expected:?}"
    );
}

/// `froot` / `froots` flat list `[r1,m1,r2,m2,…]`.
pub fn froot_pairs(items: &[ExprArc]) -> Vec<(&ExprArc, &ExprArc)> {
    assert!(items.len() % 2 == 0, "froot list length must be even");
    items.chunks(2).map(|c| (&c[0], &c[1])).collect()
}

pub fn assert_froot_has_root(items: &[ExprArc], root: &Expr, ctx: &Context) {
    let found = froot_pairs(items)
        .iter()
        .any(|(r, _)| assert_equiv(root, r.as_ref(), ctx).unwrap());
    assert!(found, "froot missing root {root:?}");
}

/// `realroot` list of `[root, multiplicity]`.
pub fn realroot_entries(items: &[ExprArc]) -> Vec<(&ExprArc, &ExprArc)> {
    items
        .iter()
        .map(|pair| match pair.as_ref() {
            Expr::List(v) if v.len() == 2 => (&v[0], &v[1]),
            other => panic!("expected [root,mult] pair, got {other:?}"),
        })
        .collect()
}

pub fn assert_realroot_has(items: &[ExprArc], root: &Expr, ctx: &Context) {
    let found = realroot_entries(items)
        .iter()
        .any(|(r, _)| assert_equiv(root, r.as_ref(), ctx).unwrap());
    assert!(found, "realroot missing {root:?}");
}
