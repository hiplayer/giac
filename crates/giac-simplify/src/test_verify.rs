//! Semantic verification for giac-simplify tests (`.doc/conformance-testing.md` §3).
//! See `.doc/issues/GIAC-expr-api-test-contains-cleanup.md` H4.
#![allow(clippy::expect_used)] // test helper: expect carries assert context

use std::sync::Arc;

use giac_core::{eval, format_expr, Context, Expr, ExprArc, FuncKind};

use crate::{assert_equiv, expand};

/// `expand(factors)` is mathematically equal to `orig`.
pub fn assert_factorization(orig: &Expr, factored: &Expr, ctx: &Context) {
    let expanded = expand(factored, ctx).expect("expand factors");
    assert!(
        assert_equiv(orig, expanded.as_ref(), ctx).expect("assert_equiv"),
        "factorization mismatch: expand got {}",
        format_expr(expanded.as_ref())
    );
}

/// `expand(factored)` via eval(Expand) — for plugin eval(Factor) results.
pub fn assert_factorization_eval(orig: &ExprArc, factored: &ExprArc, ctx: &Context) {
    let expanded = eval(
        Expr::func(FuncKind::Expand, vec![Arc::clone(factored)]).as_ref(),
        ctx,
    )
    .expect("eval(Expand)");
    assert!(
        assert_equiv(orig.as_ref(), expanded.as_ref(), ctx).expect("assert_equiv"),
        "factorization mismatch: expand got {}",
        format_expr(expanded.as_ref())
    );
}
