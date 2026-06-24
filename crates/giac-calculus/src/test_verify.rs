//! Semantic verification for giac-calculus unit tests (`.doc/conformance-testing.md` §3).
//! See `.doc/issues/GIAC-expr-api-test-contains-cleanup.md` H1.
#![allow(clippy::expect_used)] // test helper: expect carries assert context

use std::sync::Arc;

use giac_core::{eval, format_expr, Context, Expr, ExprArc, FuncKind, Ident, RelOp};
use giac_simplify::assert_equiv;

use crate::diff::diff;

/// `diff(F, var)` is mathematically equal to `integrand` under `ctx`.
pub fn assert_deriv_equals_integrand(
    integrand: &ExprArc,
    var: &Ident,
    antiderivative: &ExprArc,
    ctx: &Context,
) {
    let d = diff(antiderivative, var).expect("diff");
    assert!(
        assert_equiv(integrand.as_ref(), d.as_ref(), ctx).expect("assert_equiv"),
        "diff(F) != f: got diff {}",
        format_expr(d.as_ref())
    );
}

/// `series(f, var, center, order)` matches `expected` (truncated polynomial).
pub fn assert_series_equiv_at(
    f: &ExprArc,
    var: &Ident,
    center: &ExprArc,
    order: i64,
    expected: &Expr,
    ctx: &Context,
) {
    let got = eval(
        Expr::func(
            FuncKind::Series,
            vec![
                Arc::clone(f),
                Expr::sym(var.as_str()),
                Arc::clone(center),
                Expr::int(order),
            ],
        )
        .as_ref(),
        ctx,
    )
    .expect("eval(Series)");
    assert!(
        assert_equiv(got.as_ref(), expected, ctx).expect("assert_equiv"),
        "series mismatch: got {}",
        format_expr(got.as_ref())
    );
}

/// `taylor(f, var=x, center, order)` matches `expected`.
pub fn assert_taylor_equiv_at(
    f: &ExprArc,
    var: &Ident,
    center: &ExprArc,
    order: i64,
    expected: &Expr,
    ctx: &Context,
) {
    let got = eval(
        Expr::func(
            FuncKind::Taylor,
            vec![
                Arc::clone(f),
                Arc::new(Expr::Relation(
                    RelOp::Eq,
                    Expr::sym(var.as_str()),
                    Arc::clone(center),
                )),
                Expr::int(order),
            ],
        )
        .as_ref(),
        ctx,
    )
    .expect("eval(Taylor)");
    assert!(
        assert_equiv(got.as_ref(), expected, ctx).expect("assert_equiv"),
        "taylor mismatch: got {}",
        format_expr(got.as_ref())
    );
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::Expr;

    use crate::plugin::xcas_default;

    use super::*;

    #[test]
    fn assert_deriv_equals_integrand_smoke() {
        let ctx = xcas_default();
        let x = Ident::new("x");
        let integrand = Arc::new(Expr::mul(vec![Expr::int(2), Expr::sym("x")]));
        let antiderivative = Arc::new(Expr::pow(Expr::sym("x"), Expr::int(2)));
        assert_deriv_equals_integrand(&integrand, &x, &antiderivative, &ctx);
    }
}
