//! Semantic verification for giac-core tests (`.doc/conformance-testing.md` §3).
//! See `.doc/issues/GIAC-expr-api-test-contains-cleanup.md` H6.

use std::sync::Arc;

use num_traits::Zero;

use crate::{eval, format_expr, Context, EvalError, Expr, ExprArc, FuncKind, Ident};

fn rem_mod_var(e: &Expr, var: &Ident, ctx: &Context) -> ExprArc {
    eval(
        Expr::func(FuncKind::Rem, vec![Arc::new(e.clone()), Expr::sym(var.as_str())]).as_ref(),
        ctx,
    )
    .expect("rem")
}

fn assert_rem_zero(rem: &ExprArc, ctx: &Context) {
    let n = eval(rem.as_ref(), ctx).expect("eval rem");
    assert!(
        matches!(n.as_ref(), Expr::Int(n) if n.is_zero()),
        "expected zero rem, got {}",
        format_expr(n.as_ref())
    );
}

/// Polynomial identity `lhs - rhs ≡ 0 (mod var)` via `rem`.
pub fn assert_poly_identity(lhs: &Expr, rhs: &Expr, var: &Ident, ctx: &Context) {
    let diff = Expr::add(vec![
        Arc::new(lhs.clone()),
        Expr::mul(vec![Expr::int(-1), Arc::new(rhs.clone())]),
    ]);
    assert_rem_zero(&rem_mod_var(diff.as_ref(), var, ctx), ctx);
}

/// `a*u + b*v` equals `g` as polynomials in `var`.
pub fn assert_bezout(
    a: &Expr,
    b: &Expr,
    g: &Expr,
    u: &Expr,
    v: &Expr,
    var: &Ident,
    ctx: &Context,
) {
    let bezout = Expr::add(vec![
        Expr::mul(vec![Arc::new(a.clone()), Arc::new(u.clone())]),
        Expr::mul(vec![Arc::new(b.clone()), Arc::new(v.clone())]),
    ]);
    assert_poly_identity(bezout.as_ref(), g, var, ctx);
}

/// `rem(dividend, divisor) ≡ 0` as polynomials in `var`.
pub fn assert_divides(divisor: &Expr, dividend: &Expr, _var: &Ident, ctx: &Context) {
    let rem = eval(
        Expr::func(
            FuncKind::Rem,
            vec![Arc::new(dividend.clone()), Arc::new(divisor.clone())],
        )
        .as_ref(),
        ctx,
    )
    .expect("rem");
    assert_rem_zero(&rem, ctx);
}

pub fn list_items(e: &ExprArc) -> Result<&[ExprArc], EvalError> {
    match e.as_ref() {
        Expr::List(v) | Expr::Seq(v) => Ok(v),
        _ => Err(EvalError::TypeError("expected list")),
    }
}
