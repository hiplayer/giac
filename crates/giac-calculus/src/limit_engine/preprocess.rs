//! Limit preprocessing before MRV (GIAC-216d / `limit_symbolic_preprocess` subset).

use std::sync::Arc;

use giac_core::{eval, Context, EvalError, Expr, ExprArc, Ident};

use crate::risch::pow2expln;

/// `pow2expln` and light normalization before series / limit asymptotics.
pub(crate) fn limit_preprocess_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    series_preprocess(expr, var, ctx)
}

pub(crate) fn series_preprocess(
    expr: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let normalized = pow2expln(expr, var);
    eval(normalized.as_ref(), ctx)
}
