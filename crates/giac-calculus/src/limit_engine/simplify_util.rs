//! Shared normalization helpers for limit_engine (no `eval` fallback).
//!
//! **API 分层：** [`giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md)
//!
//! | 层级 | 内容 |
//! |------|------|
//! | **Pipeline** | `simplify_limit_expr`、`int_pow_growth_sub_rank` |
//! | **Pipeline private** | （本模块无私有 `fn`） |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{bigint_to_i64, Context, Expr, ExprArc};
use giac_simplify::{normal, ratnormal};
use num_traits::Signed;

/// **Pipeline** — `ratnormal` 后 `normal`；失败保留输入
pub(crate) fn simplify_limit_expr(expr: &ExprArc, ctx: &Context) -> ExprArc {
    let rat = match ratnormal(expr.as_ref(), ctx) {
        Ok(r) => r,
        Err(_) => Arc::clone(expr),
    };
    normal(rat.as_ref(), ctx).unwrap_or(rat)
}

/// **Stable** — 符号负常数谓词
pub(crate) fn is_negative_const_expr(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Int(n) => n.is_negative(),
        Expr::Rat(r) => r.is_negative(),
        Expr::Frac(num, _) => matches!(num.as_ref(), Expr::Int(n) if n.is_negative()),
        Expr::Mul(fs) => {
            let mut neg = false;
            for f in fs {
                if matches!(f.as_ref(), Expr::Int(n) if n.is_negative()) {
                    neg = !neg;
                }
            }
            neg
        }
        _ => false,
    }
}

/// **Pipeline** — 整数底 `a^var` 的精确 sub-rank
pub(crate) fn int_pow_growth_sub_rank(base: &ExprArc) -> Option<i64> {
    let n = match base.as_ref() {
        Expr::Int(i) => bigint_to_i64(i).ok()?,
        _ => return None,
    };
    if n <= 0 {
        return None;
    }
    Some(n.ilog2() as i64)
}
