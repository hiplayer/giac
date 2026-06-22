//! Shared `limit_engine` expression predicates.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Pipeline private** | `is_expr_zero`, `is_half_exponent` |

use giac_core::{Expr, ExprArc};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

/// Whether `e` is the integer constant zero.
pub(crate) fn is_expr_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

/// Whether `exp` is the rational exponent `1/2`.
pub(crate) fn is_half_exponent(exp: &ExprArc) -> bool {
    matches!(exp.as_ref(), Expr::Rat(r) if *r == Ratio::new(1.into(), 2.into()))
        || matches!(
            exp.as_ref(),
            Expr::Frac(n, d)
                if matches!(n.as_ref(), Expr::Int(nn) if nn.is_one())
                    && matches!(d.as_ref(), Expr::Int(dd) if dd == &BigInt::from(2))
        )
}
