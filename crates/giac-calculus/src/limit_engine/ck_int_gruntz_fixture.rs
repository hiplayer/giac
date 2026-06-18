//! Shared CK-INT-60/61 Gruntz limit shapes from upstream `check/testintegrate` L61–62.
//!
//! Decomposed as:
//! - **ratio** — `exp(inner)/exp(x) → 1` (Maxima `rtest_limit_gruntz.mac` L98)
//! - **CK-INT-60** — `exp(inner)/x → +infinity`
//! - **CK-INT-61** — `(exp(inner)-exp(x))/x → -exp(2)`

#[cfg(test)]
use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind};

/// `x*exp(-x) / (exp(-x) + exp(-2*x^2/(x+1)))`
#[cfg(test)]
pub(crate) fn inner() -> ExprArc {
    Arc::new(Expr::Frac(
        Expr::mul(vec![
            Expr::sym("x"),
            Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]),
        ]),
        Expr::add(vec![
            Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]),
            Expr::func(
                FuncKind::Exp,
                vec![Expr::mul(vec![
                    Expr::int(-2),
                    Arc::new(Expr::Frac(
                        Expr::pow(Expr::sym("x"), Expr::int(2)),
                        Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
                    )),
                ])],
            ),
        ]),
    ))
}

#[cfg(test)]
pub(crate) fn exp_inner() -> ExprArc {
    Expr::func(FuncKind::Exp, vec![inner()])
}

/// `limit(exp(inner)/exp(x), x, +infinity) = 1`
#[cfg(test)]
pub(crate) fn ratio() -> ExprArc {
    Arc::new(Expr::Frac(
        exp_inner(),
        Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
    ))
}

/// CK-INT-60: `limit(exp(inner)/x, x, +infinity) = +infinity`
#[cfg(test)]
pub(crate) fn ck_int_60() -> ExprArc {
    Expr::mul(vec![
        exp_inner(),
        Expr::pow(Expr::sym("x"), Expr::int(-1)),
    ])
}

/// CK-INT-61: `limit((exp(inner)-exp(x))/x, x, +infinity) = -exp(2)`
#[cfg(test)]
pub(crate) fn ck_int_61() -> ExprArc {
    Expr::mul(vec![
        Expr::add(vec![
            exp_inner(),
            Expr::mul(vec![Expr::int(-1), Expr::func(FuncKind::Exp, vec![Expr::sym("x")])]),
        ]),
        Expr::pow(Expr::sym("x"), Expr::int(-1)),
    ])
}

#[cfg(test)]
pub(crate) mod lines {
    pub const RATIO: &str =
        "limit(exp(x*exp(-x)/(exp(-x)+exp(-2*x^2/(x+1))))/exp(x),x,+infinity)";
    pub const CK_INT_60: &str =
        "limit(exp(x*exp(-x)/(exp(-x)+exp(-2*x^2/(x+1))))/x,x,+infinity)";
    pub const CK_INT_61: &str =
        "limit((exp(x*exp(-x)/(exp(-x)+exp(-2*x^2/(x+1))))-exp(x))/x,x,+infinity)";
}
