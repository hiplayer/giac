//! Quotient normalization for MRV limit preprocessing.
//!
//! **API 分层：** [`giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md)
//!
//! | 层级 | 内容 |
//! |------|------|
//! | **Pipeline** | `try_as_quotient`、`normalize_*`、`canonicalize_limit_entry` |
//! | **Pipeline private** | `normalize_inverse_products`、`normalize_inverse_with` |

use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind};
use num_traits::Signed;

/// **Pipeline** — 提取 `(numerator, denominator)`
pub(crate) fn try_as_quotient(expr: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    match expr.as_ref() {
        Expr::Frac(num, den) => Some((Arc::clone(num), Arc::clone(den))),
        Expr::Mul(factors) => {
            let mut num = Vec::new();
            let mut den = Vec::new();
            for f in factors {
                if let Expr::Pow(b, exp) = f.as_ref() {
                    if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) {
                        if let Expr::Int(n) = exp.as_ref() {
                            den.push(Expr::pow(Arc::clone(b), Arc::new(Expr::Int(-n))));
                            continue;
                        }
                    }
                }
                num.push(Arc::clone(f));
            }
            if den.is_empty() {
                if num.is_empty() {
                    return None;
                }
                let n = if num.len() == 1 {
                    Arc::clone(&num[0])
                } else {
                    Expr::mul(num)
                };
                return Some((n, Expr::int(1)));
            }
            let n = if num.is_empty() {
                Expr::int(1)
            } else if num.len() == 1 {
                Arc::clone(&num[0])
            } else {
                Expr::mul(num)
            };
            let d = if den.len() == 1 {
                Arc::clone(&den[0])
            } else {
                Expr::mul(den)
            };
            Some((n, d))
        }
        _ => None,
    }
}

/// **Pipeline** — 商式 Laurent 形态统一
pub(crate) fn normalize_expr_quotients(expr: &ExprArc) -> ExprArc {
    let normalized = match expr.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(normalize_expr_quotients).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(normalize_expr_quotients).collect()),
        Expr::Pow(b, e) => Expr::pow(normalize_expr_quotients(b), normalize_expr_quotients(e)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            normalize_expr_quotients(n),
            normalize_expr_quotients(d),
        )),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(normalize_expr_quotients).collect()),
        _ => Arc::clone(expr),
    };
    if let Some((n, d)) = try_as_quotient(&normalized) {
        Arc::new(Expr::Frac(n, d))
    } else {
        normalized
    }
}

/// **Pipeline** — 顶层 `Frac` → `Mul·den^-1`
pub(crate) fn unify_top_quotient(expr: &ExprArc) -> ExprArc {
    if let Expr::Frac(n, d) = expr.as_ref() {
        Expr::mul(vec![
            Arc::clone(n),
            Expr::pow(Arc::clone(d), Expr::int(-1)),
        ])
    } else {
        Arc::clone(expr)
    }
}

/// **Pipeline** — 极限入口规范形
pub(crate) fn canonicalize_limit_entry(expr: &ExprArc) -> ExprArc {
    unify_top_quotient(&normalize_inverse_products(expr))
}

/// **Pipeline** — `a*(b+c)^-1` → `a/(b+c)`
pub(crate) fn normalize_inverse_sums(expr: &ExprArc) -> ExprArc {
    normalize_inverse_with(expr, |b| matches!(b.as_ref(), Expr::Add(_)))
}

// **Pipeline private** — normalize inverse products
fn normalize_inverse_products(expr: &ExprArc) -> ExprArc {
    normalize_inverse_with(expr, |_| true)
}

// **Pipeline private** — normalize inverse with
fn normalize_inverse_with(
    expr: &ExprArc,
    accept_inverse_base: impl Fn(&ExprArc) -> bool + Copy,
) -> ExprArc {
    match expr.as_ref() {
        Expr::Add(ts) => Expr::add(
            ts.iter()
                .map(|t| normalize_inverse_with(t, accept_inverse_base))
                .collect(),
        ),
        Expr::Mul(fs) => {
            let parts: Vec<ExprArc> = fs
                .iter()
                .map(|f| normalize_inverse_with(f, accept_inverse_base))
                .collect();
            let mut num = Vec::new();
            let mut den = Vec::new();
            for f in &parts {
                if let Expr::Pow(b, e) = f.as_ref() {
                    if matches!(e.as_ref(), Expr::Int(n) if n.is_negative())
                        && accept_inverse_base(b)
                    {
                        if let Expr::Int(n) = e.as_ref() {
                            den.push(Arc::clone(b));
                            continue;
                        }
                    }
                }
                num.push(Arc::clone(f));
            }
            if den.is_empty() {
                return Expr::mul(parts);
            }
            let n = if num.is_empty() {
                Expr::int(1)
            } else if num.len() == 1 {
                num.into_iter().next().unwrap()
            } else {
                Expr::mul(num)
            };
            let d = if den.len() == 1 {
                den.into_iter().next().unwrap()
            } else {
                Expr::mul(den)
            };
            Arc::new(Expr::Frac(n, d))
        }
        Expr::Pow(b, e) => Expr::pow(
            normalize_inverse_with(b, accept_inverse_base),
            normalize_inverse_with(e, accept_inverse_base),
        ),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            normalize_inverse_with(n, accept_inverse_base),
            normalize_inverse_with(d, accept_inverse_base),
        )),
        Expr::Func(k, args) => Expr::func(
            *k,
            args.iter()
                .map(|a| normalize_inverse_with(a, accept_inverse_base))
                .collect(),
        ),
        _ => Arc::clone(expr),
    }
}
