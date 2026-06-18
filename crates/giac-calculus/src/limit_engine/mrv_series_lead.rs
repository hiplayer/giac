//! Quotient normalization for MRV limit preprocessing.

use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind};
use num_traits::Signed;

/// Extract `(numerator, denominator)` from `a/b` or `a * b^-1`.
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

/// `a/b` and `a*b^-1` share the same Laurent leading term.
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

/// Top-level `Frac` → `Mul·den^-1` for exp-diff preprocess (parse emits `Frac`, fixtures use `Mul`).
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

/// `a*(b+c+…)^-1` → `a/(b+c+…)` so merge/fold see fixture-shaped denominators.
pub(crate) fn normalize_inverse_sums(expr: &ExprArc) -> ExprArc {
    match expr.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(normalize_inverse_sums).collect()),
        Expr::Mul(fs) => {
            let parts: Vec<ExprArc> = fs.iter().map(normalize_inverse_sums).collect();
            let mut num = Vec::new();
            let mut den = None;
            for f in &parts {
                if let Expr::Pow(b, e) = f.as_ref() {
                    if matches!(e.as_ref(), Expr::Int(n) if n.is_negative())
                        && matches!(b.as_ref(), Expr::Add(_))
                    {
                        den = Some(Arc::clone(b));
                        continue;
                    }
                }
                num.push(Arc::clone(f));
            }
            if let Some(d) = den {
                let n = if num.is_empty() {
                    Expr::int(1)
                } else if num.len() == 1 {
                    num.into_iter().next().unwrap()
                } else {
                    Expr::mul(num)
                };
                return Arc::new(Expr::Frac(n, d));
            }
            Expr::mul(parts)
        }
        Expr::Pow(b, e) => Expr::pow(normalize_inverse_sums(b), normalize_inverse_sums(e)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            normalize_inverse_sums(n),
            normalize_inverse_sums(d),
        )),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(normalize_inverse_sums).collect()),
        _ => Arc::clone(expr),
    }
}
