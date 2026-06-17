//! Quotient normalization for MRV limit preprocessing.

use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind};
use num_traits::Signed;

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
    if let Some((n, d)) = extract_quotient(&normalized) {
        Arc::new(Expr::Frac(n, d))
    } else {
        normalized
    }
}

fn extract_quotient(expr: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    match expr.as_ref() {
        Expr::Frac(n, d) => Some((Arc::clone(n), Arc::clone(d))),
        Expr::Mul(fs) => {
            let mut num = Vec::new();
            let mut den = None;
            for f in fs {
                if let Expr::Pow(b, exp) = f.as_ref() {
                    if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) {
                        den = Some(Arc::clone(b));
                        continue;
                    }
                }
                num.push(Arc::clone(f));
            }
            den.map(|d| (if num.is_empty() { Expr::int(1) } else { Expr::mul(num) }, d))
        }
        _ => None,
    }
}
