//! Expr ↔ `Poly` conversion boundary (B-01 / P0-C).
//!
//! Normative contract: [expr-poly-conversion.md](../../../../.doc/expr-poly-conversion.md)

use crate::error::EvalError;
use crate::expr::{Expr, FuncKind};

/// `expr_to_poly`: `Expr::AlgExt` coefficient.
pub const ERR_ALG_EXT_COEFF: &str = "alg ext not allowed in rational polynomial";
/// `expr_to_poly`: `Expr::AlgExtC` coefficient.
pub const ERR_ALG_EXT_C_COEFF: &str = "alg extc not allowed in rational polynomial";
/// `expr_to_poly`: concrete `rootof(...)` subtree.
pub const ERR_ROOTOF_COEFF: &str = "rootof not allowed in rational polynomial";
/// `poly_alg_from_expr`: lift path blocked until `Poly<AlgExtC>` (P1).
pub const ERR_POLY_ALG_UNIMPL: &str =
    "poly_alg_from_expr requires Poly<AlgExtC> (P1); see .doc/expr-poly-conversion.md";
/// `poly_alg_from_expr`: expression has no algebraic coefficients.
pub const ERR_POLY_ALG_NO_ALG_COEFF: &str =
    "expression has no algebraic coefficients; use expr_to_poly for Poly over Q";

/// Whether `expr` contains an algebraic extension constant usable as a polynomial
/// coefficient (evaluated `AlgExt` / `AlgExtC`, or concrete `rootof(...)`).
///
/// **Stable** — predicate for callers choosing `expr_to_poly` vs `poly_alg_from_expr`.
pub fn expr_contains_alg_coeff(expr: &Expr) -> bool {
    match expr {
        Expr::AlgExt(_) | Expr::AlgExtC(_) => true,
        Expr::Func(FuncKind::RootOf, _) => true,
        Expr::Add(ts) | Expr::Mul(ts) | Expr::Seq(ts) | Expr::List(ts) => {
            ts.iter().any(|t| expr_contains_alg_coeff(t.as_ref()))
        }
        Expr::Frac(n, d) | Expr::Pow(n, d) | Expr::Complex(n, d) | Expr::Mod(n, d) => {
            expr_contains_alg_coeff(n.as_ref()) || expr_contains_alg_coeff(d.as_ref())
        }
        Expr::Func(_, args) => args.iter().any(|a| expr_contains_alg_coeff(a.as_ref())),
        Expr::Int(_) | Expr::Rat(_) | Expr::Symbol(_) | Expr::Str(_) | Expr::Undefined => false,
        Expr::Matrix(rows) | Expr::GiacMatrix(rows) => rows
            .iter()
            .flatten()
            .any(|e| expr_contains_alg_coeff(e.as_ref())),
        Expr::Relation(_, lhs, rhs) => {
            expr_contains_alg_coeff(lhs.as_ref()) || expr_contains_alg_coeff(rhs.as_ref())
        }
    }
}

/// Error for [`super::poly::expr_to_poly`] when [`expr_contains_alg_coeff`] holds.
pub(crate) fn rational_poly_reject(expr: &Expr) -> EvalError {
    EvalError::TypeError(first_alg_coeff_message(expr))
}

fn first_alg_coeff_message(expr: &Expr) -> &'static str {
    match expr {
        Expr::AlgExt(_) => ERR_ALG_EXT_COEFF,
        Expr::AlgExtC(_) => ERR_ALG_EXT_C_COEFF,
        Expr::Func(FuncKind::RootOf, _) => ERR_ROOTOF_COEFF,
        Expr::Add(ts) | Expr::Mul(ts) | Expr::Seq(ts) | Expr::List(ts) => ts
            .iter()
            .find_map(|t| {
                if expr_contains_alg_coeff(t.as_ref()) {
                    Some(first_alg_coeff_message(t.as_ref()))
                } else {
                    None
                }
            })
            .unwrap_or(ERR_ALG_EXT_COEFF),
        Expr::Frac(n, d) | Expr::Pow(n, d) | Expr::Complex(n, d) | Expr::Mod(n, d) => {
            if expr_contains_alg_coeff(n.as_ref()) {
                first_alg_coeff_message(n.as_ref())
            } else {
                first_alg_coeff_message(d.as_ref())
            }
        }
        Expr::Func(_, args) => args
            .iter()
            .find_map(|a| {
                if expr_contains_alg_coeff(a.as_ref()) {
                    Some(first_alg_coeff_message(a.as_ref()))
                } else {
                    None
                }
            })
            .unwrap_or(ERR_ROOTOF_COEFF),
        Expr::Matrix(rows) | Expr::GiacMatrix(rows) => rows
            .iter()
            .flatten()
            .find_map(|e| {
                if expr_contains_alg_coeff(e.as_ref()) {
                    Some(first_alg_coeff_message(e.as_ref()))
                } else {
                    None
                }
            })
            .unwrap_or(ERR_ALG_EXT_COEFF),
        Expr::Relation(_, lhs, rhs) => {
            if expr_contains_alg_coeff(lhs.as_ref()) {
                first_alg_coeff_message(lhs.as_ref())
            } else {
                first_alg_coeff_message(rhs.as_ref())
            }
        }
        _ => ERR_ALG_EXT_COEFF,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use crate::AlgExtData;

    #[test]
    fn detects_rootof_in_sum() {
        let min = Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ));
        let root = Expr::func(
            FuncKind::RootOf,
            vec![
                Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
                min,
            ],
        );
        let e = Expr::add(vec![Expr::sym("x"), root]);
        assert!(expr_contains_alg_coeff(&e));
        assert_eq!(
            rational_poly_reject(&e),
            EvalError::TypeError(ERR_ROOTOF_COEFF)
        );
    }

    #[test]
    fn detects_algext_in_product() {
        let min = Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ));
        let alpha = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap()
        .into_expr();
        let e = Expr::mul(vec![Expr::sym("x"), alpha]);
        assert!(expr_contains_alg_coeff(e.as_ref()));
    }

    #[test]
    fn rational_poly_has_no_alg_coeff() {
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::rat(1, 3),
        ]);
        assert!(!expr_contains_alg_coeff(&e));
    }
}
