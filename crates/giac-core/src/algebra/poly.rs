use std::sync::Arc;

pub use giac_poly::{Monomial, Poly, Var};
use giac_poly::PolyMod;
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::{EvalError, Expr, ExprArc, FuncKind, Ident};

use super::poly_conv::{
    rational_poly_reject, ERR_POLY_ALG_NO_ALG_COEFF, ERR_POLY_ALG_UNIMPL,
};

pub use super::poly_conv::expr_contains_alg_coeff;

fn var(id: &Ident) -> Var {
    Arc::from(id.as_str())
}

/// Convert an expression to a polynomial over **ℚ** in its variables.
///
/// **Stable** — path A in [expr-poly-conversion.md](../../../../.doc/expr-poly-conversion.md).
///
/// Coefficients must be rational (`Int` / `Rat`). Algebraic constants (`AlgExt`, `AlgExtC`,
/// concrete `rootof(...)`) are rejected with [`EvalError::TypeError`]; use
/// [`expr_contains_alg_coeff`] to detect them and [`poly_alg_from_expr`] for the future
/// `Poly<AlgExtC>` path.
pub fn expr_to_poly(expr: &Expr) -> Result<Poly, EvalError> {
    if expr_contains_alg_coeff(expr) {
        return Err(rational_poly_reject(expr));
    }
    expr_to_poly_inner(expr)
}

/// Lift an expression to a polynomial over algebraic coefficients (`Poly<AlgExtC>`).
///
/// **Stable (stub)** — path B; implementation blocked on P1-2/P1-4
/// ([GIAC-poly-algext-backlog.md](../../../../.doc/issues/GIAC-poly-algext-backlog.md)).
///
/// Returns `NotImplemented` when the expression contains algebraic coefficients.
/// Returns `TypeError` when the expression is purely rational (caller should use
/// [`expr_to_poly`] instead).
pub fn poly_alg_from_expr(expr: &Expr) -> Result<Poly, EvalError> {
    if !expr_contains_alg_coeff(expr) {
        return Err(EvalError::TypeError(ERR_POLY_ALG_NO_ALG_COEFF));
    }
    let _ = expr;
    Err(EvalError::NotImplemented(ERR_POLY_ALG_UNIMPL))
}

fn expr_to_poly_inner(expr: &Expr) -> Result<Poly, EvalError> {
    match expr {
        Expr::Int(n) => Ok(Poly::constant(Ratio::from_integer(n.clone()))),
        Expr::Rat(r) => Ok(Poly::constant(r.clone())),
        Expr::Symbol(id) => Ok(Poly::var(var(id))),
        Expr::Add(terms) => terms
            .iter()
            .map(|t| expr_to_poly_inner(t))
            .try_fold(Poly::zero(), |acc, p| Ok(acc.add(&p?))),
        Expr::Mul(factors) => factors
            .iter()
            .map(|f| expr_to_poly_inner(f))
            .try_fold(Poly::one(), |acc, p| Ok(acc.mul(&p?))),
        Expr::Pow(base, exp) => {
            let base_p = expr_to_poly_inner(base)?;
            if let Expr::Int(e) = exp.as_ref() {
                let e_u = crate::num_util::bigint_to_poly_exponent(e)?;
                return Ok(base_p.pow(e_u));
            }
            Err(EvalError::TypeError("non-polynomial power"))
        }
        _ => Err(EvalError::TypeError("not a polynomial expression")),
    }
}

/// **Stable** — `Poly` over ℚ → `Expr` (coefficients remain rational).
pub fn poly_to_expr(poly: &Poly) -> ExprArc {
    if poly.is_zero() {
        return Expr::int(0);
    }
    let mut terms: Vec<(u64, ExprArc)> = Vec::new();
    for (m, c) in &poly.terms {
        let coeff = ratio_to_expr(c);
        let term = monomial_to_expr(m, coeff);
        terms.push((m.degree(), term));
    }
    terms.sort_by(|a, b| b.0.cmp(&a.0));
    Expr::add(terms.into_iter().map(|(_, t)| t).collect())
}

/// Format a polynomial over ℤ/pℤ with per-coefficient `(c % p)` display (giac style).
pub fn poly_mod_to_expr(pm: &PolyMod) -> ExprArc {
    let modulus = pm
        .modulus
        .to_string()
        .parse::<i64>()
        .unwrap_or(0);
    if pm.is_zero() {
        return Arc::new(Expr::Mod(Expr::int(0), Expr::int(modulus)));
    }
    let mut terms: Vec<(u64, ExprArc)> = Vec::new();
    for (m, c) in &pm.terms {
        let rem = giac_poly::smod(
            c.val.to_string().parse().unwrap_or(0),
            modulus,
        );
        let coeff = Arc::new(Expr::Mod(Expr::int(rem), Expr::int(modulus)));
        let term = monomial_to_expr(m, coeff);
        terms.push((m.degree(), term));
    }
    terms.sort_by(|a, b| b.0.cmp(&a.0));
    Expr::add(terms.into_iter().map(|(_, t)| t).collect())
}

pub(crate) fn u64_to_expr_int(n: u64) -> ExprArc {
    match i64::try_from(n) {
        Ok(v) => Expr::int(v),
        Err(_) => Arc::new(Expr::Int(BigInt::from(n))),
    }
}

pub fn ratio_to_expr(r: &Ratio<BigInt>) -> ExprArc {
    if r.is_zero() {
        Expr::int(0)
    } else if r.denom() == &BigInt::one() {
        if let Ok(v) = r.numer().to_string().parse::<i64>() {
            Expr::int(v)
        } else {
            Arc::new(Expr::Rat(r.clone()))
        }
    } else {
        Arc::new(Expr::Rat(r.clone()))
    }
}

fn monomial_to_expr(m: &Monomial, coeff: ExprArc) -> ExprArc {
    if m.is_const() {
        return coeff;
    }
    let mut factors: Vec<ExprArc> = Vec::new();
    if !matches!(coeff.as_ref(), Expr::Int(n) if n.is_one()) {
        factors.push(coeff);
    }
    for (v, e) in m.iter() {
        let base = Expr::sym(v.as_ref());
        if e == 1 {
            factors.push(base);
        } else {
            factors.push(Expr::pow(base, u64_to_expr_int(e)));
        }
    }
    if factors.is_empty() {
        Expr::int(1)
    } else if factors.len() == 1 {
        Arc::clone(&factors[0])
    } else {
        Expr::mul(factors)
    }
}

pub fn vars_from_expr(expr: &Expr) -> Vec<Var> {
    let mut out = Vec::new();
    collect_vars(expr, &mut out);
    out.sort();
    out.dedup();
    out
}

fn collect_vars(expr: &Expr, out: &mut Vec<Var>) {
    match expr {
        Expr::Symbol(id) => out.push(var(id)),
        Expr::Add(ts) | Expr::Mul(ts) => ts.iter().for_each(|t| collect_vars(t, out)),
        Expr::Pow(b, e) => {
            collect_vars(b, out);
            collect_vars(e, out);
        }
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::poly_conv::{ERR_ALG_EXT_COEFF, ERR_POLY_ALG_NO_ALG_COEFF, ERR_POLY_ALG_UNIMPL, ERR_ROOTOF_COEFF};
    use crate::AlgExtData;

    #[test]
    fn gcd_linear_bridge() {
        let x = Ident::new("x");
        let p = expr_to_poly(
            &Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::int(-1),
            ]),
        )
        .unwrap();
        let q = expr_to_poly(
            &Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(3)),
                Expr::int(-1),
            ]),
        )
        .unwrap();
        let g = p.gcd(&q);
        assert_eq!(poly_to_expr(&g), Expr::add(vec![Expr::sym("x"), Expr::int(-1)]));
        let _ = x;
    }

    #[test]
    fn expr_to_poly_matches_expanded_power() {
        let expanded = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(4)),
            Expr::mul(vec![Expr::int(12), Expr::pow(Expr::sym("x"), Expr::int(3))]),
            Expr::mul(vec![Expr::int(54), Expr::pow(Expr::sym("x"), Expr::int(2))]),
            Expr::mul(vec![Expr::int(108), Expr::sym("x")]),
            Expr::int(81),
        ]);
        let p = expr_to_poly(&expanded).unwrap();
        let x = Poly::var("x");
        let expected = x
            .add(&Poly::constant(Ratio::from_integer(BigInt::from(3))))
            .pow(4);
        assert_eq!(p, expected);
        assert_eq!(giac_poly::factor_poly(&p), expected);
    }

    #[test]
    fn expr_to_poly_high_degree_within_limit() {
        let e = Expr::pow(Expr::sym("x"), Expr::int(100));
        assert!(expr_to_poly(&e).is_ok());
    }

    #[test]
    fn expr_to_poly_rejects_algext() {
        let min = Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ));
        let e = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap()
        .into_expr();
        assert!(matches!(
            expr_to_poly(e.as_ref()),
            Err(EvalError::TypeError(ERR_ALG_EXT_COEFF))
        ));
    }

    #[test]
    fn expr_to_poly_rejects_nested_rootof() {
        let min = Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ));
        let e = Expr::add(vec![
            Expr::sym("x"),
            Expr::func(
                FuncKind::RootOf,
                vec![
                    Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
                    min,
                ],
            ),
        ]);
        assert!(matches!(
            expr_to_poly(&e),
            Err(EvalError::TypeError(ERR_ROOTOF_COEFF))
        ));
    }

    #[test]
    fn poly_alg_from_expr_rejects_rational() {
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::int(1),
        ]);
        assert!(matches!(
            poly_alg_from_expr(&e),
            Err(EvalError::TypeError(ERR_POLY_ALG_NO_ALG_COEFF))
        ));
    }

    #[test]
    fn poly_alg_from_expr_stub_for_alg_coeff() {
        let min = Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ));
        let e = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap()
        .into_expr();
        assert!(matches!(
            poly_alg_from_expr(e.as_ref()),
            Err(EvalError::NotImplemented(ERR_POLY_ALG_UNIMPL))
        ));
    }
}
