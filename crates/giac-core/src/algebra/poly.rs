//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::sync::Arc;

pub use giac_poly::{Monomial, Poly, PolyCoeff, Var};
use giac_poly::PolyMod;
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::{EvalError, Expr, ExprArc, FuncKind, Ident};

use super::alg_ext_c::canonicalize_to_algext_c;
use super::poly_alg_coeff::AlgExtCPolyCoeff;
use super::poly_conv::{
    rational_poly_reject, ERR_POLY_ALG_NO_ALG_COEFF,
};

pub use super::poly_conv::expr_contains_alg_coeff;

/// Polynomial over algebraic coefficients (`Poly<AlgExtCPolyCoeff>`).
pub type PolyAlgExt = Poly<AlgExtCPolyCoeff>;

// **Stable** — Poly univariate generator
fn var(id: &Ident) -> Var {
    Arc::from(id.as_str())
}

/// Ring operations for shared [`poly_from_expr_shape`] (H1).
struct PolyExprRing<C: PolyCoeff> {
    zero: fn() -> Poly<C>,
    one: fn() -> Poly<C>,
    var: fn(&Ident) -> Poly<C>,
    add: fn(&Poly<C>, &Poly<C>) -> Result<Poly<C>, EvalError>,
    mul: fn(&Poly<C>, &Poly<C>) -> Result<Poly<C>, EvalError>,
    pow: fn(&Poly<C>, u64) -> Result<Poly<C>, EvalError>,
}

/// Shared polynomial-shaped Expr descent (path A and path B).
// **Pipeline private** — `poly_from_expr_shape`
fn poly_from_expr_shape<C: PolyCoeff>(
    expr: &Expr,
    leaf: fn(&Expr) -> Result<C, EvalError>,
    ring: &PolyExprRing<C>,
) -> Result<Poly<C>, EvalError> {
    match expr {
        Expr::Int(n) => Ok(Poly::ring_constant(leaf(&Expr::Int(n.clone()))?)),
        Expr::Rat(r) => Ok(Poly::ring_constant(leaf(&Expr::Rat(r.clone()))?)),
        Expr::AlgExt(_) | Expr::AlgExtC(_) | Expr::Func(FuncKind::RootOf, _) => {
            Ok(Poly::ring_constant(leaf(expr)?))
        }
        Expr::Symbol(id) => Ok((ring.var)(id)),
        Expr::Add(terms) => terms
            .iter()
            .map(|t| poly_from_expr_shape(t, leaf, ring))
            .try_fold((ring.zero)(), |acc, p| (ring.add)(&acc, &p?)),
        Expr::Mul(factors) => factors
            .iter()
            .map(|f| poly_from_expr_shape(f, leaf, ring))
            .try_fold((ring.one)(), |acc, p| (ring.mul)(&acc, &p?)),
        Expr::Pow(base, exp) => {
            let base_p = poly_from_expr_shape(base, leaf, ring)?;
            if let Expr::Int(e) = exp.as_ref() {
                let e_u = crate::num_util::bigint_to_poly_exponent(e)?;
                return (ring.pow)(&base_p, e_u);
            }
            Err(EvalError::TypeError("non-polynomial power"))
        }
        _ => Err(EvalError::TypeError("not a polynomial expression")),
    }
}

// **Pipeline private** — `rational_leaf`
fn rational_leaf(expr: &Expr) -> Result<Ratio<BigInt>, EvalError> {
    match expr {
        Expr::Int(n) => Ok(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Ok(r.clone()),
        _ => Err(EvalError::TypeError("not a rational polynomial coefficient")),
    }
}

// **Pipeline private** — `rational_poly_ring`
fn rational_poly_ring() -> PolyExprRing<Ratio<BigInt>> {
    PolyExprRing {
        zero: Poly::zero,
        one: Poly::one,
        var: |id| Poly::var(var(id)),
        add: |a, b| Ok(a.add(b)),
        mul: |a, b| Ok(a.mul(b)),
        pow: |p, e| Ok(p.pow(e)),
    }
}

// **Pipeline private** — `algext_leaf`
fn algext_leaf(expr: &Expr) -> Result<AlgExtCPolyCoeff, EvalError> {
    match expr {
        Expr::Int(_) | Expr::Rat(_) | Expr::AlgExt(_) | Expr::AlgExtC(_)
        | Expr::Func(FuncKind::RootOf, _) => {
            Ok(AlgExtCPolyCoeff(canonicalize_to_algext_c(expr)?))
        }
        _ => Err(EvalError::TypeError("not a polynomial coefficient")),
    }
}

// **Pipeline private** — `algext_poly_ring`
fn algext_poly_ring() -> PolyExprRing<AlgExtCPolyCoeff> {
    PolyExprRing {
        zero: PolyAlgExt::ring_zero,
        one: PolyAlgExt::ring_one,
        var: |id| PolyAlgExt::ring_var(var(id)),
        add: |a, b| a.try_add(b).map_err(Into::into),
        mul: |a, b| a.try_mul(b).map_err(Into::into),
        pow: |p, e| p.try_pow(e).map_err(Into::into),
    }
}

/// Convert an expression to a polynomial over **ℚ** in its variables.
///
/// **Stable** — path A in [expr-poly-conversion.md](../../../../.doc/expr-poly-conversion.md).
pub fn expr_to_poly(expr: &Expr) -> Result<Poly, EvalError> {
    if expr_contains_alg_coeff(expr) {
        return Err(rational_poly_reject(expr));
    }
    poly_from_expr_shape(expr, rational_leaf, &rational_poly_ring())
}

/// Lift an expression to a polynomial over algebraic coefficients.
///
/// **Stable** — path B in [expr-poly-conversion.md](../../../../.doc/expr-poly-conversion.md).
pub fn poly_alg_from_expr(expr: &Expr) -> Result<PolyAlgExt, EvalError> {
    if !expr_contains_alg_coeff(expr) {
        return Err(EvalError::TypeError(ERR_POLY_ALG_NO_ALG_COEFF));
    }
    poly_from_expr_shape(expr, algext_leaf, &algext_poly_ring())
}

/// Embed a `Poly` over ℚ into `PolyAlgExt` (coefficients in ℚ ⊂ K).
///
/// **Stable** — lift rational polynomial for `poly_algext_roots`.
pub fn poly_algext_from_poly(p: &Poly) -> Result<PolyAlgExt, EvalError> {
    use super::alg_ext::AlgExtData;
    use super::alg_ext_c::AlgExtCData;
    use super::ext_tower::ExtensionField;
    use super::field_arith::coords_to_expr;

    let field = ExtensionField::rational();
    let mut out = PolyAlgExt::ring_zero();
    for (m, r) in &p.terms {
        let coords = field.embed_rational(r);
        let a = AlgExtData::from_field_coords(Arc::clone(&field), coords_to_expr(&coords)?)?;
        let c = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&a)?);
        let mut term = PolyAlgExt::ring_constant(c);
        if !m.is_const() {
            for (v, e) in m.iter() {
                term = term.try_mul(&PolyAlgExt::ring_var(v.clone()).try_pow(e)?)?;
            }
        }
        out = out.try_add(&term)?;
    }
    Ok(out)
}

/// Assemble sparse terms into a sum expression (M1).
// **Pipeline private** — `assemble_poly_expr`
fn assemble_poly_expr<I>(terms: I, zero: ExprArc) -> ExprArc
where
    I: IntoIterator<Item = (Monomial, ExprArc)>,
{
    let mut collected: Vec<(u64, ExprArc)> = terms
        .into_iter()
        .map(|(m, coeff)| (m.degree(), monomial_to_expr(&m, coeff)))
        .collect();
    if collected.is_empty() {
        return zero;
    }
    collected.sort_by(|a, b| b.0.cmp(&a.0));
    Expr::add(collected.into_iter().map(|(_, t)| t).collect())
}

/// **Stable** — `Poly` over ℚ → `Expr` (coefficients remain rational).
pub fn poly_to_expr(poly: &Poly) -> ExprArc {
    assemble_poly_expr(
        poly.terms
            .iter()
            .map(|(m, c)| (m.clone(), ratio_to_expr(c))),
        Expr::int(0),
    )
}

/// **Stable** — `Poly<AlgExtC>` → `Expr` (coefficients as `AlgExt` / `rootof` / `AlgExtC`).
pub fn algext_poly_to_expr(poly: &PolyAlgExt) -> Result<ExprArc, EvalError> {
    Ok(assemble_poly_expr(
        poly.terms.iter().map(|(m, c)| {
            (
                m.clone(),
                c.as_inner().to_expr().into(),
            )
        }),
        Expr::int(0),
    ))
}

/// **Stable** — univariate `Poly` w.r.t. `var` → `poly1[coeffs…]` Expr (giac high-degree-first order).
pub fn univariate_poly_to_poly1_expr(poly: &Poly, var: &Var) -> ExprArc {
    let deg = giac_poly::univariate_degree(poly, var);
    let mut coeffs = Vec::with_capacity((deg + 1) as usize);
    for e in (0..=deg).rev() {
        coeffs.push(ratio_to_expr(&giac_poly::coeff_at(poly, var, e)));
    }
    Expr::func(FuncKind::Poly1, vec![Arc::new(Expr::Seq(coeffs))])
}

/// Format a polynomial over ℤ/pℤ with per-coefficient `(c % p)` display (giac style).
/// **Stable** — PolyMod → Expr
pub fn poly_mod_to_expr(pm: &PolyMod) -> ExprArc {
    let modulus = pm
        .modulus
        .to_string()
        .parse::<i64>()
        .unwrap_or(0);
    if pm.is_zero() {
        return Arc::new(Expr::Mod(Expr::int(0), Expr::int(modulus)));
    }
    assemble_poly_expr(
        pm.terms.iter().map(|(m, c)| {
            let rem = giac_poly::smod(
                c.val.to_string().parse().unwrap_or(0),
                modulus,
            );
            (
                m.clone(),
                Arc::new(Expr::Mod(Expr::int(rem), Expr::int(modulus))),
            )
        }),
        Arc::new(Expr::Mod(Expr::int(0), Expr::int(modulus))),
    )
}

// **Stable** — `u64_to_expr_int`
pub(crate) fn u64_to_expr_int(n: u64) -> ExprArc {
    match i64::try_from(n) {
        Ok(v) => Expr::int(v),
        Err(_) => Arc::new(Expr::Int(BigInt::from(n))),
    }
}

/// **Stable** — `Ratio<BigInt>` → `Expr` (`Int` / `Rat`).
pub fn ratio_to_expr(r: &Ratio<BigInt>) -> ExprArc {
    if r.is_zero() {
        Expr::int(0)
    } else if r.denom() == &BigInt::one() {
        if let Ok(v) = crate::num_util::bigint_to_i64(r.numer()) {
            Expr::int(v)
        } else {
            Arc::new(Expr::Rat(r.clone()))
        }
    } else {
        Arc::new(Expr::Rat(r.clone()))
    }
}

// **Pipeline private** — `monomial_to_expr`
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

/// **Stable** — sorted variables in Expr
pub fn vars_from_expr(expr: &Expr) -> Vec<Var> {
    let mut out = Vec::new();
    collect_vars(expr, &mut out);
    out.sort();
    out.dedup();
    out
}

// **Pipeline private** — `collect_vars`
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
    use super::super::poly_conv::{ERR_ALG_EXT_COEFF, ERR_POLY_ALG_NO_ALG_COEFF, ERR_ROOTOF_COEFF};
    use crate::algebra::test_fixtures::sqrt2_algext_expr;
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
    fn poly_algext_from_poly_embeds_rational() {
        let p = expr_to_poly(
            &Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::int(1),
            ]),
        )
        .unwrap();
        let p_alg = poly_algext_from_poly(&p).unwrap();
        assert_eq!(p_alg.degree(), p.degree());
        assert!(!p_alg.is_zero());
    }

    #[test]
    fn poly_alg_from_expr_constant_rootof() {
        let e = sqrt2_algext_expr();
        let p = poly_alg_from_expr(e.as_ref()).unwrap();
        assert_eq!(p.degree(), 0);
        let back = algext_poly_to_expr(&p).unwrap();
        assert!(matches!(back.as_ref(), Expr::AlgExt(_)));
    }

    #[test]
    fn poly_alg_from_expr_x_squared_minus_two_over_k() {
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(-1), sqrt2_algext_expr()]),
        ]);
        let p = poly_alg_from_expr(&e).unwrap();
        assert_eq!(p.degree(), 2);
        let back = algext_poly_to_expr(&p).unwrap();
        let round = poly_alg_from_expr(back.as_ref()).unwrap();
        assert_eq!(p, round);
    }

    #[test]
    fn flat_uni_algext_degree() {
        use giac_poly::{FlatUni, MainVar};

        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            sqrt2_algext_expr(),
        ]);
        let p = poly_alg_from_expr(&e).unwrap();
        let flat = FlatUni::try_new(p, MainVar::new("x")).unwrap();
        assert_eq!(flat.degree(), 2);
    }

    #[test]
    fn poly_alg_from_expr_univariate_roundtrip() {
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::int(-2),
            sqrt2_algext_expr(),
        ]);
        let p = poly_alg_from_expr(&e).unwrap();
        let back = algext_poly_to_expr(&p).unwrap();
        let p2 = poly_alg_from_expr(back.as_ref()).unwrap();
        assert_eq!(p, p2);
    }

    #[test]
    fn univariate_poly_to_poly1_expr_quadratic() {
        let x = Poly::var("x");
        let p = x.pow(2).sub(&Poly::constant(Ratio::from_integer(2.into())));
        let e = univariate_poly_to_poly1_expr(&p, &Var::from("x"));
        assert!(matches!(e.as_ref(), Expr::Func(FuncKind::Poly1, _)));
    }
}
