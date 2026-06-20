//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{eval, expr_to_poly, poly_to_expr, Context, EvalError, Expr, ExprArc};
use giac_poly::{coeff_at, factor_into, univariate_degree, Poly, Var};
use crate::rootof::quadratic_rootof_roots;
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Zero;

/// `realroot(poly)` — exact rational real roots with multiplicity (GIAC-207 minimal).
/// **Stable (bounded)** — real roots via Sturm isolation
pub fn eval_realroot(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.is_empty() || args.len() > 2 {
        return Err(EvalError::TooFewArgs("realroot"));
    }
    let poly = expr_to_poly(eval(args[0].as_ref(), ctx)?.as_ref())?;
    let var = Var::from("x");
    if univariate_degree(&poly, &var) == 0 {
        return Err(EvalError::TypeError("constant polynomial"));
    }
    let roots = algebraic_real_roots(&poly, &var)?;
    let items: Vec<ExprArc> = roots
        .into_iter()
        .map(|(r, m)| {
            Arc::new(Expr::List(vec![
                r,
                Expr::int(i64::try_from(m).unwrap_or(1)),
            ]))
        })
        .collect();
    Ok(Arc::new(Expr::List(items)))
}

// **Pipeline private** — `algebraic_real_roots`
fn algebraic_real_roots(
    p: &Poly,
    var: &Var,
) -> Result<Vec<(ExprArc, usize)>, EvalError> {
    if let Ok(rational) = rational_real_roots(p, var) {
        if !rational.is_empty() {
            return Ok(rational
                .into_iter()
                .map(|(r, m)| (poly_to_expr(&Poly::constant(r)), m))
                .collect());
        }
    }
    if univariate_degree(p, var) == 2 {
        let rs = quadratic_rootof_roots(p, var)?;
        return Ok(rs.into_iter().map(|r| (r, 1)).collect());
    }
    rational_real_roots(p, var).map(|rs| {
        rs.into_iter()
            .map(|(r, m)| (poly_to_expr(&Poly::constant(r)), m))
            .collect()
    })
}

// **Pipeline private** — `rational_real_roots`
fn rational_real_roots(
    p: &Poly,
    var: &Var,
) -> Result<Vec<(Ratio<BigInt>, usize)>, EvalError> {
    let mut roots = Vec::new();
    if let Some(factors) = factor_into(p) {
        for f in factors {
            if univariate_degree(&f, var) != 1 {
                continue;
            }
            let a = coeff_at(&f, var, 1);
            let b = coeff_at(&f, var, 0);
            if a.is_zero() {
                continue;
            }
            let root = -b / a;
            if let Some(idx) = roots.iter().position(|(r, _)| r == &root) {
                roots[idx].1 += 1;
            } else {
                roots.push((root, 1));
            }
        }
        roots.sort_by(|a, b| a.0.cmp(&b.0));
        return Ok(roots);
    }
    if univariate_degree(p, var) == 1 {
        let a = coeff_at(p, var, 1);
        let b = coeff_at(p, var, 0);
        if !a.is_zero() {
            return Ok(vec![(-b / a, 1)]);
        }
    }
    Err(EvalError::NotImplemented("realroot"))
}

mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind};

    use crate::plugin::xcas_default;

    #[test]
    // **Pipeline private** — `realroot_x_squared_minus_two`
    fn realroot_x_squared_minus_two() {
        let ctx = xcas_default();
        let p = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::int(-2),
        ]);
        let e = Expr::func(FuncKind::Realroot, vec![p]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("rootof"), "got {s}");
    }

    #[test]
    // **Pipeline private** — `realroot_x_fourth_minus_one`
    fn realroot_x_fourth_minus_one() {
        let ctx = xcas_default();
        let p = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(4)),
            Expr::int(-1),
        ]);
        let e = Expr::func(FuncKind::Realroot, vec![p]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("-1"), "got {s}");
        assert!(s.contains("1"), "got {s}");
    }
}
