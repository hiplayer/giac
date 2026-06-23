//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{
    canonicalize_to_algext_c, eval, expr_to_poly, ident_from_expr, poly_to_expr, Context,
    EvalError, Expr, ExprArc,
};
use giac_poly::{
    coeff_at, factor_into, square_free_factorization, sturmab_count, univariate_degree, Poly, Var,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Zero;

use crate::solve_poly::solve_irreducible_factor;

/// `realroot(poly)` or `realroot(poly,x)` — exact real roots with multiplicity (S6).
/// **Stable (bounded)** — rational roots + filtered exact algebraic (deg≤4 via S0)
pub fn eval_realroot(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.is_empty() || args.len() > 2 {
        return Err(EvalError::TooFewArgs("realroot"));
    }
    let poly = expr_to_poly(eval(args[0].as_ref(), ctx)?.as_ref())?;
    let var = match args.len() {
        1 => Var::from("x"),
        _ => {
            let id = ident_from_expr(&args[1])?;
            Var::from(id.as_str())
        }
    };
    if univariate_degree(&poly, &var) == 0 {
        return Err(EvalError::TypeError("constant polynomial"));
    }
    let roots = real_roots_with_multiplicity(&poly, &var, ctx)?;
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

// **Pipeline private** — merge rational + exact real algebraic roots
fn real_roots_with_multiplicity(
    p: &Poly,
    var: &Var,
    ctx: &Context,
) -> Result<Vec<(ExprArc, usize)>, EvalError> {
    let mut out = rational_real_roots(p, var)?
        .into_iter()
        .map(|(r, m)| (poly_to_expr(&Poly::constant(r)), m))
        .collect::<Vec<_>>();
    for (factor, mult) in square_free_factorization(p, var)? {
        if factor.is_zero() || univariate_degree(&factor, var) == 0 {
            continue;
        }
        if univariate_degree(&factor, var) == 1 {
            continue;
        }
        for root in real_algebraic_roots(&factor, var, ctx)? {
            merge_root(&mut out, root, mult);
        }
    }
    if out.is_empty() {
        return Err(EvalError::NotImplemented("realroot"));
    }
    out.sort_by(|a, b| format_root_key(&a.0).cmp(&format_root_key(&b.0)));
    Ok(out)
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
            merge_rational(&mut roots, root, 1);
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
    Ok(Vec::new())
}

// **Pipeline private** — S6: exact roots filtered by Sturm count / casus cubic
fn real_algebraic_roots(
    factor: &Poly,
    var: &Var,
    ctx: &Context,
) -> Result<Vec<ExprArc>, EvalError> {
    let rs = solve_irreducible_factor(factor, var, ctx)?;
    let bound = Ratio::from_integer(BigInt::from(1_000_000));
    let neg_bound = -bound.clone();
    let real_count = sturmab_count(factor, var, &neg_bound, &bound)?;
    if real_count == 0 {
        return Ok(vec![]);
    }
    if univariate_degree(factor, var) == 3 && real_count == 1 && rs.len() == 3 {
        if monic_cubic_discriminant(factor, var)
            .map(|d| d < Ratio::zero())
            .unwrap_or(false)
        {
            // ponytail: casus cubic — real root is adjoin-first in poly_roots sort
            return Ok(vec![rs[0].clone()]);
        }
    }
    let filtered: Vec<_> = rs.into_iter().filter(|r| expr_is_real(r)).collect();
    if filtered.len() >= real_count {
        return Ok(filtered.into_iter().take(real_count).collect());
    }
    Ok(filtered)
}

/// Discriminant of monic cubic `x³+bx²+cx+d` over ℚ.
fn monic_cubic_discriminant(factor: &Poly, var: &Var) -> Option<Ratio<BigInt>> {
    if univariate_degree(factor, var) != 3 {
        return None;
    }
    let b = coeff_at(factor, var, 2);
    let c = coeff_at(factor, var, 1);
    let d = coeff_at(factor, var, 0);
    let b2c2 = &b * &b * &c * &c;
    let four_c3 = Ratio::from_integer(4.into()) * &c * &c * &c;
    let four_b3d = Ratio::from_integer(4.into()) * &b * &b * &b * &d;
    let twenty_seven_d2 = Ratio::from_integer(27.into()) * &d * &d;
    let eighteen_bcd = Ratio::from_integer(18.into()) * &b * &c * &d;
    Some(b2c2 - four_c3 - four_b3d - twenty_seven_d2 + eighteen_bcd)
}

fn merge_rational(out: &mut Vec<(Ratio<BigInt>, usize)>, root: Ratio<BigInt>, mult: usize) {
    if let Some(idx) = out.iter().position(|(r, _)| r == &root) {
        out[idx].1 += mult;
    } else {
        out.push((root, mult));
    }
}

fn merge_root(out: &mut Vec<(ExprArc, usize)>, root: ExprArc, mult: usize) {
    if let Some(idx) = out.iter().position(|(r, _)| roots_eq(r, &root)) {
        out[idx].1 += mult;
    } else {
        out.push((root, mult));
    }
}

fn roots_eq(a: &ExprArc, b: &ExprArc) -> bool {
    Arc::ptr_eq(a, b) || format_root_key(a) == format_root_key(b)
}

fn format_root_key(e: &ExprArc) -> String {
    format!("{e:?}")
}

fn expr_is_real(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Int(_) | Expr::Rat(_) => true,
        _ => canonicalize_to_algext_c(e.as_ref())
            .map(|z| z.im.iter().all(|c| expr_is_zero(c)))
            .unwrap_or(false),
    }
}

fn expr_is_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
        || matches!(e.as_ref(), Expr::Rat(r) if r.is_zero())
}

#[cfg(test)]
mod tests {
    //! Test tiers — `.doc/test-writing-spec.md` · audit §2

    use giac_core::{eval, Expr, FuncKind, Ident};

    use super::*;
    use crate::plugin::xcas_default;
    use crate::test_verify::test_verify::{assert_is_algext_or_rootof, assert_realroot_has, list_items};

    // **A** — eval(Realroot); AlgExt/rootof pairs.
    #[test]
    fn realroot_x_squared_minus_two() {
        let ctx = xcas_default();
        let p = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::int(-2),
        ]);
        let e = Expr::func(FuncKind::Realroot, vec![p]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let entries = list_items(&r);
        assert_eq!(entries.len(), 2);
        for pair in entries {
            assert_is_algext_or_rootof(&list_items(pair)[0]);
        }
    }

    // **A** — eval(Realroot); assert_realroot_has.
    #[test]
    fn realroot_x_fourth_minus_one() {
        let ctx = xcas_default();
        let p = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(4)),
            Expr::int(-1),
        ]);
        let e = Expr::func(FuncKind::Realroot, vec![p]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let items = list_items(&r);
        assert_realroot_has(items, &Expr::int(-1), &ctx);
        assert_realroot_has(items, &Expr::int(1), &ctx);
    }

    // **B** — S6: one real root of x³−x−1 (two complex roots filtered out).
    #[test]
    fn realroot_x_cubed_minus_x_minus_1() {
        let ctx = xcas_default();
        let p = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(3)),
            Expr::mul(vec![Expr::int(-1), Expr::sym("x")]),
            Expr::int(-1),
        ]);
        let e = Expr::func(FuncKind::Realroot, vec![p, Expr::sym("x")]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        let items = list_items(&r);
        assert_eq!(items.len(), 1);
        assert_is_algext_or_rootof(&list_items(&items[0])[0]);
    }
}
