use std::sync::Arc;

use giac_core::{poly_to_expr, EvalError, Expr, ExprArc, FuncKind};
use giac_poly::{coeff_at, univariate_degree, Poly, Var};

/// Two `rootof` branches for quadratic irrational roots of `poly` in `var`.
pub fn quadratic_rootof_roots(poly: &Poly, var: &Var) -> Result<Vec<ExprArc>, EvalError> {
    if univariate_degree(poly, var) != 2 {
        return Err(EvalError::TypeError("expected quadratic"));
    }
    let minpoly = poly1_from_univariate(poly, var);
    Ok(vec![
        rootof_expr(&[1, 0], &minpoly),
        rootof_expr(&[-1, 0], &minpoly),
    ])
}

fn poly1_from_univariate(poly: &Poly, var: &Var) -> ExprArc {
    let deg = univariate_degree(poly, var);
    let mut coeffs = Vec::with_capacity((deg + 1) as usize);
    for e in (0..=deg).rev() {
        coeffs.push(poly_to_expr(&Poly::constant(coeff_at(poly, var, e))));
    }
    Expr::func(FuncKind::Poly1, vec![Arc::new(Expr::Seq(coeffs))])
}

fn rootof_expr(num: &[i64], minpoly: &ExprArc) -> ExprArc {
    let num_seq = Arc::new(Expr::Seq(num.iter().map(|&n| Expr::int(n)).collect()));
    giac_core::AlgExtData::from_rootof(&num_seq, minpoly)
        .expect("quadratic rootof")
        .into_expr()
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Context, FuncKind, RelOp};
    use giac_poly::{roots, Poly, Var};

    use super::*;
    use crate::plugin::xcas_default;

    fn x() -> Poly {
        Poly::var("t")
    }

    #[test]
    fn quadratic_rootof_has_two_branches() {
        let p = x().pow(2).sub(&Poly::constant(num_rational::Ratio::from_integer(
            num_bigint::BigInt::from(2),
        )));
        let rs = quadratic_rootof_roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 2);
        assert!(format_expr(rs[0].as_ref()).contains("rootof"));
    }

    #[test]
    fn solve_t_squared_minus_two_uses_rootof() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Solve,
            vec![
                Arc::new(Expr::Relation(
                    RelOp::Eq,
                    Expr::add(vec![
                        Expr::pow(Expr::sym("t"), Expr::int(2)),
                        Expr::int(-2),
                    ]),
                    Expr::int(0),
                )),
                Expr::sym("t"),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("rootof"), "got {s}");
    }

    #[test]
    fn rational_quadratic_still_uses_roots() {
        let p = x().pow(2).sub(&Poly::one());
        let rs = roots(&p, &Var::from("t")).unwrap();
        assert_eq!(rs.len(), 2);
    }
}
