//! Rewrite `base^exp` with non-constant exponent as `exp(exp*ln(base))` (GIAC `pow2expln`).
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md) §5.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `pow2expln` |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind, Ident};

use crate::expr_util::{depends_on_var, is_const_wrt};

/// **Stable** — `pow2expln(e, x)`; subset of GIAC `subst.cc::pow2expln(e, x)`.
pub fn pow2expln(expr: &ExprArc, var: &Ident) -> ExprArc {
    match expr.as_ref() {
        Expr::Pow(base, exp) => {
            let base_p = pow2expln(base, var);
            let exp_p = pow2expln(exp, var);
            if depends_on_var(exp, var) || (!is_const_wrt(exp, var) && depends_on_var(base, var)) {
                return Expr::func(
                    FuncKind::Exp,
                    vec![Expr::mul(vec![
                        exp_p,
                        Expr::func(FuncKind::Ln, vec![base_p]),
                    ])],
                );
            }
            Expr::pow(base_p, exp_p)
        }
        Expr::Add(ts) => Expr::add(ts.iter().map(|t| pow2expln(t, var)).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(|f| pow2expln(f, var)).collect()),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            pow2expln(n, var),
            pow2expln(d, var),
        )),
        Expr::Func(k, args) => {
            Expr::func(*k, args.iter().map(|a| pow2expln(a, var)).collect())
        }
        _ => Arc::clone(expr),
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{format_expr, Expr, FuncKind, Ident};

    use crate::plugin::xcas_default;

    use super::*;

    fn x() -> ExprArc {
        Expr::sym("x")
    }

    #[test]
    fn pow2expln_x_to_x() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let e = Expr::pow(x(), x());
        let r = pow2expln(&e, &var);
        let expected = Expr::func(
            FuncKind::Exp,
            vec![Expr::mul(vec![
                Arc::clone(&x()),
                Expr::func(FuncKind::Ln, vec![x()]),
            ])],
        );
        assert!(
            giac_simplify::assert_equiv(r.as_ref(), &expected, &ctx).expect("assert_equiv"),
            "got {}",
            format_expr(r.as_ref())
        );
    }

    #[test]
    fn pow2expln_integer_power_unchanged() {
        let var = Ident::new("x");
        let e = Expr::pow(x(), Expr::int(3));
        let r = pow2expln(&e, &var);
        assert_eq!(format_expr(r.as_ref()), "x^3");
    }

    #[test]
    fn pow2expln_x_to_const_exp() {
        let ctx = xcas_default();
        let var = Ident::new("x");
        let e = Expr::pow(Expr::int(2), x());
        let r = pow2expln(&e, &var);
        let expected = Expr::func(
            FuncKind::Exp,
            vec![Expr::mul(vec![
                x(),
                Expr::func(FuncKind::Ln, vec![Expr::int(2)]),
            ])],
        );
        assert!(
            giac_simplify::assert_equiv(r.as_ref(), &expected, &ctx).expect("assert_equiv"),
            "got {}",
            format_expr(r.as_ref())
        );
    }

    #[test]
    fn pow2expln_nested_in_mul() {
        let var = Ident::new("x");
        let e = Expr::mul(vec![Expr::pow(Expr::int(2), x()), Expr::func(FuncKind::Ln, vec![x()])]);
        let r = pow2expln(&e, &var);
        assert!(matches!(r.as_ref(), Expr::Mul(_)));
    }
}
