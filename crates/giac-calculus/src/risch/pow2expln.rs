//! Rewrite `base^exp` with non-constant exponent as `exp(exp*ln(base))` (GIAC `pow2expln`).

use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind, Ident};

fn depends_on_var(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id == var,
        Expr::Int(_) | Expr::Rat(_) => false,
        Expr::Add(ts) | Expr::Mul(ts) => ts.iter().any(|t| depends_on_var(t, var)),
        Expr::Pow(b, exp) => depends_on_var(b, var) || depends_on_var(exp, var),
        Expr::Frac(n, d) => depends_on_var(n, var) || depends_on_var(d, var),
        Expr::Func(_, args) => args.iter().any(|a| depends_on_var(a, var)),
        _ => false,
    }
}

fn is_const_wrt(e: &ExprArc, var: &Ident) -> bool {
    !depends_on_var(e, var)
}

/// `pow2expln(e, x)` — subset of GIAC `subst.cc::pow2expln(e, x)`.
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
    use giac_core::format_expr;

    use super::*;

    fn x() -> ExprArc {
        Expr::sym("x")
    }

    #[test]
    fn pow2expln_x_to_x() {
        let var = Ident::new("x");
        let e = Expr::pow(x(), x());
        let r = pow2expln(&e, &var);
        let s = format_expr(r.as_ref());
        assert!(s.contains("exp") && s.contains("ln"), "got {s}");
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
        let var = Ident::new("x");
        let e = Expr::pow(Expr::int(2), x());
        let r = pow2expln(&e, &var);
        let s = format_expr(r.as_ref());
        assert!(s.contains("exp") && s.contains("ln(2)"), "got {s}");
    }

    #[test]
    fn pow2expln_nested_in_mul() {
        let var = Ident::new("x");
        let e = Expr::mul(vec![Expr::pow(Expr::int(2), x()), Expr::func(FuncKind::Ln, vec![x()])]);
        let r = pow2expln(&e, &var);
        assert!(matches!(r.as_ref(), Expr::Mul(_)));
    }
}
