//! Shared expression utilities (variable dependence, etc.).
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../.doc/giac-calculus-api-stability.md) §6.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `depends_on_var`, `is_const_wrt` |

use giac_core::{Expr, ExprArc, Ident};

/// **Stable** — whether `e` syntactically depends on `var`.
pub fn depends_on_var(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id == var,
        Expr::Int(_) | Expr::Rat(_) => false,
        Expr::Add(ts) | Expr::Mul(ts) => ts.iter().any(|t| depends_on_var(t, var)),
        Expr::Pow(b, exp) => depends_on_var(b, var) || depends_on_var(exp, var),
        Expr::Frac(n, d) => depends_on_var(n, var) || depends_on_var(d, var),
        Expr::AlgExt(a) => {
            a.min_poly.iter().any(|c| depends_on_var(c, var))
                || a.coords.iter().any(|c| depends_on_var(c, var))
        }
        Expr::Func(_, args) => args.iter().any(|a| depends_on_var(a, var)),
        _ => false,
    }
}

/// **Stable** — whether `e` is constant with respect to `var` (syntactic).
pub fn is_const_wrt(e: &ExprArc, var: &Ident) -> bool {
    !depends_on_var(e, var)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depends_on_var_symbol() {
        let var = Ident::new("x");
        assert!(depends_on_var(&Expr::sym("x"), &var));
        assert!(!depends_on_var(&Expr::sym("y"), &var));
        assert!(!depends_on_var(&Expr::int(1), &var));
    }

    #[test]
    fn depends_on_var_nested() {
        let var = Ident::new("x");
        let e = Expr::add(vec![Expr::sym("x"), Expr::int(1)]);
        assert!(depends_on_var(&e, &var));
        assert!(is_const_wrt(&Expr::int(2), &var));
    }

    #[test]
    fn depends_on_var_pow_frac_func() {
        let var = Ident::new("x");
        let pow = Expr::pow(Expr::sym("x"), Expr::int(2));
        let frac: ExprArc = Expr::Frac(Expr::sym("x"), Expr::int(1)).into();
        let f = Expr::func(giac_core::FuncKind::Sin, vec![Expr::sym("x")]);
        assert!(depends_on_var(&pow, &var));
        assert!(depends_on_var(&frac, &var));
        assert!(depends_on_var(&f, &var));
        assert!(is_const_wrt(&pow, &Ident::new("y")));
    }
}
