//! Shared expression utilities (variable dependence, etc.).

use giac_core::{Expr, ExprArc, Ident};

/// Whether `e` syntactically depends on `var`.
pub fn depends_on_var(e: &ExprArc, var: &Ident) -> bool {
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

/// Whether `e` is constant with respect to `var` (syntactic).
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
}
