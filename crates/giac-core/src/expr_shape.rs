//! Syntactic expression shape predicates shared across giac-rs crates.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Pipeline private** | `is_var`, `var_to_expr`, `is_sin_of_var`, `is_cos_of_var`, `is_ln_of_var` |

use crate::{Expr, ExprArc, FuncKind, Ident};

/// Whether `e` is exactly the symbol `var`.
pub fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

/// Whether `e` is exactly the symbol `var` (by-value tree reference).
pub fn is_var_expr(e: &Expr, var: &Ident) -> bool {
    matches!(e, Expr::Symbol(id) if id == var)
}

/// Integration/limit variable as a symbol expression.
pub fn var_to_expr(var: &Ident) -> ExprArc {
    Expr::sym(var.as_str())
}

/// Whether `e` is `sin(var)`.
pub fn is_sin_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_var(&args[0], var)
    )
}

/// Whether `e` is `sin(var)` (by-value tree reference).
pub fn is_sin_of_var_expr(e: &Expr, var: &Ident) -> bool {
    matches!(
        e,
        Expr::Func(FuncKind::Sin, args)
            if args.len() == 1 && is_var_expr(args[0].as_ref(), var)
    )
}

/// Whether `e` is `cos(var)`.
pub fn is_cos_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Cos, args) if args.len() == 1 && is_var(&args[0], var)
    )
}

/// Whether `e` is `ln(var)`.
pub fn is_ln_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && is_var(&args[0], var)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_var_and_trig_shapes() {
        let var = Ident::new("x");
        assert!(is_var_expr(Expr::sym("x").as_ref(), &var));
        assert!(!is_var_expr(Expr::sym("y").as_ref(), &var));
        let sin_x = Expr::func(FuncKind::Sin, vec![Expr::sym("x")]);
        assert!(is_sin_of_var(&sin_x, &var));
        assert!(is_sin_of_var_expr(sin_x.as_ref(), &var));
        assert!(is_cos_of_var(
            &Expr::func(FuncKind::Cos, vec![Expr::sym("x")]),
            &var
        ));
        assert!(is_ln_of_var(
            &Expr::func(FuncKind::Ln, vec![Expr::sym("x")]),
            &var
        ));
    }

    #[test]
    fn var_to_expr_builds_symbol() {
        let var = Ident::new("t");
        assert!(matches!(var_to_expr(&var).as_ref(), Expr::Symbol(id) if id == &var));
    }
}
