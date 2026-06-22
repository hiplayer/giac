//! GIAC-226: elementary extension tower (`risch_tower` / `rlvarx` subset).
//!
//! See [`.doc/giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md) §5.
//!
//! | Tier | 函数 |
//! |------|------|
//! | **Stable** | `rlvarx`, `risch_tower` |
//! | **Pipeline private** | `collect_rlvarx`, `is_exp_or_ln`, `contains_non_elementary_transcendental`, `extension_rank` |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind, Ident};

use crate::expr_util::{depends_on_var, is_var};
use super::pow2expln::pow2expln;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RischTowerError {
    NotElementary,
}

/// **Stable** — logarithmic / exponential extension variables in `expr` that depend on `var`.
pub fn rlvarx(expr: &ExprArc, var: &Ident) -> Vec<ExprArc> {
    let mut out = Vec::new();
    collect_rlvarx(expr, var, &mut out);
    out.sort_by(|a, b| extension_rank(a).cmp(&extension_rank(b)));
    out.dedup_by(|a, b| a == b);
    out
}

/// **Stable** — returns the tower (most complex extension first) when `expr` is elementary over `var`.
pub fn risch_tower(expr: &ExprArc, var: &Ident) -> Result<Vec<ExprArc>, RischTowerError> {
    let normalized = pow2expln(expr, var);
    if contains_non_elementary_transcendental(&normalized, var) {
        return Err(RischTowerError::NotElementary);
    }
    let mut tower = rlvarx(&normalized, var);
    tower.retain(|e| !is_var(e, var));
    for ext in &tower {
        if !is_exp_or_ln(ext) {
            return Err(RischTowerError::NotElementary);
        }
    }
    tower.reverse();
    Ok(tower)
}

// **Pipeline private** — collect `exp`/`ln` extension atoms depending on `var`.
fn collect_rlvarx(expr: &ExprArc, var: &Ident, out: &mut Vec<ExprArc>) {
    if !depends_on_var(expr, var) {
        return;
    }
    match expr.as_ref() {
        Expr::Func(FuncKind::Exp | FuncKind::Ln, _) => {
            push_unique(out, Arc::clone(expr));
            if let Expr::Func(_, args) = expr.as_ref() {
                for a in args {
                    collect_rlvarx(a, var, out);
                }
            }
        }
        Expr::Pow(base, exp) => {
            if depends_on_var(exp, var) {
                let ln_base = Expr::func(FuncKind::Ln, vec![Arc::clone(base)]);
                push_unique(out, ln_base);
            }
            collect_rlvarx(base, var, out);
            collect_rlvarx(exp, var, out);
        }
        Expr::Add(ts) | Expr::Mul(ts) => {
            for t in ts {
                collect_rlvarx(t, var, out);
            }
        }
        Expr::Frac(n, d) => {
            collect_rlvarx(n, var, out);
            collect_rlvarx(d, var, out);
        }
        Expr::Func(_, args) => {
            for a in args {
                collect_rlvarx(a, var, out);
            }
        }
        Expr::Symbol(_) | Expr::Int(_) | Expr::Rat(_) => {}
        _ => {}
    }
}

// **Pipeline private** — `exp` or `ln` top-level form.
fn is_exp_or_ln(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp | FuncKind::Ln, _)
    )
}

// **Pipeline private** — detect non-elementary transcendentals (e.g. trig).
fn contains_non_elementary_transcendental(e: &ExprArc, var: &Ident) -> bool {
    if !depends_on_var(e, var) {
        return false;
    }
    match e.as_ref() {
        Expr::Func(FuncKind::Exp | FuncKind::Ln, args) => {
            args.iter().any(|a| contains_non_elementary_transcendental(a, var))
        }
        Expr::Func(_, _) => true,
        Expr::Pow(b, exp) => {
            contains_non_elementary_transcendental(b, var)
                || contains_non_elementary_transcendental(exp, var)
        }
        Expr::Add(ts) | Expr::Mul(ts) => {
            ts.iter()
                .any(|t| contains_non_elementary_transcendental(t, var))
        }
        Expr::Frac(n, d) => {
            contains_non_elementary_transcendental(n, var)
                || contains_non_elementary_transcendental(d, var)
        }
        Expr::Symbol(_) | Expr::Int(_) | Expr::Rat(_) => false,
        _ => false,
    }
}

// **Pipeline private** — append extension atom if not already present.
fn push_unique(out: &mut Vec<ExprArc>, e: ExprArc) {
    if !out.iter().any(|x| x == &e) {
        out.push(e);
    }
}

// **Pipeline private** — nesting depth for tower ordering.
fn extension_rank(e: &ExprArc) -> usize {
    match e.as_ref() {
        Expr::Symbol(_) => 0,
        Expr::Func(FuncKind::Exp | FuncKind::Ln, args) => {
            1 + args.iter().map(extension_rank).max().unwrap_or(0)
        }
        Expr::Pow(b, exp) => extension_rank(b) + extension_rank(exp),
        Expr::Add(ts) | Expr::Mul(ts) => ts.iter().map(extension_rank).max().unwrap_or(0),
        Expr::Frac(n, d) => extension_rank(n).max(extension_rank(d)),
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn x() -> ExprArc {
        Expr::sym("x")
    }

    #[test]
    fn rlvarx_exp_and_ln() {
        let var = Ident::new("x");
        let e = Expr::mul(vec![
            Expr::func(FuncKind::Exp, vec![x()]),
            Expr::func(FuncKind::Ln, vec![x()]),
        ]);
        let v = rlvarx(&e, &var);
        assert_eq!(v.len(), 2);
        assert!(risch_tower(&e, &var).is_ok());
    }

    #[test]
    fn risch_tower_rejects_trig() {
        let var = Ident::new("x");
        let e = Expr::func(FuncKind::Sin, vec![x()]);
        assert_eq!(risch_tower(&e, &var), Err(RischTowerError::NotElementary));
    }

    #[test]
    fn risch_tower_polynomial_is_empty() {
        let var = Ident::new("x");
        let e = Expr::pow(x(), Expr::int(2));
        let tower = risch_tower(&e, &var).unwrap();
        assert!(tower.is_empty());
    }

    #[test]
    fn rlvarx_nested_exp() {
        let var = Ident::new("x");
        let e = Expr::func(FuncKind::Exp, vec![Expr::func(FuncKind::Exp, vec![x()])]);
        let v = rlvarx(&e, &var);
        assert!(!v.is_empty());
        assert!(risch_tower(&e, &var).is_ok());
    }

    #[test]
    fn risch_tower_pow2expln_exp_x() {
        let var = Ident::new("x");
        let e = Expr::pow(Expr::int(2), x());
        let tower = risch_tower(&e, &var).unwrap();
        assert_eq!(tower.len(), 1);
        assert!(matches!(tower[0].as_ref(), Expr::Func(FuncKind::Exp, _)));
    }

    #[test]
    fn risch_tower_rejects_nested_trig() {
        let var = Ident::new("x");
        let e = Expr::func(FuncKind::Sin, vec![Expr::func(FuncKind::Cos, vec![x()])]);
        assert_eq!(risch_tower(&e, &var), Err(RischTowerError::NotElementary));
    }

    #[test]
    fn rlvarx_independent_of_var() {
        let var = Ident::new("x");
        let e = Expr::func(FuncKind::Exp, vec![Expr::sym("y")]);
        assert!(rlvarx(&e, &var).is_empty());
    }
}
