//! GIAC-226: elementary extension tower (`risch_tower` / `rlvarx` subset).

use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind, Ident};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RischTowerError {
    NotElementary,
}

/// Logarithmic / exponential extension variables in `expr` that depend on `var`.
pub fn rlvarx(expr: &ExprArc, var: &Ident) -> Vec<ExprArc> {
    let mut out = Vec::new();
    collect_rlvarx(expr, var, &mut out);
    out.sort_by(|a, b| extension_rank(a).cmp(&extension_rank(b)));
    out.dedup_by(|a, b| a == b);
    out
}

/// Returns the tower (most complex extension first) when `expr` is elementary over `var`.
pub fn risch_tower(expr: &ExprArc, var: &Ident) -> Result<Vec<ExprArc>, RischTowerError> {
    if contains_non_elementary_transcendental(expr, var) {
        return Err(RischTowerError::NotElementary);
    }
    let mut tower = rlvarx(expr, var);
    tower.retain(|e| !is_var(e, var));
    for ext in &tower {
        if !is_exp_or_ln(ext) {
            return Err(RischTowerError::NotElementary);
        }
    }
    tower.reverse();
    Ok(tower)
}

fn collect_rlvarx(expr: &ExprArc, var: &Ident, out: &mut Vec<ExprArc>) {
    if !depends_on(expr, var) {
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
            if depends_on(exp, var) {
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

fn is_exp_or_ln(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Exp | FuncKind::Ln, _)
    )
}

fn contains_non_elementary_transcendental(e: &ExprArc, var: &Ident) -> bool {
    if !depends_on(e, var) {
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

fn depends_on(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => id == var,
        Expr::Int(_) | Expr::Rat(_) => false,
        Expr::Add(ts) | Expr::Mul(ts) => ts.iter().any(|t| depends_on(t, var)),
        Expr::Pow(b, exp) => depends_on(b, var) || depends_on(exp, var),
        Expr::Frac(n, d) => depends_on(n, var) || depends_on(d, var),
        Expr::Func(_, args) => args.iter().any(|a| depends_on(a, var)),
        _ => false,
    }
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn push_unique(out: &mut Vec<ExprArc>, e: ExprArc) {
    if !out.iter().any(|x| x == &e) {
        out.push(e);
    }
}

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
}
