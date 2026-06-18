//! Expression size / depth guards — avoid expand/series blowups on heavy exp forms.

use giac_core::{Expr, ExprArc, FuncKind, Ident};

pub(crate) const MAX_SERIES_ORDER: usize = 10;
/// Upper bound for MRV `mrv_lead_term` ordre escalation (`series.cc` `max_series_expansion_order`).
pub(crate) const MAX_SERIES_EXPANSION_ORDER: usize = 24;
pub(crate) const MAX_SERIES_TERMS: usize = 24;
pub(crate) const MAX_SERIES_DEPTH: usize = 32;
pub(crate) const MAX_EXPAND_NODES: usize = 256;

pub(crate) fn expr_nodes(e: &ExprArc) -> usize {
    match e.as_ref() {
        Expr::Int(_) | Expr::Rat(_) | Expr::Symbol(_) => 1,
        Expr::Add(ts) | Expr::Mul(ts) => 1 + ts.iter().map(expr_nodes).sum::<usize>(),
        Expr::Pow(b, exp) => 1 + expr_nodes(b) + expr_nodes(exp),
        Expr::Frac(n, d) => 1 + expr_nodes(n) + expr_nodes(d),
        Expr::Func(_, args) => 1 + args.iter().map(expr_nodes).sum::<usize>(),
        _ => 1,
    }
}

pub(crate) fn expr_depth(e: &ExprArc) -> usize {
    match e.as_ref() {
        Expr::Int(_) | Expr::Rat(_) | Expr::Symbol(_) => 0,
        Expr::Add(ts) | Expr::Mul(ts) => {
            1 + ts.iter().map(expr_depth).max().unwrap_or(0)
        }
        Expr::Pow(b, exp) => 1 + expr_depth(b).max(expr_depth(exp)),
        Expr::Frac(n, d) => 1 + expr_depth(n).max(expr_depth(d)),
        Expr::Func(_, args) => 1 + args.iter().map(expr_depth).max().unwrap_or(0),
        _ => 0,
    }
}

pub(crate) fn too_heavy_for_expand(e: &ExprArc) -> bool {
    expr_nodes(e) > MAX_EXPAND_NODES || expr_contains_nested_exp(e)
}

pub(crate) fn mrv_series_eligible(e: &ExprArc) -> bool {
    expr_contains_exp(e)
        && expr_nodes(e) <= MAX_EXPAND_NODES
        && expr_depth(e) <= MAX_SERIES_DEPTH
}

/// Limit at `+infinity` via `mrv_lead_term` (upstream `unidirectional_limit`); broader than series-only gate.
pub(crate) fn mrv_limit_eligible(e: &ExprArc) -> bool {
    expr_nodes(e) <= MAX_EXPAND_NODES && expr_depth(e) <= MAX_SERIES_DEPTH
}

/// MRV rewrite is bounded when the rewritten tree stays under the node cap.
pub(crate) fn mrv_rewrite_bounded(e: &ExprArc) -> bool {
    expr_nodes(e) <= MAX_EXPAND_NODES
}

pub(crate) fn expr_contains_exp(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Func(FuncKind::Exp, _) => true,
        Expr::Add(ts) => ts.iter().any(expr_contains_exp),
        Expr::Mul(fs) => fs.iter().any(expr_contains_exp),
        Expr::Pow(b, exp) => expr_contains_exp(b) || expr_contains_exp(exp),
        Expr::Frac(n, d) => expr_contains_exp(n) || expr_contains_exp(d),
        Expr::Func(_, args) => args.iter().any(expr_contains_exp),
        _ => false,
    }
}

pub(crate) fn expr_contains_nested_exp(e: &ExprArc) -> bool {
    fn walk(e: &ExprArc, in_exp: bool) -> bool {
        match e.as_ref() {
            Expr::Func(FuncKind::Exp, args) => {
                if in_exp {
                    return true;
                }
                args.iter().any(|a| walk(a, true))
            }
            Expr::Add(ts) => ts.iter().any(|t| walk(t, in_exp)),
            Expr::Mul(fs) => fs.iter().any(|t| walk(t, in_exp)),
            Expr::Pow(b, exp) => walk(b, in_exp) || walk(exp, in_exp),
            Expr::Frac(n, d) => walk(n, in_exp) || walk(d, in_exp),
            Expr::Func(_, args) => args.iter().any(|a| walk(a, in_exp)),
            _ => false,
        }
    }
    walk(e, false)
}

pub(crate) fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

#[cfg(test)]
mod tests {
    use giac_core::{Expr, ExprArc, FuncKind};

    use super::*;

    #[test]
    fn expr_nodes_and_depth() {
        let e = Expr::add(vec![Expr::sym("x"), Expr::int(1)]);
        assert!(expr_nodes(&e) >= 3);
        assert!(expr_depth(&e) >= 1);
    }

    #[test]
    fn mrv_eligibility_flags() {
        let light = Expr::func(FuncKind::Exp, vec![Expr::sym("x")]);
        assert!(mrv_series_eligible(&light));
        assert!(mrv_limit_eligible(&light));
        assert!(mrv_rewrite_bounded(&light));
        assert!(!too_heavy_for_expand(&light));
    }

    #[test]
    fn nested_exp_detected() {
        let nested = Expr::func(
            FuncKind::Exp,
            vec![Expr::func(FuncKind::Exp, vec![Expr::sym("x")])],
        );
        assert!(expr_contains_nested_exp(&nested));
        assert!(too_heavy_for_expand(&nested) || expr_contains_nested_exp(&nested));
    }

    #[test]
    fn is_var_helper() {
        let var = Ident::new("x");
        assert!(is_var(&Expr::sym("x"), &var));
        assert!(!is_var(&Expr::sym("y"), &var));
    }

    #[test]
    fn expr_contains_exp_paths() {
        let exp_x = Expr::func(FuncKind::Exp, vec![Expr::sym("x")]);
        let add = Expr::add(vec![exp_x.clone(), Expr::int(1)]);
        let mul = Expr::mul(vec![exp_x.clone(), Expr::sym("y")]);
        let frac: ExprArc = Expr::Frac(Expr::int(1), exp_x.clone()).into();
        let pow = Expr::pow(exp_x.clone(), Expr::int(2));
        assert!(expr_contains_exp(&add));
        assert!(expr_contains_exp(&mul));
        assert!(expr_contains_exp(&frac));
        assert!(expr_contains_exp(&pow));
        assert!(!expr_contains_exp(&Expr::sym("x")));
    }

    #[test]
    fn expr_nodes_relation_branch() {
        use std::sync::Arc;

        let e = Arc::new(Expr::Relation(
            giac_core::RelOp::Eq,
            Expr::sym("x"),
            Expr::int(0),
        ));
        assert_eq!(expr_nodes(&e), 1);
        assert_eq!(expr_depth(&e), 0);
    }
}
