//! Most rapidly varying (MRV) set at `+infinity` (GIAC-216b / `series.cc` subset).

use std::cmp::Ordering;
use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind, Ident};
use num_traits::Signed;

use crate::risch::{depends_on_var, rlvarx};

/// giac `faster_var`, `coeff_ln`, `slower_var`.
#[derive(Clone, Debug, Default)]
pub(crate) struct MrvSet {
    pub faster: Vec<ExprArc>,
    pub coeff_ln: Vec<ExprArc>,
    pub slower: Vec<ExprArc>,
}

impl MrvSet {
    pub(crate) fn is_empty(&self) -> bool {
        self.faster.is_empty()
    }

    pub(crate) fn merge(&mut self, other: &MrvSet, var: &Ident) {
        for (f, c) in other.faster.iter().zip(other.coeff_ln.iter()) {
            merge_mrv_pair(self, Arc::clone(f), Arc::clone(c), var);
        }
        for s in &other.slower {
            if !self.slower.iter().any(|x| x == s) && !self.faster.iter().any(|x| x == s) {
                self.slower.push(Arc::clone(s));
            }
        }
    }
}

/// Build the MRV set for `expr` as `var → +infinity`.
pub(crate) fn mrv_at_plus_infinity(expr: &ExprArc, var: &Ident) -> MrvSet {
    let mut set = MrvSet::default();
    if !depends_on_var(expr, var) {
        return set;
    }
    collect_mrv(expr, var, &mut set);
    if depends_on_var(expr, var) && !set.faster.iter().any(|e| is_var(e, var)) && !set.slower.iter().any(|e| is_var(e, var)) {
        merge_mrv_pair(&mut set, var_to_expr(var), Expr::int(1), var);
    }
    set
}

fn collect_mrv(expr: &ExprArc, var: &Ident, set: &mut MrvSet) {
    if !depends_on_var(expr, var) {
        return;
    }
    match expr.as_ref() {
        Expr::Symbol(id) if id == var => {
            merge_mrv_pair(set, var_to_expr(var), Expr::int(1), var);
        }
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            let inner = &args[0];
            collect_mrv(inner, var, set);
            if depends_on_var(inner, var) {
                merge_mrv_pair(set, Arc::clone(expr), Expr::int(1), var);
            }
        }
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => {
            collect_mrv(&args[0], var, set);
        }
        Expr::Pow(base, exp) => {
            collect_mrv(base, var, set);
            collect_mrv(exp, var, set);
        }
        Expr::Add(ts) | Expr::Mul(ts) => {
            for t in ts {
                collect_mrv(t, var, set);
            }
        }
        Expr::Frac(n, d) => {
            collect_mrv(n, var, set);
            collect_mrv(d, var, set);
        }
        Expr::Func(_, args) => {
            for a in args {
                collect_mrv(a, var, set);
            }
        }
        _ => {}
    }
}

fn merge_mrv_pair(set: &mut MrvSet, elem: ExprArc, coeff_ln: ExprArc, var: &Ident) {
    if set.faster.is_empty() {
        set.faster.push(elem);
        set.coeff_ln.push(coeff_ln);
        return;
    }
    let elem_rank = growth_rank(&elem, var);
    let top_rank = growth_rank(&set.faster[0], var);
    match elem_rank.cmp(&top_rank) {
        Ordering::Greater => {
            set.slower.append(&mut set.faster);
            set.faster = vec![elem];
            set.coeff_ln = vec![coeff_ln];
        }
        Ordering::Equal => {
            if !set.faster.iter().any(|e| e == &elem) {
                set.faster.push(elem);
                set.coeff_ln.push(coeff_ln);
            }
        }
        Ordering::Less => {
            if !set.slower.iter().any(|e| e == &elem) && !set.faster.iter().any(|e| e == &elem) {
                set.slower.push(elem);
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Growth {
    Const,
    Log,
    Poly,
    Exp,
    ExpExp,
}

fn growth_rank(e: &ExprArc, var: &Ident) -> Growth {
    if !depends_on_var(e, var) {
        return Growth::Const;
    }
    match e.as_ref() {
        Expr::Symbol(id) if id == var => Growth::Poly,
        Expr::Pow(b, exp) => {
            if is_var(b, var) {
                if matches!(exp.as_ref(), Expr::Int(n) if n.is_positive()) {
                    return Growth::Poly;
                }
                if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) {
                    return Growth::Poly;
                }
            }
            growth_rank(b, var).max(growth_rank(exp, var))
        }
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => Growth::Log.max(growth_rank(&args[0], var)),
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            let inner = growth_rank(&args[0], var);
            if inner >= Growth::Exp {
                Growth::ExpExp
            } else {
                Growth::Exp
            }
        }
        Expr::Add(ts) | Expr::Mul(ts) => ts.iter().map(|t| growth_rank(t, var)).max().unwrap_or(Growth::Const),
        Expr::Frac(n, d) => growth_rank(n, var).max(growth_rank(d, var)),
        _ => Growth::Poly,
    }
}

/// Pick `w = exp(g)` going to `0` at `+infinity` and express `x` via `w`.
pub(crate) fn choose_mrv_w(set: &MrvSet, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    let mut candidates: Vec<&ExprArc> = set.faster.iter().collect();
    candidates.extend(set.slower.iter());
    let mut best: Option<(&ExprArc, ExprArc)> = None;
    for c in candidates {
        if let Expr::Func(FuncKind::Exp, args) = c.as_ref() {
            if args.len() != 1 {
                continue;
            }
            let inner = &args[0];
            if is_neg_linear(inner, var) {
                let x_expr = Expr::mul(vec![Expr::int(-1), Expr::func(FuncKind::Ln, vec![Arc::clone(c)])]);
                let score = expr_size(c);
                if best.as_ref().map_or(true, |(b, _)| score < expr_size(b)) {
                    best = Some((c, x_expr));
                }
            } else if is_neg_scaled_linear(inner, var) {
                if let Some((a, _)) = linear_coeff(inner, var) {
                    if !a.is_zero() {
                        let x_expr = Expr::mul(vec![
                            Expr::pow(a.clone(), Expr::int(-1)),
                            Expr::func(FuncKind::Ln, vec![Arc::clone(c)]),
                        ]);
                        let score = expr_size(c);
                        if best.as_ref().map_or(true, |(b, _)| score < expr_size(b)) {
                            best = Some((c, x_expr));
                        }
                    }
                }
            }
        }
    }
    best.map(|(w, x)| (Arc::clone(w), x))
}

fn is_neg_linear(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Mul(fs) if fs.len() == 2
        && matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative())
        && is_var(&fs[1], var))
        || matches!(e.as_ref(), Expr::Mul(fs) if fs.len() == 2
        && matches!(fs[1].as_ref(), Expr::Int(n) if n.is_negative())
        && is_var(&fs[0], var))
}

fn is_neg_scaled_linear(e: &ExprArc, var: &Ident) -> bool {
    linear_coeff(e, var).is_some_and(|(a, _)| {
        matches!(a.as_ref(), Expr::Int(n) if n.is_negative())
            || matches!(a.as_ref(), Expr::Frac(num, _) if matches!(num.as_ref(), Expr::Int(n) if n.is_negative()))
    })
}

fn linear_coeff(e: &ExprArc, var: &Ident) -> Option<(ExprArc, ExprArc)> {
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            if is_var(&fs[0], var) {
                Some((Arc::clone(&fs[1]), Expr::int(0)))
            } else if is_var(&fs[1], var) {
                Some((Arc::clone(&fs[0]), Expr::int(0)))
            } else {
                None
            }
        }
        Expr::Symbol(id) if id == var => Some((Expr::int(1), Expr::int(0))),
        Expr::Mul(fs) if fs.iter().any(|f| is_var(f, var)) => {
            let mut coeff = Expr::int(1);
            let mut found = false;
            for f in fs {
                if is_var(f, var) {
                    found = true;
                } else {
                    coeff = Expr::mul(vec![coeff, Arc::clone(f)]);
                }
            }
            if found {
                Some((coeff, Expr::int(0)))
            } else {
                None
            }
        }
        _ => None,
    }
}

fn expr_size(e: &ExprArc) -> usize {
    match e.as_ref() {
        Expr::Int(_) | Expr::Rat(_) | Expr::Symbol(_) => 1,
        Expr::Add(ts) | Expr::Mul(ts) => 1 + ts.iter().map(expr_size).sum::<usize>(),
        Expr::Pow(b, exp) => 1 + expr_size(b) + expr_size(exp),
        Expr::Frac(n, d) => 1 + expr_size(n) + expr_size(d),
        Expr::Func(_, args) => 1 + args.iter().map(expr_size).sum::<usize>(),
        _ => 1,
    }
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn var_to_expr(var: &Ident) -> ExprArc {
    Expr::sym(var.as_str())
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, Expr, FuncKind};

    use super::*;

    #[test]
    fn mrv_exp_neg_x() {
        let var = Ident::new("x");
        let e = Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]);
        let set = mrv_at_plus_infinity(&e, &var);
        assert!(!set.faster.is_empty());
        let (w, x_sub) = choose_mrv_w(&set, &var).expect("w");
        assert!(matches!(w.as_ref(), Expr::Func(FuncKind::Exp, _)));
        assert!(format_expr(x_sub.as_ref()).contains("ln"));
    }
}
