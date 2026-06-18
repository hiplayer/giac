//! Most rapidly varying (MRV) set at `+infinity` (GIAC-216b / `series.cc` subset).

use std::cmp::Ordering;
use std::sync::Arc;

use giac_core::{bigint_to_i64, eval, Expr, ExprArc, FuncKind, Ident, Context};
use num_traits::{Signed, ToPrimitive};

use crate::expr_util::depends_on_var;
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

    pub(crate) fn merge(&mut self, other: &MrvSet, var: &Ident, ctx: &Context) {
        for (f, c) in other.faster.iter().zip(other.coeff_ln.iter()) {
            merge_mrv_pair(self, Arc::clone(f), Arc::clone(c), var, ctx);
        }
        for s in &other.slower {
            if !self.slower.iter().any(|x| x == s) && !self.faster.iter().any(|x| x == s) {
                self.slower.push(Arc::clone(s));
            }
        }
    }
}

/// Build the MRV set for `expr` as `var → +infinity`.
pub(crate) fn mrv_at_plus_infinity(expr: &ExprArc, var: &Ident, ctx: &Context) -> MrvSet {
    let mut set = MrvSet::default();
    if !depends_on_var(expr, var) {
        return set;
    }
    collect_mrv(expr, var, &mut set, ctx);
    if depends_on_var(expr, var) && !set.faster.iter().any(|e| is_var(e, var)) && !set.slower.iter().any(|e| is_var(e, var)) {
        merge_mrv_pair(&mut set, var_to_expr(var), Expr::int(1), var, ctx);
    }
    set
}

fn collect_mrv(expr: &ExprArc, var: &Ident, set: &mut MrvSet, ctx: &Context) {
    if !depends_on_var(expr, var) {
        return;
    }
    match expr.as_ref() {
        Expr::Symbol(id) if id == var => {
            merge_mrv_pair(set, var_to_expr(var), Expr::int(1), var, ctx);
        }
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            let inner = &args[0];
            collect_mrv(inner, var, set, ctx);
            if depends_on_var(inner, var) {
                merge_mrv_pair(set, Arc::clone(expr), Expr::int(1), var, ctx);
            }
        }
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => {
            collect_mrv(&args[0], var, set, ctx);
        }
        Expr::Pow(base, exp) => {
            collect_mrv(base, var, set, ctx);
            collect_mrv(exp, var, set, ctx);
            if !depends_on_var(base, var) && depends_on_var(exp, var) {
                merge_mrv_pair(
                    set,
                    Expr::func(
                        FuncKind::Exp,
                        vec![Expr::mul(vec![Arc::clone(exp), Expr::func(FuncKind::Ln, vec![Arc::clone(base)])])],
                    ),
                    Expr::int(1),
                    var,
                    ctx,
                );
            }
        }
        Expr::Add(ts) | Expr::Mul(ts) => {
            for t in ts {
                collect_mrv(t, var, set, ctx);
            }
        }
        Expr::Frac(n, d) => {
            collect_mrv(n, var, set, ctx);
            collect_mrv(d, var, set, ctx);
        }
        Expr::Func(_, args) => {
            for a in args {
                collect_mrv(a, var, set, ctx);
            }
        }
        _ => {}
    }
}

/// Compare growth at `+∞` (upstream `mrv_compare` subset).
pub(crate) fn mrv_compare(a: &ExprArc, b: &ExprArc, var: &Ident, ctx: &Context) -> Ordering {
    let (ga, sa) = growth_rank_detailed(a, var, ctx);
    let (gb, sb) = growth_rank_detailed(b, var, ctx);
    match ga.cmp(&gb) {
        Ordering::Equal => sa.partial_cmp(&sb).unwrap_or(Ordering::Equal),
        other => other,
    }
}

fn growth_rank_detailed(e: &ExprArc, var: &Ident, ctx: &Context) -> (Growth, f64) {
    let g = growth_rank(e, var);
    let sub = match e.as_ref() {
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            linear_coeff_in_var(&args[0], var)
                .and_then(|c| try_const_f64(&c))
                .or_else(|| try_const_f64(&args[0]))
                .unwrap_or_else(|| {
                    eval(args[0].as_ref(), ctx)
                        .ok()
                        .and_then(|v| try_const_f64(&v))
                        .unwrap_or(0.0)
                })
        }
        Expr::Pow(base, exp) if !depends_on_var(base, var) && is_var(exp, var) => {
            try_const_f64(base).map(|v| v.ln()).unwrap_or(0.0)
        }
        Expr::Pow(base, exp) if is_var(base, var) => match exp.as_ref() {
            Expr::Func(FuncKind::Ln, args) if args.len() == 1 && is_var(&args[0], var) => {
                f64::INFINITY
            }
            Expr::Int(n) => bigint_to_i64(n).ok().map(|k| k as f64).unwrap_or(1.0),
            _ => 1.0,
        },
        Expr::Symbol(id) if id == var => 1.0,
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && is_var(&args[0], var) => -1.0,
        _ => 0.0,
    };
    (g, sub)
}

fn merge_mrv_pair(set: &mut MrvSet, elem: ExprArc, coeff_ln: ExprArc, var: &Ident, ctx: &Context) {
    if set.faster.is_empty() {
        set.faster.push(elem);
        set.coeff_ln.push(coeff_ln);
        return;
    }
    match mrv_compare(&elem, &set.faster[0], var, ctx) {
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
            if !depends_on_var(b, var) && depends_on_var(exp, var) {
                return Growth::Exp;
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

pub(crate) fn linear_coeff_in_var(e: &ExprArc, var: &Ident) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Symbol(id) if id == var => Some(Expr::int(1)),
        Expr::Int(_) | Expr::Rat(_) => Some(Expr::int(0)),
        Expr::Mul(fs) => {
            let mut coeff = Expr::int(1);
            let mut has_var = false;
            for f in fs {
                if is_var(f, var) {
                    has_var = true;
                } else if depends_on_var(f, var) {
                    let c = linear_coeff_in_var(f, var)?;
                    coeff = Expr::mul(vec![coeff, c]);
                    has_var = true;
                } else {
                    coeff = Expr::mul(vec![coeff, Arc::clone(f)]);
                }
            }
            if has_var { Some(coeff) } else { None }
        }
        Expr::Add(ts) => {
            let mut sum = Expr::int(0);
            for t in ts {
                let c = linear_coeff_in_var(t, var)?;
                sum = Expr::add(vec![sum, c]);
            }
            Some(sum)
        }
        _ => None,
    }
}

pub(crate) fn is_negative_const_expr(e: &ExprArc, ctx: &Context) -> bool {
    if eval(e.as_ref(), ctx).ok().is_some_and(|v| match v.as_ref() {
        Expr::Int(n) => n.is_negative(),
        Expr::Rat(r) => r.is_negative(),
        Expr::Frac(num, _) => matches!(num.as_ref(), Expr::Int(n) if n.is_negative()),
        _ => false,
    }) {
        return true;
    }
    try_const_f64(e).is_some_and(|x| x < 0.0)
}

pub(crate) fn try_const_f64(e: &ExprArc) -> Option<f64> {
    match e.as_ref() {
        Expr::Int(n) => bigint_to_i64(n).ok().map(|x| x as f64),
        Expr::Rat(r) => Some(r.to_f64()?),
        Expr::Frac(num, den) => {
            let n = try_const_f64(num)?;
            let d = try_const_f64(den)?;
            if d == 0.0 { None } else { Some(n / d) }
        }
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => {
            let v = try_const_f64(&args[0])?;
            if v > 0.0 { Some(v.ln()) } else { None }
        }
        Expr::Add(ts) => {
            let mut sum = 0.0;
            for t in ts {
                sum += try_const_f64(t)?;
            }
            Some(sum)
        }
        Expr::Mul(fs) => {
            let mut prod = 1.0;
            for f in fs {
                prod *= try_const_f64(f)?;
            }
            Some(prod)
        }
        Expr::Pow(b, exp) if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) => {
            let v = try_const_f64(b)?;
            if v == 0.0 { None } else { Some(1.0 / v) }
        }
        _ => None,
    }
}

/// Pick `w = exp(g)` going to `0` at `+infinity` and express `x` via `w`.
pub(crate) fn choose_mrv_w(
    set: &MrvSet,
    var: &Ident,
    ctx: &Context,
) -> Option<(ExprArc, ExprArc)> {
    let mut candidates: Vec<&ExprArc> = set.faster.iter().collect();
    candidates.extend(set.slower.iter());
    let mut best: Option<(&ExprArc, ExprArc)> = None;
    for c in candidates {
        if let Expr::Func(FuncKind::Exp, args) = c.as_ref() {
            if args.len() != 1 {
                continue;
            }
            let inner = &args[0];
            if let Some(coeff) = linear_coeff_in_var(inner, var) {
                if is_negative_const_expr(&coeff, ctx) {
                    let x_expr = Expr::mul(vec![
                        Expr::pow(Arc::clone(&coeff), Expr::int(-1)),
                        Expr::func(FuncKind::Ln, vec![Arc::clone(c)]),
                    ]);
                    let score = expr_size(c);
                    if best.as_ref().map_or(true, |(b, _)| score < expr_size(b)) {
                        best = Some((c, x_expr));
                    }
                    continue;
                }
            }
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
    use std::sync::Arc;

    use giac_core::{format_expr, Expr, FuncKind};

    use crate::limit_engine::preprocess::limit_preprocess_plus_infinity;
    use crate::plugin::xcas_default;

    use super::*;
    #[test]
    fn mrv_nested_exp_inner() {
        let var = Ident::new("x");
        let inner = Arc::new(Expr::Frac(
            Expr::mul(vec![
                Expr::sym("x"),
                Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]),
            ]),
            Expr::add(vec![
                Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]),
                Expr::func(
                    FuncKind::Exp,
                    vec![Expr::mul(vec![
                        Expr::int(-2),
                        Arc::new(Expr::Frac(
                            Expr::pow(Expr::sym("x"), Expr::int(2)),
                            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
                        )),
                    ])],
                ),
            ]),
        ));
        let e = Expr::func(FuncKind::Exp, vec![inner]);
        let ctx = xcas_default();
        let set = mrv_at_plus_infinity(&e, &var, &ctx);
        assert!(!set.faster.is_empty());
        assert!(choose_mrv_w(&set, &var, &ctx).is_some());
    }

    #[test]
    fn choose_mrv_seven_over_eight_pow_n() {
        let ctx = xcas_default();
        let var = Ident::new("n");
        let e = Arc::new(Expr::Frac(
            Expr::pow(Expr::int(7), Expr::sym("n")),
            Expr::pow(Expr::int(8), Expr::sym("n")),
        ));
        let pre = limit_preprocess_plus_infinity(&e, &var, &ctx).unwrap();
        let set = mrv_at_plus_infinity(&pre, &var, &ctx);
        eprintln!(
            "set faster: {:?}",
            set.faster.iter().map(|e| format_expr(e.as_ref())).collect::<Vec<_>>()
        );
        assert!(choose_mrv_w(&set, &var, &ctx).is_some(), "choose_mrv_w failed");
    }

    #[test]
    fn mrv_exp_neg_x() {
        let var = Ident::new("x");
        let e = Expr::func(FuncKind::Exp, vec![Expr::mul(vec![Expr::int(-1), Expr::sym("x")])]);
        let ctx = xcas_default();
        let set = mrv_at_plus_infinity(&e, &var, &ctx);
        assert!(!set.faster.is_empty());
        let (w, x_sub) = choose_mrv_w(&set, &var, &ctx).expect("w");
        assert!(matches!(w.as_ref(), Expr::Func(FuncKind::Exp, _)));
        assert!(format_expr(x_sub.as_ref()).contains("ln"));
    }
}
