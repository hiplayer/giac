//! Most rapidly varying (MRV) set at `+infinity` (GIAC-216b / `series.cc` subset).
//!
//! **API 分层：** [`giac-calculus-api-stability.md`](../../../../../.doc/giac-calculus-api-stability.md)
//! **专项契约：** [`limit-engine-expr-api.md`](../../../../../.doc/limit-engine-expr-api.md)
//!
//! | 层级 | 内容 |
//! |------|------|
//! | **Stable** | `mrv_at_plus_infinity`、`mrv_compare`、`vanishes_faster_than_*`、`choose_mrv_w` |
//! | **Pipeline private** | `collect_mrv`、`growth_rank*`、`merge_mrv_pair` 等遍历/合并 |

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::cmp::Ordering;
use std::sync::Arc;

use giac_core::{bigint_to_i64, Expr, ExprArc, FuncKind, Ident, Context};
use num_traits::{Signed, ToPrimitive};

use crate::expr_util::depends_on_var;
use super::simplify_util::{int_pow_growth_sub_rank, is_negative_const_expr as is_neg_const};
#[derive(Clone, Debug, Default)]
pub(crate) struct MrvSet {
    pub faster: Vec<ExprArc>,
    pub coeff_ln: Vec<ExprArc>,
    pub slower: Vec<ExprArc>,
}

impl MrvSet {
    // **Pipeline private** — MRV 集合是否为空
    pub(crate) fn is_empty(&self) -> bool {
        self.faster.is_empty()
    }

    // **Pipeline private** — 合并两个 MRV 集合
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

/// **Stable** — 构建 `+∞` MRV 集合
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

// **Pipeline private** — collect mrv
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

/// **Stable** — 谓词：`n/d→0` at `+∞`
pub(crate) fn ratio_tends_to_zero_at_plus_infinity(
    n: &ExprArc,
    d: &ExprArc,
    var: &Ident,
) -> bool {
    vanishes_faster_than_at_plus_infinity(n, d, var, &Context::default())
}

/// **Stable** — 谓词：`a/b→0` at `+∞`
pub(crate) fn vanishes_faster_than_at_plus_infinity(
    a: &ExprArc,
    b: &ExprArc,
    var: &Ident,
    ctx: &Context,
) -> bool {
    match mrv_compare(a, b, var, ctx) {
        Ordering::Less => true,
        Ordering::Equal => {
            let (_, sa) = growth_rank_detailed(a, var, ctx);
            let (_, sb) = growth_rank_detailed(b, var, ctx);
            sa < sb
        }
        Ordering::Greater => false,
    }
}

/// **Stable** — `+∞` 增长比较（上游 `mrv_compare` 子集）
pub(crate) fn mrv_compare(a: &ExprArc, b: &ExprArc, var: &Ident, ctx: &Context) -> Ordering {
    let (ga, sa) = growth_rank_detailed(a, var, ctx);
    let (gb, sb) = growth_rank_detailed(b, var, ctx);
    match ga.cmp(&gb) {
        Ordering::Equal => sa.partial_cmp(&sb).unwrap_or(Ordering::Equal),
        other => other,
    }
}

// **Pipeline private** — growth rank detailed
fn growth_rank_detailed(e: &ExprArc, var: &Ident, ctx: &Context) -> (Growth, f64) {
    match e.as_ref() {
        Expr::Add(ts) => ts
            .iter()
            .filter(|t| depends_on_var(t, var))
            .map(|t| growth_rank_detailed(t, var, ctx))
            .max_by(|a, b| match a.0.cmp(&b.0) {
                Ordering::Equal => a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal),
                other => other,
            })
            .unwrap_or((Growth::Const, 0.0)),
        Expr::Mul(fs) => {
            let mut best = (Growth::Const, 0.0);
            let mut exp_sub_sum = 0.0;
            let mut exp_count = 0usize;
            for f in fs {
                let g = growth_rank_detailed(f, var, ctx);
                if g.0 == Growth::Exp {
                    exp_sub_sum += g.1;
                    exp_count += 1;
                }
                if g.0 > best.0 || (g.0 == best.0 && g.1 > best.1) {
                    best = g;
                }
            }
            if exp_count > 0 && best.0 == Growth::Exp {
                (Growth::Exp, exp_sub_sum)
            } else {
                best
            }
        }
        _ => {
            let g = growth_rank(e, var);
            let sub = match e.as_ref() {
                Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
                    asymptotic_linear_coeff_at_plus_infinity(&args[0], var, ctx)
                        .or_else(|| {
                            linear_coeff_in_var(&args[0], var).and_then(|c| try_const_f64(&c))
                        })
                        .or_else(|| try_const_f64(&args[0]))
                        .unwrap_or(0.0)
                }
                Expr::Pow(base, exp) if !depends_on_var(base, var) && is_var(exp, var) => {
                    int_pow_growth_sub_rank(base)
                        .map(|r| r as f64)
                        .or_else(|| try_const_f64(base).map(|v| v.ln()))
                        .unwrap_or(0.0)
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
    }
}

/// **Stable** — 常数底幂商在 `+∞` 的极限
pub(crate) fn limit_const_pow_quotient_at_plus_infinity(
    expr: &ExprArc,
    var: &Ident,
) -> Option<ExprArc> {
    let (num, den) = super::mrv_series_lead::try_as_quotient(expr)?;
    let Expr::Pow(nb, ne) = num.as_ref() else {
        return None;
    };
    let Expr::Pow(db, de) = den.as_ref() else {
        return None;
    };
    if !matches!(ne.as_ref(), Expr::Symbol(id) if id == var) {
        return None;
    }
    if !matches!(de.as_ref(), Expr::Symbol(id) if id == var) {
        return None;
    }
    let nc = try_const_f64(nb)?;
    let dc = try_const_f64(db)?;
    if nc <= 0.0 || dc <= 0.0 {
        return None;
    }
    if (nc - dc).abs() < f64::EPSILON {
        return Some(Expr::int(1));
    }
    if nc < dc {
        return Some(Expr::int(0));
    }
    Some(Expr::sym("+infinity"))
}

// **Pipeline private** — asymptotic linear coeff at plus infinity
fn asymptotic_linear_coeff_at_plus_infinity(
    e: &ExprArc,
    var: &Ident,
    _ctx: &Context,
) -> Option<f64> {
    match e.as_ref() {
        Expr::Frac(n, d) => {
            let nd = add_degree_in_var(n, var);
            let dd = add_degree_in_var(d, var);
            if nd == dd + 1 {
                let nc = leading_coeff_in_var(n, var)?;
                let dc = leading_coeff_in_var(d, var)?;
                if dc.abs() > f64::EPSILON {
                    return Some(nc / dc);
                }
            }
            None
        }
        Expr::Mul(fs) => {
            let mut scalar = 1.0;
            let mut var_coeff = None;
            for f in fs {
                if let Some(c) = asymptotic_linear_coeff_at_plus_infinity(f, var, _ctx) {
                    if var_coeff.is_some() {
                        return None;
                    }
                    var_coeff = Some(c);
                } else if let Some(c) = try_const_f64(f) {
                    scalar *= c;
                } else if !depends_on_var(f, var) {
                    scalar *= try_const_f64(f).unwrap_or(1.0);
                } else {
                    return None;
                }
            }
            var_coeff.map(|c| c * scalar)
        }
        _ => None,
    }
}

// **Pipeline private** — leading coeff in var
fn leading_coeff_in_var(e: &ExprArc, var: &Ident) -> Option<f64> {
    let deg = add_degree_in_var(e, var);
    match e.as_ref() {
        Expr::Symbol(id) if id == var => Some(1.0),
        Expr::Pow(b, exp) if is_var(b, var) => {
            let k = match exp.as_ref() {
                Expr::Int(n) => bigint_to_i64(n).ok()? as isize,
                _ => return None,
            };
            (k == deg).then_some(1.0)
        }
        Expr::Mul(fs) => {
            let mut p = 1.0;
            for f in fs {
                if depends_on_var(f, var) {
                    p *= leading_coeff_in_var(f, var)?;
                } else {
                    p *= try_const_f64(f).unwrap_or(1.0);
                }
            }
            Some(p)
        }
        Expr::Add(ts) => {
            let mut sum = 0.0;
            let mut any = false;
            for t in ts {
                if add_degree_in_var(t, var) == deg {
                    sum += leading_coeff_in_var(t, var).unwrap_or(0.0);
                    any = true;
                }
            }
            any.then_some(sum)
        }
        _ => try_const_f64(e),
    }
}

/// **Stable** — 有理/常数 lead 在 `+∞`
pub(crate) fn limit_rational_const_at_plus_infinity(
    n: &ExprArc,
    d: &ExprArc,
    var: &Ident,
) -> Option<f64> {
    let nd = add_degree_in_var(n, var);
    let dd = add_degree_in_var(d, var);
    if nd == dd {
        let nc = leading_coeff_in_var(n, var)?;
        let dc = leading_coeff_in_var(d, var)?;
        if dc.abs() > f64::EPSILON {
            return Some(nc / dc);
        }
    }
    if nd < dd {
        return Some(0.0);
    }
    None
}

// **Pipeline private** — add degree in var
fn add_degree_in_var(e: &ExprArc, var: &Ident) -> isize {
    match e.as_ref() {
        Expr::Pow(b, exp) if is_var(b, var) => match exp.as_ref() {
            Expr::Int(n) => bigint_to_i64(n).ok().unwrap_or(1) as isize,
            _ => 1,
        },
        Expr::Mul(fs) => fs.iter().map(|f| add_degree_in_var(f, var)).sum(),
        Expr::Symbol(id) if id == var => 1,
        Expr::Add(ts) => ts.iter().map(|t| add_degree_in_var(t, var)).max().unwrap_or(0),
        _ if !depends_on_var(e, var) => 0,
        _ => 0,
    }
}

// **Pipeline private** — merge mrv pair
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

// **Pipeline private** — growth rank
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

/// **Stable** — 提取 `var` 的线性系数（含 Frac 指数）
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
        Expr::Frac(n, d) => {
            if depends_on_var(d, var) {
                return None;
            }
            let cn = linear_coeff_in_var(n, var)?;
            if super::mrv_w::is_expr_one(d) {
                Some(cn)
            } else if !depends_on_var(d, var) {
                Some(Arc::new(Expr::Frac(cn, Arc::clone(d))))
            } else {
                None
            }
        }
        _ => None,
    }
}

/// **Stable** — 符号负常数谓词
pub(crate) fn is_negative_const_expr(e: &ExprArc, _ctx: &Context) -> bool {
    if is_neg_const(e) {
        return true;
    }
    try_const_f64(e).is_some_and(|x| x < 0.0)
}

/// **Stable** — 尝试将表达式读为 `f64` 常数
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

/// **Stable** — 从 MRV 集选取换元 `w`
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

// **Pipeline private** — is neg linear
fn is_neg_linear(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Mul(fs) if fs.len() == 2
        && matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative())
        && is_var(&fs[1], var))
        || matches!(e.as_ref(), Expr::Mul(fs) if fs.len() == 2
        && matches!(fs[1].as_ref(), Expr::Int(n) if n.is_negative())
        && is_var(&fs[0], var))
}

// **Pipeline private** — is neg scaled linear
fn is_neg_scaled_linear(e: &ExprArc, var: &Ident) -> bool {
    linear_coeff(e, var).is_some_and(|(a, _)| {
        matches!(a.as_ref(), Expr::Int(n) if n.is_negative())
            || matches!(a.as_ref(), Expr::Frac(num, _) if matches!(num.as_ref(), Expr::Int(n) if n.is_negative()))
    })
}

// **Pipeline private** — linear coeff
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

// **Pipeline private** — expr size
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

// **Pipeline private** — is var
fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

// **Pipeline private** — var to expr
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
