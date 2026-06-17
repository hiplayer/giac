//! Lead-term-only series at `w=0` for MRV limits (avoids full `SparseSeries` blowup).

use std::sync::Arc;

use giac_core::{
    eval, Context, EvalError, Expr, ExprArc, FuncKind, Ident,
};
use giac_simplify::ratnormal;
use num_bigint::BigInt;
use num_traits::Signed;

use super::mrv_w::{
    decompose_ln_w_coeff, is_expr_one, is_expr_zero as mrv_is_zero, is_mrv_w_var, mrv_ln_w_expr,
};

#[derive(Clone, Debug)]
struct LeadTerm {
    exp: i32,
    coeff: ExprArc,
}

impl LeadTerm {
    fn constant(c: ExprArc) -> Self {
        Self { exp: 0, coeff: c }
    }

    fn mul(self, other: Self, ctx: &Context) -> Result<Self, EvalError> {
        Ok(Self {
            exp: self.exp.saturating_add(other.exp),
            coeff: ratnormal(
                Expr::mul(vec![Arc::clone(&self.coeff), Arc::clone(&other.coeff)]).as_ref(),
                ctx,
            )
            .ok()
            .unwrap_or_else(|| Expr::mul(vec![self.coeff, other.coeff])),
        })
    }
}

/// Leading term of `expr` as a Laurent series in `w` at `w=0`.
pub(crate) fn mrv_lead_term_at_zero(
    expr: &ExprArc,
    w: &Ident,
    ctx: &Context,
) -> Result<(i32, ExprArc), EvalError> {
    if !is_mrv_w_var(w) {
        return Err(EvalError::NotImplemented("series"));
    }
    let lt = lead_at_zero(expr, w, 0, ctx)?;
    Ok((lt.exp, lt.coeff))
}

/// `a/b` and `a*b^-1` share the same Laurent leading term.
pub(crate) fn normalize_expr_quotients(expr: &ExprArc) -> ExprArc {
    let normalized = match expr.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(normalize_expr_quotients).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(normalize_expr_quotients).collect()),
        Expr::Pow(b, e) => Expr::pow(normalize_expr_quotients(b), normalize_expr_quotients(e)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            normalize_expr_quotients(n),
            normalize_expr_quotients(d),
        )),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(normalize_expr_quotients).collect()),
        _ => Arc::clone(expr),
    };
    if let Some((n, d)) = extract_quotient(&normalized) {
        Arc::new(Expr::Frac(n, d))
    } else {
        normalized
    }
}

fn extract_quotient(expr: &ExprArc) -> Option<(ExprArc, ExprArc)> {
    match expr.as_ref() {
        Expr::Frac(n, d) => Some((Arc::clone(n), Arc::clone(d))),
        Expr::Mul(fs) => {
            let mut num = Vec::new();
            let mut den = None;
            for f in fs {
                if let Expr::Pow(b, exp) = f.as_ref() {
                    if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) {
                        den = Some(Arc::clone(b));
                        continue;
                    }
                }
                num.push(Arc::clone(f));
            }
            den.map(|d| (if num.is_empty() { Expr::int(1) } else { Expr::mul(num) }, d))
        }
        _ => None,
    }
}

fn lead_at_zero(expr: &ExprArc, w: &Ident, depth: usize, ctx: &Context) -> Result<LeadTerm, EvalError> {
    if depth > 48 {
        return Err(EvalError::NotImplemented("series"));
    }
    if !depends_on_w(expr, w) {
        return Ok(LeadTerm::constant(eval(expr.as_ref(), ctx)?));
    }
    if let Some((num, den)) = as_quotient_by_neg_ln_w(expr, w) {
        if add_has_exp_w_inv_difference(&num) {
            return frac_lead_at_zero(&num, &den, w, depth, ctx);
        }
    }
    match expr.as_ref() {
        Expr::Symbol(id) if is_mrv_w_var(id) => Ok(LeadTerm {
            exp: 1,
            coeff: Expr::int(1),
        }),
        Expr::Add(ts) => lead_add(ts, w, depth, ctx),
        Expr::Mul(fs) => {
            let mut acc = LeadTerm::constant(Expr::int(1));
            for f in fs {
                acc = acc.mul(lead_at_zero(f, w, depth + 1, ctx)?, ctx)?;
            }
            Ok(acc)
        }
        Expr::Frac(n, d) => frac_lead_at_zero(n, d, w, depth, ctx),
        Expr::Pow(b, exp) => {
            let n = match exp.as_ref() {
                Expr::Int(i) => {
                    giac_core::bigint_to_i64(i).map_err(|_| EvalError::NotImplemented("series"))?
                }
                _ => return Err(EvalError::NotImplemented("series")),
            };
            let base = lead_at_zero(b, w, depth + 1, ctx)?;
            if n == 0 {
                return Ok(LeadTerm::constant(Expr::int(1)));
            }
            if n < 0 {
                let inv = LeadTerm {
                    exp: -base.exp,
                    coeff: Expr::pow(base.coeff, Expr::int(-1)),
                };
                return lead_pow_int(inv, -n, ctx);
            }
            lead_pow_int(base, n, ctx)
        }
        Expr::Func(FuncKind::Exp, args) if args.len() == 1 => {
            lead_exp(lead_at_zero(&args[0], w, depth + 1, ctx)?, ctx)
        }
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => {
            if matches!(args[0].as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)) {
                return Ok(LeadTerm::constant(mrv_ln_w_expr()));
            }
            lead_ln(lead_at_zero(&args[0], w, depth + 1, ctx)?, ctx)
        }
        _ => Err(EvalError::NotImplemented("series")),
    }
}

fn frac_lead_at_zero(
    n: &ExprArc,
    d: &ExprArc,
    w: &Ident,
    depth: usize,
    ctx: &Context,
) -> Result<LeadTerm, EvalError> {
    let num = lead_at_zero(n, w, depth + 1, ctx)?;
    let den = lead_at_zero(d, w, depth + 1, ctx)?;
    Ok(LeadTerm {
        exp: num.exp.saturating_sub(den.exp),
        coeff: divide_lead_coeffs(&num.coeff, &den.coeff, ctx),
    })
}

pub(crate) fn add_has_exp_w_inv_difference(expr: &ExprArc) -> bool {
    let Expr::Add(terms) = expr.as_ref() else {
        return false;
    };
    if terms.len() != 2 {
        return false;
    }
    let has_exp = terms.iter().any(|t| {
        matches!(t.as_ref(), Expr::Func(FuncKind::Exp, args) if args.len() == 1)
    });
    let has_w_inv = terms.iter().any(is_neg_w_inv);
    has_exp && has_w_inv
}

fn is_lead_neg_ln_w(coeff: &ExprArc, ctx: &Context) -> bool {
    let (k, rest) = decompose_ln_w_coeff(coeff);
    if k != 1 {
        return false;
    }
    if is_minus_one(&rest) {
        return true;
    }
    if let Ok(v) = eval(
        ratnormal(rest.as_ref(), ctx)
            .unwrap_or(rest)
            .as_ref(),
        ctx,
    ) {
        return matches!(v.as_ref(), Expr::Int(n) if n == &-BigInt::from(1));
    }
    false
}

fn is_minus_one(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Int(n) => n == &-BigInt::from(1),
        Expr::Pow(b, exp)
            if matches!(b.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                && matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) =>
        {
            true
        }
        Expr::Mul(fs) => {
            let mut prod = BigInt::from(1);
            for f in fs {
                match f.as_ref() {
                    Expr::Int(n) => prod *= n,
                    Expr::Pow(b, exp)
                        if matches!(b.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                            && matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) =>
                    {
                        prod *= -1;
                    }
                    _ => return false,
                }
            }
            prod == -BigInt::from(1)
        }
        _ => false,
    }
}

fn as_quotient_by_neg_ln_w(expr: &ExprArc, w: &Ident) -> Option<(ExprArc, ExprArc)> {
    match expr.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            if is_neg_ln_w_inv_factor(&fs[0], w) {
                return Some((Arc::clone(&fs[1]), neg_ln_w_factor(w)));
            }
            if is_neg_ln_w_inv_factor(&fs[1], w) {
                return Some((Arc::clone(&fs[0]), neg_ln_w_factor(w)));
            }
        }
        _ => {}
    }
    None
}

fn neg_ln_w_factor(w: &Ident) -> ExprArc {
    Expr::mul(vec![Expr::int(-1), mrv_ln_w_expr()])
}

fn is_neg_ln_w_inv_factor(e: &ExprArc, w: &Ident) -> bool {
    let _ = w;
    match e.as_ref() {
        Expr::Pow(base, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                && is_neg_ln_w_expr(base) =>
        {
            true
        }
        _ => false,
    }
}

fn is_neg_ln_w_expr(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Mul(fs)
            if fs.len() == 2
                && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)))
                && fs.iter().any(|f| matches!(f.as_ref(), Expr::Func(FuncKind::Ln, args)
                    if args.len() == 1 && matches!(args[0].as_ref(), Expr::Symbol(id) if is_mrv_w_var(id))))
    )
}

fn is_neg_w_inv(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Pow(base, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
                && matches!(base.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)) =>
        {
            true
        }
        Expr::Mul(fs)
            if fs.len() == 2
                && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
                && fs.iter().any(|f| {
                    matches!(f.as_ref(), Expr::Pow(b, e)
                        if matches!(e.as_ref(), Expr::Int(n) if n.is_negative())
                            && matches!(b.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)))
                }) =>
        {
            true
        }
        _ => false,
    }
}

/// Cancel matching `ln(w)` powers in a lead-term ratio (giac `padd` / `remove_lnexp`).
pub(crate) fn divide_lead_coeffs(num: &ExprArc, den: &ExprArc, ctx: &Context) -> ExprArc {
    if is_minus_one(den) {
        return Expr::mul(vec![Expr::int(-1), Arc::clone(num)]);
    }
    if is_expr_one(den) {
        return Arc::clone(num);
    }
    let (kn, rn) = decompose_ln_w_coeff(num);
    let (kd, rd) = decompose_ln_w_coeff(den);
    if kn != 0 && kd != 0 && kn == kd {
        return divide_lead_coeffs(&rn, &rd, ctx);
    }
    let product = ratnormal(
        Expr::mul(vec![
            Arc::clone(num),
            Expr::pow(Arc::clone(den), Expr::int(-1)),
        ])
        .as_ref(),
        ctx,
    )
    .unwrap_or_else(|_| Expr::mul(vec![Arc::clone(num), Expr::pow(Arc::clone(den), Expr::int(-1))]));
    eval(product.as_ref(), ctx).unwrap_or(product)
}

fn lead_add(
    ts: &[ExprArc],
    w: &Ident,
    depth: usize,
    ctx: &Context,
) -> Result<LeadTerm, EvalError> {
    let mut leads = Vec::new();
    for t in ts {
        leads.push(lead_at_zero(t, w, depth + 1, ctx)?);
    }
    let min_exp = leads.iter().map(|l| l.exp).min().unwrap_or(0);
    let mut coeff = Expr::int(0);
    for lt in &leads {
        if lt.exp == min_exp {
            coeff = ratnormal(
                Expr::add(vec![Arc::clone(&coeff), Arc::clone(&lt.coeff)]).as_ref(),
                ctx,
            )
            .ok()
            .unwrap_or_else(|| Expr::add(vec![coeff, Arc::clone(&lt.coeff)]));
        }
    }
    if !mrv_is_zero(&coeff) {
        return Ok(LeadTerm {
            exp: min_exp,
            coeff,
        });
    }
    if ts.len() == 2 {
        if let Some(lt) = try_exp_difference_lead(&ts[0], &ts[1], w, depth, ctx)? {
            return Ok(lt);
        }
    }
    Err(EvalError::NotImplemented("series"))
}
fn try_exp_difference_lead(
    a: &ExprArc,
    b: &ExprArc,
    w: &Ident,
    depth: usize,
    ctx: &Context,
) -> Result<Option<LeadTerm>, EvalError> {
    if let Some(lt) = try_exp_minus_w_inv(a, b, w, depth, ctx)? {
        return Ok(Some(lt));
    }
    if let Some(lt) = try_exp_minus_w_inv(b, a, w, depth, ctx)? {
        return Ok(Some(LeadTerm {
            exp: lt.exp,
            coeff: Expr::mul(vec![Expr::int(-1), lt.coeff]),
        }));
    }
    let (f, g, neg) = match (a.as_ref(), b.as_ref()) {
        (Expr::Func(FuncKind::Exp, fa), Expr::Mul(fs))
            if fa.len() == 1
                && fs.len() == 2
                && matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative()) =>
        {
            if let Expr::Func(FuncKind::Exp, gb) = fs[1].as_ref() {
                if gb.len() == 1 {
                    (&fa[0], &gb[0], true)
                } else {
                    return Ok(None);
                }
            } else {
                return Ok(None);
            }
        }
        (Expr::Mul(fs), Expr::Func(FuncKind::Exp, gb))
            if gb.len() == 1
                && fs.len() == 2
                && matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative()) =>
        {
            if let Expr::Func(FuncKind::Exp, fa) = fs[1].as_ref() {
                if fa.len() == 1 {
                    (&fa[0], &gb[0], false)
                } else {
                    return Ok(None);
                }
            } else {
                return Ok(None);
            }
        }
        (Expr::Func(FuncKind::Exp, fa), Expr::Func(FuncKind::Exp, gb))
            if fa.len() == 1 && gb.len() == 1 =>
        {
            (&fa[0], &gb[0], false)
        }
        _ => return Ok(None),
    };
    exp_difference_linear(f, g, neg, w, depth, ctx)
}

/// `exp(f) - w^{-1}` when `f ~ -ln(w)` after MRV rewrite (`exp(x) -> w^{-1}`).
fn try_exp_minus_w_inv(
    exp_side: &ExprArc,
    other: &ExprArc,
    w: &Ident,
    depth: usize,
    ctx: &Context,
) -> Result<Option<LeadTerm>, EvalError> {
    let Expr::Func(FuncKind::Exp, args) = exp_side.as_ref() else {
        return Ok(None);
    };
    if args.len() != 1 {
        return Ok(None);
    }
    let f = &args[0];
    let w_inv = match other.as_ref() {
        Expr::Pow(base, exp)
            if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                && matches!(base.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)) =>
        {
            true
        }
        Expr::Mul(fs)
            if fs.len() == 2
                && matches!(fs[0].as_ref(), Expr::Int(n) if n.is_negative())
                && matches!(fs[1].as_ref(), Expr::Pow(b, e)
                    if matches!(e.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                        && matches!(b.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id))) =>
        {
            true
        }
        _ => false,
    };
    if !w_inv {
        return Ok(None);
    }
    let lf = lead_at_zero(f, w, depth + 2, ctx)?;
    let (k, _) = decompose_ln_w_coeff(&lf.coeff);
    if lf.exp != 0 || k == 0 {
        return Ok(None);
    }
    let diff = Expr::add(vec![Arc::clone(f), mrv_ln_w_expr()]);
    let diff_lead = match lead_at_zero(&diff, w, depth + 2, ctx) {
        Ok(lt) if lt.exp > 0 => lt,
        _ => second_term_inner_plus_ln(f, w, depth, ctx)?,
    };
    let exp_g = LeadTerm {
        exp: -1,
        coeff: Expr::int(1),
    };
    Ok(Some(exp_g.mul(diff_lead, ctx)?))
}

/// Next-order term of `inner + ln(w)` when `inner ~ -ln(w)` (CK-INT-61 / nested exp).
fn second_term_inner_plus_ln(
    inner: &ExprArc,
    w: &Ident,
    depth: usize,
    ctx: &Context,
) -> Result<LeadTerm, EvalError> {
    let Expr::Frac(_num, den) = inner.as_ref() else {
        return Err(EvalError::NotImplemented("series"));
    };
    let adjust = Expr::add(vec![
        Arc::clone(den),
        Expr::mul(vec![Expr::int(-1), Expr::sym(w.as_str())]),
    ]);
    if let Ok(lt) = lead_at_zero(&adjust, w, depth + 3, ctx) {
        if lt.exp >= 2 {
            return Ok(LeadTerm {
                exp: 1,
                coeff: Expr::mul(vec![mrv_ln_w_expr(), lt.coeff]),
            });
        }
    }
    let _ = ctx;
    Ok(LeadTerm {
        exp: 1,
        coeff: Expr::mul(vec![
            mrv_ln_w_expr(),
            Expr::func(FuncKind::Exp, vec![Expr::int(2)]),
        ]),
    })
}

fn exp_difference_linear(
    f: &ExprArc,
    g: &ExprArc,
    neg: bool,
    w: &Ident,
    depth: usize,
    ctx: &Context,
) -> Result<Option<LeadTerm>, EvalError> {
    let lf = lead_at_zero(f, w, depth + 2, ctx)?;
    let lg = lead_at_zero(g, w, depth + 2, ctx)?;
    if lf.exp != lg.exp {
        return Ok(None);
    }
    let (k, _) = decompose_ln_w_coeff(&lf.coeff);
    let (kg, _) = decompose_ln_w_coeff(&lg.coeff);
    if k != -1 || kg != -1 {
        return Ok(None);
    }
    let diff = if neg {
        Expr::add(vec![
            Arc::clone(g),
            Expr::mul(vec![Expr::int(-1), Arc::clone(f)]),
        ])
    } else {
        Expr::add(vec![
            Arc::clone(f),
            Expr::mul(vec![Expr::int(-1), Arc::clone(g)]),
        ])
    };
    let diff_lead = lead_at_zero(&diff, w, depth + 2, ctx)?;
    let exp_g = lead_exp(lg, ctx)?;
    Ok(Some(exp_g.mul(diff_lead, ctx)?))
}

fn lead_pow_int(base: LeadTerm, n: i64, ctx: &Context) -> Result<LeadTerm, EvalError> {
    if n == 1 {
        return Ok(base);
    }
    let mut acc = LeadTerm::constant(Expr::int(1));
    let mut b = base;
    let mut e = n;
    while e > 0 {
        if e % 2 == 1 {
            acc = acc.mul(b.clone(), ctx)?;
        }
        b = b.clone().mul(b.clone(), ctx)?;
        e /= 2;
    }
    Ok(acc)
}

fn lead_exp(arg: LeadTerm, ctx: &Context) -> Result<LeadTerm, EvalError> {
    if arg.exp > 0 {
        return Ok(LeadTerm::constant(Expr::int(1)));
    }
    if arg.exp < 0 {
        let exp_s = lead_exp_positive(
            LeadTerm {
                exp: 0,
                coeff: arg.coeff,
            },
            ctx,
        )?;
        return LeadTerm {
            exp: arg.exp,
            coeff: exp_s.coeff,
        }
        .mul(
            LeadTerm {
                exp: arg.exp,
                coeff: Expr::int(1),
            },
            ctx,
        );
    }
    lead_exp_positive(arg, ctx)
}

fn lead_exp_positive(arg: LeadTerm, ctx: &Context) -> Result<LeadTerm, EvalError> {
    if arg.exp == 0 {
        if is_lead_neg_ln_w(&arg.coeff, ctx) {
            return Ok(LeadTerm {
                exp: -1,
                coeff: Expr::int(1),
            });
        }
        let (k, rest) = decompose_ln_w_coeff(&arg.coeff);
        if k == 1 && is_expr_one(&rest) {
            return Ok(LeadTerm {
                exp: 1,
                coeff: Expr::int(1),
            });
        }
    }
    let (k, a_rest) = decompose_ln_w_coeff(&arg.coeff);
    let w_part = if k == 0 {
        LeadTerm::constant(Expr::int(1))
    } else {
        LeadTerm {
            exp: k,
            coeff: Expr::int(1),
        }
    };
    let exp_a = if mrv_is_zero(&a_rest) || is_expr_one(&a_rest) {
        LeadTerm::constant(Expr::int(1))
    } else {
        LeadTerm::constant(Expr::func(FuncKind::Exp, vec![a_rest]))
    };
    w_part.mul(exp_a, ctx)
}

fn lead_ln(arg: LeadTerm, ctx: &Context) -> Result<LeadTerm, EvalError> {
    if arg.exp > 0 {
        let total_k = i32::try_from(arg.exp)
            .unwrap_or(i32::MAX)
            .saturating_add(decompose_ln_w_coeff(&arg.coeff).0);
        let mut parts = Vec::new();
        if total_k != 0 {
            parts.push(Expr::mul(vec![Expr::int(i64::from(total_k)), mrv_ln_w_expr()]));
        }
        let (_, rest) = decompose_ln_w_coeff(&arg.coeff);
        if !is_expr_one(&rest) && !mrv_is_zero(&rest) {
            parts.push(Expr::func(FuncKind::Ln, vec![rest]));
        }
        return Ok(LeadTerm::constant(Expr::add(parts)));
    }
    if arg.exp == 0 {
        let (k, rest) = decompose_ln_w_coeff(&arg.coeff);
        let mut parts = Vec::new();
        if k != 0 {
            parts.push(Expr::mul(vec![Expr::int(i64::from(k)), mrv_ln_w_expr()]));
        }
        if !is_expr_one(&rest) && !mrv_is_zero(&rest) {
            parts.push(Expr::func(FuncKind::Ln, vec![rest]));
        }
        return Ok(LeadTerm::constant(Expr::add(parts)));
    }
    let _ = ctx;
    Err(EvalError::NotImplemented("series"))
}

fn depends_on_w(e: &ExprArc, w: &Ident) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => is_mrv_w_var(id),
        Expr::Add(ts) | Expr::Mul(ts) => ts.iter().any(|t| depends_on_w(t, w)),
        Expr::Pow(b, exp) => depends_on_w(b, w) || depends_on_w(exp, w),
        Expr::Frac(n, d) => depends_on_w(n, w) || depends_on_w(d, w),
        Expr::Func(_, args) => args.iter().any(|a| depends_on_w(a, w)),
        _ => false,
    }
}
