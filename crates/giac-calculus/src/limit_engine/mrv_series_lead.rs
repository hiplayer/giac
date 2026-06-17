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

fn lead_at_zero(expr: &ExprArc, w: &Ident, depth: usize, ctx: &Context) -> Result<LeadTerm, EvalError> {
    if depth > 48 {
        return Err(EvalError::NotImplemented("series"));
    }
    if !depends_on_w(expr, w) {
        return Ok(LeadTerm::constant(eval(expr.as_ref(), ctx)?));
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
    let frac = Arc::new(Expr::Frac(Arc::clone(n), Arc::clone(d)));
    if let Ok(series) = super::sparse_series::series_at_zero(&frac, w, 3, ctx) {
        if let Some((exp, coeff)) = series.lead() {
            return Ok(LeadTerm { exp, coeff });
        }
    }
    let num = lead_at_zero(n, w, depth + 1, ctx)?;
    let den = lead_at_zero(d, w, depth + 1, ctx)?;
    Ok(LeadTerm {
        exp: num.exp.saturating_sub(den.exp),
        coeff: ratnormal(
            Expr::mul(vec![
                Arc::clone(&num.coeff),
                Expr::pow(Arc::clone(&den.coeff), Expr::int(-1)),
            ])
            .as_ref(),
            ctx,
        )
        .ok()
        .unwrap_or_else(|| Expr::mul(vec![num.coeff, Expr::pow(den.coeff, Expr::int(-1))])),
    })
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

/// `exp(f) - exp(g) ~ exp(g)*(f-g)` when the `w^{-1}` parts cancel (CK-INT-61).
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
    if lf.exp != 0 || k != -1 {
        return Ok(None);
    }
    let diff = Expr::add(vec![Arc::clone(f), mrv_ln_w_expr()]);
    let diff_lead = match lead_at_zero(&diff, w, depth + 2, ctx) {
        Ok(lt) => lt,
        Err(_) => second_term_inner_plus_ln(f, w, depth, ctx)?,
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
