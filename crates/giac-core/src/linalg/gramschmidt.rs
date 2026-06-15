//! Gram–Schmidt orthogonalization.

use std::sync::Arc;

use crate::context::Context;
use crate::error::EvalError;
use crate::eval::eval;
use crate::expr::{Expr, ExprArc, FuncKind};
use crate::integrate::integrate;

pub fn eval_gramschmidt(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(if args.is_empty() {
            EvalError::TooFewArgs("gramschmidt")
        } else {
            EvalError::TooManyArgs("gramschmidt")
        });
    }
    let vectors = match args[0].as_ref() {
        Expr::List(v) | Expr::Seq(v) => v.clone(),
        _ => return Err(EvalError::TypeError("gramschmidt expects vector list")),
    };
    let inner = parse_inner_product(&args[1])?;
    gramschmidt_vectors(&vectors, &inner, ctx)
}

type InnerProductFn = Box<dyn Fn(&ExprArc, &ExprArc, &Context) -> Result<ExprArc, EvalError>>;

fn parse_inner_product(arg: &ExprArc) -> Result<InnerProductFn, EvalError> {
    if let Expr::Func(FuncKind::Lambda, parts) = arg.as_ref() {
        if parts.len() == 2 {
            if let Expr::List(params) = parts[0].as_ref() {
                if params.len() == 2 {
                    let body = Arc::clone(&parts[1]);
                    let p_id = match params[0].as_ref() {
                        Expr::Symbol(id) => id.clone(),
                        _ => return Err(EvalError::TypeError("lambda parameter")),
                    };
                    let q_id = match params[1].as_ref() {
                        Expr::Symbol(id) => id.clone(),
                        _ => return Err(EvalError::TypeError("lambda parameter")),
                    };
                    return Ok(Box::new(move |p, q, ctx| {
                        let mut subs = std::collections::HashMap::new();
                        subs.insert(p_id.clone(), Arc::clone(p));
                        subs.insert(q_id.clone(), Arc::clone(q));
                        let substituted = crate::eval::eval_subst_map(&body, &subs)?;
                        eval(substituted.as_ref(), ctx)
                    }));
                }
            }
        }
    }
    // Fallback: recognize integrate inner product pattern without full lambda parse
    Err(EvalError::TypeError("gramschmidt expects (p,q)->inner product"))
}

fn gramschmidt_vectors(
    vectors: &[ExprArc],
    inner: &InnerProductFn,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let n = vectors.len();
    if n == 0 {
        return Ok(Arc::new(Expr::List(vec![])));
    }
    let mut lv = vectors.to_vec();
    let mut norms = Vec::with_capacity(n);
    norms.push(inner(&lv[0], &lv[0], ctx)?);
    for i in 1..n {
        let mut proj = Expr::int(0);
        for j in 0..i {
            let ip = inner(&lv[i], &lv[j], ctx)?;
            let coeff = div_expr(&ip, &norms[j])?;
            proj = Expr::add(vec![proj, Expr::mul(vec![coeff, Arc::clone(&lv[j])])]);
        }
        lv[i] = Expr::add(vec![Arc::clone(&lv[i]), neg_expr(&proj)]);
        norms.push(inner(&lv[i], &lv[i], ctx)?);
    }
    let mut out = Vec::with_capacity(n);
    for (i, v) in lv.into_iter().enumerate() {
        let norm = eval(
            Expr::func(FuncKind::Sqrt, vec![norms[i].clone()]).as_ref(),
            ctx,
        )?;
        out.push(div_expr(&v, &norm)?);
    }
    let list = Arc::new(Expr::List(out));
    eval(list.as_ref(), ctx)
}

fn div_expr(a: &ExprArc, b: &ExprArc) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![Arc::clone(a), Expr::pow(Arc::clone(b), Expr::int(-1))]))
}

fn neg_expr(e: &ExprArc) -> ExprArc {
    Expr::mul(vec![Expr::int(-1), Arc::clone(e)])
}

#[allow(dead_code)]
pub fn default_poly_inner_product(
    p: &ExprArc,
    q: &ExprArc,
    _ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let prod = Expr::mul(vec![Arc::clone(p), Arc::clone(q)]);
    let x = crate::ident::Ident::new("x");
    integrate(&prod, &x)
}
