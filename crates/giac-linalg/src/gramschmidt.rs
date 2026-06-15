//! Gram–Schmidt orthogonalization.

use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use giac_core::Context;
use giac_core::EvalError;
use giac_core::eval;
use giac_core::{Expr, ExprArc, FuncKind};

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
                        let substituted = giac_core::eval_subst_map(&body, &subs)?;
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
    let n0 = inner(&lv[0], &lv[0], ctx)?;
    ensure_non_negative_inner(&n0, ctx)?;
    norms.push(n0);
    for i in 1..n {
        let mut proj: ExprArc = Expr::int(0);
        for j in 0..i {
            let ip = inner(&lv[i], &lv[j], ctx)?;
            let coeff = scalar_div(&ip, &norms[j], ctx)?;
            let term = eval(
                Expr::mul(vec![coeff, Arc::clone(&lv[j])]).as_ref(),
                ctx,
            )?;
            proj = Expr::add(vec![proj, term]);
        }
        lv[i] = eval(
            Expr::add(vec![Arc::clone(&lv[i]), neg_expr(&proj)]).as_ref(),
            ctx,
        )?;
        let ni = inner(&lv[i], &lv[i], ctx)?;
        ensure_non_negative_inner(&ni, ctx)?;
        norms.push(ni);
    }
    let mut out = Vec::with_capacity(n);
    for (i, v) in lv.into_iter().enumerate() {
        let norm = eval(
            Expr::func(FuncKind::Sqrt, vec![norms[i].clone()]).as_ref(),
            ctx,
        )?;
        out.push(eval(div_expr(&v, &norm)?.as_ref(), ctx)?);
    }
    let list = Arc::new(Expr::List(out));
    eval(list.as_ref(), ctx)
}

fn div_expr(a: &ExprArc, b: &ExprArc) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![Arc::clone(a), Expr::pow(Arc::clone(b), Expr::int(-1))]))
}

fn scalar_div(num: &ExprArc, den: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let n = eval(num.as_ref(), ctx)?;
    let d = eval(den.as_ref(), ctx)?;
    if let (Some(a), Some(b)) = (as_scalar(n.as_ref()), as_scalar(d.as_ref())) {
        if b.is_zero() {
            return Err(EvalError::TypeError("division by zero"));
        }
        return scalar_to_expr(a / b);
    }
    eval(div_expr(num, den)?.as_ref(), ctx)
}

fn as_scalar(e: &Expr) -> Option<Ratio<BigInt>> {
    match e {
        Expr::Int(n) => Some(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Some(r.clone()),
        _ => None,
    }
}

fn scalar_to_expr(r: Ratio<BigInt>) -> Result<ExprArc, EvalError> {
    if r.is_zero() {
        return Ok(Expr::int(0));
    }
    if r.denom().is_one() {
        return Ok(Expr::int(
            r.numer()
                .to_string()
                .parse()
                .map_err(|_| EvalError::TypeError("scalar overflow"))?,
        ));
    }
    Ok(Arc::new(Expr::Rat(r)))
}

fn neg_expr(e: &ExprArc) -> ExprArc {
    Expr::mul(vec![Expr::int(-1), Arc::clone(e)])
}

fn ensure_non_negative_inner(norm: &ExprArc, ctx: &Context) -> Result<(), EvalError> {
    let v = eval(norm.as_ref(), ctx)?;
    match v.as_ref() {
        Expr::Int(n) if n.is_negative() => {
            Err(EvalError::TypeError("inner product must be non-negative"))
        }
        Expr::Rat(r) if r.is_negative() => {
            Err(EvalError::TypeError("inner product must be non-negative"))
        }
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use giac_core::eval_subst_map;
    use giac_core::format_expr;

    fn integrate_inner_lambda() -> ExprArc {
        let body = Expr::func(
            FuncKind::Integrate,
            vec![
                Expr::mul(vec![Expr::sym("p"), Expr::sym("q")]),
                Expr::sym("x"),
                Expr::int(-1),
                Expr::int(1),
            ],
        );
        Expr::func(
            FuncKind::Lambda,
            vec![
                Arc::new(Expr::List(vec![Expr::sym("p"), Expr::sym("q")])),
                body,
            ],
        )
    }

    fn call_inner(inner: &ExprArc, p: &ExprArc, q: &ExprArc, ctx: &Context) -> ExprArc {
        let lambda = inner;
        let Expr::Func(FuncKind::Lambda, parts) = lambda.as_ref() else {
            panic!("not lambda");
        };
        let Expr::List(params) = parts[0].as_ref() else {
            panic!("bad params");
        };
        let p_id = match params[0].as_ref() {
            Expr::Symbol(id) => id.clone(),
            _ => panic!("bad p"),
        };
        let q_id = match params[1].as_ref() {
            Expr::Symbol(id) => id.clone(),
            _ => panic!("bad q"),
        };
        let mut subs = HashMap::new();
        subs.insert(p_id, Arc::clone(p));
        subs.insert(q_id, Arc::clone(q));
        let substituted = eval_subst_map(&parts[1], &subs).unwrap();
        eval(substituted.as_ref(), ctx).unwrap()
    }

    #[test]
    fn inner_products_for_gramschmidt_basis() {
        let ctx = crate::plugin::xcas_default();
        let v0 = Expr::int(1);
        let v1 = Expr::add(vec![Expr::int(1), Expr::sym("x")]);
        let lambda = integrate_inner_lambda();
        let n0 = call_inner(&lambda, &v0, &v0, &ctx);
        let cross = call_inner(&lambda, &v1, &v0, &ctx);
        assert_eq!(format_expr(n0.as_ref()), "2");
        assert_eq!(format_expr(cross.as_ref()), "2");
    }

    #[test]
    fn gramschmidt_poly_orthonormal() {
        let ctx = crate::plugin::xcas_default();
        let vectors = vec![
            Expr::int(1),
            Expr::add(vec![Expr::int(1), Expr::sym("x")]),
        ];
        let inner = parse_inner_product(&integrate_inner_lambda()).unwrap();
        let r = gramschmidt_vectors(&vectors, &inner, &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("sqrt"), "got {s}");
        assert!(s.contains("x"), "got {s}");
    }

    #[test]
    fn gramschmidt_rejects_negative_inner_product() {
        let ctx = crate::plugin::xcas_default();
        let inner: InnerProductFn = Box::new(|_, _, _| Ok(Expr::int(-1)));
        let err = gramschmidt_vectors(&[Expr::int(1)], &inner, &ctx).unwrap_err();
        assert!(matches!(err, EvalError::TypeError(_)));
    }
}
