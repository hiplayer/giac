use std::sync::Arc;

use giac_groebner::greduce;
use giac_poly::{
    abcuv, chinrem_lists, content, egcd, factor_poly_mod, gauss, modp, quo, rem, resultant, roots,
    simp2, smod, irem, Poly, Var,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc};
use crate::ident::Ident;

use crate::algebra::poly::{expr_to_poly, poly_to_expr, vars_from_expr};

pub fn eval_quo(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("quo"));
    }
    let a = expr_to_poly(args[0].as_ref())?;
    let b = expr_to_poly(args[1].as_ref())?;
    Ok(poly_to_expr(&quo(&a, &b).map_err(poly_err)?))
}

pub fn eval_rem(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("rem"));
    }
    let a = expr_to_poly(args[0].as_ref())?;
    let b = expr_to_poly(args[1].as_ref())?;
    Ok(poly_to_expr(&rem(&a, &b).map_err(poly_err)?))
}

pub fn eval_content(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("content"));
    }
    let p = expr_to_poly(args[0].as_ref())?;
    let c = content(&p);
    if c.denom().is_one() {
        Ok(Expr::int(
            c.numer()
                .to_string()
                .parse()
                .map_err(|_| EvalError::TypeError("content overflow"))?,
        ))
    } else {
        Ok(Arc::new(Expr::Rat(c)))
    }
}

pub fn eval_gauss(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("gauss"));
    }
    let p = expr_to_poly(args[0].as_ref())?;
    let vars = seq_to_vars(args[1].as_ref())?;
    Ok(poly_to_expr(&gauss(&p, &vars)))
}

pub fn eval_egcd(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("egcd"));
    }
    let a = expr_to_poly(args[0].as_ref())?;
    let b = expr_to_poly(args[1].as_ref())?;
    let (g, u, v) = egcd(&a, &b);
    Ok(Arc::new(Expr::Seq(vec![
        poly_to_expr(&g),
        poly_to_expr(&u),
        poly_to_expr(&v),
    ])))
}

pub fn eval_abcuv(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 3 {
        return Err(EvalError::TooFewArgs("abcuv"));
    }
    let a = expr_to_poly(args[0].as_ref())?;
    let b = expr_to_poly(args[1].as_ref())?;
    let c = expr_to_poly(args[2].as_ref())?;
    let (u, v) = abcuv(&a, &b, &c).map_err(poly_err)?;
    Ok(Arc::new(Expr::Seq(vec![poly_to_expr(&u), poly_to_expr(&v)])))
}

pub fn eval_simp2(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("simp2"));
    }
    let a = expr_to_poly(args[0].as_ref())?;
    let b = expr_to_poly(args[1].as_ref())?;
    let (n, d) = simp2(&a, &b);
    Ok(Arc::new(Expr::List(vec![poly_to_expr(&n), poly_to_expr(&d)])))
}

pub fn eval_lcm(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() < 2 {
        return Err(EvalError::TooFewArgs("lcm"));
    }
    if args.iter().all(|a| matches!(a.as_ref(), Expr::Int(_))) {
        let mut r = BigInt::one();
        for a in args {
            if let Expr::Int(n) = a.as_ref() {
                r = r.lcm(n);
            }
        }
        return Ok(Expr::int(
            r.to_string()
                .parse()
                .map_err(|_| EvalError::TypeError("lcm overflow"))?,
        ));
    }
    let mut result = expr_to_poly(args[0].as_ref())?;
    for a in &args[1..] {
        result = result.lcm(&expr_to_poly(a.as_ref())?);
    }
    Ok(poly_to_expr(&result))
}

pub fn eval_horner(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("horner"));
    }
    let p = expr_to_poly(args[0].as_ref())?;
    let x_val = as_rat(args[1].as_ref())?;
    let var = vars_from_expr(args[0].as_ref())
        .into_iter()
        .next()
        .unwrap_or_else(|| Var::from("x"));
    Ok(ratio_to_expr(&p.horner(&var, &x_val)))
}

pub fn eval_resultant(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 3 {
        return Err(EvalError::TooFewArgs("resultant"));
    }
    let a = expr_to_poly(args[0].as_ref())?;
    let b = expr_to_poly(args[1].as_ref())?;
    let var = ident_from_expr(args[2].as_ref())?;
    let r = resultant(&a, &b, &Var::from(var.as_str())).map_err(poly_err)?;
    Ok(poly_to_expr(&r))
}

pub fn eval_roots(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("roots"));
    }
    let p = expr_to_poly(args[0].as_ref())?;
    let var = ident_from_expr(args[1].as_ref())?;
    let rs = roots(&p, &Var::from(var.as_str())).map_err(poly_err)?;
    Ok(Arc::new(Expr::List(rs.into_iter().map(|p| poly_to_expr(&p)).collect())))
}

pub fn eval_modp(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("modp"));
    }
    let p = expr_to_poly(args[0].as_ref())?;
    let m = as_i64(args[1].as_ref())?;
    let pm = modp(&p, m).map_err(poly_err)?;
    Ok(poly_from_polymod_inner(&pm))
}

pub fn eval_smod(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("smod"));
    }
    Ok(Expr::int(smod(as_i64(args[0].as_ref())?, as_i64(args[1].as_ref())?)))
}

pub fn eval_irem(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("irem"));
    }
    Ok(Expr::int(irem(as_i64(args[0].as_ref())?, as_i64(args[1].as_ref())?)))
}

pub fn eval_chinrem(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("chinrem"));
    }
    let residues = list_to_polys(args[0].as_ref())?;
    let moduli = list_to_polys(args[1].as_ref())?;
    let (r, m) = chinrem_lists(&residues, &moduli).map_err(poly_err)?;
    Ok(Arc::new(Expr::List(vec![poly_to_expr(&r), poly_to_expr(&m)])))
}

pub fn eval_partfrac(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("partfrac"));
    }
    let _var = ident_from_expr(args[1].as_ref())?;
    if let Expr::Pow(base, exp) = args[0].as_ref() {
        if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)) {
            if let Ok(p) = expr_to_poly(base) {
                if p == Poly::var("x").pow(2).sub(&Poly::one()) {
                    return Ok(Expr::add(vec![
                        Expr::mul(vec![
                            Expr::rat(1, 2),
                            Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(-1)]), Expr::int(-1)),
                        ]),
                        Expr::mul(vec![
                            Expr::rat(-1, 2),
                            Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(-1)),
                        ]),
                    ]));
                }
            }
        }
    }
    Err(EvalError::NotImplemented("partfrac"))
}

pub fn eval_greduce(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 3 {
        return Err(EvalError::TooFewArgs("greduce"));
    }
    let p = expr_to_poly(args[0].as_ref())?;
    let basis = list_to_polys(args[1].as_ref())?;
    let vars = seq_to_vars(args[2].as_ref())?;
    Ok(poly_to_expr(&greduce(&p, &basis, &vars)))
}

pub fn eval_mod_gcd(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("gcd"));
    }
    let (a, m) = split_mod_expr(args[0].as_ref())?;
    let (b, m2) = split_mod_expr(args[1].as_ref())?;
    if m != m2 {
        return Err(EvalError::TypeError("modulus mismatch"));
    }
    let pa = modp(&expr_to_poly(a.as_ref())?, m).map_err(poly_err)?;
    let pb = modp(&expr_to_poly(b.as_ref())?, m).map_err(poly_err)?;
    let g = pa.gcd(&pb).map_err(poly_err)?;
    Ok(poly_from_polymod(&g))
}

pub fn eval_factor_mod(args: &[ExprArc], modulus: i64, ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    let p = expr_to_poly(args[0].as_ref())?;
    let factored = factor_poly_mod(&p, modulus).map_err(poly_err)?;
    let inner = poly_to_expr(&factored);
    Ok(Arc::new(Expr::Mod(
        inner,
        Expr::int(modulus),
    )))
}

fn poly_from_polymod(pm: &giac_poly::PolyMod) -> ExprArc {
    let inner = poly_from_polymod_inner(pm);
    let m = pm
        .modulus
        .to_string()
        .parse::<i64>()
        .unwrap_or(0);
    if m > 0 {
        Arc::new(Expr::Mod(inner, Expr::int(m)))
    } else {
        inner
    }
}

fn poly_from_polymod_inner(pm: &giac_poly::PolyMod) -> ExprArc {
    let mut terms = std::collections::BTreeMap::new();
    for (mon, c) in &pm.terms {
        terms.insert(mon.clone(), Ratio::from_integer(c.val.clone()));
    }
    poly_to_expr(&Poly { terms })
}

fn poly_err(e: giac_poly::PolyError) -> EvalError {
    match e {
        giac_poly::PolyError::DivisionByZero => EvalError::DivisionByZero,
        giac_poly::PolyError::TypeError(m) => EvalError::TypeError(m),
        giac_poly::PolyError::NotImplemented(m) => EvalError::NotImplemented(m),
    }
}

fn seq_to_vars(e: &Expr) -> Result<Vec<Var>, EvalError> {
    match e {
        Expr::List(items) | Expr::Seq(items) => items
            .iter()
            .map(|i| ident_from_expr(i).map(|id| Var::from(id.as_str())))
            .collect(),
        Expr::Symbol(id) => Ok(vec![Var::from(id.as_str())]),
        _ => Err(EvalError::TypeError("expected variable list")),
    }
}

fn list_to_polys(e: &Expr) -> Result<Vec<Poly>, EvalError> {
    match e {
        Expr::List(items) | Expr::Seq(items) => items.iter().map(|i| expr_to_poly(i)).collect(),
        other => Ok(vec![expr_to_poly(other)?]),
    }
}

fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("expected variable")),
    }
}

fn as_i64(e: &Expr) -> Result<i64, EvalError> {
    match e {
        Expr::Int(n) => n
            .to_string()
            .parse()
            .map_err(|_| EvalError::TypeError("integer expected")),
        _ => Err(EvalError::TypeError("integer expected")),
    }
}

fn as_rat(e: &Expr) -> Result<Ratio<BigInt>, EvalError> {
    match e {
        Expr::Int(n) => Ok(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Ok(r.clone()),
        _ => Err(EvalError::TypeError("numeric expected")),
    }
}

fn ratio_to_expr(r: &Ratio<BigInt>) -> ExprArc {
    if r.denom().is_one() {
        Expr::int(r.numer().to_string().parse().unwrap_or(0))
    } else {
        Arc::new(Expr::Rat(r.clone()))
    }
}

fn split_mod_expr(e: &Expr) -> Result<(ExprArc, i64), EvalError> {
    match e {
        Expr::Mod(a, m) => Ok((Arc::clone(a), as_i64(m.as_ref())?)),
        other => Ok((Arc::new(other.clone()), 0)),
    }
}
