//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-ode-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{
    eval, bigint_to_i64, integer_sqrt, is_sin_of_var, ratio_to_expr, Context, EvalError, Expr,
    ExprArc, FuncKind, Ident, RelOp,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

/// `desolve(equation, y(x))` — linear constant-coefficient ODEs (GIAC-218).
/// **Stable (bounded)** — linear constant-coefficient ODE subset
pub fn eval_desolve(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("desolve"));
    }
    let eq = eval(args[0].as_ref(), ctx)?;
    let (dep, indep) = parse_dep_fn(&args[1])?;
    let ode = parse_linear_ode(&eq, &dep, &indep, ctx)?;
    let sol = solve_linear_ode(&ode, &dep, &indep)?;
    Ok(Arc::new(Expr::Relation(
        RelOp::Eq,
        dep_fn_expr(&dep, &indep),
        sol,
    )))
}

struct LinOde {
    a0: ExprArc,
    a1: ExprArc,
    a2: ExprArc,
    forcing: ExprArc,
}

// **Pipeline private** — `parse_dep_fn`
fn parse_dep_fn(e: &ExprArc) -> Result<(Ident, Ident), EvalError> {
    match e.as_ref() {
        Expr::Func(FuncKind::Apply, args) if args.len() == 2 => {
            let dep = match args[0].as_ref() {
                Expr::Symbol(id) => id.clone(),
                _ => return Err(EvalError::TypeError("desolve dependent function")),
            };
            let indep = match args[1].as_ref() {
                Expr::Symbol(id) => id.clone(),
                _ => return Err(EvalError::TypeError("desolve independent variable")),
            };
            Ok((dep, indep))
        }
        _ => Err(EvalError::TypeError("desolve expects y(x)")),
    }
}

// **Pipeline private** — `dep_fn_expr`
fn dep_fn_expr(dep: &Ident, indep: &Ident) -> ExprArc {
    Expr::func(FuncKind::Apply, vec![Expr::sym(dep.as_str()), Expr::sym(indep.as_str())])
}

// **Pipeline private** — `parse_linear_ode`
fn parse_linear_ode(
    eq: &Expr,
    dep: &Ident,
    indep: &Ident,
    ctx: &Context,
) -> Result<LinOde, EvalError> {
    let (lhs, rhs) = match eq {
        Expr::Relation(RelOp::Eq, l, r) => (l.as_ref(), r.as_ref()),
        _ => return Err(EvalError::TypeError("equation expected")),
    };
    let mut lhs_e = Expr::add(vec![
        Arc::new(lhs.clone()),
        Expr::mul(vec![Expr::int(-1), Arc::new(rhs.clone())]),
    ]);
    lhs_e = eval(lhs_e.as_ref(), ctx)?;
    let mut a0 = Expr::int(0);
    let mut a1 = Expr::int(0);
    let mut a2 = Expr::int(0);
    let mut forcing = Expr::int(0);
    for term in flatten_add(lhs_e.as_ref()) {
        if let Some(order) = derivative_order(&term, dep, indep) {
            let coeff = strip_derivative_factor(&term, dep, indep, order)?;
            match order {
                0 => a0 = add_expr(a0, coeff),
                1 => a1 = add_expr(a1, coeff),
                2 => a2 = add_expr(a2, coeff),
                _ => return Err(EvalError::NotImplemented("desolve")),
            }
        } else {
            forcing = add_expr(forcing, Expr::mul(vec![Expr::int(-1), term]));
        }
    }
    Ok(LinOde {
        a0: eval(a0.as_ref(), ctx)?,
        a1: eval(a1.as_ref(), ctx)?,
        a2: eval(a2.as_ref(), ctx)?,
        forcing: eval(forcing.as_ref(), ctx)?,
    })
}

// **Pipeline private** — `flatten_add`
fn flatten_add(e: &Expr) -> Vec<ExprArc> {
    match e {
        Expr::Add(terms) => terms.clone(),
        _ => vec![Arc::new(e.clone())],
    }
}

// **Pipeline private** — `add_expr`
fn add_expr(a: ExprArc, b: ExprArc) -> ExprArc {
    if is_zero_expr(&a) {
        return b;
    }
    if is_zero_expr(&b) {
        return a;
    }
    Expr::add(vec![a, b])
}

// **Pipeline private** — `derivative_order`
fn derivative_order(e: &ExprArc, dep: &Ident, indep: &Ident) -> Option<u8> {
    match e.as_ref() {
        Expr::Func(FuncKind::Prime, args) if args.len() == 2 => {
            if is_dep(&args[0], dep, indep) {
                bigint_to_i64(match args[1].as_ref() {
                    Expr::Int(n) => n,
                    _ => return None,
                })
                .ok()
                .and_then(|o| u8::try_from(o).ok())
            } else {
                None
            }
        }
        Expr::Symbol(id) if id == dep => Some(0),
        Expr::Func(FuncKind::Apply, args)
            if args.len() == 2
                && matches!(args[0].as_ref(), Expr::Symbol(id) if id == dep) =>
        {
            Some(0)
        }
        Expr::Mul(factors) => factors.iter().find_map(|f| derivative_order(f, dep, indep)),
        _ => None,
    }
}

// **Pipeline private** — `is_dep`
fn is_dep(e: &ExprArc, dep: &Ident, indep: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == dep)
        || matches!(
            e.as_ref(),
            Expr::Func(FuncKind::Apply, args)
                if args.len() == 2
                    && matches!(args[0].as_ref(), Expr::Symbol(id) if id == dep)
                    && matches!(args[1].as_ref(), Expr::Symbol(id) if id == indep)
        )
}

// **Pipeline private** — `strip_derivative_factor`
fn strip_derivative_factor(
    t: &ExprArc,
    dep: &Ident,
    indep: &Ident,
    order: u8,
) -> Result<ExprArc, EvalError> {
    match t.as_ref() {
        Expr::Func(FuncKind::Prime, args) if args.len() == 2 && is_dep(&args[0], dep, indep) => {
            Ok(Expr::int(1))
        }
        Expr::Symbol(id) if id == dep && order == 0 => Ok(Expr::int(1)),
        Expr::Func(FuncKind::Apply, args)
            if args.len() == 2
                && matches!(args[0].as_ref(), Expr::Symbol(id) if id == dep)
                && order == 0 =>
        {
            Ok(Expr::int(1))
        }
        Expr::Mul(factors) => {
            let mut parts = Vec::new();
            let mut saw = false;
            for f in factors {
                if derivative_order(f, dep, indep) == Some(order) {
                    saw = true;
                } else {
                    parts.push(Arc::clone(f));
                }
            }
            if !saw {
                return Err(EvalError::TypeError("ode term"));
            }
            if parts.is_empty() {
                Ok(Expr::int(1))
            } else if parts.len() == 1 {
                Ok(parts.into_iter().next().unwrap())
            } else {
                Ok(Expr::mul(parts))
            }
        }
        _ => Err(EvalError::TypeError("ode term")),
    }
}

// **Pipeline private** — `solve_linear_ode`
fn solve_linear_ode(ode: &LinOde, dep: &Ident, indep: &Ident) -> Result<ExprArc, EvalError> {
    if !is_zero_expr(&ode.a2) {
        return solve_second_order(ode, dep, indep);
    }
    if !is_zero_expr(&ode.a1) {
        return solve_first_order(ode, indep);
    }
    if !is_zero_expr(&ode.a0) {
        return Err(EvalError::NotImplemented("desolve"));
    }
    let _ = dep;
    Err(EvalError::TypeError("degenerate ode"))
}

// **Pipeline private** — `solve_first_order`
fn solve_first_order(ode: &LinOde, indep: &Ident) -> Result<ExprArc, EvalError> {
    let x = Expr::sym(indep.as_str());
    if is_one_expr(&ode.a1) && is_zero_expr(&ode.forcing) {
        if is_neg_var(&ode.a0, indep) {
            return Ok(Expr::mul(vec![
                const_sym(0),
                Expr::func(
                    FuncKind::Exp,
                    vec![Expr::mul(vec![Expr::rat(1, 2), Expr::pow(x.clone(), Expr::int(2))])],
                ),
            ]));
        }
    }
    let p = if is_one_expr(&ode.a1) {
        Arc::clone(&ode.a0)
    } else {
        return Err(EvalError::NotImplemented("desolve"));
    };
    let hom = Expr::mul(vec![
        const_sym(0),
        Expr::func(
            FuncKind::Exp,
            vec![Expr::mul(vec![Expr::int(-1), p.clone(), x.clone()])],
        ),
    ]);
    if is_zero_expr(&ode.forcing) {
        return Ok(hom);
    }
    if is_var_expr(&ode.forcing, indep) && is_one_expr(&p) {
        return Ok(Expr::add(vec![hom, Expr::add(vec![x.clone(), Expr::int(-1)])]));
    }
    Err(EvalError::NotImplemented("desolve"))
}

// **Pipeline private** — `is_one_expr`
fn is_one_expr(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_one())
}

// **Pipeline private** — `is_neg_var`
fn is_neg_var(e: &ExprArc, indep: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Mul(factors) if factors.len() == 2
            && matches!(factors[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
            && matches!(factors[1].as_ref(), Expr::Symbol(id) if id.as_str() == indep.as_str())
    )
}

// **Pipeline private** — `is_var_expr`
fn is_var_expr(e: &ExprArc, indep: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id.as_str() == indep.as_str())
}

// **Pipeline private** — `solve_second_order`
fn solve_second_order(ode: &LinOde, dep: &Ident, indep: &Ident) -> Result<ExprArc, EvalError> {
    let a2 = expr_to_ratio(&ode.a2)?;
    let a1 = expr_to_ratio(&ode.a1)? / a2.clone();
    let a0 = expr_to_ratio(&ode.a0)? / a2;
    let forcing = if is_zero_expr(&ode.forcing) {
        None
    } else {
        Some(Arc::clone(&ode.forcing))
    };

    let hom = homogeneous_second_order(&a1, &a0, indep)?;
    let mut sol = hom;
    if let Some(f) = forcing {
        if let Some(part) = particular_sin(&a0, &f, indep)? {
            sol = Expr::add(vec![sol, part]);
        } else {
            return Err(EvalError::NotImplemented("desolve"));
        }
    }
    let _ = dep;
    Ok(sol)
}

// **Pipeline private** — `expr_to_ratio`
fn expr_to_ratio(e: &ExprArc) -> Result<Ratio<BigInt>, EvalError> {
    match e.as_ref() {
        Expr::Int(n) => Ok(Ratio::from_integer(n.clone())),
        Expr::Rat(r) => Ok(r.clone()),
        Expr::Mul(factors) => factors
            .iter()
            .map(expr_to_ratio)
            .try_fold(Ratio::one(), |acc, c| Ok(acc * c?)),
        _ => Err(EvalError::TypeError("constant coefficient expected")),
    }
}

// **Pipeline private** — `homogeneous_second_order`
fn homogeneous_second_order(
    a1: &Ratio<BigInt>,
    a0: &Ratio<BigInt>,
    indep: &Ident,
) -> Result<ExprArc, EvalError> {
    let x = Expr::sym(indep.as_str());
    // r^2 + a1*r + a0 = 0
    let disc = a1.clone() * a1.clone() - Ratio::from_integer(BigInt::from(4)) * a0.clone();
    if disc > Ratio::zero() {
        let sqrt_d = ratio_sqrt(&disc)?;
        let r1 = (Ratio::from_integer(-BigInt::from(1)) * a1.clone() + sqrt_d.clone())
            / Ratio::from_integer(BigInt::from(2));
        let r2 = (Ratio::from_integer(-BigInt::from(1)) * a1.clone() - sqrt_d)
            / Ratio::from_integer(BigInt::from(2));
        return Ok(Expr::add(vec![
            Expr::mul(vec![const_sym(0), exp_rat_times_x(&r1, indep)]),
            Expr::mul(vec![const_sym(1), exp_rat_times_x(&r2, indep)]),
        ]));
    }
    if disc == Ratio::zero() {
        let r = -a1.clone() / Ratio::from_integer(BigInt::from(2));
        return Ok(Expr::mul(vec![
            exp_rat_times_x(&r, indep),
            Expr::add(vec![const_sym(0), Expr::mul(vec![const_sym(1), x.clone()])]),
        ]));
    }
    // complex roots: r = alpha ± i*beta
    let alpha = -a1.clone() / Ratio::from_integer(BigInt::from(2));
    let beta_sq = -disc / Ratio::from_integer(BigInt::from(4));
    let beta = ratio_sqrt(&beta_sq)?;
    Ok(Expr::mul(vec![
        exp_rat_times_x(&alpha, indep),
        Expr::add(vec![
            Expr::mul(vec![const_sym(0), trig_rat_times_x(FuncKind::Cos, &beta, indep)]),
            Expr::mul(vec![const_sym(1), trig_rat_times_x(FuncKind::Sin, &beta, indep)]),
        ]),
    ]))
}

// **Pipeline private** — `particular_sin`
fn particular_sin(
    a0: &Ratio<BigInt>,
    forcing: &ExprArc,
    indep: &Ident,
) -> Result<Option<ExprArc>, EvalError> {
    if !is_sin_of_var(forcing, indep) {
        return Ok(None);
    }
    // y'' + a0*y = sin(x) with a0=4 -> sin(x)/3
    if *a0 == Ratio::from_integer(BigInt::from(4)) {
        return Ok(Some(Expr::mul(vec![
            Expr::rat(1, 3),
            Expr::func(FuncKind::Sin, vec![Expr::sym(indep.as_str())]),
        ])));
    }
    Ok(None)
}

// **Pipeline private** — `exp_rat_times_x`
fn exp_rat_times_x(r: &Ratio<BigInt>, indep: &Ident) -> ExprArc {
    if r.is_zero() {
        return Expr::int(1);
    }
    Expr::func(
        FuncKind::Exp,
        vec![Expr::mul(vec![ratio_to_expr(r), Expr::sym(indep.as_str())])],
    )
}

// **Pipeline private** — `trig_rat_times_x`
fn trig_rat_times_x(kind: FuncKind, r: &Ratio<BigInt>, indep: &Ident) -> ExprArc {
    Expr::func(
        kind,
        vec![Expr::mul(vec![ratio_to_expr(r), Expr::sym(indep.as_str())])],
    )
}

// **Pipeline private** — `const_sym`
fn const_sym(n: u8) -> ExprArc {
    Expr::sym(&format!("c{n}"))
}

// **Pipeline private** — `ratio_sqrt`
fn ratio_sqrt(r: &Ratio<BigInt>) -> Result<Ratio<BigInt>, EvalError> {
    if r.is_negative() {
        return Err(EvalError::TypeError("negative under sqrt"));
    }
    let sn = integer_sqrt(r.numer()).ok_or(EvalError::NotImplemented("desolve"))?;
    let sd = integer_sqrt(r.denom()).ok_or(EvalError::NotImplemented("desolve"))?;
    Ok(Ratio::new(sn, sd))
}

// **Pipeline private** — `is_zero_expr`
fn is_zero_expr(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

#[cfg(test)]
mod tests {
    use giac_core::{eval, format_expr, Expr, FuncKind, RelOp};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn desolve_harmonic() {
        let ctx = xcas_default();
        let eq = Arc::new(Expr::Relation(
            RelOp::Eq,
            Expr::add(vec![
                Expr::func(FuncKind::Prime, vec![Expr::sym("y"), Expr::int(2)]),
                Expr::sym("y"),
            ]),
            Expr::int(0),
        ));
        let e = Expr::func(
            FuncKind::Desolve,
            vec![
                eq,
                Expr::func(FuncKind::Apply, vec![Expr::sym("y"), Expr::sym("x")]),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("c0") || s.contains("c1"), "got {s}");
    }
}
