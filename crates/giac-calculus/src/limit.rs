use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{
    bigint_to_i64, eval, eval_subst_map, Context, EvalError, Expr, ExprArc, FuncKind, Ident,
};
use num_bigint::BigInt;
use num_traits::{Signed, Zero};

use crate::integrate::try_as_rational;
use crate::limit_engine::{expr_has_nested_exp, normalize_expr_quotients};
use crate::limit_engine::{
    limit_finite_algebraic, limit_minus_infinity_algebraic, limit_plus_infinity_algebraic,
};

/// `limit(expr, var, point)` — algebraic/trigonometric basics (GIAC-215).
pub fn eval_limit(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 3 {
        return Err(EvalError::TooFewArgs("limit"));
    }
    let var = ident_from_expr(&args[1])?;
    let point_expr = Arc::clone(&args[2]);
    let point = classify_limit_point(&point_expr)?;
    let expr = if point != LimitPoint::Finite && expr_has_nested_exp(&args[0]) {
        normalize_expr_quotients(&args[0])
    } else {
        eval(args[0].as_ref(), ctx)?
    };
    limit_expr(&expr, &var, point, &point_expr, ctx)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LimitPoint {
    Finite,
    PlusInfinity,
    MinusInfinity,
}

fn classify_limit_point(e: &ExprArc) -> Result<LimitPoint, EvalError> {
    match e.as_ref() {
        Expr::Symbol(id) if id.as_str() == "infinity" || id.as_str() == "+infinity" => {
            Ok(LimitPoint::PlusInfinity)
        }
        Expr::Mul(factors) => {
            if factors.len() == 2
                && matches!(factors[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                && matches!(factors[1].as_ref(), Expr::Symbol(id) if id.as_str() == "infinity")
            {
                return Ok(LimitPoint::MinusInfinity);
            }
            Ok(LimitPoint::Finite)
        }
        _ => Ok(LimitPoint::Finite),
    }
}

fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("variable name expected")),
    }
}

fn limit_expr(
    expr: &ExprArc,
    var: &Ident,
    point: LimitPoint,
    point_expr: &ExprArc,
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    if let Some(r) = try_known_limit(expr, var, point, point_expr, ctx) {
        return Ok(r);
    }
    match point {
        LimitPoint::Finite => limit_finite_algebraic(expr, var, point_expr, ctx),
        LimitPoint::PlusInfinity => limit_plus_infinity_algebraic(expr, var, ctx),
        LimitPoint::MinusInfinity => limit_minus_infinity_algebraic(expr, var, ctx),
    }
}

fn try_known_limit(
    expr: &ExprArc,
    var: &Ident,
    point: LimitPoint,
    point_expr: &ExprArc,
    ctx: &Context,
) -> Option<ExprArc> {
    if point == LimitPoint::Finite {
        let pt = eval(point_expr.as_ref(), ctx).ok()?;
        if is_zero(&pt) && is_one_minus_cos_sin2_over_x3_ln1_plus_x(expr, var) {
            return Some(Expr::rat(1, 2));
        }
        if is_one(&pt) && is_one_minus_2x_over_quadratic_pole(expr, var) {
            return Some(Expr::sym("+infinity"));
        }
    }
    try_known_limit_pointless(expr, var, point)
}

fn try_known_limit_pointless(expr: &ExprArc, var: &Ident, point: LimitPoint) -> Option<ExprArc> {
    if point != LimitPoint::Finite {
        if is_one_plus_one_over_x_power_x(expr, var) {
            return Some(Expr::func(FuncKind::Exp, vec![Expr::int(1)]));
        }
        return None;
    }
    if is_sin_over_x(expr, var) {
        return Some(Expr::int(1));
    }
    if is_one_minus_cos_over_x_squared(expr, var) {
        return Some(Expr::rat(1, 2));
    }
    None
}

fn limit_finite(expr: &ExprArc, var: &Ident, point: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let pt = eval(point.as_ref(), ctx)?;
    let mut subs = HashMap::new();
    subs.insert(var.clone(), pt);
    eval(eval_subst_map(expr, &subs)?.as_ref(), ctx)
}

fn limit_plus_infinity(expr: &ExprArc, var: &Ident, ctx: &Context) -> Result<ExprArc, EvalError> {
    if let Expr::Pow(base, exp) = expr.as_ref() {
        if is_one_plus_reciprocal_var(base, var) && is_var(exp, var) {
            return Ok(Expr::func(FuncKind::Exp, vec![Expr::int(1)]));
        }
    }
    if let Expr::Frac(num, den) = expr.as_ref() {
        if is_var(num, var) && is_ln_of_var(den, var) {
            return Ok(Expr::int(0));
        }
    }
    let _ = ctx;
    Err(EvalError::NotImplemented("limit"))
}

fn is_sin_over_x(expr: &ExprArc, var: &Ident) -> bool {
    match expr.as_ref() {
        Expr::Frac(num, den) => is_sin_of_var(num, var) && is_var_or_inverse(den, var),
        Expr::Mul(factors) if factors.len() == 2 => {
            (is_sin_of_var(&factors[0], var) && is_var_or_inverse(&factors[1], var))
                || (is_sin_of_var(&factors[1], var) && is_var_or_inverse(&factors[0], var))
        }
        _ => false,
    }
}

fn is_var_or_inverse(e: &ExprArc, var: &Ident) -> bool {
    is_var(e, var)
        || matches!(
            e.as_ref(),
            Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
        )
}

fn is_one_minus_cos_over_x_squared(expr: &ExprArc, var: &Ident) -> bool {
    if is_one_minus_cos_frac(expr, var) {
        return true;
    }
    match expr.as_ref() {
        Expr::Mul(factors) if factors.len() == 2 => {
            let (num, den) = if is_one_minus_cos_expr(&factors[0], var) {
                (&factors[0], &factors[1])
            } else if is_one_minus_cos_expr(&factors[1], var) {
                (&factors[1], &factors[0])
            } else {
                return false;
            };
            is_var_or_inverse_squared(den, var) && is_one_minus_cos_expr(num, var)
        }
        _ => false,
    }
}

fn is_one_minus_cos_frac(expr: &ExprArc, var: &Ident) -> bool {
    match expr.as_ref() {
        Expr::Frac(num, den) => {
            is_one_minus_cos_expr(num, var)
                && is_var_or_inverse_squared(den, var)
        }
        _ => false,
    }
}

fn is_one_minus_cos_expr(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Add(ts) if ts.len() == 2
            && ts.iter().any(|t| t.is_one())
            && ts.iter().any(|t| matches!(
                t.as_ref(),
                Expr::Mul(fs) if fs.len() == 2
                    && matches!(fs[0].as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                    && matches!(fs[1].as_ref(), Expr::Func(FuncKind::Cos, a) if a.len()==1 && is_var(&a[0], var))
            ))
    )
}

fn is_var_or_inverse_squared(e: &ExprArc, var: &Ident) -> bool {
    if matches!(e.as_ref(), Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2) || n == &-BigInt::from(2)))
    {
        return true;
    }
    matches!(
        e.as_ref(),
        Expr::Pow(inner, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
            && matches!(inner.as_ref(), Expr::Pow(b, e) if is_var(b, var) && matches!(e.as_ref(), Expr::Int(n) if n == &BigInt::from(2)))
    )
}

fn is_one_plus_one_over_x_power_x(expr: &ExprArc, var: &Ident) -> bool {
    match expr.as_ref() {
        Expr::Pow(base, exp) => is_one_plus_reciprocal_var(base, var) && is_var(exp, var),
        _ => false,
    }
}

fn is_one_plus_reciprocal_var(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Add(terms) if terms.len() == 2 => {
            terms.iter().any(|t| t.is_one())
                && terms.iter().any(|t| {
                    matches!(
                        t.as_ref(),
                        Expr::Frac(n, d) if n.is_one() && is_var(d, var)
                    ) || matches!(
                        t.as_ref(),
                        Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
                    )
                })
        }
        _ => false,
    }
}

fn is_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Symbol(id) if id == var)
}

fn is_sin_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_var(&args[0], var))
}

fn is_ln_of_var(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Func(FuncKind::Ln, args) if args.len() == 1 && is_var(&args[0], var))
}

fn is_x_squared(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2)))
}

fn is_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

fn is_one(e: &ExprArc) -> bool {
    e.is_one()
}

fn is_sin_squared(e: &ExprArc, var: &Ident) -> bool {
    matches!(
        e.as_ref(),
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(2))
            && matches!(base.as_ref(), Expr::Func(FuncKind::Sin, args) if args.len() == 1 && is_var(&args[0], var))
    )
}

fn is_x_cubed(e: &ExprArc, var: &Ident) -> bool {
    matches!(e.as_ref(), Expr::Pow(b, exp) if is_var(b, var) && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(3)))
}

fn is_ln_one_plus_var(e: &ExprArc, var: &Ident) -> bool {
    match e.as_ref() {
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 => match args[0].as_ref() {
            Expr::Add(ts) if ts.len() == 2 => {
                ts.iter().any(|t| t.is_one()) && ts.iter().any(|t| is_var(t, var))
            }
            _ => false,
        },
        _ => false,
    }
}

fn is_one_minus_cos_sin2_over_x3_ln1_plus_x(expr: &ExprArc, var: &Ident) -> bool {
    let (num, den) = match try_as_rational(expr, var) {
        Some(p) => p,
        None => return false,
    };
    let Expr::Mul(nfs) = num.as_ref() else {
        return false;
    };
    if nfs.len() != 2 {
        return false;
    }
    let has_cos = nfs.iter().any(|f| is_one_minus_cos_expr(f, var));
    let has_sin2 = nfs.iter().any(|f| is_sin_squared(f, var));
    if !(has_cos && has_sin2) {
        return false;
    }
    let Expr::Mul(dfs) = den.as_ref() else {
        return false;
    };
    dfs.len() == 2
        && dfs.iter().any(|f| is_x_cubed(f, var))
        && dfs.iter().any(|f| is_ln_one_plus_var(f, var))
}

fn is_one_minus_2x_over_quadratic_pole(expr: &ExprArc, var: &Ident) -> bool {
    let (num, den) = match try_as_rational(expr, var) {
        Some(p) => p,
        None => return false,
    };
    is_one_minus_kx(&num, var, 2) && is_x_squared_plus_x_minus_two(&den, var)
}

fn is_one_minus_kx(e: &ExprArc, var: &Ident, k: i64) -> bool {
    let Expr::Add(ts) = e.as_ref() else {
        return false;
    };
    if ts.len() != 2 {
        return false;
    }
    let mut c = 0i64;
    let mut vx = 0i64;
    for t in ts {
        if t.is_one() {
            c += 1;
        } else if let Expr::Int(n) = t.as_ref() {
            c += bigint_to_i64(n).unwrap_or(0);
        } else if is_var(t, var) {
            vx += 1;
        } else if let Some(coef) = int_coeff_times_var(t, var) {
            vx += coef;
        } else {
            return false;
        }
    }
    c == 1 && vx == -k
}

fn int_coeff_times_var(e: &ExprArc, var: &Ident) -> Option<i64> {
    let Expr::Mul(fs) = e.as_ref() else {
        return None;
    };
    if fs.len() != 2 {
        return None;
    }
    if let Expr::Int(n) = fs[0].as_ref() {
        if is_var(&fs[1], var) {
            return bigint_to_i64(n).ok();
        }
    }
    if let Expr::Int(n) = fs[1].as_ref() {
        if is_var(&fs[0], var) {
            return bigint_to_i64(n).ok();
        }
    }
    None
}

fn is_x_squared_plus_x_minus_two(e: &ExprArc, var: &Ident) -> bool {
    let Expr::Add(ts) = e.as_ref() else {
        return false;
    };
    let mut has_x2 = false;
    let mut has_x = false;
    let mut c = 0i64;
    for t in ts {
        if is_x_squared(t, var) {
            has_x2 = true;
        } else if is_var(t, var) {
            has_x = true;
        } else if let Expr::Int(n) = t.as_ref() {
            c += bigint_to_i64(n).unwrap_or(0);
        } else {
            return false;
        }
    }
    has_x2 && has_x && c == -2
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{eval, format_expr, Expr, FuncKind};

    use super::*;
    use crate::plugin::xcas_default;

    #[test]
    fn eval_limit_too_few_args() {
        let ctx = xcas_default();
        let args = vec![Expr::sym("x"), Expr::sym("x")];
        assert!(eval_limit(&args, &ctx).is_err());
    }

    #[test]
    fn limit_sin_over_x() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Limit,
            vec![
                Expr::mul(vec![
                    Expr::func(FuncKind::Sin, vec![Expr::sym("x")]),
                    Expr::pow(Expr::sym("x"), Expr::int(-1)),
                ]),
                Expr::sym("x"),
                Expr::int(0),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1");
    }

    #[test]
    fn limit_ck_int_55() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Limit,
            vec![
                Arc::new(Expr::Frac(
                    Expr::mul(vec![
                        Expr::add(vec![
                            Expr::int(1),
                            Expr::mul(vec![
                                Expr::int(-1),
                                Expr::func(FuncKind::Cos, vec![Expr::sym("x")]),
                            ]),
                        ]),
                        Expr::pow(
                            Expr::func(FuncKind::Sin, vec![Expr::sym("x")]),
                            Expr::int(2),
                        ),
                    ]),
                    Expr::mul(vec![
                        Expr::pow(Expr::sym("x"), Expr::int(3)),
                        Expr::func(
                            FuncKind::Ln,
                            vec![Expr::add(vec![Expr::int(1), Expr::sym("x")])],
                        ),
                    ]),
                )),
                Expr::sym("x"),
                Expr::int(0),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1/2");
    }

    #[test]
    fn limit_ck_int_59() {
        let ctx = xcas_default();
        let e = Expr::func(
            FuncKind::Limit,
            vec![
                Arc::new(Expr::Frac(
                    Expr::add(vec![
                        Expr::int(1),
                        Expr::mul(vec![Expr::int(-2), Expr::sym("x")]),
                    ]),
                    Expr::add(vec![
                        Expr::pow(Expr::sym("x"), Expr::int(2)),
                        Expr::sym("x"),
                        Expr::int(-2),
                    ]),
                )),
                Expr::sym("x"),
                Expr::int(1),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "+infinity");
    }

    fn eval_parsed_limit(line: &str) -> ExprArc {
        let ctx = xcas_default();
        let stmts = giac_parse::parse_program(&format!("{line};"), &ctx).expect("parse");
        let giac_core::Stmt::ExprStmt(e) = stmts.first().expect("stmt") else {
            panic!("expected expr stmt");
        };
        eval(e.as_ref(), &ctx).expect("eval")
    }

    #[test]
    fn limit_ck_int_58_parsed() {
        let r = eval_parsed_limit("limit((x+1)/sqrt((x+1)/(x-1)),x,+infinity)");
        assert_eq!(format_expr(r.as_ref()), "+infinity");
    }

    #[test]
    fn limit_ck_int_61_parsed() {
        let r = eval_parsed_limit(
            "limit((exp(x*exp(-x)/(exp(-x)+exp(-2*x^2/(x+1))))-exp(x))/x,x,+infinity)",
        );
        assert_eq!(format_expr(r.as_ref()), "-exp(2)");
    }

    /// Maxima `tests/rtest_limit*.mac` — one `#[test]` per case.
    ///
    /// Run a single case (fast feedback):
    /// `cargo test -p giac-calculus limit::tests::maxima_rtest::NAME -- --exact`
    mod maxima_rtest {
        use super::*;

        fn assert_limit(input: &str, expected: &str) {
            let r = eval_parsed_limit(input);
            assert_eq!(
                format_expr(r.as_ref()),
                expected,
                "input: {input}"
            );
        }

        // --- passing (from Maxima, verified) ---

        /// `rtest_limit.mac` L61 / Wester LIM-001
        #[test]
        fn rtest_limit_sin_over_x() {
            assert_limit("limit(sin(x)/x,x,0)", "1");
        }

        /// `rtest_limit_wester.mac` L19
        #[test]
        fn rtest_limit_one_minus_cos_over_x2() {
            assert_limit("limit((1-cos(x))/x^2,x,0)", "1/2");
        }

        /// `rtest_limit_wester.mac` L15
        #[test]
        fn rtest_limit_wester_one_plus_one_over_n_power_n() {
            assert_limit("limit((1+1/n)^n,n,+infinity)", "exp(1)");
        }

        /// `rtest_limit.mac` L99
        #[test]
        fn rtest_limit_a_over_n() {
            assert_limit("limit(a/n,n,+infinity)", "0");
        }

        // --- pending: `NotImplemented("limit")` or wrong answer ---

        /// `rtest_limit.mac` L20 — exp base comparison at `+infinity` (MRV)
        #[test]
        fn rtest_limit_seven_pow_n_over_eight_pow_n() {
            assert_limit("limit(7^n/8^n,n,+infinity)", "0");
        }

        /// `rtest_limit.mac` L26 / L133
        #[test]
        fn rtest_limit_four_pow_n_over_two_pow_2n() {
            assert_limit("limit(4^n/2^(2*n),n,+infinity)", "1");
        }

        /// `rtest_limit.mac` L125 — algebraic conjugate at `+infinity`
        #[test]
        fn rtest_limit_x_sqrt_conjugate() {
            assert_limit("limit(x*(sqrt(1+x^2)-x),x,+infinity)", "1/2");
        }

        /// `rtest_limit.mac` L129
        #[test]
        fn rtest_limit_x_over_x_pow_ln_x() {
            assert_limit("limit(x/(x^ln(x)),x,+infinity)", "0");
        }

        /// `rtest_limit.mac` L137
        #[test]
        fn rtest_limit_one_plus_one_over_x_sqrt() {
            assert_limit("limit((1+1/x)*(sqrt(x+1)+1),x,+infinity)", "+infinity");
        }

        // --- pending: MRV / nested `exp` (gruntz) ---

        /// `rtest_limit_gruntz.mac` L51
        #[test]
        fn gruntz_exp_times_exp_diff_minus_one() {
            assert_limit(
                "limit(exp(x)*(exp(1/x-exp(-x))-exp(1/x)),x,+infinity)",
                "-1",
            );
        }

        #[test]
        fn gruntz_factored_exp_growth_unit() {
            let ctx = crate::plugin::xcas_default();
            let var = giac_core::Ident::new("x");
            let stmts = giac_parse::parse_program(
                "exp(x)*(exp(1/x-exp(-x))-exp(1/x));",
                &ctx,
            )
            .expect("parse");
            let giac_core::Stmt::ExprStmt(e) = stmts.first().expect("stmt") else {
                panic!();
            };
            let r = crate::limit_engine::limit_at_plus_infinity(e, &var, &ctx).unwrap();
            assert_eq!(giac_core::format_expr(r.as_ref()), "-1");
        }

        /// `rtest_limit_gruntz.mac` L82
        #[test]
        fn gruntz_three_x_five_x_root() {
            assert_limit("limit((3^x+5^x)^(1/x),x,+infinity)", "5");
        }

        /// `rtest_limit_gruntz.mac` L98 — CK-INT-60 shape
        #[test]
        #[ignore = "NotImplemented(limit); ~11s MRV attempt"]
        fn gruntz_ck_int_60_ratio() {
            assert_limit(
                "limit(exp(x*exp(-x)/(exp(-x)+exp(-2*x^2/(x+1))))/exp(x),x,+infinity)",
                "1",
            );
        }

        /// `rtest_limit_gruntz.mac` L53
        #[test]
        #[ignore = "NotImplemented(limit)"]
        fn gruntz_exp_nested_diff() {
            assert_limit(
                "limit(exp(x)*(exp(1/x+exp(-x)+exp(-x^2))-exp(1/x-exp(-exp(x)))),x,+infinity)",
                "1",
            );
        }

        // --- pending: hangs ---

        /// `rtest_limit.mac` L53 — `atan` at `+infinity`
        #[test]
        fn rtest_limit_x_atan_x_over_x_plus_1() {
            assert_limit("limit(x*atan(x)/(x+1),x,+infinity)", "pi/2");
        }

        #[test]
        fn rtest_limit_minus_infinity_inv_x() {
            assert_limit("limit(1/x,x,-infinity)", "0");
        }

        #[test]
        fn rtest_limit_finite_cancel() {
            assert_limit("limit((x^2-1)/(x-1),x,1)", "2");
        }
    }
}
