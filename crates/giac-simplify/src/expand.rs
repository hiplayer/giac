//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-simplify-api-stability.md`.
//!
//!
use std::sync::Arc;

use giac_poly::modp;
use num_bigint::BigInt;
use num_traits::Zero;

use giac_core::{bigint_to_i64, bigint_to_nonneg_u32, Context, EvalError, Expr, ExprArc, FuncKind};

use giac_core::{expr_to_poly, poly_mod_to_expr, poly_to_expr};

/// Max exponent for generic symbolic power expansion (binomial / repeated multiply).
pub(crate) const MAX_EXPAND_POWER: u32 = 20;
/// Max exponent for modular polynomial power expansion.
pub(crate) const MAX_MOD_EXPAND_POWER: u32 = 50;

/// Controls how [`expand_with_policy`] distributes products and expands powers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ExpandPolicy {
    /// Distribute over all `Add`, expand binomial powers of sums.
    #[default]
    Full,
    /// Skip `Mul`×`Add` distribution and binomial expansion when any involved
    /// subtree contains `exp`/`ln` (preserves MRV / `exp(·)-1` shapes).
    NoExpDistribute,
}

/// **Stable** — distribute products over sums and expand powers of sums (`ExpandPolicy::Full`).
pub fn expand(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    expand_with_policy(expr, ctx, ExpandPolicy::Full)
}

/// **Stable** — [`expand`] with an explicit policy.
pub fn expand_with_policy(
    expr: &Expr,
    ctx: &Context,
    policy: ExpandPolicy,
) -> Result<ExprArc, EvalError> {
    match expr {
        Expr::Add(terms) => {
            let expanded: Result<Vec<_>, _> = terms
                .iter()
                .map(|t| expand_with_policy(t.as_ref(), ctx, policy))
                .collect();
            Ok(Expr::add(expanded?))
        }
        Expr::Mul(factors) => {
            let mut acc = expand_with_policy(&factors[0], ctx, policy)?;
            for f in &factors[1..] {
                acc = expand_mul_pair(
                    acc.as_ref(),
                    expand_with_policy(f.as_ref(), ctx, policy)?.as_ref(),
                    ctx,
                    policy,
                )?;
            }
            Ok(acc)
        }
        Expr::Pow(base, exp) => expand_pow(base, exp, ctx, policy),
        Expr::Frac(num, den) => Ok(Arc::new(Expr::Frac(
            expand_with_policy(num.as_ref(), ctx, policy)?,
            expand_with_policy(den.as_ref(), ctx, policy)?,
        ))),
        Expr::Mod(base, modulus) => Ok(Arc::new(Expr::Mod(
            expand_with_policy(base.as_ref(), ctx, policy)?,
            expand_with_policy(modulus.as_ref(), ctx, policy)?,
        ))),
        Expr::Complex(re, im) => Ok(Arc::new(Expr::Complex(
            expand_with_policy(re.as_ref(), ctx, policy)?,
            expand_with_policy(im.as_ref(), ctx, policy)?,
        ))),
        Expr::Func(_, _) | Expr::Symbol(_) | Expr::Int(_) | Expr::Rat(_) => {
            Ok(Arc::new(expr.clone()))
        }
        other => Ok(Arc::new(other.clone())),
    }
}

/// **Stable** — polynomial-oriented expand: no distribution through `exp`/`ln` subtrees.
pub fn expand_polynomial(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    expand_with_policy(expr, ctx, ExpandPolicy::NoExpDistribute)
}

/// True when `e` or any descendant is `exp(...)` or `ln(...)`.
///
/// **Pipeline private** — used by `ExpandPolicy::NoExpDistribute`.
pub(crate) fn expr_contains_exp_ln(e: &Expr) -> bool {
    match e {
        Expr::Func(FuncKind::Exp | FuncKind::Ln, _) => true,
        Expr::Add(terms) => terms.iter().any(|t| expr_contains_exp_ln(t.as_ref())),
        Expr::Mul(factors) => factors.iter().any(|f| expr_contains_exp_ln(f.as_ref())),
        Expr::Pow(base, exp) => {
            expr_contains_exp_ln(base.as_ref()) || expr_contains_exp_ln(exp.as_ref())
        }
        Expr::Frac(num, den) => {
            expr_contains_exp_ln(num.as_ref()) || expr_contains_exp_ln(den.as_ref())
        }
        Expr::Func(_, args) => args.iter().any(|a| expr_contains_exp_ln(a.as_ref())),
        Expr::Mod(base, modulus) => {
            expr_contains_exp_ln(base.as_ref()) || expr_contains_exp_ln(modulus.as_ref())
        }
        Expr::Complex(re, im) => {
            expr_contains_exp_ln(re.as_ref()) || expr_contains_exp_ln(im.as_ref())
        }
        _ => false,
    }
}

// **Pipeline private** — distribute one mul factor over add
fn expand_mul_pair(
    lhs: &Expr,
    rhs: &Expr,
    ctx: &Context,
    policy: ExpandPolicy,
) -> Result<ExprArc, EvalError> {
    match (lhs, rhs) {
        (Expr::Add(terms), other) | (other, Expr::Add(terms)) => {
            if policy == ExpandPolicy::NoExpDistribute {
                let touches_exp = expr_contains_exp_ln(other)
                    || terms.iter().any(|t| expr_contains_exp_ln(t.as_ref()));
                if touches_exp {
                    return Ok(Expr::mul(vec![
                        Arc::new(lhs.clone()),
                        Arc::new(rhs.clone()),
                    ]));
                }
            }
            let (terms, other) = if matches!(lhs, Expr::Add(_)) {
                (terms, other)
            } else {
                match rhs {
                    Expr::Add(t) => (t, lhs),
                    _ => unreachable!(),
                }
            };
            let parts: Result<Vec<_>, _> = terms
                .iter()
                .map(|t| expand_mul_pair(t.as_ref(), other, ctx, policy))
                .collect();
            Ok(Expr::add(parts?))
        }
        _ => Ok(Expr::mul(vec![Arc::new(lhs.clone()), Arc::new(rhs.clone())])),
    }
}

// **Pipeline private** — expand integer powers and mod-poly powers
fn expand_pow(
    base: &ExprArc,
    exp: &ExprArc,
    ctx: &Context,
    policy: ExpandPolicy,
) -> Result<ExprArc, EvalError> {
    if let Expr::Mod(b, m) = base.as_ref() {
        if let (Ok(p), Ok(mod_i)) = (expr_to_poly(b), modulus_from_expr(m)) {
            if let Expr::Int(n) = exp.as_ref() {
                if n >= &BigInt::zero() && bigint_to_nonneg_u32(n).ok() <= Some(MAX_MOD_EXPAND_POWER)
                {
                    let e = bigint_to_nonneg_u32(n)?;
                    let pm = modp(&p.pow(u64::from(e)), mod_i)?;
                    return Ok(poly_mod_to_expr(&pm));
                }
            }
        }
    }
    let base_e = expand_with_policy(base.as_ref(), ctx, policy)?;
    if let Expr::Int(n) = exp.as_ref() {
        if n >= &BigInt::zero() {
            let e = bigint_to_nonneg_u32(n)?;
            if e > MAX_EXPAND_POWER {
                return Ok(Expr::pow(base_e, Arc::clone(exp)));
            }
            if e == 0 {
                return Ok(Expr::int(1));
            }
            if let Expr::Add(terms) = base_e.as_ref() {
                if policy == ExpandPolicy::NoExpDistribute
                    && terms.iter().any(|t| expr_contains_exp_ln(t.as_ref()))
                {
                    return Ok(Expr::pow(base_e, Arc::clone(exp)));
                }
                return Ok(expand_binomial(terms, e));
            }
            if e == 1 {
                return Ok(base_e);
            }
            if matches!(
                base_e.as_ref(),
                Expr::Func(
                    FuncKind::Sin
                        | FuncKind::Cos
                        | FuncKind::Exp
                        | FuncKind::Ln
                        | FuncKind::Atan
                        | FuncKind::Tan,
                    _,
                )
            ) {
                return Ok(Expr::pow(base_e, Arc::clone(exp)));
            }
            return Ok(repeated_mul(base_e, e, ctx, policy)?);
        }
    }
    Ok(Expr::pow(base_e, Arc::clone(exp)))
}

// **Pipeline private** — repeated multiply for small integer power
fn repeated_mul(
    base: ExprArc,
    exp: u32,
    ctx: &Context,
    policy: ExpandPolicy,
) -> Result<ExprArc, EvalError> {
    let mut result = Arc::clone(&base);
    for _ in 1..exp {
        result = expand_mul_pair(result.as_ref(), base.as_ref(), ctx, policy)?;
    }
    Ok(result)
}

// **Pipeline private** — binomial power via poly or Expr::pow
fn expand_binomial(terms: &[ExprArc], n: u32) -> ExprArc {
    let sum = Expr::add(terms.to_vec());
    if let Ok(p) = expr_to_poly(sum.as_ref()) {
        return poly_to_expr(&p.pow(u64::from(n)));
    }
    Expr::pow(sum, Expr::int(n as i64))
}

/// **Stable** — expand then collect into polynomial normal form.
///
/// `AlgExt` subexpressions are passed through `eval`; transcendental leaves are
/// unchanged when `expr_to_poly` fails.
pub fn normal(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    if giac_core::contains_algext(expr) {
        let expanded = expand(expr, ctx)?;
        return giac_core::eval(expanded.as_ref(), ctx);
    }
    if let Ok(p) = expr_to_poly(expr) {
        return Ok(poly_to_expr(&p));
    }
    let expanded = expand(expr, ctx)?;
    if let Ok(p) = expr_to_poly(expanded.as_ref()) {
        return Ok(poly_to_expr(&p));
    }
    Ok(expanded)
}

// **Pipeline private** — coerce Expr modulus to i64
fn modulus_from_expr(m: &Expr) -> Result<i64, EvalError> {
    match m {
        Expr::Int(n) => bigint_to_i64(n),
        _ => Err(EvalError::TypeError("integer modulus expected")),
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::{format_expr, Context, Expr};

    #[test]
    fn normal_mod_power_displays_giac_style() {
        let ctx = crate::plugin::xcas_default();
        let e = Expr::func(
            giac_core::FuncKind::Normal,
            vec![Expr::pow(
                Arc::new(Expr::Mod(
                    Expr::add(vec![
                        Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
                        Expr::int(1),
                    ]),
                    Expr::int(13),
                )),
                Expr::int(5),
            )],
        );
        let r = giac_core::eval(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(!s.contains(" mod 13*"), "nested mod display: {s}");
        assert_eq!(
            s,
            "(6 % 13)*x^5+(2 % 13)*x^4+(2 % 13)*x^3+(1 % 13)*x^2+(10 % 13)*x+(1 % 13)"
        );
    }

    #[test]
    fn expand_square_of_sum() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(3)]), Expr::int(4));
        let r = normal(e.as_ref(), &ctx).unwrap();
        assert_eq!(
            format_expr(r.as_ref()),
            "x^4+12*x^3+54*x^2+108*x+81"
        );
    }

    #[test]
    fn normal_cancels_distributed_polynomial_terms() {
        let ctx = Context::default();
        let e = Expr::add(vec![
            Expr::mul(vec![
                Expr::int(2),
                Expr::sym("x"),
                Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            ]),
            Expr::mul(vec![Expr::int(-2), Expr::pow(Expr::sym("x"), Expr::int(2))]),
        ]);
        let r = normal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2*x");
    }

    #[test]
    fn expand_distribute_mul_over_add() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::sym("y"),
        ]);
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x*y+1*y");
    }

    #[test]
    fn expand_polynomial_distribute_rational_product() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(-1)]), Expr::int(-1)),
        ]);
        let r = expand_polynomial(e.as_ref(), &ctx).unwrap();
        assert!(
            matches!(r.as_ref(), Expr::Add(terms) if terms.len() == 2),
            "expected distributed sum, got {}",
            format_expr(r.as_ref())
        );
    }

    #[test]
    fn expand_polynomial_skips_exp_subtree() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::add(vec![
                Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
                Expr::int(1),
            ]),
            Expr::sym("y"),
        ]);
        let r = expand_polynomial(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(exp(x)+1)*y");
    }

    #[test]
    fn expand_polynomial_skips_when_other_factor_has_exp() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::func(FuncKind::Exp, vec![Expr::sym("x")]),
        ]);
        let r = expand_polynomial(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(x+1)*exp(x)");
    }

    #[test]
    fn expand_pow_single_base() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::sym("x"), Expr::int(1));
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(r, Expr::sym("x"));
    }

    #[test]
    fn expand_pow_cube_of_symbol() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::sym("x"), Expr::int(3));
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(x*x)*x");
    }

    #[test]
    fn expand_complex_and_non_poly_normal() {
        let ctx = Context::default();
        let c = Expr::Complex(Expr::sym("a"), Expr::sym("b"));
        let r = expand(&c, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "a+b*i");

        let trig = Expr::func(giac_core::FuncKind::Sin, vec![Expr::sym("x")]);
        let n = normal(trig.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(n.as_ref()), "sin(x)");
    }

    #[test]
    fn expand_rhs_add_and_zero_power() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::sym("y"),
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
        ]);
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x*y+1*y");

        let z = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(0));
        assert_eq!(expand(z.as_ref(), &ctx).unwrap(), Expr::int(1));
    }

    #[test]
    fn normal_fast_path_already_polynomial() {
        let ctx = Context::default();
        let e = Expr::add(vec![
            Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
            Expr::int(3),
        ]);
        let r = normal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2*x+3");
    }

    #[test]
    // smoke-until B-EXPAND-BINOM: delete when `expand_binomial_fallback_semantic` green
    fn expand_binomial_fallback_for_non_poly() {
        let ctx = Context::default();
        let e = Expr::pow(
            Expr::add(vec![Expr::func(giac_core::FuncKind::Sin, vec![Expr::sym("x")]), Expr::int(1)]),
            Expr::int(2),
        );
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(sin(x)+1)^2");
    }

    #[test]
    #[ignore = "B-EXPAND-BINOM: expand 应收敛 (sin(x)+1)^2 → sin(x)^2+2*sin(x)+1"]
    fn expand_binomial_fallback_semantic() {
        let ctx = Context::default();
        let e = Expr::pow(
            Expr::add(vec![Expr::func(giac_core::FuncKind::Sin, vec![Expr::sym("x")]), Expr::int(1)]),
            Expr::int(2),
        );
        let r = expand(e.as_ref(), &ctx).unwrap();
        let expected = Expr::add(vec![
            Expr::pow(Expr::func(giac_core::FuncKind::Sin, vec![Expr::sym("x")]), Expr::int(2)),
            Expr::mul(vec![Expr::int(2), Expr::func(giac_core::FuncKind::Sin, vec![Expr::sym("x")])]),
            Expr::int(1),
        ]);
        assert!(
            crate::assert_equiv(r.as_ref(), &expected, &ctx).unwrap(),
            "got {}",
            format_expr(r.as_ref())
        );
    }

    #[test]
    fn expand_recurses_into_frac() {
        let ctx = Context::default();
        let e = Expr::Frac(
            Expr::mul(vec![
                Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
                Expr::sym("y"),
            ]),
            Expr::sym("z"),
        );
        let r = expand(&e, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(x*y+1*y)/(z)");
    }
}
