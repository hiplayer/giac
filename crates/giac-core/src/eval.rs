use std::sync::Arc;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::context::Context;
use crate::error::EvalError;
use crate::expr::{Expr, ExprArc, FuncKind};
use crate::ident::Ident;
use crate::simplify::simplify;

/// Evaluate an expression in the given context.
pub fn eval(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> {
    let result = match expr {
        Expr::Int(_) | Expr::Rat(_) | Expr::Str(_) | Expr::Undefined => Arc::new(expr.clone()),
        Expr::Symbol(id) => eval_symbol(id, ctx)?,
        Expr::Add(terms) => eval_add(terms, ctx)?,
        Expr::Mul(factors) => eval_mul(factors, ctx)?,
        Expr::Pow(base, exp) => eval_pow(base, exp, ctx)?,
        Expr::Func(kind, args) => eval_func(*kind, args, ctx)?,
        Expr::Frac(num, den) => eval_frac(num, den, ctx)?,
        Expr::Complex(re, im) => {
            let re = eval(re, ctx)?;
            let im = eval(im, ctx)?;
            let c_re = try_as_complex(re.as_ref(), ctx)
                .map(|c| c.re)
                .unwrap_or(Ratio::zero());
            let c_im = try_as_complex(im.as_ref(), ctx)
                .map(|c| c.re)
                .unwrap_or(Ratio::zero());
            complex_to_expr(c_re, c_im)?
        }
        Expr::Seq(items) => {
            let ev: Result<Vec<_>, _> = items.iter().map(|e| eval(e, ctx)).collect();
            Arc::new(Expr::Seq(ev?))
        }
        Expr::List(items) => {
            let ev: Result<Vec<_>, _> = items.iter().map(|e| eval(e, ctx)).collect();
            Arc::new(Expr::List(ev?))
        }
        other => Arc::new(other.clone()),
    };
    simplify(result.as_ref(), ctx)
}

fn eval_symbol(id: &Ident, ctx: &Context) -> Result<ExprArc, EvalError> {
    if let Some(val) = ctx.get(id) {
        return Ok(Arc::clone(val));
    }
    if ctx.complex_mode && id.is_imaginary_unit() {
        return Ok(Expr::sym("i"));
    }
    Ok(Expr::sym(id.as_str()))
}

fn eval_add(terms: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    let mut complex_sum = ComplexVal::zero();
    let mut symbolic = Vec::new();
    let mut all_numeric = true;

    for t in terms {
        let ev = eval(t, ctx)?;
        if let Some(c) = try_as_complex(ev.as_ref(), ctx) {
            complex_sum = complex_sum + c;
        } else {
            all_numeric = false;
            symbolic.push(ev);
        }
    }

    if all_numeric && !complex_sum.is_zero() {
        return complex_to_expr(complex_sum.re, complex_sum.im);
    }

    if !complex_sum.is_zero() {
        symbolic.push(complex_to_expr(complex_sum.re, complex_sum.im)?);
    }

    match symbolic.len() {
        0 => Ok(Expr::int(0)),
        1 => Ok(Arc::clone(&symbolic[0])),
        _ => Ok(Expr::add(symbolic)),
    }
}

fn eval_mul(factors: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    let mut complex_prod = ComplexVal::one();
    let mut symbolic = Vec::new();
    let mut all_numeric = true;

    for f in factors {
        let ev = eval(f, ctx)?;
        if let Some(c) = try_as_complex(ev.as_ref(), ctx) {
            complex_prod = complex_prod * c;
        } else {
            all_numeric = false;
            symbolic.push(ev);
        }
    }

    if all_numeric {
        return complex_to_expr(complex_prod.re, complex_prod.im);
    }

    if complex_prod != ComplexVal::one() {
        symbolic.insert(0, complex_to_expr(complex_prod.re, complex_prod.im)?);
    }

    match symbolic.len() {
        0 => Ok(Expr::int(1)),
        1 => Ok(Arc::clone(&symbolic[0])),
        _ => Ok(Expr::mul(symbolic)),
    }
}

fn eval_pow(base: &ExprArc, exp: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let base = eval(base, ctx)?;
    let exp = eval(exp, ctx)?;

    if let (Some(b), Some(e)) = (try_as_complex(base.as_ref(), ctx), as_nonneg_int(exp.as_ref())) {
        if e == 0 {
            return Ok(Expr::int(1));
        }
        let mut result = ComplexVal::one();
        let mut base_c = b;
        let mut n = e;
        while n > 0 {
            if n % 2 == 1 {
                result = result * base_c.clone();
            }
            base_c = base_c.clone() * base_c;
            n /= 2;
        }
        return complex_to_expr(result.re, result.im);
    }

    if let (Expr::Int(b), Expr::Int(e)) = (base.as_ref(), exp.as_ref()) {
        if e >= &BigInt::zero() && e <= &BigInt::from(30) {
            let e_u = e.to_string().parse::<u32>().unwrap();
            return Ok(Expr::int(int_to_i64(&(b.pow(e_u)))?));
        }
    }

    Ok(Expr::pow(base, exp))
}

fn eval_frac(num: &ExprArc, den: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let num = eval(num, ctx)?;
    let den = eval(den, ctx)?;
    if let (Some(n), Some(d)) = (as_int(num.as_ref()), as_int(den.as_ref())) {
        if d.is_zero() {
            return Err(EvalError::DivisionByZero);
        }
        return Ok(Expr::rat(
            int_to_i64(n)?,
            int_to_i64(d)?,
        ));
    }
    Ok(Arc::new(Expr::Frac(num, den)))
}

fn eval_func(kind: FuncKind, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    let ev: Result<Vec<_>, _> = args.iter().map(|a| eval(a, ctx)).collect();
    let args = ev?;
    match kind {
        FuncKind::Abs => eval_abs(&args),
        FuncKind::Gcd => eval_gcd(&args),
        FuncKind::Conj => eval_conj(&args, ctx),
        FuncKind::Sqrt => eval_sqrt(&args),
        FuncKind::Re => eval_component(&args, true, ctx),
        FuncKind::Im => eval_component(&args, false, ctx),
        other => Err(EvalError::NotImplemented(func_name(other))),
    }
}

fn eval_abs(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(if args.is_empty() {
            EvalError::TooFewArgs("abs")
        } else {
            EvalError::TooManyArgs("abs")
        });
    }
    if let Some(c) = try_as_complex(args[0].as_ref(), &Context::new()) {
        let norm_sq = &c.re * &c.re + &c.im * &c.im;
        if norm_sq.is_zero() {
            return Ok(Expr::int(0));
        }
        if norm_sq.is_one() {
            return Ok(Expr::int(1));
        }
        if norm_sq.denom() == &BigInt::one() {
            let n = norm_sq.numer();
            if is_perfect_square(n) {
                return Ok(Expr::int(int_to_i64(&integer_sqrt(n))?));
            }
            return Ok(Expr::func(FuncKind::Sqrt, vec![Expr::int(int_to_i64(n)?)]));
        }
    }
    Ok(Expr::func(FuncKind::Abs, args.to_vec()))
}

fn eval_gcd(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() < 2 {
        return Err(EvalError::TooFewArgs("gcd"));
    }
    let mut result: Option<BigInt> = None;
    for a in args {
        let n = as_int(a.as_ref()).ok_or(EvalError::TypeError("gcd expects integers"))?;
        result = Some(match result {
            None => n.clone(),
            Some(r) => r.gcd(n),
        });
    }
    Ok(Expr::int(int_to_i64(result.as_ref().unwrap())?))
}

fn eval_conj(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(if args.is_empty() {
            EvalError::TooFewArgs("conj")
        } else {
            EvalError::TooManyArgs("conj")
        });
    }
    if let Some(c) = try_as_complex(args[0].as_ref(), ctx) {
        return complex_to_expr(c.re, -c.im);
    }
    Ok(Expr::func(FuncKind::Conj, args.to_vec()))
}

fn eval_sqrt(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(if args.is_empty() {
            EvalError::TooFewArgs("sqrt")
        } else {
            EvalError::TooManyArgs("sqrt")
        });
    }
    if let Some(n) = as_int(args[0].as_ref()) {
        if n.is_negative() {
            return Err(EvalError::TypeError("sqrt of negative integer"));
        }
        if is_perfect_square(n) {
            return Ok(Expr::int(int_to_i64(&integer_sqrt(n))?));
        }
    }
    Ok(Expr::func(FuncKind::Sqrt, args.to_vec()))
}

fn eval_component(args: &[ExprArc], real_part: bool, ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs(if real_part { "re" } else { "im" }));
    }
    if let Some(c) = try_as_complex(args[0].as_ref(), ctx) {
        let v = if real_part { c.re } else { c.im };
        return ratio_to_expr(&v);
    }
    if real_part {
        Ok(Arc::clone(&args[0]))
    } else {
        Ok(Expr::int(0))
    }
}

#[derive(Clone, Debug, PartialEq)]
struct ComplexVal {
    re: Ratio<BigInt>,
    im: Ratio<BigInt>,
}

impl ComplexVal {
    fn zero() -> Self {
        Self {
            re: Ratio::zero(),
            im: Ratio::zero(),
        }
    }

    fn one() -> Self {
        Self {
            re: Ratio::one(),
            im: Ratio::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }
}

impl std::ops::Add for ComplexVal {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl std::ops::Mul for ComplexVal {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: &self.re * &rhs.re - &self.im * &rhs.im,
            im: &self.re * &rhs.im + &self.im * &rhs.re,
        }
    }
}

fn try_as_complex(expr: &Expr, ctx: &Context) -> Option<ComplexVal> {
    match expr {
        Expr::Int(n) => Some(ComplexVal {
            re: Ratio::from_integer(n.clone()),
            im: Ratio::zero(),
        }),
        Expr::Rat(r) => Some(ComplexVal {
            re: r.clone(),
            im: Ratio::zero(),
        }),
        Expr::Symbol(id) if ctx.complex_mode && id.is_imaginary_unit() => Some(ComplexVal {
            re: Ratio::zero(),
            im: Ratio::one(),
        }),
        Expr::Add(terms) => {
            let mut sum = ComplexVal::zero();
            for t in terms {
                sum = sum + try_as_complex(t, ctx)?;
            }
            Some(sum)
        }
        Expr::Mul(factors) => {
            let mut prod = ComplexVal::one();
            for f in factors {
                prod = prod * try_as_complex(f, ctx)?;
            }
            Some(prod)
        }
        Expr::Pow(base, exp) => {
            let b = try_as_complex(base, ctx)?;
            let e = as_nonneg_int(exp)?;
            if e == 0 {
                return Some(ComplexVal::one());
            }
            let mut result = ComplexVal::one();
            let mut base_c = b;
            let mut n = e;
            while n > 0 {
                if n % 2 == 1 {
                    result = result * base_c.clone();
                }
                base_c = base_c.clone() * base_c.clone();
                n /= 2;
            }
            Some(result)
        }
        _ => None,
    }
}

fn complex_to_expr(re: Ratio<BigInt>, im: Ratio<BigInt>) -> Result<ExprArc, EvalError> {
    if im.is_zero() {
        return ratio_to_expr(&re);
    }
    if re.is_zero() {
        return ratio_to_expr(&im).map(|im_e| {
            if matches!(im_e.as_ref(), Expr::Int(n) if n.is_one()) {
                Expr::sym("i")
            } else if matches!(im_e.as_ref(), Expr::Int(n) if n == &-BigInt::one()) {
                Expr::mul(vec![Expr::int(-1), Expr::sym("i")])
            } else {
                Expr::mul(vec![im_e, Expr::sym("i")])
            }
        });
    }

    let mut terms = vec![ratio_to_expr(&re)?];
    let im_term = if im == Ratio::one() {
        Expr::sym("i")
    } else if im == Ratio::from_integer(-BigInt::one()) {
        Expr::mul(vec![Expr::int(-1), Expr::sym("i")])
    } else {
        Expr::mul(vec![ratio_to_expr(&im)?, Expr::sym("i")])
    };
    terms.push(im_term);
    Ok(Expr::add(terms))
}

fn ratio_to_expr(r: &Ratio<BigInt>) -> Result<ExprArc, EvalError> {
    if r.is_zero() {
        Ok(Expr::int(0))
    } else if r.denom() == &BigInt::one() {
        Ok(Expr::int(int_to_i64(r.numer())?))
    } else {
        Ok(Arc::new(Expr::Rat(r.clone())))
    }
}

fn as_int(expr: &Expr) -> Option<&BigInt> {
    match expr {
        Expr::Int(n) => Some(n),
        _ => None,
    }
}

fn as_nonneg_int(expr: &Expr) -> Option<u32> {
    match expr {
        Expr::Int(n) if n >= &BigInt::zero() => n.to_string().parse().ok(),
        _ => None,
    }
}

fn int_to_i64(n: &BigInt) -> Result<i64, EvalError> {
    n.to_string()
        .parse()
        .map_err(|_| EvalError::TypeError("integer out of i64 range"))
}

fn is_perfect_square(n: &BigInt) -> bool {
    if n.is_negative() {
        return false;
    }
    let root = integer_sqrt(n);
    &root * &root == *n
}

fn integer_sqrt(n: &BigInt) -> BigInt {
    n.sqrt()
}

fn func_name(kind: FuncKind) -> &'static str {
    match kind {
        FuncKind::Abs => "abs",
        FuncKind::Gcd => "gcd",
        FuncKind::Conj => "conj",
        FuncKind::Sqrt => "sqrt",
        FuncKind::Sin => "sin",
        FuncKind::Cos => "cos",
        FuncKind::Exp => "exp",
        FuncKind::Ln => "ln",
        FuncKind::Re => "re",
        FuncKind::Im => "im",
        FuncKind::Arg => "arg",
        FuncKind::Sign => "sign",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::display::format_expr;
    use crate::expr::Expr;

    fn ctx() -> Context {
        Context::xcas_default()
    }

    #[test]
    fn eval_integer_add_mul() {
        let e = Expr::add(vec![Expr::int(1), Expr::mul(vec![Expr::int(2), Expr::int(3)])]);
        let r = eval(e.as_ref(), &ctx()).unwrap();
        assert_eq!(r, Expr::int(7));
    }

    #[test]
    fn eval_gcd() {
        let e = Expr::func(FuncKind::Gcd, vec![Expr::int(45), Expr::int(75)]);
        let r = eval(e.as_ref(), &ctx()).unwrap();
        assert_eq!(format_expr(r.as_ref()), "15");
    }

    #[test]
    fn eval_abs_complex() {
        let i = Expr::sym("i");
        let e = Expr::func(
            FuncKind::Abs,
            vec![Expr::add(vec![Expr::int(1), Expr::mul(vec![Expr::int(2), i])])],
        );
        let r = eval(e.as_ref(), &ctx()).unwrap();
        assert_eq!(format_expr(r.as_ref()), "sqrt(5)");
    }

    #[test]
    fn eval_conj_power() {
        let i = Expr::sym("i");
        let inner = Expr::pow(
            Expr::add(vec![Expr::int(1), Expr::mul(vec![Expr::int(2), i])]),
            Expr::int(2),
        );
        let e = Expr::func(FuncKind::Conj, vec![inner]);
        let r = eval(e.as_ref(), &ctx()).unwrap();
        assert_eq!(format_expr(r.as_ref()), "-3-4*i");
    }
}
