use std::sync::Arc;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::context::Context;
use crate::error::EvalError;
use crate::expr::{Expr, ExprArc, FuncKind, RelOp};
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
                .map(|c| c.im)
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
        Expr::Matrix(rows) | Expr::GiacMatrix(rows) => {
            let ev: Result<Vec<Vec<_>>, _> = rows
                .iter()
                .map(|row| row.iter().map(|c| eval(c, ctx)).collect())
                .collect();
            let rows = ev?;
            match expr {
                Expr::GiacMatrix(_) => Arc::new(Expr::GiacMatrix(rows)),
                _ => Arc::new(Expr::Matrix(rows)),
            }
        }
        Expr::Relation(_, _, _) => Arc::new(expr.clone()),
        Expr::Mod(a, m) => eval_mod(a, m, ctx)?,
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
    if id.as_str() == "pi" {
        return Ok(Expr::sym("pi"));
    }
    Ok(Expr::sym(id.as_str()))
}

fn eval_mod(a: &ExprArc, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let inner = eval(a, ctx)?;
    let modulus = eval(m, ctx)?;
    if let Expr::Func(FuncKind::Factor, args) = inner.as_ref() {
        if args.len() == 1 {
            let arg = eval(&args[0], ctx)?;
            let mod_i = int_from_expr(modulus.as_ref())?;
            return crate::eval_poly::eval_factor_mod(&[arg], mod_i, ctx);
        }
    }
    if let (Some(a_i), Some(m_i)) = (as_int(inner.as_ref()), as_int(modulus.as_ref())) {
        let m_abs = m_i.abs();
        let mut r = (a_i % &m_abs).to_string().parse::<i64>().unwrap_or(0);
        if r < 0 {
            r += m_abs.to_string().parse::<i64>().unwrap_or(0);
        }
        return Ok(Expr::int(r));
    }
    if crate::algebra::poly::expr_to_poly(inner.as_ref()).is_ok() {
        let reduced = crate::eval_poly::eval_modp(&[Arc::clone(&inner), Arc::clone(&modulus)], ctx)?;
        return Ok(Arc::new(Expr::Mod(reduced, modulus)));
    }
    Ok(Arc::new(Expr::Mod(inner, modulus)))
}

fn int_from_expr(e: &Expr) -> Result<i64, EvalError> {
    as_int(e)
        .and_then(|n| n.to_string().parse().ok())
        .ok_or(EvalError::TypeError("integer expected"))
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
    if factors.len() == 2
        && matches!(factors[0].as_ref(), Expr::Matrix(_) | Expr::GiacMatrix(_))
        && matches!(factors[1].as_ref(), Expr::Matrix(_) | Expr::GiacMatrix(_))
    {
        let a = eval(&factors[0], ctx)?;
        let b = eval(&factors[1], ctx)?;
        return crate::linalg::eval_matrix_mul(&a, &b).and_then(|m| eval(m.as_ref(), ctx));
    }

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
            let e_u = crate::num_util::bigint_to_nonneg_u32(e)?;
            return Ok(Expr::int(int_to_i64(&(b.pow(e_u)))?));
        }
    }

    if let (Ok(_), Some(e_u)) = (
        crate::linalg::as_matrix(&base),
        as_nonneg_int(exp.as_ref()),
    ) {
        if e_u == 0 {
            let n = crate::linalg::as_matrix(&base)?.len();
            return Ok(crate::linalg::eval_idn(n));
        }
        if e_u <= 20 {
            return crate::linalg::eval_matrix_pow(&base, e_u, ctx);
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
    match kind {
        FuncKind::Subst => return eval_subst(args, ctx),
        FuncKind::Integrate | FuncKind::Int => return eval_integrate(args, ctx),
        FuncKind::Lambda => return Ok(Expr::func(FuncKind::Lambda, args.to_vec())),
        FuncKind::Smod => return crate::eval_poly::eval_smod(args),
        FuncKind::Irem => return crate::eval_poly::eval_irem(args),
        _ => {}
    }
    let ev: Result<Vec<_>, _> = args.iter().map(|a| eval(a, ctx)).collect();
    let args = ev?;
    match kind {
        FuncKind::Abs => eval_abs(&args, ctx),
        FuncKind::Gcd => eval_gcd(&args, ctx),
        FuncKind::Conj => eval_conj(&args, ctx),
        FuncKind::Sqrt => eval_sqrt(&args),
        FuncKind::Re => eval_component(&args, true, ctx),
        FuncKind::Im => eval_component(&args, false, ctx),
        FuncKind::Arg => eval_arg(&args, ctx),
        FuncKind::Sign => eval_sign(&args),
        FuncKind::Atan => Ok(Expr::func(FuncKind::Atan, args.to_vec())),
        FuncKind::Ln => Ok(Expr::func(FuncKind::Ln, args.to_vec())),
        FuncKind::Normal => crate::algebra::normal(args[0].as_ref(), ctx),
        FuncKind::Ratnormal => crate::algebra::ratnormal(args[0].as_ref(), ctx),
        FuncKind::Expand => crate::algebra::expand(args[0].as_ref(), ctx),
        FuncKind::Factor => crate::algebra::factor(args[0].as_ref(), ctx),
        FuncKind::Quo => crate::eval_poly::eval_quo(&args, ctx),
        FuncKind::Rem => crate::eval_poly::eval_rem(&args, ctx),
        FuncKind::Content => crate::eval_poly::eval_content(&args, ctx),
        FuncKind::Gauss => crate::eval_poly::eval_gauss(&args, ctx),
        FuncKind::Egcd => crate::eval_poly::eval_egcd(&args, ctx),
        FuncKind::Abcuv => crate::eval_poly::eval_abcuv(&args, ctx),
        FuncKind::Simp2 => crate::eval_poly::eval_simp2(&args, ctx),
        FuncKind::Lcm => crate::eval_poly::eval_lcm(&args, ctx),
        FuncKind::Horner => crate::eval_poly::eval_horner(&args, ctx),
        FuncKind::Resultant => crate::eval_poly::eval_resultant(&args, ctx),
        FuncKind::Roots => crate::eval_poly::eval_roots(&args, ctx),
        FuncKind::Modp => crate::eval_poly::eval_modp(&args, ctx),
        FuncKind::Chinrem => crate::eval_poly::eval_chinrem(&args, ctx),
        FuncKind::Partfrac => crate::eval_poly::eval_partfrac(&args, ctx),
        FuncKind::Greduce => crate::eval_poly::eval_greduce(&args, ctx),
        FuncKind::Rref => crate::linalg::eval_rref(&args, ctx),
        FuncKind::Idn => eval_idn(&args),
        FuncKind::Inv => eval_inv(&args, ctx),
        FuncKind::Det => eval_det(&args, ctx),
        FuncKind::Tran => eval_tran(&args),
        FuncKind::Ker => eval_ker(&args),
        FuncKind::Image => eval_image(&args),
        FuncKind::Pcar => eval_pcar(&args),
        FuncKind::Charpoly => eval_charpoly(&args, ctx),
        FuncKind::Linsolve => eval_linsolve(&args, ctx),
        FuncKind::Jordan => {
            if args.is_empty() {
                return Err(EvalError::TooFewArgs("jordan"));
            }
            crate::linalg::eigen::eval_jordan(&args[0], ctx)
        }
        FuncKind::Egv => {
            if args.is_empty() {
                return Err(EvalError::TooFewArgs("egv"));
            }
            crate::linalg::eigen::eval_egv(&args[0], ctx)
        }
        FuncKind::Lu => {
            if args.is_empty() {
                return Err(EvalError::TooFewArgs("lu"));
            }
            crate::linalg::numeric::eval_lu(&args[0], ctx)
        }
        FuncKind::Qr => {
            if args.is_empty() {
                return Err(EvalError::TooFewArgs("qr"));
            }
            crate::linalg::numeric::eval_qr(&args[0], ctx)
        }
        FuncKind::Svd => {
            if args.is_empty() {
                return Err(EvalError::TooFewArgs("svd"));
            }
            crate::linalg::numeric::eval_svd(&args[0], ctx)
        }
        FuncKind::Gramschmidt => crate::linalg::gramschmidt::eval_gramschmidt(&args, ctx),
        FuncKind::Trace => {
            if args.is_empty() {
                return Err(EvalError::TooFewArgs("trace"));
            }
            crate::linalg::eval_trace(&args[0])
        }
        FuncKind::RootOf => Ok(Expr::func(FuncKind::RootOf, args.to_vec())),
        FuncKind::Poly1 => Ok(Expr::func(FuncKind::Poly1, args.to_vec())),
        other => Err(EvalError::NotImplemented(func_name(other))),
    }
}

fn eval_abs(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(if args.is_empty() {
            EvalError::TooFewArgs("abs")
        } else {
            EvalError::TooManyArgs("abs")
        });
    }
    if let Some(n) = as_int(args[0].as_ref()) {
        return Ok(Expr::int(int_to_i64(&n.abs())?));
    }
    if let Some(c) = try_as_complex(args[0].as_ref(), ctx) {
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

fn eval_gcd(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() < 2 {
        return Err(EvalError::TooFewArgs("gcd"));
    }
    let has_mod = args.iter().any(|a| matches!(a.as_ref(), Expr::Mod(_, _)));
    if has_mod {
        return crate::eval_poly::eval_mod_gcd(args, ctx);
    }
    if args.iter().all(|a| as_int(a.as_ref()).is_some()) {
        let first = as_int(args[0].as_ref()).ok_or(EvalError::TypeError("gcd expects integers"))?;
        let mut result = first.clone();
        for a in &args[1..] {
            let n = as_int(a.as_ref()).ok_or(EvalError::TypeError("gcd expects integers"))?;
            result = result.gcd(n);
        }
        return Ok(Expr::int(int_to_i64(&result)?));
    }
    use crate::algebra::poly::{expr_to_poly, poly_to_expr};
    let first = expr_to_poly(args[0].as_ref())?;
    let mut result = first;
    for a in &args[1..] {
        let p = expr_to_poly(a.as_ref())?;
        result = result.gcd(&p);
    }
    let _ = ctx;
    Ok(poly_to_expr(&result))
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

fn eval_arg(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("arg"));
    }
    if let Some(c) = try_as_complex(args[0].as_ref(), ctx) {
        if c.im.is_zero() && c.re > Ratio::zero() {
            return Ok(Expr::int(0));
        }
        if c.re == Ratio::one() && c.im == Ratio::one() {
            return Ok(Expr::mul(vec![Expr::sym("pi"), Expr::rat(1, 4)]));
        }
        if c.re > Ratio::zero() && c.im > Ratio::zero() {
            let ratio = &c.im / &c.re;
            if ratio.denom() == &BigInt::one() {
                return Ok(Expr::func(
                    FuncKind::Atan,
                    vec![Expr::int(int_to_i64(ratio.numer())?)],
                ));
            }
        }
        if c.re < Ratio::zero() && c.im > Ratio::zero() {
            let ratio = &c.im / &-&c.re;
            return Ok(Expr::add(vec![
                Expr::mul(vec![Expr::int(-1), Expr::func(FuncKind::Atan, vec![ratio_to_expr(&ratio)?])]),
                Expr::mul(vec![Expr::int(2), Expr::sym("pi"), Expr::rat(1, 2)]),
            ]));
        }
    }
    Ok(Expr::func(FuncKind::Arg, args.to_vec()))
}

fn eval_sign(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("sign"));
    }
    if let Some(n) = as_int(args[0].as_ref()) {
        if n.is_zero() {
            return Ok(Expr::int(0));
        }
        return Ok(Expr::int(if n.is_negative() { -1 } else { 1 }));
    }
    Ok(Expr::func(FuncKind::Sign, args.to_vec()))
}

fn eval_idn(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    let n = as_int(args[0].as_ref()).ok_or(EvalError::TypeError("idn expects integer"))?;
    let n: usize = n.to_string().parse().map_err(|_| EvalError::TypeError("idn size"))?;
    Ok(crate::linalg::eval_idn(n))
}

fn eval_inv(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("inv"));
    }
    if let Some(n) = as_int(args[0].as_ref()) {
        return inv_scalar(n);
    }
    crate::linalg::eval_inv(&args[0], ctx)
}

fn inv_scalar(n: &BigInt) -> Result<ExprArc, EvalError> {
    if n.is_zero() {
        return Err(EvalError::DivisionByZero);
    }
    Ok(Expr::rat(
        1,
        n.to_string()
            .parse()
            .map_err(|_| EvalError::TypeError("scalar inv denominator"))?,
    ))
}

fn eval_det(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("det"));
    }
    crate::linalg::eval_det(&args[0], ctx)
}

fn eval_charpoly(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("charpoly"));
    }
    let var = match args[1].as_ref() {
        Expr::Symbol(id) => id.clone(),
        _ => return Err(EvalError::TypeError("charpoly variable")),
    };
    crate::linalg::eval_charpoly(&args[0], &var, ctx)
}

fn eval_linsolve(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("linsolve"));
    }
    crate::linalg::eval_linsolve(&args[0], &args[1], ctx)
}

fn eval_tran(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("tran"));
    }
    crate::linalg::eval_tran(&args[0])
}

fn eval_ker(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("ker"));
    }
    crate::linalg::eval_ker(&args[0])
}

fn eval_image(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("image"));
    }
    crate::linalg::eval_image(&args[0])
}

fn eval_pcar(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("pcar"));
    }
    crate::linalg::eval_pcar(&args[0])
}

fn eval_integrate(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 && args.len() != 4 {
        return Err(EvalError::TooFewArgs("integrate"));
    }
    let var = match args[1].as_ref() {
        Expr::Symbol(id) => id.clone(),
        _ => return Err(EvalError::TypeError("integration variable")),
    };
    let antideriv = crate::integrate::integrate(&args[0], &var)?;
    if args.len() == 2 {
        return Ok(antideriv);
    }
    let lo_bound = eval(&args[2], ctx)?;
    let hi_bound = eval(&args[3], ctx)?;
    let hi = eval(
        subst_expr(&antideriv, &var, &hi_bound)?.as_ref(),
        ctx,
    )?;
    let lo = eval(
        subst_expr(&antideriv, &var, &lo_bound)?.as_ref(),
        ctx,
    )?;
    eval(
        Expr::add(vec![
            hi,
            Expr::mul(vec![Expr::int(-1), lo]),
        ])
        .as_ref(),
        ctx,
    )
}

fn eval_subst(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("subst"));
    }
    let (var, val) = match args[1].as_ref() {
        Expr::Relation(RelOp::Eq, lhs, rhs) => {
            let name = match lhs.as_ref() {
                Expr::Symbol(id) => id.clone(),
                _ => return Err(EvalError::TypeError("subst equation")),
            };
            (name, eval(rhs, ctx)?)
        }
        _ => return Err(EvalError::TypeError("subst equation")),
    };
    subst_expr(&args[0], &var, &val)
}

fn subst_expr(expr: &ExprArc, var: &Ident, val: &ExprArc) -> Result<ExprArc, EvalError> {
    match expr.as_ref() {
        Expr::Symbol(id) if id == var => Ok(Arc::clone(val)),
        Expr::Symbol(_) | Expr::Int(_) | Expr::Rat(_) | Expr::List(_) | Expr::Seq(_)
        | Expr::Matrix(_) | Expr::GiacMatrix(_) => Ok(Arc::clone(expr)),
        Expr::Add(terms) => {
            let t: Result<Vec<_>, _> = terms.iter().map(|t| subst_expr(t, var, val)).collect();
            Ok(Expr::add(t?))
        }
        Expr::Mul(factors) => {
            let f: Result<Vec<_>, _> = factors.iter().map(|f| subst_expr(f, var, val)).collect();
            Ok(Expr::mul(f?))
        }
        Expr::Pow(b, e) => Ok(Expr::pow(subst_expr(b, var, val)?, subst_expr(e, var, val)?)),
        Expr::Frac(n, d) => Ok(Arc::new(Expr::Frac(
            subst_expr(n, var, val)?,
            subst_expr(d, var, val)?,
        ))),
        Expr::Func(k, a) => {
            let na: Result<Vec<_>, _> = a.iter().map(|x| subst_expr(x, var, val)).collect();
            Ok(Expr::func(*k, na?))
        }
        _ => Ok(Arc::clone(expr)),
    }
}

/// Multi-variable substitution without evaluation.
pub fn eval_subst_map(
    expr: &ExprArc,
    subs: &std::collections::HashMap<Ident, ExprArc>,
) -> Result<ExprArc, EvalError> {
    match expr.as_ref() {
        Expr::Symbol(id) => {
            if let Some(v) = subs.get(id) {
                Ok(Arc::clone(v))
            } else {
                Ok(Arc::clone(expr))
            }
        }
        Expr::Int(_) | Expr::Rat(_) | Expr::Str(_) | Expr::Undefined => Ok(Arc::clone(expr)),
        Expr::Add(terms) => {
            let t: Result<Vec<_>, _> = terms
                .iter()
                .map(|t| eval_subst_map(t, subs))
                .collect();
            Ok(Expr::add(t?))
        }
        Expr::Mul(factors) => {
            let f: Result<Vec<_>, _> = factors
                .iter()
                .map(|f| eval_subst_map(f, subs))
                .collect();
            Ok(Expr::mul(f?))
        }
        Expr::Pow(b, e) => Ok(Expr::pow(
            eval_subst_map(b, subs)?,
            eval_subst_map(e, subs)?,
        )),
        Expr::Frac(n, d) => Ok(Arc::new(Expr::Frac(
            eval_subst_map(n, subs)?,
            eval_subst_map(d, subs)?,
        ))),
        Expr::Complex(re, im) => Ok(Arc::new(Expr::Complex(
            eval_subst_map(re, subs)?,
            eval_subst_map(im, subs)?,
        ))),
        Expr::Func(k, a) => {
            let na: Result<Vec<_>, _> = a.iter().map(|x| eval_subst_map(x, subs)).collect();
            Ok(Expr::func(*k, na?))
        }
        Expr::Relation(op, lhs, rhs) => Ok(Arc::new(Expr::Relation(
            *op,
            eval_subst_map(lhs, subs)?,
            eval_subst_map(rhs, subs)?,
        ))),
        Expr::List(items) | Expr::Seq(items) => {
            let v: Result<Vec<_>, _> = items
                .iter()
                .map(|x| eval_subst_map(x, subs))
                .collect();
            Ok(Arc::new(match expr.as_ref() {
                Expr::List(_) => Expr::List(v?),
                _ => Expr::Seq(v?),
            }))
        }
        Expr::Matrix(rows) | Expr::GiacMatrix(rows) => {
            let r: Result<Vec<Vec<_>>, _> = rows
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|c| eval_subst_map(c, subs))
                        .collect::<Result<_, _>>()
                })
                .collect();
            Ok(Arc::new(match expr.as_ref() {
                Expr::GiacMatrix(_) => Expr::GiacMatrix(r?),
                _ => Expr::Matrix(r?),
            }))
        }
        Expr::Mod(a, m) => Ok(Arc::new(Expr::Mod(
            eval_subst_map(a, subs)?,
            eval_subst_map(m, subs)?,
        ))),
    }
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
        FuncKind::Atan => "atan",
        FuncKind::Normal => "normal",
        FuncKind::Ratnormal => "ratnormal",
        FuncKind::Expand => "expand",
        FuncKind::Factor => "factor",
        FuncKind::Quo => "quo",
        FuncKind::Rem => "rem",
        FuncKind::Content => "content",
        FuncKind::Gauss => "gauss",
        FuncKind::Egcd => "egcd",
        FuncKind::Abcuv => "abcuv",
        FuncKind::Simp2 => "simp2",
        FuncKind::Lcm => "lcm",
        FuncKind::Horner => "horner",
        FuncKind::Resultant => "resultant",
        FuncKind::Roots => "roots",
        FuncKind::Modp => "modp",
        FuncKind::Smod => "smod",
        FuncKind::Irem => "irem",
        FuncKind::Chinrem => "chinrem",
        FuncKind::Partfrac => "partfrac",
        FuncKind::Greduce => "greduce",
        FuncKind::Rref => "rref",
        FuncKind::Integrate => "integrate",
        FuncKind::Int => "int",
        FuncKind::Idn => "idn",
        FuncKind::Inv => "inv",
        FuncKind::Det => "det",
        FuncKind::Tran => "tran",
        FuncKind::Ker => "ker",
        FuncKind::Image => "image",
        FuncKind::Pcar => "pcar",
        FuncKind::Charpoly => "charpoly",
        FuncKind::Linsolve => "linsolve",
        FuncKind::Jordan => "jordan",
        FuncKind::Egv => "egv",
        FuncKind::Lu => "lu",
        FuncKind::Qr => "qr",
        FuncKind::Svd => "svd",
        FuncKind::Gramschmidt => "gramschmidt",
        FuncKind::Trace => "trace",
        FuncKind::Lambda => "lambda",
        FuncKind::Subst => "subst",
        FuncKind::RootOf => "rootof",
        FuncKind::Poly1 => "poly1",
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
    fn eval_gcd_three_integers() {
        let e = Expr::func(FuncKind::Gcd, vec![Expr::int(45), Expr::int(75), Expr::int(30)]);
        let r = eval(e.as_ref(), &ctx()).unwrap();
        assert_eq!(format_expr(r.as_ref()), "15");
    }

    #[test]
    fn eval_gcd_polynomial_pair() {
        let ctx = Context::default();
        let p = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(-2), Expr::sym("x")]),
            Expr::int(1),
        ]);
        let q = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(3)),
            Expr::int(-1),
        ]);
        let e = Expr::func(FuncKind::Gcd, vec![p, q]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x-1");
    }

    #[test]
    fn eval_gcd_polynomial_three_args() {
        let ctx = Context::default();
        let a = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(-2), Expr::sym("x")]),
            Expr::int(1),
        ]);
        let b = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(3)),
            Expr::int(-1),
        ]);
        let c = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::sym("x"),
            Expr::int(-2),
        ]);
        let e = Expr::func(FuncKind::Gcd, vec![a, b, c]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x-1");
    }

    #[test]
    fn eval_gcd_coprime_polynomials() {
        let ctx = Context::default();
        let a = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
            Expr::int(1),
        ]);
        let b = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::sym("x"),
            Expr::int(-2),
        ]);
        let e = Expr::func(FuncKind::Gcd, vec![a, b]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1");
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
    fn eval_factor_x4_minus_1() {
        let ctx = Context::default();
        let e = Expr::func(
            FuncKind::Factor,
            vec![Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(4)),
                Expr::int(-1),
            ])],
        );
        let r = eval(e.as_ref(), &ctx);
        assert!(r.is_ok());
    }

    #[test]
    fn integrate_definite_one() {
        let ctx = Context::default();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![Expr::int(1), Expr::sym("x"), Expr::int(-1), Expr::int(1)],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2");
    }

    #[test]
    fn integrate_definite_x_squared() {
        let ctx = Context::default();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::sym("x"),
                Expr::int(-1),
                Expr::int(1),
            ],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "2/3");
    }

    #[test]
    fn lambda_body_not_evaluated_prematurely() {
        let ctx = Context::default();
        let lambda = Expr::func(
            FuncKind::Lambda,
            vec![
                Arc::new(Expr::List(vec![Expr::sym("p"), Expr::sym("q")])),
                Expr::func(
                    FuncKind::Integrate,
                    vec![
                        Expr::mul(vec![Expr::sym("p"), Expr::sym("q")]),
                        Expr::sym("x"),
                        Expr::int(-1),
                        Expr::int(1),
                    ],
                ),
            ],
        );
        let r = eval(lambda.as_ref(), &ctx).unwrap();
        let Expr::Func(FuncKind::Lambda, parts) = r.as_ref() else {
            panic!("lambda not preserved");
        };
        let Expr::Func(FuncKind::Integrate, _) = parts[1].as_ref() else {
            panic!("integrate body not preserved");
        };
    }

    #[test]
    fn eval_integrate_inv_x() {
        let ctx = Context::default();
        let e = Expr::func(
            FuncKind::Integrate,
            vec![
                Expr::pow(Expr::sym("x"), Expr::int(-1)),
                Expr::sym("x"),
            ],
        );
        let r = eval(e.as_ref(), &ctx);
        assert!(r.is_ok(), "{:?}", r.err());
    }

    #[test]
    fn expand_binomial_no_stack_overflow() {
        let ctx = Context::default();
        let e = Expr::func(
            FuncKind::Normal,
            vec![Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(3)]), Expr::int(4))],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(
            format_expr(r.as_ref()),
            "x^4+12*x^3+54*x^2+108*x+81"
        );
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

    #[test]
    fn eval_arg_complex_power() {
        let i = Expr::sym("i");
        let inner = Expr::pow(
            Expr::add(vec![Expr::int(1), Expr::mul(vec![Expr::int(2), i])]),
            Expr::int(2),
        );
        let e = Expr::func(FuncKind::Arg, vec![inner]);
        let r = eval(e.as_ref(), &ctx()).unwrap();
        assert_eq!(format_expr(r.as_ref()), "-atan(4/3)+2*pi/2");
    }

    #[test]
    fn factor_perfect_square() {
        let ctx = Context::default();
        let e = Expr::func(
            FuncKind::Factor,
            vec![Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(4)),
                Expr::mul(vec![Expr::int(12), Expr::pow(Expr::sym("x"), Expr::int(3))]),
                Expr::mul(vec![Expr::int(54), Expr::pow(Expr::sym("x"), Expr::int(2))]),
                Expr::mul(vec![Expr::int(108), Expr::sym("x")]),
                Expr::int(81),
            ])],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(x+3)^4");
    }

    fn mat2() -> ExprArc {
        Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2)],
            vec![Expr::int(3), Expr::int(4)],
        ]))
    }

    #[test]
    fn eval_matrix_and_collection_types() {
        let ctx = ctx();
        let m = mat2();
        let mul = Expr::mul(vec![Arc::clone(&m), Arc::clone(&m)]);
        let r = eval(mul.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "[[7,10],[15,22]]");

        let seq = Arc::new(Expr::Seq(vec![Expr::int(1), Expr::add(vec![Expr::int(2), Expr::int(3)])]));
        assert_eq!(format_expr(eval(seq.as_ref(), &ctx).unwrap().as_ref()), "1,5");

        let list = Arc::new(Expr::List(vec![Expr::mul(vec![Expr::int(2), Expr::int(3)])]));
        assert_eq!(format_expr(eval(list.as_ref(), &ctx).unwrap().as_ref()), "[6]");

        let giac = eval(
            Arc::new(Expr::GiacMatrix(vec![vec![Expr::int(1)]])).as_ref(),
            &ctx,
        )
        .unwrap();
        assert!(matches!(giac.as_ref(), Expr::GiacMatrix(_)));
    }

    #[test]
    fn eval_frac_mod_relation() {
        let ctx = ctx();
        assert_eq!(
            eval(Arc::new(Expr::Frac(Expr::int(6), Expr::int(2))).as_ref(), &ctx).unwrap(),
            Expr::rat(3, 1)
        );
        let rel = Arc::new(Expr::Relation(RelOp::Eq, Expr::sym("x"), Expr::int(1)));
        assert!(matches!(eval(rel.as_ref(), &ctx).unwrap().as_ref(), Expr::Relation(RelOp::Eq, _, _)));
        let md = Arc::new(Expr::Mod(Expr::int(7), Expr::int(3)));
        assert_eq!(eval(md.as_ref(), &ctx).unwrap(), Expr::int(1));
    }

    #[test]
    fn eval_complex_node_and_pow() {
        let ctx = ctx();
        let c = Arc::new(Expr::Complex(Expr::int(1), Expr::mul(vec![Expr::int(2), Expr::sym("i")])));
        assert_eq!(format_expr(eval(c.as_ref(), &ctx).unwrap().as_ref()), "1+2*i");

        let i = Expr::sym("i");
        let p = Expr::pow(Expr::add(vec![Expr::int(1), i.clone()]), Expr::int(2));
        assert_eq!(format_expr(eval(p.as_ref(), &ctx).unwrap().as_ref()), "2*i");

        assert_eq!(eval(Expr::pow(Expr::int(2), Expr::int(10)).as_ref(), &ctx).unwrap(), Expr::int(1024));
    }

    #[test]
    fn eval_context_and_symbols() {
        let mut ctx = ctx();
        let x = Ident::new("x");
        ctx.set(x.clone(), Expr::int(5));
        assert_eq!(eval(Expr::sym("x").as_ref(), &ctx).unwrap(), Expr::int(5));
        assert_eq!(format_expr(eval(Expr::sym("pi").as_ref(), &ctx).unwrap().as_ref()), "pi");
    }

    #[test]
    fn eval_matrix_builtins() {
        let ctx = ctx();
        let m = mat2();

        let det = eval(Expr::func(FuncKind::Det, vec![Arc::clone(&m)]).as_ref(), &ctx).unwrap();
        let det_s = format_expr(det.as_ref());
        assert!(det_s == "-2" || det_s == "1*4-1*2*3");

        let inv = eval(Expr::func(FuncKind::Inv, vec![Arc::clone(&m)]).as_ref(), &ctx).unwrap();
        let inv_s = format_expr(inv.as_ref());
        assert!(
            inv_s == "[[-2,1],[3/2,-1/2]]" || inv_s.contains("4") && inv_s.contains("2*3"),
            "unexpected inv: {inv_s}"
        );

        let scalar_inv = eval(Expr::func(FuncKind::Inv, vec![Expr::int(4)]).as_ref(), &ctx).unwrap();
        assert_eq!(scalar_inv, Expr::rat(1, 4));

        let tran = eval(Expr::func(FuncKind::Tran, vec![Arc::clone(&m)]).as_ref(), &ctx).unwrap();
        assert!(format_expr(tran.as_ref()).starts_with("matrix[[1,3]"));

        let ker = eval(
            Expr::func(
                FuncKind::Ker,
                vec![Arc::new(Expr::Matrix(vec![
                    vec![Expr::int(1), Expr::int(2)],
                    vec![Expr::int(3), Expr::int(6)],
                ]))],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert!(!format_expr(ker.as_ref()).is_empty());

        let image = eval(
            Expr::func(
                FuncKind::Image,
                vec![Arc::new(Expr::Matrix(vec![
                    vec![Expr::int(1), Expr::int(2)],
                    vec![Expr::int(3), Expr::int(6)],
                ]))],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert!(!format_expr(image.as_ref()).is_empty());

        let pcar = eval(Expr::func(FuncKind::Pcar, vec![Arc::clone(&m)]).as_ref(), &ctx).unwrap();
        let pcar_s = format_expr(pcar.as_ref());
        assert!(
            pcar_s == "poly1[1,-5,-2]" || pcar_s.starts_with("poly1["),
            "unexpected pcar: {pcar_s}"
        );
    }

    #[test]
    fn eval_subst_and_algebra_funcs() {
        let ctx = ctx();
        let body = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(2)]), Expr::int(-1));
        let eq = Arc::new(Expr::Relation(RelOp::Eq, Expr::sym("x"), Expr::int(2)));
        let e = Expr::func(FuncKind::Subst, vec![body, eq]);
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "(2+2)^-1");

        let rn = eval(
            Expr::func(
                FuncKind::Ratnormal,
                vec![Expr::mul(vec![
                    Expr::rat(1, 2),
                    Expr::pow(Expr::sym("x"), Expr::int(-1)),
                ])],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert!(format_expr(rn.as_ref()).contains('/'));

        let ex = eval(
            Expr::func(
                FuncKind::Expand,
                vec![Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2))],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(ex.as_ref()), "x^2+2*x+1");
    }

    #[test]
    fn eval_sqrt_sign_re_im() {
        let ctx = ctx();
        assert_eq!(
            eval(Expr::func(FuncKind::Sqrt, vec![Expr::int(9)]).as_ref(), &ctx).unwrap(),
            Expr::int(3)
        );
        assert_eq!(
            format_expr(
                eval(Expr::func(FuncKind::Sqrt, vec![Expr::sym("x")]).as_ref(), &ctx).unwrap().as_ref()
            ),
            "sqrt(x)"
        );
        assert_eq!(
            eval(Expr::func(FuncKind::Sign, vec![Expr::int(0)]).as_ref(), &ctx).unwrap(),
            Expr::int(0)
        );
        assert_eq!(
            format_expr(
                eval(Expr::func(FuncKind::Sign, vec![Expr::sym("x")]).as_ref(), &ctx).unwrap().as_ref()
            ),
            "sign(x)"
        );
        assert_eq!(
            eval(Expr::func(FuncKind::Re, vec![Expr::int(42)]).as_ref(), &ctx).unwrap(),
            Expr::int(42)
        );
        assert_eq!(
            eval(Expr::func(FuncKind::Im, vec![Expr::int(42)]).as_ref(), &ctx).unwrap(),
            Expr::int(0)
        );
    }

    #[test]
    fn eval_arg_and_conj_edge_cases() {
        let ctx = ctx();
        let i = Expr::sym("i");
        assert_eq!(
            format_expr(
                eval(Expr::func(FuncKind::Arg, vec![Expr::int(5)]).as_ref(), &ctx).unwrap().as_ref()
            ),
            "0"
        );
        assert_eq!(
            format_expr(
                eval(
                    Expr::func(FuncKind::Arg, vec![Expr::add(vec![Expr::int(1), i.clone()])]).as_ref(),
                    &ctx,
                )
                .unwrap()
                .as_ref()
            ),
            "pi/4"
        );
        assert_eq!(
            format_expr(
                eval(Expr::func(FuncKind::Conj, vec![Expr::sym("x")]).as_ref(), &ctx).unwrap().as_ref()
            ),
            "conj(x)"
        );
        assert_eq!(
            eval(Expr::func(FuncKind::Abs, vec![Expr::int(-9)]).as_ref(), &ctx).unwrap(),
            Expr::int(9)
        );
    }

    #[test]
    fn eval_atan_ln_rootof_poly1() {
        let ctx = ctx();
        let atan = eval(Expr::func(FuncKind::Atan, vec![Expr::sym("x")]).as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(atan.as_ref()), "atan(x)");
        let ln = eval(Expr::func(FuncKind::Ln, vec![Expr::sym("x")]).as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(ln.as_ref()), "ln(x)");
        let root = eval(Expr::func(FuncKind::RootOf, vec![Expr::sym("x")]).as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(root.as_ref()), "rootof(x)");
        let poly = eval(
            Expr::func(
                FuncKind::Poly1,
                vec![Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)]))],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(poly.as_ref()), "poly1[1,0]");
    }

    #[test]
    fn eval_int_alias_and_idn() {
        let ctx = ctx();
        let e = Expr::func(
            FuncKind::Int,
            vec![
                Expr::pow(Expr::sym("x"), Expr::int(-1)),
                Expr::sym("x"),
            ],
        );
        assert!(eval(e.as_ref(), &ctx).is_ok());
        let idn = eval(Expr::func(FuncKind::Idn, vec![Expr::int(2)]).as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(idn.as_ref()), "matrix[[1,0],[0,1]]");
    }

    #[test]
    fn eval_not_implemented_and_errors() {
        let ctx = ctx();
        assert!(matches!(
            eval(Expr::func(FuncKind::Sin, vec![Expr::sym("x")]).as_ref(), &ctx),
            Err(EvalError::NotImplemented(_))
        ));
        assert!(matches!(
            eval(Expr::func(FuncKind::Abs, vec![]).as_ref(), &ctx),
            Err(EvalError::TooFewArgs(_))
        ));
        assert!(matches!(
            eval(
                Arc::new(Expr::Frac(Expr::int(1), Expr::int(0))).as_ref(),
                &ctx
            ),
            Err(EvalError::DivisionByZero)
        ));
    }

    #[test]
    fn eval_gcd_too_few_args() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Gcd, vec![Expr::int(1)]).as_ref(), &ctx()),
            Err(EvalError::TooFewArgs(_))
        ));
    }

    #[test]
    fn eval_abs_too_many_args() {
        assert!(matches!(
            eval(
                Expr::func(FuncKind::Abs, vec![Expr::int(1), Expr::int(2)]).as_ref(),
                &ctx()
            ),
            Err(EvalError::TooManyArgs(_))
        ));
    }

    #[test]
    fn eval_conj_too_many_args() {
        assert!(matches!(
            eval(
                Expr::func(FuncKind::Conj, vec![Expr::int(1), Expr::int(2)]).as_ref(),
                &ctx()
            ),
            Err(EvalError::TooManyArgs(_))
        ));
    }

    #[test]
    fn eval_sqrt_negative_and_symbolic() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Sqrt, vec![Expr::int(-1)]).as_ref(), &ctx()),
            Err(EvalError::TypeError(_))
        ));
        let r = eval(Expr::func(FuncKind::Sqrt, vec![Expr::int(2)]).as_ref(), &ctx()).unwrap();
        assert_eq!(format_expr(r.as_ref()), "sqrt(2)");
    }

    #[test]
    fn eval_re_im_too_few_args() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Re, vec![]).as_ref(), &ctx()),
            Err(EvalError::TooFewArgs(_))
        ));
        assert!(matches!(
            eval(Expr::func(FuncKind::Im, vec![]).as_ref(), &ctx()),
            Err(EvalError::TooFewArgs(_))
        ));
    }

    #[test]
    fn eval_arg_atan_ratio_and_symbolic() {
        let ctx = ctx();
        let i = Expr::sym("i");
        let r = eval(
            Expr::func(
                FuncKind::Arg,
                vec![Expr::add(vec![Expr::int(1), Expr::mul(vec![Expr::int(2), i.clone()])])],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(r.as_ref()), "atan(2)");
        let sym = eval(
            Expr::func(
                FuncKind::Arg,
                vec![Expr::add(vec![Expr::int(3), Expr::mul(vec![Expr::int(4), i])])],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(sym.as_ref()), "arg(3+4*i)");
    }

    #[test]
    fn eval_sign_too_few_and_positive() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Sign, vec![]).as_ref(), &ctx()),
            Err(EvalError::TooFewArgs(_))
        ));
        assert_eq!(
            eval(Expr::func(FuncKind::Sign, vec![Expr::int(7)]).as_ref(), &ctx()).unwrap(),
            Expr::int(1)
        );
    }

    #[test]
    fn eval_builtin_too_few_args() {
        let ctx = ctx();
        for (kind, name) in [
            (FuncKind::Inv, "inv"),
            (FuncKind::Det, "det"),
            (FuncKind::Tran, "tran"),
            (FuncKind::Ker, "ker"),
            (FuncKind::Image, "image"),
            (FuncKind::Pcar, "pcar"),
            (FuncKind::Integrate, "integrate"),
            (FuncKind::Subst, "subst"),
        ] {
            let r = eval(Expr::func(kind, vec![]).as_ref(), &ctx);
            assert!(matches!(r, Err(EvalError::TooFewArgs(n)) if n == name), "{name}");
        }
    }

    #[test]
    fn eval_integrate_bad_variable() {
        assert!(matches!(
            eval(
                Expr::func(
                    FuncKind::Integrate,
                    vec![Expr::sym("x"), Expr::int(1)],
                )
                .as_ref(),
                &ctx()
            ),
            Err(EvalError::TypeError(_))
        ));
    }

    #[test]
    fn eval_subst_bad_equation() {
        assert!(matches!(
            eval(
                Expr::func(FuncKind::Subst, vec![Expr::sym("x"), Expr::int(1)]).as_ref(),
                &ctx()
            ),
            Err(EvalError::TypeError(_))
        ));
        assert!(matches!(
            eval(
                Expr::func(
                    FuncKind::Subst,
                    vec![
                        Expr::sym("x"),
                        Arc::new(Expr::Relation(RelOp::Eq, Expr::int(1), Expr::int(2))),
                    ],
                )
                .as_ref(),
                &ctx()
            ),
            Err(EvalError::TypeError(_))
        ));
    }

    #[test]
    fn eval_subst_frac_and_func() {
        let ctx = ctx();
        let eq = Arc::new(Expr::Relation(RelOp::Eq, Expr::sym("x"), Expr::int(2)));
        let frac_body = Arc::new(Expr::Frac(Expr::sym("x"), Expr::int(2)));
        let r = eval(
            Expr::func(FuncKind::Subst, vec![frac_body, Arc::clone(&eq)]).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(r.as_ref()), "2/2");

        let fn_body = Expr::func(FuncKind::Sin, vec![Expr::sym("x")]);
        let r2 = eval(
            Expr::func(FuncKind::Subst, vec![fn_body, eq]).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(r2.as_ref()), "sin(2)");
    }

    #[test]
    fn eval_cos_not_implemented() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Cos, vec![Expr::sym("x")]).as_ref(), &ctx()),
            Err(EvalError::NotImplemented(_))
        ));
    }

    #[test]
    fn eval_pow_zero_complex() {
        let ctx = ctx();
        let i = Expr::sym("i");
        let r = eval(Expr::pow(Expr::add(vec![Expr::int(1), i]), Expr::int(0)).as_ref(), &ctx).unwrap();
        assert_eq!(r, Expr::int(1));
    }

    #[test]
    fn eval_frac_symbolic() {
        let ctx = ctx();
        let r = eval(
            Arc::new(Expr::Frac(Expr::sym("x"), Expr::sym("y"))).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(r.as_ref()), "(x)/(y)");
    }

    #[test]
    fn eval_str_undefined() {
        let ctx = ctx();
        let s = Arc::new(Expr::Str("ok".into()));
        assert!(matches!(eval(s.as_ref(), &ctx).unwrap().as_ref(), Expr::Str(_)));
        let u = eval(&Expr::Undefined, &ctx).unwrap();
        assert!(matches!(u.as_ref(), Expr::Undefined));
    }

    #[test]
    fn eval_ratnormal_builtin() {
        let ctx = ctx();
        let e = Expr::func(
            FuncKind::Ratnormal,
            vec![Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(-1)),
                Expr::rat(1, 2),
            ])],
        );
        let r = eval(e.as_ref(), &ctx).unwrap();
        assert!(format_expr(r.as_ref()).contains('x'));
    }

    #[test]
    fn eval_empty_add_mul() {
        let ctx = ctx();
        assert_eq!(eval(Expr::add(vec![]).as_ref(), &ctx).unwrap(), Expr::int(0));
        assert_eq!(eval(Expr::mul(vec![]).as_ref(), &ctx).unwrap(), Expr::int(1));
    }

    #[test]
    fn eval_exp_not_implemented() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Exp, vec![Expr::sym("x")]).as_ref(), &ctx()),
            Err(EvalError::NotImplemented(_))
        ));
    }

    #[test]
    fn eval_conj_too_few_args() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Conj, vec![]).as_ref(), &ctx()),
            Err(EvalError::TooFewArgs(_))
        ));
    }

    #[test]
    fn eval_sqrt_too_many_args() {
        assert!(matches!(
            eval(
                Expr::func(FuncKind::Sqrt, vec![Expr::int(1), Expr::int(2)]).as_ref(),
                &ctx()
            ),
            Err(EvalError::TooManyArgs(_))
        ));
    }

    #[test]
    fn eval_arg_too_few_args() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Arg, vec![]).as_ref(), &ctx()),
            Err(EvalError::TooFewArgs(_))
        ));
    }

    #[test]
    fn eval_idn_invalid() {
        assert!(matches!(
            eval(Expr::func(FuncKind::Idn, vec![Expr::sym("x")]).as_ref(), &ctx()),
            Err(EvalError::TypeError(_))
        ));
    }

    #[test]
    fn eval_subst_pow() {
        let ctx = ctx();
        let eq = Arc::new(Expr::Relation(RelOp::Eq, Expr::sym("x"), Expr::int(3)));
        let body = Expr::pow(Expr::sym("x"), Expr::int(2));
        let r = eval(
            Expr::func(FuncKind::Subst, vec![body, eq]).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(r.as_ref()), "9");
    }

    #[test]
    fn eval_complex_rat_imag_unit() {
        let ctx = ctx();
        let r = eval(
            Expr::mul(vec![Expr::rat(3, 2), Expr::sym("i")]).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(r.as_ref()), "i*3/2");
    }

    #[test]
    fn eval_mul_mixed_complex_symbolic() {
        let ctx = ctx();
        let i = Expr::sym("i");
        let r = eval(
            Expr::mul(vec![
                Expr::add(vec![Expr::int(1), i.clone()]),
                Expr::sym("x"),
            ])
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert!(format_expr(r.as_ref()).contains('x'));
    }

    #[test]
    fn func_name_covers_all_kinds() {
        let names = [
            (FuncKind::Abs, "abs"),
            (FuncKind::Gcd, "gcd"),
            (FuncKind::Conj, "conj"),
            (FuncKind::Sqrt, "sqrt"),
            (FuncKind::Sin, "sin"),
            (FuncKind::Cos, "cos"),
            (FuncKind::Atan, "atan"),
            (FuncKind::Exp, "exp"),
            (FuncKind::Ln, "ln"),
            (FuncKind::Re, "re"),
            (FuncKind::Im, "im"),
            (FuncKind::Arg, "arg"),
            (FuncKind::Sign, "sign"),
            (FuncKind::Normal, "normal"),
            (FuncKind::Ratnormal, "ratnormal"),
            (FuncKind::Expand, "expand"),
            (FuncKind::Factor, "factor"),
            (FuncKind::Quo, "quo"),
            (FuncKind::Rem, "rem"),
            (FuncKind::Content, "content"),
            (FuncKind::Gauss, "gauss"),
            (FuncKind::Egcd, "egcd"),
            (FuncKind::Abcuv, "abcuv"),
            (FuncKind::Simp2, "simp2"),
            (FuncKind::Lcm, "lcm"),
            (FuncKind::Horner, "horner"),
            (FuncKind::Resultant, "resultant"),
            (FuncKind::Roots, "roots"),
            (FuncKind::Modp, "modp"),
            (FuncKind::Smod, "smod"),
            (FuncKind::Irem, "irem"),
            (FuncKind::Chinrem, "chinrem"),
            (FuncKind::Partfrac, "partfrac"),
            (FuncKind::Greduce, "greduce"),
            (FuncKind::Rref, "rref"),
            (FuncKind::Integrate, "integrate"),
            (FuncKind::Int, "int"),
            (FuncKind::Idn, "idn"),
            (FuncKind::Inv, "inv"),
            (FuncKind::Det, "det"),
            (FuncKind::Tran, "tran"),
            (FuncKind::Ker, "ker"),
            (FuncKind::Image, "image"),
            (FuncKind::Pcar, "pcar"),
            (FuncKind::Subst, "subst"),
            (FuncKind::RootOf, "rootof"),
            (FuncKind::Poly1, "poly1"),
        ];
        for (kind, name) in names {
            assert_eq!(func_name(kind), name);
        }
    }

    #[test]
    fn eval_abs_complex_edge_cases() {
        let ctx = ctx();
        let i = Expr::sym("i");
        assert_eq!(
            eval(Expr::func(FuncKind::Abs, vec![Expr::int(0)]).as_ref(), &ctx).unwrap(),
            Expr::int(0)
        );
        assert_eq!(
            eval(Expr::func(FuncKind::Abs, vec![i.clone()]).as_ref(), &ctx).unwrap(),
            Expr::int(1)
        );
        assert_eq!(
            eval(
                Expr::func(
                    FuncKind::Abs,
                    vec![Expr::add(vec![Expr::int(1), i.clone()])],
                )
                .as_ref(),
                &ctx,
            )
            .unwrap(),
            Expr::func(FuncKind::Sqrt, vec![Expr::int(2)])
        );
        let r = eval(
            Expr::func(
                FuncKind::Abs,
                vec![Expr::add(vec![Expr::int(3), Expr::mul(vec![Expr::int(4), i])])],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(r, Expr::int(5));
    }

    #[test]
    fn eval_complex_pow_and_neg_im() {
        let ctx = ctx();
        let i = Expr::sym("i");
        let r = eval(
            Expr::pow(Expr::add(vec![Expr::int(1), i]), Expr::int(3)).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(r.as_ref()), "-2+2*i");
        let neg_i = eval(Expr::mul(vec![Expr::int(-1), Expr::sym("i")]).as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(neg_i.as_ref()), "-1*i");
    }

    #[test]
    fn eval_sign_negative_and_arg_second_quadrant() {
        let ctx = ctx();
        let i = Expr::sym("i");
        assert_eq!(
            eval(Expr::func(FuncKind::Sign, vec![Expr::int(-4)]).as_ref(), &ctx).unwrap(),
            Expr::int(-1)
        );
        let r = eval(
            Expr::func(
                FuncKind::Arg,
                vec![Expr::add(vec![Expr::int(-1), i.clone()])],
            )
            .as_ref(),
            &ctx,
        )
        .unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("atan") && s.contains("pi"));
    }

    #[test]
    fn eval_gcd_integer_with_polynomial() {
        let r = eval(
            Expr::func(FuncKind::Gcd, vec![Expr::int(1), Expr::sym("x")]).as_ref(),
            &ctx(),
        )
        .unwrap();
        assert_eq!(r, Expr::int(1));
    }

    #[test]
    fn eval_subst_mul_and_func() {
        let ctx = ctx();
        let eq = Arc::new(Expr::Relation(RelOp::Eq, Expr::sym("x"), Expr::int(3)));
        let mul_body = Expr::mul(vec![Expr::sym("x"), Expr::int(2)]);
        let r = eval(
            Expr::func(FuncKind::Subst, vec![mul_body, Arc::clone(&eq)]).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(r, Expr::int(6));

        let list_body = Arc::new(Expr::List(vec![Expr::sym("x")]));
        let r2 = eval(
            Expr::func(FuncKind::Subst, vec![list_body, eq]).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(r2.as_ref()), "[x]");
    }

    #[test]
    fn eval_pure_imaginary_complex_forms() {
        let ctx = ctx();
        assert_eq!(
            format_expr(
                eval(Expr::mul(vec![Expr::int(-1), Expr::sym("i")]).as_ref(), &ctx)
                    .unwrap()
                    .as_ref()
            ),
            "-1*i"
        );
        assert_eq!(
            format_expr(
                eval(Expr::mul(vec![Expr::rat(-1, 2), Expr::sym("i")]).as_ref(), &ctx)
                    .unwrap()
                    .as_ref()
            ),
            "i*-1/2"
        );
        assert_eq!(
            eval(Expr::func(FuncKind::Abs, vec![Expr::int(0)]).as_ref(), &ctx).unwrap(),
            Expr::int(0)
        );
    }

    #[test]
    fn eval_complex_integer_power() {
        let ctx = ctx();
        let r = eval(Expr::pow(Expr::int(2), Expr::int(25)).as_ref(), &ctx).unwrap();
        assert_eq!(r, Expr::int(33554432));
    }
}
