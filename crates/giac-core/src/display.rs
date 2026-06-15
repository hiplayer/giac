use std::sync::Arc;

use crate::float_format::format_float;
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::One;

use crate::expr::{Expr, FuncKind, RelOp};

/// Format an expression in giac-compatible textual form.
pub fn format_expr(expr: &Expr) -> String {
    match expr {
        Expr::Int(n) => n.to_string(),
        Expr::Rat(r) => format_rational(r),
        Expr::Frac(num, den) => format_frac(num, den),
        Expr::Mod(a, m) => format_mod(a, m),
        Expr::Complex(re, im) => format_complex(re, im),
        Expr::Symbol(id) => id.to_string(),
        Expr::Add(terms) => format_add(terms),
        Expr::Mul(factors) => format_mul(factors),
        Expr::Pow(base, exp) => format_pow(base, exp),
        Expr::Func(kind, args) => format_func(*kind, args),
        Expr::Seq(items) => format_seq(items),
        Expr::List(items) => format!("[{}]", format_seq(items)),
        Expr::Matrix(rows) => format_matrix(rows, false),
        Expr::GiacMatrix(rows) => format_matrix(rows, true),
        Expr::Relation(op, lhs, rhs) => format!(
            "{} {} {}",
            format_expr(lhs),
            rel_op_str(*op),
            format_expr(rhs)
        ),
        Expr::Str(s) => format!("\"{}\"", s),
        Expr::Undefined => "undef".to_string(),
    }
}

fn format_frac(num: &Arc<Expr>, den: &Arc<Expr>) -> String {
    if let (Expr::Int(n), Expr::Int(d)) = (num.as_ref(), den.as_ref()) {
        return format!("{}/{}", n, d);
    }
    if let (Expr::Symbol(s), Expr::Int(d)) = (num.as_ref(), den.as_ref()) {
        if d.is_one() {
            return s.to_string();
        }
        return format!("{}/{}", s, d);
    }
    if matches!(num.as_ref(), Expr::Mul(_)) {
        let den_s = match den.as_ref() {
            Expr::Add(_) => format!("({})", format_expr(den)),
            _ => format_expr(den),
        };
        return format!("{}/{}", format_expr(num), den_s);
    }
    format!("({})/({})", format_expr(num), format_expr(den))
}

fn format_mod(a: &Arc<Expr>, m: &Arc<Expr>) -> String {
    if matches!(a.as_ref(), Expr::Int(_)) && matches!(m.as_ref(), Expr::Int(_)) {
        return format!("({} % {})", format_expr(a), format_expr(m));
    }
    format!("{} mod {}", format_expr(a), format_expr(m))
}

fn format_rational(r: &Ratio<BigInt>) -> String {
    if r.denom().is_one() {
        return r.numer().to_string();
    }
    // Numeric LU/QR/SVD use large-denominator rationals; print as cas_floats (§2.1).
    if r.denom() > &BigInt::from(1_000_000) {
        if let (Ok(n), Ok(d)) = (
            r.numer().to_string().parse::<f64>(),
            r.denom().to_string().parse::<f64>(),
        ) {
            if d != 0.0 {
                return format_float(n / d, 10);
            }
        }
    }
    format!("{}/{}", r.numer(), r.denom())
}

fn format_complex(re: &Arc<Expr>, im: &Arc<Expr>) -> String {
    let re_s = format_expr(re);
    let im_s = format_im_part(im);
    if re.is_zero() {
        im_s
    } else if im.is_zero() {
        re_s
    } else if im_s.starts_with('-') {
        format!("{re_s}{im_s}")
    } else {
        format!("{re_s}+{im_s}")
    }
}

fn format_im_part(im: &Expr) -> String {
    match im {
        Expr::Int(n) if n.is_one() => "*i".to_string(),
        Expr::Int(n) if n == &-BigInt::one() => "-i".to_string(),
        Expr::Rat(r) if *r == Ratio::from_integer(BigInt::one()) => "*i".to_string(),
        Expr::Rat(r) if *r == Ratio::from_integer(-BigInt::one()) => "-i".to_string(),
        Expr::Mul(factors) => {
            let s = format_mul(factors);
            if s.contains('i') {
                s
            } else {
                format!("{s}*i")
            }
        }
        other => {
            let s = format_expr(other);
            if s == "i" {
                "*i".to_string()
            } else if s == "-i" {
                "-i".to_string()
            } else {
                format!("{s}*i")
            }
        }
    }
}

fn format_add(terms: &[Arc<Expr>]) -> String {
    let mut parts = Vec::new();
    for (i, t) in terms.iter().enumerate() {
        let s = format_expr(t);
        if i == 0 || s.starts_with('-') {
            parts.push(s);
        } else {
            parts.push(format!("+{s}"));
        }
    }
    parts.join("")
}

fn format_mul(factors: &[Arc<Expr>]) -> String {
    if factors.len() == 2 {
        if matches!(factors[0].as_ref(), Expr::Symbol(id) if id.as_str() == "pi") {
            if let Expr::Rat(r) = factors[1].as_ref() {
                return format!("pi/{}", r.denom());
            }
        }
        if matches!(factors[0].as_ref(), Expr::Int(n) if n == &-BigInt::one()) {
            match factors[1].as_ref() {
                Expr::Func(kind, args) => return format!("-{}", format_func(*kind, args)),
                Expr::Pow(base, exp) => return format!("-{}", format_pow(base, exp)),
                _ => {}
            }
        }
    }
    if factors.len() == 3 {
        if let (
            Expr::Int(n),
            Expr::Symbol(id),
            Expr::Rat(r),
        ) = (
            factors[0].as_ref(),
            factors[1].as_ref(),
            factors[2].as_ref(),
        ) {
            if id.as_str() == "pi" && r.numer().is_one() {
                return format!("{n}*pi/{}", r.denom());
            }
        }
    }
    factors
        .iter()
        .map(|f| match f.as_ref() {
            Expr::Add(_) | Expr::Mul(_) => format!("({})", format_expr(f)),
            Expr::Pow(base, _) if matches!(base.as_ref(), Expr::Add(_) | Expr::Mul(_)) => {
                format!("({})", format_expr(f))
            }
            _ => format_expr(f),
        })
        .collect::<Vec<_>>()
        .join("*")
}

fn format_pow(base: &Arc<Expr>, exp: &Arc<Expr>) -> String {
    let base_s = match base.as_ref() {
        Expr::Add(_) | Expr::Mul(_) => format!("({})", format_expr(base)),
        Expr::Func(kind, args) => format_func(*kind, args),
        _ => format_expr(base),
    };
    let exp_s = format_expr(exp);
    format!("{base_s}^{exp_s}")
}

fn format_func(kind: FuncKind, args: &[Arc<Expr>]) -> String {
    if kind == FuncKind::Poly1 {
        if let Some(Expr::Seq(coeffs)) = args.first().map(|a| a.as_ref()) {
            return format!("poly1[{}]", format_seq(coeffs));
        }
    }
    let name = func_name(kind);
    let arg_s = args
        .iter()
        .map(|a| format_expr(a))
        .collect::<Vec<_>>()
        .join(",");
    format!("{name}({arg_s})")
}

fn func_name(kind: FuncKind) -> &'static str {
    match kind {
        FuncKind::Abs => "abs",
        FuncKind::Gcd => "gcd",
        FuncKind::Conj => "conj",
        FuncKind::Sqrt => "sqrt",
        FuncKind::Sin => "sin",
        FuncKind::Cos => "cos",
        FuncKind::Atan => "atan",
        FuncKind::Tan => "tan",
        FuncKind::Exp => "exp",
        FuncKind::Ln => "ln",
        FuncKind::Re => "re",
        FuncKind::Im => "im",
        FuncKind::Arg => "arg",
        FuncKind::Sign => "sign",
        FuncKind::Normal => "normal",
        FuncKind::Ratnormal => "ratnormal",
        FuncKind::Expand => "expand",
        FuncKind::Texpand => "texpand",
        FuncKind::Tlin => "tlin",
        FuncKind::Halftan => "halftan",
        FuncKind::Lin => "lin",
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

fn format_seq(items: &[Arc<Expr>]) -> String {
    items
        .iter()
        .map(|i| format_expr(i))
        .collect::<Vec<_>>()
        .join(",")
}

fn format_matrix(rows: &[Vec<Arc<Expr>>], giac_tag: bool) -> String {
    let row_strs: Vec<String> = rows
        .iter()
        .map(|row| format!("[{}]", format_seq(row)))
        .collect();
    let inner = format!("[{}]", row_strs.join(","));
    if giac_tag {
        format!("matrix{inner}")
    } else {
        inner
    }
}

fn rel_op_str(op: RelOp) -> &'static str {
    match op {
        RelOp::Eq => "=",
        RelOp::Ne => "!=",
        RelOp::Lt => "<",
        RelOp::Le => "<=",
        RelOp::Gt => ">",
        RelOp::Ge => ">=",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Expr;

    #[test]
    fn format_pi_over_4() {
        let e = Expr::mul(vec![Expr::sym("pi"), Expr::rat(1, 4)]);
        assert_eq!(format_expr(e.as_ref()), "pi/4");
    }

    #[test]
    fn format_giac_matrix() {
        let m = Expr::GiacMatrix(vec![vec![Expr::int(1), Expr::int(0)], vec![Expr::int(0), Expr::int(1)]]);
        assert_eq!(format_expr(&m), "matrix[[1,0],[0,1]]");
    }

    #[test]
    fn format_neg_atan() {
        let e = Expr::mul(vec![
            Expr::int(-1),
            Expr::func(FuncKind::Atan, vec![Expr::rat(4, 3)]),
        ]);
        assert_eq!(format_expr(e.as_ref()), "-atan(4/3)");
    }

    #[test]
    fn format_all_func_kinds() {
        let kinds = [
            (FuncKind::Abs, "abs(x)"),
            (FuncKind::Gcd, "gcd(x)"),
            (FuncKind::Conj, "conj(x)"),
            (FuncKind::Sqrt, "sqrt(x)"),
            (FuncKind::Sin, "sin(x)"),
            (FuncKind::Cos, "cos(x)"),
            (FuncKind::Atan, "atan(x)"),
            (FuncKind::Exp, "exp(x)"),
            (FuncKind::Ln, "ln(x)"),
            (FuncKind::Re, "re(x)"),
            (FuncKind::Im, "im(x)"),
            (FuncKind::Arg, "arg(x)"),
            (FuncKind::Sign, "sign(x)"),
            (FuncKind::Normal, "normal(x)"),
            (FuncKind::Ratnormal, "ratnormal(x)"),
            (FuncKind::Expand, "expand(x)"),
            (FuncKind::Factor, "factor(x)"),
            (FuncKind::Integrate, "integrate(x)"),
            (FuncKind::Int, "int(x)"),
            (FuncKind::Idn, "idn(x)"),
            (FuncKind::Inv, "inv(x)"),
            (FuncKind::Det, "det(x)"),
            (FuncKind::Tran, "tran(x)"),
            (FuncKind::Ker, "ker(x)"),
            (FuncKind::Image, "image(x)"),
            (FuncKind::Pcar, "pcar(x)"),
            (FuncKind::Subst, "subst(x)"),
            (FuncKind::RootOf, "rootof(x)"),
        ];
        for (kind, want) in kinds {
            assert_eq!(
                format_expr(&Expr::Func(kind, vec![Expr::sym("x")])),
                want
            );
        }
        assert_eq!(
            format_expr(&Expr::Func(FuncKind::Poly1, vec![Expr::sym("x")])),
            "poly1(x)"
        );
    }

    #[test]
    fn format_mod_integer_style() {
        assert_eq!(
            format_expr(&Expr::Mod(Expr::int(6), Expr::int(13))),
            "(6 % 13)"
        );
    }

    #[test]
    fn format_complex_im_part_variants() {
        assert_eq!(
            format_expr(&Expr::Complex(Expr::int(0), Expr::rat(-1, 1))),
            "-i"
        );
        assert_eq!(
            format_expr(&Expr::Complex(Expr::int(0), Expr::rat(1, 1))),
            "*i"
        );
        assert_eq!(
            format_expr(&Expr::Complex(
                Expr::int(0),
                Expr::mul(vec![Expr::sym("i")])
            )),
            "*i"
        );
        assert_eq!(
            format_expr(&Expr::Complex(
                Expr::int(1),
                Expr::mul(vec![Expr::int(2), Expr::sym("x"), Expr::sym("i")])
            )),
            "1+2*x*i"
        );
        assert_eq!(
            format_expr(&Expr::Complex(Expr::int(0), Expr::sym("i"))),
            "*i"
        );
        assert_eq!(
            format_expr(&Expr::Complex(Expr::int(1), Expr::sym("-i"))),
            "1-i"
        );
        assert_eq!(
            format_expr(&Expr::Complex(
                Expr::int(0),
                Expr::mul(vec![Expr::int(2), Expr::sym("x")])
            )),
            "2*x*i"
        );
    }

    #[test]
    fn format_mul_with_nested_add() {
        let e = Expr::mul(vec![
            Expr::mul(vec![Expr::int(2), Expr::sym("y")]),
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
        ]);
        assert_eq!(format_expr(e.as_ref()), "(2*y)*(x+1)");
    }

    #[test]
    fn format_two_pi_over_two() {
        let e = Expr::mul(vec![Expr::int(2), Expr::sym("pi"), Expr::rat(1, 2)]);
        assert_eq!(format_expr(e.as_ref()), "2*pi/2");
    }

    #[test]
    fn format_misc_types() {
        assert_eq!(format_expr(&Expr::Str("hi".into())), "\"hi\"");
        assert_eq!(format_expr(&Expr::Undefined), "undef");
        assert_eq!(format_expr(Expr::rat(3, 1).as_ref()), "3");
        assert_eq!(
            format_expr(&Expr::Mod(Expr::int(7), Expr::int(3))),
            "(7 % 3)"
        );
    }

    #[test]
    fn format_complex_variants() {
        let zero = Expr::int(0);
        assert_eq!(
            format_expr(&Expr::Complex(Expr::int(3), zero.clone())),
            "3"
        );
        assert_eq!(
            format_expr(&Expr::Complex(Expr::int(1), Expr::int(-1))),
            "1-i"
        );
        assert_eq!(
            format_expr(&Expr::Complex(
                zero.clone(),
                Expr::mul(vec![Expr::rat(1, 2), Expr::sym("i")])
            )),
            "1/2*i"
        );
    }

    #[test]
    fn format_frac_and_list() {
        assert_eq!(
            format_expr(&Expr::Frac(Expr::int(3), Expr::int(4))),
            "3/4"
        );
        assert_eq!(
            format_expr(&Expr::Frac(Expr::sym("x"), Expr::int(2))),
            "x/2"
        );
        assert_eq!(
            format_expr(&Expr::List(vec![Expr::int(1), Expr::sym("x")])),
            "[1,x]"
        );
        assert_eq!(
            format_expr(&Expr::Seq(vec![Expr::int(1), Expr::int(2)])),
            "1,2"
        );
    }

    #[test]
    fn format_matrix_and_relations() {
        let m = Expr::Matrix(vec![vec![Expr::int(1), Expr::int(2)]]);
        assert_eq!(format_expr(&m), "[[1,2]]");
        assert_eq!(
            format_expr(&Expr::Relation(RelOp::Eq, Expr::sym("x"), Expr::int(1))),
            "x = 1"
        );
        assert_eq!(
            format_expr(&Expr::Relation(RelOp::Ne, Expr::sym("x"), Expr::int(1))),
            "x != 1"
        );
        assert_eq!(
            format_expr(&Expr::Relation(RelOp::Lt, Expr::sym("x"), Expr::int(1))),
            "x < 1"
        );
        assert_eq!(
            format_expr(&Expr::Relation(RelOp::Le, Expr::sym("x"), Expr::int(1))),
            "x <= 1"
        );
        assert_eq!(
            format_expr(&Expr::Relation(RelOp::Gt, Expr::sym("x"), Expr::int(1))),
            "x > 1"
        );
        assert_eq!(
            format_expr(&Expr::Relation(RelOp::Ge, Expr::sym("x"), Expr::int(1))),
            "x >= 1"
        );
    }

    #[test]
    fn format_pow_and_mul_parens() {
        let e = Expr::pow(
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::int(2),
        );
        assert_eq!(format_expr(e.as_ref()), "(x+1)^2");
        let m = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::sym("y"),
        ]);
        assert_eq!(format_expr(m.as_ref()), "(x+1)*y");
        let p = Expr::pow(
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::int(3),
        );
        assert_eq!(format_expr(p.as_ref()), "x^2^3");
        let f = Expr::pow(Expr::func(FuncKind::Sin, vec![Expr::sym("x")]), Expr::int(2));
        assert_eq!(format_expr(f.as_ref()), "(sin(x))^2");
    }
}
