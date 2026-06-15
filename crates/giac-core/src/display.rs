use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::One;

use crate::expr::{Expr, FuncKind, RelOp};

/// Format an expression in giac-compatible textual form.
pub fn format_expr(expr: &Expr) -> String {
    match expr {
        Expr::Int(n) => n.to_string(),
        Expr::Rat(r) => format_rational(r),
        Expr::Frac(num, den) => format!("({})/({})", format_expr(num), format_expr(den)),
        Expr::Mod(a, m) => format!("{} mod {}", format_expr(a), format_expr(m)),
        Expr::Complex(re, im) => format_complex(re, im),
        Expr::Symbol(id) => id.to_string(),
        Expr::Add(terms) => format_add(terms),
        Expr::Mul(factors) => format_mul(factors),
        Expr::Pow(base, exp) => format_pow(base, exp),
        Expr::Func(kind, args) => format_func(*kind, args),
        Expr::Seq(items) => format_seq(items),
        Expr::List(items) => format!("[{}]", format_seq(items)),
        Expr::Matrix(rows) => format_matrix(rows),
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

fn format_rational(r: &Ratio<BigInt>) -> String {
    if r.denom().is_one() {
        r.numer().to_string()
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
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
        if i == 0 {
            parts.push(s);
        } else if s.starts_with('-') {
            parts.push(s);
        } else {
            parts.push(format!("+{s}"));
        }
    }
    parts.join("")
}

fn format_mul(factors: &[Arc<Expr>]) -> String {
    factors
        .iter()
        .map(|f| match f.as_ref() {
            Expr::Add(_) | Expr::Pow(_, _) => format!("({})", format_expr(f)),
            _ => format_expr(f),
        })
        .collect::<Vec<_>>()
        .join("*")
}

fn format_pow(base: &Arc<Expr>, exp: &Arc<Expr>) -> String {
    let base_s = match base.as_ref() {
        Expr::Add(_) | Expr::Mul(_) | Expr::Func(_, _) => format!("({})", format_expr(base)),
        _ => format_expr(base),
    };
    let exp_s = format_expr(exp);
    format!("{base_s}^{exp_s}")
}

fn format_func(kind: FuncKind, args: &[Arc<Expr>]) -> String {
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
        FuncKind::Exp => "exp",
        FuncKind::Ln => "ln",
        FuncKind::Re => "re",
        FuncKind::Im => "im",
        FuncKind::Arg => "arg",
        FuncKind::Sign => "sign",
    }
}

fn format_seq(items: &[Arc<Expr>]) -> String {
    items
        .iter()
        .map(|i| format_expr(i))
        .collect::<Vec<_>>()
        .join(",")
}

fn format_matrix(rows: &[Vec<Arc<Expr>>]) -> String {
    let row_strs: Vec<String> = rows
        .iter()
        .map(|row| format!("[{}]", format_seq(row)))
        .collect();
    format!("[{}]", row_strs.join(","))
}

fn rel_op_str(op: RelOp) -> &'static str {
    match op {
        RelOp::Eq => "==",
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
    fn format_integer_and_sqrt() {
        assert_eq!(format_expr(Expr::int(15).as_ref()), "15");
        assert_eq!(
            format_expr(Expr::func(FuncKind::Sqrt, vec![Expr::int(5)]).as_ref()),
            "sqrt(5)"
        );
    }

    #[test]
    fn format_complex_negative_imag() {
        let e = Expr::add(vec![
            Expr::int(-3),
            Expr::mul(vec![Expr::int(-4), Expr::sym("i")]),
        ]);
        assert_eq!(format_expr(e.as_ref()), "-3-4*i");
    }
}
