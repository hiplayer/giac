//! Quadratic factorization with sqrt display (Expr string helpers for giac-simplify).
//!
//! **Partial:** `quadratic_sqrt_factor_exprs`.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::poly::Poly;
use crate::resultant::{quadratic_abc};

use super::util::{ratio_perfect_sqrt, vars_in};

/// **Partial** — Display factors `(x - (-b ± sqrt(disc))/(2a))` for a univariate quadratic.
pub fn quadratic_sqrt_factor_exprs(
    p: &Poly,
    var_name: &str,
) -> Option<Vec<(String, u32)>> {
    let vars = vars_in(p);
    if vars.len() != 1 {
        return None;
    }
    let var = &vars[0];
    let (a, b, c) = quadratic_abc(p, var)?;
    let disc = &b * &b - Ratio::from_integer(BigInt::from(4)) * &a * &c;
    if disc.is_zero() || ratio_perfect_sqrt(&disc).is_some() {
        return None;
    }
    let two_a = Ratio::from_integer(BigInt::from(2)) * &a;
    let neg_b = -&b;
    let disc_str = format_sqrt_ratio(&disc);
    let r1 = format!("({} + {disc_str})/({})", format_ratio(&neg_b), format_ratio(&two_a));
    let r2 = format!("({} - {disc_str})/({})", format_ratio(&neg_b), format_ratio(&two_a));
    Some(vec![
        (format!("{var_name} - ({r1})"), 1),
        (format!("{var_name} - ({r2})"), 1),
    ])
}

// **Pipeline private** — `format_ratio`
fn format_ratio(r: &Ratio<BigInt>) -> String {
    if r.denom().is_one() {
        r.numer().to_string()
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
}

// **Pipeline private** — `format_sqrt_ratio`
fn format_sqrt_ratio(r: &Ratio<BigInt>) -> String {
    if r.is_one() {
        return "sqrt(1)".to_string();
    }
    if r.denom().is_one() {
        return format!("sqrt({})", r.numer());
    }
    format!("sqrt({}/{})", r.numer(), r.denom())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::monomial::Var;

    #[test]
    fn sqrt_factor_x2_minus_2() {
        let x = Var::from("x");
        let mut p = Poly::var(x.clone()).pow(2);
        p = p.sub(&Poly::constant(Ratio::from_integer(BigInt::from(2))));
        let facs = quadratic_sqrt_factor_exprs(&p, "x").expect("sqrt factors");
        assert_eq!(facs.len(), 2);
        assert_eq!(facs[0].0, "x - ((0 + sqrt(8))/(2))");
    }
}
