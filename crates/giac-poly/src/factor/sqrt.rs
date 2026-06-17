//! Factorization over Q(sqrt(d)) for quadratics with non-square discriminant.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{univariate_degree};

use super::util::{ratio_perfect_sqrt, vars_in};

/// Display factors `(x - (-b ± sqrt(disc))/(2a))` for a univariate quadratic.
pub fn quadratic_sqrt_factor_exprs(
    p: &Poly,
    var_name: &str,
) -> Option<Vec<(String, u32)>> {
    let vars = vars_in(p);
    if vars.len() != 1 {
        return None;
    }
    let var = &vars[0];
    if univariate_degree(p, var) != 2 {
        return None;
    }
    let mut a = Ratio::zero();
    let mut b = Ratio::zero();
    let mut c = Ratio::zero();
    for (m, coeff) in &p.terms {
        match m.exp_of(var) {
            2 => a += coeff,
            1 => b += coeff,
            0 => c += coeff,
            _ => return None,
        }
    }
    if a.is_zero() {
        return None;
    }
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

fn format_ratio(r: &Ratio<BigInt>) -> String {
    if r.denom().is_one() {
        r.numer().to_string()
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
}

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
    use num_traits::One;

    #[test]
    fn sqrt_factor_x2_minus_2() {
        let x = Var::from("x");
        let mut p = Poly::var(x.clone()).pow(2);
        p = p.sub(&Poly::constant(Ratio::from_integer(BigInt::from(2))));
        let facs = quadratic_sqrt_factor_exprs(&p, "x").expect("sqrt factors");
        assert_eq!(facs.len(), 2);
        assert!(facs[0].0.contains("sqrt(8)"));
    }
}
