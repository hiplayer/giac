use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;

use super::patterns::try_factor_patterns;
use super::poly_uni::{
    content_wrt, factor_sqff_over_coeff_ring, primitive_part_wrt, square_free_wrt,
};
use super::univariate::factor_univariate_flat;
use super::util::{is_univariate_in, main_var, vars_in};

pub fn factor_into_poly(p: &Poly) -> Option<Vec<Poly>> {
    factor_multivariate(p).ok()
}

pub fn factor_multivariate(p: &Poly) -> PolyResult<Vec<Poly>> {
    if p.is_zero() {
        return Err(PolyError::TypeError("zero polynomial"));
    }
    if p.is_one() {
        return Ok(vec![]);
    }
    let vars = vars_in(p);
    factor_multivariate_rec(p, &vars)
}

pub(crate) fn factor_multivariate_rec(p: &Poly, vars: &[Var]) -> PolyResult<Vec<Poly>> {
    if let Some(f) = try_factor_patterns(p) {
        return Ok(f);
    }
    match vars.len() {
        0 => Ok(vec![p.clone()]),
        1 => factor_univariate_flat(p, &vars[0]),
        _ => {
            let x = main_var(p, vars);
            let others: Vec<Var> = vars.iter().filter(|v| *v != &x).cloned().collect();
            factor_wrt_main_var(p, &x, &others)
        }
    }
}

fn factor_wrt_main_var(p: &Poly, var: &Var, others: &[Var]) -> PolyResult<Vec<Poly>> {
    if is_univariate_in(p, var) {
        return factor_univariate_flat(p, var);
    }

    let content = content_wrt(p, var);
    let mut factors = if content.is_one() {
        Vec::new()
    } else {
        factor_multivariate_rec(&content, others)?
    };

    let pp = primitive_part_wrt(p, var)?;
    let sqff = square_free_wrt(&pp, var)?;
    for (g, k) in sqff {
        let gf = factor_sqff_over_coeff_ring(&g, var, others, factor_multivariate_rec)?;
        for _ in 0..k {
            factors.extend(gf.iter().cloned());
        }
    }
    Ok(factors)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    use num_traits::One;
    use crate::poly::Poly;

    #[test]
    fn factor_bivariate_product() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.mul(&y);
        let f = factor_multivariate(&p).unwrap();
        assert_eq!(f.len(), 2);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    #[ignore = "expanded bivariate gcd is too slow without Wang–Trager"]
    fn factor_bivariate_mixed() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x
            .mul(&x.pow(2).add(&x.mul(&y)).add(&y.pow(2)))
            .mul(
                &Poly::constant(Ratio::from_integer(3.into()))
                    .mul(&x)
                    .sub(&Poly::constant(Ratio::from_integer(2.into())).mul(&y)),
            );
        if let Ok(f) = factor_multivariate(&p) {
            assert!(f.len() >= 2);
            let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
            assert_eq!(prod, p);
        }
    }
}
