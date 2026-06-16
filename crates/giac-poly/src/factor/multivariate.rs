use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;

use super::patterns::try_factor_patterns;
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
    if let Some(f) = try_factor_patterns(p) {
        return Ok(f);
    }
    let vars = vars_in(p);
    match vars.len() {
        0 => Ok(vec![p.clone()]),
        1 => factor_univariate_flat(p, &vars[0]),
        _ => {
            let v = main_var(p, &vars);
            if is_univariate_in(p, &v) {
                factor_univariate_flat(p, &v)
            } else {
                Err(PolyError::NotImplemented("factor"))
            }
        }
    }
}
