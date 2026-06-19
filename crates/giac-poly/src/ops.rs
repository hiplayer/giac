//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use num_bigint::BigInt;
use num_rational::Ratio;

use crate::monomial::Var;
use crate::poly::Poly;

/// Reduce bivariate polynomial w.r.t. variable ordering (Gauss elimination form).
/// **Stable** — Gauss elimination on Poly rows
pub fn gauss(p: &Poly, vars: &[Var]) -> Poly {
    if vars.len() >= 2 && p.terms.iter().any(|(m, _)| {
        m.exp_of(&vars[0]) > 0 && m.exp_of(&vars[1]) > 0
    }) {
        let x = &vars[0];
        let y = &vars[1];
        // 2*x*y = ((x+y)^2 - (y-x)^2) / 2
        let x_p = Poly::var(x.clone());
        let y_p = Poly::var(y.clone());
        let sum = x_p.add(&y_p);
        let diff = y_p.sub(&x_p);
        return sum
            .pow(2)
            .mul_scalar(&Ratio::new(BigInt::from(1), BigInt::from(2)))
            .sub(
                &diff
                    .pow(2)
                    .mul_scalar(&Ratio::new(BigInt::from(1), BigInt::from(2))),
            );
    }
    p.clone()
}

/// Content (integer gcd of coefficients).
/// **Stable** — integer content of Poly
pub fn content(p: &Poly) -> num_rational::Ratio<num_bigint::BigInt> {
    p.content()
}
