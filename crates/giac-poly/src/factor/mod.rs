//! Polynomial factorization over ℚ (and ℤ/pℤ for `factor_poly_mod`).
//!
//! **Nested-ring gate:** ℚ[others][main] paths in `sparse` / `hensel` / `poly_uni` must use
//! [`crate::nested::UnivariateIn`] (`.divides`, `.div_rem_wrt_aux_indep`) — not [`Poly::div_rem`].
//!
//! **Upstream:** `gausspol.cc` / `ezgcd.cc` — `do_factor_hensel`, `factor`.
//! **缺口:** FAC-G2 一般参系数塔；`unitaryfactor` 二元 MVP 已接入 `poly_uni`（`factor/unitary.rs`）。
//!
//! | Tier | 入口 |
//! |------|------|
//! | **Stable (bounded)** | `factor_into`, `factor_poly`, `factor_into_by_rational_roots` |
//! | **Stable** | `factor_poly_mod`, `factor_mod_irreducibles` |
//!
//! 子模块：`multivariate` → `univariate` / `hensel` / `sparse` / `zassenhaus` / `patterns`。
//! 全函数 tier 见 `.doc/giac-poly-api-stability.md` § `factor/mod.rs`。

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
mod tower;
mod ctx;
mod eval;
mod unitary;
mod cyclotomic;
mod fpx;
mod fpx_uni;
mod hensel;
mod zassenhaus;
mod modular;
mod multivariate;
mod patterns;
mod poly_uni;
mod power;
mod sparse;
mod sqrt;
#[cfg(test)]
mod tracer;
mod univariate;
mod util;

use crate::poly::Poly;

pub use util::{ratio_perfect_sqrt, vars_in};
pub use eval::GoodEval;
pub use power::{as_perfect_power, try_linear_power};
pub use sqrt::quadratic_sqrt_factor_exprs;
pub use univariate::factor_power_pairs;
pub use univariate::{factor_univariate_flat, factor_univariate_pairs};
pub(crate) use univariate::find_rational_root;

use multivariate::factor_into_poly;
use patterns::{factor_xn_minus_one_display, try_factor_patterns};
use power::as_perfect_power as perfect_power;

/// **Stable (bounded)** — Factor into irreducible polynomial factors over ℚ when possible.
pub fn factor_into(p: &Poly) -> Option<Vec<Poly>> {
    factor_into_poly(p)
}

/// **Stable (bounded)** — Integer-style factorization display (legacy `factor_poly`).
pub fn factor_poly(p: &Poly) -> Poly {
    if let Some(f) = factor_xn_minus_one_display(p) {
        return f;
    }
    if let Some((base, exp)) = perfect_power(p) {
        if let Some(factors) = factor_into_poly(&base) {
            if factors.len() > 1 {
                return factors
                    .into_iter()
                    .fold(Poly::one(), |acc, f| acc.mul(&f))
                    .pow(exp);
            }
        }
        return base.pow(exp);
    }
    if let Some(factors) = factor_into_poly(p) {
        if factors.len() > 1 {
            return factors
                .into_iter()
                .fold(Poly::one(), |acc, f| acc.mul(&f));
        }
    }
    p.clone()
}

/// **Stable (bounded)** — univariate via rational roots
pub fn factor_into_by_rational_roots(p: &Poly, var: &crate::monomial::Var) -> crate::error::PolyResult<Vec<Poly>> {
    factor_univariate_flat(p, var)
}

/// **Stable** — factor mod p display
pub fn factor_poly_mod(p: &Poly, modulus: i64) -> Result<Poly, crate::error::EvalError> {
    modular::factor_poly_mod(p, modulus)
}

/// **Stable** — Irreducible factors over F_p (monic, with repetition).
pub fn factor_mod_irreducibles(p: &Poly, modulus: i64) -> crate::error::PolyResult<Vec<crate::modular::PolyMod>> {
    let pm = crate::modp(p, modulus)?;
    fpx::factor_fpx(&pm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_rational::Ratio;
    use num_traits::One;
    use crate::monomial::Var;
    use crate::poly::Poly;

    #[test]
    fn as_perfect_power_quadratic_squared() {
        let x = Poly::var("x");
        let p = x.pow(2).add(&Poly::one()).pow(2);
        let (base, exp) = as_perfect_power(&p).unwrap();
        assert_eq!(exp, 2);
        assert_eq!(base.pow(2), p);
    }

    #[test]
    fn factor_multivariate_returns_none_without_hanging() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.mul(&y);
        let f = factor_into(&p).expect("factor x*y");
        assert_eq!(f.len(), 2);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    fn factor_x_fourth_minus_one() {
        let p = Poly::var("x").pow(4).sub(&Poly::one());
        let f = factor_into(&p).expect("factor x^4-1");
        assert_eq!(f.len(), 3);
    }

    #[test]
    fn factor_x_cubed_plus_one() {
        let p = Poly::var("x").pow(3).add(&Poly::one());
        let f = factor_into(&p).expect("factor x^3+1");
        assert_eq!(f.len(), 2);
    }

    #[test]
    fn factor_x_times_x_squared_plus_one() {
        let p = Poly::var("x").mul(&Poly::var("x").pow(2).add(&Poly::one()));
        let f = factor_into(&p).expect("factor x*(x^2+1)");
        assert_eq!(f.len(), 2);
    }

    #[test]
    fn factor_x6_minus_y6() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.pow(6).sub(&y.pow(6));
        let f = factor_into(&p).expect("x^6-y^6");
        assert!(f.len() >= 4);
    }

    #[test]
    fn factor_x100_plus_x50_plus_1() {
        let x = Poly::var("x");
        let p = x.pow(100).add(&x.pow(50)).add(&Poly::one());
        let f = factor_into(&p).expect("factor x^100+x^50+1");
        assert_eq!(f.len(), 6);
        let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
        assert_eq!(prod, p);
    }

    #[test]
    fn try_linear_power_detects_square() {
        let p = Poly::var("x").sub(&Poly::one()).pow(2);
        let (base, exp) = try_linear_power(&p, &Var::from("x")).unwrap();
        assert_eq!(exp, 2);
        assert_eq!(base.pow(2), p);
    }
}
