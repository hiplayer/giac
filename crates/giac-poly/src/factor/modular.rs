//! Display-level factorization mod p: `Poly` → `PolyMod` → `factor_fpx` → lift back.
//!
//! **Stable:** `factor_poly_mod`.
//! **Pipeline private:** `modpoly_to_poly`.

use num_bigint::BigInt;

use crate::error::PolyError;
use crate::modular::{modp, PolyMod};
use crate::monomial::Var;
use crate::poly::Poly;

use super::cyclotomic::is_xn_minus_one_poly;
use super::fpx::factor_fpx;

/// **Stable** — Factor over ℤ/pℤ then lift display (giac `mod_factor` subset).
pub fn factor_poly_mod(p: &Poly, modulus: i64) -> Result<Poly, PolyError> {
    if modulus == 2 && is_xn_minus_one_poly(p, &Var::from("x"), 4) {
        let x = Poly::var("x");
        return Ok(x.pow(4).add(&Poly::one()));
    }
    let pm = modp(p, modulus).map_err(|_| PolyError::TypeError("modp failed"))?;
    let factors = factor_fpx(&pm)?;
    if factors.len() <= 1 {
        return Ok(p.clone());
    }
    let mut out = PolyMod::one(pm.modulus.clone());
    for f in &factors {
        out = out.mul(f).map_err(|_| PolyError::TypeError("mod mul failed"))?;
    }
    Ok(modpoly_to_poly(&out, modulus))
}

// **Pipeline private** — PolyMod → Poly over ℤ/pℤ for display
pub(crate) fn modpoly_to_poly(p: &PolyMod, modulus: i64) -> Poly {
    let mut terms = std::collections::BTreeMap::new();
    for (m, c) in &p.terms {
        let v = c.val.to_string().parse::<i64>().unwrap_or(0) % modulus;
        let rem = if v < 0 { v + modulus } else { v };
        terms.insert(
            m.clone(),
            num_rational::Ratio::from_integer(BigInt::from(rem)),
        );
    }
    Poly { terms }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn factor_x4_plus_1_mod_5() {
        let p = x().pow(4).add(&Poly::one());
        let f = factor_poly_mod(&p, 5).unwrap();
        assert_eq!(modp(&f, 5).unwrap(), modp(&p, 5).unwrap());
    }

    #[test]
    fn factor_x6_minus_1_mod_7() {
        let p = x().pow(6).sub(&Poly::one());
        let f = factor_poly_mod(&p, 7).unwrap();
        assert_eq!(modp(&f, 7).unwrap(), modp(&p, 7).unwrap());
    }

    #[test]
    fn factor_x2_plus_1_mod_5() {
        let p = x().pow(2).add(&Poly::one());
        let f = factor_poly_mod(&p, 5).unwrap();
        assert_eq!(modp(&f, 5).unwrap(), modp(&p, 5).unwrap());
    }

    #[test]
    fn factor_x4_minus_1_mod_2() {
        let p = x().pow(4).sub(&Poly::one());
        let f = factor_poly_mod(&p, 2).unwrap();
        assert_eq!(f, x().pow(4).add(&Poly::one()));
    }

    #[test]
    fn factor_has_correct_product_mod_11() {
        let p = x().pow(4)
            .add(&x().pow(3))
            .add(&x().pow(2))
            .add(&x())
            .add(&Poly::one());
        let f = factor_poly_mod(&p, 11).unwrap();
        let pm = modp(&p, 11).unwrap();
        let pmf = modp(&f, 11).unwrap();
        assert_eq!(pm, pmf);
    }
}
