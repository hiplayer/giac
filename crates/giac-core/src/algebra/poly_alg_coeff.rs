//! [`PolyCoeff`](giac_poly::PolyCoeff) bridge for [`AlgExtCData`](super::alg_ext_c::AlgExtCData).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::fmt::{self, Debug};

use giac_poly::{FieldCoeff, PolyCoeff, PolyResult};

use super::alg_ext_c::AlgExtCData;

/// Coefficient wrapper for `Poly<AlgExtCPolyCoeff>` (path B in expr-poly-conversion).
#[derive(Clone, PartialEq)]
pub struct AlgExtCPolyCoeff(pub AlgExtCData);

impl AlgExtCPolyCoeff {
    /// **Stable** — `new`
    pub fn new(data: AlgExtCData) -> Self {
        Self(data)
    }

    /// **Stable** — `Poly::as_inner`
    pub fn as_inner(&self) -> &AlgExtCData {
        &self.0
    }

    /// **Stable** — `Poly::into_inner`
    pub fn into_inner(self) -> AlgExtCData {
        self.0
    }
}

impl Debug for AlgExtCPolyCoeff {
    // **Stable** — `Poly::fmt`
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("AlgExtCPolyCoeff").field(&self.0).finish()
    }
}

impl From<AlgExtCData> for AlgExtCPolyCoeff {
    // **Stable** — `Poly::from`
    fn from(data: AlgExtCData) -> Self {
        Self(data)
    }
}

impl FieldCoeff for AlgExtCPolyCoeff {}

impl PolyCoeff for AlgExtCPolyCoeff {
    // **Stable** — `Poly::coeff_zero`
    fn coeff_zero() -> Self {
        Self(
            AlgExtCData::zero(
                super::ext_tower::ExtensionField::rational(),
            )
            .expect("rational field zero"),
        )
    }

    // **Stable** — `Poly::coeff_one`
    fn coeff_one() -> Self {
        Self(
            AlgExtCData::one(
                super::ext_tower::ExtensionField::rational(),
            )
            .expect("rational field one"),
        )
    }

    // **Pipeline private** — `coeff_is_zero`
    fn coeff_is_zero(&self) -> bool {
        self.0.is_zero().unwrap_or(true)
    }

    // **Pipeline private** — `coeff_is_one`
    fn coeff_is_one(&self) -> bool {
        self.0.is_one().unwrap_or(false)
    }

    // **Pipeline private** — `coeff_add`
    fn coeff_add(&self, rhs: &Self) -> PolyResult<Self> {
        self.0.add(rhs.as_inner()).map(Self)
    }

    // **Pipeline private** — `coeff_sub`
    fn coeff_sub(&self, rhs: &Self) -> PolyResult<Self> {
        self.0.sub(rhs.as_inner()).map(Self)
    }

    // **Pipeline private** — `coeff_neg`
    fn coeff_neg(&self) -> PolyResult<Self> {
        self.0.neg().map(Self)
    }

    // **Pipeline private** — `coeff_mul`
    fn coeff_mul(&self, rhs: &Self) -> PolyResult<Self> {
        self.0.mul(rhs.as_inner()).map(Self)
    }

    // **Pipeline private** — `coeff_div`
    fn coeff_div(&self, rhs: &Self) -> PolyResult<Self> {
        let inv = rhs.0.inv()?;
        self.coeff_mul(&Self(inv))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_poly::PolyCoeff;

    use crate::algebra::test_fixtures::sqrt2_algext;
    use crate::expr::Expr;

    use super::*;

    #[test]
    fn algext_c_coeff_mul_sqrt2() {
        let alpha = sqrt2_algext();
        let z = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&alpha).unwrap());
        let sq = z.coeff_mul(&z).unwrap();
        let two = AlgExtCPolyCoeff::from(
            AlgExtCData::from_complex_parts(&Expr::int(2), &Expr::int(0)).unwrap(),
        );
        assert!(sq.as_inner().eq_mod(two.as_inner()).unwrap());
    }
}
