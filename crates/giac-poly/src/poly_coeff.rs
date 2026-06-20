//! Coefficient ring abstraction for [`super::poly::Poly<C>`](super::poly::Poly).
//!
//! `Ratio<BigInt>` is implemented here; `AlgExtC` lives in `giac-core` (avoids circular deps).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use std::fmt::Debug;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{EvalError, PolyResult};

/// Ring operations for sparse polynomial coefficients.
///
/// **Stable** — P1-1; implementations: `Ratio<BigInt>` (giac-poly), `AlgExtCPolyCoeff` (giac-core).
pub trait PolyCoeff: Clone + PartialEq + Debug {
    // **Pipeline private** — `coeff_zero`
    fn coeff_zero() -> Self;
    // **Pipeline private** — `coeff_one`
    fn coeff_one() -> Self;
    // **Pipeline private** — `coeff_is_zero`
    fn coeff_is_zero(&self) -> bool;
    // **Pipeline private** — `coeff_is_one`
    fn coeff_is_one(&self) -> bool;
    // **Pipeline private** — `coeff_add`
    fn coeff_add(&self, rhs: &Self) -> PolyResult<Self>;
    // **Pipeline private** — `coeff_sub`
    fn coeff_sub(&self, rhs: &Self) -> PolyResult<Self>;
    // **Pipeline private** — `coeff_neg`
    fn coeff_neg(&self) -> PolyResult<Self>;
    // **Pipeline private** — `coeff_mul`
    fn coeff_mul(&self, rhs: &Self) -> PolyResult<Self>;
    // **Pipeline private** — `coeff_div`
    fn coeff_div(&self, rhs: &Self) -> PolyResult<Self>;
}

impl PolyCoeff for Ratio<BigInt> {
    // **Stable** — `Poly::coeff_zero`
    fn coeff_zero() -> Self {
        Ratio::zero()
    }

    // **Stable** — `Poly::coeff_one`
    fn coeff_one() -> Self {
        Ratio::one()
    }

    // **Stable** — `Poly::coeff_is_zero`
    fn coeff_is_zero(&self) -> bool {
        self.is_zero()
    }

    // **Stable** — `Poly::coeff_is_one`
    fn coeff_is_one(&self) -> bool {
        self.is_one()
    }

    // **Pipeline private** — `coeff_add`
    fn coeff_add(&self, rhs: &Self) -> PolyResult<Self> {
        Ok(self + rhs)
    }

    // **Pipeline private** — `coeff_sub`
    fn coeff_sub(&self, rhs: &Self) -> PolyResult<Self> {
        Ok(self - rhs)
    }

    // **Pipeline private** — `coeff_neg`
    fn coeff_neg(&self) -> PolyResult<Self> {
        Ok(-self.clone())
    }

    // **Pipeline private** — `coeff_mul`
    fn coeff_mul(&self, rhs: &Self) -> PolyResult<Self> {
        Ok(self * rhs)
    }

    // **Pipeline private** — `coeff_div`
    fn coeff_div(&self, rhs: &Self) -> PolyResult<Self> {
        if rhs.is_zero() {
            Err(EvalError::DivisionByZero)
        } else {
            Ok(self / rhs)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ratio_coeff_ring() {
        let a = Ratio::from_integer(3.into());
        let b = Ratio::from_integer(2.into());
        assert_eq!(a.coeff_add(&b).unwrap(), Ratio::from_integer(5.into()));
        assert_eq!(a.coeff_mul(&b).unwrap(), Ratio::from_integer(6.into()));
        assert_eq!(b.coeff_div(&a).unwrap(), Ratio::new(2.into(), 3.into()));
    }
}
