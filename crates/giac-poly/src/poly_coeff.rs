//! Coefficient ring abstraction for [`super::poly::Poly<C>`](super::poly::Poly).
//!
//! `Ratio<BigInt>` is implemented here; `AlgExtC` lives in `giac-core` (avoids circular deps).

use std::fmt::Debug;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{EvalError, PolyResult};

/// Ring operations for sparse polynomial coefficients.
///
/// **Stable** — P1-1; implementations: `Ratio<BigInt>` (giac-poly), `AlgExtCPolyCoeff` (giac-core).
pub trait PolyCoeff: Clone + PartialEq + Debug {
    fn coeff_zero() -> Self;
    fn coeff_one() -> Self;
    fn coeff_is_zero(&self) -> bool;
    fn coeff_is_one(&self) -> bool;
    fn coeff_add(&self, rhs: &Self) -> PolyResult<Self>;
    fn coeff_sub(&self, rhs: &Self) -> PolyResult<Self>;
    fn coeff_neg(&self) -> PolyResult<Self>;
    fn coeff_mul(&self, rhs: &Self) -> PolyResult<Self>;
    fn coeff_div(&self, rhs: &Self) -> PolyResult<Self>;
}

impl PolyCoeff for Ratio<BigInt> {
    fn coeff_zero() -> Self {
        Ratio::zero()
    }

    fn coeff_one() -> Self {
        Ratio::one()
    }

    fn coeff_is_zero(&self) -> bool {
        self.is_zero()
    }

    fn coeff_is_one(&self) -> bool {
        self.is_one()
    }

    fn coeff_add(&self, rhs: &Self) -> PolyResult<Self> {
        Ok(self + rhs)
    }

    fn coeff_sub(&self, rhs: &Self) -> PolyResult<Self> {
        Ok(self - rhs)
    }

    fn coeff_neg(&self) -> PolyResult<Self> {
        Ok(-self.clone())
    }

    fn coeff_mul(&self, rhs: &Self) -> PolyResult<Self> {
        Ok(self * rhs)
    }

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
