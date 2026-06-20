//! ℚ coefficient ring adapter for [`super::poly1::Poly1RingCtx`].
//!
//! **Tier:** Stable (crate-internal) — used by `giac-core::field_arith` facade (D2–D3).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{EvalError, PolyResult};

use super::poly1::Poly1RingCtx;

/// Marker for [`Ratio<BigInt>`] as dense poly1 coefficients.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RatioRingOps;

/// Ring context for dense poly1 over `Ratio<BigInt>`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RatioRingCtx;

impl Poly1RingCtx for RatioRingCtx {
    type Coeff = Ratio<BigInt>;

    // **Stable** — Poly zero
    fn zero(&self) -> Self::Coeff {
        Zero::zero()
    }

    // **Stable** — Poly one
    fn one(&self) -> Self::Coeff {
        One::one()
    }

    // **Stable** — Poly is zero
    fn is_zero(&self, c: &Self::Coeff) -> bool {
        Zero::is_zero(c)
    }

    // **Stable** — Poly addition
    fn add(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff> {
        Ok(a + b)
    }

    // **Stable** — Poly subtraction
    fn sub(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff> {
        Ok(a - b)
    }

    // **Stable** — Poly negation
    fn neg(&self, c: &Self::Coeff) -> PolyResult<Self::Coeff> {
        Ok(-c.clone())
    }

    // **Stable** — Poly multiplication
    fn mul(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff> {
        Ok(a * b)
    }

    // **Pipeline private** — `inv`
    fn inv(&self, c: &Self::Coeff) -> PolyResult<Self::Coeff> {
        if Zero::is_zero(c) {
            Err(EvalError::DivisionByZero)
        } else {
            Ok(<Ratio<BigInt> as One>::one() / c)
        }
    }
}
