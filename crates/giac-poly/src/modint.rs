//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::error::{PolyError, PolyResult};

/// Element of ℤ/mℤ.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModInt {
    pub val: BigInt,
    pub modulus: BigInt,
}

impl ModInt {
    /// **Stable** — `new`
    pub fn new(val: BigInt, modulus: BigInt) -> PolyResult<Self> {
        if modulus.is_zero() {
            return Err(PolyError::DivisionByZero);
        }
        Ok(Self {
            val: val.mod_floor(&modulus),
            modulus,
        })
    }

    /// **Stable** — `Poly::from_i64`
    pub fn from_i64(val: i64, modulus: i64) -> PolyResult<Self> {
        Self::new(BigInt::from(val), BigInt::from(modulus))
    }

    /// **Stable** — Poly addition
    pub fn add(&self, other: &Self) -> PolyResult<Self> {
        if self.modulus != other.modulus {
            return Err(PolyError::TypeError("modulus mismatch"));
        }
        Self::new(&self.val + &other.val, self.modulus.clone())
    }

    /// **Stable** — Poly subtraction
    pub fn sub(&self, other: &Self) -> PolyResult<Self> {
        if self.modulus != other.modulus {
            return Err(PolyError::TypeError("modulus mismatch"));
        }
        Self::new(&self.val - &other.val, self.modulus.clone())
    }

    /// **Stable** — Poly multiplication
    pub fn mul(&self, other: &Self) -> PolyResult<Self> {
        if self.modulus != other.modulus {
            return Err(PolyError::TypeError("modulus mismatch"));
        }
        Self::new(&self.val * &other.val, self.modulus.clone())
    }

    /// **Stable** — `inv`
    pub fn inv(&self) -> PolyResult<Self> {
        let eg = self.val.extended_gcd(&self.modulus);
        if eg.gcd != BigInt::one() {
            return Err(PolyError::TypeError("not invertible"));
        }
        Self::new(eg.x, self.modulus.clone())
    }

    /// **Stable** — Poly is zero
    pub fn is_zero(&self) -> bool {
        self.val.is_zero()
    }

    /// **Stable** — Poly is one
    pub fn is_one(&self) -> bool {
        self.val == BigInt::one()
    }
}

/// **Stable** — symmetric mod for i64
pub fn smod(a: i64, m: i64) -> i64 {
    let m = m.abs();
    let mut r = a % m;
    if r < 0 {
        r += m;
    }
    r
}

/// **Stable** — integer remainder
pub fn irem(a: i64, b: i64) -> i64 {
    a % b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smod_positive() {
        assert_eq!(smod(17, 5), 2);
        assert_eq!(irem(17, 5), 2);
    }
}
