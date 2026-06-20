//! [`PolyCoeff`](giac_poly::PolyCoeff) bridge for [`AlgExtCData`](super::alg_ext_c::AlgExtCData).

use std::fmt::{self, Debug};

use giac_poly::{PolyCoeff, PolyResult};

use super::alg_ext_c::AlgExtCData;

/// Coefficient wrapper for `Poly<AlgExtCPolyCoeff>` (path B in expr-poly-conversion).
#[derive(Clone, PartialEq)]
pub struct AlgExtCPolyCoeff(pub AlgExtCData);

impl AlgExtCPolyCoeff {
    pub fn new(data: AlgExtCData) -> Self {
        Self(data)
    }

    pub fn as_inner(&self) -> &AlgExtCData {
        &self.0
    }

    pub fn into_inner(self) -> AlgExtCData {
        self.0
    }
}

impl Debug for AlgExtCPolyCoeff {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("AlgExtCPolyCoeff").field(&self.0).finish()
    }
}

impl From<AlgExtCData> for AlgExtCPolyCoeff {
    fn from(data: AlgExtCData) -> Self {
        Self(data)
    }
}

impl PolyCoeff for AlgExtCPolyCoeff {
    fn coeff_zero() -> Self {
        Self(
            AlgExtCData::zero(
                super::ext_tower::ExtensionField::rational(),
            )
            .expect("rational field zero"),
        )
    }

    fn coeff_one() -> Self {
        Self(
            AlgExtCData::one(
                super::ext_tower::ExtensionField::rational(),
            )
            .expect("rational field one"),
        )
    }

    fn coeff_is_zero(&self) -> bool {
        self.0.is_zero().unwrap_or(true)
    }

    fn coeff_is_one(&self) -> bool {
        self.0.is_one().unwrap_or(false)
    }

    fn coeff_add(&self, rhs: &Self) -> PolyResult<Self> {
        self.0.add(rhs.as_inner()).map(Self)
    }

    fn coeff_sub(&self, rhs: &Self) -> PolyResult<Self> {
        self.0.sub(rhs.as_inner()).map(Self)
    }

    fn coeff_neg(&self) -> PolyResult<Self> {
        self.0.neg().map(Self)
    }

    fn coeff_mul(&self, rhs: &Self) -> PolyResult<Self> {
        self.0.mul(rhs.as_inner()).map(Self)
    }

    fn coeff_div(&self, rhs: &Self) -> PolyResult<Self> {
        let inv = rhs.0.inv()?;
        self.coeff_mul(&Self(inv))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_poly::PolyCoeff;

    use crate::expr::{Expr, FuncKind};
    use crate::AlgExtData;

    use super::*;

    #[test]
    fn algext_c_coeff_mul_sqrt2() {
        let min = Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ));
        let alpha = AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap();
        let z = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&alpha).unwrap());
        let sq = z.coeff_mul(&z).unwrap();
        let two = AlgExtCPolyCoeff::from(
            AlgExtCData::from_complex_parts(&Expr::int(2), &Expr::int(0)).unwrap(),
        );
        assert!(sq.as_inner().eq_mod(two.as_inner()).unwrap());
    }
}
