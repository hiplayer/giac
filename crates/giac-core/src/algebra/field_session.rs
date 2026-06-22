//! Ambient **K** + working **L** for `poly_algext_roots` (G5 / PR-B′).
//!
//! Coefficients are lifted and aligned on **L**; `ambient` stays fixed for normalized
//! polynomial input. See [FieldSession plan](../../../../.doc/issues/GIAC-poly-roots-field-session-plan.md).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::sync::Arc;

use giac_poly::PolyCoeff;
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::One;

use crate::error::EvalError;

use super::alg_ext::{algext_cube_root, algext_square_roots, AlgExtData};
use super::alg_ext_c::AlgExtCData;
use super::ext_tower::ExtensionField;
use super::field_arith::coords_to_expr;
use super::poly_alg_coeff::AlgExtCPolyCoeff;

/// Session carrying ambient coefficient field **K** and monotonically growing working field **L**.
/// **Stable (bounded)** — roots pipeline field session
pub struct FieldSession {
    ambient: Arc<ExtensionField>,
    working: Arc<ExtensionField>,
}

impl FieldSession {
    /// **Stable** — `FieldSession::new`
    pub fn new(ambient: Arc<ExtensionField>) -> Self {
        let working = Arc::clone(&ambient);
        Self { ambient, working }
    }

    /// **Stable** — ambient K (normalized poly coefficients)
    pub fn ambient(&self) -> &Arc<ExtensionField> {
        &self.ambient
    }

    /// **Stable** — current working L
    pub fn working(&self) -> &Arc<ExtensionField> {
        &self.working
    }

    /// **Stable** — zero in L
    pub fn zero(&self) -> AlgExtCPolyCoeff {
        coeff_in_field(&self.working, |f| AlgExtCData::zero(Arc::clone(f)).expect("zero"))
    }

    /// **Stable** — one in L
    pub fn one(&self) -> AlgExtCPolyCoeff {
        coeff_in_field(&self.working, |f| AlgExtCData::one(Arc::clone(f)).expect("one"))
    }

    /// **Stable** — integer constant in L
    pub fn int(&self, n: i64) -> Result<AlgExtCPolyCoeff, EvalError> {
        rat(
            &self.working,
            Ratio::from_integer(BigInt::from(n)),
        )
    }

    /// **Stable** — 1/2 in L
    pub fn half(&self) -> Result<AlgExtCPolyCoeff, EvalError> {
        rat(&self.working, Ratio::new(1.into(), 2.into()))
    }

    /// Embed `c` into current **L** (may bump working).
    /// **Stable** — lift coeff to working field
    pub fn lift(&mut self, c: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        if Arc::ptr_eq(&c.as_inner().field, &self.working) {
            return Ok(c.clone());
        }
        if ExtensionField::is_subfield_of(&c.as_inner().field, &self.working) {
            let z = self.zero();
            let (lifted, _) = self.align(c, &z)?;
            return Ok(lifted);
        }
        if ExtensionField::is_subfield_of(&self.working, &c.as_inner().field) {
            self.bump_to(&c.as_inner().field);
            return Ok(c.clone());
        }
        let z = self.zero();
        let (lifted, _) = self.align(c, &z)?;
        Ok(lifted)
    }

    /// Align two coefficients to a common field; bumps **L** to the align target.
    /// **Stable** — align pair on working field
    pub fn align(
        &mut self,
        a: &AlgExtCPolyCoeff,
        b: &AlgExtCPolyCoeff,
    ) -> Result<(AlgExtCPolyCoeff, AlgExtCPolyCoeff), EvalError> {
        let (aa, bb) = AlgExtCData::align_pair(a.as_inner(), b.as_inner())?;
        self.bump_to(&aa.field);
        Ok((AlgExtCPolyCoeff::from(aa), AlgExtCPolyCoeff::from(bb)))
    }

    /// Formal i·z with z ∈ K embedded in current **L** (i² = −1, not adjoined to the tower).
    /// **Stable (bounded)** — multiply by formal i
    pub fn mul_formal_i(&mut self, z: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let z = self.lift(z)?;
        let inner = z.as_inner();
        Ok(AlgExtCPolyCoeff::from(AlgExtCData {
            field: Arc::clone(&inner.field),
            re: coords_to_expr(&inner.field.zero_coords())?,
            im: inner.re.clone(),
            root_index: None,
        }))
    }

    /// Canonical principal square root on **L** (may extend the tower).
    ///
    /// ℚ negative → i√|u|; real u → adjoin x²−u (or u²−u over K); see `sqrt_euler_options`
    /// when Euler needs both real and i·√(−u) candidates.
    /// **Stable (bounded)** — principal sqrt on working field
    pub fn sqrt_principal(&mut self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let u = self.lift(u)?;
        if u.coeff_is_zero() {
            return Ok(self.zero());
        }
        if is_negative_rational(&u) {
            let abs = self.neg(&u)?;
            let beta = self.adjoin_sqrt(&abs)?;
            return self.mul_formal_i(&beta);
        }
        self.adjoin_sqrt(&u)
    }

    /// Adjoin √u (real part) and return a square root generator in **L**.
    /// **Stable (bounded)** — adjoin sqrt primitive
    pub fn adjoin_sqrt(&mut self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let u = self.lift(u)?;
        let inner = u.as_inner();
        if !inner.im.iter().all(|e| e.is_zero()) {
            return Err(EvalError::NotImplemented("adjoin_sqrt complex"));
        }
        let re = AlgExtData::from_field_coords(Arc::clone(&inner.field), inner.re.clone())?;
        let mut roots = algext_square_roots(&re)?;
        let beta = roots
            .pop()
            .ok_or(EvalError::NotImplemented("adjoin_sqrt empty"))?;
        let out = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&beta)?);
        self.bump_to(&out.as_inner().field);
        Ok(out)
    }

    /// Adjoin ∛u (real part) and return a cube root generator in **L**.
    /// **Stable (bounded)** — adjoin cbrt primitive
    pub fn adjoin_cbrt(&mut self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        let u = self.lift(u)?;
        let inner = u.as_inner();
        if !inner.im.iter().all(|e| e.is_zero()) {
            return Err(EvalError::NotImplemented("adjoin_cbrt complex"));
        }
        let re = AlgExtData::from_field_coords(Arc::clone(&inner.field), inner.re.clone())?;
        let beta = algext_cube_root(&re)?;
        let out = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&beta)?);
        self.bump_to(&out.as_inner().field);
        Ok(out)
    }

    /// Adjoin a primitive cube root of unity ω (ω²+ω+1=0) over **L**.
    /// **Stable (bounded)** — adjoin ω for pure cubic roots
    pub fn adjoin_primitive_cube_root_of_unity(
        &mut self,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let parent = Arc::clone(self.working());
        let one = parent.one_coords();
        let layer = vec![one.clone(), one.clone(), one];
        let omega_field = if parent.dimension() == 1 {
            let one = Ratio::one();
            ExtensionField::adjoin_irreducible_over_q(vec![one.clone(), one.clone(), one])?
        } else {
            ExtensionField::adjoin_irreducible_parent_coeffs(&parent, layer)?
        };
        let coords = omega_field.generator_coords();
        let omega = AlgExtData::from_field_coords(
            Arc::clone(&omega_field),
            coords_to_expr(&coords)?,
        )?;
        let out = AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&omega)?);
        self.bump_to(&out.as_inner().field);
        Ok(out)
    }

    /// Reset **L** to a prior checkpoint (Euler branch search).
    /// **Pipeline private** — `set_working`
    pub(crate) fn set_working(&mut self, field: &Arc<ExtensionField>) {
        self.working = Arc::clone(field);
    }

    /// **Stable** — add with auto-align
    pub fn add(
        &mut self,
        a: &AlgExtCPolyCoeff,
        b: &AlgExtCPolyCoeff,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let (a, b) = self.align(a, b)?;
        a.coeff_add(&b)
    }

    /// **Stable** — multiply with auto-align
    pub fn mul(
        &mut self,
        a: &AlgExtCPolyCoeff,
        b: &AlgExtCPolyCoeff,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let (a, b) = self.align(a, b)?;
        a.coeff_mul(&b)
    }

    /// **Stable** — divide with auto-align
    pub fn div(
        &mut self,
        a: &AlgExtCPolyCoeff,
        b: &AlgExtCPolyCoeff,
    ) -> Result<AlgExtCPolyCoeff, EvalError> {
        let (a, b) = self.align(a, b)?;
        a.coeff_div(&b)
    }

    /// **Stable** — negate (no align)
    pub fn neg(&self, a: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError> {
        a.coeff_neg()
    }

    /// Grow **L** toward `new` when `new` is a superfield or align common field.
    /// **Pipeline private** — `bump_to`
    pub(crate) fn bump_to(&mut self, new: &Arc<ExtensionField>) {
        if Arc::ptr_eq(&self.working, new) {
            return;
        }
        if ExtensionField::is_subfield_of(&self.working, new) {
            self.working = Arc::clone(new);
            return;
        }
        if ExtensionField::is_subfield_of(new, &self.working) {
            return;
        }
        if new.dimension() >= self.working.dimension() {
            self.working = Arc::clone(new);
        }
    }
}

// **Pipeline private** — rational constant in field
fn rat(field: &Arc<ExtensionField>, r: Ratio<BigInt>) -> Result<AlgExtCPolyCoeff, EvalError> {
    let coords = field.embed_rational(&r);
    let a = AlgExtData::from_field_coords(Arc::clone(field), coords_to_expr(&coords)?)?;
    Ok(AlgExtCPolyCoeff::from(AlgExtCData::from_alg_ext(&a)?))
}

fn coeff_in_field(
    field: &Arc<ExtensionField>,
    f: impl FnOnce(&Arc<ExtensionField>) -> AlgExtCData,
) -> AlgExtCPolyCoeff {
    AlgExtCPolyCoeff::from(f(field))
}

// **Pipeline private** — negative constant in ℚ ⊂ K
fn is_negative_rational(c: &AlgExtCPolyCoeff) -> bool {
    use super::field_arith::{pad_to_len, rationalize_poly1};
    use num_traits::Zero;

    let inner = c.as_inner();
    if !inner.im.iter().all(|e| e.is_zero()) {
        return false;
    }
    let Ok(re) = rationalize_poly1(&inner.re) else {
        return false;
    };
    if inner.field.dimension() != 1 {
        return false;
    }
    pad_to_len(&re, 1)[0] < Ratio::zero()
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_poly::PolyCoeff;
    use num_rational::Ratio;
    use serial_test::serial;

    use crate::algebra::test_fixtures::{k1_adjoin_sqrt2, sqrt2_algext};

    use super::*;

    #[serial]
    #[test]
    fn int_and_one_on_rational() {
        let k = ExtensionField::rational();
        let mut session = FieldSession::new(k);
        let two = session.int(2).unwrap();
        let one = session.one();
        assert!(!two.coeff_is_zero());
        assert!(!one.coeff_is_zero());
        assert!(Arc::ptr_eq(session.working(), session.ambient()));
    }

    #[serial]
    #[test]
    fn lift_and_align_bump_working_on_k1() {
        let k1 = k1_adjoin_sqrt2();
        let mut session = FieldSession::new(Arc::clone(&k1));
        let sqrt2 = AlgExtCPolyCoeff::from(
            AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap(),
        );
        let four = session.int(4).unwrap();
        let disc = session.mul(&sqrt2, &four).unwrap();
        let dim_before = session.working().dimension();
        let beta = session.adjoin_sqrt(&disc).unwrap();
        assert!(
            session.working().dimension() > dim_before,
            "adjoin_sqrt should grow working field"
        );
        let sq = session.mul(&beta, &beta).unwrap();
        let (sq_a, disc_a) = session.align(&sq, &disc).unwrap();
        assert!(sq_a.coeff_sub(&disc_a).unwrap().coeff_is_zero());
    }

    #[serial]
    #[test]
    fn adjoin_sqrt_matches_algext_square_roots() {
        let k1 = k1_adjoin_sqrt2();
        let mut session = FieldSession::new(Arc::clone(&k1));
        let sqrt2 = AlgExtCPolyCoeff::from(
            AlgExtCData::from_alg_ext(&sqrt2_algext()).unwrap(),
        );
        let four = session.int(4).unwrap();
        let u = session.mul(&sqrt2, &four).unwrap();
        let beta = session.adjoin_sqrt(&u).unwrap();
        let sq = session.mul(&beta, &beta).unwrap();
        let (sq_a, u_a) = session.align(&sq, &u).unwrap();
        assert!(sq_a.coeff_sub(&u_a).unwrap().coeff_is_zero());
    }
}
