//! Complex algebraic numbers z = re + im·i over an [`ExtensionField`](super::ext_tower::ExtensionField).
//!
//! `i` is formal with i² = −1; it is **not** adjoined to the tower until explicitly
//! required (Phase 2b+). Normative: [GIAC-algext-adoption.md](../../../../.doc/issues/GIAC-algext-adoption.md) §8.3.

use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc, FuncKind};

use super::alg_ext::AlgExtData;
use super::ext_tower::ExtensionField;
use super::field_arith::{
    coords_to_expr, min_poly_exprs_to_q, pad_to_len, poly1_coeffs, rationalize_poly1,
    ratio_to_expr_arc, CoordsQ,
};

/// z ∈ K[i]/(i²+1) where K is the tower-top [`ExtensionField`].
#[derive(Clone, Debug, PartialEq)]
pub struct AlgExtCData {
    pub field: Arc<ExtensionField>,
    /// Real part coords in K (length = field.dimension()).
    pub re: Vec<ExprArc>,
    /// Imaginary part coords in K.
    pub im: Vec<ExprArc>,
    pub root_index: Option<u32>,
}

impl AlgExtCData {
    pub fn zero(field: Arc<ExtensionField>) -> Result<Self, EvalError> {
        let z = field.zero_coords();
        Self::from_coords_q(field, &z, &z, None)
    }

    pub fn one(field: Arc<ExtensionField>) -> Result<Self, EvalError> {
        let re = field.one_coords();
        let im = field.zero_coords();
        Self::from_coords_q(field, &re, &im, None)
    }

    pub fn from_alg_ext(a: &AlgExtData) -> Result<Self, EvalError> {
        let re = rationalize_poly1(&a.coords)?;
        let im = a.field.zero_coords();
        Self::from_coords_q(Arc::clone(&a.field), &re, &im, a.root_index)
    }

    pub fn from_complex_parts(re: &ExprArc, im: &ExprArc) -> Result<Self, EvalError> {
        let (re_field, re_q) = expr_to_field_element(re)?;
        let (im_field, im_q) = expr_to_field_element(im)?;
        let aligned = ExtensionField::align_elements(&re_field, &re_q, &im_field, &im_q)?;
        Self::from_coords_q(
            aligned.field,
            &aligned.left,
            &aligned.right,
            None,
        )
    }

    fn from_coords_q(
        field: Arc<ExtensionField>,
        re: &CoordsQ,
        im: &CoordsQ,
        root_index: Option<u32>,
    ) -> Result<Self, EvalError> {
        let dim = field.dimension();
        Ok(Self {
            field,
            re: coords_to_expr(&pad_to_len(re, dim))?,
            im: coords_to_expr(&pad_to_len(im, dim))?,
            root_index,
        })
    }

    fn re_q(&self) -> Result<CoordsQ, EvalError> {
        Ok(pad_to_len(
            &rationalize_poly1(&self.re)?,
            self.field.dimension(),
        ))
    }

    fn im_q(&self) -> Result<CoordsQ, EvalError> {
        Ok(pad_to_len(
            &rationalize_poly1(&self.im)?,
            self.field.dimension(),
        ))
    }

    fn ensure_same_field(&self, other: &Self) -> Result<(), EvalError> {
        if !Arc::ptr_eq(&self.field, &other.field) && self.field != other.field {
            return Err(EvalError::TypeError("AlgExtC field mismatch"));
        }
        Ok(())
    }

    pub fn add(&self, other: &Self) -> Result<Self, EvalError> {
        let pair = self.align_with(other)?;
        let re = pair
            .field
            .element_add(&pair.self_re, &pair.other_re)?;
        let im = pair
            .field
            .element_add(&pair.self_im, &pair.other_im)?;
        Self::from_coords_q(Arc::clone(&pair.field), &re, &im, None)
    }

    pub fn sub(&self, other: &Self) -> Result<Self, EvalError> {
        let pair = self.align_with(other)?;
        let re = pair
            .field
            .element_sub(&pair.self_re, &pair.other_re)?;
        let im = pair
            .field
            .element_sub(&pair.self_im, &pair.other_im)?;
        Self::from_coords_q(Arc::clone(&pair.field), &re, &im, None)
    }

    pub fn mul(&self, other: &Self) -> Result<Self, EvalError> {
        let pair = self.align_with(other)?;
        let f = &pair.field;
        let ac = f.element_mul(&pair.self_re, &pair.other_re)?;
        let bd = f.element_mul(&pair.self_im, &pair.other_im)?;
        let ad = f.element_mul(&pair.self_re, &pair.other_im)?;
        let bc = f.element_mul(&pair.self_im, &pair.other_re)?;
        let re = f.element_sub(&ac, &bd)?;
        let im = f.element_add(&ad, &bc)?;
        Self::from_coords_q(Arc::clone(f), &re, &im, None)
    }

    pub fn neg(&self) -> Result<Self, EvalError> {
        let re = self.field.element_neg(&self.re_q()?)?;
        let im = self.field.element_neg(&self.im_q()?)?;
        Self::from_coords_q(Arc::clone(&self.field), &re, &im, self.root_index)
    }

    /// Multiplicative inverse using N(z) = re² + im² ∈ K.
    pub fn inv(&self) -> Result<Self, EvalError> {
        let re = self.re_q()?;
        let im = self.im_q()?;
        let re2 = self.field.element_mul(&re, &re)?;
        let im2 = self.field.element_mul(&im, &im)?;
        let norm = self.field.element_add(&re2, &im2)?;
        let norm_inv = self.field.element_inv(&norm)?;
        let re_part = self.field.element_mul(&re, &norm_inv)?;
        let im_part = self.field.element_mul(&im, &norm_inv)?;
        let im_neg = self.field.element_neg(&im_part)?;
        Self::from_coords_q(Arc::clone(&self.field), &re_part, &im_neg, self.root_index)
    }

    pub fn eq_mod(&self, other: &Self) -> Result<bool, EvalError> {
        let pair = self.align_with(other)?;
        Ok(pair.field.element_eq_mod(&pair.self_re, &pair.other_re)?
            && pair.field.element_eq_mod(&pair.self_im, &pair.other_im)?)
    }

    pub fn is_zero(&self) -> Result<bool, EvalError> {
        Ok(self.field.element_is_zero(&self.re_q()?)
            && self.field.element_is_zero(&self.im_q()?))
    }

    pub fn is_one(&self) -> Result<bool, EvalError> {
        Ok(self.field.element_eq_mod(&self.re_q()?, &self.field.one_coords())?
            && self.field.element_is_zero(&self.im_q()?))
    }

    /// Convert to legacy `Expr` shapes for display / eval compatibility.
    pub fn into_expr(self) -> ExprArc {
        Arc::new(self.to_expr())
    }

    pub fn to_expr(&self) -> Expr {
        if self.im.iter().all(|c| c.is_zero()) {
            if self.re.iter().all(|c| c.is_zero()) {
                return Expr::Int(0.into());
            }
            return self.re_to_legacy_expr();
        }
        if self.re.iter().all(|c| c.is_zero()) {
            return Expr::Complex(Expr::int(0), self.im_to_legacy_expr());
        }
        Expr::AlgExtC(Arc::new(self.clone()))
    }

    fn re_to_legacy_expr(&self) -> Expr {
        match AlgExtData::from_field_coords(Arc::clone(&self.field), self.re.clone()) {
            Ok(a) => Expr::AlgExt(Arc::new(a)),
            Err(_) => Expr::Add(self.re.clone()),
        }
    }

    fn im_to_legacy_expr(&self) -> ExprArc {
        match AlgExtData::from_field_coords(Arc::clone(&self.field), self.im.clone()) {
            Ok(a) => a.into_expr(),
            Err(_) => Arc::new(Expr::Add(self.im.clone())),
        }
    }

    fn align_with(&self, other: &Self) -> Result<AlignedPair, EvalError> {
        if Arc::ptr_eq(&self.field, &other.field) || self.field == other.field {
            return Ok(AlignedPair {
                field: Arc::clone(&self.field),
                self_re: self.re_q()?,
                self_im: self.im_q()?,
                other_re: other.re_q()?,
                other_im: other.im_q()?,
            });
        }
        let aligned_re =
            ExtensionField::align_elements(&self.field, &self.re_q()?, &other.field, &other.re_q()?)?;
        let aligned_im =
            ExtensionField::align_elements(&self.field, &self.im_q()?, &other.field, &other.im_q()?)?;
        debug_assert_eq!(aligned_re.field.id(), aligned_im.field.id());
        Ok(AlignedPair {
            field: aligned_re.field,
            self_re: aligned_re.left,
            self_im: aligned_im.left,
            other_re: aligned_re.right,
            other_im: aligned_im.right,
        })
    }
}

struct AlignedPair {
    field: Arc<ExtensionField>,
    self_re: CoordsQ,
    self_im: CoordsQ,
    other_re: CoordsQ,
    other_im: CoordsQ,
}

/// Canonicalize `Expr` leaves to [`AlgExtCData`].
pub fn canonicalize_to_algext_c(e: &Expr) -> Result<AlgExtCData, EvalError> {
    match e {
        Expr::AlgExtC(z) => Ok((**z).clone()),
        Expr::AlgExt(a) => AlgExtCData::from_alg_ext(a),
        Expr::Complex(re, im) => AlgExtCData::from_complex_parts(re, im),
        Expr::Int(n) => {
            let q = ExtensionField::rational();
            let mut re = q.zero_coords();
            if let Some(c) = re.last_mut() {
                *c = Ratio::from_integer(n.clone());
            }
            let im = q.zero_coords();
            AlgExtCData::from_coords_q(q, &re, &im, None)
        }
        Expr::Rat(r) => {
            let q = ExtensionField::rational();
            let mut re = q.zero_coords();
            if let Some(c) = re.last_mut() {
                *c = r.clone();
            }
            let im = q.zero_coords();
            AlgExtCData::from_coords_q(q, &re, &im, None)
        }
        Expr::Func(FuncKind::RootOf, args) if args.len() == 2 => {
            let a = AlgExtData::from_rootof(&args[0], &args[1])?;
            AlgExtCData::from_alg_ext(&a)
        }
        _ => Err(EvalError::TypeError("cannot canonicalize to AlgExtC")),
    }
}

fn expr_to_field_element(e: &ExprArc) -> Result<(Arc<ExtensionField>, CoordsQ), EvalError> {
    match e.as_ref() {
        Expr::AlgExt(a) => {
            let q = rationalize_poly1(&a.coords)?;
            Ok((Arc::clone(&a.field), q))
        }
        Expr::AlgExtC(z) => Ok((Arc::clone(&z.field), z.re_q()?)),
        Expr::Int(n) => {
            let q = ExtensionField::rational();
            let mut re = q.zero_coords();
            if let Some(c) = re.last_mut() {
                *c = Ratio::from_integer(n.clone());
            }
            Ok((q, re))
        }
        Expr::Rat(r) => {
            let q = ExtensionField::rational();
            let mut re = q.zero_coords();
            if let Some(c) = re.last_mut() {
                *c = r.clone();
            }
            Ok((q, re))
        }
        Expr::Func(FuncKind::RootOf, args) if args.len() == 2 => {
            let a = AlgExtData::from_rootof(&args[0], &args[1])?;
            let q = rationalize_poly1(&a.coords)?;
            Ok((Arc::clone(&a.field), q))
        }
        _ => Err(EvalError::TypeError("field element expected")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Expr;

    fn sqrt2_algext() -> AlgExtData {
        let min = Arc::new(Expr::Func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                Expr::int(0),
                Expr::int(-2),
            ]))],
        ));
        AlgExtData::from_rootof(
            &Arc::new(Expr::Seq(vec![Expr::int(1), Expr::int(0)])),
            &min,
        )
        .unwrap()
    }

    #[test]
    fn algext_c_from_algext_is_real() {
        let a = sqrt2_algext();
        let z = AlgExtCData::from_alg_ext(&a).unwrap();
        assert!(z.im.iter().all(|c| c.is_zero()));
        assert_eq!(z.field.dimension(), 2);
    }

    #[test]
    fn i_times_sqrt2_squared_is_minus_two() {
        let a = sqrt2_algext();
        let sqrt2 = AlgExtCData::from_alg_ext(&a).unwrap();
        let i_sqrt2 = AlgExtCData::from_complex_parts(&Expr::int(0), &a.into_expr()).unwrap();
        let sq = i_sqrt2.mul(&i_sqrt2).unwrap();
        let minus_two = AlgExtCData::from_complex_parts(
            &ratio_to_expr_arc(&Ratio::from_integer((-2).into())),
            &Expr::int(0),
        )
        .unwrap();
        // Align fields: sqrt2 lives in ℚ(√2), −2 in ℚ — common merges them.
        assert!(sq.eq_mod(&minus_two).unwrap());
    }

    #[test]
    fn canonicalize_algext_roundtrip() {
        let a = sqrt2_algext();
        let z = canonicalize_to_algext_c(&Expr::AlgExt(Arc::new(a.clone()))).unwrap();
        let back = z.to_expr();
        assert!(matches!(back, Expr::AlgExt(_)));
    }

    #[test]
    fn canonicalize_complex_with_algext_im() {
        let a = sqrt2_algext();
        let e = Expr::Complex(Expr::int(0), a.into_expr());
        let z = canonicalize_to_algext_c(&e).unwrap();
        assert!(!z.im.iter().all(|c| c.is_zero()));
    }
}
