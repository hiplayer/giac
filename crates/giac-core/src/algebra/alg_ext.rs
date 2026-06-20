//! Algebraic extension elements (upstream `_EXT` / `ref_algext`).
//!
//! An [`AlgExtData`] is an element of K ≅ K[x]/(min_poly) where K is an [`ExtensionField`].
//! Complex algebraic numbers use [`AlgExtCData`](super::alg_ext_c::AlgExtCData) or legacy
//! `Expr::Complex(re, im)`.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-core-algebra-api-stability.md`.
//!
use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc, FuncKind};

use super::ext_tower::ExtensionField;
use super::field_arith::{
    canonical_poly1_expr, coords_to_expr, embed_in_cube_extension, embed_in_square_extension,
    generator_coords, min_poly_exprs_to_q, minpoly_at_cube, minpoly_at_square,
    mult_matrix_of_element, pad_to_len, poly1_coeffs,
    poly_degree, rationalize_poly1, ratio_to_expr_arc, CoordsQ,
};

/// Element of an algebraic extension field.
#[derive(Clone, Debug, PartialEq)]
pub struct AlgExtData {
    /// Ambient field K (tower-top).
    pub field: Arc<ExtensionField>,
    /// Element as polynomial in the primitive root, giac `poly1` order.
    pub coords: Vec<ExprArc>,
    pub root_index: Option<u32>,
}

impl AlgExtData {
    /// Backward-compatible `min_poly` view (defining polynomial over ℚ).
    /// **Stable** — `Poly::min_poly`
    pub fn min_poly(&self) -> Vec<ExprArc> {
        self.field.top_min_poly_exprs().unwrap_or_default()
    }

    /// **Stable** — Poly zero
    pub fn zero(field: Arc<ExtensionField>) -> Self {
        let coords = field.zero_coords();
        Self::from_coords_q(field, coords, None).unwrap()
    }

    /// **Stable** — Poly one
    pub fn one(field: Arc<ExtensionField>) -> Self {
        let coords = field.one_coords();
        Self::from_coords_q(field, coords, None).unwrap()
    }

    /// **Stable** — `from_field_coords`
    pub fn from_field_coords(field: Arc<ExtensionField>, coords: Vec<ExprArc>) -> Result<Self, EvalError> {
        let q = rationalize_poly1(&coords)?;
        let dim = field.dimension();
        Self::from_coords_q(field, pad_to_len(&q, dim), None)
    }

    // **Pipeline private** — `from_coords_q`
    fn from_coords_q(
        field: Arc<ExtensionField>,
        coords: CoordsQ,
        root_index: Option<u32>,
    ) -> Result<Self, EvalError> {
        let dim = field.dimension();
        // Coords must already be in this field's operational basis (primitive over ℚ when
        // parent is Base; tensor u^0..u^{e-1} blocks when parent is a nontrivial extension).
        // No mod reduction here — use field.element_* after construction.
        Ok(Self {
            field,
            coords: coords_to_expr(&pad_to_len(&coords, dim))?,
            root_index,
        })
    }

    /// **Stable** — `into_expr`
    pub fn into_expr(self) -> ExprArc {
        Arc::new(Expr::AlgExt(Arc::new(self)))
    }

    /// **Stable** — `from_rootof`
    pub fn from_rootof(num: &ExprArc, min_poly: &ExprArc) -> Result<Self, EvalError> {
        Self::from_rootof_over(num, min_poly, &ExtensionField::rational())
    }

    /// Like [`Self::from_rootof`] but adjoin over `parent` (default `parent = ℚ` in [`Self::from_rootof`]).
    /// **Stable** — `from_rootof_over`
    pub fn from_rootof_over(
        num: &ExprArc,
        min_poly: &ExprArc,
        parent: &Arc<ExtensionField>,
    ) -> Result<Self, EvalError> {
        let min_q = min_poly_exprs_to_q(&poly1_coeffs(min_poly)?)?;
        let field = ExtensionField::adjoin_irreducible(parent, min_q)?;
        let coords = match num.as_ref() {
            Expr::Seq(items) => rationalize_poly1(items)?,
            Expr::Int(_) | Expr::Rat(_) => rationalize_poly1(std::slice::from_ref(num))?,
            _ => return Err(EvalError::TypeError("rootof numerator")),
        };
        Self::from_coords_q(field, coords, None)
    }

    /// **Stable** — Poly addition
    pub fn add(&self, other: &Self) -> Result<Self, EvalError> {
        let (field, a, b) = align_pair(self, other)?;
        let sum = field.element_add(&a, &b)?;
        Self::from_coords_q(field, sum, None)
    }

    /// **Stable** — Poly subtraction
    pub fn sub(&self, other: &Self) -> Result<Self, EvalError> {
        let (field, a, b) = align_pair(self, other)?;
        let diff = field.element_sub(&a, &b)?;
        Self::from_coords_q(field, diff, None)
    }

    /// **Stable** — Poly multiplication
    pub fn mul(&self, other: &Self) -> Result<Self, EvalError> {
        let (field, a, b) = align_pair(self, other)?;
        let prod = field.element_mul(&a, &b)?;
        Self::from_coords_q(field, prod, None)
    }

    /// **Stable** — `inv`
    pub fn inv(&self) -> Result<Self, EvalError> {
        let a = self.coords_q()?;
        let inv = self.field.element_inv(&a)?;
        Self::from_coords_q(Arc::clone(&self.field), inv, self.root_index)
    }

    /// **Stable** — `eq_mod`
    pub fn eq_mod(&self, other: &Self) -> Result<bool, EvalError> {
        let (field, a, b) = align_pair(self, other)?;
        field.element_eq_mod(&a, &b)
    }

    /// **Stable** — Poly negation
    pub fn neg(&self) -> Result<Self, EvalError> {
        self.mul_rational(&Ratio::from_integer(-BigInt::one()))
    }

    /// **Stable** — `mul_rational`
    pub fn mul_rational(&self, r: &Ratio<BigInt>) -> Result<Self, EvalError> {
        if r.is_zero() {
            return Ok(Self::zero(Arc::clone(&self.field)));
        }
        if r.is_one() {
            return Ok(self.clone());
        }
        let coords = self.coords_q()?;
        let scaled = self
            .field
            .element_mul(&coords, &self.field.embed_rational(r))?;
        Self::from_coords_q(Arc::clone(&self.field), scaled, self.root_index)
    }

    /// **Stable** — Poly is zero
    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|c| c.is_zero())
    }

    /// **Stable** — Poly is one
    pub fn is_one(&self) -> bool {
        matches!(self.coords.as_slice(), [c] if c.is_one())
            && self.field.dimension() == 1
    }

    /// Build `rootof` display form from current coords + minpoly (L0).
    ///
    /// Does **not** call `common` / `align_elements` — display stays on the element's
    /// minimal ambient field. See [GIAC-lazy-common-tower-plan.md] §11.1 L0/L1.
    /// **Stable** — `to_rootof_expr`
    pub fn to_rootof_expr(&self) -> ExprArc {
        let coords = canonical_poly1_expr(&self.coords);
        let min_poly = canonical_poly1_expr(&self.min_poly());
        let num = Arc::new(Expr::Seq(coords));
        let min = if min_poly.len() == 1
            && matches!(min_poly[0].as_ref(), Expr::Func(FuncKind::Poly1, _))
        {
            Arc::clone(&min_poly[0])
        } else {
            Arc::new(Expr::Func(
                FuncKind::Poly1,
                vec![Arc::new(Expr::Seq(min_poly))],
            ))
        };
        Expr::func(FuncKind::RootOf, vec![num, min])
    }

    // **Pipeline private** — `coords_q`
    fn coords_q(&self) -> Result<CoordsQ, EvalError> {
        Ok(pad_to_len(
            &rationalize_poly1(&self.coords)?,
            self.field.dimension(),
        ))
    }
}

// **Pipeline private** — `fields_same`
fn fields_same(a: &Arc<ExtensionField>, b: &Arc<ExtensionField>) -> bool {
    Arc::ptr_eq(a, b) || **a == **b
}

// **Pipeline private** — `align_pair`
fn align_pair(a: &AlgExtData, b: &AlgExtData) -> Result<(Arc<ExtensionField>, CoordsQ, CoordsQ), EvalError> {
    let aligned = ExtensionField::align_elements(&a.field, &a.coords_q()?, &b.field, &b.coords_q()?)?;
    Ok((aligned.field, aligned.left, aligned.right))
}

/// How [`fold_algext_sum_mode`] merges unlike [`ExtensionField`] groups.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FoldAlgExtMode {
    /// Default: Expr traversal order; unlike fields may stay an `Add` tree (S1).
    #[default]
    Split,
    /// Sort groups by `field.id()` then pairwise `add` (S1-opt); order-independent.
    Canonical,
}

/// Fold a sum of `AlgExt` / rational leaves (eval `Add` patch).
///
/// **Lazy common (Split):** group by same [`ExtensionField`] (`Arc::ptr_eq` or `field ==`);
/// within-group merge uses [`AlgExtData::add`] without cross-domain `common`. When exactly
/// one group exists and a new term lies in a different field, that pairwise `add` may
/// `common`; otherwise unlike fields stay as an `Add` tree. See
/// [GIAC-lazy-common-tower-plan.md] §11.6 E.
/// **Stable** — canonical sum of AlgExt terms
pub fn fold_algext_sum(terms: &[ExprArc]) -> Result<ExprArc, EvalError> {
    fold_algext_sum_mode(terms, FoldAlgExtMode::Split)
}

/// Like [`fold_algext_sum`] with explicit merge mode (plan S1-opt).
/// **Stable** — fold_algext_sum with mode
pub fn fold_algext_sum_mode(terms: &[ExprArc], mode: FoldAlgExtMode) -> Result<ExprArc, EvalError> {
    let mut groups: Vec<(Arc<ExtensionField>, AlgExtData)> = Vec::new();
    let mut rat_sum = Ratio::<BigInt>::zero();
    let mut rest = Vec::new();
    for t in terms {
        match t.as_ref() {
            Expr::AlgExt(a) => {
                if let Some((_, acc)) = groups.iter_mut().find(|(f, _)| fields_same(f, &a.field)) {
                    *acc = acc.add(a)?;
                } else if groups.len() == 1 && mode == FoldAlgExtMode::Split {
                    let (f, acc) = groups.remove(0);
                    let sum = acc.add(a)?;
                    groups.push((Arc::clone(&sum.field), sum));
                    let _ = f;
                } else {
                    groups.push((Arc::clone(&a.field), (**a).clone()));
                }
            }
            Expr::Int(n) => rat_sum += Ratio::from_integer(n.clone()),
            Expr::Rat(r) => rat_sum += r.clone(),
            _ => rest.push(Arc::clone(t)),
        }
    }
    if mode == FoldAlgExtMode::Canonical && groups.len() > 1 {
        groups.sort_by_key(|(f, _)| f.id());
        let mut acc = groups.remove(0).1;
        for (_, next) in groups {
            acc = acc.add(&next)?;
        }
        groups = vec![(Arc::clone(&acc.field), acc)];
    }
    if !rat_sum.is_zero() {
        if groups.len() == 1 {
            let (field, acc) = &mut groups[0];
            let coords = acc.coords_q()?;
            let merged = field.element_add(&coords, &field.embed_rational(&rat_sum))?;
            *acc = AlgExtData::from_coords_q(Arc::clone(field), merged, acc.root_index)?;
        } else {
            rest.push(ratio_to_expr_arc(&rat_sum));
        }
    }
    for (_, a) in groups {
        if !a.is_zero() {
            rest.push(a.into_expr());
        }
    }
    match rest.len() {
        0 => Ok(Expr::int(0)),
        1 => Ok(Arc::clone(&rest[0])),
        _ => Ok(Expr::add(rest)),
    }
}

// **Pipeline private** — `complex_algext_parts`
fn complex_algext_parts(e: &Expr) -> Option<(ExprArc, ExprArc)> {
    match e {
        Expr::Complex(re, im) => Some((Arc::clone(re), Arc::clone(im))),
        Expr::AlgExt(_) | Expr::AlgExtC(_) => Some((Arc::new(e.clone()), Expr::int(0))),
        Expr::Int(n) => Some((Arc::new(Expr::Int(n.clone())), Expr::int(0))),
        Expr::Rat(r) => Some((Arc::new(Expr::Rat(r.clone())), Expr::int(0))),
        _ => None,
    }
}

// **Pipeline private** — `complex_algext_to_expr`
fn complex_algext_to_expr(re: ExprArc, im: ExprArc) -> Result<ExprArc, EvalError> {
    if im.is_zero() {
        Ok(re)
    } else if re.is_zero() {
        Ok(Arc::new(Expr::Complex(Expr::int(0), im)))
    } else {
        Ok(Arc::new(Expr::Complex(re, im)))
    }
}

/// **Stable** — sum with complex AlgExtC parts
pub fn fold_complex_algext_sum(terms: &[ExprArc]) -> Result<ExprArc, EvalError> {
    let mut res_re = Vec::new();
    let mut res_im = Vec::new();
    for t in terms {
        let (re, im) = complex_algext_parts(t.as_ref())
            .ok_or(EvalError::TypeError("not complex algext sum"))?;
        if !re.is_zero() {
            res_re.push(re);
        }
        if !im.is_zero() {
            res_im.push(im);
        }
    }
    let re_sum = match res_re.len() {
        0 => Expr::int(0),
        1 => res_re.remove(0),
        _ => fold_algext_sum(&res_re)?,
    };
    let im_sum = match res_im.len() {
        0 => Expr::int(0),
        1 => res_im.remove(0),
        _ => fold_algext_sum(&res_im)?,
    };
    complex_algext_to_expr(re_sum, im_sum)
}

/// **Stable** — product with complex AlgExtC parts
pub fn fold_complex_algext_product(factors: &[ExprArc]) -> Result<ExprArc, EvalError> {
    let mut acc_re = Expr::int(1);
    let mut acc_im = Expr::int(0);
    for f in factors {
        let (re, im) = complex_algext_parts(f.as_ref())
            .ok_or(EvalError::TypeError("not complex algext product"))?;
        let (new_re, new_im) = complex_algext_mul_parts(acc_re, acc_im, re, im)?;
        acc_re = new_re;
        acc_im = new_im;
    }
    complex_algext_to_expr(acc_re, acc_im)
}

// **Pipeline private** — `complex_algext_mul_parts`
fn complex_algext_mul_parts(
    ar: ExprArc,
    ai: ExprArc,
    br: ExprArc,
    bi: ExprArc,
) -> Result<(ExprArc, ExprArc), EvalError> {
    let ac = fold_algext_product(&[Arc::clone(&ar), Arc::clone(&br)])?;
    let bd = fold_algext_product(&[Arc::clone(&ai), bi.clone()])?;
    let ad = fold_algext_product(&[Arc::clone(&ar), bi.clone()])?;
    let bc = fold_algext_product(&[ai, br])?;
    let re = fold_algext_sum(&[ac, fold_algext_product(&[Expr::int(-1), bd])?])?;
    let im = fold_algext_sum(&[ad, bc])?;
    Ok((re, im))
}

/// **Stable** — canonical product of AlgExt terms
pub fn fold_algext_product(factors: &[ExprArc]) -> Result<ExprArc, EvalError> {
    let mut acc_ext: Option<AlgExtData> = None;
    let mut rat_prod = Ratio::<BigInt>::one();
    let mut rest = Vec::new();
    for f in factors {
        match f.as_ref() {
            Expr::AlgExt(a) => {
                acc_ext = Some(match acc_ext {
                    None => (**a).clone(),
                    Some(e) => e.mul(a)?,
                });
            }
            Expr::Int(n) => {
                if n.is_zero() {
                    return Ok(Expr::int(0));
                }
                rat_prod *= Ratio::from_integer(n.clone());
            }
            Expr::Rat(r) => {
                if r.is_zero() {
                    return Ok(Expr::int(0));
                }
                rat_prod *= r.clone();
            }
            _ => rest.push(Arc::clone(f)),
        }
    }
    if let Some(mut e) = acc_ext {
        if !rat_prod.is_one() {
            e = e.mul_rational(&rat_prod)?;
        }
        if !e.is_one() {
            rest.push(e.into_expr());
        }
    } else if !rat_prod.is_one() {
        rest.push(ratio_to_expr_arc(&rat_prod));
    }
    match rest.len() {
        0 => Ok(Expr::int(1)),
        1 => Ok(Arc::clone(&rest[0])),
        _ => Ok(Expr::mul(rest)),
    }
}

/// **Stable** — Func(RootOf) → AlgExt Expr
pub fn try_rootof_to_algext(args: &[ExprArc]) -> Result<ExprArc, EvalError> {
    if args.len() != 2 {
        return Err(EvalError::TooFewArgs("rootof"));
    }
    Ok(AlgExtData::from_rootof(&args[0], &args[1])?.into_expr())
}

/// **Stable** — subtree contains AlgExt or rootof
pub fn contains_algext(e: &Expr) -> bool {
    match e {
        Expr::AlgExt(_) | Expr::AlgExtC(_) => true,
        Expr::Add(terms) | Expr::Mul(terms) | Expr::Seq(terms) | Expr::List(terms) => {
            terms.iter().any(|t| contains_algext(t.as_ref()))
        }
        Expr::Frac(n, d) | Expr::Pow(n, d) | Expr::Complex(n, d) | Expr::Mod(n, d) => {
            contains_algext(n.as_ref()) || contains_algext(d.as_ref())
        }
        Expr::Func(_, args) => args.iter().any(|a| contains_algext(a.as_ref())),
        _ => false,
    }
}

/// **Stable** — view Expr as AlgExtData if present
pub fn try_as_algext_data(e: &Expr) -> Option<AlgExtData> {
    match e {
        Expr::AlgExt(a) => Some((**a).clone()),
        Expr::AlgExtC(z) if z.im.iter().all(|c| c.is_zero()) => {
            AlgExtData::from_field_coords(Arc::clone(&z.field), z.re.clone()).ok()
        }
        Expr::Func(FuncKind::RootOf, args) if args.len() == 2 => {
            AlgExtData::from_rootof(&args[0], &args[1]).ok()
        }
        _ => None,
    }
}

/// **Stable (bounded)** — square roots in extension field
pub fn algext_square_roots(u: &AlgExtData) -> Result<Vec<AlgExtData>, EvalError> {
    if u.is_zero() {
        return Ok(vec![AlgExtData::zero(Arc::clone(&u.field))]);
    }
    let f = u.coords_q()?;
    if u.field.dimension() == 1 {
        let u_val = pad_to_len(&f, 1)[0].clone();
        let min_s = vec![Ratio::one(), Ratio::zero(), -u_val];
        let sqrt_field = ExtensionField::adjoin_irreducible_over_q(min_s)?;
        let beta = AlgExtData::from_coords_q(
            Arc::clone(&sqrt_field),
            sqrt_field.generator_coords(),
            None,
        )?;
        return Ok(vec![beta.clone(), beta.neg()?]);
    }
    {
        let field = Arc::clone(&u.field);
        let one = field.one_coords();
        let zero = field.zero_coords();
        let neg = field.element_neg(&f)?;
        let sqrt_field =
            ExtensionField::adjoin_irreducible_parent_coeffs(&field, vec![one, zero, neg])?;
        let beta = AlgExtData::from_coords_q(
            Arc::clone(&sqrt_field),
            sqrt_field.generator_coords(),
            None,
        )?;
        return Ok(vec![beta.clone(), beta.neg()?]);
    }
}

/// One cube root of `u` in an adjoin extension (Cardano / resolvent helper).
/// **Stable (bounded)** — cube root in extension field
pub fn algext_cube_root(u: &AlgExtData) -> Result<AlgExtData, EvalError> {
    if u.is_zero() {
        return Ok(AlgExtData::zero(Arc::clone(&u.field)));
    }
    let f = u.coords_q()?;
    if u.field.dimension() == 1 {
        let u_val = pad_to_len(&f, 1)[0].clone();
        let min_c = vec![Ratio::one(), Ratio::zero(), Ratio::zero(), -u_val];
        let cbrt_field = ExtensionField::adjoin_irreducible_over_q(min_c)?;
        return AlgExtData::from_coords_q(
            Arc::clone(&cbrt_field),
            cbrt_field.generator_coords(),
            None,
        );
    }
    {
        let field = Arc::clone(&u.field);
        let one = field.one_coords();
        let zero = field.zero_coords();
        let neg = field.element_neg(&f)?;
        let cbrt_field = ExtensionField::adjoin_irreducible_parent_coeffs(
            &field,
            vec![one, zero.clone(), zero, neg],
        )?;
        return AlgExtData::from_coords_q(
            Arc::clone(&cbrt_field),
            cbrt_field.generator_coords(),
            None,
        );
    }
}

// **Pipeline private** — `cube_minpoly_via_matrix`
fn cube_minpoly_via_matrix(f: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> Result<CoordsQ, EvalError> {
    let n = poly_degree(m);
    let mat_u = mult_matrix_of_element(f, m);
    let mut block = vec![vec![Ratio::zero(); 3 * n]; 3 * n];
    for i in 0..n {
        block[i][n + i] = Ratio::one();
        block[n + i][2 * n + i] = Ratio::one();
        for j in 0..n {
            block[2 * n + i][j] = mat_u[i][j].clone();
        }
    }
    Ok(super::field_arith::char_poly_matrix(&block))
}

/// **Stable (bounded)** — sqrt branches as Expr list
pub fn algext_sqrt_branches(u: &AlgExtData) -> Result<Vec<ExprArc>, EvalError> {
    if let Ok(real) = algext_square_roots(u) {
        return Ok(real.into_iter().map(|a| a.into_expr()).collect());
    }
    let neg_u = u.neg()?;
    let betas = algext_square_roots(&neg_u)?;
    let beta = betas
        .into_iter()
        .next()
        .ok_or(EvalError::NotImplemented("alg ext square root"))?;
    let neg_beta = beta.neg()?;
    Ok(vec![
        Arc::new(Expr::Complex(Expr::int(0), beta.into_expr())),
        Arc::new(Expr::Complex(Expr::int(0), neg_beta.into_expr())),
    ])
}

/// Merge two extension elements into a common field (`common_EXT`).
/// **Stable** — common extension for two AlgExt values
pub fn common_ext(
    a: &AlgExtData,
    b: &AlgExtData,
) -> Result<(AlgExtData, AlgExtData, AlgExtData), EvalError> {
    let aligned = ExtensionField::align_elements(&a.field, &a.coords_q()?, &b.field, &b.coords_q()?)?;
    let gamma = AlgExtData::from_coords_q(
        Arc::clone(&aligned.field),
        aligned.field.generator_coords(),
        None,
    )?;
    let ae = AlgExtData::from_coords_q(
        Arc::clone(&aligned.field),
        aligned.left,
        a.root_index,
    )?;
    let be = AlgExtData::from_coords_q(
        Arc::clone(&aligned.field),
        aligned.right,
        b.root_index,
    )?;
    Ok((gamma, ae, be))
}

// **Pipeline private** — `sqrt_minpoly_via_matrix`
fn sqrt_minpoly_via_matrix(f: &[Ratio<BigInt>], m: &[Ratio<BigInt>]) -> Result<CoordsQ, EvalError> {
    let n = poly_degree(m);
    let mat_u = mult_matrix_of_element(f, m);
    let mut block = vec![vec![Ratio::zero(); 2 * n]; 2 * n];
    for i in 0..n {
        block[i][n + i] = Ratio::one();
        for j in 0..n {
            block[n + i][j] = mat_u[i][j].clone();
        }
    }
    Ok(super::field_arith::char_poly_matrix(&block))
}

#[cfg(test)]
mod tests {
    //! Test tiers: **A** / **B** — `.doc/test-writing-spec.md`
    //! Audit: `.doc/issues/GIAC-expr-api-test-audit.md` §6

    use super::*;
    use crate::algebra::test_fixtures::{
        algext_on_t1b_k2, algext_with_coords, cbrt2_algext, duplicate_field_arc,
        neg_sqrt2_algext, sqrt2_algext, sqrt3_algext, t1b_k2_adjoin_sqrt3_over_k1,
    };

    // **B** — AlgExt::mul; eq_mod α² = 2.
    #[test]
    fn algext_mul_squares_to_two() {
        let alpha = sqrt2_algext();
        let sq = alpha.mul(&alpha).unwrap();
        let two = AlgExtData::from_field_coords(
            Arc::clone(&alpha.field),
            vec![Expr::int(2)],
        )
        .unwrap();
        assert!(sq.eq_mod(&two).unwrap());
    }

    // **B** — AlgExt add cancellation.
    #[test]
    fn algext_add_neg_cancels() {
        let pos = sqrt2_algext();
        let neg = neg_sqrt2_algext();
        let z = pos.add(&neg).unwrap();
        assert!(z.is_zero());
    }

    // **B** — algext_square_roots; subfield embedding β² = √2.
    #[test]
    fn algext_sqrt_of_sqrt2() {
        let sqrt2 = sqrt2_algext();
        let k1 = Arc::clone(&sqrt2.field);
        let roots = algext_square_roots(&sqrt2).unwrap();
        assert_eq!(roots.len(), 2);
        let beta = &roots[0];
        let k4 = Arc::clone(&beta.field);
        let emb = ExtensionField::try_subfield_embedding(&k1, &k4)
            .unwrap()
            .expect("K1 embeds in sqrt adjoin");
        let sqrt2_in_k4 = emb.apply(&sqrt2.coords_q().unwrap());
        let beta_sq = beta.mul(beta).unwrap();
        assert!(
            k4.element_eq_mod(&beta_sq.coords_q().unwrap(), &sqrt2_in_k4).unwrap(),
            "beta^2 should equal sqrt2"
        );
    }

    // **B** — algext_cube_root; eq_mod β³ = 2.
    #[test]
    fn algext_cube_root_of_two() {
        let q = ExtensionField::rational();
        let two = AlgExtData::from_field_coords(Arc::clone(&q), vec![Expr::int(2)]).unwrap();
        let beta = algext_cube_root(&two).unwrap();
        let b3 = beta.mul(&beta).unwrap().mul(&beta).unwrap();
        assert!(b3.eq_mod(&two).unwrap(), "cbrt(2)^3 should equal 2");
    }

    // **B** — algext_sqrt_branches on -√2; eq_mod branch².
    #[test]
    fn algext_sqrt_of_neg_sqrt2_is_complex() {
        let neg_sqrt2 = neg_sqrt2_algext();
        let branches = algext_sqrt_branches(&neg_sqrt2).unwrap();
        assert_eq!(branches.len(), 2);
        for b in &branches {
            let data = try_as_algext_data(b.as_ref()).expect("sqrt branch as AlgExt");
            let sq = data.mul(&data).unwrap();
            assert!(sq.eq_mod(&neg_sqrt2).unwrap(), "branch^2 should be -sqrt(2)");
        }
    }

    // **B** — fold_complex_algext_sum cancellation.
    #[test]
    fn complex_algext_sum_cancels() {
        let alpha = sqrt2_algext();
        let pos = Arc::new(Expr::Complex(Expr::int(0), alpha.clone().into_expr()));
        let neg = Arc::new(Expr::Complex(
            Expr::int(0),
            alpha.neg().unwrap().into_expr(),
        ));
        let sum = fold_complex_algext_sum(&[pos, neg]).unwrap();
        assert!(sum.is_zero());
    }

    // **B** — to_rootof_expr / try_as_algext_data roundtrip; eq_mod.
    #[test]
    fn algext_to_rootof_roundtrip_display() {
        let alpha = sqrt2_algext();
        let e = alpha.to_rootof_expr();
        let back = try_as_algext_data(e.as_ref()).expect("rootof roundtrip");
        assert!(back.eq_mod(&alpha).unwrap());
    }

    // **B** — AlgExt::inv; α·α⁻¹ = 1.
    #[test]
    fn algext_inv_divides_to_one() {
        let alpha = sqrt2_algext();
        let inv = alpha.inv().unwrap();
        let one = alpha.mul(&inv).unwrap();
        assert!(one.is_one() || one.eq_mod(&AlgExtData::one(Arc::clone(&alpha.field))).unwrap());
    }

    // **B** — ℚ ↪ Q(√2) embed; dimension after add.
    #[test]
    fn subfield_embed_rational_into_sqrt2() {
        let sqrt2 = sqrt2_algext();
        let three = AlgExtData::from_field_coords(
            ExtensionField::rational(),
            vec![Expr::int(3)],
        )
        .unwrap();
        let sum = three.add(&sqrt2).unwrap();
        assert_eq!(sum.field.dimension(), 2);
    }

    // **B** — common_ext(√2,∛2); slow/ignored perf case.
    #[test]
    #[ignore = "primitive-element common(√2,∛2) char poly is slow; see ext_tower perf follow-up"]
    fn common_ext_sqrt2_cbrt2() {
        let sqrt2 = sqrt2_algext();
        let cbrt2 = cbrt2_algext();
        let (_gamma, ae, be) = common_ext(&sqrt2, &cbrt2).unwrap();
        let sum = ae.add(&be).unwrap();
        assert!(!sum.is_zero());
        assert_eq!(sum.field.dimension(), 6);
    }

    // **B** — common_over_q cache; add commutativity in compositum.
    #[test]
    fn algext_add_reverse_order_after_common_cache() {
        let sqrt2 = sqrt2_algext();
        let one = AlgExtData::from_field_coords(
            ExtensionField::rational(),
            vec![Expr::int(1)],
        )
        .unwrap();
        let _ = ExtensionField::common_over_q(&one.field, &sqrt2.field).unwrap();
        let sum_ab = one.add(&sqrt2).unwrap();
        let sum_ba = sqrt2.add(&one).unwrap();
        assert!(sum_ab.eq_mod(&sum_ba).unwrap());
        assert_eq!(sum_ab.field.id(), sqrt2.field.id());
        assert!(!sum_ab.is_zero());
    }

    // **B** — fold_algext_sum merges equal fields without Arc::ptr_eq.
    #[test]
    fn fold_algext_sum_merges_equal_fields_without_ptr_eq() {
        let alpha = sqrt2_algext();
        let dup_field = duplicate_field_arc(&alpha.field);
        assert!(!Arc::ptr_eq(&alpha.field, &dup_field));
        assert_eq!(*alpha.field, *dup_field);

        let beta = algext_with_coords(dup_field, alpha.coords_q().unwrap());
        let expected = alpha.add(&beta).unwrap();
        let folded = fold_algext_sum(&[alpha.into_expr(), beta.into_expr()]).unwrap();
        match folded.as_ref() {
            Expr::AlgExt(a) => assert!(a.eq_mod(&expected).unwrap()),
            other => panic!("expected single AlgExt, got {other:?}"),
        }
    }

    // **B** — from_coords_q roundtrip on K2 tensor fixture.
    #[test]
    fn from_coords_q_preserves_k2_tensor_coords_roundtrip() {
        let fix = t1b_k2_adjoin_sqrt3_over_k1();
        assert_eq!(fix.k2.dimension(), 4);
        for coords in [&fix.sqrt2_in_k2, &fix.beta_in_k2] {
            let sample = algext_on_t1b_k2(&fix, coords);
            let round = algext_with_coords(Arc::clone(&fix.k2), sample.coords_q().unwrap());
            assert_eq!(round.coords_q().unwrap(), sample.coords_q().unwrap());
            assert!(round.eq_mod(&sample).unwrap());
        }
    }

    // **B** — fold_algext_sum merges ℚ(5) into K2 via embed_rational.
    #[test]
    fn fold_algext_sum_rat_on_k2_merges_via_embed_rational() {
        let fix = t1b_k2_adjoin_sqrt3_over_k1();
        let sqrt2_k2 = algext_on_t1b_k2(&fix, &fix.sqrt2_in_k2);
        let five = Ratio::from_integer(5.into());
        let expected = sqrt2_k2
            .add(&algext_with_coords(
                Arc::clone(&fix.k2),
                fix.k2.embed_rational(&five),
            ))
            .unwrap();
        let folded =
            fold_algext_sum(&[sqrt2_k2.clone().into_expr(), Expr::int(5)]).unwrap();
        match folded.as_ref() {
            Expr::AlgExt(a) => assert!(a.eq_mod(&expected).unwrap()),
            other => panic!("expected single AlgExt, got {other:?}"),
        }
        let mut wrong_coords = sqrt2_k2.coords_q().unwrap();
        *wrong_coords.last_mut().unwrap() += five;
        let wrong = algext_on_t1b_k2(&fix, &wrong_coords);
        assert!(!wrong.eq_mod(&expected).unwrap());
    }

    // **B** — embed_rational block layout in K2.
    #[test]
    fn embed_rational_on_k2_places_constant_in_block_u0() {
        let fix = t1b_k2_adjoin_sqrt3_over_k1();
        assert_eq!(
            fix.k2.embed_rational(&Ratio::from_integer(5.into())),
            vec![
                Ratio::zero(),
                Ratio::from_integer(5.into()),
                Ratio::zero(),
                Ratio::zero(),
            ]
        );
    }

    // **B** — fold_algext_sum_mode Canonical order independence.
    #[test]
    fn fold_algext_sum_canonical_order_independent() {
        let sqrt2 = sqrt2_algext();
        let sqrt3 = sqrt3_algext();
        let terms_a = [sqrt2.clone().into_expr(), sqrt3.clone().into_expr()];
        let terms_b = [sqrt3.clone().into_expr(), sqrt2.clone().into_expr()];
        let folded_a = fold_algext_sum_mode(&terms_a, FoldAlgExtMode::Canonical).unwrap();
        let folded_b = fold_algext_sum_mode(&terms_b, FoldAlgExtMode::Canonical).unwrap();
        match (folded_a.as_ref(), folded_b.as_ref()) {
            (Expr::AlgExt(a), Expr::AlgExt(b)) => assert!(a.eq_mod(b).unwrap()),
            _ => panic!("expected single AlgExt from canonical fold"),
        }
    }

    // **B** — fold_algext_sum lazy common(√2,√3).
    #[test]
    fn fold_algext_sum_split_two_fields_merges_via_lazy_common() {
        let sqrt2 = sqrt2_algext();
        let sqrt3 = sqrt3_algext();
        let expected = sqrt2.add(&sqrt3).unwrap();
        let folded = fold_algext_sum(&[sqrt2.into_expr(), sqrt3.into_expr()]).unwrap();
        match folded.as_ref() {
            Expr::AlgExt(a) => assert!(a.eq_mod(&expected).unwrap()),
            other => panic!("expected merged AlgExt, got {other:?}"),
        }
    }

    // **B** — fold_algext_sum rat + AlgExt on K1.
    #[test]
    fn fold_algext_sum_rat_on_k1_still_merges_into_algext() {
        let alpha = sqrt2_algext();
        let expected = alpha
            .add(&algext_with_coords(
                Arc::clone(&alpha.field),
                alpha.field.embed_rational(&Ratio::from_integer(3.into())),
            ))
            .unwrap();
        let folded = fold_algext_sum(&[alpha.clone().into_expr(), Expr::int(3)]).unwrap();
        match folded.as_ref() {
            Expr::AlgExt(a) => assert!(a.eq_mod(&expected).unwrap()),
            other => panic!("expected single AlgExt, got {other:?}"),
        }
    }

    // **A** — eval(Expr::Frac) on AlgExt; division in extension field.
    #[test]
    fn algext_frac_via_eval() {
        use crate::{eval, Context};
        let ctx = Context::default();
        let alpha_data = sqrt2_algext();
        let alpha = alpha_data.clone().into_expr();
        let frac = Arc::new(Expr::Frac(Expr::int(2), alpha));
        let r = eval(frac.as_ref(), &ctx).unwrap();
        let quotient = match r.as_ref() {
            Expr::AlgExt(a) => (**a).clone(),
            other => panic!("expected AlgExt, got {other:?}"),
        };
        let prod = quotient.mul(&alpha_data).unwrap();
        let two = AlgExtData::from_field_coords(
            Arc::clone(&alpha_data.field),
            vec![Expr::int(2)],
        )
        .unwrap();
        assert!(prod.eq_mod(&two).unwrap());
    }
}
