//! F4′ — Galois automorphisms on splitting fields.
//!
//! Used by `poly_roots::sqrt_in_field_euler_second` to compute √β = σ(√α) without
//! a second blind adjoin when β = σ(α) for σ ∈ Gal(P/ℚ).

use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::EvalError;

use super::ext_tower::{self, ExtensionField};
use super::field_arith::{pad_to_len, CoordsQ};

/// **Pipeline private** — solve `A·c = b` over ℚ; `A` is `rows × cols`, returns `c` or `None`.
fn solve_linear_system(
    rows: &[Vec<Ratio<BigInt>>],
    b: &[Ratio<BigInt>],
) -> Option<Vec<Ratio<BigInt>>> {
    let n = rows.len();
    if n == 0 {
        return None;
    }
    let cols = rows[0].len();
    let mut aug = vec![vec![Ratio::zero(); cols + 1]; n];
    for i in 0..n {
        for j in 0..cols {
            aug[i][j] = rows[i][j].clone();
        }
        aug[i][cols] = pad_to_len(b, n)[i].clone();
    }
    let (pivot_cols, inconsistent) = ext_tower::gauss_elim_rref(&mut aug);
    if inconsistent {
        return None;
    }
    let mut c = vec![Ratio::zero(); cols];
    for (row, &pc) in pivot_cols.iter().enumerate() {
        if pc < cols {
            c[pc] = aug[row][cols].clone();
        }
    }
    Some(c)
}

/// **Pipeline private** — embed parent coords into child (lower block).
fn embed_parent_coords(
    field: &Arc<ExtensionField>,
    parent: &Arc<ExtensionField>,
    a: &CoordsQ,
) -> Result<CoordsQ, EvalError> {
    let pd = parent.dimension();
    let cd = field.dimension();
    if cd == pd {
        return Ok(pad_to_len(a, pd));
    }
    let mut out = field.zero_coords();
    let a = pad_to_len(a, pd);
    #[allow(clippy::manual_memcpy)] // Ratio<BigInt> is not Copy
    for i in 0..pd.min(cd) {
        out[i] = a[i].clone();
    }
    Ok(out)
}

/// **Pipeline private** — if `x = Σ c_k src^k` in `field`, return `Σ c_k dst^k`.
fn try_conjugate_by_power_basis(
    field: &Arc<ExtensionField>,
    x: &CoordsQ,
    src: &CoordsQ,
    dst: &CoordsQ,
    max_deg: usize,
) -> Result<Option<CoordsQ>, EvalError> {
    let n = field.dimension();
    let x = pad_to_len(x, n);
    let src = pad_to_len(src, n);
    let dst = pad_to_len(dst, n);

    let mut pows = vec![field.one_coords()];
    let mut cur = src.clone();
    for _ in 0..max_deg {
        pows.push(cur.clone());
        cur = field.element_mul(&cur, &src)?;
    }
    let cols = pows.len();
    let mut rows = vec![vec![Ratio::zero(); cols]; n];
    for i in 0..n {
        for k in 0..cols {
            rows[i][k] = pad_to_len(&pows[k], n)[i].clone();
        }
    }
    let Some(coeffs) = solve_linear_system(&rows, &x) else {
        return Ok(None);
    };

    let mut out = field.zero_coords();
    let mut dst_pow = field.one_coords();
    for ck in coeffs {
        if !ck.is_zero() {
            let term = field.element_mul(&field.embed_rational(&ck), &dst_pow)?;
            out = field.element_add(&out, &term)?;
        }
        dst_pow = field.element_mul(&dst_pow, &dst)?;
    }
    Ok(Some(out))
}

/// **Pipeline private** — σ(x) via conjugation in a cubic (or smaller) subfield, embed back.
fn try_apply_sigma_on_parent(
    parent: &Arc<ExtensionField>,
    x: &CoordsQ,
    src: &CoordsQ,
    dst: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let mut cur: Option<Arc<ExtensionField>> = Some(Arc::clone(parent));
    while let Some(sub) = cur {
        let emb = if Arc::ptr_eq(&sub, parent) {
            None
        } else {
            ExtensionField::try_subfield_embedding(&sub, parent)?
        };
        let (xs, ss, sd) = match &emb {
            Some(e) => (
                ext_tower::try_preimage_under_embedding(e, x),
                ext_tower::try_preimage_under_embedding(e, src),
                ext_tower::try_preimage_under_embedding(e, dst),
            ),
            None => (Some(x.clone()), Some(src.clone()), Some(dst.clone())),
        };
        if let (Some(xs), Some(ss), Some(sd)) = (xs, ss, sd) {
            let deg = sub.dimension();
            if let Some(img) = try_conjugate_by_power_basis(&sub, &xs, &ss, &sd, deg)? {
                return Ok(match &emb {
                    Some(e) => Some(e.apply(&img)),
                    None => Some(img),
                });
            }
        }
        cur = sub.parent_field().map(Arc::clone);
    }
    Ok(None)
}

/// **Pipeline private** — verify σ(src)=dst and σ(u)=v.
pub(crate) fn conjugate_map_sends(
    field: &Arc<ExtensionField>,
    src: &CoordsQ,
    dst: &CoordsQ,
    u: &CoordsQ,
    v: &CoordsQ,
    _z_roots: &[CoordsQ],
) -> Result<bool, EvalError> {
    let Some(img_src) = try_apply_sigma_on_parent(field, src, src, dst)? else {
        return Ok(false);
    };
    if !field.element_eq_mod(&img_src, dst)? {
        return Ok(false);
    }
    let Some(img_u) = try_apply_sigma_on_parent(field, u, src, dst)? else {
        return Ok(false);
    };
    field.element_eq_mod(&img_u, v)
}

fn coords_square_eq_mod(
    field: &ExtensionField,
    cand: &CoordsQ,
    u: &CoordsQ,
) -> Result<bool, EvalError> {
    let sq = field.element_mul(cand, cand)?;
    field.element_eq_mod(&sq, u)
}

/// **Pipeline private** — σ(κ) with κ²=u, σ(u)=v, σ|_P from subfield conjugation.
fn try_sigma_kappa(
    field: &Arc<ExtensionField>,
    parent: &Arc<ExtensionField>,
    src: &CoordsQ,
    dst: &CoordsQ,
    v_emb: &CoordsQ,
    kappa: &CoordsQ,
) -> Result<Option<CoordsQ>, EvalError> {
    let pd = parent.dimension();
    let try_pair = |a: &CoordsQ, b: &CoordsQ| -> Result<Option<CoordsQ>, EvalError> {
        let a_emb = embed_parent_coords(field, parent, a)?;
        let b_emb = embed_parent_coords(field, parent, b)?;
        let h = field.element_add(&a_emb, &field.element_mul(&b_emb, kappa)?)?;
        for cand in [h.clone(), field.element_neg(&h)?] {
            if coords_square_eq_mod(field, &cand, v_emb)? {
                return Ok(Some(cand));
            }
        }
        Ok(None)
    };

    let zero = parent.zero_coords();
    let one = parent.one_coords();
    let mut sparse = vec![zero.clone(), one.clone()];
    for i in 0..pd {
        let mut e = parent.zero_coords();
        e[i] = Ratio::one();
        sparse.push(e.clone());
        let mut ne = parent.zero_coords();
        ne[i] = Ratio::from_integer((-1).into());
        sparse.push(ne);
    }

    for s in &sparse {
        let a = try_apply_sigma_on_parent(parent, s, src, dst)?
            .unwrap_or_else(|| parent.zero_coords());
        for t in &sparse {
            let b = try_apply_sigma_on_parent(parent, t, src, dst)?
                .unwrap_or_else(|| parent.zero_coords());
            if let Some(hit) = try_pair(&a, &b)? {
                return Ok(Some(hit));
            }
        }
    }
    Ok(None)
}

/// **Pipeline private** — F4′: √v = σ(√u) when σ(α)=β on resolvent splitting field.
pub(crate) fn try_galois_sqrt_second(
    field: &Arc<ExtensionField>,
    u: &CoordsQ,
    v: &CoordsQ,
    sqrt_u: &CoordsQ,
    z_roots: &[CoordsQ],
) -> Result<Option<CoordsQ>, EvalError> {
    let parent = match field.parent_field() {
        Some(p) => Arc::clone(p),
        None => return Ok(None),
    };
    let pd = parent.dimension();
    if field.dimension() != 2 * pd {
        return Ok(None);
    }
    let u = pad_to_len(u, pd);
    let v = pad_to_len(v, pd);
    let sqrt_u = pad_to_len(sqrt_u, field.dimension());
    let kappa = field.generator_coords();
    let v_emb = embed_parent_coords(field, &parent, &v)?;

    for src in z_roots {
        let src = pad_to_len(src, pd);
        if !parent.element_eq_mod(&src, &u)? {
            continue;
        }
        for dst in z_roots {
            let dst = pad_to_len(dst, pd);
            if !parent.element_eq_mod(&dst, &v)? {
                continue;
            }
            if !conjugate_map_sends(&parent, &src, &dst, &u, &v, z_roots)? {
                continue;
            }

            if let Some(h) = ext_tower::ExtensionField::try_square_root_in_field(field, &v_emb)? {
                return Ok(Some(h));
            }

            if let Some(h) = try_sigma_kappa(field, &parent, &src, &dst, &v_emb, &kappa)? {
                return Ok(Some(h));
            }

            if let Ok(inv_u) = parent.element_inv(&u) {
                let ratio = parent.element_mul(&v, &inv_u)?;
                if let Some(sigma_ratio) = try_apply_sigma_on_parent(&parent, &ratio, &src, &dst)? {
                    let ratio_emb = embed_parent_coords(field, &parent, &sigma_ratio)?;
                    let cand = field.element_mul(&sqrt_u, &ratio_emb)?;
                    for h in [cand.clone(), field.element_neg(&cand)?] {
                        if coords_square_eq_mod(field, &h, &v_emb)? {
                            return Ok(Some(h));
                        }
                    }
                }
            }

            // ponytail: σ(√u) = √u · σ(zj/zi) for resolvent root ratios (S₃).
            for zi in z_roots {
                let zi = pad_to_len(zi, pd);
                let Ok(inv_zi) = parent.element_inv(&zi) else {
                    continue;
                };
                for zj in z_roots {
                    let zj = pad_to_len(zj, pd);
                    let r = parent.element_mul(&zj, &inv_zi)?;
                    if let Some(sr) = try_apply_sigma_on_parent(&parent, &r, &src, &dst)? {
                        let r_emb = embed_parent_coords(field, &parent, &sr)?;
                        let cand = field.element_mul(&sqrt_u, &r_emb)?;
                        for h in [cand.clone(), field.element_neg(&cand)?] {
                            if coords_square_eq_mod(field, &h, &v_emb)? {
                                return Ok(Some(h));
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(None)
}
