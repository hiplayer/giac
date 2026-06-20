//! Dense poly1 ↔ sparse univariate conversion (ℚ coefficients).
//!
//! Locks giac **`poly1` HighFirst** vs giac-poly **ascending** layout per
//! [GIAC-dense-poly1-refactor](.doc/issues/GIAC-dense-poly1-refactor.md) §4.2 and
//! [GIAC-lazy-common-tower-plan](.doc/issues/GIAC-lazy-common-tower-plan.md) §12.9.
//!
//! **Tier:** Stable (crate-internal).

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Zero;

use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::univariate_coeffs_ascending;

/// **Stable (crate-internal)** — reverse coefficient order (HighFirst ↔ Ascending).
pub fn reverse_coeffs(c: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>> {
    c.iter().rev().cloned().collect()
}

/// **Stable (crate-internal)** — ascending dense → giac `poly1` HighFirst.
pub fn ascending_to_dense_high_first(asc: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>> {
    reverse_coeffs(asc)
}

/// **Stable (crate-internal)** — giac `poly1` HighFirst → ascending dense.
pub fn dense_high_first_to_ascending(hi: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>> {
    reverse_coeffs(hi)
}

/// **Stable (crate-internal)** — sparse `Poly` univariate in `var` → dense HighFirst.
pub fn sparse_ascending_to_dense_high_first(p: &Poly, var: &Var) -> Vec<Ratio<BigInt>> {
    ascending_to_dense_high_first(&univariate_coeffs_ascending(p, var))
}

/// **Stable (crate-internal)** — dense HighFirst → sparse univariate in `var`.
pub fn dense_high_first_to_sparse(p: &[Ratio<BigInt>], var: &Var) -> Poly {
    poly_from_ascending_coeffs(var, &dense_high_first_to_ascending(p))
}

fn poly_from_ascending_coeffs(var: &Var, coeffs: &[Ratio<BigInt>]) -> Poly {
    let mut out = Poly::zero();
    for (e, c) in coeffs.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let term = if e == 0 {
            Poly::constant(c.clone())
        } else {
            Poly::constant(c.clone()).mul(&Poly::var(var.clone()).pow(e as u64))
        };
        out = out.add(&term);
    }
    out
}
