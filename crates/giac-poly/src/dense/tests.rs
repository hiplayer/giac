//! Dense poly1 unit tests — locks **HighFirst** semantics (giac `poly1`).
//!
//! Cases extracted from `giac-core::algebra::field_arith` behaviour, not ported tests
//! (core has no `#[cfg(test)]` on that module). See
//! [GIAC-dense-poly1-refactor](.doc/issues/GIAC-dense-poly1-refactor.md) §7.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_bigint::BigInt;
use num_rational::Ratio;

use super::poly1::{
    Poly1Order, add, div_rem, ext_gcd, inv_mod, mul, neg, poly_degree, reduce_mod_monic, scale, sub,
    trim,
};
use super::ratio_ring::RatioRingCtx;

// **Pipeline private** — `q`
fn q(n: i64) -> Ratio<BigInt> {
    Ratio::from_integer(BigInt::from(n))
}

// **Pipeline private** — `ctx`
fn ctx() -> RatioRingCtx {
    RatioRingCtx
}

#[test]
// **Pipeline private** — `high_first_degree_and_trim`
fn high_first_degree_and_trim() {
    let c = ctx();
    // x^2 - 2
    let mut p = vec![q(1), q(0), q(-2)];
    assert_eq!(poly_degree(&p), 2);
    trim(&c, &mut p, Poly1Order::HighFirst);
    assert_eq!(p, vec![q(1), q(0), q(-2)]);

    // 0*x^2 + x + 1 → x + 1
    let mut qpoly = vec![q(0), q(1), q(1)];
    trim(&c, &mut qpoly, Poly1Order::HighFirst);
    assert_eq!(qpoly, vec![q(1), q(1)]);

    // all zero → [0]
    let mut zero = vec![q(0), q(0)];
    trim(&c, &mut zero, Poly1Order::HighFirst);
    assert_eq!(zero, vec![q(0)]);
}

#[test]
// **Pipeline private** — `high_first_add_sub`
fn high_first_add_sub() {
    let c = ctx();
    // (x+1) + (x+2) = 2x + 3 → [2, 3]
    let a = vec![q(1), q(1)];
    let b = vec![q(1), q(2)];
    assert_eq!(add(&c, &a, &b, Poly1Order::HighFirst).unwrap(), vec![q(2), q(3)]);
    // (x+1)-(x+2) = -1 → trim drops explicit zero high term
    assert_eq!(sub(&c, &a, &b, Poly1Order::HighFirst).unwrap(), vec![q(-1)]);
}

#[test]
// **Pipeline private** — `high_first_mul_x_plus_1_times_x_plus_2`
fn high_first_mul_x_plus_1_times_x_plus_2() {
    let c = ctx();
    // (x+1)(x+2) = x^2 + 3x + 2
    let a = vec![q(1), q(1)];
    let b = vec![q(1), q(2)];
    assert_eq!(
        mul(&c, &a, &b, Poly1Order::HighFirst).unwrap(),
        vec![q(1), q(3), q(2)]
    );
}

#[test]
// **Pipeline private** — `high_first_neg_and_scale`
fn high_first_neg_and_scale() {
    let c = ctx();
    let p = vec![q(1), q(0), q(-2)];
    assert_eq!(neg(&c, &p, Poly1Order::HighFirst).unwrap(), vec![q(-1), q(0), q(2)]);
    assert_eq!(scale(&c, &p, &q(2), Poly1Order::HighFirst).unwrap(), vec![q(2), q(0), q(-4)]);
    assert_eq!(scale(&c, &p, &q(0), Poly1Order::HighFirst).unwrap(), vec![q(0)]);
}

#[test]
// **Pipeline private** — `high_first_div_rem`
fn high_first_div_rem() {
    let c = ctx();
    // x^2 + 3x + 2 = (x+1) * (x+2) + 0
    let dividend = vec![q(1), q(3), q(2)];
    let divisor = vec![q(1), q(1)];
    let (quot, rem) = div_rem(&c, &dividend, &divisor, Poly1Order::HighFirst).unwrap();
    assert_eq!(quot, vec![q(1), q(2)]);
    assert_eq!(rem, vec![q(0)]);

    // x^2 + 1 = (x+1)(x-1) + 2
    let a = vec![q(1), q(0), q(1)];
    let b = vec![q(1), q(1)];
    let (q_out, r) = div_rem(&c, &a, &b, Poly1Order::HighFirst).unwrap();
    assert_eq!(q_out, vec![q(1), q(-1)]);
    assert_eq!(r, vec![q(2)]);
}

#[test]
// **Pipeline private** — `high_first_reduce_mod_x_squared_minus_2`
fn high_first_reduce_mod_x_squared_minus_2() {
    let c = ctx();
    // x^3 mod (x^2 - 2) = 2x  →  [2, 0]
    let x3 = vec![q(1), q(0), q(0), q(0)];
    let m = vec![q(1), q(0), q(-2)];
    assert_eq!(
        reduce_mod_monic(&c, &x3, &m, Poly1Order::HighFirst).unwrap(),
        vec![q(2), q(0)]
    );

    // x^4 mod (x^2 - 2) = 4  →  [4]
    let x4 = vec![q(1), q(0), q(0), q(0), q(0)];
    assert_eq!(
        reduce_mod_monic(&c, &x4, &m, Poly1Order::HighFirst).unwrap(),
        vec![q(4)]
    );
}

#[test]
// **Pipeline private** — `high_first_reduce_skips_non_monic_modulus`
fn high_first_reduce_skips_non_monic_modulus() {
    let c = ctx();
    let p = vec![q(1), q(0), q(0), q(0)];
    // 2x^2 - 4 — leading coeff not ±1
    let m = vec![q(2), q(0), q(-4)];
    assert_eq!(
        reduce_mod_monic(&c, &p, &m, Poly1Order::HighFirst).unwrap(),
        p
    );
}

#[test]
// **Pipeline private** — `high_first_inv_mod_and_ext_gcd`
fn high_first_inv_mod_and_ext_gcd() {
    let c = ctx();
    // gcd(x, x^2-2) = 1 (coprime); ext_gcd returns constant gcd
    let x = vec![q(1), q(0)];
    let m = vec![q(1), q(0), q(-2)];
    let (g, _s) = ext_gcd(&c, &x, &m, Poly1Order::HighFirst).unwrap();
    assert_eq!(g, vec![q(1)]);

    // 2^{-1} = 1/2 in Q[x]/(x^2-2)
    let two = vec![q(2)];
    let inv = inv_mod(&c, &two, &m, Poly1Order::HighFirst).unwrap();
    let prod = mul(&c, &two, &inv, Poly1Order::HighFirst).unwrap();
    assert_eq!(
        reduce_mod_monic(&c, &prod, &m, Poly1Order::HighFirst).unwrap(),
        vec![q(1)]
    );
}

#[test]
// **Pipeline private** — `ascending_order_matches_reversed_high_first_mul`
fn ascending_order_matches_reversed_high_first_mul() {
    let c = ctx();
    // (1+x)(2+x) = 2 + 3x + x^2
    let a = vec![q(1), q(1)];
    let b = vec![q(2), q(1)];
    assert_eq!(
        mul(&c, &a, &b, Poly1Order::Ascending).unwrap(),
        vec![q(2), q(3), q(1)]
    );
}

#[test]
// **Pipeline private** — `reverse_coeffs_involution`
fn reverse_coeffs_involution() {
    let hi = vec![q(1), q(3), q(2)];
    assert_eq!(super::convert::reverse_coeffs(&super::convert::reverse_coeffs(&hi)), hi);
    let asc = vec![q(2), q(3), q(1)];
    assert_eq!(super::convert::ascending_to_dense_high_first(&asc), hi);
    assert_eq!(super::convert::dense_high_first_to_ascending(&hi), asc);
}

#[test]
// **Pipeline private** — `sparse_dense_high_first_roundtrip`
fn sparse_dense_high_first_roundtrip() {
    use std::sync::Arc;

    use crate::resultant::univariate_coeffs_ascending;

    let var: Arc<str> = Arc::from("x");
    // x^2 + 3x + 2
    let hi = vec![q(1), q(3), q(2)];
    let sparse = super::convert::dense_high_first_to_sparse(&hi, &var);
    assert_eq!(
        super::convert::sparse_ascending_to_dense_high_first(&sparse, &var),
        hi
    );
    assert_eq!(
        univariate_coeffs_ascending(&sparse, &var),
        super::convert::dense_high_first_to_ascending(&hi)
    );
}

#[test]
// **Pipeline private** — `dense_div_rem_matches_univariate_div_rem_wrt`
fn dense_div_rem_matches_univariate_div_rem_wrt() {
    use std::sync::Arc;

    use crate::resultant::univariate_coeffs_ascending;
    use crate::subresultant::univariate_div_rem_wrt;

    let var: Arc<str> = Arc::from("x");
    let a = vec![q(1), q(3), q(2)];
    let b = vec![q(1), q(1)];
    let (dq, dr) = div_rem(&ctx(), &a, &b, Poly1Order::HighFirst).unwrap();
    let pa = super::convert::dense_high_first_to_sparse(&a, &var);
    let pb = super::convert::dense_high_first_to_sparse(&b, &var);
    let (sq, sr) = univariate_div_rem_wrt(&pa, &pb, &var);
    assert_eq!(
        univariate_coeffs_ascending(&super::convert::dense_high_first_to_sparse(&dq, &var), &var),
        univariate_coeffs_ascending(&sq, &var)
    );
    assert_eq!(
        univariate_coeffs_ascending(&super::convert::dense_high_first_to_sparse(&dr, &var), &var),
        univariate_coeffs_ascending(&sr, &var)
    );
}
