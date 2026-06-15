use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::One;

use crate::{chinrem_lists, factor_into, gauss, irem, modp, resultant, smod, ModInt, Monomial, Poly, Var};

fn x() -> Poly {
    Poly::var("x")
}

#[test]
fn modp_reduces_coefficients() {
    let p = x().pow(2).add(&Poly::constant(Ratio::from_integer(BigInt::from(5))));
    let pm = modp(&p, 13).unwrap();
    assert_eq!(pm.terms.len(), 2);
}

#[test]
fn polymod_gcd_linear_mod_13() {
    let a = x()
        .pow(2)
        .mul_scalar(&Ratio::from_integer(BigInt::from(2)))
        .add(&Poly::constant(Ratio::from_integer(BigInt::from(5))));
    let b = x()
        .pow(2)
        .mul_scalar(&Ratio::from_integer(BigInt::from(5)))
        .add(&x().mul_scalar(&Ratio::from_integer(BigInt::from(2))))
        .add(&Poly::constant(Ratio::from_integer(BigInt::from(-3))));
    let pa = modp(&a, 13).unwrap();
    let pb = modp(&b, 13).unwrap();
    let g = pa.gcd(&pb).unwrap();
    assert!(g.terms.len() >= 2, "expected non-constant gcd");
}

#[test]
fn polymod_div_rem() {
    let p = x().pow(2).sub(&Poly::one());
    let d = x().sub(&Poly::one());
    let pm = modp(&p, 13).unwrap();
    let pd = modp(&d, 13).unwrap();
    let (_, r) = pm.div_rem(&pd).unwrap();
    assert!(r.is_zero());
}

#[test]
fn modint_inv() {
    let two = ModInt::from_i64(2, 13).unwrap();
    let inv = two.inv().unwrap();
    assert_eq!(two.mul(&inv).unwrap().val, BigInt::one());
}

#[test]
fn chinrem_lists_two_moduli() {
    let a1 = x().add(&Poly::constant(Ratio::from_integer(BigInt::from(2))));
    let a2 = x().pow(2).add(&Poly::one());
    let m1 = x().add(&Poly::one());
    let m2 = x().pow(2).add(&x()).add(&Poly::one());
    let (r, combined) = chinrem_lists(&[a1.clone(), a2.clone()], &[m1.clone(), m2.clone()]).unwrap();
    let (_, rem1) = r.sub(&a1).div_rem(&m1);
    let (_, rem2) = r.sub(&a2).div_rem(&m2);
    assert!(rem1.is_zero());
    assert!(rem2.is_zero());
    assert_eq!(combined, m1.mul(&m2));
}

#[test]
fn factor_into_x4_minus_1() {
    let p = x().pow(4).sub(&Poly::one());
    let f = factor_into(&p).unwrap();
    assert_eq!(f.len(), 3);
    let product = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
    assert_eq!(product, p);
}

#[test]
fn gauss_2xy() {
    let p = x()
        .mul(&Poly::var("y"))
        .mul_scalar(&Ratio::from_integer(BigInt::from(2)));
    let g = gauss(&p, &[Var::from("x"), Var::from("y")]);
    let xy = x().mul(&Poly::var("y")).mul_scalar(&Ratio::from_integer(BigInt::from(2)));
    assert_eq!(g, xy);
}

#[test]
fn resultant_zero_shared_factor() {
    let a = x().pow(2).sub(&Poly::one());
    let b = x().pow(3).sub(&Poly::one());
    let r = resultant(&a, &b, &Var::from("x")).unwrap();
    assert!(r.is_zero());
}

#[test]
fn smod_irem_values() {
    assert_eq!(smod(17, 5), 2);
    assert_eq!(irem(17, 5), 2);
}

#[test]
fn monomial_lex_order() {
    let order = [Var::from("x"), Var::from("y")];
    let xy = Monomial::var("x").mul(&Monomial::var("y"));
    let x2 = Monomial::var("x").mul(&Monomial::var("x"));
    assert!(x2.cmp_lex(&xy, &order) == std::cmp::Ordering::Greater);
    assert!(Monomial::var("x").divides(&xy));
    assert!(!xy.divides(&x2));
}

#[test]
fn leading_term_lex_picks_x_over_y() {
    let order = [Var::from("x"), Var::from("y")];
    let p = x()
        .add(&Poly::var("y").mul_scalar(&Ratio::from_integer(BigInt::from(3))));
    let (lt, _) = p.leading_term_lex(&order).unwrap();
    assert_eq!(lt, &Monomial::var("x"));
}
