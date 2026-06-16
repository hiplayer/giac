//! Univariate factorization over ℚ via Zassenhaus (GIAC `modfactor.cc` Phase B MVP).

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::error::{PolyError, PolyResult};
use crate::modint::ModInt;
use crate::modular::PolyMod;
use crate::monomial::{Monomial, Var};
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

use super::fpx::{self, factor_fpx};

const PRIMES: &[i64] = &[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];

fn x_var() -> Var {
    Var::from("x")
}

fn monomial_x_pow(e: u64) -> Monomial {
    if e == 0 {
        Monomial::one()
    } else {
        let x = x_var();
        let mut m = Monomial::var(x.clone());
        for _ in 1..e {
            m = m.mul(&Monomial::var(x.clone()));
        }
        m
    }
}

fn is_int_poly(p: &Poly) -> bool {
    p.terms.values().all(|c| c.denom().is_one())
}

fn integer_coeffs(p: &Poly, var: &Var) -> Option<Vec<BigInt>> {
    if !is_int_poly(p) {
        return None;
    }
    let d = univariate_degree(p, var);
    let mut out = vec![BigInt::zero(); d as usize + 1];
    for e in 0..=d {
        let c = coeff_at(p, var, e);
        if !c.denom().is_one() {
            return None;
        }
        out[e as usize] = c.numer().clone();
    }
    Some(out)
}

fn mignotte_bound(p: &Poly, var: &Var) -> BigInt {
    let d = univariate_degree(p, var) as i64;
    if d <= 0 {
        return BigInt::one();
    }
    let mut norm = BigInt::zero();
    for e in 0..=d as u64 {
        let c = coeff_at(p, var, e);
        if c.is_zero() {
            continue;
        }
        let a = c.numer().abs();
        if a > norm {
            norm = a;
        }
    }
    let mut n = BigInt::from(d + 1);
    if d % 2 != 0 {
        n *= 2;
    }
    let sqrt_n = n.sqrt();
    let pow2 = BigInt::from(2).pow((1 + d / 2) as u32);
    (sqrt_n + 1) * norm * pow2
}

fn coeff_mod(p: &PolyMod, e: u64) -> ModInt {
    let var = x_var();
    for (m, c) in &p.terms {
        if m.exp_of(&var) == e && m.iter().all(|(v, _)| v == &var) {
            return c.clone();
        }
    }
    ModInt::new(BigInt::zero(), p.modulus.clone()).unwrap()
}

fn poly_mod_from_coeffs(coeffs: &[BigInt], modulus: &BigInt) -> PolyResult<PolyMod> {
    let mut terms = std::collections::BTreeMap::new();
    for (e, c) in coeffs.iter().enumerate() {
        let r = c.mod_floor(modulus);
        if !r.is_zero() {
            terms.insert(
                monomial_x_pow(e as u64),
                ModInt::new(r, modulus.clone())?,
            );
        }
    }
    Ok(PolyMod {
        terms,
        modulus: modulus.clone(),
    })
}

fn poly_mod_from_poly(p: &Poly, var: &Var, modulus: &BigInt) -> PolyResult<PolyMod> {
    let coeffs = integer_coeffs(p, var).ok_or(PolyError::TypeError("non-integer poly"))?;
    poly_mod_from_coeffs(&coeffs, modulus)
}

fn make_monic_mod(p: &PolyMod) -> PolyResult<PolyMod> {
    let d = fpx::degree(p);
    if d == 0 {
        return Ok(p.clone());
    }
    let lc = coeff_mod(p, d);
    if lc.is_one() {
        return Ok(p.clone());
    }
    let inv = lc.inv()?;
    let mut coeffs = Vec::new();
    for e in 0..=d {
        coeffs.push(coeff_mod(p, e).mul(&inv)?.val);
    }
    poly_mod_from_coeffs(&coeffs, &p.modulus)
}

fn at_modulus(p: &PolyMod, modulus: &BigInt) -> PolyResult<PolyMod> {
    let d = fpx::degree(p);
    let mut coeffs = Vec::new();
    for e in 0..=d {
        coeffs.push(coeff_mod(p, e).val.clone());
    }
    poly_mod_from_coeffs(&coeffs, modulus)
}

fn derivative_mod(p: &PolyMod) -> PolyMod {
    let d = fpx::degree(p);
    let mut terms = std::collections::BTreeMap::new();
    for e in 1..=d {
        let c = coeff_mod(p, e);
        let scaled = c
            .mul(&ModInt::new(BigInt::from(e as i64), p.modulus.clone()).unwrap())
            .unwrap();
        if !scaled.is_zero() {
            terms.insert(monomial_x_pow(e - 1), scaled);
        }
    }
    PolyMod {
        terms,
        modulus: p.modulus.clone(),
    }
}

fn is_square_free_mod(p: &Poly, var: &Var, prime: i64) -> bool {
    let pm = match crate::modp(p, prime) {
        Ok(m) => m,
        Err(_) => return false,
    };
    let dp = derivative_mod(&pm);
    pm.gcd(&dp).map(|g| fpx::degree(&g) == 0).unwrap_or(false)
}

fn extgcd_mod(a: &PolyMod, b: &PolyMod) -> Option<(PolyMod, PolyMod, PolyMod)> {
    let mut old_r = a.clone();
    let mut r = b.clone();
    let mut old_s = PolyMod::one(a.modulus.clone());
    let mut s = PolyMod::zero(a.modulus.clone());
    let mut old_t = PolyMod::zero(a.modulus.clone());
    let mut t = PolyMod::one(a.modulus.clone());
    while !r.is_zero() {
        let (q, rem) = old_r.div_rem(&r).ok()?;
        old_r = r;
        r = rem;
        let new_s = old_s.sub(&q.mul(&s).ok()?).ok()?;
        let new_t = old_t.sub(&q.mul(&t).ok()?).ok()?;
        old_s = s;
        s = new_s;
        old_t = t;
        t = new_t;
    }
    Some((old_r, old_s, old_t))
}

fn polymod_to_int_poly(p: &PolyMod, var: &Var) -> Poly {
    let d = fpx::degree(p);
    let half = &p.modulus / 2;
    let mut out = Poly::zero();
    for e in 0..=d {
        let mut c = coeff_mod(p, e).val.clone();
        if c > half {
            c -= &p.modulus;
        }
        if c.is_zero() {
            continue;
        }
        let term = if e == 0 {
            Poly::constant(Ratio::from_integer(c))
        } else {
            Poly::constant(Ratio::from_integer(c)).mul(&Poly::var(var.clone()).pow(e))
        };
        out = out.add(&term);
    }
    out
}

fn divides_exact_quotient(num: &Poly, den: &Poly) -> Option<Poly> {
    let (q, r) = num.div_rem(den);
    if r.is_zero() {
        Some(q)
    } else {
        None
    }
}

/// Linear Hensel lift for two monic coprime factors mod `p` (GIAC `liftl`, n=2).
fn hensel_lift_two(
    q: &Poly,
    var: &Var,
    prime: i64,
    f0: PolyMod,
    g0: PolyMod,
    bound: &BigInt,
) -> Option<(PolyMod, PolyMod)> {
    let p = BigInt::from(prime);
    let (g, s, t) = extgcd_mod(&f0, &g0)?;
    if fpx::degree(&g) != 0 {
        return None;
    }
    let mut f = f0;
    let mut g = g0;
    let mut mod_k = p.clone();
    while &mod_k < bound {
        let mod_next = &mod_k * &p;
        let q_mod = poly_mod_from_poly(q, var, &mod_next).ok()?;
        let pi = f.mul(&g).ok()?;
        let pi = at_modulus(&pi, &mod_next).ok()?;
        let mut diff = q_mod.sub(&pi).ok()?;
        diff = coeff_div_mod(&diff, &mod_k).ok()?;
        diff = at_modulus(&diff, &p).ok()?;
        let f_orig = at_modulus(&f, &p).ok()?;
        let g_orig = at_modulus(&g, &p).ok()?;
        let s_p = at_modulus(&s, &p).ok()?;
        let t_p = at_modulus(&t, &p).ok()?;
        let df = lift_correction(&diff, &t_p, &f_orig, &mod_k).ok()?;
        let dg = lift_correction(&diff, &s_p, &g_orig, &mod_k).ok()?;
        f = f.add(&df).ok()?;
        g = g.add(&dg).ok()?;
        f.modulus = mod_next.clone();
        g.modulus = mod_next.clone();
        mod_k = mod_next;
    }
    Some((f, g))
}

fn coeff_div_mod(p: &PolyMod, d: &BigInt) -> PolyResult<PolyMod> {
    let deg = fpx::degree(p);
    let mut coeffs = Vec::new();
    for e in 0..=deg {
        let c = coeff_mod(p, e).val.clone();
        if &c % d != BigInt::zero() {
            return Err(PolyError::NotImplemented("non-exact coeff division"));
        }
        coeffs.push(c / d);
    }
    poly_mod_from_coeffs(&coeffs, &p.modulus)
}

fn lift_correction(
    q: &PolyMod,
    u: &PolyMod,
    denom: &PolyMod,
    scale: &BigInt,
) -> PolyResult<PolyMod> {
    let prod = q.mul(u)?;
    let (_, rem) = prod.div_rem(denom)?;
    let mut coeffs = Vec::new();
    for e in 0..=fpx::degree(&rem) {
        coeffs.push(coeff_mod(&rem, e).val.clone() * scale);
    }
    poly_mod_from_coeffs(&coeffs, &q.modulus)
}

pub fn try_zassenhaus_factor(g: &Poly, var: &Var) -> Option<Vec<Poly>> {
    if !is_int_poly(g) {
        return None;
    }
    let d = univariate_degree(g, var);
    if d < 3 {
        return None;
    }
    let bound = mignotte_bound(g, var) * BigInt::from(2);
    for &prime in PRIMES {
        let lc = coeff_at(g, var, d);
        if lc.is_zero() || (lc.numer() % BigInt::from(prime)) == BigInt::zero() {
            continue;
        }
        if !is_square_free_mod(g, var, prime) {
            continue;
        }
        let pm = crate::modp(g, prime).ok()?;
        let mut facs = factor_fpx(&pm).ok()?;
        if facs.len() != 2 {
            continue;
        }
        facs[0] = make_monic_mod(&facs[0]).ok()?;
        facs[1] = make_monic_mod(&facs[1]).ok()?;
        let lifted = match hensel_lift_two(g, var, prime, facs[0].clone(), facs[1].clone(), &bound) {
            Some(v) => v,
            None => continue,
        };
        let f_int = polymod_to_int_poly(&lifted.0, var);
        let g_int = polymod_to_int_poly(&lifted.1, var);
        if let Some(quo) = divides_exact_quotient(g, &f_int) {
            if univariate_degree(&quo, var) > 0 {
                return Some(vec![f_int, quo]);
            }
        }
        if let Some(quo) = divides_exact_quotient(g, &g_int) {
            if univariate_degree(&quo, var) > 0 {
                return Some(vec![g_int, quo]);
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic1() -> Poly {
        let x = Poly::var("x");
        x.pow(3).sub(&x).add(&Poly::one())
    }

    fn cubic2() -> Poly {
        let x = Poly::var("x");
        x.pow(3).add(&x).add(&Poly::one())
    }

    #[test]
    fn zassenhaus_at_41() {
        let p = cubic1().mul(&cubic2());
        let pm = crate::modp(&p, 41).unwrap();
        let n = factor_fpx(&pm).unwrap().len();
        assert_eq!(n, 2, "factor count mod 41");
        let f0 = make_monic_mod(&factor_fpx(&pm).unwrap()[0]).unwrap();
        let g0 = make_monic_mod(&factor_fpx(&pm).unwrap()[1]).unwrap();
        let bound = mignotte_bound(&p, &Var::from("x")) * BigInt::from(2);
        let lifted = hensel_lift_two(&p, &Var::from("x"), 41, f0.clone(), g0.clone(), &bound);
        assert!(lifted.is_some(), "hensel lift");
        let (fl, gl) = lifted.unwrap();
        let fi = polymod_to_int_poly(&fl, &Var::from("x"));
        assert!(divides_exact_quotient(&p, &fi).is_some(), "div fi");
        let f = try_zassenhaus_factor(&p, &Var::from("x"));
        assert!(f.is_some(), "zassenhaus at p=41");
    }

    #[test]
    fn zassenhaus_two_cubics() {
        let p = cubic1().mul(&cubic2());
        let f = try_zassenhaus_factor(&p, &Var::from("x")).expect("zassenhaus");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }
}
