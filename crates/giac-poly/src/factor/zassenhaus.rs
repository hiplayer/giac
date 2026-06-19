//! Univariate Zassenhaus + modular Hensel lifting (integer poly → factors over ℚ).
//!
//! **Partial:** `try_zassenhaus_factor`.
//! **Pipeline private:** modular egcd, Hensel lift, factor combination recovery.

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
const MAX_COMBINE_FACTORS: usize = 12;

// **Pipeline private** — `x_var`
fn x_var() -> Var {
    Var::from("x")
}

// **Pipeline private** — `monomial_x_pow`
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

// **Pipeline private** — `is_int_poly`
fn is_int_poly(p: &Poly) -> bool {
    p.terms.values().all(|c| c.denom().is_one())
}

// **Pipeline private** — `integer_coeffs`
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

// **Pipeline private** — `mignotte_bound`
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

// **Pipeline private** — `coeff_mod`
fn coeff_mod(p: &PolyMod, e: u64) -> ModInt {
    let var = x_var();
    for (m, c) in &p.terms {
        if m.exp_of(&var) == e && m.iter().all(|(v, _)| v == &var) {
            return c.clone();
        }
    }
    ModInt::new(BigInt::zero(), p.modulus.clone()).unwrap()
}

// **Pipeline private** — `poly_mod_from_coeffs`
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

// **Pipeline private** — `poly_mod_from_poly`
fn poly_mod_from_poly(p: &Poly, var: &Var, modulus: &BigInt) -> PolyResult<PolyMod> {
    let coeffs = integer_coeffs(p, var).ok_or(PolyError::TypeError("non-integer poly"))?;
    poly_mod_from_coeffs(&coeffs, modulus)
}

// **Pipeline private** — `make_monic_mod`
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

// **Pipeline private** — `at_modulus`
fn at_modulus(p: &PolyMod, modulus: &BigInt) -> PolyResult<PolyMod> {
    let d = fpx::degree(p);
    let mut coeffs = Vec::new();
    for e in 0..=d {
        coeffs.push(coeff_mod(p, e).val.clone());
    }
    poly_mod_from_coeffs(&coeffs, modulus)
}

// **Pipeline private** — `derivative_mod`
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

// **Pipeline private** — `is_square_free_mod`
fn is_square_free_mod(p: &Poly, var: &Var, prime: i64) -> bool {
    let pm = match crate::modp(p, prime) {
        Ok(m) => m,
        Err(_) => return false,
    };
    let dp = derivative_mod(&pm);
    pm.gcd(&dp)
        .map(|g| fpx::degree(&g) == 0)
        .unwrap_or(false)
}

// **Pipeline private** — `extgcd_mod`
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

/// `Σ u[i] * Π_{j≠i} f[j] = 1` for pairwise coprime univariate `f[i]` (GIAC `egcd`).
// **Pipeline private** — `egcd_factor_list_mod`
fn egcd_factor_list_mod(factors: &[PolyMod]) -> Option<Vec<PolyMod>> {
    let n = factors.len();
    if n == 0 {
        return None;
    }
    if n == 1 {
        return Some(vec![PolyMod::one(factors[0].modulus.clone())]);
    }
    let p = factors[0].modulus.clone();
    let mut pi: Vec<PolyMod> = vec![factors[n - 1].clone()];
    for k in 1..n - 1 {
        pi.push(pi[k - 1].mul(&factors[n - k - 1]).ok()?);
    }
    let mut u: Vec<PolyMod> = Vec::with_capacity(n);
    let mut c = PolyMod::one(p.clone());
    for k in 0..n - 1 {
        let (g, v, big_u) = extgcd_mod(&factors[k], &pi[n - k - 2])?;
        if fpx::degree(&g) != 0 {
            return None;
        }
        let prod = big_u.mul(&c).ok()?;
        let (_, u_k) = prod.div_rem(&factors[k]).ok()?;
        u.push(u_k);
        let prod2 = v.mul(&c).ok()?;
        let (_, new_c) = prod2.div_rem(&pi[n - k - 2]).ok()?;
        c = new_c;
    }
    u.push(c);
    Some(u)
}

// **Pipeline private** — `polymod_to_int_poly`
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

// **Pipeline private** — `divides_exact`
fn divides_exact(num: &Poly, den: &Poly) -> bool {
    num.div_rem(den).1.is_zero()
}

// **Pipeline private** — `divides_exact_quotient`
fn divides_exact_quotient(num: &Poly, den: &Poly) -> Option<Poly> {
    let (quo, rem) = num.div_rem(den);
    if rem.is_zero() {
        Some(quo)
    } else {
        None
    }
}

// **Pipeline private** — `coeff_div_mod`
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

// **Pipeline private** — `lift_correction`
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

/// Linear Hensel lift for two monic coprime factors mod `p` (GIAC `liftl`, n=2).
// **Pipeline private** — `hensel_lift_two`
fn hensel_lift_two(
    q: &Poly,
    var: &Var,
    prime: i64,
    mut f: PolyMod,
    mut g: PolyMod,
    bound: &BigInt,
) -> Option<(PolyMod, PolyMod)> {
    let p = BigInt::from(prime);
    let (gcd, s, t) = extgcd_mod(&f, &g)?;
    if fpx::degree(&gcd) != 0 {
        return None;
    }
    let mut mod_k = p.clone();
    while &mod_k < bound {
        let mod_next = &mod_k * &p;
        let q_mod = poly_mod_from_poly(q, var, &mod_next).ok()?;
        let pi = f.mul(&g).ok()?;
        let pi = at_modulus(&pi, &mod_next).ok()?;
        let mut diff = q_mod.sub(&pi).ok()?;
        diff = coeff_div_mod(&diff, &mod_k).ok()?;
        diff = at_modulus(&diff, &p).ok()?;
        if diff.is_zero() {
            f.modulus = mod_next.clone();
            g.modulus = mod_next.clone();
            mod_k = mod_next;
            continue;
        }
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

/// Linear Hensel lift for n mod-p factors to mod p^k (GIAC `liftl`).
// **Pipeline private** — `hensel_lift_n`
fn hensel_lift_n(
    q: &Poly,
    var: &Var,
    prime: i64,
    mut factors: Vec<PolyMod>,
    bound: &BigInt,
) -> Option<Vec<PolyMod>> {
    let n = factors.len();
    if n == 0 {
        return None;
    }
    if n == 1 {
        return Some(factors);
    }
    if n == 2 {
        let (a, b) = hensel_lift_two(q, var, prime, factors[0].clone(), factors[1].clone(), bound)?;
        return Some(vec![a, b]);
    }

    let p = BigInt::from(prime);
    let orig: Vec<PolyMod> = factors
        .iter()
        .map(|f| at_modulus(f, &p))
        .collect::<Result<_, _>>()
        .ok()?;
    let bezout = egcd_factor_list_mod(&orig)?;
    let mut mod_k = p.clone();
    while &mod_k < bound {
        let mod_next = &mod_k * &p;
        let mut pi = factors[0].clone();
        for f in factors.iter().skip(1) {
            pi = pi.mul(f).ok()?;
            pi = at_modulus(&pi, &mod_next).ok()?;
        }
        let q_mod = poly_mod_from_poly(q, var, &mod_next).ok()?;
        let mut diff = q_mod.sub(&pi).ok()?;
        diff = coeff_div_mod(&diff, &mod_k).ok()?;
        diff = at_modulus(&diff, &p).ok()?;
        if diff.is_zero() {
            for f in &mut factors {
                f.modulus = mod_next.clone();
            }
            mod_k = mod_next;
            continue;
        }
        for (f, (u, o)) in factors.iter_mut().zip(bezout.iter().zip(orig.iter())) {
            let u_p = at_modulus(u, &p).ok()?;
            let df = lift_correction(&diff, &u_p, o, &mod_k).ok()?;
            *f = f.add(&df).ok()?;
            f.modulus = mod_next.clone();
        }
        mod_k = mod_next;
    }
    Some(factors)
}

// **Pipeline private** — `polymod_product`
fn polymod_product(factors: &[PolyMod], indices: &[usize]) -> Option<PolyMod> {
    if indices.is_empty() {
        return None;
    }
    let mut out = factors[indices[0]].clone();
    for &i in &indices[1..] {
        out = out.mul(&factors[i]).ok()?;
    }
    Some(out)
}

// **Pipeline private** — `recover_factors_from_lifted`
fn recover_factors_from_lifted(g: &Poly, var: &Var, lifted: &[PolyMod]) -> Option<Vec<Poly>> {
    if lifted.is_empty() {
        return None;
    }
    if lifted.len() == 1 {
        let f = polymod_to_int_poly(&lifted[0], var);
        return divides_exact_quotient(g, &f).map(|q| vec![f, q]);
    }
    if let Some(out) = extract_factors_via_combine(g, var, lifted) {
        return Some(out);
    }
    // Single lifted factor may be exact while the paired one is not (asymmetric rounding).
    for f in lifted {
        let fi = polymod_to_int_poly(f, var);
        if let Some(q) = divides_exact_quotient(g, &fi) {
            return Some(vec![fi, q]);
        }
    }
    None
}

/// Subset search on lifted mod-p^k factors (GIAC `combine` MVP).
// **Pipeline private** — `extract_factors_via_combine`
fn extract_factors_via_combine(g: &Poly, var: &Var, lifted: &[PolyMod]) -> Option<Vec<Poly>> {
    let n = lifted.len();
    if n == 0 {
        return None;
    }
    if n == 1 {
        let f = polymod_to_int_poly(&lifted[0], var);
        return divides_exact_quotient(g, &f).map(|_| vec![f]);
    }
    extract_combine_rec(g, var, lifted, &mut Vec::new())
}

// **Pipeline private** — `extract_combine_rec`
fn extract_combine_rec(
    g: &Poly,
    var: &Var,
    lifted: &[PolyMod],
    acc: &mut Vec<Poly>,
) -> Option<Vec<Poly>> {
    if lifted.is_empty() {
        if g.is_one() || univariate_degree(g, var) == 0 {
            return Some(std::mem::take(acc));
        }
        acc.push(g.clone());
        return Some(std::mem::take(acc));
    }
    if lifted.len() == 1 {
        let f = polymod_to_int_poly(&lifted[0], var);
        if let Some(quo) = divides_exact_quotient(g, &f) {
            acc.push(f);
            return extract_combine_rec(&quo, var, &[], acc);
        }
        acc.push(g.clone());
        return Some(std::mem::take(acc));
    }

    let max_k = lifted.len() / 2 + lifted.len() % 2;
    for k in 1..=max_k {
        let mut idx: Vec<usize> = (0..k).collect();
        loop {
            let prod = polymod_product(lifted, &idx)?;
            let cand = polymod_to_int_poly(&prod, var);
            if let Some(quo) = divides_exact_quotient(g, &cand) {
                let rest: Vec<PolyMod> = lifted
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| !idx.contains(i))
                    .map(|(_, f)| f.clone())
                    .collect();
                acc.push(cand);
                if let Some(out) = extract_combine_rec(&quo, var, &rest, acc) {
                    return Some(out);
                }
                acc.pop();
            }
            if k == lifted.len() {
                break;
            }
            let mut i = k;
            while i > 0 && idx[i - 1] == lifted.len() - k + i - 1 {
                i -= 1;
            }
            if i == 0 {
                break;
            }
            idx[i - 1] += 1;
            for j in i..k {
                idx[j] = idx[j - 1] + 1;
            }
        }
    }
    None
}

// **Pipeline private** — `verify_product`
fn verify_product(factors: &[Poly], g: &Poly) -> bool {
    factors.iter().fold(Poly::one(), |a, b| a.mul(b)) == *g
}

/// **Partial** — Zassenhaus+Hensel lift
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
        if facs.len() < 2 {
            continue;
        }
        if facs.len() > MAX_COMBINE_FACTORS {
            continue;
        }
        for f in &mut facs {
            *f = make_monic_mod(f).ok()?;
        }
        let lifted = match hensel_lift_n(g, var, prime, facs, &bound) {
            Some(v) => v,
            None => continue,
        };
        if let Some(out) = recover_factors_from_lifted(g, var, &lifted) {
            if verify_product(&out, g) {
                return Some(out);
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
    fn zassenhaus_two_cubics() {
        let p = cubic1().mul(&cubic2());
        let f = try_zassenhaus_factor(&p, &Var::from("x")).expect("zassenhaus");
        assert_eq!(f.len(), 2);
        assert_eq!(f.iter().fold(Poly::one(), |acc, q| acc.mul(q)), p);
    }

    #[test]
    fn zassenhaus_sophie_germain_quartic() {
        let x = Poly::var("x");
        let g = x.pow(4).add(&Poly::constant(Ratio::from_integer(BigInt::from(4))));
        let facs = try_zassenhaus_factor(&g, &Var::from("x")).expect("quartic");
        assert_eq!(facs.len(), 2);
        assert!(verify_product(&facs, &g));
        let degs: Vec<u64> = facs
            .iter()
            .map(|f| univariate_degree(f, &Var::from("x")))
            .collect();
        assert_eq!(degs, vec![2, 2]);
    }

    #[test]
    fn zassenhaus_quadratic_times_cubic() {
        // (x^2+1)(x^3-x+1): no rational roots, degree 5 → Zassenhaus path
        let x = Poly::var("x");
        let f1 = x.pow(2).add(&Poly::one());
        let f2 = x.pow(3).sub(&x).add(&Poly::one());
        let g = f1.mul(&f2);
        let facs = try_zassenhaus_factor(&g, &Var::from("x")).expect("deg-5 product");
        assert_eq!(facs.len(), 2);
        assert!(verify_product(&facs, &g));
    }
}
