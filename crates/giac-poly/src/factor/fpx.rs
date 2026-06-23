//! Irreducible factorization of polynomials over finite fields (Cantor–Zassenhaus).
//!
//! **Stable:** `factor_fpx`, `degree`.
//! **Pipeline private:** Yun square-free, distinct-degree, CZ block split, Berlekamp-style linear.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::error::{EvalError, PolyResult};
use crate::modint::ModInt;
use crate::modular::PolyMod;
use crate::monomial::{Monomial, Var};

use super::fpx_uni::{self, mod_int, var_poly};

// **Pipeline private** — `x_var`
fn x_var() -> Var {
    Var::from("x")
}

// **Pipeline private** — `mi`
fn mi(val: i64, modulus: &BigInt) -> ModInt {
    mod_int(val, modulus)
}

// **Pipeline private** — `is_poly_one`
fn is_poly_one(p: &PolyMod) -> bool {
    p.terms.len() == 1
        && p
            .terms
            .get(&Monomial::one())
            .is_some_and(|c| c.is_one())
}

/// **Stable** — total degree in `x`
pub fn degree(p: &PolyMod) -> u64 {
    fpx_uni::univariate_degree(p, &x_var())
}

// **Pipeline private** — `coeff`
fn coeff(p: &PolyMod, exp: u64) -> ModInt {
    fpx_uni::coeff_at(p, &x_var(), exp)
}

// **Pipeline private** — `from_coeffs`
fn from_coeffs(modulus: &BigInt, coeffs: &[ModInt]) -> PolyMod {
    fpx_uni::from_modint_coeffs(&x_var(), modulus, coeffs)
}

// **Pipeline private** — `x_poly`
fn x_poly(modulus: &BigInt) -> PolyMod {
    var_poly(&x_var(), modulus)
}

// **Pipeline private** — `one_poly`
fn one_poly(modulus: &BigInt) -> PolyMod {
    PolyMod::one(modulus.clone())
}

// **Pipeline private** — `mod_poly`
fn mod_poly(a: &PolyMod, m: &PolyMod) -> PolyResult<PolyMod> {
    let (_, r) = a.div_rem(m)?;
    Ok(r)
}

// **Stable** — exact division if remainder zero
fn div_exact(a: &PolyMod, b: &PolyMod) -> PolyResult<PolyMod> {
    let (q, r) = a.div_rem(b)?;
    if !r.is_zero() {
        return Err(EvalError::NotImplemented("poly division"));
    }
    Ok(q)
}

// **Pipeline private** — `make_monic`
fn make_monic(p: &PolyMod) -> PolyResult<PolyMod> {
    fpx_uni::make_monic(p, &x_var())
}

// **Pipeline private** — `derivative`
fn derivative(p: &PolyMod) -> PolyResult<PolyMod> {
    fpx_uni::derivative(p, &x_var())
}

// **Pipeline private** — `eval`
fn eval(p: &PolyMod, x: i64) -> PolyResult<ModInt> {
    let d = degree(p);
    let mut acc = mi(0, &p.modulus);
    let xv = mi(x, &p.modulus);
    for e in (0..=d).rev() {
        acc = acc.mul(&xv)?;
        acc = acc.add(&coeff(p, e))?;
    }
    Ok(acc)
}

// **Pipeline private** — `powmod`
fn powmod(base: &PolyMod, exp: &BigInt, modulus: &PolyMod) -> PolyResult<PolyMod> {
    let mut result = one_poly(&modulus.modulus);
    let mut b = mod_poly(base, modulus)?;
    let mut e = exp.clone();
    let two = BigInt::from(2);
    while !e.is_zero() {
        if (&e % &two).is_one() {
            result = mod_poly(&result.mul(&b)?, modulus)?;
        }
        b = mod_poly(&b.mul(&b)?, modulus)?;
        e >>= 1;
    }
    Ok(result)
}

/// Evaluate `f(x)` at `x = g` modulo `modulus` (polynomial composition).
// **Pipeline private** — `compose`
fn compose(f: &PolyMod, g: &PolyMod, modulus: &PolyMod) -> PolyResult<PolyMod> {
    let d = degree(f);
    if d == 0 {
        return Ok(from_coeffs(&modulus.modulus, &[coeff(f, 0)]));
    }
    let mut acc = from_coeffs(&modulus.modulus, &[coeff(f, d)]);
    for e in (0..d).rev() {
        acc = mod_poly(&acc.mul(g)?, modulus)?;
        let c = coeff(f, e);
        if !c.is_zero() {
            acc = acc.add(&from_coeffs(&modulus.modulus, &[c]))?;
        }
    }
    Ok(acc)
}

/// `f(x^p) mod modulus` (GIAC `xtoxpowerpn` without qmatrix).
// **Pipeline private** — `subst_x_to_xp`
fn subst_x_to_xp(f: &PolyMod, prime: &BigInt, modulus: &PolyMod) -> PolyResult<PolyMod> {
    let xp = powmod(&x_poly(&modulus.modulus), prime, modulus)?;
    compose(f, &xp, modulus)
}

// **Pipeline private** — `linear_factor`
fn linear_factor(root: i64, modulus: &BigInt) -> PolyMod {
    from_coeffs(modulus, &[mi(-root, modulus), mi(1, modulus)])
}

struct Lcg {
    state: u64,
}

impl Lcg {
    // **Pipeline private** — `new`
    fn new(seed: u64) -> Self {
        Self { state: seed.max(1) }
    }

    // **Stable** — `Poly::next_u64`
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }

    // **Stable** — `Poly::next_i64`
    fn next_i64(&mut self, modulus: &BigInt) -> i64 {
        let p = modulus.to_string().parse::<i64>().unwrap_or(65537);
        (self.next_u64() % p as u64) as i64
    }
}

// **Pipeline private** — `random_poly`
fn random_poly(deg: u64, modulus: &BigInt, rng: &mut Lcg) -> PolyMod {
    let mut coeffs = Vec::with_capacity(deg as usize + 1);
    for e in 0..=deg {
        if e == deg {
            coeffs.push(mi(1, modulus));
        } else {
            coeffs.push(mi(rng.next_i64(modulus), modulus));
        }
    }
    from_coeffs(modulus, &coeffs)
}

// **Pipeline private** — `FpxModRing`
struct FpxModRing;

impl crate::square_free::SquareFreeRing for FpxModRing {
    type Poly = PolyMod;

    fn is_zero(&self, p: &PolyMod) -> bool {
        p.is_zero()
    }

    fn is_one(&self, p: &PolyMod) -> bool {
        is_poly_one(p)
    }

    fn derivative(&self, p: &PolyMod) -> PolyResult<PolyMod> {
        derivative(p)
    }

    fn gcd(&self, a: &PolyMod, b: &PolyMod) -> PolyResult<PolyMod> {
        a.gcd(b)
    }

    fn div_exact(&self, a: &PolyMod, b: &PolyMod) -> PolyResult<PolyMod> {
        div_exact(a, b)
    }

    fn sub(&self, a: &PolyMod, b: &PolyMod) -> PolyResult<PolyMod> {
        a.sub(b)
    }

    fn max_exponent(&self, p: &PolyMod) -> usize {
        degree(p) as usize
    }
}

// **Pipeline private** — `square_free_yun`
fn square_free_yun(p: &PolyMod) -> PolyResult<Vec<(PolyMod, usize)>> {
    crate::square_free::square_free_yun_mod(&FpxModRing, p, degree)
}

/// GIAC `ddf`: distinct-degree factorization into blocks of fixed irreducible degree.
// **Pipeline private** — `distinct_degree_factorization`
fn distinct_degree_factorization(q: &PolyMod) -> PolyResult<Vec<(PolyMod, u64)>> {
    let prime = q.modulus.clone();
    let x = x_poly(&prime);
    let mut blocks = Vec::new();

    let x_p = powmod(&x, &prime, q)?;
    let mut ddfactor = q.gcd(&x_p.sub(&x)?)?;
    let mut qrem = div_exact(q, &ddfactor)?;
    if degree(&ddfactor) > 0 {
        blocks.push((ddfactor, 1));
    }

    let qdeg = degree(q);
    let mut i = 2u64;
    while i <= qdeg {
        let deg_rem = degree(&qrem);
        if deg_rem == 0 {
            break;
        }
        if deg_rem < 2 * i {
            blocks.push((qrem, deg_rem));
            break;
        }
        let exp = prime.pow(i as u32);
        let x_pi = powmod(&x, &exp, &qrem)?;
        ddfactor = qrem.gcd(&x_pi.sub(&x)?)?;
        let k = degree(&ddfactor);
        if k > 0 {
            blocks.push((ddfactor.clone(), i));
        }
        if k == deg_rem {
            break;
        }
        if k > 0 {
            qrem = div_exact(&qrem, &ddfactor)?;
        }
        i += 1;
    }
    Ok(blocks)
}

// **Pipeline private** — `extract_linear_factors`
fn extract_linear_factors(p: &PolyMod) -> PolyResult<Vec<PolyMod>> {
    let pval = p
        .modulus
        .to_string()
        .parse::<i64>()
        .map_err(|_| EvalError::TypeError("modulus too large"))?;
    let mut rest = p.clone();
    let mut out = Vec::new();
    for r in 0..pval {
        let lin = linear_factor(r, &p.modulus);
        loop {
            let (_, rem) = rest.div_rem(&lin)?;
            if !rem.is_zero() {
                break;
            }
            out.push(lin.clone());
            rest = div_exact(&rest, &lin)?;
        }
    }
    if !rest.is_zero() && !is_poly_one(&rest) {
        out.push(rest);
    }
    Ok(out)
}

/// GIAC `cantor_zassenhaus` for one DDF block of irreducible degree `i`.
// **Pipeline private** — `cantor_zassenhaus_block`
fn cantor_zassenhaus_block(block: &PolyMod, i: u64, rng: &mut Lcg) -> PolyResult<Vec<PolyMod>> {
    let k = degree(block);
    if k == 0 {
        return Ok(vec![]);
    }
    if k == i {
        // ponytail: DDF block already irreducible of degree i (4B §2.3.1)
        return Ok(vec![block.clone()]);
    }
    if i == 1 {
        return extract_linear_factors(block);
    }

    let prime = &block.modulus;
    let p_i64 = prime
        .to_string()
        .parse::<i64>()
        .map_err(|_| EvalError::TypeError("modulus too large"))?;

    for _attempt in 2..=50 {
        let pp = random_poly(2 * i - 1, prime, rng);
        let mut pp_acc = pp.clone();

        if p_i64 == 2 {
            let mut somme = pp.clone();
            let m = (prime.bits() as u64).saturating_mul(i).max(1);
            for _ in 1..m {
                pp_acc = mod_poly(&pp_acc.mul(&pp_acc)?, block)?;
                somme = somme.add(&pp_acc)?;
            }
            pp_acc = somme;
        } else {
            let mut ppp = pp.clone();
            for _ in 1..i {
                ppp = subst_x_to_xp(&ppp, prime, block)?;
                pp_acc = mod_poly(&pp_acc.mul(&ppp)?, block)?;
            }
            let half = (prime - 1) / 2;
            pp_acc = powmod(&pp_acc, &half, block)?;
            pp_acc = pp_acc.sub(&one_poly(prime))?;
        }

        let fact1 = block.gcd(&pp_acc)?;
        let deg = degree(&fact1);
        if deg == 0 || deg == k {
            continue;
        }
        let fact2 = div_exact(block, &fact1)?;
        let mut out = cantor_zassenhaus_block(&fact1, i, rng)?;
        out.extend(cantor_zassenhaus_block(&fact2, i, rng)?);
        return Ok(out);
    }
    Err(EvalError::NotImplemented("cantor-zassenhaus split"))
}

// **Pipeline private** — `factor_square_free`
fn factor_square_free(p: &PolyMod) -> PolyResult<Vec<PolyMod>> {
    let q = make_monic(p)?;
    let blocks = distinct_degree_factorization(&q)?;
    let mut rng = Lcg::new(0xC0FFEE);
    let mut out = Vec::new();
    for (block, deg) in blocks {
        out.extend(cantor_zassenhaus_block(&block, deg, &mut rng)?);
    }
    Ok(out)
}

/// **Stable** — Full factorization in F_p[x] into monic irreducible factors (with repetition).
pub fn factor_fpx(p: &PolyMod) -> PolyResult<Vec<PolyMod>> {
    if p.is_zero() {
        return Err(EvalError::TypeError("zero polynomial"));
    }
    if degree(p) == 0 {
        return Ok(vec![]);
    }
    let sqff = square_free_yun(p)?;
    let mut out = Vec::new();
    for (g, exp) in sqff {
        let mut facs = factor_square_free(&g)?;
        for f in &mut facs {
            *f = make_monic(f)?;
        }
        for _ in 0..exp {
            out.extend(facs.iter().cloned());
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modular::modp;
    use crate::poly::Poly;

    fn poly(s: &str) -> Poly {
        // minimal parser for tests: "x^4+1", "x^6-1"
        let x = Poly::var("x");
        match s {
            "x^4+1" => x.pow(4).add(&Poly::one()),
            "x^6-1" => x.pow(6).sub(&Poly::one()),
            "x^2+1" => x.pow(2).add(&Poly::one()),
            "x^3+x+1" => x.pow(3).add(&Poly::var("x")).add(&Poly::one()),
            _ => panic!("unknown test poly {s}"),
        }
    }

    fn assert_product(facs: &[PolyMod], p: i64, expected: &Poly) {
        let prod = facs.iter().fold(Poly::one(), |acc, f| {
            acc.mul(&super::super::modular::modpoly_to_poly(f, p))
        });
        assert_eq!(modp(&prod, p).unwrap(), modp(expected, p).unwrap());
    }

    #[test]
    fn fpx_x4_plus_1_mod_5() {
        let pm = modp(&poly("x^4+1"), 5).unwrap();
        let f = factor_fpx(&pm).unwrap();
        assert_eq!(f.len(), 2);
        assert_eq!(degree(&f[0]), 2);
        assert_eq!(degree(&f[1]), 2);
        assert_product(&f, 5, &poly("x^4+1"));
    }

    #[test]
    fn fpx_x6_minus_1_mod_7() {
        let pm = modp(&poly("x^6-1"), 7).unwrap();
        let f = factor_fpx(&pm).unwrap();
        assert_eq!(f.len(), 6);
        assert!(f.iter().all(|p| degree(p) == 1));
        assert_product(&f, 7, &poly("x^6-1"));
    }

    #[test]
    fn fpx_x2_plus_1_mod_5() {
        let pm = modp(&poly("x^2+1"), 5).unwrap();
        let f = factor_fpx(&pm).unwrap();
        assert_eq!(f.len(), 2);
        assert_product(&f, 5, &poly("x^2+1"));
    }

    #[test]
    fn fpx_cubic_irreducible_mod_7() {
        let pm = modp(&poly("x^3+x+1"), 7).unwrap();
        let f = factor_fpx(&pm).unwrap();
        assert_eq!(f.len(), 1);
        assert_eq!(degree(&f[0]), 3);
    }
}
