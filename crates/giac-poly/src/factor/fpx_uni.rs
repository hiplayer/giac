//! Univariate view of `PolyMod` in one variable (F_p[x] / ℤ/p^kℤ[x]).
//!
//! Shared by `fpx` (mod-p factorization) and `zassenhaus` (Hensel lift).
//! **Pipeline private (`pub(crate)`):** callers must ensure `p` is univariate in `var`.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::Zero;

use crate::error::PolyResult;
use crate::modint::ModInt;
use crate::modular::PolyMod;
use crate::monomial::{Monomial, Var};
use crate::poly::Poly;

// **Pipeline private** — `monomial_pow`
pub(crate) fn monomial_pow(var: &Var, exp: u64) -> Monomial {
    if exp == 0 {
        Monomial::one()
    } else {
        let mut m = Monomial::var(var.clone());
        for _ in 1..exp {
            m = m.mul(&Monomial::var(var.clone()));
        }
        m
    }
}

// **Pipeline private** — `mod_int`
pub(crate) fn mod_int(val: i64, modulus: &BigInt) -> ModInt {
    ModInt::new(BigInt::from(val), modulus.clone()).unwrap()
}

/// Degree in `var`; missing terms count as exponent 0.
// **Pipeline private** — `univariate_degree`
pub(crate) fn univariate_degree(p: &PolyMod, var: &Var) -> u64 {
    p.terms
        .keys()
        .map(|m| m.exp_of(var))
        .max()
        .unwrap_or(0)
}

/// Coefficient of `var^e`; 0 if absent.
// **Pipeline private** — `coeff_at`
pub(crate) fn coeff_at(p: &PolyMod, var: &Var, e: u64) -> ModInt {
    for (m, c) in &p.terms {
        if m.exp_of(var) == e && m.iter().all(|(v, _)| v == var) {
            return c.clone();
        }
    }
    mod_int(0, &p.modulus)
}

// **Pipeline private** — `set_coeff`
pub(crate) fn set_coeff(
    terms: &mut BTreeMap<Monomial, ModInt>,
    var: &Var,
    exp: u64,
    c: ModInt,
) {
    let m = monomial_pow(var, exp);
    if c.is_zero() {
        terms.remove(&m);
    } else {
        terms.insert(m, c);
    }
}

/// Build univariate `PolyMod` from dense `ModInt` coefficients (index = exponent).
// **Pipeline private** — `from_modint_coeffs`
pub(crate) fn from_modint_coeffs(var: &Var, modulus: &BigInt, coeffs: &[ModInt]) -> PolyMod {
    let mut terms = BTreeMap::new();
    for (e, c) in coeffs.iter().enumerate() {
        if !c.is_zero() {
            set_coeff(&mut terms, var, e as u64, c.clone());
        }
    }
    PolyMod {
        terms,
        modulus: modulus.clone(),
    }
}

/// Build univariate `PolyMod` from integer coefficients, reducing mod `modulus`.
// **Pipeline private** — `from_bigint_coeffs`
pub(crate) fn from_bigint_coeffs(
    var: &Var,
    modulus: &BigInt,
    coeffs: &[BigInt],
) -> PolyResult<PolyMod> {
    let mut terms = BTreeMap::new();
    for (e, c) in coeffs.iter().enumerate() {
        let r = c.mod_floor(modulus);
        if !r.is_zero() {
            set_coeff(
                &mut terms,
                var,
                e as u64,
                ModInt::new(r, modulus.clone())?,
            );
        }
    }
    Ok(PolyMod {
        terms,
        modulus: modulus.clone(),
    })
}

/// Monic generator `var` over F_p.
// **Pipeline private** — `var_poly`
pub(crate) fn var_poly(var: &Var, modulus: &BigInt) -> PolyMod {
    let mut terms = BTreeMap::new();
    terms.insert(Monomial::var(var.clone()), mod_int(1, modulus));
    PolyMod {
        terms,
        modulus: modulus.clone(),
    }
}

// **Pipeline private** — `make_monic`
pub(crate) fn make_monic(p: &PolyMod, var: &Var) -> PolyResult<PolyMod> {
    let d = univariate_degree(p, var);
    if d == 0 {
        return Ok(p.clone());
    }
    let lc = coeff_at(p, var, d);
    if lc.is_one() {
        return Ok(p.clone());
    }
    let inv = lc.inv()?;
    let mut coeffs = Vec::with_capacity(d as usize + 1);
    for e in 0..=d {
        coeffs.push(coeff_at(p, var, e).mul(&inv)?);
    }
    Ok(from_modint_coeffs(var, &p.modulus, &coeffs))
}

// **Pipeline private** — `derivative`
pub(crate) fn derivative(p: &PolyMod, var: &Var) -> PolyResult<PolyMod> {
    let d = univariate_degree(p, var);
    if d == 0 {
        return Ok(PolyMod::zero(p.modulus.clone()));
    }
    let mut coeffs = Vec::new();
    for e in 1..=d {
        let c = coeff_at(p, var, e);
        let scaled = c.mul(&mod_int(e as i64, &p.modulus))?;
        if !scaled.is_zero() {
            coeffs.push((e - 1, scaled));
        }
    }
    let mut terms = BTreeMap::new();
    for (e, c) in coeffs {
        set_coeff(&mut terms, var, e, c);
    }
    Ok(PolyMod {
        terms,
        modulus: p.modulus.clone(),
    })
}

/// Re-embed coefficients under a new modulus (Hensel lift).
// **Pipeline private** — `at_modulus`
pub(crate) fn at_modulus(p: &PolyMod, var: &Var, modulus: &BigInt) -> PolyResult<PolyMod> {
    let d = univariate_degree(p, var);
    let coeffs: Vec<BigInt> = (0..=d)
        .map(|e| coeff_at(p, var, e).val.clone())
        .collect();
    from_bigint_coeffs(var, modulus, &coeffs)
}

/// Lift `PolyMod` to integer `Poly` with symmetric residue in `[-modulus/2, modulus/2]`.
// **Pipeline private** — `to_centered_int_poly`
pub(crate) fn to_centered_int_poly(p: &PolyMod, var: &Var) -> Poly {
    let d = univariate_degree(p, var);
    let half = &p.modulus / 2;
    let mut out = Poly::zero();
    for e in 0..=d {
        let mut c = coeff_at(p, var, e).val.clone();
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
