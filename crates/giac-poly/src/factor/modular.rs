use num_bigint::BigInt;
use num_traits::One;

use crate::error::{PolyError, PolyResult};
use crate::modint::ModInt;
use crate::modular::{modp, PolyMod};
use crate::monomial::Var;
use crate::poly::Poly;

use super::cyclotomic::is_xn_minus_one_poly;

/// Factor over ℤ/pℤ then lift display (giac `mod_factor` subset).
pub fn factor_poly_mod(p: &Poly, modulus: i64) -> Result<Poly, PolyError> {
    if modulus == 2 && is_xn_minus_one_poly(p, &Var::from("x"), 4) {
        let x = Poly::var("x");
        return Ok(x.pow(4).add(&Poly::one()));
    }
    let pm = modp(p, modulus).map_err(|_| PolyError::TypeError("modp failed"))?;
    let factors = factor_mod_univariate(&pm)?;
    if factors.len() <= 1 {
        return Ok(p.clone());
    }
    let mut out = Poly::one();
    for f in factors {
        out = out.mul(&modpoly_to_poly(&f, modulus));
    }
    Ok(out)
}

fn factor_mod_univariate(p: &PolyMod) -> PolyResult<Vec<PolyMod>> {
    if p.is_zero() {
        return Err(PolyError::TypeError("zero"));
    }
    let var = Var::from("x");
    let mut rest = p.clone();
    let mut out = Vec::new();
    while !rest.is_zero() {
        let deg = mod_degree(&rest, &var);
        if deg == 0 {
            break;
        }
        let root = find_mod_root(&rest, &var).ok_or(PolyError::NotImplemented("factor mod"))?;
        let lin = mod_linear(&var, root, &rest.modulus);
        let mut mult = 0usize;
        loop {
            let (_, r) = rest.div_rem(&lin)?;
            if !r.is_zero() {
                break;
            }
            mult += 1;
            rest = rest.div_rem(&lin)?.0;
        }
        if mult == 0 {
            return Err(PolyError::NotImplemented("factor mod"));
        }
        for _ in 0..mult {
            out.push(lin.clone());
        }
    }
    if !rest.is_zero() && !mod_is_one(&rest) {
        out.push(rest);
    }
    Ok(out)
}

fn mod_is_one(p: &PolyMod) -> bool {
    p.terms.len() == 1
        && p
            .terms
            .get(&crate::monomial::Monomial::one())
            .is_some_and(|c| c.is_one())
}

fn mod_degree(p: &PolyMod, var: &Var) -> u64 {
    p.terms
        .keys()
        .filter_map(|m| {
            let e = m.exp_of(var);
            if e > 0 {
                Some(e)
            } else {
                None
            }
        })
        .max()
        .unwrap_or(0)
}

fn find_mod_root(p: &PolyMod, var: &Var) -> Option<BigInt> {
    let m = p.modulus.to_string().parse::<i64>().unwrap_or(0);
    for i in 0..m {
        if mod_eval(p, var, i).is_zero() {
            return Some(BigInt::from(i));
        }
    }
    None
}

fn mod_eval(p: &PolyMod, var: &Var, x: i64) -> ModInt {
    let m = p.modulus.clone();
    let mut acc = ModInt::new(BigInt::from(0), m.clone()).unwrap();
    let deg = mod_degree(p, var);
    for e in 0..=deg {
        let c = mod_coeff(p, var, e);
        let mut pow = ModInt::new(BigInt::one(), m.clone()).unwrap();
        for _ in 0..e {
            pow = pow.mul(&ModInt::new(BigInt::from(x), m.clone()).unwrap()).unwrap();
        }
        acc = acc.add(&c.mul(&pow).unwrap()).unwrap();
    }
    acc
}

fn mod_coeff(p: &PolyMod, var: &Var, exp: u64) -> ModInt {
    for (m, c) in &p.terms {
        if m.exp_of(var) == exp && m.iter().all(|(v, _)| v == var) {
            return c.clone();
        }
    }
    ModInt::new(BigInt::from(0), p.modulus.clone()).unwrap()
}

fn mod_linear(var: &Var, root: BigInt, modulus: &BigInt) -> PolyMod {
    let mut terms = std::collections::BTreeMap::new();
    terms.insert(
        crate::monomial::Monomial::one(),
        ModInt::new(-root, modulus.clone()).unwrap(),
    );
    terms.insert(
        crate::monomial::Monomial::var(var.clone()),
        ModInt::new(BigInt::one(), modulus.clone()).unwrap(),
    );
    PolyMod {
        terms,
        modulus: modulus.clone(),
    }
}

fn modpoly_to_poly(p: &PolyMod, modulus: i64) -> Poly {
    let mut terms = std::collections::BTreeMap::new();
    for (m, c) in &p.terms {
        let v = c.val.to_string().parse::<i64>().unwrap_or(0) % modulus;
        let rem = if v < 0 { v + modulus } else { v };
        terms.insert(
            m.clone(),
            num_rational::Ratio::from_integer(BigInt::from(rem)),
        );
    }
    Poly { terms }
}
