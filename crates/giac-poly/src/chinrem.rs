//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use crate::error::{EvalError, PolyResult};
use crate::poly::{egcd, Poly};

/// Polynomial CRT: find `r` with `r ≡ a1 (mod m1)` and `r ≡ a2 (mod m2)`.
/// **Stable** — Chinese remainder two residues
pub fn chinrem(a1: &Poly, a2: &Poly, m1: &Poly, m2: &Poly) -> PolyResult<Poly> {
    let (g, s, _t) = egcd(m1, m2);
    if !g.is_one() {
        return Err(EvalError::TypeError("moduli not coprime"));
    }
    let diff = a2.sub(a1);
    let (_, rem) = diff.div_rem(m2);
    let q = rem.mul(&s);
    let (_, q_mod) = q.div_rem(m2);
    Ok(a1.add(&m1.mul(&q_mod)))
}

/// Combine parallel lists of residues/moduli; returns `[solution, m1*m2*...]`.
/// **Stable** — CRT fold over lists
pub fn chinrem_lists(residues: &[Poly], moduli: &[Poly]) -> PolyResult<(Poly, Poly)> {
    if residues.len() != moduli.len() || residues.is_empty() {
        return Err(EvalError::TypeError("chinrem list length"));
    }
    let mut r = residues[0].clone();
    let mut m = moduli[0].clone();
    for i in 1..residues.len() {
        r = chinrem(&r, &residues[i], &m, &moduli[i])?;
        m = m.mul(&moduli[i]);
    }
    Ok((r, m))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_rational::Ratio;

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn chinrem_two_linear_moduli() {
        let a1 = x().add(&Poly::constant(Ratio::from_integer(BigInt::from(2))));
        let a2 = x().pow(2).add(&Poly::one());
        let m1 = x().add(&Poly::one());
        let m2 = x().pow(2).add(&x()).add(&Poly::one());
        let r = chinrem(&a1, &a2, &m1, &m2).unwrap();
        let (_, rem1) = r.sub(&a1).div_rem(&m1);
        let (_, rem2) = r.sub(&a2).div_rem(&m2);
        assert!(rem1.is_zero());
        assert!(rem2.is_zero());
    }
}
