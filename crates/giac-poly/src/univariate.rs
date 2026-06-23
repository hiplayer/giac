//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::error::{EvalError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

// **Pipeline private** — `univariate_coeffs`
fn univariate_coeffs(p: &Poly, var: &Var) -> Vec<Ratio<BigInt>> {
    crate::resultant::univariate_coeffs_ascending(p, var)
}

// **Pipeline private** — `trim_coeffs`
fn trim_coeffs(coeffs: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>> {
    let mut out = coeffs.to_vec();
    while out.len() > 1 && out.last().is_some_and(|c| c.is_zero()) {
        out.pop();
    }
    if out.is_empty() {
        vec![Ratio::zero()]
    } else {
        out
    }
}

// **Pipeline private** — `poly_from_coeffs`
fn poly_from_coeffs(var: &Var, coeffs: &[Ratio<BigInt>]) -> Poly {
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

// **Pipeline private** — `univariate_div_rem`
fn univariate_div_rem(
    a: &[Ratio<BigInt>],
    b: &[Ratio<BigInt>],
) -> (Vec<Ratio<BigInt>>, Vec<Ratio<BigInt>>) {
    let mut r = trim_coeffs(a);
    let b = trim_coeffs(b);
    if b.len() == 1 && b[0].is_zero() {
        return (vec![Ratio::zero()], r);
    }
    if b.len() == 1 {
        let q: Vec<_> = r.iter().map(|c| c / &b[0]).collect();
        return (trim_coeffs(&q), vec![Ratio::zero()]);
    }
    let db = b.len() - 1;
    let mut q = vec![Ratio::zero(); r.len().saturating_sub(db).max(1)];
    while r.len() > db {
        let da = r.len() - 1;
        if r[da].is_zero() {
            r.pop();
            if r.is_empty() {
                r = vec![Ratio::zero()];
            }
            continue;
        }
        let coeff = r[da].clone() / b[db].clone();
        let shift = da - db;
        if shift >= q.len() {
            q.resize(shift + 1, Ratio::zero());
        }
        q[shift] += coeff.clone();
        for j in 0..=db {
            r[j + shift] -= coeff.clone() * b[j].clone();
        }
        r = trim_coeffs(&r);
    }
    (trim_coeffs(&q), r)
}

// **Pipeline private** — `coeffs_to_integer_primitive`
fn coeffs_to_integer_primitive(coeffs: &[Ratio<BigInt>]) -> Vec<BigInt> {
    if coeffs.is_empty() {
        return vec![BigInt::zero()];
    }
    let mut den = BigInt::one();
    for c in coeffs {
        den = den.lcm(c.denom());
    }
    let mut ic: Vec<BigInt> = coeffs
        .iter()
        .map(|c| (c * Ratio::from_integer(den.clone())).numer().clone())
        .collect();
    let mut g = BigInt::zero();
    for c in &ic {
        if c.is_zero() {
            continue;
        }
        g = if g.is_zero() { c.abs() } else { g.gcd(c) };
    }
    if g.is_zero() {
        return vec![BigInt::zero()];
    }
    for c in &mut ic {
        *c /= &g;
    }
    while ic.len() > 1 && ic.last().is_some_and(|c| c.is_zero()) {
        ic.pop();
    }
    ic
}

// **Pipeline private** — `is_zero_int`
fn is_zero_int(c: &[BigInt]) -> bool {
    trim_int(c) == vec![BigInt::zero()]
}

// **Pipeline private** — `int_exact_div_rem`
fn int_exact_div_rem(a: &[BigInt], b: &[BigInt]) -> (Vec<BigInt>, Vec<BigInt>) {
    let mut r = trim_int(a);
    let b = trim_int(b);
    if b.len() == 1 {
        if b[0].is_zero() {
            return (vec![BigInt::zero()], r);
        }
        let q: Vec<BigInt> = r.iter().map(|c| c / &b[0]).collect();
        return (trim_int(&q), vec![BigInt::zero()]);
    }
    let db = b.len() - 1;
    let mut q = vec![BigInt::zero(); r.len().saturating_sub(db).max(1)];
    while r.len() > db {
        let da = r.len() - 1;
        if r[da].is_zero() {
            r.pop();
            if r.is_empty() {
                r = vec![BigInt::zero()];
            }
            continue;
        }
        let coeff = &r[da] / &b[db];
        let shift = da - db;
        if shift >= q.len() {
            q.resize(shift + 1, BigInt::zero());
        }
        q[shift] += &coeff;
        for j in 0..=db {
            r[j + shift] -= &coeff * &b[j];
        }
        r = trim_int(&r);
    }
    (trim_int(&q), r)
}

// **Pipeline private** — `pseudo_remainder`
fn pseudo_remainder(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let a = trim_int(a);
    let b = trim_int(b);
    if is_zero_int(&b) || a.len() < b.len() {
        return a;
    }
    let da = a.len() - 1;
    let db = b.len() - 1;
    let exp = (da - db + 1) as u32;
    let scale = b[db].pow(exp);
    let scaled: Vec<BigInt> = a.iter().map(|c| c * &scale).collect();
    let (_, r) = int_exact_div_rem(&scaled, &b);
    coeffs_to_integer_primitive(
        &r.iter()
            .map(|c| Ratio::from_integer(c.clone()))
            .collect::<Vec<_>>(),
    )
}


// **Pipeline private** — `trim_int`
fn trim_int(c: &[BigInt]) -> Vec<BigInt> {
    let mut out = c.to_vec();
    while out.len() > 1 && out.last().is_some_and(|v| v.is_zero()) {
        out.pop();
    }
    if out.is_empty() {
        vec![BigInt::zero()]
    } else {
        out
    }
}

// **Pipeline private** — `poly_from_int_coeffs`
fn poly_from_int_coeffs(var: &Var, coeffs: &[BigInt]) -> Poly {
    let trimmed = trim_int(coeffs);
    let ratios: Vec<Ratio<BigInt>> = trimmed
        .iter()
        .map(|c| Ratio::from_integer(c.clone()))
        .collect();
    monic_univariate(var, &ratios)
}

// **Pipeline private** — `monic_univariate`
fn monic_univariate(var: &Var, coeffs: &[Ratio<BigInt>]) -> Poly {
    let trimmed = trim_coeffs(coeffs);
    let lc = trimmed.last().cloned().unwrap_or_else(Ratio::zero);
    if lc.is_zero() || lc.is_one() {
        return poly_from_coeffs(var, &trimmed);
    }
    let inv = Ratio::one() / lc;
    let monic: Vec<_> = trimmed.iter().map(|c| c * &inv).collect();
    poly_from_coeffs(var, &monic)
}

/// ∂p/∂x for univariate `p` in `var`.
/// **Stable** — derivative w.r.t. var
pub fn univariate_derivative(p: &Poly, var: &Var) -> Poly {
    let deg = univariate_degree(p, var);
    let mut out = Poly::zero();
    for e in 1..=deg {
        let c = coeff_at(p, var, e);
        if c.is_zero() {
            continue;
        }
        let coef = c * Ratio::from_integer(BigInt::from(e));
        let term = if e == 1 {
            Poly::constant(coef)
        } else {
            Poly::constant(coef).mul(&Poly::var(var.clone()).pow(e - 1))
        };
        out = out.add(&term);
    }
    out
}

/// Square-free factorization `p = ∏ f_k^k` (giac `Tsqff_char0`).
/// **Stable** — Yun square-free factors
pub fn square_free_factorization(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    let ring = crate::square_free::UniVarRing { var };
    crate::square_free::square_free_yun(&ring, p)
}

/// Product of distinct square-free factors (`p = ∏ f_k^k` → `∏ f_k`).
/// **Stable** — product of square-free factors
pub fn square_free_part(p: &Poly, var: &Var) -> PolyResult<Poly> {
    if p.is_zero() {
        return Err(EvalError::TypeError("zero polynomial"));
    }
    let mut prod = Poly::one();
    for (g, _) in square_free_factorization(p, var)? {
        prod = prod.mul(&g);
    }
    Ok(prod)
}

/// Substitute `var -> sub` in univariate polynomial `p`.
/// **Stable** — substitute var → Poly
pub fn substitute_univariate(p: &Poly, var: &Var, sub: &Poly) -> Poly {
    let d = univariate_degree(p, var);
    let mut out = Poly::zero();
    for e in 0..=d {
        let c = coeff_at(p, var, e);
        if c.is_zero() {
            continue;
        }
        out = out.add(&Poly::constant(c).mul(&sub.pow(e)));
    }
    out
}

/// Product of square-free factors with odd multiplicity (giac `sturm` / `sturmab` convention).
/// **Stable** — odd multiplicity factor
pub fn odd_multiplicity_part(p: &Poly, var: &Var) -> PolyResult<Poly> {
    if p.is_zero() {
        return Err(EvalError::TypeError("zero polynomial"));
    }
    let xv = Poly::var(var.clone());
    let mut w = p.clone();
    let mut x_power = 0usize;
    while coeff_at(&w, var, 0).is_zero() {
        let Some(q) = univariate_div_exact(&w, &xv, var) else {
            break;
        };
        w = q;
        x_power += 1;
    }
    let mut odd = odd_part_core(&w, var)?;
    if x_power % 2 == 1 {
        odd = xv.mul(&odd);
    }
    Ok(odd)
}

// **Pipeline private** — `odd_part_core`
fn odd_part_core(p: &Poly, var: &Var) -> PolyResult<Poly> {
    let ring = crate::square_free::UniVarRing { var };
    let mut w = p.clone();
    let mut y = univariate_derivative(p, var);
    crate::square_free::gcd_reduce(&ring, &mut w, &mut y)?;
    y = y.sub(&univariate_derivative(&w, var));

    let mut odd = Poly::one();
    let mut k = 1usize;
    let max_k = univariate_degree(p, var) as usize + 2;
    while !y.is_zero() && k <= max_k {
        let g = crate::square_free::gcd_reduce(&ring, &mut w, &mut y)?;
        if !g.is_one() && k % 2 == 1 {
            odd = odd.mul(&g);
        }
        y = y.sub(&univariate_derivative(&w, var));
        k += 1;
    }
    if k % 2 == 1 && !w.is_one() {
        odd = odd.mul(&w);
    }
    Ok(odd)
}

/// Subresultant gcd for univariate polynomials in `var`.
/// **Stable** — univariate gcd
pub fn gcd_univariate(p: &Poly, q: &Poly, var: &Var) -> Poly {
    univariate_gcd(p, q, var)
}

// **Pipeline private** — `univariate_gcd`
pub(crate) fn univariate_gcd(p: &Poly, q: &Poly, var: &Var) -> Poly {
    let mut a = coeffs_to_integer_primitive(&univariate_coeffs(p, var));
    let mut b = coeffs_to_integer_primitive(&univariate_coeffs(q, var));
    if is_zero_int(&a) {
        return poly_from_int_coeffs(var, &b);
    }
    if is_zero_int(&b) {
        return poly_from_int_coeffs(var, &a);
    }
    loop {
        if is_zero_int(&b) {
            return poly_from_int_coeffs(var, &a);
        }
        let r = pseudo_remainder(&a, &b);
        if is_zero_int(&r) {
            return poly_from_int_coeffs(var, &b);
        }
        a = b;
        b = r;
    }
}

// **Pipeline private** — `univariate_div_exact`
pub(crate) fn univariate_div_exact(p: &Poly, d: &Poly, var: &Var) -> Option<Poly> {
    let (q, r) = univariate_div_rem(&univariate_coeffs(p, var), &univariate_coeffs(d, var));
    if r.len() == 1 && r[0].is_zero() {
        Some(poly_from_coeffs(var, &q))
    } else {
        None
    }
}

// **Pipeline private** — `univariate_rem`
fn univariate_rem(a: &Poly, b: &Poly, var: &Var) -> Poly {
    let (_, r) = univariate_div_rem(&univariate_coeffs(a, var), &univariate_coeffs(b, var));
    poly_from_coeffs(var, &r)
}

/// Sturm sequence for univariate `p` (classical: P0=squarefree(p), P1=p', P_{i+1} = -rem(P_{i-1}, P_i)).
/// **Stable** — Sturm chain
pub fn sturm_sequence(p: &Poly, var: &Var) -> PolyResult<Vec<Poly>> {
    if univariate_degree(p, var) == 0 {
        return Err(EvalError::TypeError("constant polynomial"));
    }
    let deriv = univariate_derivative(p, var);
    let g = univariate_gcd(p, &deriv, var);
    let sq = univariate_div_exact(p, &g, var).unwrap_or_else(|| p.clone());
    let mut seq = vec![sq.clone(), univariate_derivative(&sq, var)];
    while seq.len() >= 2 {
        let n = seq.len();
        let r = univariate_rem(&seq[n - 2], &seq[n - 1], var);
        if r.is_zero() {
            break;
        }
        seq.push(r.mul_scalar(&Ratio::from_integer(-BigInt::one())));
    }
    Ok(seq)
}

/// Evaluate univariate polynomial at a rational point.
/// **Stable** — Horner eval
pub fn eval_univariate_at(p: &Poly, var: &Var, x: &Ratio<BigInt>) -> Ratio<BigInt> {
    p.horner(var, x)
}

// **Pipeline private** — `sign_of_ratio`
fn sign_of_ratio(r: &Ratio<BigInt>) -> i8 {
    if r.is_zero() {
        0
    } else if r > &Ratio::zero() {
        1
    } else {
        -1
    }
}

/// Count sign changes in `values`, skipping zeros (Sturm's theorem convention).
/// **Stable** — sign change count in sequence
pub fn sign_variations(values: &[Ratio<BigInt>]) -> usize {
    let signs: Vec<i8> = values
        .iter()
        .map(sign_of_ratio)
        .filter(|&s| s != 0)
        .collect();
    if signs.len() < 2 {
        return 0;
    }
    signs.windows(2).filter(|w| w[0] != w[1]).count()
}

/// Sturm sign-variation count V(a) for sequence `seq` at point `a`.
/// **Stable** — Sturm sign count at point
pub fn sturm_sign_variations_at(seq: &[Poly], var: &Var, a: &Ratio<BigInt>) -> usize {
    let values: Vec<_> = seq.iter().map(|p| eval_univariate_at(p, var, a)).collect();
    sign_variations(&values)
}

/// Root count in `(a, b]` per giac `sturmab` (odd-multiplicity factors only).
/// **Stable** — root count in (a,b)
pub fn sturmab_count(
    p: &Poly,
    var: &Var,
    a: &Ratio<BigInt>,
    b: &Ratio<BigInt>,
) -> PolyResult<usize> {
    let q = odd_multiplicity_part(p, var)?;
    if univariate_degree(&q, var) == 0 {
        return Ok(0);
    }
    let seq = sturm_sequence(&q, var)?;
    let va = sturm_sign_variations_at(&seq, var, a);
    let b_eval = if b.is_zero() {
        Ratio::new(BigInt::from(-1), BigInt::from(1_000_000))
    } else {
        b.clone()
    };
    let vb = sturm_sign_variations_at(&seq, var, &b_eval);
    Ok(va.saturating_sub(vb))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn sturm_x_cubed_plus_one_has_three_polys() {
        let p = x().pow(3).add(&Poly::one());
        let seq = sturm_sequence(&p, &Var::from("x")).unwrap();
        assert_eq!(seq.len(), 3);
    }

    #[test]
    fn square_free_x_squared_times_cubic() {
        let p = x().pow(2).mul(&x().pow(3).add(&Poly::constant(Ratio::from_integer(
            BigInt::from(2),
        ))));
        let factors = square_free_factorization(&p, &Var::from("x")).unwrap();
        assert!(
            factors.iter().any(|(_, m)| *m == 1),
            "expected a multiplicity-1 factor, got {factors:?}"
        );
        let odd = odd_multiplicity_part(&p, &Var::from("x")).unwrap();
        assert_eq!(odd, x().pow(3).add(&Poly::constant(Ratio::from_integer(BigInt::from(2)))));
    }

    #[test]
    fn sturmab_x_squared_times_cubic() {
        let p = x().pow(2).mul(&x().pow(3).add(&Poly::constant(Ratio::from_integer(
            BigInt::from(2),
        ))));
        let a = Ratio::from_integer(-BigInt::from(2));
        let b = Ratio::zero();
        let count = sturmab_count(&p, &Var::from("x"), &a, &b).unwrap();
        assert_eq!(count, 1);
    }
}
