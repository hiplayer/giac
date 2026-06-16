//! Parametric resultant Res_x(P(x,t), Q(x)) eliminating `x`, yielding a polynomial in `t`.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, resultant, univariate_degree};
use crate::univariate::univariate_derivative;

/// Build `num(x) - t * den'(x)` as a polynomial in `(x, t)`.
pub fn num_minus_t_derivative(num: &Poly, den: &Poly, x: &Var, t: &Var) -> Poly {
    let dp = univariate_derivative(den, x);
    let mut p1 = embed_univariate_x(num, x);
    for (m, c) in &dp.terms {
        let xe = m.exp_of(x);
        let mut term = Poly::var(t.clone()).mul_scalar(c);
        if xe > 0 {
            term = term.mul(&Poly::var(x.clone()).pow(xe));
        }
        p1 = p1.sub(&term);
    }
    p1
}

fn embed_univariate_x(p: &Poly, x: &Var) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        if m.iter().all(|(v, _)| v == x) || m.is_const() {
            out.terms.insert(m.clone(), c.clone());
        }
    }
    out
}

/// Resultant eliminating `x`; result is univariate in `t` (via interpolation).
pub fn tresultant_eliminate_x(a: &Poly, b: &Poly, x: &Var, t: &Var) -> PolyResult<Poly> {
    let da = univariate_degree(a, x);
    let db = univariate_degree(b, x);
    if da == 0 || db == 0 {
        return Err(PolyError::TypeError("degenerate tresultant"));
    }
    let max_deg = da.max(db) as usize + 1;
    let mut samples = Vec::new();
    for k in 1..=(max_deg + 1) {
        let alpha = Ratio::from_integer(BigInt::from(k as i64));
        let ae = eval_param_poly(a, t, &alpha, x);
        let val = resultant(&ae, b, x)?;
        let c = val
            .terms
            .get(&crate::monomial::Monomial::one())
            .cloned()
            .unwrap_or_else(|| coeff_at(&val, x, 0));
        samples.push((alpha, c));
    }
    Ok(lagrange_poly(&samples, t))
}

fn lagrange_poly(samples: &[(Ratio<BigInt>, Ratio<BigInt>)], t: &Var) -> Poly {
    let mut out = Poly::zero();
    let n = samples.len();
    for i in 0..n {
        let (xi, yi) = &samples[i];
        let mut basis = Poly::constant(yi.clone());
        for j in 0..n {
            if i == j {
                continue;
            }
            let (xj, _) = &samples[j];
            let linear = Poly::var(t.clone()).sub(&Poly::constant(xj.clone()));
            let scale = Ratio::one() / (xi.clone() - xj.clone());
            basis = basis.mul(&linear.mul_scalar(&scale));
        }
        out = out.add(&basis);
    }
    out
}

/// Evaluate `p(x,t)` at `t = alpha`, returning a univariate polynomial in `x`.
pub fn eval_param_poly(p: &Poly, t: &Var, alpha: &Ratio<BigInt>, x: &Var) -> Poly {
    let mut out = Poly::zero();
    for (m, c) in &p.terms {
        let te = m.exp_of(t);
        let mut scale = c.clone();
        for _ in 0..te {
            scale *= alpha.clone();
        }
        let xe = m.exp_of(x);
        let xp = if xe == 0 {
            Poly::one()
        } else {
            Poly::var(x.clone()).pow(xe)
        };
        out = out.add(&xp.mul_scalar(&scale));
    }
    out
}

/// Rational roots of a univariate polynomial in `t` (degree ≤ 4).
pub fn rational_roots_in_t(p: &Poly, t: &Var) -> PolyResult<Vec<Ratio<BigInt>>> {
    let d = univariate_degree(p, t);
    if d == 0 {
        return Ok(vec![]);
    }
    if d == 1 {
        let a = coeff_at(p, t, 1);
        let b = coeff_at(p, t, 0);
        if a.is_zero() {
            return Ok(vec![]);
        }
        return Ok(vec![-b / a]);
    }
    let candidates = rational_root_candidates(p, t);
    let mut roots = Vec::new();
    for r in candidates {
        if eval_univariate_at_t(p, t, &r).is_zero() && !roots.contains(&r) {
            roots.push(r);
        }
    }
    Ok(roots)
}

fn eval_univariate_at_t(p: &Poly, t: &Var, val: &Ratio<BigInt>) -> Ratio<BigInt> {
    let mut sum = Ratio::zero();
    for (m, c) in &p.terms {
        if !m.iter().all(|(v, _)| v == t) && !m.is_const() {
            continue;
        }
        let e = m.exp_of(t);
        let mut pow = Ratio::one();
        for _ in 0..e {
            pow *= val.clone();
        }
        sum += c.clone() * pow;
    }
    sum
}

fn rational_root_candidates(p: &Poly, t: &Var) -> Vec<Ratio<BigInt>> {
    let d = univariate_degree(p, t);
    if d == 0 {
        return vec![];
    }
    let lc = coeff_at(p, t, d);
    let ac = coeff_at(p, t, 0);
    let mut numer_divs = vec![BigInt::from(1)];
    let mut denom_divs = vec![BigInt::from(1)];
    push_divisors(ac.numer(), &mut numer_divs);
    push_divisors(ac.denom(), &mut denom_divs);
    push_divisors(lc.numer(), &mut denom_divs);
    push_divisors(lc.denom(), &mut numer_divs);
    let mut out = Vec::new();
    for p in &numer_divs {
        for q in &denom_divs {
            if q.is_zero() {
                continue;
            }
            for sign in [1i64, -1] {
                let r = Ratio::new(p.clone() * sign, q.clone());
                if !out.contains(&r) {
                    out.push(r);
                }
            }
        }
    }
    out.sort();
    out
}

fn push_divisors(n: &BigInt, out: &mut Vec<BigInt>) {
    if n.is_zero() {
        return;
    }
    let a = n.abs();
    let mut i = BigInt::one();
    while &i * &i <= a {
        if (&a % &i).is_zero() {
            if !out.contains(&i) {
                out.push(i.clone());
            }
            let q = &a / &i;
            if !out.contains(&q) {
                out.push(q);
            }
        }
        i += BigInt::one();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn x() -> Var {
        Var::from("x")
    }

    fn t() -> Var {
        Var::from("__rt")
    }

    #[test]
    fn tresultant_one_over_x_squared_plus_one() {
        let xv = x();
        let tv = t();
        let den = Poly::var("x").pow(2).add(&Poly::one());
        let num = Poly::one();
        let p1 = num_minus_t_derivative(&num, &den, &xv, &tv);
        let r = tresultant_eliminate_x(&p1, &den, &xv, &tv).unwrap();
        let roots = rational_roots_in_t(&r, &tv).unwrap();
        assert!(!roots.is_empty());
    }

    #[test]
    fn tresultant_one_over_x_fourth_plus_one() {
        let xv = x();
        let tv = t();
        let den = Poly::var("x").pow(4).add(&Poly::one());
        let num = Poly::one();
        let p1 = num_minus_t_derivative(&num, &den, &xv, &tv);
        let r = tresultant_eliminate_x(&p1, &den, &xv, &tv).unwrap();
        assert_eq!(univariate_degree(&r, &tv), 4);
        assert!(rational_roots_in_t(&r, &tv).unwrap().is_empty());
    }
}
