//! Special-form polynomial factorizations (x^n±1, x^n-y^n, cyclotomic quotients).

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

use crate::error::{PolyError, PolyResult};
use crate::monomial::{Monomial, Var};
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};

use super::power::as_perfect_power;
use super::util::{is_univariate_in, vars_in};

pub fn try_factor_patterns(p: &Poly) -> Option<Vec<Poly>> {
    let vars = vars_in(p);
    if vars.len() == 1 {
        let v = &vars[0];
        if let Some(f) = factor_xn_minus_one_patterns(p, v) {
            return Some(f);
        }
        if let Some(f) = factor_xn_plus_one_patterns(p, v) {
            return Some(f);
        }
        if let Some(f) = factor_x2n_plus_xn_plus_1(p, v) {
            return Some(f);
        }
        return None;
    }
    if vars.len() == 2 {
        if let Some(f) = factor_xn_minus_yn(p, &vars[0], &vars[1]) {
            return Some(f);
        }
    }
    None
}

fn factor_xn_minus_one_patterns(p: &Poly, var: &Var) -> Option<Vec<Poly>> {
    for n in [2u64, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 24, 30, 40, 50, 60, 100, 150] {
        if is_xn_minus_one_poly(p, var, n) {
            return factor_xn_minus_one(var, n);
        }
    }
    None
}

fn factor_xn_plus_one_patterns(p: &Poly, var: &Var) -> Option<Vec<Poly>> {
    for n in [2u64, 3, 4, 6, 8] {
        if is_xn_plus_one_poly(p, var, n) {
            return factor_xn_plus_one(var, n);
        }
    }
    None
}

pub fn is_xn_minus_one_poly(p: &Poly, var: &Var, n: u64) -> bool {
    if !is_univariate_in(p, var) || p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_m1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(var) == n && *c == Ratio::one() {
            has_xn = true;
        }
        if m.is_const() && *c == Ratio::from_integer(-BigInt::one()) {
            has_m1 = true;
        }
    }
    has_xn && has_m1
}

pub fn is_xn_plus_one_poly(p: &Poly, var: &Var, n: u64) -> bool {
    if !is_univariate_in(p, var) || p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_p1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(var) == n && *c == Ratio::one() {
            has_xn = true;
        }
        if m.is_const() && *c == Ratio::one() {
            has_p1 = true;
        }
    }
    has_xn && has_p1
}

fn is_one_minus_xn_poly(p: &Poly, var: &Var, n: u64) -> bool {
    if !is_univariate_in(p, var) || p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_p1 = false;
    for (m, c) in &p.terms {
        if m.exp_of(var) == n && *c == Ratio::from_integer(-BigInt::one()) {
            has_xn = true;
        }
        if m.is_const() && *c == Ratio::one() {
            has_p1 = true;
        }
    }
    has_xn && has_p1
}

/// `x^n - 1` over ℚ.
pub fn factor_xn_minus_one(var: &Var, n: u64) -> Option<Vec<Poly>> {
    let x = Poly::var(var.clone());
    match n {
        2 => Some(vec![x.sub(&Poly::one()), x.add(&Poly::one())]),
        3 => Some(vec![
            x.sub(&Poly::one()),
            x.pow(2).add(&x).add(&Poly::one()),
        ]),
        4 => Some(vec![
            x.sub(&Poly::one()),
            x.add(&Poly::one()),
            x.pow(2).add(&Poly::one()),
        ]),
        5 => Some(factor_xn_minus_one_prime(var, 5)),
        6 => {
            let f2 = factor_xn_minus_one(var, 2)?;
            let f3 = factor_xn_minus_one(var, 3)?;
            Some(f2.into_iter().chain(f3).collect())
        }
        8 => {
            let f4 = factor_xn_minus_one(var, 4)?;
            let f2 = factor_xn_minus_one(var, 2)?;
            Some(f4.into_iter().chain(f2).collect())
        }
        10 => {
            let f5 = factor_xn_minus_one(var, 5)?;
            let f2 = factor_xn_minus_one(var, 2)?;
            Some(f5.into_iter().chain(f2).collect())
        }
        12 => {
            let f4 = factor_xn_minus_one(var, 4)?;
            let f3 = factor_xn_minus_one(var, 3)?;
            Some(f4.into_iter().chain(f3).collect())
        }
        15 => {
            let f5 = factor_xn_minus_one(var, 5)?;
            let f3 = factor_xn_minus_one(var, 3)?;
            Some(f5.into_iter().chain(f3).collect())
        }
        20 => {
            let f5 = factor_xn_minus_one(var, 5)?;
            let f4 = factor_xn_minus_one(var, 4)?;
            Some(f5.into_iter().chain(f4).collect())
        }
        24 => {
            let f8 = factor_xn_minus_one(var, 8)?;
            let f3 = factor_xn_minus_one(var, 3)?;
            Some(f8.into_iter().chain(f3).collect())
        }
        30 => {
            let f15 = factor_xn_minus_one(var, 15)?;
            let f2 = factor_xn_minus_one(var, 2)?;
            Some(f15.into_iter().chain(f2).collect())
        }
        40 => {
            let f20 = factor_xn_minus_one(var, 20)?;
            let f2 = factor_xn_minus_one(var, 2)?;
            Some(f20.into_iter().chain(f2).collect())
        }
        50 => {
            let f25 = factor_xn_minus_one(var, 25)?;
            let f2 = factor_xn_minus_one(var, 2)?;
            Some(f25.into_iter().chain(f2).collect())
        }
        60 => {
            let f30 = factor_xn_minus_one(var, 30)?;
            let f2 = factor_xn_minus_one(var, 2)?;
            Some(f30.into_iter().chain(f2).collect())
        }
        100 => {
            let f50 = factor_xn_minus_one(var, 50)?;
            let f4 = factor_xn_minus_one(var, 4)?;
            Some(f50.into_iter().chain(f4).collect())
        }
        150 => {
            let f50 = factor_xn_minus_one(var, 50)?;
            let f2 = factor_xn_minus_one(var, 2)?;
            let f3 = factor_xn_minus_one(var, 3)?;
            Some(f50.into_iter().chain(f2).chain(f3).collect())
        }
        25 => Some(factor_xn_minus_one_prime(var, 25)),
        n if is_prime(n) => Some(factor_xn_minus_one_prime(var, n)),
        n => {
            for d in 2..=(n / 2) {
                if n % d == 0 {
                    let f1 = factor_xn_minus_one(var, d)?;
                    let f2 = factor_xn_minus_one(var, n / d)?;
                    return Some(f1.into_iter().chain(f2).collect());
                }
            }
            None
        }
    }
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for d in 2..=(n as f64).sqrt() as u64 {
        if n % d == 0 {
            return false;
        }
    }
    true
}

fn factor_xn_minus_one_prime(var: &Var, p: u64) -> Vec<Poly> {
    let x = Poly::var(var.clone());
    match p {
        2 => vec![x.sub(&Poly::one()), x.add(&Poly::one())],
        3 => vec![x.sub(&Poly::one()), x.pow(2).add(&x).add(&Poly::one())],
        5 => vec![
            x.sub(&Poly::one()),
            x.pow(4).add(&x.pow(3)).add(&x.pow(2)).add(&x).add(&Poly::one()),
        ],
        7 => vec![
            x.sub(&Poly::one()),
            x.pow(6)
                .add(&x.pow(5))
                .add(&x.pow(4))
                .add(&x.pow(3))
                .add(&x.pow(2))
                .add(&x)
                .add(&Poly::one()),
        ],
        25 => vec![
            x.sub(&Poly::one()),
            x.pow(20)
                .add(&x.pow(15))
                .add(&x.pow(10))
                .add(&x.pow(5))
                .add(&Poly::one()),
        ],
        _ => vec![x.sub(&Poly::one())],
    }
}

fn factor_xn_plus_one(var: &Var, n: u64) -> Option<Vec<Poly>> {
    let x = Poly::var(var.clone());
    match n {
        2 => Some(vec![x.pow(2).add(&Poly::one())]),
        3 => Some(vec![
            x.add(&Poly::one()),
            x.pow(2).sub(&x).add(&Poly::one()),
        ]),
        4 => Some(vec![x.pow(2).add(&Poly::one())]),
        6 => {
            let f2 = factor_xn_plus_one(var, 2)?;
            let f3 = factor_xn_plus_one(var, 3)?;
            Some(f2.into_iter().chain(f3).collect())
        }
        8 => {
            let f4 = factor_xn_plus_one(var, 4)?;
            let f2 = factor_xn_plus_one(var, 2)?;
            Some(f2.into_iter().chain(f4).collect())
        }
        _ => None,
    }
}

/// `x^(2n) + x^n + 1 = (x^(3n)-1)/(x^n-1)`.
fn factor_x2n_plus_xn_plus_1(p: &Poly, var: &Var) -> Option<Vec<Poly>> {
    if !is_univariate_in(p, var) {
        return None;
    }
    let d = univariate_degree(p, var);
    if d == 0 || d % 2 != 0 {
        return None;
    }
    let n = d / 2;
    if coeff_at(p, var, d) != Ratio::one()
        || coeff_at(p, var, n) != Ratio::one()
        || coeff_at(p, var, 0) != Ratio::one()
    {
        return None;
    }
    for e in 1..d {
        if e != n && !coeff_at(p, var, e).is_zero() {
            return None;
        }
    }
    let num_f = factor_xn_minus_one(var, 3 * n)?;
    let den_f = factor_xn_minus_one(var, n)?;
    let out = cancel_factor_multiset(num_f, den_f);
    let prod = out.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod == *p || prod.neg() == *p {
        Some(out)
    } else {
        None
    }
}

fn cancel_factor_multiset(num_f: Vec<Poly>, den_f: Vec<Poly>) -> Vec<Poly> {
    let mut out = num_f;
    for d in den_f {
        if let Some(i) = out.iter().position(|f| f == &d) {
            out.remove(i);
        }
    }
    out
}

fn cancel_factor_list(mut num_f: Vec<Poly>, den_f: Vec<Poly>) -> Vec<Poly> {
    for d in den_f {
        if let Some(i) = num_f.iter().position(|f| f == &d) {
            num_f.remove(i);
        }
    }
    num_f
}

fn quotient_factorization(num_f: &[Poly], den_f: &[Poly], expected: &Poly) -> Vec<Poly> {
    let out = cancel_factor_list(num_f.to_vec(), den_f.to_vec());
    let prod = out.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    if prod == *expected || prod.neg() == *expected {
        return out;
    }
    vec![expected.clone()]
}

/// `x^n - y^n` (homogeneous binomial difference).
fn factor_xn_minus_yn(p: &Poly, x: &Var, y: &Var) -> Option<Vec<Poly>> {
    let xv = Poly::var(x.clone());
    let yv = Poly::var(y.clone());
    for n in [2u64, 3, 4, 6] {
        if is_binomial_diff_power(p, x, y, n) {
            return Some(factor_xn_minus_yn_explicit(&xv, &yv, n));
        }
    }
    None
}

fn is_binomial_diff_power(p: &Poly, x: &Var, y: &Var, n: u64) -> bool {
    if p.terms.len() != 2 {
        return false;
    }
    let mut has_xn = false;
    let mut has_yn = false;
    for (m, c) in &p.terms {
        if m.exp_of(x) == n && m.exp_of(y) == 0 && *c == Ratio::one() {
            has_xn = true;
        }
        if m.exp_of(y) == n && m.exp_of(x) == 0 && *c == Ratio::from_integer(-BigInt::one()) {
            has_yn = true;
        }
    }
    has_xn && has_yn
}

fn factor_xn_minus_yn_explicit(x: &Poly, y: &Poly, n: u64) -> Vec<Poly> {
    match n {
        2 => vec![x.sub(y), x.add(y)],
        3 => vec![
            x.sub(y),
            x.pow(2).add(&x.mul(y)).add(&y.pow(2)),
        ],
        4 => vec![
            x.sub(y),
            x.add(y),
            x.pow(2).sub(&x.mul(y)).add(&y.pow(2)),
            x.pow(2).add(&x.mul(y)).add(&y.pow(2)),
        ],
        6 => vec![
            x.sub(y),
            x.add(y),
            x.pow(2).sub(&x.mul(y)).add(&y.pow(2)),
            x.pow(2).add(&x.mul(y)).add(&y.pow(2)),
        ],
        _ => vec![x.pow(n).sub(&y.pow(n))],
    }
}

/// Legacy wrapper used by `factor_poly`.
pub fn factor_xn_minus_one_display(p: &Poly) -> Option<Poly> {
    let vars = vars_in(p);
    if vars.len() != 1 {
        return None;
    }
    let v = &vars[0];
    if is_xn_minus_one_poly(p, v, 4) {
        let x = Poly::var(v.clone());
        return Some(
            x.sub(&Poly::one())
                .mul(&x.add(&Poly::one()))
                .mul(&x.pow(2).add(&Poly::one())),
        );
    }
    None
}