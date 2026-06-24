//! Integer factorization (`ifactor`).

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-simplify-api-stability.md`.
//!
//!
use std::sync::Arc;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use giac_core::{Expr, ExprArc};

const SMALL_PRIMES: &[u64] = &[
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
];

/// **Stable** — factor `n` into primes: `p1^e1 * p2^e2 * ...`.
pub fn ifactor(n: &BigInt) -> ExprArc {
    if n.is_zero() {
        return Expr::int(0);
    }
    if n.is_one() || n == &-BigInt::from(1) {
        return Arc::new(Expr::Int(n.clone()));
    }
    let negative = n.is_negative();
    let mut n = n.abs();
    let mut factors: Vec<(BigInt, u32)> = Vec::new();

    for &p in SMALL_PRIMES {
        let p = BigInt::from(p);
        let mut e = 0u32;
        while (&n % &p).is_zero() {
            n /= &p;
            e += 1;
        }
        if e > 0 {
            factors.push((p, e));
        }
        if n.is_one() {
            break;
        }
    }

    while !n.is_one() {
        if n.bits() as usize <= 64 {
            let ni = n.to_u64_digits().1.first().copied().unwrap_or(0);
            if is_probable_prime_u64(ni) {
                factors.push((n, 1));
                break;
            }
        }
        let Some(d) = pollard_rho(&n) else {
            factors.push((n, 1));
            break;
        };
        if d.is_one() || d == n {
            factors.push((n, 1));
            break;
        }
        let mut e = 0u32;
        while (&n % &d).is_zero() {
            n /= &d;
            e += 1;
        }
        factors.push((d, e));
    }

    let mut parts = Vec::new();
    if negative {
        parts.push(Expr::int(-1));
    }
    for (p, e) in factors {
        if e == 1 {
            parts.push(Arc::new(Expr::Int(p)));
        } else {
            parts.push(Expr::pow(Arc::new(Expr::Int(p)), Expr::int(i64::from(e))));
        }
    }
    if parts.is_empty() {
        Expr::int(1)
    } else if parts.len() == 1 {
        parts.remove(0)
    } else {
        Expr::mul(parts)
    }
}

// **Pipeline private** — Miller-Rabin for u64
fn is_probable_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &p in SMALL_PRIMES {
        if n == p {
            return true;
        }
        if n % p == 0 {
            return false;
        }
    }
    let mut d = n - 1;
    let mut s = 0u32;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    for a in [2u64, 3, 5, 7, 11] {
        if a >= n {
            continue;
        }
        let mut x = mod_pow(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut ok = false;
        for _ in 1..s {
            x = mod_mul(x, x, n);
            if x == n - 1 {
                ok = true;
                break;
            }
        }
        if !ok {
            return false;
        }
    }
    true
}

// **Pipeline private** — u64 modular multiply
fn mod_mul(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

// **Pipeline private** — u64 modular exponentiation
fn mod_pow(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let mut out = 1u64;
    base %= m;
    while exp > 0 {
        if exp % 2 == 1 {
            out = mod_mul(out, base, m);
        }
        base = mod_mul(base, base, m);
        exp /= 2;
    }
    out
}

// **Pipeline private** — Pollard rho on BigInt
fn pollard_rho(n: &BigInt) -> Option<BigInt> {
    if n.bits() as usize <= 64 {
        let ni = n.to_u64_digits().1.first().copied().unwrap_or(0);
        if ni < 4 {
            return None;
        }
        return pollard_rho_u64(ni).map(BigInt::from);
    }
    let mut x = BigInt::from(2);
    let mut y = BigInt::from(2);
    let mut c = BigInt::from(1);
    let mut d = BigInt::one();
    while d.is_one() {
        x = (&x * &x + &c) % n;
        y = (&(&y * &y + &c) % n * &y + &c) % n;
        d = (if x > y { &x - &y } else { &y - &x }).gcd(n);
        if d == *n {
            c += 1;
            x = BigInt::from(2);
            y = BigInt::from(2);
            d = BigInt::one();
        }
    }
    if d.is_one() || d == *n {
        None
    } else {
        Some(d)
    }
}

// **Pipeline private** — Pollard rho on u64
fn pollard_rho_u64(n: u64) -> Option<u64> {
    if n % 2 == 0 {
        return Some(2);
    }
    let mut x = 2u64;
    let mut y = 2u64;
    let mut c = 1u64;
    let mut d = 1u64;
    while d == 1 {
        x = (mod_mul(x, x, n) + c) % n;
        y = (mod_mul(y, y, n) + c) % n;
        y = (mod_mul(y, y, n) + c) % n;
        d = x.abs_diff(y).gcd(&n);
        if d == n {
            c += 1;
            x = 2;
            y = 2;
            d = 1;
        }
    }
    if d == 1 || d == n {
        None
    } else {
        Some(d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::format_expr;

    #[test]
    fn ifactor_small() {
        assert_eq!(format_expr(ifactor(&BigInt::from(360)).as_ref()), "2^3*3^2*5");
        assert_eq!(format_expr(ifactor(&BigInt::from(17)).as_ref()), "17");
    }
}
