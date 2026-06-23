//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

use crate::error::{EvalError, PolyResult};
use crate::exp::bigint_pow;
use crate::monomial::Var;
use crate::poly::Poly;

/// Sylvester resultant of univariate polynomials in `var`.
/// **Stable** — univariate resultant
pub fn resultant(a: &Poly, b: &Poly, var: &Var) -> PolyResult<Poly> {
    let g = a.gcd(b);
    if univariate_degree(&g, var) > 0 {
        return Ok(Poly::zero());
    }
    let da = univariate_degree(a, var);
    let db = univariate_degree(b, var);
    if da == 0 || db == 0 {
        let c = if da == 0 {
            univariate_leading_coeff(a, var)
        } else {
            univariate_leading_coeff(b, var)
        };
        let exp = if da == 0 { db } else { da };
        let other = if da == 0 { b } else { a };
        let c_pow = bigint_pow(&c, exp).ok_or(EvalError::TypeError("exponent too large"))?;
        return Ok(other.pow(exp).mul_scalar(&Ratio::from_integer(c_pow)));
    }
    if da == 1 || db == 1 {
        return Ok(sylvester_det2(a, b, var));
    }
    let ac = univariate_coefficients(a, var, da);
    let bc = univariate_coefficients(b, var, db);
    let det = sylvester_det(&ac, da as usize, &bc, db as usize);
    Ok(Poly::constant(det))
}

/// **Stable** — ascending univariate coefficient vector in ℚ[var].
pub fn univariate_coeffs_ascending(p: &Poly, var: &Var) -> Vec<Ratio<BigInt>> {
    let deg = univariate_degree(p, var);
    (0..=deg).map(|e| coeff_at(p, var, e)).collect()
}

// **Pipeline private** — `univariate_coefficients`
fn univariate_coefficients(p: &Poly, var: &Var, degree: u64) -> Vec<Ratio<BigInt>> {
    let mut coeffs = univariate_coeffs_ascending(p, var);
    coeffs.resize((degree + 1) as usize, Ratio::zero());
    coeffs.truncate((degree + 1) as usize);
    coeffs
}

// **Pipeline private** — `sylvester_det`
fn sylvester_det(
    a: &[Ratio<BigInt>],
    deg_a: usize,
    b: &[Ratio<BigInt>],
    deg_b: usize,
) -> Ratio<BigInt> {
    let n = deg_b;
    let m = deg_a;
    let size = m + n;
    let mut mat = vec![vec![Ratio::zero(); size]; size];
    for i in 0..n {
        for j in 0..=m {
            mat[i][i + j] = a[m - j].clone();
        }
    }
    for j in 0..m {
        for k in 0..=n {
            mat[n + j][j + k] = b[n - k].clone();
        }
    }
    det_rational(&mut mat)
}

// **Pipeline private** — `det_rational`
fn det_rational(mat: &mut [Vec<Ratio<BigInt>>]) -> Ratio<BigInt> {
    let n = mat.len();
    let mut det = Ratio::one();
    for k in 0..n {
        let mut pivot_row = k;
        while pivot_row < n && mat[pivot_row][k].is_zero() {
            pivot_row += 1;
        }
        if pivot_row == n {
            return Ratio::zero();
        }
        if pivot_row != k {
            mat.swap(k, pivot_row);
            det = -det;
        }
        let pivot = mat[k][k].clone();
        det *= &pivot;
        for i in (k + 1)..n {
            if mat[i][k].is_zero() {
                continue;
            }
            let factor = mat[i][k].clone() / pivot.clone();
            for j in (k + 1)..n {
                let pivot_val = mat[k][j].clone();
                mat[i][j] -= &pivot_val * &factor;
            }
        }
    }
    det
}

// **Pipeline private** — `sylvester_det2`
fn sylvester_det2(a: &Poly, b: &Poly, var: &Var) -> Poly {
    let a1 = coeff_at(a, var, 1);
    let a0 = coeff_at(a, var, 0);
    let b1 = coeff_at(b, var, 1);
    let b0 = coeff_at(b, var, 0);
    Poly::constant(a1 * b0 - a0 * b1)
}

/// **Stable** — univariate coefficient at exponent
pub fn coeff_at(p: &Poly, var: &Var, exp: u64) -> Ratio<BigInt> {
    for (m, c) in &p.terms {
        if exp == 0 && m.is_const() {
            return c.clone();
        }
        if m.exp_of(var) == exp && m.iter().all(|(v, _)| v == var) {
            return c.clone();
        }
    }
    Ratio::zero()
}

/// **Stable** — degree w.r.t. var
pub fn univariate_degree(p: &Poly, var: &Var) -> u64 {
    p.degree_wrt(var)
}

// **Pipeline private** — `univariate_leading_coeff`
fn univariate_leading_coeff(p: &Poly, var: &Var) -> BigInt {
    let deg = univariate_degree(p, var);
    for (m, c) in p.terms.iter().rev() {
        if m.exp_of(var) == deg && c.denom().is_one() {
            return c.numer().clone();
        }
    }
    BigInt::one()
}

/// Roots of univariate polynomial (low-degree exact cases).
/// **Stable (bounded)** — low-degree exact roots as Poly factors
pub fn roots(p: &Poly, var: &Var) -> PolyResult<Vec<Poly>> {
    let d = univariate_degree(p, var);
    match d {
        0 => {
            if p.is_zero() {
                Ok(vec![])
            } else {
                Err(EvalError::TypeError("constant has no roots"))
            }
        }
        1 => {
            let mut a = Ratio::zero();
            let mut b = Ratio::zero();
            for (m, c) in &p.terms {
                let exp = m.exp_of(var);
                if exp == 1 {
                    a += c;
                } else if exp == 0 {
                    b += c;
                }
            }
            if a.is_zero() {
                return Err(EvalError::TypeError("not linear"));
            }
            Ok(vec![Poly::constant(-&b / &a)])
        }
        2 => quadratic_roots(p, var),
        3 if is_xn_minus_one(p, var, 3) => Ok(vec![Poly::constant(Ratio::one())]),
        _ => Err(EvalError::NotImplemented("roots")),
    }
}

pub use crate::quadratic::quadratic_abc;

// **Pipeline private** — `quadratic_roots`
fn quadratic_roots(p: &Poly, var: &Var) -> PolyResult<Vec<Poly>> {
    let q = crate::quadratic::quadratic_coeffs(p, var)
        .ok_or(EvalError::TypeError("not quadratic"))?;
    Ok(crate::quadratic::quadratic_rational_roots(&q)?
        .into_iter()
        .map(Poly::constant)
        .collect())
}

// **Pipeline private** — `is_xn_minus_one`
fn is_xn_minus_one(p: &Poly, var: &Var, n: u64) -> bool {
    if p.terms.len() != 2 {
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

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn resultant_shared_factor_is_zero() {
        let a = x().pow(2).sub(&Poly::one());
        let b = x().pow(3).sub(&Poly::one());
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        assert!(r.is_zero());
    }

    #[test]
    fn resultant_two_linear_polys() {
        let a = x().sub(&Poly::one());
        let b = x().add(&Poly::one());
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        assert_eq!(r, Poly::constant(Ratio::from_integer(BigInt::from(2))));
    }

    #[test]
    fn resultant_constant_times_linear() {
        let a = Poly::constant(Ratio::from_integer(BigInt::from(5)));
        let b = x().sub(&Poly::one());
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        let expected = b.mul_scalar(&Ratio::from_integer(BigInt::from(5)));
        assert_eq!(r, expected);
    }

    #[test]
    fn resultant_linear_times_constant() {
        let a = x().sub(&Poly::one());
        let b = Poly::constant(Ratio::from_integer(BigInt::from(7)));
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        let expected = a.mul_scalar(&Ratio::from_integer(BigInt::from(7)));
        assert_eq!(r, expected);
    }

    #[test]
    fn resultant_quadratic() {
        let a = x().pow(2).add(&Poly::one());
        let b = x().pow(2).sub(&Poly::one());
        let r = resultant(&a, &b, &Var::from("x")).unwrap();
        assert_eq!(r, Poly::constant(Ratio::from_integer(BigInt::from(4))));
    }

    #[test]
    fn quadratic_abc_x2_minus_2() {
        let x = Var::from("x");
        let mut p = Poly::var(x.clone()).pow(2);
        p = p.sub(&Poly::constant(Ratio::from_integer(BigInt::from(2))));
        let (a, b, c) = quadratic_abc(&p, &x).unwrap();
        assert_eq!(a, Ratio::one());
        assert!(b.is_zero());
        assert_eq!(c, Ratio::from_integer(BigInt::from(-2)));
    }

    #[test]
    fn roots_linear() {
        let p = x().sub(&Poly::one());
        let rs = roots(&p, &Var::from("x")).unwrap();
        assert_eq!(rs.len(), 1);
        assert_eq!(rs[0], Poly::constant(Ratio::one()));
    }

    #[test]
    fn roots_zero_polynomial() {
        let rs = roots(&Poly::zero(), &Var::from("x")).unwrap();
        assert!(rs.is_empty());
    }

    #[test]
    fn roots_constant_nonzero_errors() {
        let err = roots(&Poly::one(), &Var::from("x")).unwrap_err();
        assert!(matches!(err, EvalError::TypeError(_)));
    }

    #[test]
    fn roots_x3_minus_one_real_root() {
        let p = x().pow(3).sub(&Poly::one());
        let rs = roots(&p, &Var::from("x")).unwrap();
        assert_eq!(rs, vec![Poly::constant(Ratio::one())]);
    }

    #[test]
    fn roots_quadratic_perfect_square() {
        let p = x().pow(2).sub(&x().mul_scalar(&Ratio::from_integer(BigInt::from(2)))).add(&Poly::one());
        let rs = roots(&p, &Var::from("x")).unwrap();
        assert_eq!(rs.len(), 1);
        assert_eq!(rs[0], Poly::constant(Ratio::one()));
    }

    #[test]
    fn roots_quadratic_irrational_discriminant() {
        let p = x().pow(2).add(&Poly::one());
        let err = roots(&p, &Var::from("x")).unwrap_err();
        assert!(matches!(err, EvalError::NotImplemented(_)));
    }
}
