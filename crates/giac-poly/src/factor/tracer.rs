//! Regression tests mapped from upstream `check/testfactor` (Issue 2.3/2.4).
//!
//! L20/L22 已绿。测试辅助为 Pipeline private。

use crate::monomial::Var;
use crate::poly::Poly;

use super::multivariate::factor_multivariate;

// **Pipeline private** — `x`
fn x() -> Poly {
    Poly::var("x")
}

// **Pipeline private** — `y`
fn y() -> Poly {
    Poly::var("y")
}

// **Pipeline private** — `z`
fn z() -> Poly {
    Poly::var("z")
}

// **Pipeline private** — `b_var`
fn b_var() -> Poly {
    Poly::var("b")
}

// **Pipeline private** — `c_var`
fn c_var() -> Poly {
    Poly::var("c")
}

// **Pipeline private** — `assert_factors`
fn assert_factors(p: &Poly, min_count: usize) {
    let f = factor_multivariate(p).unwrap_or_else(|e| panic!("factor failed: {e:?}"));
    assert!(
        f.len() >= min_count,
        "expected >= {min_count} factors, got {} for p={p:?}",
        f.len()
    );
    let prod = f.iter().fold(Poly::one(), |acc, q| acc.mul(q));
    assert_eq!(prod, *p, "product mismatch");
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    use num_traits::One;

    /// testfactor line 16: `(x-y+1)*(x-y)*(x-y-1)`
    #[test]
    fn testfactor_line16_three_shifted_linears() {
        let p = x()
            .sub(&y())
            .add(&Poly::one())
            .mul(&x().sub(&y()))
            .mul(&x().sub(&y()).sub(&Poly::one()));
        assert_factors(&p, 3);
    }

    /// testfactor line 17: `(x^2-3*x+1)*(x^2+x+1)`
    #[test]
    fn testfactor_line17_two_quadratics() {
        let q1 = x().pow(2).sub(&x().mul_scalar(&Ratio::from_integer(3.into()))).add(&Poly::one());
        let q2 = x().pow(2).add(&x()).add(&Poly::one());
        assert_factors(&q1.mul(&q2), 2);
    }

    /// testfactor line 20: `(x+b+c)*(x^2-x*b-x*c+b^2-b*c+c^2)`
    #[test]
    fn testfactor_line20_parametric_cubic_factor() {
        let p = x()
            .add(&b_var())
            .add(&c_var())
            .mul(
                &x()
                    .pow(2)
                    .sub(&x().mul(&b_var()))
                    .sub(&x().mul(&c_var()))
                    .add(&b_var().pow(2))
                    .sub(&b_var().mul(&c_var()))
                    .add(&c_var().pow(2)),
            );
        assert_factors(&p, 2);
    }

    /// testfactor line 21: `(x-y-z)*(x-y+z)*(x+y+z)`
    #[test]
    fn testfactor_line21_three_linear_ternary() {
        let p = x()
            .sub(&y())
            .sub(&z())
            .mul(&x().sub(&y()).add(&z()))
            .mul(&x().add(&y()).add(&z()));
        assert_factors(&p, 3);
    }

    /// testfactor line 22: `(3*x-y^2+y-5)*(x*y+3*x-y^2-1)`
    #[test]
    fn testfactor_line22_bivariate_mixed_degree() {
        let p = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x())
            .sub(&y().pow(2))
            .add(&y())
            .sub(&Poly::constant(Ratio::from_integer(5.into())))
            .mul(
                &x()
                    .mul(&y())
                    .add(&Poly::constant(Ratio::from_integer(3.into())).mul(&x()))
                    .sub(&y().pow(2))
                    .sub(&Poly::one()),
            );
        assert_factors(&p, 2);
    }

    /// testfactor line 12: `(x^3-x+1)*(x^3+x+1)` — cubics without rational roots
    #[test]
    fn testfactor_line12_two_cubics() {
        let c1 = x().pow(3).sub(&x()).add(&Poly::one());
        let c2 = x().pow(3).add(&x()).add(&Poly::one());
        assert_factors(&c1.mul(&c2), 2);
    }

    /// testfactor line 24: `x^6-y^6`
    #[test]
    fn testfactor_line24_x6_minus_y6() {
        assert_factors(&x().pow(6).sub(&y().pow(6)), 4);
    }

    /// Synthetic gate: bilinear product — unitaryfactor path (sparse/hensel also succeed).
    #[test]
    fn testfactor_unitary_bilinear_gate() {
        let p = x()
            .add(&y())
            .sub(&Poly::one())
            .mul(&x().sub(&y()).sub(&Poly::one()));
        assert_factors(&p, 2);
    }

    /// Hensel-fail gate (line 25 target): L22+y^3 — multi-point coeff interp (P2a).
    #[test]
    fn testfactor_line25_unitaryfactor_gate() {
        let f1 = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x())
            .sub(&y().pow(2))
            .add(&y())
            .sub(&Poly::constant(Ratio::from_integer(5.into())));
        let f2 = x()
            .mul(&y())
            .add(&Poly::constant(Ratio::from_integer(3.into())).mul(&x()))
            .sub(&y().pow(2))
            .sub(&Poly::one())
            .add(&y().pow(3));
        let p = f1.mul(&f2);
        assert_factors(&p, 2);
    }

    /// line 26: line25 f2 + (-x^2*y) cross term — harder mixed-degree gate.
    #[test]
    fn testfactor_line26_y3_x2y_gate() {
        let f1 = Poly::constant(Ratio::from_integer(3.into()))
            .mul(&x())
            .sub(&y().pow(2))
            .add(&y())
            .sub(&Poly::constant(Ratio::from_integer(5.into())));
        let f2 = x()
            .mul(&y())
            .add(&Poly::constant(Ratio::from_integer(3.into())).mul(&x()))
            .sub(&y().pow(2))
            .sub(&Poly::one())
            .add(&y().pow(3))
            .sub(&x().pow(2).mul(&y()));
        let p = f1.mul(&f2);
        assert_factors(&p, 2);
    }
}
