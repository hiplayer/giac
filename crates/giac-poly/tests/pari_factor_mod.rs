//! Compare `factor_poly_mod` against PARI `factor(Mod(...))` when gp is available.

use giac_poly::{factor_mod_irreducibles, factor_poly_mod, modp, Poly};

fn gp_factor_count(coeffs: &[(i64, u64)], prime: i64) -> Option<usize> {
    let gp = "/home/kanli.hu/upstream/pari/gp";
    if !std::path::Path::new(gp).exists() {
        return None;
    }
    let mut terms = String::new();
    for (c, e) in coeffs {
        if *c == 0 {
            continue;
        }
        if !terms.is_empty() {
            terms.push_str(" + ");
        }
        terms.push_str(&format!("Mod({c},{prime})*x^{e}"));
    }
    if terms.is_empty() {
        terms = format!("Mod(0,{prime})");
    }
    let output = std::process::Command::new(gp)
        .arg("-q")
        .arg("-s")
        .arg(format!("print(#factor({terms}));"))
        .output()
        .ok()?;
    let s = String::from_utf8_lossy(&output.stdout);
    s.trim().parse().ok()
}

fn poly_from_coeffs(coeffs: &[(i64, u64)]) -> Poly {
    let x = Poly::var("x");
    let mut p = Poly::zero();
    for (c, e) in coeffs {
        if *c == 0 {
            continue;
        }
        let term = if *e == 0 {
            Poly::constant(num_rational::Ratio::from_integer((*c).into()))
        } else {
            Poly::constant(num_rational::Ratio::from_integer((*c).into())).mul(&x.pow(*e))
        };
        p = p.add(&term);
    }
    p
}

#[test]
fn pari_matches_giac_factor_mod() {
    let cases: &[(&[(i64, u64)], i64)] = &[
        (&[(1, 4), (1, 0)], 5),
        (&[(1, 6), (-1, 0)], 7),
        (&[(1, 2), (1, 0)], 5),
        (&[(1, 3), (1, 1), (1, 0)], 7),
        (&[(1, 4), (1, 3), (1, 2), (1, 1), (1, 0)], 11),
    ];
    for (coeffs, p) in cases {
        let poly = poly_from_coeffs(coeffs);
        let giac_count = factor_mod_irreducibles(&poly, *p)
            .expect("factor_mod_irreducibles")
            .len();
        if let Some(pari_count) = gp_factor_count(coeffs, *p) {
            assert_eq!(
                giac_count, pari_count,
                "factor count mismatch for poly={poly:?} mod {p}"
            );
        }
        let prod = factor_poly_mod(&poly, *p).unwrap();
        assert_eq!(modp(&prod, *p).unwrap(), modp(&poly, *p).unwrap());
    }
}
