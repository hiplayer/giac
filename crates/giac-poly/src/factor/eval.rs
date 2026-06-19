//! Good evaluation points for multivariate factorization.
//!
//! **Upstream:** `ezgcd.cc` `find_good_eval`, `peval_1`; used in `do_factor_hensel`
//! for irreducibility probes and Hensel/sparse seeds.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
pub use super::ctx::GoodEval;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Zero;

use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::univariate_degree;

use super::poly_uni::substitute_poly;
use super::univariate::factor_univariate_flat;

/// **Pipeline private** — substitute auxiliary vars with scalars; keep `main` univariate.
pub fn peval_at_main(p: &Poly, main: &Var, auxes: &[&Var], values: &[Ratio<BigInt>]) -> Poly {
    debug_assert_eq!(auxes.len(), values.len());
    let mut out = p.clone();
    for (v, c) in auxes.iter().zip(values.iter()) {
        out = substitute_poly(&out, v, &Poly::constant(c.clone()));
    }
    out
}

/// **Pipeline private** — find evaluation preserving `main`-degree (upstream `find_good_eval`).
pub fn find_good_eval(
    p: &Poly,
    main: &Var,
    auxes: &[&Var],
    start: &[i64],
) -> Option<GoodEval> {
    let target_deg = univariate_degree(p, main);
    if target_deg == 0 {
        return None;
    }
    for pts in eval_point_candidates(auxes.len(), start) {
        let vals: Vec<Ratio<BigInt>> = pts
            .iter()
            .map(|k| Ratio::from_integer(BigInt::from(*k)))
            .collect();
        let ev = peval_at_main(p, main, auxes, &vals);
        if univariate_degree(&ev, main) == target_deg && !ev.is_zero() {
            return Some(GoodEval::new(ev, vals, target_deg));
        }
    }
    None
}

// **Pipeline private** — trial points: `start` first, then small integers / shifts.
fn eval_point_candidates(n_aux: usize, start: &[i64]) -> Vec<Vec<i64>> {
    if n_aux == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut push = |v: Vec<i64>| {
        if v.len() == n_aux && !out.iter().any(|u| u == &v) {
            out.push(v);
        }
    };
    if start.len() == n_aux {
        push(start.to_vec());
    } else {
        push(vec![0; n_aux]);
    }
    if n_aux == 1 {
        for k in 1i64..=12 {
            push(vec![k]);
            push(vec![-k]);
        }
    } else {
        push(vec![1; n_aux]);
        push(vec![-1; n_aux]);
        for i in 0..n_aux {
            let mut v = vec![0; n_aux];
            v[i] = 1;
            push(v.clone());
            v[i] = 2;
            push(v);
        }
        push((1..=n_aux as i64).collect());
        push((1..=n_aux as i64).map(|k| -k).collect());
    }
    out
}

/// **Pipeline private** — upstream `do_factor_hensel`: two good evals, single factor → irreducible.
pub fn looks_irreducible_by_good_eval(p: &Poly, main: &Var, auxes: &[&Var]) -> bool {
    let zero = vec![0i64; auxes.len()];
    let one = vec![1i64; auxes.len()];
    for start in [&zero[..], &one[..]] {
        let ge = match find_good_eval(p, main, auxes, start) {
            Some(x) => x,
            None => return false,
        };
        let facs = match factor_univariate_flat(&ge.evaluated, main) {
            Ok(f) => f,
            Err(_) => return false,
        };
        let nontrivial = facs
            .iter()
            .filter(|f| univariate_degree(f, main) > 0)
            .count();
        if nontrivial == 1 {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    #[test]
    fn find_good_eval_preserves_degree() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.pow(2).sub(&y).mul(&x.add(&Poly::one()));
        let ge = find_good_eval(&p, &Var::from("x"), &[&Var::from("y")], &[0]).unwrap();
        assert_eq!(univariate_degree(&ge.evaluated, &Var::from("x")), 3);
        assert_eq!(ge.preserved_main_degree, 3);
        assert!(!ge.evaluated.is_zero());
    }

    #[test]
    fn find_good_eval_skips_bad_start() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        // x*y: y=0 kills x-degree; later candidates still find a good point.
        let p = x.mul(&y);
        assert!(find_good_eval(&p, &Var::from("x"), &[&Var::from("y")], &[0]).is_some());
        let ge =
            find_good_eval(&p, &Var::from("x"), &[&Var::from("y")], &[0]).unwrap();
        assert_ne!(ge.values[0], Ratio::zero());
        assert_eq!(univariate_degree(&ge.evaluated, &Var::from("x")), 1);
    }

    #[test]
    fn irreducibility_probe_detects_x2_plus_y2_plus_1() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.pow(2).add(&y.pow(2)).add(&Poly::one());
        assert!(looks_irreducible_by_good_eval(
            &p,
            &Var::from("x"),
            &[&Var::from("y")]
        ));
    }

    #[test]
    fn irreducibility_probe_rejects_reducible() {
        let x = Poly::var("x");
        let y = Poly::var("y");
        let p = x.pow(2).sub(&y.pow(2));
        assert!(!looks_irreducible_by_good_eval(
            &p,
            &Var::from("x"),
            &[&Var::from("y")]
        ));
    }
}
