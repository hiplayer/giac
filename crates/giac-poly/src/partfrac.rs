//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use num_bigint::BigInt;
use num_rational::Ratio;

use num_traits::{One, Zero};

use crate::error::{EvalError, PolyResult};
use crate::factor::{
    factor_into, factor_power_pairs,
};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};
use crate::univariate::square_free_factorization;
/// Partial fraction terms `(coeff, denominator factor)` for `num/den` in `var`.
/// **Stable (bounded)** — partial fraction terms
pub fn partfrac_terms(
    num: &Poly,
    den: &Poly,
    var: &Var,
) -> PolyResult<(Option<Poly>, Vec<(Ratio<BigInt>, Poly)>)> {
    let (poly_part, terms) = partfrac_rational_terms(num, den, var)?;
    let mut linear = Vec::new();
    for (n, d) in terms {
        if univariate_degree(&n, var) == 0 && univariate_degree(&d, var) >= 1 {
            linear.push((coeff_at(&n, var, 0), d));
        } else {
            return Err(EvalError::NotImplemented("partfrac nonlinear factor"));
        }
    }
    Ok((poly_part, linear))
}

/// Partial fractions with polynomial numerators `(numer, denom_factor)`.
/// **Stable (bounded)** — partfrac with poly part
pub fn partfrac_rational_terms(
    num: &Poly,
    den: &Poly,
    var: &Var,
) -> PolyResult<(Option<Poly>, Vec<(Poly, Poly)>)> {
    if den.is_zero() {
        return Err(EvalError::DivisionByZero);
    }
    let (poly_part, rem) = if univariate_degree(num, var) >= univariate_degree(den, var) {
        let (q, r) = num.div_rem(den);
        (Some(q), r)
    } else {
        (None, num.clone())
    };
    if rem.is_zero() {
        return Ok((poly_part, vec![]));
    }
    if univariate_degree(&rem, var) >= univariate_degree(den, var) {
        return Err(EvalError::TypeError("improper rational remainder"));
    }

    if let Ok(terms) = partfrac_by_square_free(&rem, den, var) {
        return Ok((poly_part, drop_zero_numerators(terms)));
    }

    let factors = factor_into(den).ok_or(EvalError::NotImplemented("partfrac factor"))?;
    if factors.is_empty() {
        return Err(EvalError::TypeError("empty factorization"));
    }
    if factors.iter().all(|f| univariate_degree(f, var) == 1) {
        let mut terms = Vec::with_capacity(factors.len());
        for f in &factors {
            let root = linear_root(f, var)?;
            let mut denom_prod = Ratio::one();
            for g in &factors {
                if g != f {
                    denom_prod *= g.horner(var, &root);
                }
            }
            if denom_prod.is_zero() {
                return Err(EvalError::TypeError("repeated linear factor"));
            }
            let coeff = rem.horner(var, &root) / denom_prod;
            terms.push((Poly::constant(coeff), f.clone()));
        }
        return Ok((poly_part, drop_zero_numerators(terms)));
    }
    let terms = partfrac_mixed_affine(&rem, &factors, var)?;
    Ok((poly_part, drop_zero_numerators(terms)))
}

// **Pipeline private** — `drop_zero_numerators`
fn drop_zero_numerators(terms: Vec<(Poly, Poly)>) -> Vec<(Poly, Poly)> {
    terms.into_iter().filter(|(n, _)| !n.is_zero()).collect()
}

/// Partial fractions via square-free factorization (GIAC-224).
// **Pipeline private** — `partfrac_by_square_free`
fn partfrac_by_square_free(
    num: &Poly,
    den: &Poly,
    var: &Var,
) -> PolyResult<Vec<(Poly, Poly)>> {
    let sqff = expand_sqff_factors(&denominator_power_factors(den, var)?, var)?;
    if sqff.is_empty() {
        return Err(EvalError::TypeError("empty denominator"));
    }

    let mut denom_powers = Vec::new();
    for (g, mult) in &sqff {
        if univariate_degree(g, var) > 3 {
            return Err(EvalError::NotImplemented("partfrac nonlinear factor"));
        }
        for j in 1..=*mult {
            denom_powers.push(g.pow(j as u64));
        }
    }

    if denom_powers.len() == 1 && univariate_degree(&denom_powers[0], var) <= 2 {
        let g = &denom_powers[0];
        if univariate_degree(g, var) == 1 {
            let root = linear_root(g, var)?;
            let slope = coeff_at(g, var, 1);
            if slope.is_zero() {
                return Err(EvalError::TypeError("degenerate linear factor"));
            }
            let coeff = num.horner(var, &root) / slope;
            return Ok(vec![(Poly::constant(coeff), g.clone())]);
        }
        if univariate_degree(g, var) == 2 && univariate_degree(num, var) == 0 {
            return partfrac_one_quadratic(num, g, var);
        }
    }

    if sqff.iter().all(|(_, m)| *m == 1)
        && sqff
            .iter()
            .any(|(g, _)| univariate_degree(g, var) == 2)
    {
        return partfrac_square_free_affine_numerators(num, den, var, &sqff);
    }

    partfrac_affine_power_system(num, den, var, &sqff)
}

/// Partial fractions with numerators up to `deg(g)-1` for each `g^j` term.
// **Pipeline private** — `partfrac_affine_power_system`
fn partfrac_affine_power_system(
    num: &Poly,
    den: &Poly,
    var: &Var,
    sqff: &[(Poly, usize)],
) -> PolyResult<Vec<(Poly, Poly)>> {
    let mut term_specs = Vec::new();
    for (g, mult) in sqff {
        let gdeg = univariate_degree(g, var) as usize;
        if gdeg == 0 || gdeg > 3 {
            return Err(EvalError::NotImplemented("partfrac nonlinear factor"));
        }
        for j in 1..=*mult {
            term_specs.push((g.pow(j as u64), gdeg - 1));
        }
    }
    let n_unknowns: usize = term_specs.iter().map(|(_, d)| d + 1).sum();
    let nden = univariate_degree(den, var) as usize;
    let mut matrix = vec![vec![Ratio::zero(); n_unknowns]; nden];
    let mut rhs = vec![Ratio::zero(); nden];
    let mut col = 0usize;
    for (d_k, max_pow) in &term_specs {
        let cofactor = den.div_rem(d_k).0;
        for k in 0..=*max_pow {
            let scaled = cofactor.mul(&Poly::var(var.clone()).pow(k as u64));
            for i in 0..nden {
                matrix[i][col] = coeff_at(&scaled, var, i as u64);
            }
            col += 1;
        }
    }
    for i in 0..nden {
        rhs[i] = coeff_at(num, var, i as u64);
    }
    let coeffs = solve_linear_system(&matrix, &rhs)
        .ok_or(EvalError::NotImplemented("partfrac linear system"))?;
    let mut col = 0usize;
    let mut out = Vec::new();
    for (d_k, max_pow) in term_specs {
        let mut numer = Poly::zero();
        for k in 0..=max_pow {
            numer = numer.add(
                &Poly::var(var.clone())
                    .pow(k as u64)
                    .mul_scalar(&coeffs[col]),
            );
            col += 1;
        }
        out.push((numer, d_k));
    }
    Ok(out)
}

/// Square-free factors with multiplicity; rational roots when Yun sqff stalls.
// **Pipeline private** — `denominator_power_factors`
fn denominator_power_factors(den: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    if let Ok(factors) = square_free_factorization(den, var) {
        let deg_ok = factors
            .iter()
            .all(|(g, _)| univariate_degree(g, var) <= univariate_degree(den, var).saturating_sub(1).max(1));
        let split = factors.len() > 1
            || factors
                .iter()
                .any(|(g, _)| univariate_degree(g, var) < univariate_degree(den, var));
        if deg_ok && split {
            return Ok(factors);
        }
    }
    factor_power_pairs(den, var)
}

/// Split square-free factors of degree > 2 into linear/quadratic pieces.
// **Pipeline private** — `expand_sqff_factors`
fn expand_sqff_factors(
    sqff: &[(Poly, usize)],
    var: &Var,
) -> PolyResult<Vec<(Poly, usize)>> {
    let mut out = Vec::new();
    for (g, mult) in sqff {
        if univariate_degree(g, var) <= 1 {
            out.push((g.clone(), *mult));
            continue;
        }
        if let Some(factors) = factor_into(g) {
            if factors.len() > 1
                || univariate_degree(&factors[0], var) < univariate_degree(g, var)
            {
                for f in factors {
                    out.push((f, *mult));
                }
                continue;
            }
        }
        if univariate_degree(g, var) <= 3 {
            out.push((g.clone(), *mult));
            continue;
        }
        if let Some(factors) = factor_into(g) {
            for f in factors {
                out.push((f, *mult));
            }
            continue;
        }
        let sub = factor_power_pairs(g, var)?;
        for (f, m) in sub {
            out.push((f, *mult * m));
        }
    }
    Ok(out)
}

// **Pipeline private** — `partfrac_square_free_affine_numerators`
fn partfrac_square_free_affine_numerators(
    num: &Poly,
    den: &Poly,
    var: &Var,
    sqff: &[(Poly, usize)],
) -> PolyResult<Vec<(Poly, Poly)>> {
    let mut term_specs = Vec::new();
    for (g, mult) in sqff {
        if *mult != 1 || univariate_degree(g, var) > 3 {
            return Err(EvalError::NotImplemented("partfrac nonlinear factor"));
        }
        term_specs.push((g.clone(), univariate_degree(g, var)));
    }
    let n_unknowns: usize = term_specs
        .iter()
        .map(|(_, d)| *d as usize)
        .sum();
    let nden = univariate_degree(den, var) as usize;
    let mut matrix = vec![vec![Ratio::zero(); n_unknowns]; nden];
    let mut rhs = vec![Ratio::zero(); nden];
    let mut col = 0usize;
    for (g, g_deg) in &term_specs {
        let d_k = g.clone();
        let cofactor = den.div_rem(&d_k).0;
        for j in 0..*g_deg {
            let t_pow = if j == 0 {
                Poly::one()
            } else {
                Poly::var(var.clone()).pow(j)
            };
            let scaled = cofactor.mul(&t_pow);
            for i in 0..nden {
                matrix[i][col] = coeff_at(&scaled, var, i as u64);
            }
            col += 1;
        }
    }
    for i in 0..nden {
        rhs[i] = coeff_at(num, var, i as u64);
    }
    let coeffs = solve_linear_system(&matrix, &rhs)
        .ok_or(EvalError::NotImplemented("partfrac linear system"))?;
    let mut col = 0usize;
    let mut out = Vec::new();
    for (g, g_deg) in term_specs {
        let mut numer = Poly::zero();
        for j in 0..g_deg {
            numer = numer.add(&Poly::var(var.clone()).pow(j).mul_scalar(&coeffs[col]));
            col += 1;
        }
        out.push((numer, g));
    }
    Ok(out)
}

// **Pipeline private** — `partfrac_one_quadratic`
fn partfrac_one_quadratic(num: &Poly, quad: &Poly, var: &Var) -> PolyResult<Vec<(Poly, Poly)>> {
    let c = coeff_at(num, var, 0);
    let a = coeff_at(quad, var, 2);
    let b = coeff_at(quad, var, 1);
    let d = coeff_at(quad, var, 0);
    if a.is_zero() {
        return Err(EvalError::TypeError("not quadratic"));
    }
    let disc = b.clone() * b.clone() - Ratio::from_integer(BigInt::from(4)) * a.clone() * d;
    if disc > Ratio::zero() {
        if let Some(factors) = factor_into(quad) {
            if factors.len() > 1 {
                return partfrac_by_square_free(num, quad, var);
            }
        }
        return Err(EvalError::NotImplemented("partfrac real quadratic split"));
    }
    if disc == Ratio::zero() {
        if let Some(factors) = factor_into(quad) {
            if factors.len() == 1 && univariate_degree(&factors[0], var) == 1 {
                return partfrac_affine_power_system(num, quad, var, &[(factors[0].clone(), 2)]);
            }
        }
        return partfrac_affine_power_system(num, quad, var, &[(quad.clone(), 2)]);
    }
    let quad_numer = Poly::constant(c / a);
    Ok(vec![(quad_numer, quad.clone())])
}

// **Pipeline private** — `solve_linear_system`
fn solve_linear_system(
    matrix: &[Vec<Ratio<BigInt>>],
    rhs: &[Ratio<BigInt>],
) -> Option<Vec<Ratio<BigInt>>> {
    let n = rhs.len();
    if matrix.is_empty() || matrix[0].len() != n {
        return None;
    }
    let mut a: Vec<Vec<Ratio<BigInt>>> = matrix.to_vec();
    let mut b = rhs.to_vec();
    let mut pivot_row = 0usize;
    for col in 0..n {
        let mut pivot = None;
        for row in pivot_row..n {
            if !a[row][col].is_zero() {
                pivot = Some(row);
                break;
            }
        }
        let p = pivot?;
        a.swap(pivot_row, p);
        b.swap(pivot_row, p);
        let pivot_val = a[pivot_row][col].clone();
        for j in col..n {
            a[pivot_row][j] /= &pivot_val;
        }
        b[pivot_row] /= &pivot_val;
        for row in 0..n {
            if row == pivot_row || a[row][col].is_zero() {
                continue;
            }
            let factor = a[row][col].clone();
            let pivot_vals: Vec<_> = a[pivot_row].to_vec();
            let pivot_b = b[pivot_row].clone();
            for j in col..n {
                a[row][j] -= &factor * &pivot_vals[j];
            }
            b[row] -= factor * pivot_b;
        }
        pivot_row += 1;
    }
    Some(b)
}

// **Pipeline private** — `partfrac_mixed_affine`
fn partfrac_mixed_affine(
    num: &Poly,
    factors: &[Poly],
    var: &Var,
) -> PolyResult<Vec<(Poly, Poly)>> {
    let den = factors.iter().fold(Poly::one(), |acc, f| acc.mul(f));
    let sqff: Vec<(Poly, usize)> = factors.iter().map(|f| (f.clone(), 1)).collect();
    partfrac_square_free_affine_numerators(num, &den, var, &sqff)
}

// **Pipeline private** — `linear_root`
fn linear_root(f: &Poly, var: &Var) -> PolyResult<Ratio<BigInt>> {
    let a = coeff_at(f, var, 1);
    if a.is_zero() {
        return Err(EvalError::TypeError("not linear"));
    }
    Ok(-coeff_at(f, var, 0) / a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::factor::{factor_into_by_rational_roots, find_rational_root};

    fn x() -> Poly {
        Poly::var("x")
    }

    #[test]
    fn find_rational_root_on_x_plus_one_times_x_fourth_minus_one() {
        let den = x().add(&Poly::one()).mul(&x().pow(4).sub(&Poly::one()));
        let r = find_rational_root(&den, &Var::from("x"));
        assert!(r.is_some());
    }

    #[test]
    fn find_rational_root_on_x_plus_one_sq_times_x_sq_plus_one() {
        let den = x()
            .add(&Poly::one())
            .pow(2)
            .mul(&x().pow(2).add(&Poly::one()));
        let r = find_rational_root(&den, &Var::from("x"));
        assert_eq!(r, Some(Ratio::from_integer((-1).into())));
    }

    #[test]
    fn factor_by_roots_x_plus_one_times_x_fourth_minus_one() {
        let den = x().add(&Poly::one()).mul(&x().pow(4).sub(&Poly::one()));
        let factors = factor_into_by_rational_roots(&den, &Var::from("x")).unwrap();
        assert!(factors.len() >= 3);
    }

    #[test]
    fn denominator_factors_x_plus_one_times_x_fourth_minus_one() {
        let den = x().add(&Poly::one()).mul(&x().pow(4).sub(&Poly::one()));
        let factors = denominator_power_factors(&den, &Var::from("x")).unwrap();
        assert!(factors.len() >= 2);
        let expanded = expand_sqff_factors(&factors, &Var::from("x")).unwrap();
        assert!(expanded.iter().all(|(g, _)| univariate_degree(g, &Var::from("x")) <= 2));
    }

    #[test]
    fn partfrac_cubic_irreducible_denominator() {
        let num = x();
        let den = x().pow(3).add(&Poly::constant(Ratio::from_integer(2.into())));
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 1);
        assert_eq!(univariate_degree(&terms[0].1, &Var::from("x")), 3);
        let rebuilt = terms.iter().fold(Poly::zero(), |acc, (n, d)| {
            acc.add(&n.mul(&den.div_rem(d).0))
        });
        assert_eq!(rebuilt, num);
    }

    #[test]
    fn partfrac_three_quarters_over_x_fourth_minus_one() {
        let num = Poly::constant(Ratio::new(3.into(), 4.into()));
        let den = x().pow(4).sub(&Poly::one());
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 3);
    }

    #[test]
    fn partfrac_by_sqff_x_over_x_plus_one_times_x_fourth_minus_one() {
        let num = x();
        let den = x().add(&Poly::one()).mul(&x().pow(4).sub(&Poly::one()));
        let terms = partfrac_by_square_free(&num, &den, &Var::from("x")).unwrap();
        assert!(terms.len() >= 3);
    }

    #[test]
    fn partfrac_x_over_x_plus_one_times_x_fourth_minus_one() {
        let num = x();
        let den = x().add(&Poly::one()).mul(&x().pow(4).sub(&Poly::one()));
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert!(terms.len() >= 3);
    }

    #[test]
    fn partfrac_one_over_x_squared_minus_one() {
        let num = Poly::one();
        let den = x().pow(2).sub(&Poly::one());
        let (_, terms) = partfrac_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 2);
        let sum: Ratio<BigInt> = terms.iter().map(|(c, _)| c).sum();
        assert!(sum.is_zero());
    }

    #[test]
    fn partfrac_x_over_repeated_linear() {
        let num = x();
        let den = x()
            .sub(&Poly::one())
            .mul(&x().add(&Poly::one()).pow(2));
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 3);
    }

    #[test]
    fn partfrac_one_over_x_squared_minus_one_squared() {
        let num = Poly::one();
        let den = x().pow(2).sub(&Poly::one()).pow(2);
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 4);
    }

    #[test]
    fn partfrac_biquadratic_half_angle_denominator() {
        let t = Poly::var("t");
        let den = t
            .pow(4)
            .mul_scalar(&Ratio::from_integer((-1).into()))
            .add(&t.pow(3).mul_scalar(&Ratio::from_integer(4.into())))
            .add(&t.pow(2).mul_scalar(&Ratio::from_integer((-2).into())))
            .add(&t.mul_scalar(&Ratio::from_integer(4.into())))
            .add(&Poly::constant(Ratio::from_integer((-1).into())));
        let num = t
            .pow(2)
            .mul_scalar(&Ratio::from_integer(2.into()))
            .add(&t.mul_scalar(&Ratio::from_integer(8.into())))
            .add(&Poly::constant(Ratio::from_integer(2.into())));
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("t")).unwrap();
        assert_eq!(terms.len(), 2);
        let recomposed = terms.iter().fold(Poly::zero(), |acc, (n, d)| {
            acc.add(&n.mul(&den.div_rem(d).0))
        });
        assert_eq!(recomposed, num);
    }

    #[test]
    fn partfrac_one_over_x_times_x_squared_plus_one() {
        let num = Poly::one();
        let den = x().mul(&x().pow(2).add(&Poly::one()));
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 2);
        let recomposed = terms.iter().fold(Poly::zero(), |acc, (n, d)| {
            acc.add(&n.mul(&den.div_rem(d).0))
        });
        assert_eq!(recomposed, num);
    }

    #[test]
    fn partfrac_xplus1_over_x_squared_minus_one() {
        let num = x().add(&Poly::one());
        let den = x().pow(2).sub(&Poly::one());
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        let recomposed = terms.iter().fold(Poly::zero(), |acc, (n, d)| {
            acc.add(&n.mul(&den.div_rem(d).0))
        });
        assert_eq!(recomposed, num);
    }

    #[test]
    fn partfrac_ck_int_05_denominator() {
        let three = Poly::constant(Ratio::from_integer(3.into()));
        let den = three
            .mul(&x())
            .mul(&x().pow(2).add(&x()).add(&Poly::one()))
            .mul(&x().sub(&Poly::one()).pow(3));
        let num = Poly::one();
        let r = partfrac_rational_terms(&num, &den, &Var::from("x"));
        eprintln!("ck05 partfrac: {:?}", r.as_ref().map(|(_, t)| t.len()));
        let (_, terms) = r.unwrap();
        assert!(terms.len() >= 4, "got {} terms", terms.len());
    }

    #[test]
    fn partfrac_one_over_x_minus_one_squared() {
        let num = Poly::one();
        let den = x().sub(&Poly::one()).pow(2);
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        let recomposed = terms.iter().fold(Poly::zero(), |acc, (n, d)| {
            acc.add(&n.mul(&den.div_rem(d).0))
        });
        assert_eq!(recomposed, num);
        assert!(!terms.is_empty());
    }

    #[test]
    fn partfrac_x_over_x_minus_one_squared() {
        let num = x();
        let den = x().sub(&Poly::one()).pow(2);
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 2);
        let recomposed = terms.iter().fold(Poly::zero(), |acc, (n, d)| {
            acc.add(&n.mul(&den.div_rem(d).0))
        });
        assert_eq!(recomposed, num);
    }

    #[test]
    fn partfrac_one_over_x_squared_minus_four() {
        let num = Poly::one();
        let den = x().pow(2).sub(&Poly::constant(Ratio::from_integer(4.into())));
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 2);
    }

    #[test]
    fn partfrac_one_over_one_minus_x_squared() {
        let num = Poly::one();
        let den = Poly::one().sub(&x().pow(2));
        let (_, terms) = partfrac_rational_terms(&num, &den, &Var::from("x")).unwrap();
        assert_eq!(terms.len(), 2);
    }
}
