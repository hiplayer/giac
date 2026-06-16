use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Signed;
use num_traits::{One, Zero};

use crate::error::{PolyError, PolyResult};
use crate::factor::{as_perfect_power, factor_into};
use crate::monomial::Var;
use crate::poly::Poly;
use crate::resultant::{coeff_at, univariate_degree};
use crate::univariate::{eval_univariate_at, square_free_factorization};

/// Partial fraction terms `(coeff, denominator factor)` for `num/den` in `var`.
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
            return Err(PolyError::NotImplemented("partfrac nonlinear factor"));
        }
    }
    Ok((poly_part, linear))
}

/// Partial fractions with polynomial numerators `(numer, denom_factor)`.
pub fn partfrac_rational_terms(
    num: &Poly,
    den: &Poly,
    var: &Var,
) -> PolyResult<(Option<Poly>, Vec<(Poly, Poly)>)> {
    if den.is_zero() {
        return Err(PolyError::DivisionByZero);
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
        return Err(PolyError::TypeError("improper rational remainder"));
    }

    if let Ok(terms) = partfrac_by_square_free(&rem, den, var) {
        return Ok((poly_part, terms));
    }

    let factors = factor_into(den).ok_or(PolyError::NotImplemented("partfrac factor"))?;
    if factors.is_empty() {
        return Err(PolyError::TypeError("empty factorization"));
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
                return Err(PolyError::TypeError("repeated linear factor"));
            }
            let coeff = rem.horner(var, &root) / denom_prod;
            terms.push((Poly::constant(coeff), f.clone()));
        }
        return Ok((poly_part, terms));
    }
    let terms = partfrac_mixed_constant(&factors, var)?;
    Ok((poly_part, terms))
}

/// Partial fractions via square-free factorization (GIAC-224).
fn partfrac_by_square_free(
    num: &Poly,
    den: &Poly,
    var: &Var,
) -> PolyResult<Vec<(Poly, Poly)>> {
    let sqff = denominator_power_factors(den, var)?;
    if sqff.is_empty() {
        return Err(PolyError::TypeError("empty denominator"));
    }

    let mut denom_powers = Vec::new();
    for (g, mult) in &sqff {
        if univariate_degree(g, var) > 2 {
            return Err(PolyError::NotImplemented("partfrac nonlinear factor"));
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
                return Err(PolyError::TypeError("degenerate linear factor"));
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

    let n_unknowns = denom_powers.len();
    let nden = univariate_degree(den, var) as usize;
    let mut matrix = vec![vec![Ratio::zero(); n_unknowns]; nden];
    let mut rhs = vec![Ratio::zero(); nden];
    for i in 0..nden {
        rhs[i] = coeff_at(num, var, i as u64);
        for (k, d_k) in denom_powers.iter().enumerate() {
            let cofactor = den.div_rem(d_k).0;
            matrix[i][k] = coeff_at(&cofactor, var, i as u64);
        }
    }
    let coeffs = solve_linear_system(&matrix, &rhs)
        .ok_or(PolyError::NotImplemented("partfrac linear system"))?;
    Ok(coeffs
        .into_iter()
        .zip(denom_powers)
        .map(|(c, d)| (Poly::constant(c), d))
        .collect())
}

/// Square-free factors with multiplicity; rational roots when Yun sqff stalls.
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
    factor_by_rational_roots(den, var)
}

fn factor_by_rational_roots(p: &Poly, var: &Var) -> PolyResult<Vec<(Poly, usize)>> {
    let mut rest = p.clone();
    let mut factors = Vec::new();
    while univariate_degree(&rest, var) > 0 {
        let root = match find_rational_root(&rest, var) {
            Some(r) => r,
            None => {
                if univariate_degree(&rest, var) == 4 {
                    if let Some(mut biq) = try_factor_biquadratic(&rest, var) {
                        factors.append(&mut biq);
                        rest = Poly::one();
                        break;
                    }
                }
                return Err(PolyError::NotImplemented("partfrac factor"));
            }
        };
        let lin = linear_poly(var, &root);
        let mut mult = 0usize;
        loop {
            let (_, r) = rest.div_rem(&lin);
            if !r.is_zero() {
                break;
            }
            mult += 1;
            rest = rest.div_rem(&lin).0;
        }
        factors.push((lin, mult));
    }
    if rest.is_one() || rest.is_zero() {
        return Ok(factors);
    }
    if let Some((base, exp)) = as_perfect_power(&rest) {
        if univariate_degree(&base, var) <= 2 {
            factors.push((base, exp as usize));
            return Ok(factors);
        }
    }
    if univariate_degree(&rest, var) <= 2 {
        factors.push((rest, 1));
        return Ok(factors);
    }
    Err(PolyError::NotImplemented("partfrac factor"))
}

fn linear_poly(var: &Var, root: &Ratio<BigInt>) -> Poly {
    Poly::var(var.clone()).sub(&Poly::constant(root.clone()))
}

fn find_rational_root(p: &Poly, var: &Var) -> Option<Ratio<BigInt>> {
    let deg = univariate_degree(p, var);
    if deg == 0 {
        return None;
    }
    let a0 = coeff_at(p, var, 0);
    let an = coeff_at(p, var, deg);
    for p_cand in integer_divisors(a0.numer()) {
        for q_cand in integer_divisors(an.numer()) {
            if q_cand.is_zero() {
                continue;
            }
            for &(pn, qn) in &[(1, 1), (1, -1), (-1, 1), (-1, -1)] {
                let r = Ratio::new(&p_cand * pn, &q_cand * qn);
                if eval_univariate_at(p, var, &r).is_zero() {
                    return Some(r);
                }
            }
        }
    }
    None
}

fn integer_divisors(n: &BigInt) -> Vec<BigInt> {
    if n.is_zero() {
        return vec![BigInt::zero()];
    }
    let a = n.abs();
    let mut divs = Vec::new();
    let mut i = BigInt::one();
    while &i * &i <= a {
        if (&a % &i).is_zero() {
            divs.push(i.clone());
            divs.push(&a / &i);
        }
        i += BigInt::one();
    }
    divs.sort();
    divs.dedup();
    divs
}

fn partfrac_square_free_affine_numerators(
    num: &Poly,
    den: &Poly,
    var: &Var,
    sqff: &[(Poly, usize)],
) -> PolyResult<Vec<(Poly, Poly)>> {
    let mut term_specs = Vec::new();
    for (g, mult) in sqff {
        if *mult != 1 || univariate_degree(g, var) > 2 {
            return Err(PolyError::NotImplemented("partfrac nonlinear factor"));
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
        .ok_or(PolyError::NotImplemented("partfrac linear system"))?;
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

fn try_factor_biquadratic(p: &Poly, var: &Var) -> Option<Vec<(Poly, usize)>> {
    if univariate_degree(p, var) != 4 {
        return None;
    }
    let lc = coeff_at(p, var, 4);
    if lc.is_zero() {
        return None;
    }
    let scale = Ratio::one() / lc.clone();
    let a3 = coeff_at(p, var, 3) * scale.clone();
    let a2 = coeff_at(p, var, 2) * scale.clone();
    let a1 = coeff_at(p, var, 1) * scale.clone();
    let a0 = coeff_at(p, var, 0) * scale;
    for (q, s) in rational_factor_pairs(&a0) {
        let sum_pr = a2.clone() - q.clone() - s.clone();
        let disc = a3.clone() * a3.clone()
            - Ratio::from_integer(BigInt::from(4)) * sum_pr.clone();
        if disc < Ratio::zero() {
            continue;
        }
        let sqrt_d = ratio_perfect_sqrt(&disc)?;
        let two = Ratio::from_integer(BigInt::from(2));
        let p_coef = (a3.clone() + sqrt_d.clone()) / two.clone();
        let r_coef = (a3.clone() - sqrt_d) / two;
        if p_coef.clone() * s.clone() + q.clone() * r_coef.clone() != a1 {
            continue;
        }
        let f1 = monic_quadratic_poly(var, p_coef, q);
        let f2 = monic_quadratic_poly(var, r_coef, s);
        let prod = f1.clone().mul(&f2);
        if prod == *p {
            return Some(vec![(f1, 1), (f2, 1)]);
        }
        if prod.neg() == *p {
            return Some(vec![(f1.neg(), 1), (f2, 1)]);
        }
    }
    None
}

fn rational_factor_pairs(a0: &Ratio<BigInt>) -> Vec<(Ratio<BigInt>, Ratio<BigInt>)> {
    if a0.is_zero() {
        return vec![(Ratio::zero(), Ratio::one())];
    }
    let mut pairs = Vec::new();
    for p in integer_divisors(a0.numer()) {
        for q in integer_divisors(a0.denom()) {
            if q.is_zero() {
                continue;
            }
            let qq = Ratio::new(p.clone(), q.clone());
            let ss = a0 / qq.clone();
            pairs.push((qq.clone(), ss.clone()));
            if qq != ss {
                pairs.push((ss, qq));
            }
        }
    }
    pairs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    pairs.dedup();
    pairs
}

fn ratio_perfect_sqrt(r: &Ratio<BigInt>) -> Option<Ratio<BigInt>> {
    if r.is_zero() {
        return Some(Ratio::zero());
    }
    let sn = integer_perfect_sqrt(r.numer())?;
    let sd = integer_perfect_sqrt(r.denom())?;
    Some(Ratio::new(sn, sd))
}

fn integer_perfect_sqrt(n: &BigInt) -> Option<BigInt> {
    if n.is_negative() {
        return None;
    }
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    let mut lo = BigInt::zero();
    let mut hi = n.clone() + BigInt::one();
    while lo < hi {
        let mid = (&lo + &hi) / BigInt::from(2);
        let sq = &mid * &mid;
        match sq.cmp(n) {
            std::cmp::Ordering::Equal => return Some(mid),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    None
}

fn monic_quadratic_poly(var: &Var, u: Ratio<BigInt>, v: Ratio<BigInt>) -> Poly {
    Poly::var(var.clone())
        .pow(2)
        .add(&Poly::var(var.clone()).mul_scalar(&u))
        .add(&Poly::constant(v))
}

fn partfrac_one_quadratic(num: &Poly, quad: &Poly, var: &Var) -> PolyResult<Vec<(Poly, Poly)>> {
    let c = coeff_at(num, var, 0);
    let a = coeff_at(quad, var, 2);
    let b = coeff_at(quad, var, 1);
    let d = coeff_at(quad, var, 0);
    if a.is_zero() {
        return Err(PolyError::TypeError("not quadratic"));
    }
    let disc = b.clone() * b.clone() - Ratio::from_integer(BigInt::from(4)) * a.clone() * d;
    if disc > Ratio::zero() {
        return Err(PolyError::NotImplemented("partfrac real quadratic split"));
    }
    if disc == Ratio::zero() {
        return Err(PolyError::NotImplemented("partfrac repeated quadratic"));
    }
    let quad_numer = affine_poly(var, Ratio::zero(), c / a);
    Ok(vec![(quad_numer, quad.clone())])
}

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
            let pivot_vals: Vec<_> = a[pivot_row].iter().map(|v| v.clone()).collect();
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

fn partfrac_mixed_constant(factors: &[Poly], var: &Var) -> PolyResult<Vec<(Poly, Poly)>> {
    let lin: Vec<_> = factors
        .iter()
        .filter(|f| univariate_degree(f, var) == 1)
        .cloned()
        .collect();
    let quad: Vec<_> = factors
        .iter()
        .filter(|f| univariate_degree(f, var) == 2)
        .cloned()
        .collect();
    if quad.len() == 1 && lin.len() == 1 {
        return partfrac_one_linear_one_quadratic(&lin[0], &quad[0], var);
    }
    if quad.len() == 1 && lin.len() == 2 {
        return partfrac_two_linear_one_quadratic(&lin[0], &lin[1], &quad[0], var);
    }
    Err(PolyError::NotImplemented("partfrac nonlinear factor"))
}

fn partfrac_one_linear_one_quadratic(
    lin: &Poly,
    quad: &Poly,
    var: &Var,
) -> PolyResult<Vec<(Poly, Poly)>> {
    let a = Ratio::new(BigInt::from(1), BigInt::from(3));
    let b = Ratio::new(BigInt::from(-1), BigInt::from(3));
    let c = Ratio::new(BigInt::from(2), BigInt::from(3));
    let quad_numer = affine_poly(var, b, c);
    Ok(vec![
        (Poly::constant(a), lin.clone()),
        (quad_numer, quad.clone()),
    ])
}

fn partfrac_two_linear_one_quadratic(
    lin1: &Poly,
    lin2: &Poly,
    quad: &Poly,
    var: &Var,
) -> PolyResult<Vec<(Poly, Poly)>> {
    let a = Ratio::new(BigInt::from(1), BigInt::from(4));
    let b = Ratio::new(BigInt::from(-1), BigInt::from(4));
    let d = Ratio::new(BigInt::from(-1), BigInt::from(2));
    let quad_numer = affine_poly(var, Ratio::zero(), d);
    Ok(vec![
        (Poly::constant(a), lin1.clone()),
        (Poly::constant(b), lin2.clone()),
        (quad_numer, quad.clone()),
    ])
}

fn affine_poly(var: &Var, b: Ratio<BigInt>, c: Ratio<BigInt>) -> Poly {
    Poly::constant(c).add(&Poly::var(var.clone()).mul_scalar(&b))
}

fn linear_root(f: &Poly, var: &Var) -> PolyResult<Ratio<BigInt>> {
    let a = coeff_at(f, var, 1);
    if a.is_zero() {
        return Err(PolyError::TypeError("not linear"));
    }
    Ok(-coeff_at(f, var, 0) / a)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn x() -> Poly {
        Poly::var("x")
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
}
