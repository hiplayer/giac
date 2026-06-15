use super::to_dmatrix;

/// Real eigenvalues of a square matrix (when all eigenvalues are real).
pub fn real_eigenvalues(a: &[Vec<f64>]) -> Option<Vec<f64>> {
    let m = to_dmatrix(a)?;
    let n = m.nrows();
    if n == 0 {
        return Some(vec![]);
    }
    let schur = m.schur();
    if let Some(evals) = schur.eigenvalues() {
        return Some(evals.as_slice().to_vec());
    }
    let complex = schur.complex_eigenvalues();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let z = complex[i];
        if z.im.abs() > 1e-8 {
            return None;
        }
        out.push(z.re);
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eigenvalues_empty() {
        assert_eq!(real_eigenvalues(&[]), None);
    }

    #[test]
    fn eigenvalues_complex_returns_none() {
        // rotation matrix: eigenvalues are ±i
        let a = vec![vec![0.0, 1.0], vec![-1.0, 0.0]];
        assert!(real_eigenvalues(&a).is_none());
    }

    #[test]
    fn eigenvalues_diagonal_3x3() {
        let a = vec![vec![1.0, 0.0, 0.0], vec![0.0, 2.0, 0.0], vec![0.0, 0.0, 3.0]];
        let ev = real_eigenvalues(&a).unwrap();
        let mut ev = ev;
        ev.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((ev[0] - 1.0).abs() < 1e-9);
        assert!((ev[1] - 2.0).abs() < 1e-9);
        assert!((ev[2] - 3.0).abs() < 1e-9);
    }
}
