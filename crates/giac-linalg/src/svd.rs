use super::{from_dmatrix, to_dmatrix};

/// SVD decomposition via nalgebra: `A = U Σ Vᵀ`.
/// Returns `(U, sigma, Vᵀ)` where `sigma` is a vector of singular values (sorted descending).
pub fn svd_decomp(a: &[Vec<f64>]) -> Option<(Vec<Vec<f64>>, Vec<f64>, Vec<Vec<f64>>)> {
    let mat = to_dmatrix(a)?;
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    if nrows == 0 || ncols == 0 {
        return None;
    }

    let svd = mat.svd(true, true);

    let u = svd.u?;
    let vt = svd.v_t?;

    // Singular values: diagonal of Σ, sorted descending (nalgebra guarantees this)
    let k = nrows.min(ncols);
    let sigma: Vec<f64> = (0..k).map(|i| svd.singular_values[i]).collect();

    Some((from_dmatrix(&u), sigma, from_dmatrix(&vt)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn svd_2x2() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let (_u, s, _vt) = svd_decomp(&a).unwrap();
        assert!(s[0] > s[1], "singular values should be descending");
        assert!(s[0] > 5.0, "largest SV should be > 5");
    }

    #[test]
    fn svd_3x2_rectangular() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
        let (u, s, vt) = svd_decomp(&a).unwrap();
        assert_eq!(u.len(), 3);
        assert_eq!(s.len(), 2);
        assert_eq!(vt.len(), 2);
        // Verify A ≈ U Σ Vᵀ
        for i in 0..3 {
            for j in 0..2 {
                let mut sum = 0.0;
                for k in 0..2 {
                    sum += u[i][k] * s[k] * vt[k][j];
                }
                assert!((sum - a[i][j]).abs() < 1e-9, "A[{i}][{j}] mismatch");
            }
        }
    }

    #[test]
    fn svd_3x3_reconstructs() {
        let a = vec![
            vec![1.0, 2.0, 1.0],
            vec![3.0, 4.0, 1.0],
            vec![1.0, 5.0, 6.0],
        ];
        let (u, s, vt) = svd_decomp(&a).unwrap();
        assert_eq!(s.len(), 3);
        for i in 0..3 {
            for j in 0..3 {
                let mut sum = 0.0;
                for k in 0..3 {
                    sum += u[i][k] * s[k] * vt[k][j];
                }
                assert!(
                    (sum - a[i][j]).abs() < 1e-9,
                    "A[{i}][{j}] got {sum} want {}",
                    a[i][j]
                );
            }
        }
    }
}
