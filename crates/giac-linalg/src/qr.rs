use super::{from_dmatrix, to_dmatrix};

/// Gram–Schmidt QR decomposition via nalgebra.
/// Returns `(Q, R)` where Q has orthonormal columns and R is upper triangular.
pub fn qr_decomp(a: &[Vec<f64>]) -> Option<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    let mat = to_dmatrix(a)?;
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    if nrows == 0 || ncols == 0 {
        return None;
    }

    let qr = mat.qr();

    // Extract Q (nrows × ncols) — qr.q() gives the orthogonal factor
    let q = qr.q();
    // Extract R (ncols × ncols) — upper triangular
    let r = qr.r();

    Some((from_dmatrix(&q), from_dmatrix(&r)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qr_2x2_reconstructs() {
        let a = vec![vec![3.0, 5.0], vec![4.0, 5.0]];
        let (q, r) = qr_decomp(&a).unwrap();
        // Verify A ≈ Q × R
        for i in 0..2 {
            for j in 0..2 {
                let mut sum = 0.0;
                for k in 0..2 {
                    sum += q[i][k] * r[k][j];
                }
                assert!((sum - a[i][j]).abs() < 1e-9, "A[{i}][{j}] mismatch");
            }
        }
    }

    #[test]
    fn qr_3x2_tall_matrix() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
        let (q, r) = qr_decomp(&a).unwrap();
        assert_eq!(q.len(), 3); // 3 rows
        assert_eq!(q[0].len(), 2); // 2 cols (thin Q)
        assert_eq!(r.len(), 2); // 2 rows
        assert_eq!(r[0].len(), 2); // 2 cols
        // Verify reconstruction
        for i in 0..3 {
            for j in 0..2 {
                let mut sum = 0.0;
                for k in 0..2 {
                    sum += q[i][k] * r[k][j];
                }
                assert!((sum - a[i][j]).abs() < 1e-9, "A[{i}][{j}] mismatch");
            }
        }
    }
}
