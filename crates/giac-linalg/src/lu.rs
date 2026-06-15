use nalgebra::DMatrix;

use super::{from_dmatrix, to_dmatrix};

/// Doolittle LU with partial pivoting via nalgebra.
/// Returns `(perm, L, U)` where `perm` is 0-based row permutation indices.
///
/// The permutation satisfies: P*A = L*U, where P is constructed from `perm`
/// such that row `i` of P*A comes from row `perm[i]` of A.
pub fn lu_decomp(a: &[Vec<f64>]) -> Option<(Vec<usize>, Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    let mat = to_dmatrix(a)?;
    let n = mat.nrows();
    if n != mat.ncols() || n == 0 {
        return None;
    }

    let lu = mat.lu();

    // Build permutation vector by applying permutation to identity rows
    let p = lu.p();

    // Apply permutation to identity matrix to discover mapping
    let mut identity: DMatrix<f64> = DMatrix::identity(n, n);
    p.permute_rows(&mut identity);
    // identity[(i, j)] == 1 means original row j maps to permuted row i
    let mut perm_vec = vec![0usize; n];
    for i in 0..n {
        for j in 0..n {
            if identity[(i, j)] == 1.0 {
                perm_vec[i] = j;
                break;
            }
        }
    }

    let l = lu.l();
    let u = lu.u();

    Some((perm_vec, from_dmatrix(&l), from_dmatrix(&u)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lu_empty_fails() {
        assert!(lu_decomp(&[]).is_none());
    }

    #[test]
    fn lu_non_square_fails() {
        assert!(lu_decomp(&[vec![1.0, 2.0], vec![3.0, 4.0, 5.0]]).is_none());
    }

    #[test]
    fn lu_2x2() {
        let a = vec![vec![3.0, 5.0], vec![4.0, 5.0]];
        let (p, l, u) = lu_decomp(&a).unwrap();
        assert_eq!(p.len(), 2);
        // Verify PA = LU
        for i in 0..2 {
            for j in 0..2 {
                let pa_ij = a[p[i]][j];
                let mut lu_ij = 0.0;
                for k in 0..2 {
                    lu_ij += l[i][k] * u[k][j];
                }
                assert!((lu_ij - pa_ij).abs() < 1e-9, "PA[{i}][{j}]={pa_ij} != LU={lu_ij}");
            }
        }
    }

    #[test]
    fn lu_3x3_identity() {
        let a = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 1.0]];
        let (_, l, u) = lu_decomp(&a).unwrap();
        assert!((l[0][0] - 1.0).abs() < 1e-12);
        assert!((u[0][0] - 1.0).abs() < 1e-12);
    }
}
