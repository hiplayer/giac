#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod f64_eigen;
mod gramschmidt;
mod lu;
mod numeric;
mod plugin;
mod qr;
mod symbolic;
mod symbolic_eigen;
mod svd;

pub use f64_eigen::real_eigenvalues;
pub use gramschmidt::eval_gramschmidt;
pub use lu::lu_decomp;
pub use numeric::{eval_lu, eval_qr, eval_svd};
pub use plugin::{install_linalg, xcas_default, DefaultLinalgPlugin};
pub use qr::qr_decomp;
pub use symbolic::{
    as_matrix, eval_charpoly, eval_det, eval_idn, eval_image, eval_inv, eval_ker,
    eval_linsolve, eval_matrix_mul, eval_matrix_pow, eval_pcar, eval_rref, eval_trace,
    eval_tran, f64_to_expr_numeric, try_to_f64_matrix,
};
pub use symbolic_eigen::{eval_egv, eval_jordan};
pub use svd::svd_decomp;

use nalgebra::DMatrix;

/// Convert a row-major Vec<Vec<f64>> to nalgebra DMatrix.
pub fn to_dmatrix(m: &[Vec<f64>]) -> Option<DMatrix<f64>> {
    if m.is_empty() {
        return None;
    }
    let nrows = m.len();
    let ncols = m[0].len();
    if !m.iter().all(|r| r.len() == ncols) {
        return None;
    }
    let mut entries = Vec::with_capacity(nrows * ncols);
    for row in m {
        entries.extend_from_slice(row);
    }
    Some(DMatrix::from_row_slice(nrows, ncols, &entries))
}

/// Convert a nalgebra DMatrix back to row-major Vec<Vec<f64>>.
pub fn from_dmatrix(m: &DMatrix<f64>) -> Vec<Vec<f64>> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    (0..nrows)
        .map(|i| (0..ncols).map(|j| m[(i, j)]).collect())
        .collect()
}

/// Format a float with `digits` significant figures (giac `evalf` style).
pub fn format_float(v: f64, digits: u32) -> String {
    giac_core::float_format::format_float(v, digits)
}
