#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod eigen;
mod lu;
mod qr;
mod svd;

pub use eigen::real_eigenvalues;
pub use lu::lu_decomp;
pub use qr::qr_decomp;
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
    // nalgebra stores data column-major; we provide row-major data
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
    if v.is_nan() {
        return "undef".to_string();
    }
    if v.is_infinite() {
        return if v.is_sign_positive() {
            "inf".to_string()
        } else {
            "-inf".to_string()
        };
    }
    if v == 0.0 {
        return "0".to_string();
    }
    let exp = v.abs().log10().floor() as i32;
    let scale = 10_f64.powi(digits as i32 - 1 - exp);
    let rounded = (v * scale).round() / scale;
    let s = format!("{rounded:.12}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}
