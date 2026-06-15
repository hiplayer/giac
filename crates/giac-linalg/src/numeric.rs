//! Numeric matrix decompositions (bridges to `giac-linalg`).

use std::sync::Arc;

use crate::{lu_decomp, qr_decomp, svd_decomp};

use giac_core::Context;
use giac_core::EvalError;
use giac_core::{Expr, ExprArc};
use crate::symbolic::{f64_to_expr_numeric, try_to_f64_matrix};

pub fn eval_lu(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let a = try_to_f64_matrix(m, ctx)?;
    let (perm, l, u) = lu_decomp(&a).ok_or(EvalError::TypeError("LU failed"))?;
    let p_expr: Vec<ExprArc> = perm.into_iter().map(|i| Expr::int(i as i64)).collect();
    Ok(Arc::new(Expr::Seq(vec![
        Arc::new(Expr::List(p_expr)),
        f64_matrix_to_giac(&l),
        f64_matrix_to_giac(&u),
    ])))
}

pub fn eval_qr(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let a = try_to_f64_matrix(m, ctx)?;
    let (q, r) = qr_decomp(&a).ok_or(EvalError::TypeError("QR failed"))?;
    Ok(Arc::new(Expr::Seq(vec![
        f64_matrix_to_giac(&q),
        f64_matrix_to_giac(&r),
    ])))
}

pub fn eval_svd(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let a = try_to_f64_matrix(m, ctx)?;
    let (u, sigma, vt) = svd_decomp(&a).ok_or(EvalError::TypeError("SVD failed"))?;
    let u_e = f64_matrix_to_giac(&u);
    let s_e: Vec<ExprArc> = sigma.into_iter().map(f64_to_expr_numeric).collect();
    let vt_e = f64_matrix_to_giac(&vt);
    Ok(Arc::new(Expr::Seq(vec![
        u_e,
        Arc::new(Expr::List(s_e)),
        vt_e,
    ])))
}

fn f64_matrix_to_giac(m: &[Vec<f64>]) -> ExprArc {
    let rows: Vec<Vec<ExprArc>> = m
        .iter()
        .map(|row| row.iter().map(|&v| f64_to_expr_numeric(v)).collect())
        .collect();
    Arc::new(Expr::GiacMatrix(rows))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin::xcas_default;

    fn mat2() -> ExprArc {
        Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2)],
            vec![Expr::int(3), Expr::int(4)],
        ]))
    }

    #[test]
    fn eval_lu_qr_svd_2x2() {
        let ctx = xcas_default();
        let m = mat2();
        assert!(eval_lu(&m, &ctx).is_ok());
        assert!(eval_qr(&m, &ctx).is_ok());
        let svd = eval_svd(&m, &ctx).unwrap();
        assert!(matches!(svd.as_ref(), Expr::Seq(items) if items.len() == 3));
    }
}
