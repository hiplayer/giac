//! giac-linalg plugin: wires symbolic/numeric matrix ops into `giac-core::Context`.

use std::sync::Arc;

use giac_core::{Context, EvalError, ExprArc, Ident, LinalgPlugin};

use crate::gramschmidt::eval_gramschmidt;
use crate::numeric::{eval_lu, eval_qr, eval_svd};
use crate::symbolic::{
    eval_charpoly, eval_det, eval_idn, eval_image, eval_inv, eval_ker, eval_linsolve,
    eval_matrix_mul, eval_matrix_pow, eval_pcar, eval_rref, eval_trace, eval_tran,
};
use crate::symbolic_eigen::{eval_egv, eval_jordan};

/// Default implementation of [`LinalgPlugin`].
pub struct DefaultLinalgPlugin;

impl LinalgPlugin for DefaultLinalgPlugin {
    fn eval_matrix_mul(&self, a: &ExprArc, b: &ExprArc) -> Result<ExprArc, EvalError> {
        eval_matrix_mul(a, b)
    }

    fn eval_matrix_pow(&self, m: &ExprArc, exp: u32, ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_matrix_pow(m, exp, ctx)
    }

    fn eval_idn(&self, n: usize) -> ExprArc {
        eval_idn(n)
    }

    fn eval_rref(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_rref(args, ctx)
    }

    fn eval_inv(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_inv(m, ctx)
    }

    fn eval_det(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_det(m, ctx)
    }

    fn eval_tran(&self, m: &ExprArc) -> Result<ExprArc, EvalError> {
        eval_tran(m)
    }

    fn eval_ker(&self, m: &ExprArc) -> Result<ExprArc, EvalError> {
        eval_ker(m)
    }

    fn eval_image(&self, m: &ExprArc) -> Result<ExprArc, EvalError> {
        eval_image(m)
    }

    fn eval_pcar(&self, m: &ExprArc) -> Result<ExprArc, EvalError> {
        eval_pcar(m)
    }

    fn eval_charpoly(
        &self,
        m: &ExprArc,
        var: &Ident,
        ctx: &Context,
    ) -> Result<ExprArc, EvalError> {
        eval_charpoly(m, var, ctx)
    }

    fn eval_linsolve(
        &self,
        eqs: &ExprArc,
        vars: &ExprArc,
        ctx: &Context,
    ) -> Result<ExprArc, EvalError> {
        eval_linsolve(eqs, vars, ctx)
    }

    fn eval_jordan(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_jordan(m, ctx)
    }

    fn eval_egv(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_egv(m, ctx)
    }

    fn eval_lu(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_lu(m, ctx)
    }

    fn eval_qr(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_qr(m, ctx)
    }

    fn eval_svd(&self, m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_svd(m, ctx)
    }

    fn eval_gramschmidt(&self, args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
        eval_gramschmidt(args, ctx)
    }

    fn eval_trace(&self, m: &ExprArc) -> Result<ExprArc, EvalError> {
        eval_trace(m)
    }
}

/// Install the default linear algebra plugin on `ctx`.
pub fn install_linalg(ctx: &mut Context) {
    ctx.set_linalg_plugin(Arc::new(DefaultLinalgPlugin));
}

/// `Context::xcas_default()` with linear algebra enabled.
pub fn xcas_default() -> Context {
    let mut ctx = Context::xcas_default();
    install_linalg(&mut ctx);
    ctx
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{eval, format_expr, Expr, ExprArc, FuncKind};

    use super::xcas_default;

    fn mat2() -> ExprArc {
        Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2)],
            vec![Expr::int(3), Expr::int(4)],
        ]))
    }

    #[test]
    fn eval_matrix_builtins_via_plugin() {
        let ctx = xcas_default();
        let m = mat2();

        let det = eval(Expr::func(FuncKind::Det, vec![Arc::clone(&m)]).as_ref(), &ctx).unwrap();
        let det_s = format_expr(det.as_ref());
        assert!(det_s == "-2" || det_s == "1*4-1*2*3");

        let inv = eval(Expr::func(FuncKind::Inv, vec![Arc::clone(&m)]).as_ref(), &ctx).unwrap();
        let inv_s = format_expr(inv.as_ref());
        assert!(
            inv_s == "[[-2,1],[3/2,-1/2]]" || inv_s.contains("4") && inv_s.contains("2*3"),
            "unexpected inv: {inv_s}"
        );

        let tran = eval(Expr::func(FuncKind::Tran, vec![Arc::clone(&m)]).as_ref(), &ctx).unwrap();
        assert!(format_expr(tran.as_ref()).starts_with("matrix[[1,3]"));

        let mul = eval(
            Expr::mul(vec![Arc::clone(&m), Arc::clone(&m)]).as_ref(),
            &ctx,
        )
        .unwrap();
        assert_eq!(format_expr(mul.as_ref()), "[[7,10],[15,22]]");

        let idn = eval(Expr::func(FuncKind::Idn, vec![Expr::int(2)]).as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(idn.as_ref()), "matrix[[1,0],[0,1]]");
    }
}
