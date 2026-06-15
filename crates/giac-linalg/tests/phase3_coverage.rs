//! Phase 3 linear-algebra coverage: exercise giac-linalg APIs not hit by conformance alone.

use std::sync::Arc;

use giac_core::{eval, format_expr, Expr, ExprArc, FuncKind, Ident, LinalgPlugin};
use giac_linalg::{
    as_matrix, eval_charpoly, eval_det, eval_egv, eval_idn, eval_image, eval_inv, eval_jordan,
    eval_ker, eval_linsolve, eval_lu, eval_matrix_mul, eval_matrix_pow, eval_pcar, eval_qr,
    eval_rref, eval_svd, eval_trace, eval_tran, f64_to_expr_numeric, format_float,
    is_identity_matrix, real_eigenvalues, to_dmatrix, DefaultLinalgPlugin,
};

fn mat2() -> ExprArc {
    Arc::new(Expr::Matrix(vec![
        vec![Expr::int(1), Expr::int(2)],
        vec![Expr::int(3), Expr::int(4)],
    ]))
}

fn ctx() -> giac_core::Context {
    giac_calculus::xcas_default()
}

#[test]
fn numeric_lu_qr_svd_via_expr() {
    let ctx = ctx();
    let m = mat2();
    for f in [
        eval_lu(&m, &ctx),
        eval_qr(&m, &ctx),
        eval_svd(&m, &ctx),
    ] {
        assert!(f.is_ok(), "{f:?}");
    }
}

#[test]
fn symbolic_rref_inv_ker_image() {
    let ctx = ctx();
    let m = mat2();
    let rref = eval_rref(&[Arc::clone(&m)], &ctx).unwrap();
    assert!(format_expr(rref.as_ref()).contains('1'));

    let inv = eval_inv(&m, &ctx).unwrap();
    assert!(!format_expr(inv.as_ref()).is_empty());

    let singular = Arc::new(Expr::Matrix(vec![
        vec![Expr::int(1), Expr::int(2)],
        vec![Expr::int(3), Expr::int(6)],
    ]));
    let ker = eval_ker(&singular).unwrap();
    assert!(!format_expr(ker.as_ref()).is_empty());
    let image = eval_image(&singular).unwrap();
    assert!(!format_expr(image.as_ref()).is_empty());
}

#[test]
fn symbolic_linsolve_charpoly_pcar_trace() {
    let ctx = ctx();
    let m = mat2();
    let x = Ident::new("x");
    let _y = Ident::new("y");

    let eqs = Arc::new(Expr::List(vec![
        Arc::new(Expr::Relation(
            giac_core::RelOp::Eq,
            Expr::add(vec![
                Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
                Expr::sym("y"),
            ]),
            Expr::int(3),
        )),
        Arc::new(Expr::Relation(
            giac_core::RelOp::Eq,
            Expr::add(vec![Expr::sym("x"), Expr::mul(vec![Expr::int(-1), Expr::sym("y")])]),
            Expr::int(1),
        )),
    ]));
    let vars = Arc::new(Expr::List(vec![Expr::sym("x"), Expr::sym("y")]));
    let sol = eval_linsolve(&eqs, &vars, &ctx).unwrap();
    assert!(format_expr(sol.as_ref()).contains('x') || format_expr(sol.as_ref()).contains('1'));

    let cp = eval_charpoly(&m, &x, &ctx).unwrap();
    assert!(format_expr(cp.as_ref()).contains('x'));

    let pcar = eval_pcar(&m).unwrap();
    assert!(format_expr(pcar.as_ref()).starts_with("poly1["));

    let tr_s = format_expr(eval_trace(&m).unwrap().as_ref());
    assert!(tr_s == "5" || tr_s == "1+4", "trace got {tr_s}");
}

#[test]
fn matrix_mul_pow_tran_idn() {
    let ctx = ctx();
    let m = mat2();
    let prod = eval_matrix_mul(&m, &m).unwrap();
    let prod = eval(prod.as_ref(), &ctx).unwrap();
    assert_eq!(format_expr(prod.as_ref()), "[[7,10],[15,22]]");

    let sq = eval_matrix_pow(&m, 2, &ctx).unwrap();
    assert_eq!(format_expr(sq.as_ref()), "[[7,10],[15,22]]");

    let id = eval_idn(2);
    assert!(is_identity_matrix(&id));

    let t = eval_tran(&m).unwrap();
    assert!(format_expr(t.as_ref()).starts_with("matrix[[1,3]"));
}

#[test]
fn eigen_egv_jordan_2x2() {
    let ctx = ctx();
    let m = Arc::new(Expr::Matrix(vec![
        vec![Expr::int(1), Expr::int(2)],
        vec![Expr::int(3), Expr::int(4)],
    ]));
    let ev = eval_egv(&m, &ctx).unwrap();
    assert!(format_expr(ev.as_ref()).contains('(') || format_expr(ev.as_ref()).contains('/'));

    let jblock = Arc::new(Expr::Matrix(vec![
        vec![Expr::int(1), Expr::int(1)],
        vec![Expr::int(0), Expr::int(1)],
    ]));
    let j = eval_jordan(&jblock, &ctx).unwrap();
    assert!(format_expr(j.as_ref()).contains("[[1,1]"));
}

#[test]
fn linalg_plugin_dispatches_all() {
    let ctx = ctx();
    let plugin = DefaultLinalgPlugin;
    let m = mat2();
    let x = Ident::new("x");

    plugin.eval_matrix_mul(&m, &m).unwrap();
    plugin.eval_matrix_pow(&m, 2, &ctx).unwrap();
    plugin.eval_idn(2);
    plugin.eval_rref(&[Arc::clone(&m)], &ctx).unwrap();
    plugin.eval_inv(&m, &ctx).unwrap();
    plugin.eval_det(&m, &ctx).unwrap();
    plugin.eval_tran(&m).unwrap();
    plugin.eval_ker(&m).unwrap();
    plugin.eval_image(&m).unwrap();
    plugin.eval_pcar(&m).unwrap();
    plugin.eval_charpoly(&m, &x, &ctx).unwrap();
    plugin.eval_jordan(&m, &ctx).ok();
    plugin.eval_egv(&m, &ctx).unwrap();
    plugin.eval_lu(&m, &ctx).unwrap();
    plugin.eval_qr(&m, &ctx).unwrap();
    plugin.eval_svd(&m, &ctx).unwrap();
    plugin.eval_trace(&m).unwrap();

    let eqs = Arc::new(Expr::List(vec![Arc::new(Expr::Relation(
        giac_core::RelOp::Eq,
        Expr::add(vec![Expr::sym("x"), Expr::sym("y")]),
        Expr::int(3),
    ))]));
    let vars = Arc::new(Expr::List(vec![Expr::sym("x"), Expr::sym("y")]));
    plugin.eval_linsolve(&eqs, &vars, &ctx).unwrap();
}

#[test]
fn f64_helpers_and_eigenvalues() {
    assert_eq!(format_float(0.0, 10), "0");
    assert_eq!(format_float(f64::NAN, 10), "undef");
    assert!(to_dmatrix(&[]).is_none());

    let ev = real_eigenvalues(&[vec![1.0, 0.0], vec![0.0, 2.0]]).unwrap();
    assert_eq!(ev.len(), 2);

    assert!(real_eigenvalues(&[vec![0.0, 1.0], vec![-1.0, 0.0]]).is_none());

    let rat = f64_to_expr_numeric(0.5);
    assert!(format_expr(rat.as_ref()).contains('/') || format_expr(rat.as_ref()) == "1/2");
}

#[test]
fn as_matrix_and_det_numeric() {
    let ctx = ctx();
    let m = mat2();
    assert_eq!(as_matrix(&m).unwrap().len(), 2);
    assert_eq!(eval_det(&m, &ctx).unwrap(), Expr::int(-2));
}

#[test]
fn eval_builtin_matrix_ops_via_core() {
    let ctx = ctx();
    let m = mat2();
    let det = eval(Expr::func(FuncKind::Det, vec![Arc::clone(&m)]).as_ref(), &ctx).unwrap();
    assert_eq!(det, Expr::int(-2));
}

#[test]
fn matrix_pow_zero_and_mod_rref() {
    let ctx = ctx();
    let m = mat2();
    let id = eval_matrix_pow(&m, 0, &ctx).unwrap();
    assert!(is_identity_matrix(&id));

    let row1 = Arc::new(Expr::Mod(
        Arc::new(Expr::List(vec![Expr::int(1), Expr::int(2), Expr::int(3)])),
        Expr::int(5),
    ));
    let row2 = Arc::new(Expr::Mod(
        Arc::new(Expr::List(vec![Expr::int(4), Expr::int(5), Expr::int(6)])),
        Expr::int(5),
    ));
    let r = eval_rref(&[row1, row2], &ctx).unwrap();
    let s = format_expr(r.as_ref());
    assert!(s.contains('1'), "mod rref got {s}");
}

#[test]
fn inv_3x3_and_singular_errors() {
    let ctx = ctx();
    let m3 = Arc::new(Expr::Matrix(vec![
        vec![Expr::int(1), Expr::int(0), Expr::int(0)],
        vec![Expr::int(0), Expr::int(2), Expr::int(0)],
        vec![Expr::int(0), Expr::int(0), Expr::int(3)],
    ]));
    eval_inv(&m3, &ctx).unwrap();

    let singular3 = Arc::new(Expr::Matrix(vec![
        vec![Expr::int(1), Expr::int(2), Expr::int(3)],
        vec![Expr::int(2), Expr::int(4), Expr::int(6)],
        vec![Expr::int(1), Expr::int(1), Expr::int(1)],
    ]));
    assert!(eval_inv(&singular3, &ctx).is_err());

    let bad_mul = eval_matrix_mul(&mat2(), &Arc::new(Expr::Matrix(vec![vec![Expr::int(1)]])));
    assert!(bad_mul.is_err());
}

#[test]
fn charpoly_3x3_and_egv_jordan_limits() {
    let ctx = ctx();
    let m3 = Arc::new(Expr::Matrix(vec![
        vec![Expr::int(2), Expr::int(0), Expr::int(0)],
        vec![Expr::int(0), Expr::int(3), Expr::int(0)],
        vec![Expr::int(0), Expr::int(0), Expr::int(4)],
    ]));
    let x = Ident::new("x");
    let cp = eval_charpoly(&m3, &x, &ctx).unwrap();
    assert!(format_expr(cp.as_ref()).contains('x'));

    let ev = eval_egv(&m3, &ctx).unwrap();
    assert_eq!(format_expr(ev.as_ref()), "2,3,4");

    let j = eval_jordan(&m3, &ctx).unwrap();
    assert!(format_expr(j.as_ref()).contains("[[2,0,0]"));

    let big = Arc::new(Expr::GiacMatrix(vec![
        vec![Expr::int(1); 4],
        vec![Expr::int(1); 4],
        vec![Expr::int(1); 4],
        vec![Expr::int(1); 4],
    ]));
    assert!(eval_egv(&big, &ctx).is_err());
    assert!(eval_jordan(&big, &ctx).is_err());
}

#[test]
fn gramschmidt_via_eval_and_plugin() {
    let ctx = ctx();
    let vectors = Arc::new(Expr::List(vec![Expr::int(1), Expr::add(vec![Expr::int(1), Expr::sym("x")])]));
    let inner = Expr::func(
        FuncKind::Lambda,
        vec![
            Arc::new(Expr::List(vec![Expr::sym("p"), Expr::sym("q")])),
            Expr::func(
                FuncKind::Integrate,
                vec![
                    Expr::mul(vec![Expr::sym("p"), Expr::sym("q")]),
                    Expr::sym("x"),
                    Expr::int(-1),
                    Expr::int(1),
                ],
            ),
        ],
    );
    let r = eval(
        Expr::func(FuncKind::Gramschmidt, vec![vectors, inner]).as_ref(),
        &ctx,
    )
    .unwrap();
    let s = format_expr(r.as_ref());
    assert!(s.contains("sqrt"), "gramschmidt got {s}");

    assert!(DefaultLinalgPlugin.eval_gramschmidt(&[], &ctx).is_err());
}

#[test]
fn det_singular_and_trace_errors() {
    let ctx = ctx();
    let singular = Arc::new(Expr::Matrix(vec![
        vec![Expr::int(1), Expr::int(2)],
        vec![Expr::int(2), Expr::int(4)],
    ]));
    assert_eq!(eval_det(&singular, &ctx).unwrap(), Expr::int(0));

    let nonsquare = Arc::new(Expr::Matrix(vec![
        vec![Expr::int(1), Expr::int(2)],
        vec![Expr::int(3), Expr::int(4)],
        vec![Expr::int(5), Expr::int(6)],
    ]));
    assert!(eval_trace(&nonsquare).is_err());
    assert!(eval_tran(&Arc::new(Expr::Matrix(vec![]))).is_err());
}
