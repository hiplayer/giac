use std::sync::Arc;

use crate::{eval, format_expr, Context, Expr, FuncKind};

#[test]
fn eval_quo_direct() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Quo,
        vec![
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(3)), Expr::int(-1)]),
            Expr::add(vec![Expr::sym("x"), Expr::int(-1)]),
        ],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    assert_eq!(format_expr(r.as_ref()), "x^2+x+1");
}

#[test]
fn eval_rem_direct() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Rem,
        vec![
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(3)), Expr::int(-1)]),
            Expr::add(vec![Expr::sym("x"), Expr::int(-1)]),
        ],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    assert_eq!(format_expr(r.as_ref()), "0");
}

#[test]
fn eval_content_direct() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Content,
        vec![Expr::add(vec![
            Expr::mul(vec![Expr::int(6), Expr::pow(Expr::sym("x"), Expr::int(2))]),
            Expr::int(12),
        ])],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    assert_eq!(format_expr(r.as_ref()), "6");
}

#[test]
fn eval_smod_irem_direct() {
    let ctx = Context::xcas_default();
    assert_eq!(
        format_expr(
            eval(
                Expr::func(FuncKind::Smod, vec![Expr::int(17), Expr::int(5)]).as_ref(),
                &ctx
            )
            .unwrap()
            .as_ref()
        ),
        "2"
    );
    assert_eq!(
        format_expr(
            eval(
                Expr::func(FuncKind::Irem, vec![Expr::int(17), Expr::int(5)]).as_ref(),
                &ctx
            )
            .unwrap()
            .as_ref()
        ),
        "2"
    );
}

#[test]
fn eval_factor_direct() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Factor,
        vec![Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(4)),
            Expr::int(-1),
        ])],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    assert_eq!(format_expr(r.as_ref()), "(x-1)*(x+1)*(x^2+1)");
}

#[test]
fn eval_mod_gcd_direct() {
    let ctx = Context::xcas_default();
    let a = Arc::new(Expr::Mod(
        Expr::add(vec![
            Expr::mul(vec![Expr::int(2), Expr::pow(Expr::sym("x"), Expr::int(2))]),
            Expr::int(5),
        ]),
        Expr::int(13),
    ));
    let b = Arc::new(Expr::Mod(
        Expr::add(vec![
            Expr::mul(vec![Expr::int(5), Expr::pow(Expr::sym("x"), Expr::int(2))]),
            Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
            Expr::int(-3),
        ]),
        Expr::int(13),
    ));
    let e = Expr::func(FuncKind::Gcd, vec![a, b]);
    let r = eval(e.as_ref(), &ctx).unwrap();
    assert!(format_expr(r.as_ref()).contains('x'));
}

#[test]
fn eval_horner_resultant() {
    let ctx = Context::xcas_default();
    let h = Expr::func(
        FuncKind::Horner,
        vec![
            Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(4)),
                Expr::mul(vec![Expr::int(2), Expr::pow(Expr::sym("x"), Expr::int(3))]),
                Expr::mul(vec![Expr::int(-3), Expr::pow(Expr::sym("x"), Expr::int(2))]),
                Expr::sym("x"),
                Expr::int(-2),
            ]),
            Expr::int(1),
        ],
    );
    assert_eq!(format_expr(eval(h.as_ref(), &ctx).unwrap().as_ref()), "-1");

    let res = Expr::func(
        FuncKind::Resultant,
        vec![
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(-1)]),
            Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(3)), Expr::int(-1)]),
            Expr::sym("x"),
        ],
    );
    assert_eq!(format_expr(eval(res.as_ref(), &ctx).unwrap().as_ref()), "0");
}

#[test]
fn eval_partfrac_direct() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Partfrac,
        vec![
            Expr::pow(
                Expr::add(vec![
                    Expr::pow(Expr::sym("x"), Expr::int(2)),
                    Expr::int(-1),
                ]),
                Expr::int(-1),
            ),
            Expr::sym("x"),
        ],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    let s = format_expr(r.as_ref());
    assert!(s.contains("x-1") && s.contains("x+1"));
}

#[test]
fn eval_chinrem_direct() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Chinrem,
        vec![
            Arc::new(Expr::List(vec![
                Expr::add(vec![Expr::sym("x"), Expr::int(2)]),
                Expr::add(vec![Expr::pow(Expr::sym("x"), Expr::int(2)), Expr::int(1)]),
            ])),
            Arc::new(Expr::List(vec![
                Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
                Expr::add(vec![
                    Expr::pow(Expr::sym("x"), Expr::int(2)),
                    Expr::sym("x"),
                    Expr::int(1),
                ]),
            ])),
        ],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    assert!(format_expr(r.as_ref()).starts_with('['));
}

#[test]
fn eval_greduce_direct() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Greduce,
        vec![
            Expr::add(vec![
                Expr::mul(vec![Expr::sym("x"), Expr::sym("y")]),
                Expr::int(-1),
            ]),
            Arc::new(Expr::List(vec![
                Expr::add(vec![
                    Expr::pow(Expr::sym("x"), Expr::int(2)),
                    Expr::mul(vec![Expr::int(-1), Expr::pow(Expr::sym("y"), Expr::int(2))]),
                ]),
                Expr::add(vec![
                    Expr::mul(vec![Expr::int(2), Expr::sym("x"), Expr::sym("y")]),
                    Expr::mul(vec![Expr::int(-1), Expr::pow(Expr::sym("y"), Expr::int(2))]),
                ]),
                Expr::pow(Expr::sym("y"), Expr::int(3)),
            ])),
            Arc::new(Expr::List(vec![Expr::sym("x"), Expr::sym("y"), Expr::sym("z")])),
        ],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    let s = format_expr(r.as_ref());
    assert_eq!(s, "1/2*y^2-1");
}

#[test]
fn eval_greduce_circle() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Greduce,
        vec![
            Expr::add(vec![
                Expr::pow(Expr::sym("x"), Expr::int(2)),
                Expr::pow(Expr::sym("y"), Expr::int(2)),
                Expr::int(-1),
            ]),
            Arc::new(Expr::List(vec![
                Expr::add(vec![
                    Expr::pow(Expr::sym("x"), Expr::int(2)),
                    Expr::mul(vec![Expr::int(-1), Expr::pow(Expr::sym("y"), Expr::int(2))]),
                ]),
                Expr::add(vec![
                    Expr::mul(vec![Expr::int(2), Expr::sym("x"), Expr::sym("y")]),
                    Expr::int(-1),
                ]),
            ])),
            Arc::new(Expr::List(vec![Expr::sym("x"), Expr::sym("y")])),
        ],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    assert_eq!(format_expr(r.as_ref()), "2*y^2-1");
}

#[test]
fn eval_normal_mod_power() {
    let ctx = Context::xcas_default();
    let e = Expr::func(
        FuncKind::Normal,
        vec![Expr::pow(
            Arc::new(Expr::Mod(
                Expr::add(vec![
                    Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
                    Expr::int(1),
                ]),
                Expr::int(13),
            )),
            Expr::int(5),
        )],
    );
    let r = eval(e.as_ref(), &ctx).unwrap();
    let s = format_expr(r.as_ref());
    assert!(s.contains("(6 % 13)*x^5") || s.contains("x^5"), "got {s}");
    verify_sympy_style_mod(&s);
}

fn verify_sympy_style_mod(s: &str) {
    assert!(s.contains("% 13"));
}
