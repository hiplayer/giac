//! Re-exports symbolic simplification from `giac-core` (facade crate per migration plan).

#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

pub use giac_core::{assert_equiv, expand, factor, is_zero, normal, ratnormal, sub};

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, normal, expand, factor, assert_equiv, Context, Expr};

    #[test]
    fn facade_normal_matches_core() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2));
        let r = normal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x^2+2*x+1");
    }

    #[test]
    fn facade_expand_distributes() {
        let ctx = Context::default();
        let e = Expr::mul(vec![
            Expr::add(vec![Expr::sym("x"), Expr::int(1)]),
            Expr::sym("y"),
        ]);
        let r = expand(e.as_ref(), &ctx).unwrap();
        assert!(format_expr(r.as_ref()).contains("x*y"));
    }

    #[test]
    fn facade_factor_x4_minus_1() {
        let ctx = Context::default();
        let e = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(4)),
            Expr::int(-1),
        ]);
        let r = factor(e.as_ref(), &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("(x-1)") && s.contains("(x+1)"));
    }

    #[test]
    fn facade_assert_equiv_commutative() {
        let ctx = Context::default();
        let a = Expr::add(vec![Expr::sym("x"), Expr::int(1)]);
        let b = Expr::add(vec![Expr::int(1), Expr::sym("x")]);
        assert!(assert_equiv(a.as_ref(), b.as_ref(), &ctx).unwrap());
    }

    #[test]
    fn facade_normal_mod_power() {
        let ctx = Context::default();
        let e = Expr::pow(
            std::sync::Arc::new(Expr::Mod(
                Expr::add(vec![
                    Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
                    Expr::int(1),
                ]),
                Expr::int(13),
            )),
            Expr::int(3),
        );
        let r = normal(e.as_ref(), &ctx).unwrap();
        assert!(format_expr(r.as_ref()).contains("% 13"));
    }
}
