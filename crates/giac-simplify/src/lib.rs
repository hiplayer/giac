//! Re-exports symbolic simplification from `giac-core` (facade crate per migration plan).

#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

pub use giac_core::{expand, factor, normal, ratnormal};

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, normal, Context, Expr};

    #[test]
    fn facade_normal_matches_core() {
        let ctx = Context::default();
        let e = Expr::pow(Expr::add(vec![Expr::sym("x"), Expr::int(1)]), Expr::int(2));
        let r = normal(e.as_ref(), &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "x^2+2*x+1");
    }
}
