//! Semantic verification for giac-linalg unit tests (`.doc/conformance-testing.md` §3).
//! See `.doc/issues/GIAC-expr-api-test-contains-cleanup.md` H2.
#![allow(clippy::expect_used)] // test helper: expect carries assert context

use std::collections::HashMap;
use std::sync::Arc;

use giac_core::{eval, eval_subst_map, format_expr, Context, Expr, ExprArc, FuncKind, Ident, RelOp};
use giac_simplify::{assert_equiv, expand, is_zero};

use crate::symbolic::as_matrix;

pub fn matrix_rows(e: &ExprArc) -> Vec<Vec<ExprArc>> {
    as_matrix(e).expect("expected matrix")
}

pub fn list_items(e: &ExprArc) -> &[ExprArc] {
    match e.as_ref() {
        Expr::List(v) | Expr::Seq(v) => v,
        other => panic!("expected list, got {other:?}"),
    }
}

pub fn var_idents(vars: &ExprArc) -> Vec<Ident> {
    list_items(vars)
        .iter()
        .map(|e| match e.as_ref() {
            Expr::Symbol(id) => id.clone(),
            other => panic!("expected variable symbol, got {other:?}"),
        })
        .collect()
}

pub fn ctx_with_bindings(ctx: &Context, vars: &[Ident], values: &[ExprArc]) -> Context {
    assert_eq!(vars.len(), values.len(), "vars/values length mismatch");
    let mut c = ctx.clone();
    for (v, val) in vars.iter().zip(values) {
        c.set(v.clone(), Arc::clone(val));
    }
    c
}

pub fn assert_matrix_equiv(a: &ExprArc, b: &ExprArc, ctx: &Context) {
    let a_rows = matrix_rows(a);
    let b_rows = matrix_rows(b);
    assert_eq!(
        a_rows.len(),
        b_rows.len(),
        "row count mismatch: {} vs {}",
        format_expr(a.as_ref()),
        format_expr(b.as_ref())
    );
    for (i, (ar, br)) in a_rows.iter().zip(b_rows.iter()).enumerate() {
        assert_eq!(ar.len(), br.len(), "column count mismatch in row {i}");
        for (j, (ac, bc)) in ar.iter().zip(br.iter()).enumerate() {
            assert!(
                assert_equiv(ac.as_ref(), bc.as_ref(), ctx).expect("assert_equiv"),
                "matrix[{i},{j}] mismatch: {} vs {}",
                format_expr(ac.as_ref()),
                format_expr(bc.as_ref())
            );
        }
    }
}

pub fn assert_charpoly_equiv(m: &ExprArc, var: &Ident, expected: &Expr, ctx: &Context) {
    let cp = crate::symbolic::eval_charpoly(m, var, ctx).expect("charpoly");
    assert!(
        assert_equiv(cp.as_ref(), expected, ctx).expect("assert_equiv"),
        "charpoly mismatch: got {}",
        format_expr(cp.as_ref())
    );
}

pub fn assert_linsolve_satisfies(eqs: &ExprArc, vars: &ExprArc, sol: &ExprArc, ctx: &Context) {
    let var_list = var_idents(vars);
    let values = list_items(sol).to_vec();
    assert_eq!(
        var_list.len(),
        values.len(),
        "solution arity mismatch: {} vars, {} values",
        var_list.len(),
        values.len()
    );
    let mut subs = HashMap::new();
    for (v, val) in var_list.iter().zip(&values) {
        subs.insert(v.clone(), Arc::clone(val));
    }
    for eq in list_items(eqs) {
        let zero = equation_to_zero(eq);
        let subbed = eval_subst_map(&zero, &subs).expect("subst");
        let ev = eval(subbed.as_ref(), ctx).expect("eval");
        let n = eval(
            Expr::func(FuncKind::Normal, vec![ev]).as_ref(),
            ctx,
        )
        .expect("normal");
        let expanded = expand(n.as_ref(), ctx).expect("expand");
        assert!(
            assert_equiv(expanded.as_ref(), &Expr::int(0), ctx).expect("assert_equiv"),
            "equation not satisfied: residual {}",
            format_expr(expanded.as_ref())
        );
    }
}

fn equation_to_zero(eq: &ExprArc) -> ExprArc {
    match eq.as_ref() {
        Expr::Relation(RelOp::Eq, lhs, rhs) => Expr::add(vec![
            Arc::clone(lhs),
            Expr::mul(vec![Expr::int(-1), Arc::clone(rhs)]),
        ]),
        other => panic!("expected equation, got {other:?}"),
    }
}

pub fn assert_gramschmidt_orthogonal<F>(
    basis: &[ExprArc],
    inner: F,
    ctx: &Context,
) where
    F: Fn(&ExprArc, &ExprArc, &Context) -> Result<ExprArc, giac_core::EvalError>,
{
    for i in 0..basis.len() {
        for j in (i + 1)..basis.len() {
            let ip = inner(&basis[i], &basis[j], ctx).expect("inner product");
            assert!(
                is_zero(ip.as_ref(), ctx).expect("is_zero"),
                "basis vectors {i} and {j} not orthogonal: ⟨vi,vj⟩={}",
                format_expr(ip.as_ref())
            );
        }
    }
}

pub fn list_basis_vectors(basis_list: &ExprArc) -> Vec<ExprArc> {
    list_items(basis_list).to_vec()
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use giac_core::{Expr, Ident};

    use super::*;

    #[test]
    fn assert_matrix_equiv_smoke() {
        let ctx = crate::plugin::xcas_default();
        let a = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(0)],
            vec![Expr::int(0), Expr::int(1)],
        ]));
        let b = Arc::clone(&a);
        assert_matrix_equiv(&a, &b, &ctx);
    }

    #[test]
    fn assert_charpoly_equiv_smoke() {
        let ctx = crate::plugin::xcas_default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2)],
            vec![Expr::int(3), Expr::int(4)],
        ]));
        let x = Ident::new("x");
        let expected = Expr::add(vec![
            Expr::pow(Expr::sym("x"), Expr::int(2)),
            Expr::mul(vec![Expr::int(-5), Expr::sym("x")]),
            Expr::int(-2),
        ]);
        assert_charpoly_equiv(&m, &x, expected.as_ref(), &ctx);
    }

    #[test]
    #[ignore = "GIAC-expr-api T1/2C: normal rational sum not proved zero (linsolve residual)"]
    fn assert_linsolve_satisfies_smoke() {
        let ctx = crate::plugin::xcas_default();
        let eqs = Arc::new(Expr::List(vec![Arc::new(Expr::Relation(
            RelOp::Eq,
            Expr::add(vec![
                Expr::mul(vec![Expr::int(2), Expr::sym("x")]),
                Expr::sym("y"),
            ]),
            Expr::int(3),
        ))]));
        let vars = Arc::new(Expr::List(vec![Expr::sym("x"), Expr::sym("y")]));
        let sol = Arc::new(Expr::List(vec![
            Expr::rat(4, 3),
            Expr::rat(1, 3),
        ]));
        assert_linsolve_satisfies(&eqs, &vars, &sol, &ctx);
    }
}
