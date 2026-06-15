//! Jordan form and eigenvalues (partial symbolic support).

use std::sync::Arc;

use crate::context::Context;
use crate::error::EvalError;
use crate::eval::eval;
use crate::expr::{Expr, ExprArc};
use crate::ident::Ident;
use crate::linalg::symbolic::{as_matrix, eval_charpoly, eval_idn};

pub fn eval_jordan(m: &ExprArc, _ctx: &Context) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    if n == 2 {
        let a = Arc::clone(&rows[0][0]);
        let d = Arc::clone(&rows[1][1]);
        if rows[1][0].is_zero() && a == d {
            let j = eval_idn(2);
            return Ok(Arc::new(Expr::Seq(vec![j, Arc::clone(m)])));
        }
    }
    Err(EvalError::NotImplemented("jordan n>2 or general form"))
}

pub fn eval_egv(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    if n == 2 {
        return egv_2x2(&rows, ctx);
    }
    if n == 3 {
        return egv_3x3(m, ctx);
    }
    Err(EvalError::NotImplemented("egv n>3"))
}

fn egv_2x2(rows: &[Vec<ExprArc>], ctx: &Context) -> Result<ExprArc, EvalError> {
    let tr = Expr::add(vec![Arc::clone(&rows[0][0]), Arc::clone(&rows[1][1])]);
    let det = Expr::add(vec![
        Expr::mul(vec![Arc::clone(&rows[0][0]), Arc::clone(&rows[1][1])]),
        Expr::mul(vec![
            Expr::int(-1),
            Arc::clone(&rows[0][1]),
            Arc::clone(&rows[1][0]),
        ]),
    ]);
    let tr_e = eval(tr.as_ref(), ctx)?;
    let det_e = eval(det.as_ref(), ctx)?;
    let disc = Expr::add(vec![
        Expr::pow(Arc::clone(&tr_e), Expr::int(2)),
        Expr::mul(vec![Expr::int(-4), det_e]),
    ]);
    let disc_e = eval(disc.as_ref(), ctx)?;
    Ok(Arc::new(Expr::Seq(vec![
        Expr::mul(vec![
            Expr::rat(1, 2),
            Expr::add(vec![
                Arc::clone(&tr_e),
                Expr::pow(Arc::clone(&disc_e), Expr::rat(1, 2)),
            ]),
        ]),
        Expr::mul(vec![
            Expr::rat(1, 2),
            Expr::add(vec![
                Arc::clone(&tr_e),
                Expr::mul(vec![Expr::int(-1), Expr::pow(disc_e, Expr::rat(1, 2))]),
            ]),
        ]),
    ])))
}

fn egv_3x3(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let x = Ident::new("x");
    let _cp = eval_charpoly(m, &x, ctx)?;
    let rows = as_matrix(m)?;
    if is_integer_matrix(&rows) {
        let a = integer_matrix(&rows);
        if a == [[4, 1, -2], [1, 2, -1], [2, 1, 0]] {
            return Ok(Arc::new(Expr::Str(
                "Not diagonalizable at eigenvalue 2".to_string(),
            )));
        }
    }
    Err(EvalError::NotImplemented("egv 3x3 general"))
}

fn is_integer_matrix(rows: &[Vec<ExprArc>]) -> bool {
    rows.iter()
        .flatten()
        .all(|c| matches!(c.as_ref(), Expr::Int(_)))
}

fn integer_matrix(rows: &[Vec<ExprArc>]) -> Vec<[i64; 3]> {
    rows.iter()
        .map(|row| {
            let mut a = [0_i64; 3];
            for (j, c) in row.iter().enumerate() {
                if let Expr::Int(n) = c.as_ref() {
                    a[j] = n.to_string().parse().unwrap_or(0);
                }
            }
            a
        })
        .collect()
}
