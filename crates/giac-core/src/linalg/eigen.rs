//! Jordan form and eigenvalues.

use std::collections::HashMap;
use std::sync::Arc;

use giac_linalg::real_eigenvalues;

use crate::context::Context;
use crate::error::EvalError;
use crate::eval::eval;
use crate::expr::{Expr, ExprArc};
use crate::ident::Ident;
use crate::num_util::bigint_to_i64;
use crate::linalg::symbolic::{as_matrix, eval_charpoly, eval_idn, try_to_f64_matrix};

type Mat3 = [[i64; 3]; 3];

pub fn eval_jordan(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    match n {
        2 => jordan_2x2(m, &rows, ctx),
        3 => jordan_3x3(m, &rows, ctx),
        _ => Err(EvalError::NotImplemented("jordan n>3")),
    }
}

pub fn eval_egv(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    match n {
        1 => Ok(Arc::clone(&rows[0][0])),
        2 => egv_2x2(&rows, ctx),
        3 => egv_3x3(m, &rows, ctx),
        _ => Err(EvalError::NotImplemented("egv n>3")),
    }
}

fn jordan_2x2(m: &ExprArc, rows: &[Vec<ExprArc>], _ctx: &Context) -> Result<ExprArc, EvalError> {
    if is_integer_matrix(rows) {
        let a = int_mat2(rows);
        if is_jordan_form_2x2(&a) {
            return Ok(Arc::new(Expr::Seq(vec![Arc::clone(m), eval_idn(2)])));
        }
    }
    let a00 = &rows[0][0];
    let a11 = &rows[1][1];
    if rows[1][0].is_zero() && a00 == a11 {
        return Ok(Arc::new(Expr::Seq(vec![Arc::clone(m), eval_idn(2)])));
    }
    Err(EvalError::NotImplemented("jordan 2x2 general"))
}

fn jordan_3x3(m: &ExprArc, rows: &[Vec<ExprArc>], _ctx: &Context) -> Result<ExprArc, EvalError> {
    if is_integer_matrix(rows) {
        let a = int_mat3(rows);
        if is_jordan_form_3x3(&a) {
            return Ok(Arc::new(Expr::Seq(vec![Arc::clone(m), eval_idn(3)])));
        }
    }
    Err(EvalError::NotImplemented("jordan 3x3 general"))
}

fn is_jordan_form_2x2(a: &[[i64; 2]; 2]) -> bool {
    a[1][0] == 0 && (a[0][1] == 0 || a[0][1] == 1)
}

fn is_jordan_form_3x3(a: &Mat3) -> bool {
    for i in 0..3 {
        for j in 0..i {
            if a[i][j] != 0 {
                return false;
            }
        }
        for j in (i + 2)..3 {
            if a[i][j] != 0 {
                return false;
            }
        }
    }
    for i in 0..2 {
        let s = a[i][i + 1];
        if s != 0 && s != 1 {
            return false;
        }
        if s == 1 && a[i][i] != a[i + 1][i + 1] {
            return false;
        }
    }
    true
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

fn egv_3x3(m: &ExprArc, rows: &[Vec<ExprArc>], ctx: &Context) -> Result<ExprArc, EvalError> {
    if is_integer_matrix(rows) {
        return egv_int3(int_mat3(rows), ctx);
    }
    if let Ok(a) = try_to_f64_matrix(m, ctx) {
        if let Some(vals) = real_eigenvalues(&a) {
            let ev: Vec<ExprArc> = vals
                .into_iter()
                .map(f64_to_expr_sorted)
                .collect();
            return Ok(Arc::new(Expr::Seq(ev)));
        }
    }
    let x = Ident::new("x");
    let _cp = eval_charpoly(m, &x, ctx)?;
    Err(EvalError::NotImplemented("egv 3x3 general"))
}

fn egv_int3(a: Mat3, ctx: &Context) -> Result<ExprArc, EvalError> {
    if is_diagonal_int3(&a) {
        return Ok(Arc::new(Expr::Seq(vec![
            Expr::int(a[0][0]),
            Expr::int(a[1][1]),
            Expr::int(a[2][2]),
        ])));
    }

    let (c2, c1, c0) = charpoly_coeffs_int3(&a);
    let roots = cubic_integer_roots(c2, c1, c0)
        .ok_or(EvalError::NotImplemented("egv 3x3 roots"))?;

    if let Some(lam) = non_diagonalizable_eigenvalue(&a, &roots) {
        return Ok(Arc::new(Expr::Str(format!(
            "Not diagonalizable at eigenvalue {lam}"
        ))));
    }

    let ev: Vec<ExprArc> = roots.into_iter().map(Expr::int).collect();
    let _ = ctx;
    Ok(Arc::new(Expr::Seq(ev)))
}

fn is_diagonal_int3(a: &Mat3) -> bool {
    (0..3).all(|i| (0..3).all(|j| i == j || a[i][j] == 0))
}

fn charpoly_coeffs_int3(a: &Mat3) -> (i64, i64, i64) {
    let c2 = -(a[0][0] + a[1][1] + a[2][2]);
    let c1 = det2_int(a[0][0], a[0][1], a[1][0], a[1][1])
        + det2_int(a[0][0], a[0][2], a[2][0], a[2][2])
        + det2_int(a[1][1], a[1][2], a[2][1], a[2][2]);
    let c0 = -det3_int(a);
    (c2, c1, c0)
}

fn det2_int(a: i64, b: i64, c: i64, d: i64) -> i64 {
    a * d - b * c
}

fn det3_int(a: &Mat3) -> i64 {
    a[0][0] * det2_int(a[1][1], a[1][2], a[2][1], a[2][2])
        - a[0][1] * det2_int(a[1][0], a[1][2], a[2][0], a[2][2])
        + a[0][2] * det2_int(a[1][0], a[1][1], a[2][0], a[2][1])
}

fn cubic_eval(c2: i64, c1: i64, c0: i64, r: i64) -> i64 {
    r * r * r + c2 * r * r + c1 * r + c0
}

fn cubic_integer_roots(c2: i64, c1: i64, c0: i64) -> Option<Vec<i64>> {
    if c0 == 0 {
        let mut roots = vec![0_i64];
        let (q1, q0) = (c2, c1);
        roots.extend(quadratic_integer_roots(1, q1, q0)?);
        roots.sort();
        return Some(roots);
    }

    let divisors = positive_divisors(c0.unsigned_abs());
    let mut candidates: Vec<i64> = divisors
        .iter()
        .flat_map(|&d| {
            let d = d as i64;
            if d == 0 {
                vec![0]
            } else {
                vec![d, -d]
            }
        })
        .collect();
    candidates.sort();
    candidates.dedup();

    for &r in &candidates {
        if r != 0 && cubic_eval(c2, c1, c0, r) == 0 {
            let (q2, q1, q0) = synthetic_div_cubic(c2, c1, c0, r);
            let mut roots = vec![r];
            roots.extend(quadratic_integer_roots(q2, q1, q0)?);
            roots.sort();
            return Some(roots);
        }
    }
    None
}

fn synthetic_div_cubic(c2: i64, c1: i64, c0: i64, r: i64) -> (i64, i64, i64) {
    let q2 = 1_i64;
    let q1 = c2 + r;
    let q0 = c1 + r * q1;
    let _ = c0 + r * q0;
    (q2, q1, q0)
}

fn quadratic_integer_roots(a: i64, b: i64, c: i64) -> Option<Vec<i64>> {
    if a == 0 {
        if b == 0 {
            return if c == 0 { Some(vec![]) } else { None };
        }
        if c % b != 0 {
            return None;
        }
        return Some(vec![-c / b]);
    }
    let disc = b * b - 4 * a * c;
    if disc < 0 {
        return None;
    }
    let s = integer_sqrt(disc)?;
    if (b - s) % (2 * a) != 0 || (b + s) % (2 * a) != 0 {
        return None;
    }
    let r1 = (-b + s) / (2 * a);
    let r2 = (-b - s) / (2 * a);
    let mut roots = vec![r1, r2];
    roots.sort();
    Some(roots)
}

fn integer_sqrt(n: i64) -> Option<i64> {
    if n < 0 {
        return None;
    }
    let r = (n as f64).sqrt().round() as i64;
    if r * r == n {
        Some(r)
    } else {
        None
    }
}

fn positive_divisors(n: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut d = 1;
    while d * d <= n {
        if n % d == 0 {
            out.push(d);
            if d * d != n {
                out.push(n / d);
            }
        }
        d += 1;
    }
    out.sort();
    out
}

fn mat3_sub_lambda(a: &Mat3, lam: i64) -> Mat3 {
    let mut m = *a;
    m[0][0] -= lam;
    m[1][1] -= lam;
    m[2][2] -= lam;
    m
}

fn mat3_rank(mut m: Mat3) -> usize {
    let mut rank = 0usize;
    let mut col = 0usize;
    while rank < 3 && col < 3 {
        let mut pivot = rank;
        while pivot < 3 && m[pivot][col] == 0 {
            pivot += 1;
        }
        if pivot == 3 {
            col += 1;
            continue;
        }
        m.swap(rank, pivot);
        let pv = m[rank][col];
        for r in (rank + 1)..3 {
            if m[r][col] != 0 {
                let factor = m[r][col] / pv;
                for c in col..3 {
                    m[r][c] -= factor * m[rank][c];
                }
            }
        }
        rank += 1;
        col += 1;
    }
    rank
}

fn non_diagonalizable_eigenvalue(a: &Mat3, roots: &[i64]) -> Option<i64> {
    let mut alg: HashMap<i64, usize> = HashMap::new();
    for &r in roots {
        *alg.entry(r).or_default() += 1;
    }
    for (&lam, &multiplicity) in &alg {
        let geom = 3 - mat3_rank(mat3_sub_lambda(a, lam));
        if geom < multiplicity {
            return Some(lam);
        }
    }
    None
}

fn f64_to_expr_sorted(v: f64) -> ExprArc {
    if (v - v.round()).abs() < 1e-9 {
        return Expr::int(v.round() as i64);
    }
    crate::linalg::symbolic::f64_to_expr_numeric(v)
}

fn is_integer_matrix(rows: &[Vec<ExprArc>]) -> bool {
    rows.iter()
        .flatten()
        .all(|c| matches!(c.as_ref(), Expr::Int(_)))
}

fn int_mat2(rows: &[Vec<ExprArc>]) -> [[i64; 2]; 2] {
    let mut a = [[0_i64; 2]; 2];
    for (i, row) in rows.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if let Expr::Int(n) = c.as_ref() {
                a[i][j] = bigint_to_i64(n).unwrap_or(0);
            }
        }
    }
    a
}

fn int_mat3(rows: &[Vec<ExprArc>]) -> Mat3 {
    let mut a = [[0_i64; 3]; 3];
    for (i, row) in rows.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if let Expr::Int(n) = c.as_ref() {
                a[i][j] = bigint_to_i64(n).unwrap_or(0);
            }
        }
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format_expr;

    #[test]
    fn egv_3x3_not_diagonalizable() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(4), Expr::int(1), Expr::int(-2)],
            vec![Expr::int(1), Expr::int(2), Expr::int(-1)],
            vec![Expr::int(2), Expr::int(1), Expr::int(0)],
        ]));
        let r = eval_egv(&m, &ctx).unwrap();
        assert_eq!(
            format_expr(r.as_ref()),
            "\"Not diagonalizable at eigenvalue 2\""
        );
    }

    #[test]
    fn egv_3x3_diagonal() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(2), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(3)],
        ]));
        let r = eval_egv(&m, &ctx).unwrap();
        assert_eq!(format_expr(r.as_ref()), "1,2,3");
    }

    #[test]
    fn jordan_3x3_block() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(1), Expr::int(0)],
            vec![Expr::int(0), Expr::int(1), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(2)],
        ]));
        let r = eval_jordan(&m, &ctx).unwrap();
        let s = format_expr(r.as_ref());
        assert!(s.contains("[[1,1,0]"), "got {s}");
        assert!(s.contains("[[1,0,0]"), "got {s}");
    }

    #[test]
    fn egv_3x3_jordan_block_eigenvalue() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(1), Expr::int(0)],
            vec![Expr::int(0), Expr::int(1), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(2)],
        ]));
        let r = eval_egv(&m, &ctx).unwrap();
        assert_eq!(
            format_expr(r.as_ref()),
            "\"Not diagonalizable at eigenvalue 1\""
        );
    }

    #[test]
    fn cubic_roots_example() {
        let a = [[4, 1, -2], [1, 2, -1], [2, 1, 0]];
        let (c2, c1, c0) = charpoly_coeffs_int3(&a);
        let roots = cubic_integer_roots(c2, c1, c0).unwrap();
        assert_eq!(roots, vec![2, 2, 2]);
    }
}
