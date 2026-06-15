use std::sync::Arc;

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::error::EvalError;
use crate::expr::{Expr, ExprArc};

pub fn eval_matrix_mul(a: &ExprArc, b: &ExprArc) -> Result<ExprArc, EvalError> {
    let a_rows = as_matrix(a)?;
    let b_rows = as_matrix(b)?;
    let cols_a = a_rows[0].len();
    let rows_a = a_rows.len();
    let _rows_b = b_rows.len();
    if cols_a != b_rows.len() {
        return Err(EvalError::TypeError("incompatible matrix dimensions"));
    }
    let cols_b = b_rows[0].len();
    let mut out = Vec::with_capacity(rows_a);
    for i in 0..rows_a {
        let mut row = Vec::with_capacity(cols_b);
        for j in 0..cols_b {
            let mut sum = Expr::int(0);
            for k in 0..cols_a {
                let term = Expr::mul(vec![
                    Arc::clone(&a_rows[i][k]),
                    Arc::clone(&b_rows[k][j]),
                ]);
                sum = Expr::add(vec![sum, term]);
            }
            row.push(sum);
        }
        out.push(row);
    }
    Ok(Arc::new(Expr::Matrix(out)))
}

pub fn eval_idn(n: usize) -> ExprArc {
    let mut rows = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(n);
        for j in 0..n {
            row.push(if i == j { Expr::int(1) } else { Expr::int(0) });
        }
        rows.push(row);
    }
    Arc::new(Expr::Matrix(rows))
}

pub fn eval_tran(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    if rows.is_empty() {
        return Err(EvalError::TypeError("empty matrix"));
    }
    let nrows = rows.len();
    let ncols = rows[0].len();
    let mut out = vec![vec![Expr::int(0); nrows]; ncols];
    for i in 0..nrows {
        for j in 0..ncols {
            out[j][i] = Arc::clone(&rows[i][j]);
        }
    }
    Ok(Arc::new(Expr::Matrix(out)))
}

pub fn eval_det(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    match n {
        1 => Ok(Arc::clone(&rows[0][0])),
        2 => det2(&rows[0][0], &rows[0][1], &rows[1][0], &rows[1][1]),
        3 => det3(&rows),
        _ => Err(EvalError::NotImplemented("det for n>3")),
    }
}

fn det2(a: &ExprArc, b: &ExprArc, c: &ExprArc, d: &ExprArc) -> Result<ExprArc, EvalError> {
    Ok(Expr::add(vec![
        Expr::mul(vec![Arc::clone(a), Arc::clone(d)]),
        Expr::mul(vec![Expr::int(-1), Arc::clone(b), Arc::clone(c)]),
    ]))
}

fn det3(rows: &[Vec<ExprArc>]) -> Result<ExprArc, EvalError> {
    let mut terms = Vec::new();
    for perm in 0..6i64 {
        let (i0, i1, i2, sign) = match perm {
            0 => (0, 1, 2, 1),
            1 => (0, 2, 1, -1),
            2 => (1, 0, 2, -1),
            3 => (1, 2, 0, 1),
            4 => (2, 0, 1, 1),
            5 => (2, 1, 0, -1),
            _ => unreachable!(),
        };
        let mut f = vec![
            Arc::clone(&rows[0][i0 as usize]),
            Arc::clone(&rows[1][i1 as usize]),
            Arc::clone(&rows[2][i2 as usize]),
        ];
        if sign < 0 {
            f.insert(0, Expr::int(-1));
        }
        terms.push(Expr::mul(f));
    }
    Ok(Expr::add(terms))
}

pub fn eval_inv(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    if n == 1 {
        return Ok(Expr::pow(Arc::clone(&rows[0][0]), Expr::int(-1)));
    }
    if n == 2 {
        let det = det2(&rows[0][0], &rows[0][1], &rows[1][0], &rows[1][1])?;
        let inv_det = Expr::pow(det, Expr::int(-1));
        return Ok(Arc::new(Expr::Matrix(vec![
            vec![
                Expr::mul(vec![Arc::clone(&rows[1][1]), Arc::clone(&inv_det)]),
                Expr::mul(vec![
                    Expr::int(-1),
                    Arc::clone(&rows[0][1]),
                    Arc::clone(&inv_det),
                ]),
            ],
            vec![
                Expr::mul(vec![
                    Expr::int(-1),
                    Arc::clone(&rows[1][0]),
                    Arc::clone(&inv_det),
                ]),
                Expr::mul(vec![Arc::clone(&rows[0][0]), Arc::clone(&inv_det)]),
            ],
        ])));
    }
    Err(EvalError::NotImplemented("inv for n>2"))
}

pub fn eval_ker(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    // RREF-based null space for small rational matrices
    let (rref, pivot_cols) = rref(&rows)?;
    let ncols = rows[0].len();
    let free_cols: Vec<usize> = (0..ncols).filter(|c| !pivot_cols.contains(c)).collect();
    if free_cols.is_empty() {
        return Ok(Arc::new(Expr::Matrix(vec![])));
    }
    let mut basis = Vec::new();
    for &fc in &free_cols {
        let mut vec = vec![Expr::int(0); ncols];
        vec[fc] = Expr::int(1);
        for (row_i, piv_col) in pivot_cols.iter().enumerate() {
            if let Some(cell) = rref.get(row_i).and_then(|r| r.get(*piv_col)) {
                if let Some(c) = as_neg_rat(cell) {
                    vec[*piv_col] = c;
                }
            }
        }
        basis.push(vec);
    }
    Ok(Arc::new(Expr::Matrix(basis)))
}

pub fn eval_image(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let (rref, pivot_cols) = rref(&rows)?;
    let mut basis = Vec::new();
    for (row_i, piv_col) in pivot_cols.iter().enumerate() {
        if row_i < rref.len() {
            let mut nonzero = false;
            for c in &rref[row_i] {
                if !c.is_zero() {
                    nonzero = true;
                    break;
                }
            }
            if nonzero {
                basis.push(rref[row_i].clone());
            }
        }
        let _ = piv_col;
    }
    if basis.is_empty() {
        Ok(Arc::new(Expr::Matrix(vec![vec![Expr::int(0); rows[0].len()]])))
    } else {
        Ok(Arc::new(Expr::Matrix(basis)))
    }
}

pub fn eval_pcar(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 2 {
        let a = Arc::clone(&rows[0][0]);
        let d = Arc::clone(&rows[1][1]);
        let trace = Expr::add(vec![a, d]);
        let det = det2(&rows[0][0], &rows[0][1], &rows[1][0], &rows[1][1])?;
        let coeffs = vec![
            Expr::int(1),
            Expr::mul(vec![Expr::int(-1), trace]),
            det,
        ];
        return Ok(Expr::func(
            crate::expr::FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(coeffs))],
        ));
    }
    if n == 3 {
        return pcar_3x3(&rows);
    }
    Err(EvalError::NotImplemented("pcar n>3"))
}

fn pcar_3x3(rows: &[Vec<ExprArc>]) -> Result<ExprArc, EvalError> {
    let trace = Expr::add(vec![
        Arc::clone(&rows[0][0]),
        Arc::clone(&rows[1][1]),
        Arc::clone(&rows[2][2]),
    ]);
    let m11 = det2(&rows[1][1], &rows[1][2], &rows[2][1], &rows[2][2])?;
    let m22 = det2(&rows[0][0], &rows[0][2], &rows[2][0], &rows[2][2])?;
    let m33 = det2(&rows[0][0], &rows[0][1], &rows[1][0], &rows[1][1])?;
    let minors = Expr::add(vec![m11, m22, m33]);
    let det = det3(rows)?;
    let coeffs = vec![
        Expr::int(1),
        Expr::mul(vec![Expr::int(-1), trace]),
        minors,
        Expr::mul(vec![Expr::int(-1), det]),
    ];
    Ok(Expr::func(
        crate::expr::FuncKind::Poly1,
        vec![Arc::new(Expr::Seq(coeffs))],
    ))
}

pub fn eval_rref(args: &[ExprArc], _ctx: &crate::Context) -> Result<ExprArc, EvalError> {
    if args.is_empty() {
        return Err(EvalError::TooFewArgs("rref"));
    }
    let mut rows = Vec::new();
    let mut modulus: Option<i64> = None;
    for arg in args {
        let (row, m) = split_mod_row(arg)?;
        if let Some(prev) = modulus {
            if prev != m {
                return Err(EvalError::TypeError("modulus mismatch"));
            }
        } else if m != 0 {
            modulus = Some(m);
        }
        rows.push(row);
    }
    if let Some(m) = modulus {
        let out = mod_rref(&rows, m)?;
        return Ok(Arc::new(Expr::Matrix(out)));
    }
    let (out, _) = rref(&rows)?;
    Ok(Arc::new(Expr::Matrix(out)))
}

fn split_mod_row(e: &ExprArc) -> Result<(Vec<ExprArc>, i64), EvalError> {
    match e.as_ref() {
        Expr::Mod(inner, m) => {
            let m_i = int_from_expr(m.as_ref())?;
            Ok((row_to_vec(inner)?, m_i))
        }
        other => Ok((row_to_vec(&Arc::new(other.clone()))?, 0)),
    }
}

fn row_to_vec(e: &ExprArc) -> Result<Vec<ExprArc>, EvalError> {
    match e.as_ref() {
        Expr::List(items) | Expr::Seq(items) => Ok(items.clone()),
        Expr::Matrix(rows) if rows.len() == 1 => Ok(rows[0].clone()),
        _ => Err(EvalError::TypeError("expected row vector")),
    }
}

fn int_from_expr(e: &Expr) -> Result<i64, EvalError> {
    match e {
        Expr::Int(n) => n
            .to_string()
            .parse()
            .map_err(|_| EvalError::TypeError("integer expected")),
        _ => Err(EvalError::TypeError("integer expected")),
    }
}

fn mod_rref(rows: &[Vec<ExprArc>], modulus: i64) -> Result<Vec<Vec<ExprArc>>, EvalError> {
    let mut m: Vec<Vec<i64>> = rows
        .iter()
        .map(|row| row.iter().map(|c| int_from_expr(c.as_ref())).collect())
        .collect::<Result<_, _>>()?;
    let ncols = m[0].len();
    let mut pivot_row = 0usize;
    for col in 0..ncols {
        if pivot_row >= m.len() {
            break;
        }
        let mut sel = None;
        for r in pivot_row..m.len() {
            if m[r][col] % modulus != 0 {
                sel = Some(r);
                break;
            }
        }
        let Some(r) = sel else { continue };
        m.swap(pivot_row, r);
        let inv = mod_inv(m[pivot_row][col], modulus)?;
        for c in 0..ncols {
            m[pivot_row][c] = smod_i64(m[pivot_row][c] * inv, modulus);
        }
        for r in 0..m.len() {
            if r == pivot_row || m[r][col] == 0 {
                continue;
            }
            let factor = m[r][col];
            for c in 0..ncols {
                m[r][c] = smod_i64(m[r][c] - factor * m[pivot_row][c], modulus);
            }
        }
        pivot_row += 1;
    }
    Ok(m.into_iter()
        .map(|row| row.into_iter().map(Expr::int).collect())
        .collect())
}

fn smod_i64(a: i64, m: i64) -> i64 {
    giac_poly::smod(a, m)
}

fn mod_inv(a: i64, m: i64) -> Result<i64, EvalError> {
    let eg = BigInt::from(a).extended_gcd(&BigInt::from(m));
    if eg.gcd != BigInt::one() {
        return Err(EvalError::TypeError("not invertible"));
    }
    Ok(smod_i64(
        eg.x.to_string().parse().unwrap_or(0),
        m,
    ))
}

fn rref(rows: &[Vec<ExprArc>]) -> Result<(Vec<Vec<ExprArc>>, Vec<usize>), EvalError> {
    let mut m = rows.to_vec();
    let ncols = m[0].len();
    let mut pivot_cols = Vec::new();
    let mut pivot_row = 0usize;
    for col in 0..ncols {
        if pivot_row >= m.len() {
            break;
        }
        let mut sel = None;
        for r in pivot_row..m.len() {
            if !is_zero_expr(&m[r][col]) {
                sel = Some(r);
                break;
            }
        }
        let Some(r) = sel else { continue };
        m.swap(pivot_row, r);
        pivot_cols.push(col);
        let pivot_val = Arc::clone(&m[pivot_row][col]);
        for c in 0..ncols {
            m[pivot_row][c] = div_expr(&m[pivot_row][c], &pivot_val)?;
        }
        for r in 0..m.len() {
            if r == pivot_row {
                continue;
            }
            if !is_zero_expr(&m[r][col]) {
                let factor = Arc::clone(&m[r][col]);
                for c in 0..ncols {
                    let sub = Expr::mul(vec![factor.clone(), Arc::clone(&m[pivot_row][c])]);
                    m[r][c] = Expr::add(vec![Arc::clone(&m[r][c]), Expr::mul(vec![Expr::int(-1), sub])]);
                }
            }
        }
        pivot_row += 1;
    }
    Ok((m, pivot_cols))
}

fn div_expr(a: &ExprArc, b: &ExprArc) -> Result<ExprArc, EvalError> {
    Ok(Expr::mul(vec![Arc::clone(a), Expr::pow(Arc::clone(b), Expr::int(-1))]))
}

fn is_zero_expr(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
        || matches!(e.as_ref(), Expr::Rat(r) if r.is_zero())
}

fn as_neg_rat(e: &ExprArc) -> Option<ExprArc> {
    match e.as_ref() {
        Expr::Int(n) => Some(Expr::int(-n.to_string().parse::<i64>().ok()?)),
        Expr::Rat(r) => Some(Arc::new(Expr::Rat(-r.clone()))),
        _ => None,
    }
}

pub fn as_matrix(e: &ExprArc) -> Result<Vec<Vec<ExprArc>>, EvalError> {
    match e.as_ref() {
        Expr::Matrix(rows) | Expr::GiacMatrix(rows) => Ok(rows.clone()),
        _ => Err(EvalError::TypeError("expected matrix")),
    }
}

#[allow(dead_code)] // exercised from matrix unit tests
pub fn is_identity_matrix(m: &ExprArc) -> bool {
    let Ok(rows) = as_matrix(m) else {
        return false;
    };
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return false;
    }
    for (i, row) in rows.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            let want_one = i == j;
            match c.as_ref() {
                Expr::Int(n) if want_one && n == &BigInt::one() => {}
                Expr::Int(n) if !want_one && n.is_zero() => {}
                _ if want_one && c.is_one() => {}
                _ if !want_one && c.is_zero() => {}
                _ => return false,
            }
        }
    }
    true
}

pub fn inv_scalar(n: &BigInt) -> Result<ExprArc, EvalError> {
    if n.is_zero() {
        return Err(EvalError::DivisionByZero);
    }
    Ok(Expr::rat(
        1,
        n.to_string()
            .parse()
            .map_err(|_| EvalError::TypeError("scalar inv denominator"))?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eval::eval;
    use crate::format_expr;
    use crate::Context;

    #[test]
    fn idn_2x2() {
        let m = eval_idn(2);
        assert!(is_identity_matrix(&m));
    }

    #[test]
    fn det_2x2() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2)],
            vec![Expr::int(3), Expr::int(4)],
        ]));
        let d = eval_det(&m).unwrap();
        let r = eval(d.as_ref(), &ctx).unwrap();
        assert_eq!(r, Expr::int(-2));
    }

    #[test]
    fn pcar_3x3_singular() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2), Expr::int(3)],
            vec![Expr::int(1), Expr::int(3), Expr::int(6)],
            vec![Expr::int(2), Expr::int(5), Expr::int(9)],
        ]));
        let p = eval_pcar(&m).unwrap();
        let r = eval(p.as_ref(), &ctx).unwrap();
        assert_eq!(
            crate::format_expr(r.as_ref()),
            "poly1[1,-13,1,0]"
        );
    }

    #[test]
    fn matrix_mul_incompatible_dimensions() {
        let a = Arc::new(Expr::Matrix(vec![vec![Expr::int(1), Expr::int(2)]]));
        let b = Arc::new(Expr::Matrix(vec![vec![Expr::int(1)]]));
        assert!(eval_matrix_mul(&a, &b).is_err());
    }

    #[test]
    fn det_1x1_and_3x3() {
        let ctx = Context::default();
        let m1 = Arc::new(Expr::Matrix(vec![vec![Expr::int(7)]]));
        assert_eq!(
            eval(eval_det(&m1).unwrap().as_ref(), &ctx).unwrap(),
            Expr::int(7)
        );

        let m3 = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(2), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(3)],
        ]));
        assert_eq!(
            eval(eval_det(&m3).unwrap().as_ref(), &ctx).unwrap(),
            Expr::int(6)
        );
    }

    #[test]
    fn matrix_error_paths() {
        let bad = Arc::new(Expr::int(1));
        assert!(as_matrix(&bad).is_err());
        assert!(eval_tran(&bad).is_err());
        assert!(eval_det(&bad).is_err());
        assert!(eval_inv(&bad).is_err());

        let empty = Arc::new(Expr::Matrix(vec![]));
        assert!(eval_tran(&empty).is_err());

        let nonsquare = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2)],
            vec![Expr::int(3), Expr::int(4)],
            vec![Expr::int(5), Expr::int(6)],
        ]));
        assert!(eval_det(&nonsquare).is_err());

        let big = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2), Expr::int(3)],
            vec![Expr::int(4), Expr::int(5), Expr::int(6)],
            vec![Expr::int(7), Expr::int(8), Expr::int(9)],
            vec![Expr::int(0), Expr::int(1), Expr::int(2)],
        ]));
        assert!(matches!(
            eval_inv(&big),
            Err(EvalError::TypeError(_))
        ));
        assert!(matches!(
            eval_pcar(&big),
            Err(EvalError::NotImplemented(_))
        ));
    }

    #[test]
    fn inv_scalar_and_identity_check() {
        assert!(inv_scalar(&BigInt::zero()).is_err());
        assert_eq!(inv_scalar(&BigInt::from(5)).unwrap(), Expr::rat(1, 5));
        let m = eval_idn(3);
        assert!(is_identity_matrix(&m));
        let not_id = Arc::new(Expr::Matrix(vec![vec![Expr::int(2)]]));
        assert!(!is_identity_matrix(&not_id));
    }

    #[test]
    fn pcar_2x2() {
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2)],
            vec![Expr::int(3), Expr::int(4)],
        ]));
        let p = eval_pcar(&m).unwrap();
        assert!(matches!(p.as_ref(), Expr::Func(crate::expr::FuncKind::Poly1, _)));
    }

    #[test]
    fn ker_full_rank_and_image() {
        let full = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(0)],
            vec![Expr::int(0), Expr::int(1)],
        ]));
        let ker = eval_ker(&full).unwrap();
        assert!(matches!(ker.as_ref(), Expr::Matrix(rows) if rows.is_empty()));

        let image = eval_image(&full).unwrap();
        match image.as_ref() {
            Expr::Matrix(rows) => assert_eq!(rows.len(), 2),
            _ => panic!("expected matrix image"),
        }
    }

    #[test]
    fn ker_nontrivial_nullspace() {
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2), Expr::int(3)],
            vec![Expr::int(2), Expr::int(4), Expr::int(6)],
        ]));
        let ker = eval_ker(&m).unwrap();
        match ker.as_ref() {
            Expr::Matrix(rows) => assert!(!rows.is_empty()),
            _ => panic!("expected kernel basis"),
        }
    }

    #[test]
    fn det_and_inv_sizes() {
        let ctx = Context::default();
        let m1 = Arc::new(Expr::Matrix(vec![vec![Expr::int(4)]]));
        let inv1 = eval_inv(&m1).unwrap();
        assert_eq!(format_expr(eval(inv1.as_ref(), &ctx).unwrap().as_ref()), "4^-1");

        let m4 = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(0), Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(1), Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(1), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(0), Expr::int(1)],
        ]));
        assert!(matches!(
            eval_det(&m4),
            Err(EvalError::NotImplemented(_))
        ));
        assert!(matches!(
            eval_inv(&Arc::new(Expr::Matrix(vec![
                vec![Expr::int(1), Expr::int(0), Expr::int(0)],
                vec![Expr::int(0), Expr::int(1), Expr::int(0)],
                vec![Expr::int(0), Expr::int(0), Expr::int(1)],
            ]))),
            Err(EvalError::NotImplemented(_))
        ));
    }

    #[test]
    fn image_zero_matrix() {
        let z = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0)],
        ]));
        let image = eval_image(&z).unwrap();
        match image.as_ref() {
            Expr::Matrix(rows) => {
                assert_eq!(rows.len(), 1);
                assert!(rows[0].iter().all(|c| c.is_zero()));
            }
            _ => panic!("expected matrix"),
        }
    }

    #[test]
    fn identity_matrix_rat_one() {
        let m = Arc::new(Expr::Matrix(vec![vec![Expr::rat(2, 2)]]));
        assert!(is_identity_matrix(&m));
    }
}
