//! Symbolic matrix algorithms (`Expr::Matrix`).

use std::collections::HashMap;
use std::sync::Arc;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Zero};

use giac_core::normal;
use giac_core::Context;
use giac_core::EvalError;
use giac_core::eval;
use giac_core::{Expr, ExprArc, FuncKind, RelOp};
use giac_core::Ident;

pub fn as_matrix(e: &ExprArc) -> Result<Vec<Vec<ExprArc>>, EvalError> {
    match e.as_ref() {
        Expr::Matrix(rows) | Expr::GiacMatrix(rows) => Ok(rows.clone()),
        _ => Err(EvalError::TypeError("expected matrix")),
    }
}

pub fn eval_matrix_mul(a: &ExprArc, b: &ExprArc) -> Result<ExprArc, EvalError> {
    let a_rows = as_matrix(a)?;
    let b_rows = as_matrix(b)?;
    let cols_a = a_rows[0].len();
    let rows_a = a_rows.len();
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

pub fn eval_matrix_pow(m: &ExprArc, exp: u32, ctx: &Context) -> Result<ExprArc, EvalError> {
    if exp == 0 {
        let n = as_matrix(m)?.len();
        return Ok(eval_idn(n));
    }
    let mut result = Arc::clone(m);
    for _ in 1..exp {
        result = eval_matrix_mul(&result, m)?;
        result = eval(result.as_ref(), ctx)?;
    }
    Ok(result)
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
    Arc::new(Expr::GiacMatrix(rows))
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
    Ok(Arc::new(Expr::GiacMatrix(out)))
}

pub fn eval_trace(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    let terms: Vec<_> = (0..n).map(|i| Arc::clone(&rows[i][i])).collect();
    Ok(Expr::add(terms))
}

pub fn eval_det(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    let raw = bareiss_det(&rows)?;
    let ev = eval(raw.as_ref(), ctx)?;
    normal(ev.as_ref(), ctx)
}

fn bareiss_det(m: &[Vec<ExprArc>]) -> Result<ExprArc, EvalError> {
    let n = m.len();
    if n == 1 {
        return Ok(Arc::clone(&m[0][0]));
    }
    let mut a: Vec<Vec<ExprArc>> = m.to_vec();
    let mut prev = Expr::int(1);
    for k in 0..n - 1 {
        if is_zero_expr(&a[k][k]) {
            let mut swap_row = None;
            for i in (k + 1)..n {
                if !is_zero_expr(&a[i][k]) {
                    swap_row = Some(i);
                    break;
                }
            }
            let Some(r) = swap_row else {
                return Ok(Expr::int(0));
            };
            a.swap(k, r);
            prev = Expr::int(-1);
        }
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let num = Expr::add(vec![
                    Expr::mul(vec![Arc::clone(&a[i][j]), Arc::clone(&a[k][k])]),
                    Expr::mul(vec![
                        Expr::int(-1),
                        Arc::clone(&a[i][k]),
                        Arc::clone(&a[k][j]),
                    ]),
                ]);
                a[i][j] = if is_one_expr(&prev) {
                    num
                } else {
                    Expr::mul(vec![Arc::clone(&prev), num])
                };
            }
        }
        if k < n - 2 {
            let div = Arc::clone(&a[k][k]);
            for i in (k + 1)..n {
                for j in (k + 1)..n {
                    a[i][j] = div_expr(&a[i][j], &div)?;
                }
            }
        }
        prev = Arc::clone(&a[k][k]);
    }
    Ok(Arc::clone(&a[n - 1][n - 1]))
}

pub fn eval_inv(m: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
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
    let mut aug = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = rows[i].clone();
        for j in 0..n {
            row.push(if i == j { Expr::int(1) } else { Expr::int(0) });
        }
        aug.push(row);
    }
    let (rref, _) = rref(&aug)?;
    let det_pivot = (0..n).all(|i| !is_zero_expr(&rref[i][i]));
    if !det_pivot {
        return Err(EvalError::TypeError("singular matrix"));
    }
    let inv: Vec<Vec<ExprArc>> = (0..n).map(|i| rref[i][n..2 * n].to_vec()).collect();
    let out = Arc::new(Expr::Matrix(inv));
    eval(out.as_ref(), ctx)
}

pub fn eval_ker(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
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
        for (row_i, &piv_col) in pivot_cols.iter().enumerate() {
            if let Some(cell) = rref.get(row_i).and_then(|r| r.get(fc)) {
                if !is_zero_expr(cell) {
                    vec[piv_col] = neg_expr(cell);
                }
            }
        }
        basis.push(vec);
    }
    Ok(Arc::new(Expr::Matrix(basis)))
}

pub fn eval_image(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let (_, pivot_cols) = rref(&rows)?;
    let mut basis = Vec::new();
    for &col in &pivot_cols {
        let col_vec: Vec<ExprArc> = rows.iter().map(|row| Arc::clone(&row[col])).collect();
        basis.push(col_vec);
    }
    if basis.is_empty() {
        Ok(Arc::new(Expr::Matrix(vec![vec![Expr::int(0); rows[0].len()]])))
    } else {
        Ok(Arc::new(Expr::Matrix(basis)))
    }
}

pub fn eval_pcar(m: &ExprArc) -> Result<ExprArc, EvalError> {
    charpoly_coeffs(m)
}

pub fn eval_charpoly(m: &ExprArc, var: &Ident, ctx: &Context) -> Result<ExprArc, EvalError> {
    let coeffs = charpoly_coeffs(m)?;
    let mut terms = Vec::new();
    let Expr::Func(FuncKind::Poly1, args) = coeffs.as_ref() else {
        return Err(EvalError::TypeError("charpoly internal error"));
    };
    let Some(Expr::Seq(c)) = args.first().map(|a| a.as_ref()) else {
        return Err(EvalError::TypeError("charpoly internal error"));
    };
    for (i, coeff) in c.iter().enumerate() {
        if is_zero_expr(coeff) {
            continue;
        }
        let power = c.len() - 1 - i;
        let x = Arc::new(Expr::Symbol(var.clone()));
        let term = if power == 0 {
            Arc::clone(coeff)
        } else if power == 1 {
            Expr::mul(vec![Arc::clone(coeff), x])
        } else {
            Expr::mul(vec![Arc::clone(coeff), Expr::pow(x, Expr::int(power as i64))])
        };
        terms.push(term);
    }
    if terms.is_empty() {
        return Ok(Expr::int(0));
    }
    normal(Expr::add(terms).as_ref(), ctx)
}

fn charpoly_coeffs(m: &ExprArc) -> Result<ExprArc, EvalError> {
    let rows = as_matrix(m)?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix"));
    }
    if n == 1 {
        return Ok(Expr::func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                neg_expr(&rows[0][0]),
            ]))],
        ));
    }
    if n == 2 {
        let trace = Expr::add(vec![Arc::clone(&rows[0][0]), Arc::clone(&rows[1][1])]);
        let det = bareiss_det(&rows)?;
        return Ok(Expr::func(
            FuncKind::Poly1,
            vec![Arc::new(Expr::Seq(vec![
                Expr::int(1),
                neg_expr(&trace),
                det,
            ]))],
        ));
    }
    if n == 3 {
        return pcar_3x3(&rows);
    }
    // Berkowitz algorithm: division-free, works for arbitrary n
    berkowitz_charpoly(&rows)
}

/// Berkowitz algorithm for characteristic polynomial coefficients.
///
/// For an n×n matrix A, computes the coefficients [c_0, c_1, ..., c_n]
/// of the characteristic polynomial det(λI - A) = Σ c_i λ^i,
/// where c_n = 1 (monic).
///
/// The algorithm uses the recurrence:
///   M_k = A[1..k, 1..k]  (k×k leading principal submatrix)
///   c^{(k)}_0 = 1
///   c^{(k)}_i = -(1/i) Σ_{j=1}^{i} tr(M^{j}_k) · c^{(k)}_{i-j}   for i = 1..k
///
/// Implemented as the Toeplitz matrix-vector product (Berkowitz's formulation).
fn berkowitz_charpoly(rows: &[Vec<ExprArc>]) -> Result<ExprArc, EvalError> {
    let n = rows.len();
    // Build the Berkowitz "T" vector: T[j] = -tr(A^{j+1}) for j = 0..n-1
    // And the "C" matrix powers: C[j] = A^{j+1} for j = 0..n-1
    //
    // Simpler approach: use the Faddeev–LeVerrier recurrence directly.
    // c_n = 1
    // c_{n-k} = -(1/k) * Σ_{j=1}^{k} tr(A^j) * c_{n-k+j}   for k = 1..n
    //
    // This gives coefficients [c_0, c_1, ..., c_n] with c_n = 1.

    let mut traces: Vec<ExprArc> = Vec::with_capacity(n);
    // Compute powers of A and their traces
    // A^1
    let mut power = rows.to_vec();
    traces.push(trace_of_rows(&power)?);
    // A^2 .. A^n
    for _ in 1..n {
        power = mul_row_matrices(&power, rows)?;
        traces.push(trace_of_rows(&power)?);
    }

    // Faddeev–LeVerrier: compute coefficients from c_0 to c_n
    // c[0] = (-1)^n * det(A) — but we build from the top
    // Actually use the standard recurrence:
    //   p_0 = 1
    //   p_k = -1/k * Σ_{i=1}^{k} tr(A^i) * p_{k-i}   for k = 1..n
    // Then char poly = Σ p_{n-k} * x^k = p_n*x^0 + p_{n-1}*x^1 + ... + p_0*x^n
    // Wait, let's be precise:
    // The characteristic polynomial det(λI - A) = λ^n + p_1*λ^{n-1} + ... + p_n
    // p_0 = 1
    // p_k = -(1/k) * Σ_{i=1}^{k} tr(A^i) * p_{k-i}

    let mut p = Vec::with_capacity(n + 1);
    p.push(Expr::int(1)); // p_0

    for k in 1..=n {
        let mut sum = Expr::int(0);
        for i in 1..=k {
            // sum += tr(A^i) * p[k-i]
            let term = Expr::mul(vec![
                Arc::clone(&traces[i - 1]),
                Arc::clone(&p[k - i]),
            ]);
            sum = Expr::add(vec![sum, term]);
        }
        // p_k = -sum / k
        // Since all our traces and p's are symbolic, we can't easily divide by k.
        // Instead, use the division-free Berkowitz form by keeping the sign:
        // p_k = -(1/k) * sum — but this requires rational arithmetic.
        // Since our Expr supports Rat, we produce: p_k = -sum * (1/k)
        let coeff = Expr::mul(vec![
            Expr::int(-1),
            sum,
            Expr::rat(1, k as i64),
        ]);
        p.push(coeff);
    }

    // char poly det(λI - A) = p_0*λ^n + p_1*λ^{n-1} + ... + p_n
    // Our Poly1 representation stores [leading, ..., constant] = [p_0, p_1, ..., p_n]
    // which is [1, p_1, p_2, ..., p_n]
    Ok(Expr::func(
        FuncKind::Poly1,
        vec![Arc::new(Expr::Seq(p))],
    ))
}

fn trace_of_rows(rows: &[Vec<ExprArc>]) -> Result<ExprArc, EvalError> {
    let n = rows.len();
    if n == 0 || rows.iter().any(|r| r.len() != n) {
        return Err(EvalError::TypeError("non-square matrix for trace"));
    }
    let terms: Vec<_> = (0..n).map(|i| Arc::clone(&rows[i][i])).collect();
    Ok(Expr::add(terms))
}

/// Multiply two matrices given as row-major Vec<Vec<ExprArc>>.
fn mul_row_matrices(a: &[Vec<ExprArc>], b: &[Vec<ExprArc>]) -> Result<Vec<Vec<ExprArc>>, EvalError> {
    let nrows = a.len();
    let ncols = b[0].len();
    let inner = b.len();
    if a[0].len() != inner {
        return Err(EvalError::TypeError("incompatible matrix dimensions"));
    }
    let mut out = Vec::with_capacity(nrows);
    for i in 0..nrows {
        let mut row = Vec::with_capacity(ncols);
        for j in 0..ncols {
            let mut sum = Expr::int(0);
            for k in 0..inner {
                let term = Expr::mul(vec![
                    Arc::clone(&a[i][k]),
                    Arc::clone(&b[k][j]),
                ]);
                sum = Expr::add(vec![sum, term]);
            }
            row.push(sum);
        }
        out.push(row);
    }
    Ok(out)
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
    let det = bareiss_det(rows)?;
    Ok(Expr::func(
        FuncKind::Poly1,
        vec![Arc::new(Expr::Seq(vec![
            Expr::int(1),
            neg_expr(&trace),
            minors,
            neg_expr(&det),
        ]))],
    ))
}

fn det2(a: &ExprArc, b: &ExprArc, c: &ExprArc, d: &ExprArc) -> Result<ExprArc, EvalError> {
    Ok(Expr::add(vec![
        Expr::mul(vec![Arc::clone(a), Arc::clone(d)]),
        Expr::mul(vec![Expr::int(-1), Arc::clone(b), Arc::clone(c)]),
    ]))
}

pub fn eval_rref(args: &[ExprArc], ctx: &Context) -> Result<ExprArc, EvalError> {
    if args.is_empty() {
        return Err(EvalError::TooFewArgs("rref"));
    }
    let (rows, modulus) = parse_rref_args(args)?;
    if let Some(m) = modulus {
        let out = mod_rref(&rows, m)?;
        return Ok(Arc::new(Expr::Matrix(out)));
    }
    let (out, _) = rref(&rows)?;
    let out: Vec<Vec<ExprArc>> = out
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|c| eval(c.as_ref(), ctx))
                .collect::<Result<_, _>>()
        })
        .collect::<Result<_, _>>()?;
    eval(Arc::new(Expr::Matrix(out)).as_ref(), ctx)
}

fn parse_rref_args(args: &[ExprArc]) -> Result<(Vec<Vec<ExprArc>>, Option<i64>), EvalError> {
    if args.len() == 1 {
        if let Ok(rows) = as_matrix(&args[0]) {
            return Ok((rows, None));
        }
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
    Ok((rows, modulus))
}

pub fn eval_linsolve(eqs: &ExprArc, vars: &ExprArc, ctx: &Context) -> Result<ExprArc, EvalError> {
    let equations = match eqs.as_ref() {
        Expr::List(v) | Expr::Seq(v) => v.clone(),
        _ => return Err(EvalError::TypeError("linsolve expects equation list")),
    };
    let var_list: Vec<Ident> = match vars.as_ref() {
        Expr::List(v) | Expr::Seq(v) => v
            .iter()
            .map(|e| match e.as_ref() {
                Expr::Symbol(id) => Ok(id.clone()),
                _ => Err(EvalError::TypeError("variable name expected")),
            })
            .collect::<Result<_, _>>()?,
        _ => return Err(EvalError::TypeError("linsolve expects variable list")),
    };
    let n = var_list.len();
    let mut aug = Vec::with_capacity(equations.len());
    for eq in &equations {
        let diff = equation_to_zero(eq)?;
        let diff = eval(diff.as_ref(), ctx)?;
        let mut row = Vec::with_capacity(n + 1);
        for var in &var_list {
            row.push(linear_coeff(&diff, var, &var_list, ctx)?);
        }
        let mut all_zero = HashMap::new();
        for v in &var_list {
            all_zero.insert(v.clone(), Expr::int(0));
        }
        let constant = eval_subst(&diff, &all_zero, ctx)?;
        row.push(neg_expr(&constant));
        aug.push(row);
    }
    let aug: Vec<Vec<ExprArc>> = aug
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|c| eval(c.as_ref(), ctx))
                .collect::<Result<_, _>>()
        })
        .collect::<Result<_, _>>()?;
    let (rref, pivot_cols) = rref(&aug)?;
    let mut solution = vec![Expr::int(0); n];
    for (row_i, &piv_col) in pivot_cols.iter().enumerate() {
        if piv_col >= n {
            continue;
        }
        let mut rhs = Arc::clone(&rref[row_i][n]);
        for j in (piv_col + 1)..n {
            if !is_zero_expr(&rref[row_i][j]) {
                rhs = Expr::add(vec![
                    rhs,
                    Expr::mul(vec![
                        Expr::int(-1),
                        Arc::clone(&rref[row_i][j]),
                        Arc::clone(&solution[j]),
                    ]),
                ]);
            }
        }
        solution[piv_col] = div_expr(&rhs, &rref[row_i][piv_col])?;
    }
    let out: Vec<ExprArc> = solution.into_iter().collect();
    eval(Arc::new(Expr::List(out)).as_ref(), ctx)
}

fn equation_to_zero(eq: &ExprArc) -> Result<ExprArc, EvalError> {
    match eq.as_ref() {
        Expr::Relation(RelOp::Eq, lhs, rhs) => Ok(Expr::add(vec![
            Arc::clone(lhs),
            Expr::mul(vec![Expr::int(-1), Arc::clone(rhs)]),
        ])),
        _ => Ok(Arc::clone(eq)),
    }
}

fn linear_coeff(
    eq: &ExprArc,
    var: &Ident,
    all_vars: &[Ident],
    ctx: &Context,
) -> Result<ExprArc, EvalError> {
    let mut one = HashMap::new();
    let mut zero = HashMap::new();
    for v in all_vars {
        zero.insert(v.clone(), Expr::int(0));
        one.insert(
            v.clone(),
            if v == var {
                Expr::int(1)
            } else {
                Expr::int(0)
            },
        );
    }
    let at_one = eval_subst(eq, &one, ctx)?;
    let at_zero = eval_subst(eq, &zero, ctx)?;
    Ok(Expr::add(vec![
        at_one,
        Expr::mul(vec![Expr::int(-1), at_zero]),
    ]))
}

fn eval_subst(
    expr: &ExprArc,
    subs: &HashMap<Ident, ExprArc>,
    _ctx: &Context,
) -> Result<ExprArc, EvalError> {
    giac_core::eval_subst_map(expr, subs)
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
        Expr::Matrix(rows) | Expr::GiacMatrix(rows) if rows.len() == 1 => Ok(rows[0].clone()),
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

pub(crate) fn rref(rows: &[Vec<ExprArc>]) -> Result<(Vec<Vec<ExprArc>>, Vec<usize>), EvalError> {
    if rows.is_empty() {
        return Err(EvalError::TypeError("empty matrix"));
    }
    if rows.iter().flatten().all(|c| matches!(c.as_ref(), Expr::Int(_))) {
        return int_rref(rows);
    }
    symbolic_rref(rows)
}

fn int_rref(rows: &[Vec<ExprArc>]) -> Result<(Vec<Vec<ExprArc>>, Vec<usize>), EvalError> {
    let mut m: Vec<Vec<Ratio<BigInt>>> = rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|c| match c.as_ref() {
                    Expr::Int(n) => Ratio::from_integer(n.clone()),
                    Expr::Rat(r) => r.clone(),
                    _ => Ratio::from_integer(BigInt::zero()),
                })
                .collect()
        })
        .collect();
    let nrows = m.len();
    let ncols = m[0].len();
    let mut pivot_cols = Vec::new();
    let mut pivot_row = 0usize;
    for col in 0..ncols {
        if pivot_row >= nrows {
            break;
        }
        let mut sel = None;
        for r in pivot_row..nrows {
            if !m[r][col].is_zero() {
                sel = Some(r);
                break;
            }
        }
        let Some(r) = sel else { continue };
        m.swap(pivot_row, r);
        pivot_cols.push(col);
        let pivot = m[pivot_row][col].clone();
        for c in 0..ncols {
            m[pivot_row][c] /= &pivot;
        }
        for r in 0..nrows {
            if r == pivot_row || m[r][col].is_zero() {
                continue;
            }
            let factor = m[r][col].clone();
            let pivot_row_vals: Vec<_> = m[pivot_row].clone();
            for c in 0..ncols {
                m[r][c] -= &factor * &pivot_row_vals[c];
            }
        }
        pivot_row += 1;
    }
    let out = m
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|v| {
                    if v.denom().is_one() {
                        Expr::int(
                            v.numer()
                                .to_string()
                                .parse::<i64>()
                                .unwrap_or(0),
                        )
                    } else {
                        Expr::rat(
                            v.numer().to_string().parse().unwrap_or(0),
                            v.denom().to_string().parse().unwrap_or(1),
                        )
                    }
                })
                .collect()
        })
        .collect();
    Ok((out, pivot_cols))
}

fn symbolic_rref(rows: &[Vec<ExprArc>]) -> Result<(Vec<Vec<ExprArc>>, Vec<usize>), EvalError> {
    if rows.is_empty() {
        return Err(EvalError::TypeError("empty matrix"));
    }
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
                    m[r][c] = Expr::add(vec![Arc::clone(&m[r][c]), neg_expr(&sub)]);
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

fn neg_expr(e: &ExprArc) -> ExprArc {
    Expr::mul(vec![Expr::int(-1), Arc::clone(e)])
}

fn is_zero_expr(e: &ExprArc) -> bool {
    if matches!(e.as_ref(), Expr::Int(n) if n.is_zero()) {
        return true;
    }
    if matches!(e.as_ref(), Expr::Rat(r) if r.is_zero()) {
        return true;
    }
    if let Ok(ev) = eval(e.as_ref(), &Context::default()) {
        return ev.is_zero();
    }
    false
}

fn is_one_expr(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_one())
        || matches!(e.as_ref(), Expr::Rat(r) if *r == Ratio::from_integer(BigInt::one()))
}

#[allow(dead_code)]
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

pub fn try_to_f64_matrix(m: &ExprArc, ctx: &Context) -> Result<Vec<Vec<f64>>, EvalError> {
    let rows = as_matrix(m)?;
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|c| expr_to_f64(c, ctx))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect()
}

pub fn expr_to_f64(c: &ExprArc, ctx: &Context) -> Result<f64, EvalError> {
    let ev = eval(c.as_ref(), ctx)?;
    match ev.as_ref() {
        Expr::Int(n) => n
            .to_string()
            .parse()
            .map_err(|_| EvalError::TypeError("float conversion")),
        Expr::Rat(r) => {
            let num: f64 = r
                .numer()
                .to_string()
                .parse()
                .map_err(|_| EvalError::TypeError("float conversion"))?;
            let den: f64 = r
                .denom()
                .to_string()
                .parse()
                .map_err(|_| EvalError::TypeError("float conversion"))?;
            Ok(num / den)
        }
        _ => Err(EvalError::TypeError("numeric matrix expected")),
    }
}

pub fn f64_to_expr_numeric(v: f64) -> ExprArc {
    if !v.is_finite() {
        return Expr::int(0);
    }
    let tol = 1e-10_f64 * v.abs().max(1.0);
    if (v - v.round()).abs() < tol {
        return Expr::int(v.round() as i64);
    }
    const DEN: i64 = 1_000_000_000_000;
    let num = (v * DEN as f64).round() as i64;
    if num == 0 {
        return Expr::int(0);
    }
    let (n, d) = giac_core::reduce_rational_pair(num, DEN);
    if d == 1 {
        Expr::int(n)
    } else {
        Expr::rat(n, d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use giac_core::format_expr;
    use giac_core::Context;

    #[test]
    fn det_2x2_numeric() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(2)],
            vec![Expr::int(3), Expr::int(4)],
        ]));
        let d = eval_det(&m, &ctx).unwrap();
        assert_eq!(d, Expr::int(-2));
    }

    #[test]
    fn charpoly_4x4_identity() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(0), Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(1), Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(1), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(0), Expr::int(1)],
        ]));
        let x = Ident::new("x");
        let cp = eval_charpoly(&m, &x, &ctx).unwrap();
        let s = format_expr(cp.as_ref());
        // det(λI - I) = (λ-1)^4 => expanded: x^4-4*x^3+6*x^2-4*x+1
        assert!(
            s.contains("x^4") && s.contains("-4"),
            "charpoly 4x4 identity got: {s}"
        );
    }

    #[test]
    fn charpoly_4x4_diagonal() {
        let ctx = Context::default();
        let m = Arc::new(Expr::Matrix(vec![
            vec![Expr::int(1), Expr::int(0), Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(2), Expr::int(0), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(3), Expr::int(0)],
            vec![Expr::int(0), Expr::int(0), Expr::int(0), Expr::int(4)],
        ]));
        let x = Ident::new("x");
        let cp = eval_charpoly(&m, &x, &ctx).unwrap();
        let s = format_expr(cp.as_ref());
        // det(λI - diag(1,2,3,4)) = (x-1)(x-2)(x-3)(x-4)
        // = x^4 - 10x^3 + 35x^2 - 50x + 24
        assert!(
            s.contains("x^4") && s.contains("-10"),
            "charpoly 4x4 diagonal got: {s}"
        );
    }
}
