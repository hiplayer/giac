//! MRV 辅助变量 `w`（`_mrv_w`）——级数系数里的 `ln(w)` / `(-ln(w))^k` 语义。
//!
//! **通用表示层 API 规则：** [`algorithm-expr-api.md`](../../../../../.doc/algorithm-expr-api.md)  
//! **本模块专项契约：** [`limit-engine-expr-api.md`](../../../../../.doc/limit-engine-expr-api.md)
//!
//! | 层级 | 内容 |
//! |------|------|
//! | **Stable** | `mrv_*_expr`、`canonical_mrv_coeff`、`decompose_mrv_coeff`、`expr_contains_*` |
//! | **Temporary** | `drift_*`（仅 `canonical_mrv_coeff` 内；退役 Phase 3A） |
//! | **Pipeline private** | `decompose_neg_ln_power`、`is_ln_w` 等分解辅助 |
//!
//! ## 稳定 API（模块外只应直接使用这些）
//!
//! | 前缀 / 类型 | 用途 |
//! |-------------|------|
//! | `mrv_*_expr`, `neg_ln_w_*_expr` | 规范原子；输出形态固定，可 `==` 比较 |
//! | `canonical_mrv_coeff` | **唯一**漂移收敛入口；`ExprArc → ExprArc` |
//! | `decompose_mrv_coeff` → `MrvCoeffParts` | 语义分解；**输入任意，内部先 canonical** |
//! | `decompose_ln_w_coeff` | 仅 `k·ln(w)` 加法部分 |
//! | `expr_contains_*`, `is_expr_*`, `is_neg_w_inv`, `is_neg_ln_first_power` | 谓词（`is_neg_w_inv` = `w^-1`，≠ `(-ln w)^-1`） |
//!
//! ## 临时 API（`drift_` 前缀，私有，禁止模块外调用）
//!
//! 仅服务于 [`canonical_mrv_coeff`]；退役见 follow-up issue Phase 3A。

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-calculus-api-stability.md`.
//!
use std::sync::Arc;

use giac_core::{Expr, ExprArc, FuncKind, Ident};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

pub(crate) const MRV_W: &str = "_mrv_w";


/// **Stable** — 谓词：是否为 `_mrv_w` 符号
pub(crate) fn is_mrv_w_var(id: &Ident) -> bool {
    id.as_str() == MRV_W
}

/// **Stable** — 规范原子 `Symbol(_mrv_w)`
pub(crate) fn mrv_w_expr() -> ExprArc {
    Expr::sym(MRV_W)
}

/// **Stable** — 规范原子 `Ln(_mrv_w)`
pub(crate) fn mrv_ln_w_expr() -> ExprArc {
    Expr::func(FuncKind::Ln, vec![mrv_w_expr()])
}

/// **Stable** — 规范原子 `-ln(w)`
pub(crate) fn neg_ln_w_expr() -> ExprArc {
    Expr::mul(vec![Expr::int(-1), mrv_ln_w_expr()])
}

/// **Stable** — 规范原子 `(-ln(w))^-1`
pub(crate) fn neg_ln_w_inv_expr() -> ExprArc {
    Expr::pow(neg_ln_w_expr(), Expr::int(-1))
}

/// **Stable** — 谓词：表达式为常数 1
pub(crate) fn is_expr_one(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_one())
}

/// **Stable** — 谓词：表达式为常数 0
pub(crate) fn is_expr_zero(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_zero())
}

/// **Stable** — 谓词：Laurent 因子 `w^-1`（≠ `(-ln w)^-1`）
pub(crate) fn is_neg_w_inv(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Pow(base, exp)
            if is_negative_unit_exp(exp)
                && matches!(base.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)) =>
        {
            true
        }
        Expr::Mul(fs)
            if fs.len() == 2
                && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
                && fs.iter().any(|f| {
                    matches!(f.as_ref(), Expr::Pow(b, e)
                        if is_negative_unit_exp(e)
                            && matches!(b.as_ref(), Expr::Symbol(id) if is_mrv_w_var(id)))
                }) =>
        {
            true
        }
        _ => false,
    }
}

// **Pipeline private** — is negative unit exp
fn is_negative_unit_exp(exp: &ExprArc) -> bool {
    match exp.as_ref() {
        Expr::Int(n) => n.is_negative(),
        Expr::Rat(r) => r.is_negative(),
        Expr::Frac(n, d) => {
            matches!(n.as_ref(), Expr::Int(nn) if nn.is_one()) && is_negative_one_like(d)
        }
        _ => false,
    }
}

// **Pipeline private** — is negative one like
fn is_negative_one_like(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n.is_negative())
        || matches!(
            e.as_ref(),
            Expr::Mul(fs)
                if fs.len() == 2
                    && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n == &BigInt::from(1)))
                    && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n.is_negative()))
        )
}

/// **Stable** — 谓词：子树含 `w`
pub(crate) fn expr_contains_w_var(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Symbol(id) => is_mrv_w_var(id),
        Expr::Add(ts) => ts.iter().any(expr_contains_w_var),
        Expr::Mul(fs) => fs.iter().any(expr_contains_w_var),
        Expr::Pow(b, exp) => expr_contains_w_var(b) || expr_contains_w_var(exp),
        Expr::Frac(n, d) => expr_contains_w_var(n) || expr_contains_w_var(d),
        Expr::Func(_, args) => args.iter().any(expr_contains_w_var),
        _ => false,
    }
}

/// **Stable** — 谓词：子树含 `ln(w)`
pub(crate) fn expr_contains_ln_w(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && expr_contains_w_var(&args[0]) => true,
        Expr::Add(ts) => ts.iter().any(expr_contains_ln_w),
        Expr::Mul(fs) => fs.iter().any(expr_contains_ln_w),
        Expr::Pow(b, exp) => expr_contains_ln_w(b) || expr_contains_ln_w(exp),
        Expr::Frac(n, d) => expr_contains_ln_w(n) || expr_contains_ln_w(d),
        Expr::Func(_, args) => args.iter().any(expr_contains_ln_w),
        _ => false,
    }
}

/// **Stable** — 谓词：规范后等于 `(-ln(w))^1`
pub(crate) fn is_neg_ln_first_power(e: &ExprArc) -> bool {
    let e = canonical_mrv_coeff(e);
    if e == neg_ln_w_expr() {
        return true;
    }
    matches!(
        e.as_ref(),
        Expr::Pow(b, exp)
            if canonical_mrv_coeff(b) == neg_ln_w_expr()
                && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(1))
    )
}


/// **Stable** — MRV 系数规范化；唯一漂移收敛入口
pub(crate) fn canonical_mrv_coeff(expr: &ExprArc) -> ExprArc {
    let folded = match expr.as_ref() {
        Expr::Add(ts) => Expr::add(ts.iter().map(canonical_mrv_coeff).collect()),
        Expr::Mul(fs) => Expr::mul(fs.iter().map(canonical_mrv_coeff).collect()),
        Expr::Pow(b, e) => Expr::pow(canonical_mrv_coeff(b), canonical_mrv_coeff(e)),
        Expr::Frac(n, d) => Arc::new(Expr::Frac(
            canonical_mrv_coeff(n),
            canonical_mrv_coeff(d),
        )),
        Expr::Func(k, args) => Expr::func(*k, args.iter().map(canonical_mrv_coeff).collect()),
        _ => Arc::clone(expr),
    };
    drift_fold_ln_atoms(&folded)
}

/// `decompose_mrv_coeff` 的返回值；消费者用方法判断，勿散落 `is_neg_ln_w_*`。
#[derive(Clone, Debug)]
pub(crate) struct MrvCoeffParts {
    /// 加法意义下 `k·ln(w)` 的 `k`（padd 时与 `ln(w)` 幂相消）。
    pub ln_w_pow: i32,
    /// 乘法意义下 `(-ln(w))^neg_ln_pow`；负指数表示含 `(-ln(w))^-1` 因子。
    pub neg_ln_pow: i32,
    pub rest: ExprArc,
}

impl MrvCoeffParts {
    /// **Stable** — 系数仍含未处理 `ln(w)` / `(-ln(w))^k`
    pub fn pending_for_series(&self) -> bool {
        self.ln_w_pow != 0 || self.neg_ln_pow != 0
    }

    /// **Stable** — 含 `(-ln(w))^-1`（需 peel + padd）
    pub fn has_neg_ln_inv(&self) -> bool {
        self.neg_ln_pow < 0
    }
}

/// **Stable** — 语义分解 `(-ln(w))^k · ln(w) 加法部分 · rest`
pub(crate) fn decompose_mrv_coeff(e: &ExprArc) -> MrvCoeffParts {
    let e = canonical_mrv_coeff(e);
    let (neg_ln_pow, rest) = decompose_neg_ln_power(&e);
    let (ln_w_pow, rest) = decompose_ln_w_coeff(&rest);
    MrvCoeffParts {
        ln_w_pow,
        neg_ln_pow,
        rest,
    }
}

/// **Stable** — 仅 `k·ln(w)` 加法基分解
pub(crate) fn decompose_ln_w_coeff(e: &ExprArc) -> (i32, ExprArc) {
    match e.as_ref() {
        Expr::Add(ts) => {
            let mut k = 0i32;
            let mut rest = Vec::new();
            for t in ts {
                if let Some(lk) = ln_w_mul_coeff(t) {
                    k = k.saturating_add(lk);
                    continue;
                }
                let (tk, tr) = decompose_ln_w_coeff(t);
                k = k.saturating_add(tk);
                if !is_expr_one(&tr) && !is_expr_zero(&tr) {
                    rest.push(tr);
                }
            }
            (k, rest_expr(rest))
        }
        Expr::Mul(fs) => {
            let mut k = 0i32;
            let mut rest = Vec::new();
            for f in fs {
                if let Some(lk) = ln_w_mul_coeff(f) {
                    k = k.saturating_add(lk);
                    continue;
                }
                let (fk, fr) = decompose_ln_w_coeff(f);
                k = k.saturating_add(fk);
                if !is_expr_one(&fr) && !is_expr_zero(&fr) {
                    rest.push(fr);
                }
            }
            (k, rest_expr(rest))
        }
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && expr_contains_w_var(&args[0]) => {
            (1, Expr::int(1))
        }
        Expr::Int(n) if n.is_zero() => (0, Expr::int(0)),
        Expr::Int(n) if n.is_one() => (0, Expr::int(1)),
        Expr::Int(n) if n == &-BigInt::from(1) => (0, Expr::int(-1)),
        Expr::Frac(n, d) if is_expr_one(d) => decompose_ln_w_coeff(n),
        Expr::Frac(n, d) => {
            let (kn, rn) = decompose_ln_w_coeff(n);
            let (kd, rd) = decompose_ln_w_coeff(d);
            if kd == 0 && is_expr_one(&rd) {
                (kn, rn)
            } else if kn == 0 && is_expr_one(&rn) {
                (-kd, rd)
            } else {
                (0, Arc::clone(e))
            }
        }
        other => (0, Arc::new(other.clone())),
    }
}

// **Pipeline private** — decompose neg ln power
fn decompose_neg_ln_power(e: &ExprArc) -> (i32, ExprArc) {
    if is_neg_ln_atom(e) {
        return (1, Expr::int(1));
    }
    if e == &neg_ln_w_inv_expr() {
        return (-1, Expr::int(1));
    }
    match e.as_ref() {
        Expr::Pow(base, exp) if is_neg_ln_atom(base) => {
            if let Expr::Int(n) = exp.as_ref() {
                if let Ok(p) = giac_core::bigint_to_i64(n) {
                    return (i32::try_from(p).unwrap_or(0), Expr::int(1));
                }
            }
            (0, Arc::clone(e))
        }
        Expr::Mul(fs) => {
            let mut k = 0i32;
            let mut rest = Vec::new();
            for f in fs {
                if is_neg_ln_atom(f) {
                    k = k.saturating_add(1);
                    continue;
                }
                if f == &neg_ln_w_inv_expr() {
                    k = k.saturating_sub(1);
                    continue;
                }
                if let Expr::Pow(base, exp) = f.as_ref() {
                    if is_neg_ln_atom(base) {
                        if let Expr::Int(n) = exp.as_ref() {
                            if let Ok(p) = giac_core::bigint_to_i64(n) {
                                k = k.saturating_add(i32::try_from(p).unwrap_or(0));
                                continue;
                            }
                        }
                    }
                }
                rest.push(Arc::clone(f));
            }
            (k, rest_expr(rest))
        }
        Expr::Frac(n, d) => {
            let (kn, rn) = decompose_neg_ln_power(n);
            let (kd, rd) = decompose_neg_ln_power(d);
            (kn.saturating_sub(kd), Arc::new(Expr::Frac(rn, rd)))
        }
        _ => (0, Arc::clone(e)),
    }
}

// **Pipeline private** — is neg ln atom
fn is_neg_ln_atom(e: &ExprArc) -> bool {
    e == &neg_ln_w_expr()
}

// **Pipeline private** — is ln w
fn is_ln_w(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Ln, args) if args.len() == 1 && expr_contains_w_var(&args[0])
    )
}

// **Pipeline private** — ln w mul coeff
fn ln_w_mul_coeff(f: &ExprArc) -> Option<i32> {
    if is_ln_w(f) {
        return Some(1);
    }
    if is_neg_ln_atom(f) {
        return Some(-1);
    }
    match f.as_ref() {
        Expr::Mul(fs) => {
            let mut n = 1i32;
            let mut has_ln = false;
            for x in fs {
                if is_ln_w(x) {
                    has_ln = true;
                    continue;
                }
                if is_neg_ln_atom(x) {
                    has_ln = true;
                    n = n.saturating_mul(-1);
                    continue;
                }
                if let Expr::Int(i) = x.as_ref() {
                    if let Ok(v) = giac_core::bigint_to_i64(i) {
                        n = n.saturating_mul(i32::try_from(v).unwrap_or(0));
                        continue;
                    }
                }
                return None;
            }
            if has_ln {
                Some(n)
            } else {
                None
            }
        }
        _ => None,
    }
}

// **Pipeline private** — rest expr
fn rest_expr(mut parts: Vec<ExprArc>) -> ExprArc {
    if parts.is_empty() {
        Expr::int(1)
    } else if parts.len() == 1 {
        parts.pop().unwrap()
    } else {
        Expr::mul(parts)
    }
}


// **Temporary** — `drift_fold_ln_atoms` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A
fn drift_fold_ln_atoms(expr: &ExprArc) -> ExprArc {
    if drift_is_neg_ln_inv_shape(expr) {
        return neg_ln_w_inv_expr();
    }
    if drift_is_neg_ln_shape(expr) {
        return neg_ln_w_expr();
    }
    if matches!(
        expr.as_ref(),
        Expr::Pow(b, exp)
            if drift_is_neg_ln_shape(b)
                && matches!(exp.as_ref(), Expr::Int(n) if n == &BigInt::from(1))
    ) {
        return neg_ln_w_expr();
    }
    expr.clone()
}

// **Temporary** — `drift_is_neg_ln_shape` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A
fn drift_is_neg_ln_shape(e: &ExprArc) -> bool {
    is_neg_ln_atom(e) || drift_is_neg_ln_expanded(e)
}

// **Temporary** — `drift_is_neg_ln_inv_shape` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A
fn drift_is_neg_ln_inv_shape(e: &ExprArc) -> bool {
    e == &neg_ln_w_inv_expr()
        || matches!(
            e.as_ref(),
            Expr::Pow(base, exp)
                if is_neg_ln_atom(base)
                    && matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
        )
        || drift_is_neg_ln_inv_expanded(e)
}

// **Temporary** — `drift_is_neg_ln_expanded` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A
fn drift_is_neg_ln_expanded(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Mul(fs) if fs.len() == 2 => {
            (drift_is_ln_w_symbol(&fs[0]) && drift_is_negative_one_like(&fs[1]))
                || (drift_is_ln_w_symbol(&fs[1]) && drift_is_negative_one_like(&fs[0]))
        }
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) => {
            drift_is_negative_one_like(base)
        }
        _ => false,
    }
}

// **Temporary** — `drift_is_neg_ln_inv_expanded` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A
fn drift_is_neg_ln_inv_expanded(e: &ExprArc) -> bool {
    match e.as_ref() {
        Expr::Pow(base, exp) if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative()) => {
            drift_is_neg_ln_shape(base)
        }
        Expr::Frac(n, d) if is_expr_one(n) && drift_is_neg_ln_shape(d) => true,
        _ => false,
    }
}

// **Temporary** — `drift_is_negative_one_like` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A
fn drift_is_negative_one_like(e: &ExprArc) -> bool {
    matches!(e.as_ref(), Expr::Int(n) if n == &-BigInt::from(1))
        || matches!(
            e.as_ref(),
            Expr::Pow(base, exp)
                if matches!(exp.as_ref(), Expr::Int(n) if n.is_negative())
                    && matches!(base.as_ref(), Expr::Mul(fs)
                        if fs.len() == 2
                            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n == &-BigInt::from(1)))
                            && fs.iter().any(|f| matches!(f.as_ref(), Expr::Int(n) if n == &BigInt::from(1))))
        )
}

// **Temporary** — `drift_is_ln_w_symbol` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A
fn drift_is_ln_w_symbol(e: &ExprArc) -> bool {
    matches!(
        e.as_ref(),
        Expr::Func(FuncKind::Ln, args)
            if args.len() == 1 && matches!(args[0].as_ref(), Expr::Symbol(id) if is_mrv_w_var(id))
    )
}

#[cfg(test)]
mod tests {
    use giac_core::{format_expr, Expr};

    use super::*;

    #[test]
    fn decompose_mrv_coeff_neg_ln_inv_cancels_in_ratio() {
        let num = Expr::mul(vec![neg_ln_w_inv_expr(), Expr::sym("c")]);
        let den = Expr::mul(vec![neg_ln_w_inv_expr(), Expr::sym("d")]);
        let pn = decompose_mrv_coeff(&num);
        let pd = decompose_mrv_coeff(&den);
        assert_eq!(pn.neg_ln_pow, -1);
        assert_eq!(pd.neg_ln_pow, -1);
        assert!(format_expr(pn.rest.as_ref()).contains("c"));
        assert!(format_expr(pd.rest.as_ref()).contains("d"));
    }

    #[test]
    fn canonical_mrv_coeff_preserves_neg_ln_to_first_power() {
        let drift_ln = Expr::pow(
            Expr::mul(vec![
                Expr::pow(Expr::mul(vec![Expr::int(1), Expr::int(-1)]), Expr::int(-1)),
                mrv_ln_w_expr(),
            ]),
            Expr::int(1),
        );
        let c = canonical_mrv_coeff(&drift_ln);
        assert_eq!(c, neg_ln_w_expr());
    }

    #[test]
    fn canonical_mrv_coeff_folds_drift_neg_ln() {
        let drift = Expr::mul(vec![
            Expr::pow(Expr::mul(vec![Expr::int(1), Expr::int(-1)]), Expr::int(-1)),
            mrv_ln_w_expr(),
        ]);
        let c = canonical_mrv_coeff(&drift);
        assert_eq!(c, neg_ln_w_expr());
    }

    #[test]
    fn is_neg_ln_first_power_atom_and_drift() {
        assert!(is_neg_ln_first_power(&neg_ln_w_expr()));
        let drift = Expr::pow(
            Expr::mul(vec![
                Expr::pow(Expr::mul(vec![Expr::int(1), Expr::int(-1)]), Expr::int(-1)),
                mrv_ln_w_expr(),
            ]),
            Expr::int(1),
        );
        assert!(is_neg_ln_first_power(&drift));
        assert!(!is_neg_ln_first_power(&neg_ln_w_inv_expr()));
    }

    #[test]
    fn mrv_coeff_parts_pending_flags() {
        let p = decompose_mrv_coeff(&neg_ln_w_inv_expr());
        assert!(p.has_neg_ln_inv());
        assert!(p.pending_for_series());
        let q = decompose_mrv_coeff(&Expr::sym("a"));
        assert!(!q.pending_for_series());
    }
}
