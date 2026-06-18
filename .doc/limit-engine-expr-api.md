# limit_engine 表达式 API 规范（giac-calculus）

> **通用规则（三层 API、表示层边界、新增函数清单）见**  
> [algorithm-expr-api.md](algorithm-expr-api.md)。本文是其在 `limit_engine` / MRV 系数域的**专项实例**。

**范围：** `giac-rs/crates/giac-calculus/src/limit_engine/`  
**上游对照：** `giac-2.0.0/src/series.cc`（`mrv_lead_term`、`remove_lnexp`、`series__SPOL1`）  
**相关 issue：** [GIAC-limit-mrv-followup](issues/GIAC-limit-mrv-followup.md) Phase 1–3A

---

## 1. 类型边界：Expr 树，不是多项式

| 层 | 合法类型 | 禁止 |
|----|----------|------|
| MRV 换元、级数系数、`ln(w)` / `(-ln(w))^k` 语义 | `ExprArc`（`giac_core::Expr`） | 把系数先 `expr_to_poly` 再比形状 |
| 有理函数 lead（`+∞` 代数比较、valuation） | `ExprArc` → 局部 `Poly` | 用 `Poly` 表示含 `exp`/`ln(w)` 的系数 |
| 级数对象 | `SparseSeries`（Laurent in `w` + `ExprArc` 系数） | 系数未经 `canonical_mrv_coeff` 直接 `==` 比较 |

**形式漂移（AST drift）：** `ratnormal` / `normal` / preprocess 会把同一数学对象写成不同 AST，例如：

- `-ln(w)` → `((1*-1)^-1)*ln(w)`
- `w^-1` → `w^(1)/(-1)`（指数 `Frac(1, Mul(-1))`）
- `exp(f)-1` → `exp(f)*(-1)`（`simplify_limit_expr` 后）

**规则：** 凡涉及 `ln(w)`、`(-ln(w))^k`、`w^-1` 语义的 **判断、分解、相消**，必须先走 **稳定 API**（见 §2），禁止在调用方散落形状匹配。

---

## 2. 三层函数分类

（分类定义见 [algorithm-expr-api.md §2](algorithm-expr-api.md#2-三层函数分类全-crate-统一)；以下为 limit_engine 落地。）

---

## 3. `mrv_w.rs` — MRV 系数稳定 API（权威清单）

### 3.1 规范原子构造器

| 函数 | 输入 | 输出 | 输出形态保证 |
|------|------|------|----------------|
| `mrv_w_expr()` | — | `ExprArc` | `Symbol("_mrv_w")` |
| `mrv_ln_w_expr()` | — | `ExprArc` | `Ln(_mrv_w)` |
| `neg_ln_w_expr()` | — | `ExprArc` | `Mul(-1, Ln(_mrv_w))`（**唯一** `-ln(w)` 原子） |
| `neg_ln_w_inv_expr()` | — | `ExprArc` | `Pow(neg_ln_w_expr(), -1)`（**唯一** `(-ln(w))^-1` 原子） |

**用途：** 测试断言、构造已知系数；**不要**手写 `Mul(-1, Ln(w))` 代替构造器。

### 3.2 规范化（漂移收敛唯一入口）

| 函数 | 输入 | 输出 | 契约 |
|------|------|------|------|
| `canonical_mrv_coeff(expr)` | 任意 `ExprArc`（级数系数子树） | `ExprArc` | 递归子树后，将可识别的 `ln(w)` 漂移折为 §3.1 原子；**不**调用 `ratnormal`；**不**保证整体化为最简商式 |

**调用方规则：**

- 比较两个系数是否含相同 `(-ln(w))^k` 因子 → 先 `canonical_mrv_coeff`，再 `decompose_mrv_coeff` 或构造器 `==`。
- **禁止**在 `mrv_w` 外新增「识别 `((1*-1)^-1)*ln(w)`」类函数；漂移吸收只许在 `canonical_mrv_coeff` 内扩展 `drift_*`。

### 3.3 语义分解

| 函数 | 输入 | 输出 | 契约 |
|------|------|------|------|
| `decompose_mrv_coeff(e)` | `ExprArc` | `MrvCoeffParts` | 内部先 `canonical_mrv_coeff`；语义 `e ≡ (-ln(w))^neg_ln_pow · (k·ln(w) 加法部分) · rest`（乘法/加法意义见字段注释） |
| `decompose_ln_w_coeff(e)` | `ExprArc` | `(i32, ExprArc)` | `(k, rest)` 表示加法意义 `k·ln(w) + rest`；`Frac(n,1)` 展开到分子；**不**单独抽取 `(-ln(w))^k`（用 `decompose_mrv_coeff`） |

**`MrvCoeffParts` 方法（优先于字段散落判断）：**

| 方法 | 含义 |
|------|------|
| `pending_for_series()` | 系数仍含未处理 `ln(w)` / `(-ln(w))^k` |
| `has_neg_ln_inv()` | 含 `(-ln(w))^-1`（需 peel + padd） |

### 3.4 谓词（只认规范或漂移吸收后的树）

| 函数 | 输入 | 输出 | 说明 |
|------|------|------|------|
| `is_mrv_w_var(id)` | `&Ident` | `bool` | 是否为 `_mrv_w` |
| `is_expr_one` / `is_expr_zero` | `&ExprArc` | `bool` | 常数 0/1 |
| `is_neg_w_inv(e)` | `&ExprArc` | `bool` | **`w^-1` Laurent 因子**，含 `Pow(w, Frac(1,-1))` 等漂移；**不是** `(-ln(w))^-1` |
| `expr_contains_w_var(e)` | `&ExprArc` | `bool` | 子树含 `w` |
| `expr_contains_ln_w(e)` | `&ExprArc` | `bool` | 子树含 `ln(w)` |
| `is_neg_ln_first_power(e)` | `&ExprArc` | `bool` | 规范后等于 `(-ln(w))^1`（peel 分母识别；**不是** `(-ln(w))^-1`） |

### 3.5 临时 API（`mrv_w.rs` 私有，`drift_*`）

| 函数 | 用途 | 退役 |
|------|------|------|
| `drift_fold_ln_atoms` | `canonical_mrv_coeff` 叶节点折叠 | Phase 3A |
| `drift_is_neg_ln_shape` / `drift_is_neg_ln_inv_shape` | 识别漂移 AST | Phase 3A |
| `drift_is_neg_ln_expanded` / `drift_is_neg_ln_inv_expanded` | 展开形识别 | Phase 3A |
| `drift_is_ln_w_symbol` / `drift_is_negative_one_like` | 子形状 | Phase 3A |

**模块内私有、非 drift 但仍属实现细节（禁止外调）：**  
`decompose_neg_ln_power`, `is_neg_ln_atom`, `is_ln_w`, `ln_w_mul_coeff`, `rest_expr`, `is_negative_unit_exp`, `is_negative_one_like`。

---

## 4. 邻接模块稳定 API（消费 `mrv_w` 的契约）

### 4.1 `remove_lnexp.rs`

| 函数 | 输入 | 输出 | 前置 / 后置 |
|------|------|------|-------------|
| `remove_lnexp(expr, ctx)` | `ExprArc`, `Context` | `ExprArc` | 输入可为换元后任意树；输出仍可能含 `ln(w)`，应用 `canonical_mrv_coeff` 再分解 |
| `divide_lead_coeffs(num, den, ctx)` | 两个系数 `ExprArc` | `ExprArc` | **padd 商**；内部 `decompose_mrv_coeff` 同幂 `ln(w)` / `(-ln(w))^k` 相消；fallback `ratnormal` |
| `expr_contains_exp_or_ln(e)` | `ExprArc` | `bool` | 与 MRV 系数规范无关的粗谓词 |

**管线私有（临时改写，非稳定 API）：**  
`try_rewrite_exp_minus_w_inv`, `try_rewrite_w_inv_times_exp_minus_one`, `try_collapse_w_inv_exp_*`, `lead_after_exp_ln_cancel`, `second_term_inner_plus_ln_expr`。

### 4.2 `mrv_lead_term.rs`（MRV lead 管线入口）

| 函数 | 输入 | 输出 |
|------|------|------|
| `mrv_lead_term_plus_infinity(expr, var, ctx)` | 原式或 MRV-preprocess 后式 | `Result<MrvLeadTerm, EvalError>` |
| `limit_from_mrv_lead_term(lead, var, ctx)` | `MrvLeadTerm` | `Result<ExprArc, EvalError>` |
| `limit_unidirectional_plus_infinity` | 同上 | 递归 lead |
| `series_lead_at_zero` | `w=0` 级数 lead | `MrvLeadTerm` |

**管线顺序契约（`mrv_lead_term_plus_infinity`）：**

```text
limit_preprocess_mrv → rewrite_in_mrv_w → canonical_mrv_coeff → peel_neg_ln_w_inv
  →（命中）lead_from_peeled_core → remove_lnexp → …
  →（未命中）simplify_limit_expr → mrv_series_lead_loop
```

**关键：** `peel` 必须在 `simplify_limit_expr` **之前**（否则 `exp(f)-1` 漂成 `exp(f)*(-1)`）。

**管线私有：** `peel_neg_ln_w_inv`, `lead_from_peeled_core`, `mrv_series_lead_loop*`, `rewrite_in_mrv_w`（剥壳分母识别用 `mrv_w::is_neg_ln_first_power`）。

### 4.3 `sparse_series.rs`（级数系数）

| 函数 | 输入 | 输出 |
|------|------|------|
| `series_at_zero` / `series_at_zero_order` | `ExprArc`, `Ident` w, order | `SparseSeries` |
| `series_spdiv_one` | `SparseSeries` | 升阶除法 |

**系数契约：** 遍历系数时用 `decompose_mrv_coeff(c).pending_for_series()` 判断符号 `ln(w)`；**禁止** `is_neg_ln_w_expr` 类散落匹配。

### 4.4 明确不用 `mrv_w` 的模块

| 模块 | 类型 | 说明 |
|------|------|------|
| `mrv.rs` | `ExprArc` + 增长比较 | MRV 集合与 `choose_mrv_w`；`linear_coeff_in_var` 处理 `Frac` 指数 |
| `preprocess.rs` | `ExprArc` | `limit_preprocess_struct` 可破坏 `exp(·)-1`；MRV 门禁勿先调 struct |
| `asymptotic.rs` | `ExprArc` + 局部 `Poly` | 代数 lead / 度数比；CK-61 快路径待 Phase 2 退役 |
| `exp_diff.rs` | `ExprArc` | `exp(f)-1` / `exp(f)*L` 代数；`is_exp_minus_one_factor` 只认 `Add` 形 |

---

## 5. 反模式（审查 checklist）

| ❌ 反模式 | ✅ 应做 |
|----------|---------|
| `if format_expr(e).contains("ln(_mrv_w)")` | `expr_contains_ln_w` / `decompose_mrv_coeff` |
| `expr_to_poly` 判断 `ln(w)` 幂 | `decompose_mrv_coeff` |
| 模块外 `drift_*` 或复制其逻辑 | 扩展 `canonical_mrv_coeff` 或修 `giac-simplify` |
| 对 `ratnormal` 结果直接 `peel` 前不 `canonical` | `canonical_mrv_coeff` 再 `peel_neg_ln_w_inv` |
| 手写 `is_neg_ln_w_inv_expr` 形状表 | `decompose_mrv_coeff` + `has_neg_ln_inv()` |
| 在 `sparse_series` / `mrv_lead_term` 新增第三个 `(-ln(w))` 特例 | 回 `mrv_w` 或 `remove_lnexp` 补一条 padd 规则 |

---

## 6. 新增函数流程

1. **归类：** 稳定 / 临时（`drift_`）/ 管线私有？
2. **类型：** 签名是否全程 `ExprArc`（或 `SparseSeries`）？若需 `Poly`，是否仅限无 `exp`/`ln` 子树？
3. **I/O 契约：** 用表格写清输入形态、输出形态、是否调用 `canonical_mrv_coeff`。
4. **漂移：** 新识别的漂移形 → 只加在 `drift_*` + 单测 `canonical_mrv_coeff_*`；语义操作只加在 `decompose_*` / `divide_lead_coeffs`。
5. **单测：** 稳定 API 必须有构造器往返 + 漂移折叠用例；禁止仅 `format_expr` 子串断言。
6. **文档：** 更新本文 §3/§4 表格；临时函数注明退役条件。

---

## 7. 当前审查结论（2026-06-18）

| 状态 | 项 |
|------|-----|
| ✅ 符合 | `mrv_w` 稳定/临时分层；`drift_*` 全私有；外模块已迁 `decompose_mrv_coeff` |
| ✅ 符合 | `divide_lead_coeffs` 统一 padd；`is_neg_w_inv` 吸收 `w^-1` 指数漂移 |
| ✅ 符合 | `is_neg_ln_first_power` 已上提到 `mrv_w`；`mrv_lead_term` 谓词统一 `expr_contains_*` |
| ⚠️ Phase 3A | 全部 `drift_*` 删除后，`canonical_mrv_coeff` 改委托 `giac-simplify` |
| ⚠️ Phase 2 | `limit_preprocess_struct` 与 MRV peel 链不兼容；门禁用 `ck_int_61()` 直进 MRV |

---

## 8. 参考

- 模块实现：[`mrv_w.rs`](../giac-rs/crates/giac-calculus/src/limit_engine/mrv_w.rs)
- 算法优先： [`.cursor/rules/algorithm-before-patch.mdc`](../.cursor/rules/algorithm-before-patch.mdc)
- Phase 计划： [`GIAC-limit-mrv-followup.md`](issues/GIAC-limit-mrv-followup.md)
