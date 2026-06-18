# exp_diff 表达式 API 规范（giac-calculus / limit_engine）

> **通用规则（三层 API、表示层边界、新增函数清单）见**  
> [algorithm-expr-api.md](algorithm-expr-api.md)。本文是其在 `exp_diff` / `exp(f)-exp(g)` 差分域的**专项实例**。

**范围：** `giac-rs/crates/giac-calculus/src/limit_engine/exp_diff.rs`  
**上游对照：** `giac-2.0.0/src/series.cc`（`limit_symbolic_preprocess`、`remove_lnexp` 前的 x 层折叠）  
**相关 issue：** [GIAC-limit-exp-difference-unification](issues/GIAC-limit-exp-difference-unification.md)、[GIAC-limit-mrv-followup](issues/GIAC-limit-mrv-followup.md)

**与 MRV 层分工：** 本模块处理 **原变量 `x`** 上的 `exp(A)-exp(B)`、`exp(ε)-1` 等；含 `_mrv_w` / `ln(w)` 系数语义见 [limit-engine-expr-api.md](limit-engine-expr-api.md)（`mrv_w`、`remove_lnexp`）。

---

## 1. 类型边界：Expr 树（x 层）

| 层 | 合法类型 | 禁止 |
|----|----------|------|
| `exp` 差分、`exp(ε)-1` 因子、+∞ 增长分类 | `ExprArc` | 对含 `exp`/`ln` 子树 `expr_to_poly` 判 lead |
| 调用 `ratnormal` / `normal` / `simplify_limit_expr` 之后 | 须重新走 `canonical_exp_diff` 或 `decompose_*` | 假定 AST 与化简前 `==` |
| MRV 换元后 `w` 层 | 委托 `mrv_w` / `remove_lnexp` | 在本模块复制 `(-ln(w))^k` 形状表 |

**形式漂移（AST drift）：** preprocess / `pow2expln` / `ratnormal` 会把同一差分写成不同树，例如：

- `exp(f)-1` → `Add(exp(f), -1)` 或 `Mul(exp(f), -1)`（`simplify_limit_expr` 后）
- `exp(A)-exp(B)` → 未折叠的 `Add` 或已折叠的 `exp(s)*(exp(ε)-1)`
- `exp(-x)` ↔ `Pow(exp(x), -1)`（`unwrap_exp_arg` 吸收）

**规则：** 凡判断「是否为 `exp(ε)-1` 因子」「能否分解为 shifted difference」，必须先走 **§3 稳定 API**；禁止在调用方新增第三个 `if looks_like_exp_*` 分支。

---

## 2. 三层函数分类

（分类定义见 [algorithm-expr-api.md §2](algorithm-expr-api.md#2-三层函数分类全-crate-统一)。）

| 层 | 本模块 |
|----|--------|
| **稳定 API** | §3 |
| **临时 API** | 无 `drift_*` / `shim_*`（漂移吸收在 `detect_exp_difference` 内部，未单独导出） |
| **管线私有** | §4 |

---

## 3. `exp_diff.rs` — 稳定 API（权威清单）

### 3.1 规范入口（漂移收敛）

| 函数 | 输入 | 输出 | 契约 |
|------|------|------|------|
| `canonical_exp_diff(expr)` | 任意 `ExprArc`（x 层子树或整式） | `ExprArc` | 自底向上递归后，将可识别的 `exp(A)-exp(B)` 折为 `exp(scale_log)·(exp(ε)-1)`；**不**调用 `ratnormal`；**不**保证整体最简商式 |

**历史名：** `fold_exp_shifted_difference`（已退役，调用方须改用 `canonical_exp_diff`）。

**调用方：** `preprocess::limit_preprocess_mrv`、`limit_preprocess_struct`、`remove_lnexp`（入口折叠）、单测。

### 3.2 规范构造器

| 函数 | 输入 | 输出 | 输出形态保证 |
|------|------|------|----------------|
| `exp_scale_times_exp_minus_one(scale, ε)` | `ExprArc`, `ExprArc` | `ExprArc` | `scale * (exp(ε) - 1)`（`Add` 形 `-1`，非 `Mul(-1, exp(ε))`） |

### 3.3 语义分解

| 函数 | 输入 | 输出 | 契约 |
|------|------|------|------|
| `match_exp_times_exp_minus_one(expr)` | 已折叠或部分折叠的 `ExprArc` | `Option<(scale_log, ε)>` | 识别 `exp(scale_log)·(exp(ε)-1)`；内部 `detect_exp_difference`；**不**保证输入已 `canonical_exp_diff` |
| `exp_minus_one_epsilon(e)` | 单因子 `ExprArc` | `Option<ExprArc>` | 若 `e` 为 `exp(ε)-1`（`Add` 形），返回 `ε` |
| `unwrap_exp_arg(e)` | `ExprArc` | `Option<(negated, arg)>` | `exp(arg)` 或 `-exp(arg)` 的参数 |

**禁止：** 模块外复制 `detect_exp_difference` / `emit_exp_difference` 逻辑。

### 3.4 谓词（只认规范或检测器吸收后的树）

| 函数 | 含义 |
|------|------|
| `is_exp_minus_one_factor(e)` | `e` 为 `exp(ε)-1` 因子（**只认 `Add` 形**，不认 `Mul(-1, exp(ε))`） |
| `vanishes_at_plus_infinity(e, var)` | `e → 0` 当 `var → +∞` |
| `epsilon_vanishes_at_plus_infinity(e, var)` | `ε` 级小量（Gruntz 预处理用） |

### 3.5 +∞ 增长分类（稳定语义 API）

| 函数 | 用途 |
|------|------|
| `classify_exp_at_plus_infinity(e, var)` | 无符号 `exp` 因子在 +∞ 的相对增长 |
| `classify_signed_exp_at_plus_infinity(e, var)` | 带符号分式积 |
| `unwrap_signed_frac(e)` | 提取 `(sign, inner)` 供分类 |

### 3.6 预处理子步骤（跨模块稳定，边界见 §5）

下列函数被 `preprocess` / `asymptotic` 调用；**语义**稳定，但输出形态依赖调用顺序（须与 §5 管线表一致）：

| 函数 | 作用 |
|------|------|
| `simplify_add_sum` | 合并同类 `Add` |
| `simplify_exp_argument_adds` | 合并 `exp` 参数中的加法 |
| `balance_exp_arguments_frac_var(expr, var)` | 平衡分式指数与 `-var` |
| `first_order_exp_vanishing_epsilon(expr, var)` | Gruntz ε 展开（**破坏** MRV peel 所需的 `exp(·)-1` 形） |
| `algebraize_exp_vanishing_products(expr, var)` | `exp(f)*L` 代数化 |
| `rewrite_exp_minus_scale_inv(...)` | 特定倒数改写 |
| `fold_exp_shifted_difference` 的别名 `preprocess::factor_exp_shifted_difference` | 转发 `canonical_exp_diff` |

---

## 4. 管线私有（禁止模块外复制）

| 函数 | 用途 |
|------|------|
| `detect_exp_difference` / `emit_exp_difference` | `canonical_exp_diff` 内部 |
| `detect_exp_difference_mul` / `_add` | 检测子例程 |
| `flatten_mul_factors`, `signed_add_terms`, `peel_unit_negative_factor` | AST 遍历 |
| `canonical_add_term` | 加法项符号规范（私有，非模块外 API） |
| `try_balance_frac_minus_var` | 单次分式平衡尝试 |
| `balance_frac_minus_var`, `distribute_mul_over_add`, … | 内部改写 |

---

## 5. 邻接模块消费契约

| 模块 | 调用 | 顺序 / 注意 |
|------|------|-------------|
| `preprocess` | `canonical_exp_diff` 在 `pow2expln` 前后各一次 | MRV 路径用 `limit_preprocess_mrv`（**无** `first_order_exp_vanishing_epsilon`） |
| `preprocess` | `limit_preprocess_struct` 含 `first_order_exp_vanishing_epsilon` | **勿**在 MRV peel 前用于含 `(-ln(w))^-1` 的式子 |
| `remove_lnexp` | 入口 `canonical_exp_diff` | 之后走 `mrv_w::canonical_mrv_coeff` |
| `mrv_lead_term` | **不**直接依赖本模块 peel 谓词 | peel 前 **勿** `simplify_limit_expr` |
| `asymptotic` | `classify_*`, `simplify_add_sum` | 代数 lead 与 MRV fallback 分叉 |

---

## 6. 反模式

| ❌ | ✅ |
|----|-----|
| `format_expr(e).contains("exp(")` 判差分形态 | `canonical_exp_diff` + `match_exp_times_exp_minus_one` |
| 对 `Mul(-1, exp(ε))` 用 `is_exp_minus_one_factor` | 先 `canonical_exp_diff` 或走 MRV peel 链 |
| 在 `mrv_lead_term` / `sparse_series` 复制 `exp(ε)-1` 检测 | 调用 `is_exp_minus_one_factor` / `exp_minus_one_epsilon` |
| 第三个 per-case `looks_like_exp_*` | 扩展 `detect_exp_difference` 或 owning 算法 |

---

## 7. 单测要求

- 稳定 API：`exp_scale_times_exp_minus_one` 往返 + `canonical_exp_diff` 在漂移输入上与 `emit_exp_difference` 一致。
- 优先 `assert_eq!` / 构造器比较；`format_expr.contains` 仅作辅助断言。
- CK-INT fixture 见 `ck_int_gruntz_fixture` 与 `exp_diff` 模块内测试。

---

## 8. 当前审查结论（2026-06-18，P0 后）

| 状态 | 项 |
|------|-----|
| ✅ | 规范入口命名为 `canonical_exp_diff`；本文档与索引表已补 |
| ✅ | `mrv_w::is_neg_ln_first_power` 上提；`mrv_lead_term` 谓词统一 `expr_contains_*` |
| ⚠️ | `first_order_exp_vanishing_epsilon` 与 MRV peel 链仍冲突（Phase 2） |
| ⚠️ | `is_exp_minus_one_factor` 不认 `Mul(-1, exp(ε))` — 依赖 peel 顺序或 Phase 3A simplify |

---

## 9. 参考

- 实现：[`exp_diff.rs`](../giac-rs/crates/giac-calculus/src/limit_engine/exp_diff.rs)
- MRV 系数域：[limit-engine-expr-api.md](limit-engine-expr-api.md)
- 算法优先：[`.cursor/rules/algorithm-before-patch.mdc`](../.cursor/rules/algorithm-before-patch.mdc)
