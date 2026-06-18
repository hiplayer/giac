# GIAC-limit-exp-diff — `exp` 差分抵消的统一算法与收敛路径

**状态:** open  
**类型:** AFK（实现向）/ 部分 HITL（是否保留 `limit_factored_exp_growth` 快路径直至 MRV 收敛）  
**相关:** [GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md)、[GIAC-limit-maxima-upstream-alignment](GIAC-limit-maxima-upstream-alignment.md)、[GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md)  
**代码:** `preprocess::factor_exp_shifted_difference`、`asymptotic::limit_factored_exp_growth_at_infinity`、`remove_lnexp.rs`、`mrv_lead_term.rs`  
**验收用例:** `maxima_rtest::gruntz_exp_times_exp_diff_minus_one`（✅）；`gruntz_ck_int_60_ratio`、`gruntz_exp_nested_diff`（❌）

---

## 问题陈述

gruntz 类极限含 **嵌套 `exp` 差分**，例如：

```text
limit(exp(x)*(exp(1/x-exp(-x))-exp(1/x)), x, +infinity)  →  -1
```

giac-rs 已通过 **preprocess 因子分解 + 专用抵消** 使该用例通过，但实现分散在两层，且与 upstream **`remove_lnexp` + MRV 级数** 尚未收敛为一条管线。

本 issue 记录：**是否存在统一算法**、当前实现处于哪一层、以及如何与 GIAC-216e 合并，避免继续新增 `exp` 形状表。

---

## 设计原则（normative）

以下三条为 `exp` 差分 / gruntz 类极限的 **架构约束**；新代码与 review 均须遵守。

### P1 — 同一理论，两层规则

**`factor_exp_shifted_difference` 是 `remove_lnexp` 在原变量 `x` 层的前置规则（差分共因子），不是另一套理论。**

| | `x` 层（preprocess） | `w` 层（MRV 换元后） |
|--|----------------------|----------------------|
| 函数 | `factor_exp_shifted_difference` | `remove_lnexp` |
| 典型形状 | `exp(a)*(exp(b+c)-exp(b))` | `exp(f) - w^{-1}` |
| 目标形式 | `exp(主导) * (exp(ε) - 1)` | `w^{-1} * (exp(ε) - 1)` |

二者数学上同属 Gruntz/MRV 的「提出主导指数，余下化为 `exp(ε)-1`」。实现上应朝 **共享规则表 / 同一递归遍历** 收敛（见 Phase B），不得在 preprocess 与 `remove_lnexp` 各维护一套互不相关的特例逻辑。

### P2 — 快路径有保质期

**`limit_factored_exp_growth_at_infinity` 是临时快路径。**

- **现状：** MRV + `remove_lnexp` 未收敛前，用于避免嵌套 `exp` 在级数上挂起；仅覆盖窄形状（如 `L = x + rest`、`S = -exp(-x)`）。
- **目标：** Phase C 验收后 **删除** 生产路径中的调用，或 **仅保留** `#[cfg(test)]` 对照（与删除前 MRV 路径输出比对）。
- **禁止：** 以「再补一条形状」方式扩展 `limit_factored_exp_growth_*`；新 gruntz 用例应驱动 MRV/`remove_lnexp`，而非复制该模式。

### P3 — 禁止 per-gruntz 形状表

**不得为每条 gruntz 再写独立 `limit_*_at_infinity` 形状函数。**

允许扩展的位置 **仅**：

1. **`remove_lnexp` 规则表**（`w` 层，`exp` 差分 / `exp(f)-w^{-1}` / `ln` 展开等）
2. **`preprocess` 递归因子分解**（`x` 层，`factor_exp_shifted_difference` 及其泛化）

反模式（review 应拒绝）：

- 新增 `limit_gruntz_ck60_*`、`limit_nested_exp_diff_*` 等按用例命名的极限路由
- 在 `limit_at_plus_infinity` 主流程中插入第三条、第四条 `exp` 专用分支
- 在 `asymptotic.rs` 为单测绿而写死 pattern match，却不归入上述两表之一

**验收挂钩：** `cargo test -p giac-calculus --lib` 全绿且 **无新增** `limit_*_at_infinity` 形状函数（`sparse_series` / `remove_lnexp` / preprocess 扩展除外）。

---

## 统一数学框架（Gruntz / MRV）

核心恒等式（提取公共指数因子）：

```text
exp(A) * (exp(B) - exp(B'))  =  exp(A + B') * (exp(B - B') - 1)
```

当 `B - B' → 0`（在 `x → +∞` 换元后）时，使用一阶近似：

```text
exp(ε) - 1  ~  ε     （ε → 0）
```

**关键：** 不能对两个因子分别取极限（`lim(exp(L)) * lim(exp(S)-1)` 会错），必须在代数上先合并再取极限。  
gruntz `-1` 例：因子分解后 `exp(x+1/x) * (exp(-exp(-x))-1)`，用 `exp(S)-1 ~ S` 得 `exp(x+1/x) * (-exp(-x))`，再抵消 `exp(x)·exp(-x)` 得 `-exp(1/x) → -1`。

在 MRV 换元后，upstream 用同一思想的另一写法（`series.cc` `remove_lnexp`）：

```text
exp(f) - w^{-1}  =  w^{-1} * (exp(f + ln(w)) - 1)
```

其中 `w` 为 MRV 辅助变量（如 `w = exp(-x)`）。`exp(f + ln(w)) - 1` 中指数趋于 0，由 `SparseSeries` / `padd` 取主项。

**结论（数学上）：** 有统一算法——Gruntz/MRV 的「提取主导 `exp`，剩余化为 `exp(ε)-1` 再级数展开」；**代码上** giac-rs 尚未完全统一。

---

## giac-rs 现状：三层分工

| 模块 | 角色 |
|------|------|
| **`exp_diff.rs`** | P1 共享规则：`fold_exp_shifted_difference`、`exp_scale_times_exp_minus_one`、`try_limit_*_preprocessed`（P2 临时） |
| **预处理（`x` 上）** | `preprocess` 调用 `exp_diff::fold_exp_shifted_difference` |
| **MRV 级数（`w` 上）** | `remove_lnexp` 调用 `exp_diff` + `ln_expand` / `exp_series` |

### 已通过用例走的路径

```text
exp(x)*(exp(1/x-exp(-x))-exp(1/x))
  → factor_exp_shifted_difference
  → exp(x+1/x) * (exp(-exp(-x)) - 1)
  → limit_factored_exp_growth（L=x+1/x, S=-exp(-x), 抵消 → -exp(1/x)）
  → -1
```

**未走** `limit_unidirectional_plus_infinity` / `remove_lnexp` 主路径（此前 MRV 级数在该式上会挂起或超时）。

### 与 upstream 对照

| 能力 | upstream (`series.cc`) | giac-rs |
|------|------------------------|---------|
| `exp` 差分因子分解 | 化简 / `limit_symbolic_preprocess` 内嵌 | `factor_exp_shifted_difference`（显式） |
| `exp(ε)-1` @ MRV | `remove_lnexp` + `padd` | `remove_lnexp`（部分）；gruntz 用形状快路径绕过 |
| 取极限 | `mrv_lead_term` → exponent/coeff | `limit_factored_exp_growth` 窄特例 + MRV 回退 |

---

## 目标统一管线

与 [GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md) 及 upstream 对齐后的目标：

```text
limit_preprocess（merge_exp, pow2expln, factor_exp_shifted_difference）
    → MRV 选 ω
    → 换元 w
    → remove_lnexp（exp(f) - w^{-1}  →  w^{-1}(exp(ε)-1)）
    → SparseSeries 展开 ε（ordre 递增）
    → 主项 → limit
```

实现须遵守上文 **[设计原则](#设计原则normative)**（P1–P3）。

## 剩余缺口（2 条 Maxima gruntz）

| 用例 | 期望 | 为何快路径不够 |
|------|------|----------------|
| CK-INT-60 比值 | `1` | 无 `exp(x)` 与 `-exp(-x)` 的简单抵消；需 MRV 比值 + `remove_lnexp` peel |
| gruntz 嵌套 exp 差 | `1` | 差分结构更深，需 216e 级数在 `w=0` 上稳定主项 |

二者应通过 **扩展 `remove_lnexp` + `mrv_series_lead_loop`** 覆盖，而非复制 `limit_factored_exp_growth` 模式。

---

## 分阶段实施

### Phase A — 文档与测试基线（本 issue）

- [x] 记录统一框架与三层分工（本文档）
- [x] 在 `preprocess::tests` 保留 `factor_exp_shifted_difference` 回归
- [x] 在 `exp_diff::tests` 加入 x 层 / w 层恒等式对照

### Phase B — 规则收敛（落实 P1、P3）✅

1. [x] `exp_diff.rs`：`fold_exp_shifted_difference`（x 层）与 `exp_scale_times_exp_minus_one`（共享核心）
2. [x] `remove_lnexp` 调用 `fold_exp_shifted_difference` + `exp_scale_times_exp_minus_one`（w 层）
3. [x] `preprocess::factor_exp_shifted_difference` → `exp_diff` 别名
4. [x] `limit_factored_exp_growth_*` 从 `asymptotic.rs` 移除 → `exp_diff::try_limit_*_preprocessed`（P2 暂留 `exp_diff`，非形状表）

### Phase C — MRV 主路径接管 gruntz（落实 P2，进行中）

1. 按 [GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md) 补强 `series_lead_at_zero`：`remove_lnexp` 后 `SparseSeries` 取主项，ordre 递增有效。
2. `gruntz_ck_int_60_ratio`、`gruntz_exp_nested_diff` un-ignore 并走 MRV 路径。
3. **删除或降级** `limit_factored_exp_growth_at_infinity`（P2）；`gruntz_exp_times_exp_diff_minus_one` 仍绿。

### 验收

```bash
cargo test -p giac-calculus --lib maxima_rtest -- --include-ignored
cargo test -p giac-calculus --lib
```

- Maxima **14/14**（含 2 条仍 ignore 的 gruntz）
- 满足设计原则 P2、P3：无新增 `limit_*_at_infinity` 形状函数（仅 `remove_lnexp` / `sparse_series` / preprocess 扩展）

---

## 参考

- upstream `giac-1.5.0/src/series.cc`：`remove_lnexp`（~2261）、`padd`（~263–328）、`mrv_lead_term`
- giac-rs：`crates/giac-calculus/src/limit_engine/preprocess.rs`、`remove_lnexp.rs`、`asymptotic.rs`
- Gruntz, *On Computing Limits in a Symbolic Manipulation World*（MRV 与 `exp(ε)-1` 标准处理）
