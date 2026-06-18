# GIAC-limit-pipeline — limit 引擎分层管线落地缺口

**状态:** open  
**类型:** AFK（架构 / 实现向）  
**相关:** [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)、[GIAC-limit-exp-difference-unification](GIAC-limit-exp-difference-unification.md)、[GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md)、[GIAC-limit-maxima-upstream-alignment](GIAC-limit-maxima-upstream-alignment.md)  
**代码:** `limit_engine/{preprocess,exp_diff,asymptotic,mrv_lead_term,remove_lnexp,sparse_series}.rs`  
**验收:** `cargo test -p giac-calculus --lib` 全绿；`maxima_rtest` 含 ignore 用例 un-ignore 后仍绿

---

## 问题陈述

limit 引擎在文档与代码注释中已约定 **四层分工**（x 预处理 → 极限入口 → MRV 换元 → w 层级数），但 **层与层之间尚未串成一条可独立工作的管线**：

- 极限入口（`asymptotic`）仍用 `try_limit_*` **跳层** 绕过 MRV；
- w 层级数（`SparseSeries` + `remove_lnexp`）对嵌套 exp 仍 **不稳定**，回退到 lead-only 形状特例；
- x 层 / w 层对同一恒等式 `exp(scale)*(exp(ε)-1)` 的 **fold 触发条件未统一**，曾导致 CK-INT-61 回归。

本 issue 记录：**各层职责、已落地部分、缺口基础设施、分阶段验收**，作为 exp-diff 统一（P1–P3）与 216e 级数收敛的 **总览索引**。

---

## 目标分层（normative）

```text
┌─────────────────────────────────────────────────────────┐
│ ① x 层 preprocess                                        │
│    pow2expln / merge_exp / fold_exp_shifted_difference   │
│    职责：在原变量 x 上，保守 fold exp 差分为              │
│           exp(L) * (exp(ε) - 1)（有公共项 / mul 内嵌）   │
└──────────────────────────┬──────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────┐
│ ② 极限入口 asymptotic::limit_at_plus_infinity            │
│    职责：优先走 MRV 通用路径；fallback 仅作兜底           │
│    终态：无 try_limit_* 形状快路径（P2 退役）           │
└──────────────────────────┬──────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────┐
│ ③ MRV 换元 mrv + mrv_lead_term                           │
│    职责：mrv 集 → 选 ω → rewrite_in_mrv_w → 递归框架     │
└──────────────────────────┬──────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────┐
│ ④ w 层 remove_lnexp + SparseSeries                       │
│    职责：exp(f)-w⁻¹ → w⁻¹(exp(ε)-1)；级数取主项；        │
│           coeff 无 w/ln(w) → limit_from_mrv_lead         │
└─────────────────────────────────────────────────────────┘
```

**「落地」定义：** 任意 gruntz 类极限按 ①→④ 串行完成，**不依赖** ② 的形状表、④ 的 lead-only 回退特例。

---

## 已落地（Phase B 及近期 refactor）

| 项 | 状态 | 说明 |
|----|------|------|
| x 层 detect + emit | ✅ | `exp_diff::fold_exp_shifted_difference` 单路径；`ExpDiffForm { scale_log, epsilon }` |
| 共享 emit 核心 | ✅ | `exp_scale_times_exp_minus_one` 供 x 层与 w 层共用 |
| w 层调用 x 层 fold | ✅ | `remove_lnexp` 入口调 `fold_exp_shifted_difference` |
| MRV 骨架 | ✅ | `mrv_at_plus_infinity`、`rewrite_in_mrv_w`、`limit_unidirectional_plus_infinity` |
| CK-INT-61 | ✅ | x 层 **禁止** 无公共项的 `exp(A)-exp(B)` fold；由 w 层 / MRV 处理 |
| `limit_factored_exp_growth_*` 迁出 | ✅ | 逻辑收敛至 `exp_diff::try_limit_*_preprocessed`（仍为 P2 临时） |

---

## 未落地缺口

### G1 — ② 极限入口仍「跳层」（P2 未退役）

`limit_at_plus_infinity` 当前在 MRV **之前** 依次尝试：

| 函数 | 文件 | 性质 |
|------|------|------|
| `try_limit_exp_over_exp_via_quotient` | `exp_diff.rs` | 写死 `exp/exp` + 商 trick |
| `try_limit_exp_times_exp_minus_one_preprocessed` | `exp_diff.rs` | 含 `is_neg_exp_of_neg_var`、`exp_rest_after_unit_var` 等形状判断 |
| `try_limit_exp_finite_exponent_at_plus_infinity` | `exp_diff.rs` | 较通用，但仍 bypass MRV |

**问题：** 新 gruntz 用例倾向于再加 `try_*`，违反 [exp-diff P3](GIAC-limit-exp-difference-unification.md#p3--禁止-per-gruntz-形状表)。  
**目标：** Phase C 完成后删除或降级为 `#[cfg(test)]` 对照。

---

### G2 — x 层 / w 层 fold 策略未形式化（P1 边界）

同一恒等式，两层 **触发条件不同**：

| 层 | 何时 fold `exp(A)-exp(B)` |
|----|---------------------------|
| **x 层** | 仅 `remove_add_term` 成功（shared-subterm）；如 `exp(x+1)-exp(x)` |
| **w 层** | `exp(f) - w^{-1}`；如 CK-INT-61 的 `(exp(inner)-exp(x))/x` 换元后 |

**CK-INT-61 回归机理（已验证）：** x 层无条件 fold `exp(inner)-exp(x)` → MRV 级数主项发散，系数残留 `_mrv_w` → `NotImplemented("limit")`。恢复 shared-subterm 约束后修复。

**缺口：** 代码中尚无统一的 `ExpDiffPolicy`（或等价文档 + 单点注释），review 时易再次泛化 x 层 fold。

---

### G3 — ④ 级数主项基础设施不完整（核心阻塞）

详见 [GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md)。摘要：

| 缺口 | 现象 |
|------|------|
| `try_order` 被 `MAX_SERIES_ORDER` cap | 升阶重试无实际效果 |
| `(-ln(w))^{-1}` 语义 | 换元后的 `x^{-1}` 未纳入完整 `padd` 消去链 |
| 同阶相消丢下一阶 | `exp(inner')` 与 `-w^{-1}` 系数 1 + (-1) = 0 后 lead 错误 |
| upstream 前置改写未移植 | upscale、ln(exp^k)、ordre 迭代、spdiv 等 |
| lead-only 回退 | `second_term_inner_plus_ln_expr` 等形状特例仍扛 CK-INT-61 |

**结论：** `exp_diff` 结构侧（detect+emit）已基本到位；**缺的是 MRV 换元后的 SparseSeries 引擎**，不是再加 x 层 `try_rewrite_*`。

---

### G4 — 仍 ignore / 非 MRV 主路径的验收用例

| 用例 | 期望 | 当前 |
|------|------|------|
| `gruntz_exp_nested_diff` | `1` | `#[ignore]`；fold 可工作，完整 limit 未走通 MRV |
| `gruntz_ck_int_60_ratio` | `1` | 部分环境走 `try_limit_exp_over_exp_via_quotient` 快路径 |
| Maxima gruntz 14/14 | 全绿 | 2 条仍 ignore |

---

## 模块职责对照表

| 模块 | 设计职责 | 落地状态 |
|------|----------|----------|
| `preprocess.rs` | x 层：pow2expln、merge_exp、fold | ✅ 已接 `exp_diff` |
| `exp_diff.rs` | P1 fold + P2 临时 `try_limit_*` | ⚠️ fold ✅；`try_limit_*` 待退役 |
| `asymptotic.rs` | 极限入口、MRV 调度 | ⚠️ 仍前置 3 条快路径 |
| `mrv.rs` / `mrv_lead_term.rs` | MRV 选元、换元、递归 | ✅ 骨架完整 |
| `remove_lnexp.rs` | w 层 exp/ln 规则 + padd | ⚠️ 部分；含形状特例 |
| `sparse_series.rs` | w=0 级数、主项 | ❌ 复杂式不稳定（216e） |

---

## 与相关 issue 的分工

| Issue | 本 issue 关系 |
|-------|---------------|
| [GIAC-limit-exp-difference-unification](GIAC-limit-exp-difference-unification.md) | P1–P3 设计原则、exp 恒等式、Phase B/C 任务 |
| [GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md) | G3 级数基础设施的技术细节与 upstream 对照 |
| [GIAC-limit-maxima-upstream-alignment](GIAC-limit-maxima-upstream-alignment.md) | Maxima rtest 总体验收 |

**本 issue 不重复** 216e 的算法路线 A/B/C，仅作 **分层落地状态的总索引**。

---

## 分阶段实施

### Phase 1 — 策略文档化（本 issue）

- [x] 记录四层分工与「落地」定义
- [x] 记录 G1–G4 缺口与 CK-INT-61 回归机理
- [ ] 在 `exp_diff.rs` 模块注释链接本 issue + x/w fold 边界一句 normative 说明

### Phase 2 — 级数基础设施（依赖 216e）

1. [ ] `series_lead_at_zero`：`remove_lnexp` 后稳定主项；ordre 递增有效
2. [ ] `(-ln(w))^{-1}` / padd 消去链
3. [ ] 收缩 `mrv_lead_fallback` / `second_term_inner_plus_ln_expr` 特例

### Phase 3 — 极限入口收敛（依赖 Phase 2）

1. [ ] `gruntz_exp_nested_diff`、`gruntz_ck_int_60_ratio` un-ignore，走 MRV 主路径
2. [ ] 删除或 `#[cfg(test)]` 降级 `try_limit_exp_*`（保留输出对照测试）
3. [ ] `maxima_rtest` 14/14 无 ignore

### 验收

```bash
cargo test -p giac-calculus --lib
cargo test -p giac-calculus --lib maxima_rtest -- --include-ignored
```

- 185+ unit tests 全绿
- 无新增 `limit_*_at_infinity` per-gruntz 形状函数（P3）
- gruntz ignore 清零

---

## 反模式（review checklist）

- [ ] 在 `exp_diff.rs` 新增 `try_rewrite_*` / `try_limit_gruntz_*`
- [ ] x 层 fold 无公共项的 `exp(A)-exp(B)`（CK-INT-61 类）
- [ ] 在 `asymptotic.rs` 插入第四条、第五条 exp 专用分支而不走 MRV
- [ ] 为单条 ignore 用例写独立 `limit_ck_int_*` 形状函数

---

## 参考

- 对话结论：CK-INT-61 回归根因 = x 层 fold 越界；HEAD `exp_diff` 二分验证
- upstream `giac-1.5.0/src/series.cc`：`mrv_lead_term`、`remove_lnexp`、`padd`
- giac-rs：`crates/giac-calculus/src/limit_engine/`
