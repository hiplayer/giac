# GIAC-limit-mrv-followup — MRV 主路径收敛与 limit_engine 绕行清理

**状态:** open  
**类型:** AFK（实现向）  
**优先级:** P0（Phase 1–2）/ P1（Phase 3）  
**相关:** [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)（LIM-G2/G3/G4/G6）、[GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md)、[GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md)、[GIAC-limit-exp-difference-unification](GIAC-limit-exp-difference-unification.md)  
**上游对照:** `giac-2.0.0/src/series.cc`（`mrv_lead_term`、`remove_lnexp`、`series__SPOL1`，~2853–2975）  
**代码:** `limit_engine/{mrv_lead_term,remove_lnexp,sparse_series,asymptotic,exp_diff,simplify_util}.rs`；`giac-simplify/expand.rs`  
**验收:** `cargo nextest run -p giac-calculus` 192/192；CK-58/60/61 三角不回归

---

## 问题陈述

2026-06-18 一轮清理后，`limit_engine` 已消除 eval/ratnormal 双 fallback、商式归一化宽/窄分离、
`gruntz_exp_nested_diff`（LIM-G1）等缺口。但 **MRV 主路径仍未独立处理 CK-INT-61**，
仍依赖 `limit_exp_times_frac_quotient_at_plus_infinity` 形状特例；另有若干 **功能绕行**
待在 owning 模块补齐（`algorithm-before-patch`）。

本 issue 记录四项后续工作的 **分阶段计划、依赖、验收门禁**，作为 AFK 实现索引。

---

## 背景：目标管线（不变）

```text
preprocess → MRV → remove_lnexp → SparseSeries lead → limit
```

**禁止**新增 `try_limit_*` 或 per-CK 形状表；特例退役条件见 Phase 2。

**已完成的近期基础**（本 issue 的前置，勿回退）：

- `sparse_series` 升阶不再硬卡 `MAX_SERIES_ORDER=10`（MRV 可用 `MAX_SERIES_EXPANSION_ORDER`）
- `simplify_util::simplify_limit_expr`（`ratnormal→normal`）替代部分 `eval` 收尾
- `canonicalize_limit_entry` / `normalize_inverse_sums` 宽/窄分离（CK-58 parsed 与 CK-61 共存）

---

## 依赖总览

```text
Phase 1 (LIM-G3)  padd / (-ln w)⁻¹
        ↓
Phase 2 (LIM-G2)  MRV lead 收敛 CK-61 → 退役 limit_exp_times_frac_quotient
        ↓
Phase 3C            limit_from_mrv_lead_term 去 eval(coeff)

Phase 3A (并行)     giac-simplify 受限 expand
Phase 3B (并行)     ratnormal unwrap_or 显式化
```

---

## Phase 1 — LIM-G3：`(-ln(w))⁻¹` / `x⁻¹` padd 消去

### 问题

MRV 换元后（`w = exp(-x)`，`x → -ln(w)`）表达式常含：

```text
(exp(inner') - w^-1) * (-ln(w))^-1
```

`(-ln(w))^-1` 即 `x^-1`，不是普通 `w^k` Laurent 项。当前：

- `lead_coeff_ready` 拒绝含 `ln(w)` / `(-ln(w))^-1` 的系数
- `divide_lead_coeffs` / `remove_lnexp` padd 链不完整
- `exp(inner')` 与 `-w^-1` 同阶相消后，次阶项未被稳定提取

### 任务

| ID | 工作 | 文件 | 验收 |
|----|------|------|------|
| 1.1 | 审计 `peel_neg_ln_w_inv` / `combine_lead_with_ln_inv` / `lead_from_peeled_core` 与 upstream 差异 | `mrv_lead_term.rs` | 文档化 CK-61 换元后形态 |
| 1.2 | 扩展 `divide_lead_coeffs`：`ln(w)^k` 与 `(-ln(w))^-1` 的 padd 规则 | `remove_lnexp.rs`, `mrv_w.rs` | 单测：`c·(-ln(w))⁻¹ / d·(-ln(w))⁻¹ → c/d` |
| 1.3 | `series_at_zero` 对 `_mrv_w`：`Add` 先 `remove_lnexp` 再展开，补同阶 cancel 后次阶路径 | `sparse_series.rs` | 单测：同阶 `w^-1` 相消后露出次阶 |
| 1.4 | `lead_coeff_ready` / `is_series_coeff_undef`：区分「待处理 `(-ln(w))⁻¹`」与「已 padd 的 `ln(w)` 幂」 | `mrv_lead_term.rs` | 与 1.2 语义一致 |

### Phase 1 API 规范

通用规则：[algorithm-expr-api.md](../algorithm-expr-api.md)（Cursor：`algorithm-expr-api.mdc`、`algorithm-before-patch.mdc`）。  
limit_engine 专项契约：[limit-engine-expr-api.md](../limit-engine-expr-api.md)（仅 `.doc/`，无独立 Cursor 规则）。

### 门禁测试

- `mrv_lead_ck_int_61`：对 `ck_int_61()` 经 `mrv_lead_term_plus_infinity` 内部 MRV preprocess 必须 `Ok`，`limit_from_mrv_lead_term → -exp(2)`（勿先 `limit_preprocess_struct`，其 `algebraize` 会破坏 `exp(·)-1` 形态）
- 不破坏：`gruntz_exp_nested_diff`、`engine_ck_int_61`、`limit_ck_int_58_parsed`

### 估时

1–2 轮迭代（核心难点 1.2 / 1.3）

---

## Phase 2 — LIM-G2：MRV 主路径收敛 CK-61，退役特例

### 问题

`limit_exp_times_frac_quotient_at_plus_infinity`（`asymptotic.rs`）通过
`dominant_ln_exponent` + `float_to_expr` 处理 CK-61，属于 **形状快路径**，
与 [GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md) Phase C 退役目标冲突。

### 前置

**Phase 1 完成**（G3 是 G2 硬依赖）

### 任务

| ID | 工作 | 文件 | 验收 |
|----|------|------|------|
| 2.1 | 移植 upstream 前置改写：`rewrite_ln_exp_in_f`、`ln(exp(g)^k)` 重写（查漏补缺） | `mrv_lead_term.rs` | CK-61 换元后 `series__SPOL1` 输入对齐 upstream |
| 2.2 | 验证 `mrv_series_lead_loop_inner` 的 `spdiv` + ordre×1.5 升阶（order cap 已修） | `mrv_lead_term.rs`, `sparse_series.rs` | 扩展 `series_at_zero_order_escalates_beyond_default_cap` 为 MRV `w` 场景 |
| 2.3 | 硬化 `mrv_lead_ck_int_61`：`expect("mrv lead")` + `assert -exp(2)` | `mrv_lead_term.rs` | 测试不可静默跳过 |
| 2.4 | **删除** `limit_exp_times_frac_quotient_at_plus_infinity` 及 `limit_preprocessed` 链中的 `.or_else` 调用 | `asymptotic.rs` | 见退役条件 |
| 2.5 | 更新 [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md) LIM-G2；勾选 pipeline Phase C | `.doc/issues/` | 文档与代码一致 |

### 退役条件（全部满足才删 2.4）

```text
mrv_lead_term_plus_infinity(preprocess(ck_int_61))  → Ok
limit_from_mrv_lead_term(lead)                      → -exp(2)
limit_at_plus_infinity(preprocess(ck_int_61))         → -exp(2)   // 无特例介入
```

### 估时

1–2 轮（依赖 Phase 1 质量）

---

## Phase 3A — `giac-simplify` 受限 `expand`（替代 `distribute_mul_over_add`）

### 问题

`exp_diff::distribute_mul_over_add` 是本地 distribute 实现；全局 `giac_simplify::expand`
曾破坏 CK-58（对含 `exp` 子树过度展开）。按 `algorithm-before-patch`，能力应归 **owning 模块**
`giac-simplify`，`limit_engine` 只调用。

### 任务

| ID | 工作 | 文件 | 验收 |
|----|------|------|------|
| 3A.1 | 新增 `expand_polynomial` 或 `expand(..., ExpandPolicy::NoExpDistribute)` | `giac-simplify/expand.rs` | 单测：`(x+1)*((x-1)^-1)` OK；含 `exp` 子树不 distribute |
| 3A.2 | `exp_diff::simplify_balanced_frac` 分子改用受限 expand | `exp_diff.rs` | `fold_ck_int_61_shape`、`limit_ck_int_58_parsed` 仍绿 |
| 3A.3 | 删除或薄包装 `distribute_mul_over_add`（可保留为对照测试） | `exp_diff.rs` | 无行为回归 |

### 门禁

`cargo nextest run -p giac-calculus` 全绿，重点 CK-58/60/61

### 估时

0.5–1 轮（先定 `ExpandPolicy` API）

---

## Phase 3B — `ratnormal(...).unwrap_or` 静默降级清理

### 问题

约 15 处 `ratnormal` 失败时原样继续；含 `exp` 的式子有理化失败后，
MRV/级数建立在 **未归一** 表达式上，属隐式绕行。

### 策略

集中到 `simplify_limit_expr`，按场景分三类：

| 类 | 场景 | 处理 |
|----|------|------|
| A | 结果规范化 | `simplify_limit_expr`（`asymptotic::normalize_limit_result` 已用） |
| B | MRV 预处理必经 | `ratnormal` 失败 → `normal` → 仍失败则 `NotImplemented("series")` |
| C | 级数系数 normalize | `simplify_limit_expr`（`normalize_series_coeff` 已部分完成） |

### 替换优先级

| 优先级 | 文件 | 约处数 |
|--------|------|--------|
| P0 | `mrv_lead_term.rs` | 7 |
| P1 | `sparse_series.rs` | 2 |
| P1 | `remove_lnexp.rs` | 1 |
| P2 | `asymptotic.rs` | 4 |
| P2 | `preprocess.rs` | 1 |

### 验收

- 单测：`ratnormal` 对含 `exp` 式子失败时，路径行为明确（非静默原样）
- 192/192 不回归

### 估时

1 轮（机械替换 + 补测试）

---

## Phase 3C — `limit_from_mrv_lead_term` 去 `eval(coeff)`

### 问题

`lead.exponent == 0` 时对 coeff 调用 `eval`，损失符号信息；应等 MRV 系数纯符号化后去除。

### 前置

Phase 2（coeff 须为纯常数或已知符号形式，如 `-exp(2)`）

### 任务

| ID | 工作 | 验收 |
|----|------|------|
| 3C.1 | `rewrite_ln_w` 后断言 coeff 不含 `w`、`ln(w)`、原 `var` | `depends_on_var(coeff, var) == false` |
| 3C.2 | `simplify_limit_expr` + 符号比较替代 `eval(coeff)` | `limit_from_mrv_lead_term({coeff:-exp(2)}) → -exp(2)` |
| 3C.3 | `sign_infinity` 分支：coeff 符号用 `is_negative_const_expr` | MRV 负无穷路径 |

### 估时

0.5 轮（宜与 Phase 2 同步）

---

## 里程碑

| 里程碑 | 内容 | 目标 |
|--------|------|------|
| **M1** | Phase 1 完成，`mrv_lead_ck_int_61` 硬化通过 | 第 1 周 |
| **M2** | Phase 2 完成，删除 `limit_exp_times_frac_quotient` | 第 1–2 周 |
| **M3** | Phase 3B + 3C 完成 | 第 2 周 |
| **M4** | Phase 3A 受限 expand 落地 | 第 2–3 周（可与 M3 并行） |

---

## 风险与约束

1. **禁止**在 G2 未收敛前删除 `limit_exp_times_frac_quotient`
2. **禁止**在 `limit_engine` 复制 `giac-simplify` 的 expand/normal；3A 必须在 `giac-simplify` 补
3. **禁止**合并商式归一化宽/窄规则（`canonicalize_limit_entry` vs `normalize_inverse_sums`）
4. 每阶段结束：`cargo nextest run -p giac-calculus` + CK-58/60/61 三角验证
5. 新算法债须更新 [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)，不得只加 `#[ignore]` 无条目

---

## 验证命令

```bash
cd giac-rs
cargo nextest run -p giac-calculus
cargo nextest run -p giac-calculus -E \
  'test(engine_ck_int_61) or test(limit_ck_int_61_parsed) or test(limit_ck_int_58_parsed) or test(mrv_lead_ck_int_61) or test(gruntz_exp_nested_diff)'
```

---

## 参考

- `.cursor/rules/algorithm-before-patch.mdc`
- `.doc/module-division.md` — `series.cc` 模块映射
- `giac-2.0.0/check/testintegrate` L61–62（CK-INT-60/61）
- `giac-2.0.0/check/testlimit` L10,L21
