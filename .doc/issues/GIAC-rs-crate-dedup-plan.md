# GIAC-rs — crate 去重与临时代码退役计划

**状态:** open  
**类型:** 执行计划 / AFK 竖切  
**来源:** `giac-rs/crates` 模块审计（2026-06-22）  
**父索引:** [GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md)（本计划补 **函数合并 + crate 内整理**；不重复 0A–4C 已跟踪项）  
**规范:** [algorithm-expr-api.md](../algorithm-expr-api.md)、[algorithm-before-patch.mdc](../../.cursor/rules/algorithm-before-patch.mdc)  
**快照:** 2026-06-22

---

## 1. 问题陈述

`giac-rs/crates` 在快速移植阶段积累了四类可清理债务：

| 类别 | 典型位置 | 风险 |
|------|----------|------|
| **Temporary**（shim / drift） | `equiv::canonical_radical`、`ratnormal_algext`、`mrv_w::drift_*` | 测试通过但语义漂移；调用方复制 |
| **Partial 启发式** | `integrate` / `limit_engine` 大量 `try_*` | 与主路径双轨；按用例膨胀 |
| **拷贝工具 fn** | `ident_from_expr`×6、`is_var`×10+ | 修一处漏多处 |
| **平行算法路径** | 二次分解三处、solve 未接 `poly_algext_roots` | 行为不一致、重复维护 |

**原则（与 tech-debt 一致）：**

1. **行为不变优先** — Phase D0 只做搬运/合并，不改数学语义。  
2. **进管线，不旁路** — 退役 Temporary/Partial 时扩展主路径，禁止在调用方加第三个形状特例。  
3. **每 PR 可独立合入** — 竖切 + `cargo test-timeout` 全绿。  
4. **不重复开 parent issue** — 下列子项若已在 [GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md)、[GIAC-poly-p3-6](GIAC-poly-p3-6-quartic-roots-gaps.md)、[GIAC-limit-mrv-followup](GIAC-limit-mrv-followup.md) 跟踪，本计划仅作 **执行顺序与落点** 索引。

---

## 2. 临时代码清单（退役目标）

### 2.1 明确 **Temporary**（须写退役条件，删时连带单测）

| ID | 函数 / 模块 | 文件 | 退役条件 | 本计划子项 |
|----|-------------|------|----------|------------|
| T-01 | `ratnormal_algext` | `giac-simplify/ratnormal.rs` | [GIAC-algext-adoption](GIAC-algext-adoption.md) A-03 `ext_reduce` | → tech-debt **2B** |
| T-02 | `canonical_radical`, `inv_sqrt_to_mul` | `giac-simplify/equiv.rs` | 迁入 `normal` 或稳定 `canonical_radical` pub API | → **D4-2** |
| T-03 | `try_factor_quadratic_rootof`, `try_factor_quadratic_sqrt` | `giac-simplify/factor.rs` | 二次路径统一后经 `rootof` / giac-poly 共享核 | → **D1** |
| T-04 | `drift_*`（7 个） | `giac-calculus/limit_engine/mrv_w.rs` | [GIAC-limit-mrv-followup](GIAC-limit-mrv-followup.md) Phase 3A | → tech-debt **2A** |
| T-05 | `split_depressed_quartic` | `giac-core/algebra/poly_roots.rs` | 测试迁到 Euler 路径后删除 legacy Ferrari | → **D4-3** |
| T-06 | 空 `stubs.rs` | `giac-calculus`, `giac-ode` | 无 fn；可删模块声明或保留一行 doc | → **D4-4** |

### 2.2 **Partial** stub（范围窄，按域 issue 退役）

| ID | 函数 | 文件 | 退役方向 | 本计划子项 |
|----|------|------|----------|------------|
| P-01 | `eval_fsolve` | `giac-solve/fsolve.rs` | 完整 Newton + 区间 | 独立 backlog（本计划不展开） |
| P-02 | `biquadratic_rootof_roots`（一般四次） | `giac-solve/rootof.rs` | `poly_algext_roots` + Expr 包装 | → **D3** |
| P-03 | `try_integrate_*`（~25） | `integrate.rs`, `integrate_heuristics.rs` | Risch / partfrac 主路径 | → **D2**；长期见 calculus-api-stability §4 |
| P-04 | limit 快路径（`classify_*`, `first_order_exp_vanishing_epsilon`, …） | `limit_engine/exp_diff.rs` 等 | MRV/Gruntz 主路径 | → tech-debt **2A** 域 |
| P-05 | `giac-poly` factor `try_*` 链 | `factor/*` | FAC-G1 后 `factor_multivariate_rec` | → **D5**；见 [giac-poly-api-stability](../giac-poly-api-stability.md) §FAC |

---

## 3. 子任务总览

| ID | 标题 | 优先级 | 行为变化 | Blocked by | 状态 |
|----|------|--------|----------|------------|------|
| [D0-1](#giac-dedup-d0-1) | `ident_from_expr` 收拢到 `giac-core` | P0 | 无 | — | **done** |
| [D0-2](#giac-dedup-d0-2) | `is_var` / `var_to_expr` 收拢到 `expr_util` | P0 | 无 | — | **done** |
| [D0-3](#giac-dedup-d0-3) | 三角/对数形状检测共享 | P0 | 无 | D0-2 | **done** |
| [D0-4](#giac-dedup-d0-4) | `limit_engine` 内 `is_half_exponent` / `is_expr_zero` | P1 | 无 | — | **done** |
| [D1-1](#giac-dedup-d1-1) | `giac-poly` 二次系数提取共享核 | P1 | 无 | — | **done** |
| [D1-2](#giac-dedup-d1-2) | `giac-simplify/factor` 二次路径改调共享核 | P1 | 无（等价形） | D1-1 | **done** |
| [D1-3](#giac-dedup-d1-3) | 退役 `try_factor_quadratic_*` 重复逻辑 | P1 | 无 | D1-2, D3-2 | open |
| [D2-1](#giac-dedup-d2-1) | 全部 `try_integrate_*` 迁入 `integrate_heuristics` | P1 | 无 | — | **done** (PR-5) |
| [D2-2](#giac-dedup-d2-2) | `integrate.rs` 仅保留稳定规则 + 调度 | P2 | 无 | D2-1 | partial |
| [D3-1](#giac-dedup-d3-1) | `eval_solve` 接 `poly_algext_roots` | P0 | 有（能力扩展） | [F1–F5](GIAC-poly-quartic-roots-F1-F5.md), P3-6 | **wired** (PR-6); roots 层 `t⁴+t+1` 无 ignore |
| [D3-2](#giac-dedup-d3-2) | `quadratic_rootof_roots` 与 factor 共享二次核 | P1 | 无 | D1-1 | open |
| [D4-1](#giac-dedup-d4-1) | 按 tech-debt 2A/2B 删 drift / shim | P1 | 无 | 对应 follow-up issue | open |
| [D4-2](#giac-dedup-d4-2) | `assert_equiv` 漂移收敛（tech-debt 2C） | P2 | 可能改 golden | — | open |
| [D4-3](#giac-dedup-d4-3) | 删 `split_depressed_quartic` | P2 | 无 | D3-1 测试覆盖 | open |
| [D4-4](#giac-dedup-d4-4) | 清理空 `stubs.rs` 模块 | P3 | 无 | — | open |
| [D5-1](#giac-dedup-d5-1) | factor `try_*` 链收敛（FAC-G1 后） | P2 | 无 | giac-poly FAC-G1 | open |

---

## 4. 依赖图

```text
Phase D0（工具合并，无行为变化）
  D0-1 ident_from_expr ──────────────────────────┐
  D0-2 is_var / var_to_expr ──► D0-3 trig shapes   │
  D0-4 limit_engine 小工具（可并行）              │
                                                  │
Phase D1（二次多项式共享）                         │
  D1-1 giac-poly quadratic_coeffs ──► D1-2 simplify factor
                    └──► D3-2 rootof 共享          │
                    └──► D1-3 删重复 Temporary     │
                                                  │
Phase D2（integrate 文件整理）                     │
  D2-1 迁 try_integrate_* ──► D2-2 integrate 瘦身  │
                                                  │
Phase D3（solve ↔ poly_roots 接线）                │
  F1–F5 / P3-6 ──► D3-1 eval_solve 接 poly_algext_roots
  D1-1 ──► D3-2                                   │
  D3-1 ──► D4-3 删 legacy Ferrari                  │
                                                  │
Phase D4（临时 API 退役）                          │
  GIAC-limit-mrv-followup 3A ──► D4-1 (2A)       │
  GIAC-algext-adoption A-03 ──► D4-1 (2B)         │
  D4-2 assert_equiv (2C, HITL)                    │
  D4-4 stubs 清理（独立）                          │
                                                  │
Phase D5（giac-poly factor，独立长线）             │
  FAC-G1 ──► D5-1 try_* 收敛                      │
```

**建议抓取顺序（AFK）：** D0-1 → D0-2 → D0-3 → D1-1 → D2-1 →（等待 F1–F5）→ D3-1。

---

## 5. 子任务详情

### GIAC-dedup-D0-1

**标题:** `ident_from_expr` 收拢到 `giac-core`  
**优先级:** P0  
**落点:** `giac-core/src/ident.rs` 或新建 `giac-core/src/expr_parse.rs`（`pub(crate)`）

#### 现状

相同实现至少 6 份：

| 文件 |
|------|
| `giac-solve/solve.rs`, `sturm.rs`, `fsolve.rs` |
| `giac-core/eval_poly.rs` |
| `giac-calculus/limit.rs`, `series.rs` |

#### What to build

```rust
// giac-core — pub(crate) 或 crate 内模块
pub(crate) fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError>
```

各调用方删本地 `fn`，改 `use giac_core::...`（或 `giac_core::internal` 路径，避免扩大公开 API）。

#### Acceptance criteria

- [ ] 6 处本地定义删除，行为不变
- [ ] `cargo test-timeout` 全绿
- [ ] 不新增 `pub` 导出（除非写入 `giac-core` API stability doc）

---

### GIAC-dedup-D0-2

**标题:** `is_var` / `var_to_expr` 收拢到 `giac-calculus::expr_util`  
**优先级:** P0  
**落点:** `giac-calculus/src/expr_util.rs`

#### 现状

`is_var` 出现在 `diff`, `integrate`, `limit`, `limit_engine/{mod,mrv,asymptotic,exp_diff,mrv_lead_term,bounds}`, `risch/tower` 等。  
`var_to_expr` 出现在 `integrate`（`pub(crate)`）, `series`, `limit_engine/{mod,mrv,asymptotic,sparse_series}`。

`expr_util` 已有 **Stable** `depends_on_var` / `is_const_wrt`，缺 `is_var` / `var_to_expr`。

#### What to build

- 在 `expr_util` 增加 `pub(crate) fn is_var`、`pub(crate) fn var_to_expr`（或 **Stable** 若文档化契约）。
- `integrate.rs` 删除 `pub(crate)` 再导出，heuristics 改从 `expr_util` 引用。

#### Acceptance criteria

- [ ] `giac-calculus` 内 `fn is_var` / `fn var_to_expr` 定义 ≤1 处
- [ ] `*-api-stability.md` §2 或 §3 补一行（若标 Stable）
- [ ] 全绿测试

---

### GIAC-dedup-D0-3

**标题:** 三角/对数形状检测共享  
**优先级:** P0  
**Blocked by:** D0-2  
**落点:** `giac-calculus::expr_util`（跨 crate 需求时再下沉 `giac-core`）

#### 现状

| 函数 | 重复位置 |
|------|----------|
| `is_sin_of_var` | `solve`, `integrate`, `limit`, `ode/desolve` |
| `is_cos_of_var`, `is_ln_of_var` | `integrate`, `limit` |

#### What to build

- `pub(crate) fn is_sin_of_var(e: &ExprArc, var: &Ident) -> bool`（及 cos/ln）。
- `giac-solve`、`giac-ode` 依赖 `giac-calculus` **不合适** → 若 solve/ode 需要，下沉到 `giac-core::expr` 旁的小模块 `shape.rs`（`pub(crate)`）。

**决策（默认）：** 先放 `giac-core` 的 `pub(crate)` `expr_shape.rs`，避免 solve → calculus 依赖环。

#### Acceptance criteria

- [ ] 各 crate 本地 `is_sin_of_var` 删除
- [ ] 无新增 crate 循环依赖
- [ ] 全绿测试

---

### GIAC-dedup-D0-4

**标题:** `limit_engine` 内 `is_half_exponent` / `is_expr_zero` 合并  
**优先级:** P1  
**落点:** `giac-calculus/src/limit_engine/util.rs`（新建，或 `expr_util` 的 `pub(crate)` re-export）

#### 现状

`is_half_exponent`: `preprocess.rs`, `sparse_series.rs`, `asymptotic.rs`  
`is_expr_zero`: `sparse_series.rs`, `mrv_w.rs`

#### Acceptance criteria

- [ ] 各函数定义 1 处
- [ ] limit 相关测试全绿

---

### GIAC-dedup-D1-1

**标题:** `giac-poly` 二次系数提取共享核  
**优先级:** P1  
**落点:** `giac-poly/src/resultant.rs` 或 `giac-poly/src/univariate.rs`

#### What to build

```rust
/// (a, b, c) for a·var² + b·var + c; None if not quadratic univariate
pub fn quadratic_abc(p: &Poly, var: &Var) -> Option<(Ratio, Ratio, Ratio)>
```

供以下调用方使用：

- `factor/sqrt.rs::quadratic_sqrt_factor_exprs`
- `giac-simplify/factor.rs::try_factor_quadratic_*`
- `giac-solve/rootof.rs::quadratic_rootof_roots`（校验次数后）

#### Acceptance criteria

- [ ] `quadratic_sqrt_factor_exprs` 改调 `quadratic_abc`，单测 `sqrt_factor_x2_minus_2` 仍过
- [ ] 新 fn 标 **Stable** 或 **Partial** 并写入 `giac-poly-api-stability.md`

---

### GIAC-dedup-D1-2

**标题:** `giac-simplify/factor` 二次路径改调共享核  
**优先级:** P1  
**Blocked by:** D1-1

#### What to build

- `try_factor_quadratic_rootof` / `try_factor_quadratic_sqrt` 删重复系数循环，调用 `quadratic_abc` + `ratio_perfect_sqrt`。
- 可选：`try_factor_quadratic_rootof` 内部改调 `giac_solve::quadratic_rootof_roots` 再转 `(x - root)` 因子形（需评估 crate 依赖：simplify 是否可依赖 solve — **当前不可**；保持 simplify 内 rootof 构造或抽 `giac-core` 辅助）。

**约束:** `giac-simplify` 不得依赖 `giac-solve`；rootof 构造留在 simplify 或 `giac-core::AlgExtData`。

#### Acceptance criteria

- [ ] `factor.rs` 中系数提取循环删除
- [ ] factor 相关 golden / `assert_equiv` 全绿

---

### GIAC-dedup-D1-3

**标题:** 退役 simplify 二次 Temporary 重复逻辑  
**优先级:** P1  
**Blocked by:** D1-2, D3-2

#### What to build

当 D3-2 完成后，评估 `try_factor_quadratic_rootof` 是否可由 `eval_solve` / 统一二次 API 覆盖；若 factor 热路径仍需 display 形，保留薄包装并标 **Pipeline private**，去掉 **Temporary** 标记。

#### Acceptance criteria

- [ ] `giac-simplify-api-stability.md` §4 Temporary 表更新
- [ ] 无净增形状分支

---

### GIAC-dedup-D2-1

**标题:** 全部 `try_integrate_*` 迁入 `integrate_heuristics.rs`  
**优先级:** P1  
**落点:** `giac-calculus/src/integrate_heuristics.rs`

#### 现状

`integrate.rs` 仍含 5 个 `try_integrate_*`（`sin2x_cos`, `var_shifted_sqrt`, `sin_over_cos_squared`, `exp_trig`, `var_over_quadratic_squared`），且 `integrate()` 主函数直接调用 `try_integrate_tan_plus_tan_cubed` 等；`integrate_heuristics` 又从 `integrate` import 一批 `try_integrate_*` — **双向耦合**。

#### What to build

1. 将 `integrate.rs` 内所有 `try_integrate_*` 移至 `integrate_heuristics.rs`（或 `integrate/rules.rs`）。
2. `integrate()` 入口仅：`try_integrate_heuristic` → `try_as_rational` → 稳定 `match` 分派。
3. `integrate_heuristics` 不再从 `integrate` import `try_integrate_*`；共享 `is_const_wrt` 等从 `expr_util` 引用。

#### Acceptance criteria

- [x] `integrate.rs` 无 `try_integrate_*` 定义（迁至 `integrate_try_rules.rs`）
- [ ] `giac-calculus-api-stability.md` inventory 模块归属更新
- [x] integrate 单测全绿

---

### GIAC-dedup-D2-2

**标题:** `integrate.rs` 仅保留稳定规则 + 调度  
**优先级:** P2  
**Blocked by:** D2-1

#### Acceptance criteria

- [ ] `integrate.rs` 行数显著下降（启发式全在 heuristics 模块）
- [ ] 文档 §4 Partial 表仍准确（退役计划不变）

---

### GIAC-dedup-D3-1

**标题:** `eval_solve` 接 `poly_algext_roots`  
**优先级:** P0（能力项）  
**Blocked by:** [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md)、[P3-6](GIAC-poly-p3-6-quartic-roots-gaps.md)  
**落点:** `giac-solve/solve.rs`

#### 现状

```text
eval_solve: giac_poly::roots (ℚ) → NotImplemented → quadratic_rootof_roots / biquadratic_rootof_roots
poly_algext_roots: FieldSession 完整代数求根（giac-core，已 export）
```

二者未接线；一般四次仍 `NotImplemented("general quartic rootof")`。

#### What to build

1. `equation_to_poly` 后，若 `giac_poly::roots` 返回 `NotImplemented`：
   - `poly_alg_from_expr` → `poly_algext_roots` → `algext_poly_to_expr` 列表。
2. 删除或收窄 `biquadratic_rootof_roots` 旁路（保留二次 `rootof` 作 fast path 可选）。
3. `solve(t^4+t+1=0,t)` 四根验收（与 P3-6 normative 一致）。

#### Acceptance criteria

- [ ] `GIAC-poly-p3-6` §1 验收满足（`t⁴+t+1` 待 F1–F5）
- [ ] `giac-solve-api-stability.md` §5「一般四次」行改为 done 或链 F1–F5
- [x] 无在 solve 新写 Ferrari 特判
- [x] `poly_algext_from_poly` + `poly_algext_roots` 回退链（PR-6）

---

### GIAC-dedup-D3-2

**标题:** `quadratic_rootof_roots` 与 factor 共享二次核  
**优先级:** P1  
**Blocked by:** D1-1  
**落点:** `giac-solve/rootof.rs`

#### What to build

`quadratic_rootof_roots` 用 `quadratic_abc` 校验；与 `giac-poly::quadratic_sqrt_factor_exprs` 共享判别式逻辑。

#### Acceptance criteria

- [ ] `rootof` 测试 `quadratic_rootof_has_two_branches` 仍过
- [ ] 无重复系数循环

---

### GIAC-dedup-D4-1

**标题:** 按 tech-debt 2A/2B 删 drift / shim  
**优先级:** P1  
**落点:** 见 [GIAC-expr-api-tech-debt §2A/2B](GIAC-expr-api-tech-debt.md)

不重复 spec；本计划仅提醒：**删 `drift_*` / `ratnormal_algext` 必须与 MRV Phase 3A / ext_reduce 同 PR 或紧随其后**，避免测试绿但语义回退。

---

### GIAC-dedup-D4-2

**标题:** `assert_equiv` 漂移收敛（tech-debt 2C）  
**优先级:** P2  
**类型:** HITL（可能触及 golden）

见 tech-debt **2C**；`canonical_radical` 迁入 `normal` 后删除 T-02。

---

### GIAC-dedup-D4-3

**标题:** 删 `split_depressed_quartic`  
**优先级:** P2  
**Blocked by:** D3-1  
**落点:** `giac-core/algebra/poly_roots.rs`

#### What to build

- 将依赖 legacy Ferrari 的测试改为 Euler / `poly_algext_roots` 路径。
- 删除 `split_depressed_quartic` 及仅其使用的 pipeline private 辅助（若有）。

#### Acceptance criteria

- [ ] `poly_roots` 测试无 legacy 路径依赖
- [ ] `giac-core-algebra-api-stability.md` 更新

---

### GIAC-dedup-D4-4

**标题:** 清理空 `stubs.rs` 模块  
**优先级:** P3  
**落点:** `giac-calculus/src/lib.rs`, `giac-ode/src/lib.rs`

#### What to build

- 删除 `mod stubs;` 及空文件，或保留单文件 doc 注释说明历史（团队择一）。
- `giac-solve/stubs.rs` 保留（有 re-export）。

#### Acceptance criteria

- [ ] `cargo build -p giac-calculus -p giac-ode` 无警告
- [ ] api-stability doc 模块表同步

---

### GIAC-dedup-D5-1

**标题:** factor `try_*` 链收敛（FAC-G1 后）  
**优先级:** P2  
**Blocked by:** giac-poly FAC-G1  
**落点:** `giac-poly/src/factor/mod.rs`

见 [giac-poly-api-stability.md](../giac-poly-api-stability.md) §FAC-G1–G3 **退役目标**。本计划不拆子项；FAC-G1 PR 描述须引用本文 §2.2 P-05。

---

## 6. 明确不合并（避免误伤）

| 项 | 理由 |
|----|------|
| `poly_coeff` trait vs `poly_alg_coeff` | 不同系数类型（`Ratio` vs `AlgExtC`），有意分层 |
| 各 crate `plugin.rs` + `xcas_default()` | 组合模式，非重复 |
| `rational_num_den` 三版本 | 语义不同（静态 / eval 后 / 负幂 Mul）；合并需 `enum RationalPartsMode`，单独立项 |
| `limit_engine` 大量 pipeline private | 模块内递归；跨文件合并降低可读性 |
| `giac-poly::nested` 类型族 | 环语义显式化，见 nested-ring doc |

---

## 7. 验收门禁（父计划）

- [ ] D0 全完成：工具 fn 拷贝数较审计基线减少 ≥80%
- [ ] D3-1 完成：`solve(t^4+t+1=0,t)` 四根（P3-6 normative）
- [ ] 每个子 PR：`cargo test-timeout` + `cargo ci-clippy` 全绿
- [ ] 每个子 PR：§7.2 两行摘要（删除的临时匹配 / 新增 fn tier）
- [ ] Temporary 表（§2.1）项逐项 closed 或链到已完成 follow-up issue

---

## 8. PR 切片建议

| PR | 子项 | 预估规模 |
|----|------|----------|
| PR-1 | D0-1 | 小（~6 文件删重复） |
| PR-2 | D0-2 + D0-3 | 中（calculus + core shape） |
| PR-3 | D0-4 | 小 |
| PR-4 | D1-1 + D1-2 | 中 |
| PR-5 | D2-1 | 中（搬函数，无逻辑改） |
| PR-6 | D3-1 | 大（依赖 F1–F5） |
| PR-7 | D1-3 + D3-2 + D4-3 | 中 |
| PR-8+ | D4-1, D4-2, D5-1 | 随上游 issue 并行 |

---

## 9. 参考

- [GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md) — API 分层母索引
- [GIAC-poly-p3-6-quartic-roots-gaps](GIAC-poly-p3-6-quartic-roots-gaps.md) — 四次求根验收
- [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) — 塔修复竖切
- [giac-simplify-api-stability.md](../giac-simplify-api-stability.md) §4 Temporary
- [giac-solve-api-stability.md](../giac-solve-api-stability.md) §5 演化
- [giac-calculus-api-stability.md](../giac-calculus-api-stability.md) §4 integrate/limit Partial
- [giac-poly-api-stability.md](../giac-poly-api-stability.md) §FAC-G*
