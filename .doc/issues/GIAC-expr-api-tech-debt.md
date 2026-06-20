# GIAC-expr-api-tech-debt — 按 algorithm-expr-api 消除表示层 / API / 测试技术债

**状态:** open  
**类型:** 索引 / AFK（子项见下表）  
**依据:** [algorithm-expr-api.md](../algorithm-expr-api.md) §1–§7；配套 [algorithm-expr-api.mdc](../../.cursor/rules/algorithm-expr-api.mdc)、[algorithm-before-patch.mdc](../../.cursor/rules/algorithm-before-patch.mdc)  
**不重复:** 算法能力缺口仍归 [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)；本索引只清 **API 分层、显式上下文、临时 API 退役、测试契约、边界审计**。  
**快照:** 2026-06-20

---

## 1. 问题陈述

giac-rs 算法 crate 在 [algorithm-expr-api.md](../algorithm-expr-api.md) 落地程度上不均：

| 维度 | 现状 | 文档要求 |
|------|------|----------|
| API tier + 稳定性文档 | `giac-simplify` / `giac-poly` 有 doc + `annotate_api_tiers.py`；`giac-calculus` 有 doc 无 inventory；`giac-core/algebra` ~364 fn 仅 7 处 tier；`giac-solve` / `giac-ode` / `giac-groebner` **零 tier** | §2 三层分类；§7.2 inventory |
| 显式算法上下文 | `poly_roots.rs` 散落 `align_coeff`；**无 `FieldSession`** | §6.3 |
| 临时 API | `mrv_w` 8× `drift_*`；`equiv::canonical_radical`；`ratnormal_algext` shim | §2 须写退役条件 + issue |
| 测试契约 | 全库 ~120+ 处 `format_expr(...).contains("…")` 语义断言；`assert_equiv` ~30 处 | §3 / §5 |
| §4 索引 | `giac-solve` / `giac-groebner` / `giac-ode` / `giac-core::algebra` 未入表 | §4 |

**相对健康：** 无 `looks_like_*` 泛滥；`expr_to_poly` 对 `AlgExt` 显式拒绝；limit_engine `canonical_mrv_coeff` / `exp_diff` 契约文档齐全。

---

## 2. 子 issue 总览（tracker）

| ID | 标题 | Type | 优先级 | Blocked by | 状态 |
|----|------|------|--------|------------|------|
| [0A](#giac-expr-api-0a) | 扩展 `annotate_api_tiers.py` + inventory | AFK | P0 | — | **done** |
| [0B](#giac-expr-api-0b) | 补齐 `*-api-stability.md` + §4 索引 | AFK | P0 | 0A | **done** |
| [1A](#giac-expr-api-1a) | `FieldSession` 落地 | AFK | P0 | — | open |
| [1B](#giac-expr-api-1b) | `PolyInK` 入口规范化 | AFK | P0 | 1A | open |
| [2A](#giac-expr-api-2a) | 退役 `mrv_w::drift_*` | AFK | P1 | [GIAC-limit-mrv-followup](GIAC-limit-mrv-followup.md) Phase 3A | open |
| [2B](#giac-expr-api-2b) | 退役 `ratnormal_algext` shim | AFK | P2 | [GIAC-algext-adoption](GIAC-algext-adoption.md) A-03 | open |
| [2C](#giac-expr-api-2c) | `assert_equiv` 漂移收敛 | HITL | P2 | — | open |
| [3A](#giac-expr-api-3a) | `giac-solve` 测试去字符串语义 | AFK | P1 | 0B | **done** |
| [3B](#giac-expr-api-3b) | `giac-calculus` 积分 / limit 单测契约 | AFK | P1 | — | **partial** |
| [3D](#giac-expr-api-3d) | `giac-core/algebra` 扩域单测 | AFK | P1 | 1A | **partial** |
| [4A](#giac-expr-api-4a) | `expr_to_poly` 调用方审计 | AFK | P2 | 0B | open |
| [4B](#giac-expr-api-4b) | factor 静默 fallback 显式化 | AFK | P2 | — | open |
| [4C](#giac-expr-api-4c) | 嵌套环上下文 / `div_rem` 误用审计 | AFK | P2 | [GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md) | open |

**与已有 issue 的直接映射（不重复开项）：**

| 本索引 | 已有 issue | 关系 |
|--------|-----------|------|
| 1A / 1B | [GIAC-poly-roots-field-session-plan](GIAC-poly-roots-field-session-plan.md) | **执行该方案**，本索引仅跟踪 |
| 1A | [GIAC-poly-p3-6-quartic-roots-gaps](GIAC-poly-p3-6-quartic-roots-gaps.md) G5 | FieldSession 验收项 |
| 2A | [GIAC-limit-mrv-followup](GIAC-limit-mrv-followup.md) Phase 3A | 子任务 |
| 2B | [GIAC-algext-adoption](GIAC-algext-adoption.md) A-03 | 子任务 |
| 4C | [GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md) | 子任务 |

---

## 3. 依赖图

```text
Phase 0 (文档基线)
  0A inventory 脚本 ──► 0B 稳定性文档 + §4 索引

Phase 1 (显式上下文)
  1A FieldSession ──► 1B PolyInK
                 └──► 3D algebra 单测

Phase 2 (临时 API 退役)
  GIAC-limit-mrv-followup 3A ──► 2A drift_*
  GIAC-algext-adoption A-03 ──► 2B ratnormal_algext
  2C assert_equiv canonical (HITL，可并行)

Phase 3 (测试契约)
  0B ──► 3A solve 测试
  3A, 3B ──► 3C conformance
  1A ──► 3D algebra 测试
  3B 可与 0B 并行

Phase 4 (边界审计)
  0B ──► 4A expr_to_poly 审计
  4B factor 静默失败 (独立)
  GIAC-poly-nested-ring-types ──► 4C div_rem 审计
```

---

## 4. 验收门禁（父 issue）

- [ ] 子 issue 0A–4C 均有 tracker 条目且状态可跟踪
- [ ] `cargo test-timeout` 全绿（各子 PR 独立保证）
- [ ] 算法 crate PR 描述含 §7.2 两行摘要：「删除/合并的临时匹配」「新增 fn + tier」
- [ ] [algorithm-expr-api.md §4](../algorithm-expr-api.md#4-各-crate-规范入口索引) 索引与 `*-api-stability.md` 一致

---

# 子 issue 详情（tracker 格式）

---

## GIAC-expr-api-0A

**标题:** 扩展 `annotate_api_tiers.py` + per-crate inventory  
**Type:** AFK  
**优先级:** P0  
**Blocked by:** None — 可立即开始  
**落点:** `giac-rs/scripts/annotate_api_tiers.py`

### What to build

将 `annotate_api_tiers.py` 从仅 `giac-simplify` + `giac-poly` 扩展到：

- `giac-calculus`
- `giac-core/src/algebra`
- `giac-solve`
- `giac-ode`
- `giac-groebner`

`--inventory` 输出各 crate 函数 tier 统计；未标注私有 `fn` 自动补 `// **Pipeline private**`（或输出待补清单）。更新脚本 docstring 与贡献说明引用 [algorithm-expr-api.md §7.2](../algorithm-expr-api.md#72-测试通过后提交--合入前复审)。

### Acceptance criteria

- [ ] 五 crate 均可 `python3 scripts/annotate_api_tiers.py --inventory` 刷新
- [ ] idempotent：已有 tier 标记的 fn 不被覆盖
- [ ] inventory 计数写入对应 `*-api-stability.md` Per-file 表（或 0B 消费该输出）

---

## GIAC-expr-api-0B

**标题:** 补齐 `*-api-stability.md` + algorithm-expr-api §4 索引  
**Type:** AFK  
**优先级:** P0  
**Blocked by:** [0A](#giac-expr-api-0a)（可并行起草，0A 后对齐 inventory）  
**落点:** `.doc/giac-*-api-stability.md`；[algorithm-expr-api.md §4](../algorithm-expr-api.md)

### What to build

新增稳定性文档：

| 文档 | Crate / 模块 |
|------|-------------|
| `giac-core-algebra-api-stability.md` | `giac-core::algebra`（`poly`/`poly_conv`/`poly_roots`/`ext_tower`/`alg_ext`/`field_arith`） |
| `giac-solve-api-stability.md` | `giac-solve` |
| `giac-ode-api-stability.md` | `giac-ode` |
| `giac-groebner-api-stability.md` | `giac-groebner` |

每份含：§1 注释格式、§2 公开 API 表、§3 I/O 契约（§3 模板）。更新 [algorithm-expr-api.md §4](../algorithm-expr-api.md#4-各-crate-规范入口索引) 索引行。

### Acceptance criteria

- [ ] 四份新 doc 存在且与源码 `/// **Stable**` 等标记一致
- [ ] `poly_alg_from_expr` / `poly_algext_roots` / `solve` / `desolve` 等有契约表
- [ ] `algorithm-expr-api.md` §4 补全四 crate 行

---

## GIAC-expr-api-1A

**标题:** `FieldSession` 落地（扩域求根显式上下文）  
**Type:** AFK  
**优先级:** P0  
**Blocked by:** None  
**落点:** `giac-core/src/algebra/field_session.rs`（新）；`poly_roots.rs`  
**执行方案:** [GIAC-poly-roots-field-session-plan](GIAC-poly-roots-field-session-plan.md)

### What to build

实现 `FieldSession { ambient: K, working: L }` 及 §6.3 最小 API（`.int` / `.half` / `.zero` / `.one` / `.lift` / `.adjoin_sqrt` / `.align`）。`poly_algext_roots` 子算法签名带 `&FieldSession`；Cardano / resolvent / split 在 **working L** 上构造常数；消除 `poly_roots.rs` 内散落 `align_coeff` + `ring_int(K,…)`。

### Acceptance criteria

- [ ] 现有 7 个 `poly_roots` 单测保持绿
- [ ] `t⁴+t+1=0` 单测去 `#[ignore]`，运行 ≤10s
- [ ] `poly_roots.rs` 内 `align_coeff` 收敛为 `session.align`（或删除私有副本）
- [ ] [GIAC-poly-p3-6-quartic-roots-gaps](GIAC-poly-p3-6-quartic-roots-gaps.md) G5 满足

---

## GIAC-expr-api-1B

**标题:** `PolyInK` 入口规范化包装  
**Type:** AFK  
**优先级:** P0  
**Blocked by:** [1A](#giac-expr-api-1a)  
**落点:** `giac-core/src/algebra/poly_roots.rs`

### What to build

私有 `PolyInK { poly, ambient }` 包装类型；禁止未 normalize 的多项式进入 `roots_dispatch`。入口固定单路径：`infer_field → normalize_coeffs → monic_univariate → FieldSession::new(K) → roots_dispatch`。

### Acceptance criteria

- [ ] 外部无法构造未 normalize 的 `PolyInK`
- [ ] 文档 [algorithm-expr-api.md §6.3](../algorithm-expr-api.md#63-实例扩域求根-fieldsessionpolyc-over-k) 最小 API 草图与实现一致
- [ ] `verify` 与算法共用同一 normalize + monic 路径

---

## GIAC-expr-api-2A

**标题:** 退役 `mrv_w::drift_*`（MRV 系数 AST 漂移吸收）  
**Type:** AFK  
**优先级:** P1  
**Blocked by:** [GIAC-limit-mrv-followup](GIAC-limit-mrv-followup.md) **Phase 3A**（`giac-simplify` 受限 expand）  
**落点:** `giac-calculus/src/limit_engine/mrv_w.rs`

### What to build

Phase 3A 落地后删除 8 个 `drift_*`（`drift_fold_ln_atoms`、`drift_is_neg_ln_shape` 等）。`canonical_mrv_coeff` 单测保留「化简后漂移形态」覆盖，但只经规范入口识别，不经私有 drift 表。

### Acceptance criteria

- [ ] `mrv_w.rs` 无 `drift_*` 函数
- [ ] `limit_engine` 相关单测 + CK-INT 子集不回归
- [ ] [giac-calculus-api-stability.md](../giac-calculus-api-stability.md) Temporary 表更新

---

## GIAC-expr-api-2B

**标题:** 退役 `ratnormal_algext` shim  
**Type:** AFK  
**优先级:** P2  
**Blocked by:** [GIAC-algext-adoption](GIAC-algext-adoption.md) **A-03**（AlgExt 原生有理化）  
**落点:** `giac-simplify/src/ratnormal.rs`

### What to build

`ratnormal` 对含 `AlgExt` 式子走稳定 AlgExt 有理化路径，删除 `ratnormal_algext` eval 绕行 shim。

### Acceptance criteria

- [ ] `ratnormal_algext` 删除
- [ ] `ratnormal_algext_square_minus_two` 等单测改走新路径
- [ ] [giac-simplify-api-stability.md](../giac-simplify-api-stability.md) Temporary 表更新

---

## GIAC-expr-api-2C

**标题:** `assert_equiv` 漂移收敛（`canonical_radical` 归属决策）  
**Type:** HITL  
**优先级:** P2  
**Blocked by:** None  
**落点:** `giac-simplify/src/equiv.rs`；可能 `giac-simplify` 公开 API

### What to build

决策 `canonical_radical` / `inv_sqrt_to_mul` 归宿：

- **选项 A：** 升格为 `giac-simplify` 稳定 `canonical_radical_expr`（有 I/O 契约 + 单测）
- **选项 B：** 并入 `normal` 主路径，测试层不再依赖私有 drift

Conformance 与单元测试覆盖 radical 往返 + 化简后形态。

### Acceptance criteria

- [ ] HITL 决策记录在 issue 评论或 ADR 一行
- [ ] `equiv.rs` 无未文档化的 Temporary drift（或已升格 Stable）
- [ ] conformance `assert_equiv` 覆盖 `1/sqrt(n) ↔ sqrt(n)/n` 类形态

---

## GIAC-expr-api-3A

**标题:** `giac-solve` 测试去字符串语义断言  
**Type:** AFK  
**优先级:** P1  
**Blocked by:** [0B](#giac-expr-api-0b)  
**落点:** `giac-solve/src/{rootof,realroot,froot,solve}.rs`

### What to build

将 `format_expr(...).contains("rootof")` 等字符串语义断言替换为：

- `try_as_algext_data` / `expr_contains_alg_coeff`
- 解代入原方程（`assert_equiv` 或多项式 remainder 为零）
- `eq_mod`（扩域根验证）

在 [giac-solve-api-stability.md](../giac-solve-api-stability.md)（0B）中文档化 `solve` / `rootof` 输出契约。

### Acceptance criteria

- [ ] [test-writing-spec.md](../test-writing-spec.md) 分层验收：每个 touched 测试在 [GIAC-expr-api-test-audit.md](GIAC-expr-api-test-audit.md) 有 A/B/C 行
- [ ] 每个公开 solve **eval 入口**至少 **1 个 A**（代入验证 / `eq_mod`）
- [ ] `giac-solve` 内零 `contains` **语义**断言（**C** display 须标注）
- [ ] `cargo test-timeout -p giac-solve` 全绿

**进展：** rootof/froot/realroot 已 A 化；`solve_linear_system` ⚠️（`eval_const_expr` 待迁入 `test_verify`，见 audit §2）。

---

## GIAC-expr-api-3B

**标题:** `giac-calculus` 积分 / limit 单测契约升级  
**Type:** AFK  
**优先级:** P1  
**Blocked by:** None  
**落点:** `integrate.rs`、`exp_diff.rs`、`preprocess.rs`（优先）；其余 calculus 单测次之

### What to build

按 [algorithm-expr-api.md §3](../algorithm-expr-api.md#3-稳定-api-契约模板)：

- **积分：** `diff` 还原被积函数，或 `assert_equiv`
- **极限：** 数值探测 / `assert_equiv`，非 `format_expr.contains` 认形状
- **MRV / exp_diff：** 保留 `_mrv_w` 等 **display 快照** 时与 **语义** 断言分离

### Acceptance criteria

- [ ] [test-writing-spec.md](../test-writing-spec.md) 分层验收 + [GIAC-expr-api-test-audit.md](GIAC-expr-api-test-audit.md) 更新
- [ ] 每个 touched **Stable** API 至少 **1 个 B**；极限 fixture 至少 **1 个 A′**（数学结果，非中间 `contains`）
- [ ] `integrate.rs`：补 **A**（`eval(Integrate)`）或登记 gap；`contains` 语义清零或标 **C**
- [ ] `cargo test-timeout -p giac-calculus` 全绿

**进展（2026-06-20）：** `exp_diff.rs` / `preprocess.rs` 语义 `contains` 已清零；`integrate.rs` 仍为最大缺口（audit §5）。

---

## GIAC-expr-api-3C

**标题:** conformance 层统一等价策略  
**Type:** AFK  
**优先级:** P2  
**Blocked by:** [3A](#giac-expr-api-3a), [3B](#giac-expr-api-3b)  
**落点:** `giac-rs/tests/conformance/`；[conformance-testing.md](../conformance-testing.md)

### What to build

审计 conformance 与 `phase*_*.rs`：字面 golden 失败时优先 `assert_equiv`（§3）。在 [conformance-testing.md](../conformance-testing.md) 补充「何时字面 / 何时 assert_equiv」决策表；新 conformance 默认不用 `contains` 判语义。**单测分层**见 [test-writing-spec.md](../test-writing-spec.md)。

### Acceptance criteria

- [ ] 决策表写入 `conformance-testing.md`
- [ ] `tests/conformance/` 中新增用例遵循该表
- [ ] 现有 `contains` 断言分类为「display 快照」或改为等价判定

---

## GIAC-expr-api-3D

**标题:** `giac-core/algebra` 扩域单测契约  
**Type:** AFK  
**优先级:** P1  
**Blocked by:** [1A](#giac-expr-api-1a)  
**落点:** `alg_ext.rs`、`poly_roots.rs` 单测

### What to build

扩域相关单测用 `eq_mod`、已知根代入、`poly_algext_roots` 验证，替代 `contains("rootof")`。验证路径与 `FieldSession` 共用 normalize + monic。

### Acceptance criteria

- [ ] `alg_ext.rs` / `poly_roots.rs` 无字符串语义断言
- [ ] 四次根单测（1A 启用后）含 `eq_mod` 两两对齐
- [ ] `cargo test-timeout -p giac-core` algebra 相关全绿

---

## GIAC-expr-api-4A

**标题:** `expr_to_poly` 调用方审计  
**Type:** AFK  
**优先级:** P2  
**Blocked by:** [0B](#giac-expr-api-0b)  
**落点:** `giac-solve`、`giac-calculus`、`giac-simplify`；[expr-poly-conversion.md](../expr-poly-conversion.md)

### What to build

清单化各 crate 中 `expr_to_poly` 调用：是否先 `eval`、是否可能含超越函数、是否应改 `poly_alg_from_expr`。缺 guard 处加 `expr_contains_alg_coeff` 或文档化前置条件。禁止新增「先 `expr_to_poly` 再认 exp/ln 形状」路径。

### Acceptance criteria

- [ ] 审计表写入 `expr-poly-conversion.md` 或 `giac-solve-api-stability.md`
- [ ] 无未文档化的 transcendental → poly 误用
- [ ] grep 基线 + 回归检查纳入 PR 模板说明

---

## GIAC-expr-api-4B

**标题:** factor 静默 fallback 显式化  
**Type:** AFK  
**优先级:** P2  
**Blocked by:** None  
**落点:** `giac-poly/src/factor/`；[giac-poly-api-stability.md](../giac-poly-api-stability.md)

### What to build

审查 `factor/univariate.rs` 等 `Ok(vec![g.clone()])` 路径：区分 **Stable (bounded)** 与应 `None` / `Err` 的路径。对齐 [algorithm-expr-api.md §7.2](../algorithm-expr-api.md#72-测试通过后提交--合入前复审)「静默失败」审查项——不可约时不得假成功因子列表。

### Acceptance criteria

- [ ] 每条 fallback 在 `giac-poly-api-stability.md` 有边界说明
- [ ] 无未文档化的 `Ok(vec![g])` 当不可约
- [ ] 相关单测区分「bounded 成功」与「应失败」

---

## GIAC-expr-api-4C

**标题:** 嵌套环上下文 / `div_rem` 误用审计  
**Type:** AFK  
**优先级:** P2  
**Blocked by:** [GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)  
**落点:** `giac-poly`；[giac-poly-nested-ring-types.md](GIAC-poly-nested-ring-types.md)

### What to build

grep 审计 `Poly::div_rem` 用于验 ℚ[others][main] 整除的误用；推广 `UnivariateIn::divides` / `quo_exact_wrt` / `SqffRingCtx`。新代码禁止裸树 + 隐式主元假设。

### Acceptance criteria

- [ ] 误用清单归零或每条登记 gap issue
- [ ] [algorithm-expr-api.md §5](../algorithm-expr-api.md#5-反模式全-crate) 反模式表与代码一致
- [ ] 新 PR 在 giac-poly 改动时引用本审计

---

## 5. PR 检查清单（复制到算法 crate PR）

测试全绿后 **必须** 复审 diff（[§7.2](../algorithm-expr-api.md#72-测试通过后提交--合入前复审)）：

| 审查项 | 通过？ |
|--------|--------|
| 新增/改单测已更新 [GIAC-expr-api-test-audit.md](GIAC-expr-api-test-audit.md) | |
| Stable API 有 B；eval 入口有 A（见 [test-writing-spec.md](../test-writing-spec.md)） | |
| 临时匹配（`looks_like` / per-case / 调用方 shim）净减少或 justify + gap issue | |
| 新主路径已覆盖的旧 `drift_*` / `shim_*` 已删或标退役 | |
| 新增 / 改签名 `fn` 已标 tier 并更新 `*-api-stability.md` | |
| 算法 crate 改动时跑过 `annotate_api_tiers.py --inventory` | |
| 无静默 factor / try 假成功 | |

**PR 描述两行摘要：**

1. 删除 / 合并的临时匹配：…
2. 新增 fn + tier：…

---

## 6. 参考

- [test-writing-spec.md](../test-writing-spec.md) — A / A′ / B / C 单测编写规范
- [GIAC-expr-api-test-audit.md](GIAC-expr-api-test-audit.md) — 3A/3B/3D 审计表
- [algorithm-expr-api.md](../algorithm-expr-api.md)
- [conformance-testing.md](../conformance-testing.md) §3 `assert_equiv`
- [expr-poly-conversion.md](../expr-poly-conversion.md)
- [giac-calculus-api-stability.md](../giac-calculus-api-stability.md)
- [giac-poly-api-stability.md](../giac-poly-api-stability.md)
- [giac-simplify-api-stability.md](../giac-simplify-api-stability.md)
