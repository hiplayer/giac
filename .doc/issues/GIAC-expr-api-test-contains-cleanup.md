# GIAC-expr-api — crate 单测 `contains` 语义债清理

**状态:** open  
**类型:** AFK（测试契约 / 共享 helper）  
**父项:** [GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md) §3B 延伸；[test-writing-spec.md](../test-writing-spec.md) §2（**禁止** `contains`）  
**审计表:** [GIAC-expr-api-test-audit.md](GIAC-expr-api-test-audit.md)（逐测登记）  
**快照:** 2026-06-24（`giac-rs/crates` 内 `assert!(…contains…)` 约 **66** 处；**P0 H1/H2 已落地**）

---

## 1. 问题陈述

giac-rs 算法 crate 内大量单测用 `format_expr` + `contains` 做语义断言。这在 CAS 里等价于「输出字符串里出现过某子串」，**不能**证明数学正确，且与 [test-writing-spec.md](../test-writing-spec.md) A/B/C 分层冲突。

| 维度 | 现状 | 目标 |
|------|------|------|
| 语义断言 | `assert_equiv` / `matches!` / 代入方程 — 约 **50** 处，集中在 solve、integrate 部分 | A/B 层首选数学等价或结构匹配 |
| 弱断言 | `contains` — 约 **66** 处，linalg/series/Risch 最重 | **清零**语义 `contains`；仅 C 层完整 `assert_eq!` golden |
| 共享 helper | `giac-solve::test_verify` 已有 | calculus / linalg / ode / simplify 对齐同一模式 |
| 父项 3B | integrate/diff/limit 部分已 A 化 | 本 issue 覆盖 **剩余 crate + 基础设施** |

**不算本 issue 的 `contains`：** 生产代码 `Vec::contains`；`lexer` 的 `Token` 枚举检查；`expr_contains_*` 谓词函数名。

---

## 2. 优先级总览

| 优先级 | ID | 标题 | 处数（约） | Blocked by | 类型 |
|--------|-----|------|-----------|------------|------|
| **P0** | [H1](#h1-calculus-test_verify) | `giac-calculus::test_verify`（微分还原 + 级数截断） | 解锁 ~30 | — | **done** |
| **P0** | [H2](#h2-linalg-test_verify) | `giac-linalg::test_verify`（矩阵 / charpoly） | 解锁 ~20 | — | **done** |
| **P1** | [C1](#c1-simplify-contains) | `giac-simplify` 插件 + expand/trig/ratnormal | ~11 | — | 改测 |
| **P1** | [C2](#c2-limit-engine-contains) | `limit_engine` remove_lnexp / mrv / mrv_w | ~9 | — | 改测 |
| **P1** | [C3](#c3-core-eval-contains) | `giac-core` eval + eval_poly_tests | ~8 | H6（egcd 部分） | 改测 |
| **P1** | [H3](#h3-ode-test_verify) | `giac-ode::test_verify`（ODE 代入验算） | 解锁 2 | — | helper |
| **P2** | [C4](#c4-series-sparse-series) | `series.rs` + `sparse_series.rs` | ~13 | **H1** I4 | 改测 |
| **P2** | [C5](#c5-risch-heuristics) | `risch/*` + `integrate_heuristics.rs` | ~16 | **H1** I1 / T3 | 改测 |
| **P2** | [C6](#c6-linalg-coverage) | `phase3_coverage` + symbolic + gramschmidt + eigen | ~20 | **H2** | 改测 |
| **P2** | [H4](#h4-simplify-factor-verify) | expr 层 `assert_factorization` | 解锁 3 | — | helper |
| **P3** | [C7](#c7-ode-plugin) | `giac-ode` plugin 通路 | 1 | **H3** | 改测 |
| **P3** | [C8](#c8-wasm-smoke) | `giac-wasm` 两条 | 2 | — | 改测 |
| **P3** | [H5](#h5-poly-sqrt-factor-api) | `quadratic_sqrt_factor_exprs` 返回 `Expr` | 1 | API | 表示层 |
| **P3** | [C9](#c9-parser-debug) | `parser` Debug 子串 → `matches!` | 5 | — | 改测 |

```text
P0  helper 地基
  H1 calculus:test_verify ──┬──► C4 series
                            └──► C5 risch/heuristics (部分需 T3)
  H2 linalg:test_verify ────────► C6 linalg 全套

P1  无阻塞快改（可与 P0 并行）
  C1 simplify · C2 limit_engine · C3 core · H3 ode helper · H4 factor verify

P2  依赖 P0
  C4 · C5 · C6 · H4 落地后改 factor 测

P3  长尾
  C7 · C8 · H5 · C9
```

---

## 3. 基础设施（P0–P1）

### H1 — calculus `test_verify` {#h1-calculus-test_verify}

**文件:** `giac-calculus/src/test_verify.rs`（新建，仿 `giac-solve`）

| 函数 | 用途 | 解锁测试 |
|------|------|----------|
| `assert_deriv_equals_integrand(f, var, F, ctx)` | `diff(F,var)` 与 `f` 做 `assert_equiv` | Risch、heuristics、smoke 积分 |
| `assert_series_equiv_at(f, var, center, order, expected, ctx)` | 截断级数逐项 `assert_equiv`（或 `series` eval 后与期望多项式比） | `series.rs` 9、`sparse_series.rs` 4 |

**验收:** `cargo test-timeout -p giac-calculus test_verify`；`series_exp_at_zero`、`pow2expln_*` 已改 helper。

**已知:** `assert_linsolve_satisfies` 已实现但 **未接线** phase3（`8/3+1/3-3` 有理和 `normal` 未归零，见 T1/2C）。

**阻塞登记（T3，不阻塞 helper 本身）:** `diff(integrate(f))` 对 `ln(abs(x))`、部分代数原函数可能失败 — 失败用例登记 [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)，暂保留 `assert_equiv` 对已知闭式。

---

### H2 — linalg `test_verify` {#h2-linalg-test_verify}

**文件:** `giac-linalg/src/test_verify.rs`（新建）

| 函数 | 用途 | 解锁测试 |
|------|------|----------|
| `assert_matrix_equiv(a, b, ctx)` | `Expr::Matrix` 逐元素 `assert_equiv` | rref、inv、jordan、egv |
| `assert_charpoly_equiv(m, var, expected, ctx)` | `charpoly(m,var)` vs 展开特征多项式 | `symbolic.rs`、`phase3_coverage` |
| `assert_linsolve_satisfies(eqs, vars, sol, ctx)` | 解代入方程组（可复用 solve 的 `eval_at` 模式） | `phase3_coverage` linsolve |
| `assert_gramschmidt_orthogonal(basis, inner, ctx)` | 对内积 `⟨vi,vj⟩` 验零（非对角） | `gramschmidt.rs`、`phase3_coverage` |

**验收:** `cargo test-timeout -p giac-linalg`；`phase3_coverage` rref/charpoly + 模块内 smoke 已改 helper。

**已知:** `assert_linsolve_satisfies` 待 `normal` 有理和归零后接线。

**能力缺口:** 符号内积 + `integrate` 作内积定义时，依赖 calculus 通路稳定；正交验算失败先记 gap，勿回退 `contains`。

---

### H3 — ode `test_verify` {#h3-ode-test_verify}

**文件:** `giac-ode/src/test_verify.rs`（新建）

| 函数 | 用途 |
|------|------|
| `assert_ode_solution(eq, sol, dep_var, ctx)` | 对 `sol` 求导代入 `eq`（`y'`,`y''`…），lhs-rhs `assert_equiv` 为 0 |

**解锁:** `desolve.rs::desolve_harmonic`、`plugin.rs` eval 通路。

---

### H4 — simplify `assert_factorization` {#h4-simplify-factor-verify}

**位置:** `giac-simplify/src/test_verify.rs` 或 `giac-core` 测试模块

| 函数 | 用途 |
|------|------|
| `assert_factorization(orig, factors, ctx)` | 展开因子乘积 `assert_equiv` 原式 |

**解锁:** `plugin::eval_factor_*`、`eval_poly_tests::eval_partfrac_with_denom`。

---

### H5 — poly sqrt factor 返回 `Expr` {#h5-poly-sqrt-factor-api}

**问题:** `giac-poly/src/factor/sqrt.rs::quadratic_sqrt_factor_exprs` 返回 `(String, …)`，无法在 expr 层验收。

**方案（二选一）:**  
1. 改返回 `(ExprArc, usize)` / `Poly` + `product_equals`；或  
2. 测留在 poly 层用 `FactorSet::product_equals`（已有 `factor/ctx.rs`）。

**解锁:** `sqrt_factor_x2_minus_2` 1 处。

---

### H6 — poly 恒等式 helper（P1 局部）{#h6-poly-identity}

**函数:** `assert_poly_identity(lhs, rhs, var, ctx)` — `rem(lhs-rhs, var)==0` 或 `assert_equiv`。

**解锁:** `eval_poly_tests::eval_egcd_abcuv`（验 `a*u+b*v=gcd`）。

---

## 4. 改测任务（按 crate）

### C1 — giac-simplify {#c1-simplify-contains}

| 文件 | 测试 | 改法 | 层 |
|------|------|------|-----|
| `plugin.rs` | `eval_factor_via_plugin` | 删 L116 冗余 `contains`（L117 已有 golden） | C |
| `plugin.rs` | `eval_factor_x_squared_minus_two_rootof` | `matches!` + `try_as_algext_data` / **H4** | A |
| `expand.rs` | `normal_mod_power_displays_giac_style` | 完整 `assert_eq!` golden | C |
| `expand.rs` | `normal_binomial_fourth` 等 | `assert_equiv` | B |
| `trig.rs` | `texpand_*` | `assert_equiv` 对规范展开式 | B |
| `ratnormal.rs` | 4 个 weak 测试 | `assert_equiv` / `assert_eq!` | B |

**优先级:** P1 · **无阻塞** · 估 **0.5d**

---

### C2 — limit_engine {#c2-limit-engine-contains}

| 文件 | 测试 | 改法 | 层 |
|------|------|------|-----|
| `remove_lnexp.rs` | 4 个 divide/remove 测试 | Stable API + `assert_equiv`（已有 `expr_contains_ln_w`） | B |
| `mrv_w.rs` | decompose rest | `assert_equiv(rest, c)` | B |
| `mrv.rs` | choose_mrv_w | `expr_contains_ln_w` 或 `assert_equiv` 替换 display | B |
| `mrv_lead_term.rs` | rewrite_in_mrv_w | **A′** `limit(pre)==目标` 或 B 谓词 | A′/B |
| `sparse_series.rs` | `w_inv_cancel` lead | lead coeff `assert_equiv` | B |

**优先级:** P1 · **无阻塞** · 估 **0.5d**

---

### C3 — giac-core {#c3-core-eval-contains}

| 文件 | 测试 | 改法 |
|------|------|------|
| `eval.rs` | `eval_mul_mixed_complex_symbolic` | `assert_equiv` → `(1+i)*x` |
| `eval.rs` | `eval_sign_negative_and_arg_second_quadrant` | `assert_equiv` → `3π/4` 或 `assert_eq!` |
| `eval_poly_tests.rs` | gcd / content / partfrac | `assert_equiv` 或 **H4** / **H6** |
| `eval_poly_tests.rs` | `eval_egcd_abcuv` | **H6** |

**优先级:** P1 · 估 **0.5d**

---

### C4 — series {#c4-series-sparse-series}

| 文件 | 处数 | 改法 | 阻塞 |
|------|------|------|------|
| `series.rs` | 9 | **H1** `assert_series_equiv_at` + **A** `eval(Series/Taylor)` | **H1** |
| `sparse_series.rs` | 4 | 同上或 `to_expr` 后 `assert_equiv` | **H1** |

**优先级:** P2 · 估 **1d**

---

### C5 — Risch & heuristics {#c5-risch-heuristics}

| 文件 | 处数 | 改法 | 阻塞 |
|------|------|------|------|
| `risch/algebraic_rt.rs` | 6 | **H1** `assert_deriv_equals_integrand` | **H1** + T3 风险 |
| `risch/rothstein_trager.rs` | 2 | 同上 | 同上 |
| `risch/pow2expln.rs` | 2 | `assert_equiv`（无阻塞） | — |
| `integrate_heuristics.rs` | 4 | **H1** 或 CK golden `assert_equiv` | **H1** |

**优先级:** P2 · 估 **1d**（含 T3 失败登记）

---

### C6 — linalg {#c6-linalg-coverage}

| 文件 | 处数 | 改法 | 阻塞 |
|------|------|------|------|
| `tests/phase3_coverage.rs` | 10 | **H2** 全套 | **H2** |
| `src/symbolic.rs` | 2 | `assert_charpoly_equiv` | **H2** |
| `src/gramschmidt.rs` | 2 | 正交验算 | **H2** |
| `src/symbolic_eigen.rs` | 2 | `assert_matrix_equiv` | **H2** |
| `src/plugin.rs` | 1 | 删 contains 备选，保留 `assert_eq!` | — |

**优先级:** P2 · 估 **1.5d**

---

### C7 — ode plugin {#c7-ode-plugin}

| 测试 | 改法 | 阻塞 |
|------|------|------|
| `plugin::eval_desolve_via_plugin` | **H3** | **H3** |

**优先级:** P3 · **0.25d**

---

### C8 — wasm {#c8-wasm-smoke}

| 测试 | 改法 |
|------|------|
| `eval_integrate_x` / `eval_diff_x_squared` | 解析 + `assert_equiv`（P3 快改，也可提前） |

**优先级:** P3 · **0.25d**

---

### C9 — parser {#c9-parser-debug}

| 测试 | 改法 |
|------|------|
| `parse_stmt_debug_shape` | `matches!(expr, Expr::Pow(..))` 等，非 display |

**优先级:** P3 · **0.25d**（非 CAS 语义，但规范禁止子串）

---

## 5. 与父项 / 审计的衔接

| 父项 ID | 关系 |
|---------|------|
| [GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md) **3B** | 本 issue = 3B **剩余范围**（series / risch / linalg / ode / simplify 插件） |
| [GIAC-expr-api-test-audit.md](GIAC-expr-api-test-audit.md) | 每改一个 `#[test]` **必须**更新审计表一行 |
| 审计 **T1** | A 层 `assert_equiv` 须 `xcas_default()` 含 simplify — helper 内强制 |
| 审计 **T3** | `diff(integrate)` 未覆盖路径 — C5 失败登记 gap，禁止为绿而 `contains` |
| 审计 **T4** | 本 issue **H1–H4** 直接消 T4 |

**3B 状态建议:** tech-debt 表 3B 从 `partial` 改为「见 [GIAC-expr-api-test-contains-cleanup](GIAC-expr-api-test-contains-cleanup.md)」。

---

## 6. 验收标准（MVP）

- [ ] P0：**H1 + H2** 落地，各有 ≥1 个生产调用方
- [ ] P1：**C1 + C2 + C3** 语义 `contains` 清零；审计表更新
- [ ] P2：**C4 + C5 + C6** 清零（T3 失败项登记 gap + `#[ignore]` 或保留闭式 `assert_equiv`，**禁止**新 `contains`）
- [ ] P3：**C7–C9 + H5** 收尾
- [ ] `rg 'assert!.*\.contains\(' giac-rs/crates` 仅剩排除项（lexer / 已注释 C 层备选）
- [ ] `cargo test-timeout` 全绿

---

## 7. 执行建议（AFK 抓取顺序）

1. **H1** `assert_deriv_equals_integrand` + 一个 risch 测试改写（验证 T3）
2. **H2** `assert_matrix_equiv` + `phase3_coverage` rref 一条
3. **C1** simplify 全文件（快、无依赖）
4. **C2** limit_engine（快、无依赖）
5. **H1** `assert_series_equiv_at` + **C4**
6. **C5** risch/heuristics（按 T3 结果分叉）
7. **C6** linalg 剩余
8. **H3 + C7**、**H5**、**C8–C9**

---

## 8. 维护

- 新增 `contains` 语义断言：**禁止**；PR 须链本 issue 或审计表说明
- 完成子项后更新本文 §2 状态列 + [GIAC-expr-api-test-audit.md](GIAC-expr-api-test-audit.md) §1 汇总
