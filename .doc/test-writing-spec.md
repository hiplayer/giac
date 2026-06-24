# 算法 crate 单测编写规范

与 [conformance-testing.md](conformance-testing.md)（golden / check 层）互补：本文规范 **giac-rs 算法 crate 内 `#[test]`** 的分层、断言手段、**能力缺口双轨测例**与 PR 要求。表示层 API 契约见 [algorithm-expr-api.md](algorithm-expr-api.md) §2–§3。

**审计表（3A/3B 首批落点）：** [issues/GIAC-expr-api-test-audit.md](issues/GIAC-expr-api-test-audit.md)

---

## 1. 三层分类（A / B / C）

单测必须先归类再写断言；**禁止**为去掉 `contains` 而混用层级。

| 层级 | 名称 | 测什么 | 调用方式 | 断言什么 |
|------|------|--------|----------|----------|
| **A** | **主入口 / 完整通路** | 用户可见行为 | `eval` + `FuncKind::*`，Context 带对应 **plugin**（如 `xcas_default()`） | 数学结论：解代入方程、极限值、`assert_equiv` 与期望式、错误类型 |
| **A′** | **子通路 / 集成** | 多步管线、无单独 Stable API 的端到端片段 | 如 `parse → limit_preprocess_* → limit_at_*`；仍用真实 ctx | 与 **A** 相同：最终数学结果（如 `limit → -1`），**不断言**中间 AST 字符串 |
| **B** | **模块内部 / Stable API** | `*-api-stability.md` 或 `/// **Stable**` 标注函数的 I/O 契约 | **直接**调用该 API，不绕 `eval` 主入口 | 规范构造器 + `assert_equiv` / `match_*` / `decompose_*`；另加 **一种** 文档允许的 AST 漂移形态 |
| **C** | **Display 快照** | 打印顺序、括号、非语义排版 | 任意 | `assert_eq!(format_expr(...), "...")` 或有限备选字面；**不得**作为唯一语义依据 |

```text
用户 / golden check
        │
        ▼
   ┌─────────┐     子步骤组合、无独立 Stable 出口
   │ A 通路   │◄──────────────── A′ 子通路
   │ eval(*)  │
   └────┬────┘
        │ 依赖
        ▼
   ┌─────────┐
   │ B 单元   │  canonical_* / match_* / 规范构造器
   └────┬────┘
        │ 可选附带
        ▼
   ┌─────────┐
   │ C 快照   │  须 #[display_snapshot] 或测试名 *_display
   └─────────┘
```

### 1.1 与 Stable / Pipeline 的对应

| 源码 tier（[algorithm-expr-api.md](algorithm-expr-api.md) §2） | 应用层单测 |
|----------------------------------------------------------------|------------|
| **Stable** | **必须有 B**；每个函数至少 1 个 B（规范构造器 ± 漂移） |
| **Partial** | B 推荐；可用 A′ 覆盖行为边界 |
| **Pipeline / Pipeline private** | 不单独要求 B；由 **A 或 A′** 间接覆盖，或 Pipeline 步骤极难端到端时允许 **C** + issue |

---

## 2. 断言手段（按层级选用）

| 手段 | A / A′ | B | C |
|------|--------|---|---|
| `assert_equiv` / `is_zero` | ✅ 首选 | ✅ 首选 | ❌ |
| 代入原方程 / `eq_mod` / 多项式 remainder | ✅ solve / roots | 视 API | ❌ |
| `match_*` / `try_as_algext_data` / 结构 `matches!` | 仅当验输出类型 | ✅ | ❌ |
| `format_expr` 精确 golden | 仅极限等已有 check 对齐时 | 漂移形态备选 | ✅ 主手段 |
| `contains("…")` | ❌ | ❌ | ❌（C 用完整 golden） |
| 测试内私有 parser（如手搓 `eval_const_expr`） | ⚠️ 临时；须迁入 `test_verify` 并文档化 | ❌ | ❌ |

**Context：** A / A′ 使用与生产一致的 plugin 集（见各 crate `xcas_default()`）。若 `assert_equiv` 依赖 `normal`，Context **必须**安装 `giac-simplify` plugin，不得在单测里假设裸 `Context::default()` 可化简。

---

## 3. 编写规则

### 3.1 A — 主入口通路

1. 构造 `Expr::func(FuncKind::Solve | Integrate | Limit | …)`（或 parse 后 `eval`）。
2. 断言 **数学语义**，不断言中间步骤是否含某子串。
3. 每个 **公开 eval 入口**（`*-api-stability.md` 列出）至少 **1 个 A**。
4. 可放在实现模块的 `mod tests` 内，但文件名/模块名应让人看出是通路测（如 `eval_solve_via_plugin`）。

**示例（giac-solve）：** `eval(Solve(eq, var))` → 解集代入原方程（`assert_equation_solutions` / `eq_mod`）。

### 3.2 B — Stable API 单元

1. 只调用 **一个** Stable 入口（或规范构造器 + 被测函数）。
2. 期望用 **规范构造器** 或 stability doc 中的「输出形态」写出，不用偶然 `format_expr`。
3. 至少覆盖 **一种** 化简/排序导致的合法漂移（§3 algorithm-expr-api）。
4. **禁止** 用 A′ 整条管线（如 `limit_at_plus_infinity`）代替 B——管线通过 ≠ API 契约成立。

**示例（exp_diff）：** `canonical_exp_diff(input)` 与 `exp_scale_times_exp_minus_one(exp(x), 1)` 做 `assert_equiv`。

### 3.3 A′ — 子通路

用于 limit 预处理 + 求极限等 **无单一 Stable 函数** 的集成行为：

- 允许跳过 `eval(Limit)` 若当前 harness 未就绪，但必须在审计表注明。
- 断言 **最终极限/积分结果**；中间步仅用 B 或 C 单独覆盖。

### 3.4 C — Display 快照

1. 测试名后缀 `_display` 或模块注释 `// display_snapshot`。
2. 不与语义断言写在同一 `assert!`；C 失败时不应掩盖数学错误。
3. **Pipeline private** 步骤若只有 C，须在 audit 表标「缺 B」或登记 gap issue。

---

## 4. 反模式（禁止「为通过而测」）

| ❌ | ✅ |
|----|-----|
| 全文件 grep 替换 `contains` → 任意 `assert!` | 先贴 A/B/C 标签再改断言 |
| 用 `match_exp_times_exp_minus_one` 测 `scale*(exp(ε)-1)`（scale 非 `exp`） | 用 `exp_minus_one_epsilon` + 对 scale 的独立断言；或修正 API 文档/分解函数 |
| B 测试里跑完整 `limit_at_plus_infinity` | B 只测 Stable API；极限结果放到 A′ |
| A 测试里手搓常量 parser 且不文档化 | 使用 `test_verify` 共享 helper + stability doc 约定输出 |
| 直接 `integrate()` 却声称覆盖「用户通路」 | A 用 `eval(Integrate)`；直接 `integrate` 标 **B** |
| 仅 `tree_contains` / `contains` 证明预处理「做过事」 | 写清 preprocess 契约，或 A′ 断言 `limit(pre)==目标值` |
| 语义断言做不到时，把 A/B 改成弱 golden / `contains` 让 CI 绿 | **双轨**：保留 smoke + 新增 `#[ignore]` 语义测 + §8 登记（见 §6） |
| smoke 测例长期保留、与语义测重复维护 | smoke 标 `smoke-until`；阻塞解除后 **删除** smoke，只留语义测 |

---

## 5. 共享设施

| 设施 | 位置 | 用途 |
|------|------|------|
| `assert_equiv` / `is_zero` | `giac-simplify` | A / B 语义等价 |
| `test_verify` | 各 crate `src/test_verify.rs`（及 `giac-core/tests/semantic_pending.rs`） | A 层：代入方程、微分还原、矩阵/级数/ODE 等 |
| `xcas_default()` | 各 crate `plugin.rs` | A / A′ 的 Context；缺 simplify 时见 §6 / blocker 表 |
| Golden check | `giac-rs/tests/conformance/` | 系统级；**不替代** crate 内 A/B |

新增验证逻辑：**优先**放入 `test_verify`，禁止在多个测试里复制 ad-hoc 解析。

### 5.1 `test_verify` 索引（按 crate）

| Crate | 模块 | 典型 helper | 用于 |
|-------|------|-------------|------|
| `giac-simplify` | `test_verify` | `assert_factorization` | factor / partfrac 重组 |
| `giac-solve` | `test_verify` | `assert_equation_solutions`, `assert_roots_zero_poly` | solve / roots |
| `giac-calculus` | `test_verify` | `assert_deriv_equals_integrand`, `assert_series_equiv_at` | integrate / series |
| `giac-linalg` | `test_verify` | `assert_matrix_equiv`, `assert_linsolve_satisfies` | 矩阵 / linsolve |
| `giac-ode` | `test_verify` | `assert_desolve_lin_ode`, `subst_constants` | desolve 残差 |
| `giac-core` | `test_verify`（`#[cfg(test)]`） | `assert_poly_identity`, `assert_bezout`, `assert_divides` | rem / egcd（**不用** `giac_simplify::assert_equiv`，避免 dev-dep 双编译类型冲突） |
| `giac-core` | `tests/semantic_pending.rs` | 集成测 + `assert_equiv` | core eval 语义待解禁 |

谓词名含 `expr_contains_*`（实现函数，非 `format_expr.contains`）**允许**在 B 层作结构契约，**禁止**对 `format_expr` 输出做子串 `contains`。

---

## 6. 能力缺口：双轨测例（smoke-until + `#[ignore]` 语义）

当 **目标 A/B 断言**（`assert_equiv`、代入方程、残差归零、重组恒等式等）因算法或化简缺口 **暂时失败** 时，**不得**用弱 golden、`contains`、非空输出等「绕过去」冒充已验收。

采用 **双轨**（源自 [GIAC-expr-api-test-contains-cleanup](issues/GIAC-expr-api-test-contains-cleanup.md)）：

```text
目标契约（将来要绿）
  └── #[ignore] 语义测 ── 登记 B-*，CI 默认跳过
临时行为记录（修复后删）
  └── smoke-until 活跃测 ── 记录当前输出/结构，cargo test 默认跑
```

### 6.1 何时走双轨

| 情形 | 做法 |
|------|------|
| 语义断言 **可直接写绿** | 只写 A/B 测，**不要**再留 smoke |
| 语义断言 **已知会失败**（缺口已确认） | smoke-until + `#[ignore]` 语义测 + §8 登记 |
| 仅需 **C 层** display golden | 单轨 C 即可，无需 `#[ignore]` |
| 不确定能否 `assert_equiv` | 先写语义测跑一遍；失败则回退双轨，**禁止**静默改成弱断言 |

### 6.2 smoke-until（活跃 smoke，待删）

在 **现有** smoke/golden/结构测上方加注释（`rg 'smoke-until'` 可枚举）：

```rust
// smoke-until B-ODE: delete when `desolve_harmonic_satisfies_ode` green
#[test]
fn desolve_harmonic() { /* 当前行为：golden / 非空 / 结构 */ }
```

约定：

- `B-*` 与 [cleanup issue §8](issues/GIAC-expr-api-test-contains-cleanup.md#8-已知阻塞ignore-登记) 表内 ID **一致**。
- `delete when` 后写 **对侧** `#[ignore]` 语义测函数名（或多个，逗号分隔）。
- 测例内 **仅部分** 为 smoke 时（如混合测里的一个代码块），在块首标 `smoke-until`，解除后 **删块** 不必删整测。
- smoke **不是** 目标规格；禁止在 PR 里把「smoke 绿了」当作语义验收。

### 6.3 `#[ignore]` 语义测（目标契约）

另起 **独立** 测试函数，命名建议后缀 `_semantic` / `_satisfies_*` / `_recomposes`：

```rust
#[test]
#[ignore = "B-ODE: diff/normal 未将 c0,c1 当常数，ODE 残差无法归零"]
fn desolve_harmonic_satisfies_ode() {
    // 目标：assert_desolve_lin_ode(...) 或 assert_equiv(...)
}
```

约定：

- `ignore` 字符串 **必须以 `B-`、`L1-` 或审计 ID（如 `T3`）开头**，附一句缺口说明。
- 函数体写 **最终想要的** A/B 断言，不要用当前偶然输出。
- 默认 `cargo test` **不跑**；解除阻塞后：去 `ignore` → 全绿 → 删对应 smoke-until → 更新 §8 表。

本地验证：`cargo test -- --ignored <fn_name>`。

### 6.4 阻塞登记与解除

| 动作 | 要求 |
|------|------|
| 新增双轨 | 在 [GIAC-expr-api-test-contains-cleanup.md §8](issues/GIAC-expr-api-test-contains-cleanup.md#8-已知阻塞ignore-登记) 增一行：`B-*`、语义测名、smoke 待删、解除路径 |
| 解除阻塞 | ① 去 `#[ignore]` ② 语义测绿 ③ **删除** smoke-until 测例/代码块 ④ 删 §8 行 ⑤ 更新 [审计表](issues/GIAC-expr-api-test-audit.md) |
| 新增 `contains` 语义断言 | **禁止**；无例外 |

**不算语义 `contains`：** `Vec::contains`；lexer `Token` 检查；实现谓词 `expr_contains_*`；C 层完整 `assert_eq!(format_expr, "...")`。

### 6.5 与 A/B/C 的关系

| 轨道 | 相当于 | CI |
|------|--------|-----|
| 普通 A/B/C | 正式验收 | 默认跑，须绿 |
| smoke-until | 临时 C 或弱结构记录 | 默认跑，须绿，**修复后删除** |
| `#[ignore]` 语义 | 欠账的 A/B | 默认跳过，解除后升格为普通 A/B |

**禁止** 用 smoke-until 替代本该有的 B；Stable API 仍须另有（或计划另有）非 ignore 的 B 或明确登记 Pipeline 例外。

---

## 7. PR / issue 验收

算法 crate 单测 PR 须满足：

- [ ] 新增/修改测试在 [GIAC-expr-api-test-audit.md](issues/GIAC-expr-api-test-audit.md) 更新一行（或子表）
- [ ] 每个 touched **Stable** API 有 **B**（或说明由已有 B 覆盖）
- [ ] 每个 touched **公开 eval 入口** 有 **A**（或说明由已有 A 覆盖）
- [ ] **C** 已标注且非唯一语义依据
- [ ] 无未文档化的测试内 helper / parser
- [ ] 若语义断言暂不可绿：**smoke-until** + **`#[ignore]` 语义测** + [§8 阻塞表](issues/GIAC-expr-api-test-contains-cleanup.md#8-已知阻塞ignore-登记) 已登记（§6）
- [ ] 无新增 `format_expr` / `assert!` 的 **语义** `contains`
- [ ] `cargo test-timeout -p <crate>` 全绿（`#[ignore]` 除外）
- [ ] 动 eval 命令时已引用 [conformance-testing.md §3.5](conformance-testing.md#35-命令-io-契约l1-normative) 契约节；L1 双轨见 [§8](#8-l1-conformance-失败须人工确认)

---

## 8. L1 conformance 失败（须人工确认）

与 [conformance-testing.md §3.5–§3.6](conformance-testing.md#35-命令-io-契约l1-normative) 配套：crate 内 A/B 双轨（§6）管 **单元测**；本节管 **`giac_check_*` / SymPy 属性门禁**。

### 8.1 原则

| 允许 | 禁止 |
|------|------|
| 对照契约章节查因、修实现 | 为绿 CI 弱化 `sympy_verify.py` 属性判定 |
| `#[ignore = "L1-*: …"]` 挂起 **整条 L1** | 删除 failing 的 `factor_sympy_line_*` 等 |
| `smoke-until L1-*` 弱测维持默认 CI | 用 smoke 冒充 L1 已验收 |
| 人工确认后修订契约文档再改 L1 定义 | 未确认就改功能「迁就」测试 |

**动实现前：** PR 须引用该命令在 `*-api-stability.md`（或 `*-expr-api.md`）的 I/O 契约节。

### 8.2 双轨（L1）

```text
目标 L1     →  #[ignore = "L1-FACTOR-15: …"]  原 assert_factor_line_sympy / sympy 行测
临时 smoke  →  smoke-until L1-FACTOR-15: delete when `factor_sympy_line_15` green
登记        →  GIAC-expr-api-test-contains-cleanup.md §8.1
```

- **L1 ID 格式：** `L1-<命令>-<简述>` 或 `L1-<CHECK>-<行号>`（如 `L1-FACTOR-15`）。
- smoke 仅记录「当前不崩 / 有输出 / 结构」；**不得**替代 `expand(factor(p))=expand(p)` 等属性。
- 解除：去 `ignore` → L1 绿 → **删** smoke-until → 删 §8.1 行。

### 8.3 PR 自检（动到 eval 命令时）

- [ ] 已引用契约章节（conformance-testing §3.5 表）
- [ ] 若 L1 暂红：`#[ignore]` + `smoke-until` + §8.1 已登记
- [ ] 未改 `sympy_verify.py` 属性逻辑（除非契约文档已先改且人工确认）
- [ ] `cargo nextest run --release -p giac-conformance --test giac_check_<域>` 对应该域全跑

---

## 9. 参考

- [conformance-testing.md](conformance-testing.md) §3 `assert_equiv`、check golden
- [algorithm-expr-api.md](algorithm-expr-api.md) §3 稳定 API 契约模板
- [issues/GIAC-expr-api-tech-debt.md](issues/GIAC-expr-api-tech-debt.md) 3A / 3B / 3C
- [issues/GIAC-expr-api-test-audit.md](issues/GIAC-expr-api-test-audit.md) 审计表
- [issues/GIAC-expr-api-test-contains-cleanup.md](issues/GIAC-expr-api-test-contains-cleanup.md) contains 清理进度与 **§8 阻塞登记表**（双轨测例实例）
