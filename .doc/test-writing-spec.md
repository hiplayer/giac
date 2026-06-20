# 算法 crate 单测编写规范

与 [conformance-testing.md](conformance-testing.md)（golden / check 层）互补：本文规范 **giac-rs 算法 crate 内 `#[test]`** 的分层、断言手段与 PR 要求。表示层 API 契约见 [algorithm-expr-api.md](algorithm-expr-api.md) §2–§3。

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

---

## 5. 共享设施

| 设施 | 位置 | 用途 |
|------|------|------|
| `assert_equiv` / `is_zero` | `giac-simplify` | A / B 语义等价 |
| `test_verify` | 各 crate `src/test_verify.rs`（测试专用） | A 层：代入方程、rootof 列表、扩域谓词 |
| `xcas_default()` | 各 crate `plugin.rs` | A / A′ 的 Context；缺 simplify 时见 audit blocker |
| Golden check | `giac-rs/tests/conformance/` | 系统级；**不替代** crate 内 A/B |

新增验证逻辑：**优先**放入 `test_verify`，禁止在多个测试里复制 ad-hoc 解析。

---

## 6. PR / issue 验收（替代「contains 减少 80%」）

算法 crate 单测 PR 须满足：

- [ ] 新增/修改测试在 [GIAC-expr-api-test-audit.md](issues/GIAC-expr-api-test-audit.md) 更新一行（或子表）
- [ ] 每个 touched **Stable** API 有 **B**（或说明由已有 B 覆盖）
- [ ] 每个 touched **公开 eval 入口** 有 **A**（或说明由已有 A 覆盖）
- [ ] **C** 已标注且非唯一语义依据
- [ ] 无未文档化的测试内 helper /parser
- [ ] `cargo test-timeout -p <crate>` 全绿

---

## 7. 参考

- [conformance-testing.md](conformance-testing.md) §3 `assert_equiv`、check golden
- [algorithm-expr-api.md](algorithm-expr-api.md) §3 稳定 API 契约模板
- [issues/GIAC-expr-api-tech-debt.md](issues/GIAC-expr-api-tech-debt.md) 3A / 3B / 3C
- [issues/GIAC-expr-api-test-audit.md](issues/GIAC-expr-api-test-audit.md) 审计表
