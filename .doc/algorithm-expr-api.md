# 算法模块表达式/表示层 API 规范（giac-rs 通用）

**适用范围：** 所有以 `Expr` / `Poly` / 其他表示层互转的算法 crate（见 [module-division.md](module-division.md)）。  
**专项实例：** [limit-engine-expr-api.md](limit-engine-expr-api.md)（`limit_engine` / MRV 系数）、[exp-diff-expr-api.md](exp-diff-expr-api.md)（`exp_diff` / x 层差分）。  
**配套规则：** `.cursor/rules/algorithm-expr-api.mdc`、`algorithm-before-patch.mdc`（通用；子域细节仅保留在本文档体系，无 per-module Cursor 规则）。

---

## 1. 表示层边界（防「形式漂移」）

同一数学对象可有多种 **合法但不等价于 AST 相等** 的表示。化简、有理化、预处理会改变树形，**不能**用 `format_expr` 子串或裸 `==` 做语义判断。

| 表示 | 类型（crate） | 适用场景 | 不适用 |
|------|---------------|----------|--------|
| 符号表达式树 | `ExprArc`（`giac-core`） | 含 `exp`/`ln`/参数化/代数扩域；管线编排 | 需要 guaranteed 多项式环运算时直接当 Poly 用 |
| 多项式 | `Poly`（`giac-poly`） | 已知变量上多项式 GCD、valuation、度数比 | 子树含 `exp`、`ln`、MRV 系数语义 |
| 其他 | `SparseSeries`、`AlgExt`、矩阵… | 见各子模块文档 | 跨表示混用且不经过转换 API |

**通用规则：**

1. 在 A 表示层做语义判断，就留在 A 层完成；**禁止**「先 `expr_to_poly` 再认形状」除非上游算法明确要求且变量域无超越函数。
2. 调用 `ratnormal` / `normal` / `expand` 之后，**禁止**假定 AST 与调用前 `==`；须走该域的 **规范入口**（`canonical_*`）或 **分解 API**（`decompose_*`）。
3. **禁止**在调用方复制 owning 模块里的形状表；缺能力回 owning 模块补稳定 API（见 [algorithm-before-patch](../.cursor/rules/algorithm-before-patch.mdc)）。

---

## 2. 三层函数分类（全 crate 统一）

```text
┌──────────────────────────────────────────────────────────┐
│ 稳定 API — pub(crate) 或公开 API，有 I/O 契约，须有单测   │
│  canonical_* / decompose_* / 规范构造器 / 域内谓词        │
├──────────────────────────────────────────────────────────┤
│ 临时 API — 私有，仅规范入口内部；必须注明退役条件         │
│  drift_*（AST 漂移）  shim_*（跨模块绕行）  legacy_*     │
├──────────────────────────────────────────────────────────┤
│ 管线私有 — 私有，编排步骤；禁止它模块复制或导出谓词       │
│  try_* / peel_* / 主循环 / 单次改写                      │
└──────────────────────────────────────────────────────────┘
```

| 类别 | 可见性 | 命名约定 | 模块外 |
|------|--------|----------|--------|
| **稳定 API** | `pub(crate)` 或 `pub` | `canonical_<域>_*`、`decompose_*`、`*_expr`（原子）、`expr_contains_*` | ✅ 允许 |
| **临时 API** | `fn` 私有 | **`drift_*`**：化简导致的 AST 漂移吸收 | ❌ 禁止 |
| **临时 API** | `fn` 私有 | **`shim_*`**：短期绕行他模块缺陷 | ❌ 禁止；须写替换 issue |
| **管线私有** | `fn` 私有 | `try_*`、`peel_*`、`*_loop` 等 | ❌ 禁止复制 |

**退役：** 每个 `drift_*` / `shim_*` 须在模块头或 `.doc/issues/` 写明删除条件（如「Phase 3A simplify 落地后删除 drift」）。

---

## 3. 稳定 API 契约模板（新增/审查时用）

每个稳定函数应能填下表（可写在模块 `//!` 或 `.doc/` 子页）：

| 字段 | 说明 |
|------|------|
| **函数** | 名称 |
| **输入类型** | 如 `&ExprArc`；是否要求已 canonical |
| **输出类型** | 如 `ExprArc` / `MrvCoeffParts` / `Result<Poly, _>` |
| **输出形态** | 保证的原子/范式（或明确「不保证最简」） |
| **是否调用规范入口** | 内部是否先 `canonical_*` |
| **禁止** | 如「禁止输入含 ln(w) 未剥离」 |

**单测要求：** 规范构造器往返 + 至少一种漂移/化简后形态；避免仅 `format_expr.contains(...)`。

---

## 4. 各 crate 规范入口（索引）

| Crate / 子模块 | 主表示 | 规范 / 边界 API | 专项文档 |
|----------------|--------|-----------------|----------|
| `giac-simplify` | `ExprArc` | `normal`, `ratnormal`, `expand`（各有适用范围） | [giac-simplify-api-stability.md](giac-simplify-api-stability.md) |
| `giac-poly` | `Poly` | `factor_into`, `partfrac_rational_terms`, `Poly` 环运算；`expr_to_poly` 在 giac-core | [giac-poly-api-stability.md](giac-poly-api-stability.md) |
| `giac-calculus` / `limit_engine` | `ExprArc` + `SparseSeries` | `mrv_w::canonical_mrv_coeff`, `decompose_mrv_coeff`, `remove_lnexp` | [limit-engine-expr-api.md](limit-engine-expr-api.md) |
| `giac-calculus`（全 crate） | `ExprArc` | `integrate`, `eval_limit`, `depends_on_var`；源码 `/// **Stable**` 标记 | [giac-calculus-api-stability.md](giac-calculus-api-stability.md) |
| `giac-calculus` / `exp_diff` | `ExprArc`（x 层 `exp` 差分） | `canonical_exp_diff`, `match_exp_times_exp_minus_one`, `is_exp_minus_one_factor` | [exp-diff-expr-api.md](exp-diff-expr-api.md) |
| `giac-calculus` / `risch` | `ExprArc` | 积分塔、`transcendental` 边界 | — |
| `giac-calculus` / `intg` | `ExprArc` | 规则表入口 `_integrate` | — |
| `giac-solve` | `ExprArc` / `Poly` | `solve`、Sturm 多项式前提 | — |
| `giac-groebner` | `Poly` | Gröbner 基、变量序 | — |

未单独成文的子模块：新增 **稳定** 表示层 API 时，在本表补一行或增加 `.doc/<module>-expr-api.md` 链接。

---

## 5. 反模式（全 crate）

| ❌ | ✅ |
|----|-----|
| `format_expr(e).contains("ln")` 判语义 | 域内 `expr_contains_*` / `decompose_*` |
| 对含 `exp` 的式子 `expr_to_poly` 判 lead | 走级数 / MRV / 上游 transcendental 路径 |
| 在 B 模块复制 A 模块的 `drift_*` 逻辑 | 调用 A 的 `canonical_*` 或修 `giac-simplify` |
| 第三个 `if looks_like_ck_*` 分支 | 合并为一条算法或 owning 模块 API |
| 临时 `shim` 无 issue / 无退役注释 | 注释 + issue 链接 |

---

## 6. 新增函数检查清单

### 6.1 新增时（动刀前 / 动刀中）

1. **归类：** 稳定 / `drift_` / `shim_` / 管线私有？  
2. **表示层：** 全程 `ExprArc` 还是 `Poly`？若跨表示，转换点是否唯一？  
3. **I/O 契约：** 填 §3 表格。  
4. **漂移：** 新 AST 漂移形 → 只加 `drift_*` + `canonical_*` 单测；语义 → `decompose_*` 或 owning 算法。  
5. **文档：** 更新本表或子模块 doc；临时 API 写退役条件。  
6. **与 algorithm-before-patch 一致：** 不为单测在调用方打补丁。

### 6.2 测试通过后（提交 / 合入前复审）

**顺序：** 实现完成 → `cargo test-timeout`（及域内相关单测）**全绿** → **再 review 一遍 diff**（禁止测试一绿即停）。

| 审查项 | 做法 |
|--------|------|
| **临时匹配是否减少** | 对照 diff：`if looks_like`、per-case 形状表、重复 `try_*` 兜底、调用方 `shim` 是否 **净减少或合并**；若净增，须说明为何不能算法化（链到 gap issue） |
| **可否删旧临时** | 新主路径已覆盖的旧 `drift_*` / `shim_*` / 重复 eval 提升 → **同 PR 删除或开 issue 并注释退役** |
| **新增 fn 稳定性** | 每个 **新增或签名变更** 的 `fn` 标注 tier 并写入 crate 稳定性 doc 表 |
| **注释与 inventory** | 源码 `/// **Stable**` / `// **Temporary**` / `// **Pipeline private**`；算法 crate 跑 `python3 scripts/annotate_api_tiers.py --inventory` 刷新 Per-file 表 |
| **静默失败** | `try_*` 失败不得假成功（如 factor 静默 `Ok(vec![g])` 当不可约）；须走完整 fallback 链或明确边界 |

**新增函数 tier 判定（摘要）：**

| Tier | 条件 |
|------|------|
| **Stable** / **Stable (bounded)** | 有 I/O 契约、可跨模块调用、有单测；`pub` 或 intentional `pub(crate)` |
| **Partial** | 启发式 / 规则表子集；失败 `None`/`Err`；注释写扩展或退役计划 |
| **Temporary** | `shim_*`、`drift_*`、AlgExt stub；**必须**写退役 issue |
| **Pipeline private** | 模块内编排；禁止它模块复制 |

各 crate 细则：[giac-simplify-api-stability.md](giac-simplify-api-stability.md)、[giac-poly-api-stability.md](giac-poly-api-stability.md)、[giac-calculus-api-stability.md](giac-calculus-api-stability.md)。

**PR 描述建议：** 列出「删除 / 合并的临时匹配」与「新增 fn + tier」两行摘要。

---

## 7. 参考

- C++ 模块分工：[module-division.md](module-division.md)  
- 算法缺口：[issues/GIAC-algorithm-gaps-open.md](issues/GIAC-algorithm-gaps-open.md)  
- limit 专项：[limit-engine-expr-api.md](limit-engine-expr-api.md)
