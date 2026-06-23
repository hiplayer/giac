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
| 有理多项式 | `Poly`（`giac-poly`，系数 **ℚ**） | `expr_to_poly` → GCD、valuation、factor（ℚ） | 系数含 `AlgExt`/`rootof` |
| 代数系数多项式 | `Poly<AlgExtC>`（规划，P1） | `poly_alg_from_expr` → partfrac/factor/roots over K | P1 未实现前勿假设存在 |
| 其他 | `SparseSeries`、`AlgExt`、矩阵… | 见各子模块文档 | 跨表示混用且不经过转换 API |

**通用规则：**

1. 在 A 表示层做语义判断，就留在 A 层完成；**禁止**「先 `expr_to_poly` 再认形状」除非上游算法明确要求且变量域无超越函数。
2. 调用 `ratnormal` / `normal` / `expand` 之后，**禁止**假定 AST 与调用前 `==`；须走该域的 **规范入口**（`canonical_*`）或 **分解 API**（`decompose_*`）。
3. **禁止**在调用方复制 owning 模块里的形状表；缺能力回 owning 模块补稳定 API（见 [algorithm-before-patch](../.cursor/rules/algorithm-before-patch.mdc)）。
4. 树状表示**不承载**环/域/canonical 等算法上下文时，按 [§6](#6-算法上下文erased-tree--explicit-context) 显式携带；禁止裸树做语义判断。

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

**单测要求：** 规范构造器往返 + 至少一种漂移/化简后形态；避免仅 `format_expr.contains(...)`。**分层与 PR 规则：** [test-writing-spec.md](test-writing-spec.md)；**审计表：** [issues/GIAC-expr-api-test-audit.md](issues/GIAC-expr-api-test-audit.md)。

---

## 4. 各 crate 规范入口（索引）

| Crate / 子模块 | 主表示 | 规范 / 边界 API | 专项文档 |
|----------------|--------|-----------------|----------|
| `giac-simplify` | `ExprArc` | `normal`, `ratnormal`, `expand`（各有适用范围） | [giac-simplify-api-stability.md](giac-simplify-api-stability.md) |
| `giac-poly` | `Poly` | `factor_into`, `partfrac_rational_terms`, `Poly` 环运算；Expr 边界见下 | [giac-poly-api-stability.md](giac-poly-api-stability.md) |
| `giac-core` / `algebra::poly` | `Poly`（ℚ）/ `PolyAlgExt` / `Expr` | **`expr_to_poly`** / **`poly_to_expr`** / **`poly_alg_from_expr`** / **`algext_poly_to_expr`** / **`expr_contains_alg_coeff`** / **`poly_algext_roots`** | [giac-core-algebra-api-stability.md](giac-core-algebra-api-stability.md)、[expr-poly-conversion.md](expr-poly-conversion.md) |
| `giac-poly` | `Poly<C>` / `PolyCoeff` | 系数环泛化；ℚ 默认 `Poly`；`ring_*` 构造非 ℚ 多项式 | [giac-poly-api-stability.md](giac-poly-api-stability.md) |
| `giac-calculus` / `limit_engine` | `ExprArc` + `SparseSeries` | `mrv_w::canonical_mrv_coeff`, `decompose_mrv_coeff`, `remove_lnexp` | [limit-engine-expr-api.md](limit-engine-expr-api.md) |
| `giac-calculus`（全 crate） | `ExprArc` | `integrate`, `eval_limit`, `depends_on_var`；源码 `/// **Stable**` 标记 | [giac-calculus-api-stability.md](giac-calculus-api-stability.md) |
| `giac-calculus` / `exp_diff` | `ExprArc`（x 层 `exp` 差分） | `canonical_exp_diff`, `match_exp_times_exp_minus_one`, `is_exp_minus_one_factor` | [exp-diff-expr-api.md](exp-diff-expr-api.md) |
| `giac-calculus` / `risch` | `ExprArc` | 积分塔、`transcendental` 边界 | — |
| `giac-calculus` / `intg` | `ExprArc` | 规则表入口 `_integrate` | — |
| `giac-solve` | `ExprArc` / `Poly` | `eval_solve`, `eval_froot`, `eval_realroot`, `quadratic_rootof_roots` | [giac-solve-api-stability.md](giac-solve-api-stability.md) |
| `giac-ode` | `ExprArc` | `eval_desolve` | [giac-ode-api-stability.md](giac-ode-api-stability.md) |
| `giac-groebner` | `Poly` | `greduce`, `greduce_mod`（变量序显式） | [giac-groebner-api-stability.md](giac-groebner-api-stability.md) |

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
| 裸树做语义判断（`==`、形状匹配、`div_rem` 混用） | 显式上下文 + 域内 `canonical_*` / 专用类型 |
| `Poly::div_rem` 验 ℚ[others][main] 整除 | `UnivariateIn::divides` / `quo_exact_wrt` |
| `Poly<AlgExtCPolyCoeff>` 裸树 + 隐式 `align_coeff` | **FieldSession**（ambient K + working L） |
| `PolyCoeff::coeff_one()` 当 K 上单位元 | `session.one()` / `ring_one(&K)` |
| sqrt/cbrt 后仍 `ring_int(ambient, 2)` | `session.bump()` 后在 **L** 上造常数 |

---

## 6. 算法上下文（erased tree + explicit context）

**通用原则：** giac-rs 的树状表示层（`ExprArc`、`Poly`、`SparseSeries`、`AlgExt`…）通常是 **intentionally erased**——树只编码**数据形状**（项、子树、系数值），**不保证**承载算法所需的全部语义环境。**相关算法必须显式携带上下文**；这不是 `Poly` 或扩域专有风格，而是全表示层的默认设计模式。

| 表示 | 树里有什么 | 树里通常没有 | 上下文载体（示例） |
|------|-----------|-------------|-------------------|
| `Poly` | 项、指数、系数 | 环、主元、管线阶段 | `MainVar` / `UnivariateIn` / `SqffRingCtx` |
| `Poly<AlgExtC>` | 稀疏项 + 系数 | ambient **K**、working **L** | `FieldSession` |
| `ExprArc` | AST 节点 | 域内 canonical 形、求值绑定 | `canonical_*` / `decompose_*`、`Context` |
| `SparseSeries` | Laurent 项 + 系数子树 | `ln(w)` 等 MRV 语义 | `canonical_mrv_coeff` |
| `AlgExt` | 坐标、生成元 | 与当前 working field 对齐 | tower / `common` / `align` |

### 6.1 何时需要显式上下文

| 情形 | 做法 |
|------|------|
| **类型已足够特化** | 树 + 专用类型即可（如 `FlatUni` 上 `div_rem` 语义唯一） |
| **入口范式固定、ambient 不变** | 规范入口 + I/O 契约约束（如 `expr_to_poly` 后全程 ℚ[x]） |
| **同形多义 或 ambient 会变** | **必须**显式上下文；禁止假设「光看树就知道在哪个环/域」 |

**判据（动刀前自问）：**

1. 去掉上下文后，同一棵树是否有两种合法运算语义？
2. ambient（环、域、变量序、塔顶）是否会在算法推进中变化？

任一为是 → 必须显式携带；禁止在子步骤散落隐式对齐或默认假设。

### 6.2 三种上下文载体

不必统一成一种 `Session` 类型；按域选用：

| 形态 | 何时用 | 特征 | 示例 |
|------|--------|------|------|
| **标签 / 包装类型** | 入口固定、ambient 不变 | 不可变、编码「在哪个环/主元」 | `MainVar`、`TnEmbed`、`PolyInK` |
| **规范 API** | AST 漂移、域内原子形 | `canonical_*` 收敛形态后再比较/分解 | `canonical_mrv_coeff`、`canonical_exp_diff` |
| **可变 session** | 计算中扩大分裂域/塔 | ambient 固定 + working 单调扩大 | `FieldSession`（§6.3） |

**管线排错顺序（`Poly` / factor 等）：** 先答「在哪个环？主元是谁？哪一阶段？」，再问算术为何失败。详见 [giac-poly-nested-ring-types.md](issues/GIAC-poly-nested-ring-types.md)。

### 6.3 实例：扩域求根 `FieldSession`（`Poly<C>` over K）

`Poly<AlgExtCPolyCoeff>` 的稀疏树只保证「有哪些项、系数是哪些 `AlgExtC`」；**不保证**所有系数在同一 `ExtensionField`，更不保证与 Cardano/resolvent 步骤中的**当前分裂域 L** 一致。这是 §6「可变 session」在扩域上的特化；P3-6 实践表明：公式正确仍会因域错位而算错或 `align` 爆炸。

**规则（算法设计，非可选风格）：**

1. **Ambient K** — 多项式语义所在的系数域；入口由 `infer_field` + `normalize_coeffs` 固定，之后不变。
2. **Working L** — 求根/分解过程中单调扩大的分裂域；每次 `sqrt` / `cbrt` / T3+ `adjoin` 后 **必须** `L ← max(K, adjoin…)`，后续 `ring_int`、`ring_half`、中间结果都在 **L** 上构造。
3. **类型承载** — 子算法签名带 `&FieldSession`（或等价 `WorkingField`），**禁止**只传 `&Poly` + 裸 `&Arc<ExtensionField>` 并假设仍是 K。
4. **对齐内聚** — `add/mul/div` 仅接受 `session.align(a,b)` 后的对，或 session 上的方法；禁止在 Cardano/resolvent 各层手写 `align_coeff` 链。
5. **与表示层分离** — `poly_alg_from_expr` 可产生 ℚ/K 混合系数；**算法入口** normalize 到 K 后再进 session；`verify` 与算法共用同一 normalize+monic 路径。

**最小 API（已实现，见 [giac-core-algebra-api-stability.md](giac-core-algebra-api-stability.md) `field_session.rs` / `PolyInK`）：**

```text
FieldSession { ambient: K, working: L }
  .int / .half / .zero / .one          // 总在 working 上
  .lift / .align / .add / .mul / .div
  .adjoin_sqrt / .adjoin_cbrt
  .adjoin_primitive_cube_root_of_unity  // 纯三次 ω 分支（PR-C′）
```

**`build_resolvent_cubic`（PR-D′）：** 对 depressed 四次 \(y^4+py^2+qy+r\)，
\(R(z)=z^3-pz^2-4rz+(4pr-q^2)\)；golden：`t^4+t+1` → \((p,q,r)=(0,1,1)\) → \(z^3-4z-1\)。

**Eval / solve 接线（R5，[GIAC-ext-registry-removal-plan](issues/GIAC-ext-registry-removal-plan.md)）：**

- `Context` 内嵌 `Rc<FieldSession>`；`Clone` 共享 extension cache。独立 `Context::new()` → 独立 cache。静态 `ExtensionField::*` 无 session 时不跨调用 dedup（R5b ephemeral）。
- `poly_algext_roots_for_ctx`：`fork_ambient(K)` 保留 roots 的 K/L 语义，复用 `ctx` 的 `common_cache` / adjoin dedup。
- `eval` 中 `fold_algext_*_for_ctx`：跨域 AlgExt 合并走 `ctx.session()`，避免与 solve 分裂 cache。
- 静态 `ExtensionField::common_over_q` 在无 Context 时仍可用（同一 per-thread session）。

**适用：** `poly_algext_roots`、resolvent cubic、Cardano、双二次 split、将来 `Poly<AlgExtC>::factor/gcd` over K。  
**索引：** [GIAC-poly-roots-field-session-plan.md](issues/GIAC-poly-roots-field-session-plan.md)；[GIAC-poly-p3-6-quartic-roots-gaps.md](issues/GIAC-poly-p3-6-quartic-roots-gaps.md) G5；[expr-poly-conversion.md](expr-poly-conversion.md) path B。

---

## 7. 新增函数检查清单

### 7.1 新增时（动刀前 / 动刀中）

1. **归类：** 稳定 / `drift_` / `shim_` / 管线私有？  
2. **表示层：** 全程 `ExprArc` 还是 `Poly`？若跨表示，转换点是否唯一？  
3. **算法上下文：** 按 §6 判据——树是否同形多义？ambient 是否会变？选用标签 / `canonical_*` / session 中哪一种？  
4. **I/O 契约：** 填 §3 表格。  
5. **漂移：** 新 AST 漂移形 → 只加 `drift_*` + `canonical_*` 单测；语义 → `decompose_*` 或 owning 算法。  
6. **文档：** 更新本表或子模块 doc；临时 API 写退役条件。  
7. **与 algorithm-before-patch 一致：** 不为单测在调用方打补丁。

### 7.2 测试通过后（提交 / 合入前复审）

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

## 8. 参考

- C++ 模块分工：[module-division.md](module-division.md)  
- 算法缺口：[issues/GIAC-algorithm-gaps-open.md](issues/GIAC-algorithm-gaps-open.md)  
- 表示层 / API / 测试技术债：[issues/GIAC-expr-api-tech-debt.md](issues/GIAC-expr-api-tech-debt.md)  
- limit 专项：[limit-engine-expr-api.md](limit-engine-expr-api.md)
