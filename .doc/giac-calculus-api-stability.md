# giac-calculus API 稳定性分层

**规范来源：** [algorithm-expr-api.md](algorithm-expr-api.md)  
**专项：** [limit-engine-expr-api.md](limit-engine-expr-api.md)、[exp-diff-expr-api.md](exp-diff-expr-api.md)  
**代码注释约定：** 见下文 §1；各源文件 `//!` 模块头含本文件子表。

---

## 1. 源码注释格式（全 crate 统一）

| 标记 | 用于 | 可见性 |
|------|------|--------|
| `/// **Stable** — …` | 有 I/O 契约、可跨模块调用、须有单测 | `pub` / `pub(crate)` |
| `/// **Partial** — …` | 可用但范围窄或将与 MRV/主路径冲突；注释写 **退役** | `pub` / `pub(crate)` |
| `/// **Pipeline** — …` | crate 内编排步骤，非规范形保证 | 多为 `pub(crate)` |
| `// **Pipeline private** — …` | 模块内单次改写 / 检测 / 遍历 | `fn` 私有 |
| `// **Temporary** — …` | `drift_*` / `shim_*`；须写退役条件 | `fn` 私有 |

**未逐行标注的私有函数：** 默认 **Pipeline private**（见各模块 `//!` 清单）。

**提交前复审（测试全绿后）：** 见 [algorithm-expr-api.md §6.2](algorithm-expr-api.md#62-测试通过后提交--合入前复审) — 临时匹配净减少、新增 `fn` tier 已更新本文 §2–§5。

---

## 2. Crate 公开 API（`lib.rs` re-export）

| 函数 | 层级 | 模块 |
|------|------|------|
| `diff` | **Stable** | `diff` |
| `eval_diff` | **Stable** | `eval_diff` |
| `integrate` | **Stable** | `integrate` |
| `eval_integrate` | **Stable** | `eval_integrate` |
| `eval_limit` | **Stable** | `limit` |
| `eval_series` | **Stable** | `series` |
| `eval_risch` | **Partial** | `risch` |
| `hermite_reduce` | **Stable** | `risch/hermite` |
| `pow2expln` | **Stable** | `risch/pow2expln` |
| `risch_tower` / `rlvarx` | **Stable** | `risch/tower` |
| `rothstein_trager_integrate` | **Partial** | `risch/rothstein_trager` |
| `install_calculus` / `xcas_default` | **Stable** | `plugin` |

---

## 3. 子模块索引

| 子模块 | 稳定 API 文档 | 演化 issue |
|--------|---------------|------------|
| `limit_engine` | [limit-engine-expr-api.md](limit-engine-expr-api.md)、[exp-diff-expr-api.md](exp-diff-expr-api.md) | GIAC-limit-mrv-followup、GIAC-limit-exp-difference-unification |
| `integrate` / `integrate_heuristics` | 本文 §4 | GIAC-algorithm-gaps-open |
| `risch` | 本文 §5 | GIAC-algorithm-gaps-open |
| `expr_util` | 本文 §6 | GIAC-limit-mrv-followup Phase R1 |
| `diff` / `series` / `limit` | 本文 §7 | — |

---

## 4. `integrate` / `integrate_heuristics`

| 层级 | 函数 |
|------|------|
| **Stable** | `integrate`, `eval_integrate`, `integrate_frac`, `try_as_rational`, `ln_abs_expr`, `var_to_expr`, `is_var`, `is_const_wrt`, `is_exp_of_var` |
| **Partial** | `try_integrate_*` 规则表（启发式；逐步迁入 Risch/partfrac） |
| **Pipeline** | `integrate_heuristics::try_integrate_heuristic` |
| **Pipeline private** | `integrate.rs` 内其余 `is_*` / 分部 / 换元辅助 |

---

## 5. `risch`

| 层级 | 函数 |
|------|------|
| **Stable** | `hermite_reduce`, `pow2expln`, `rlvarx`, `risch_tower` |
| **Partial** | `eval_risch`, `rothstein_trager_integrate`, `try_algebraic_rt_*`, `try_integrate_x4_plus_one` |
| **Pipeline private** | `algebraic_rt.rs` / `tower.rs` 内部分解 |

---

## 6. `expr_util`

| 层级 | 函数 |
|------|------|
| **Stable** | `depends_on_var`, `is_const_wrt` |
| **演化** | `flatten_mul_*` / `signed_add_*` 计划迁入（GIAC-limit-mrv-followup Phase R1） |

---

## 7. `diff` / `series` / `limit` / `partfrac_integrate`

| 模块 | **Stable** |
|------|------------|
| `diff` | `diff`, `eval_diff` |
| `series` | `eval_series` |
| `limit` | `eval_limit` |
| `partfrac_integrate` | `integrate_one_over_quadratic`, `integrate_const_over_rational` |

---

## 8. 维护

新增 `pub` / `pub(crate)` 函数时：

1. 在源码加 `/// **Stable|Partial|Pipeline** — …`
2. 更新本文件对应 § 表
3. 若属 limit_engine / exp_diff，同步 [limit-engine-expr-api.md](limit-engine-expr-api.md) 或 [exp-diff-expr-api.md](exp-diff-expr-api.md)
