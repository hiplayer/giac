# giac-poly API 稳定性分层

**规范来源:** [algorithm-expr-api.md](algorithm-expr-api.md)  
**上游缺口:** [issues/GIAC-simplify-poly-upstream-gaps.md](issues/GIAC-simplify-poly-upstream-gaps.md) §2  
**Expr ↔ Poly 边界:** [expr-poly-conversion.md](expr-poly-conversion.md)（normative）。P1 表示层：[giac-poly-p1-representation.md](giac-poly-p1-representation.md)。**ℚ 路径：** `expr_to_poly`；**代数系数：** `poly_alg_from_expr` → `PolyAlgExt`。本 crate 算法管线仍全程 `Poly`（ℚ）；`Poly<C>` 骨架见 P1。

---

## 1. 源码注释格式

| 标记 | 用于 | 可见性 |
|------|------|--------|
| `/// **Stable** — …` | 环运算、分解、partfrac 主入口 | `pub` |
| `/// **Stable (bounded)** — …` | 输入次数/变元数受限 | `pub` |
| `/// **Stable (crate-internal)** — …` | 嵌套环精确除法、content/pp 实现 | `pub(crate)` |
| `/// **Partial** — …` | 启发式/模式表；失败返回 `None` 或 Err | `pub` |
| `pub(crate) fn try_*` | factor 管线步骤 | crate 内 |
| `fn try_*` | 私有 fallback / 形状检测 | 模块内 |

**除法语义（禁止混用）：**

| API | 环 | 语义 |
|-----|-----|------|
| `Poly::div_rem` | 多元展示环 | leading-monomial 除法 |
| `subresultant::nested_div_rem_wrt_in` | **ℚ[others][var]** nested | leading-term 商；`deg(d)=0` 通常不 closed |
| `subresultant::nested_exact_quo_wrt_in` | **ℚ[others][var]** | 上一元除法，余式非零 → `Err` |
| `univ_wrt::univariate_div_rem_wrt` | **K[var]** flat，K 域 | Euclidean；`deg r < deg b`；`d=0`/`lc=0` → `Err`（**crate-internal**） |
| `univ_wrt::{gcd_wrt,egcd_wrt,square_free_part_wrt}` | **K[var]** flat | 要求 `C: FieldCoeff`；**经 [`FlatUni`] 方法**（L1-3 起不再 crate 根 re-export） |
| `nested::FlatUni::div_rem` | **K[var]** flat | 绑定 `MainVar`；内部调 `univ_wrt` |

相关 issue：[GIAC-poly-flat-field-division-layering](issues/GIAC-poly-flat-field-division-layering.md)、[GIAC-poly-nested-ring-types](issues/GIAC-poly-nested-ring-types.md)。

**命名:** 公开 `try_*`（如 `try_hensel_lift_bivariate`）表示 **可选算法路径**，非 [algorithm-expr-api](algorithm-expr-api.md) 意义的临时 `shim_*`；失败时静默 `None`，调用方须处理。

**提交前复审（测试全绿后）：** 见 [algorithm-expr-api.md §7.2](algorithm-expr-api.md#72-测试通过后提交--合入前复审) — 检查 factor 等模块临时 fallback 是否净减少、新增 `fn` tier 是否已更新本文 Per-file 表。

---

## 2. Crate 公开 API（`lib.rs` re-export）

### 2.1 核心类型与环运算 — **Stable**

| 符号 | 模块 | 说明 |
|------|------|------|
| `Poly`, `PolyQ`, `Poly<C>` | `poly` | 稀疏多项式；默认 `Poly` = ℚ |
| `PolyCoeff`, `FieldCoeff` | `poly_coeff` | 系数环 / 系数域 trait（flat gcd 需后者） |
| `Monomial`, `Var` | `monomial` | 指数向量（与 `C` 无关，P1-5） |
| `UnivariateIn<C>`, `UnivariatePoly<C>`, `FlatUni<C>` | `nested` | 单变量视图；`FlatUni` = flat K[var]（`div_rem` 经 `univ_wrt`） |
| `PolyError`, `PolyResult` | `error` | |
| `quo`, `rem`, `egcd`, `simp2`, `abcuv` | `poly` | 精确除法；失败 → Err |
| `content`, `gauss` | `ops` | |
| `ModInt`, `smod`, `irem` | `modint` | |
| `PolyMod`, `modp` | `modular` | |

### 2.2 一元 / 结果式 / Sturm — **Stable**

| 符号 | 模块 | 说明 |
|------|------|------|
| `resultant`, `coeff_at`, `univariate_degree`, `roots` | `resultant` | `roots` 仅低次精确根 |
| `univariate_derivative`, `square_free_factorization`, `square_free_part` | `univariate` | |
| `substitute_univariate`, `odd_multiplicity_part` | `univariate` | |
| `eval_univariate_at`, `sign_variations` | `univariate` | |
| `sturm_sequence`, `sturm_sign_variations_at`, `sturmab_count` | `univariate` | giac-solve 用 |
| `chinrem`, `chinrem_lists` | `chinrem` | |

### 2.2.1 Flat K[var]（`univ_wrt` + `FlatUni`）— **Stable**

| 符号 | 环 | 说明 |
|------|-----|------|
| `FlatUni::{try_new,div_rem,gcd,sqff,…}` | **K[var]** | **首选** flat 一元入口（`C: FieldCoeff`） |
| `is_univariate_in`, `scalar_coeff_wrt` | — | 一元检测 / 系数抽取（crate 根 re-export） |
| `derivative_wrt`, `content_scalars`, `quadratic_coeffs_wrt` | **K[var]** | 辅助（crate 根 re-export） |
| `univ_wrt::{gcd_wrt,div_rem,…}` | **K[var]** | **crate-internal**；勿裸调，用 `FlatUni` |

### 2.3 因式分解 — **Stable (bounded)**

| 符号 | 模块 | 边界 / 缺口 |
|------|------|-------------|
| `factor_into` | `factor` | 失败 → `None`；FAC-G1–G3 外形状 |
| `factor_poly` | `factor` | 展示用乘积；非 guaranteed 不可约列表 |
| `factor_into_by_rational_roots` | `factor` | 一元有理根链 |
| `factor_poly_mod`, `factor_mod_irreducibles` | `factor` | 模 p |
| `as_perfect_power`, `try_linear_power` | `factor/power` | 幂次检测 |
| `quadratic_sqrt_factor_exprs` | `factor/sqrt` | 二次 sqrt 形（Expr 侧配合） |
| `factor_power_pairs` | `factor/univariate` | 带重数的因子对 |
| `vars_in`, `ratio_perfect_sqrt` | `factor/util` | |

#### 2.3.1 一元 `Ok(vec![g])` 单因子 fallback（4B / §7.2）

`factor_square_free` / `factor_quadratic` 等在算法未分裂时返回 **单元素列表 `[g]`**。数学上恒为合法因子分解（乘积 = 输入）；调用方须知此为 **不可约见证** 或 **低次正确因子**，非「假失败」。

| 位置 | 条件 | 语义 | 层级 |
|------|------|------|------|
| `factor_square_free` | `d == 0` | 常数项（无 `var` 次数）→ 单因子自身 | **Bounded OK** |
| `factor_square_free` | `d == 1` | 一次多项式已是不可约因子 | **Bounded OK** |
| `factor_square_free` | 末尾 fallback | Zassenhaus / 双二次 / 双三次 / 完美幂均失败后，**square-free** `g` 在 ℚ 上视为不可约单块 | **Partial** — 见 FAC-G*；乘积仍 = `g` |
| `factor_quadratic` | `quadratic_abc` 失败 | 非标准二次形 → 原式单块 | **Bounded OK** |
| `factor_quadratic` | Δ 非 ℚ 平方 | 二次在 ℚ 不可约 → `[p]` | **Bounded OK**（如 `x²+1`） |
| `factor_power_pairs_core` | `rest` 次数 ≤ 2 | 余式为一次/二次 → `(rest, 1)` | **Bounded OK** |
| `factor_power_pairs_core` | 无有理根且次数 > 2 | **`Err(NotImplemented)`** — 非静默 | — |
| `cantor_zassenhaus_block` (fpx) | `k == i` | DDF 块次数 = 目标不可约度 → 块本身 | **Bounded OK**（𝔽_p） |

**禁止误解：** `[g]` 仅当 `g` square-free 且算法链已穷尽 **已实现的** 分裂手段；不等于「giac C++ 完整 factor」对高次多项式的保证。

**单测：** `factor/univariate.rs` `tests::factor_univariate_*` — 可约二次分裂 vs 不可约二次单块 vs 乘积还原。

### 2.4 部分分式 — **Stable (bounded)**

| 符号 | 模块 | 边界 |
|------|------|------|
| `partfrac_terms` | `partfrac` | 依赖 `factor_into` |
| `partfrac_rational_terms` | `partfrac` | 非线性因子 / 重复二次 → `NotImplemented` |

### 2.5 参数结果式（Risch / RT）— **Stable**

| 符号 | 模块 | 说明 |
|------|------|------|
| `tresultant_eliminate_x`, `num_minus_t_derivative` | `tresultant` | |
| `eval_param_poly`, `rational_roots_in_t` | `tresultant` | |
| `biquadratic_res_conjugate_pairs`, `biquartic_conjugate_pairs` | `tresultant` | |
| `AlgebraicRt`, `ConjugatePair` | `tresultant` | |

### 2.6 嵌套环表示层 — **Stable**

| 符号 | 模块 | 说明 |
|------|------|------|
| `MainVar`, `UnivariateIn`, `UnivariateOver` | `nested` | ℚ[others][main] 视图；`.divides` / `.exact_quo_dividing` |
| `UnivariatePoly` | `nested` | 拥有的 ℚ[others][main] 元素 |
| `CoeffRingPoly` | `nested` | ℚ[others] 系数环；`.gcd` / `.exact_quo` |
| `TnEmbed`, `BivariateEmbed` | `nested` | sparse_bi `eval_tn`；嵌入像带 `(main,t,n,aux)` |
| `div_rem_wrt_aux_indep` | `nested` | Hensel：除子在 aux 上常系数时的 `(q,r)` |
| `eval_aux` | `nested` | 在 ℚ[others][main] 中对 aux 赋值，main 不变 |
| `GoodEval` | `factor/eval` | 好点赋值，保持 `preserved_main_degree` |
| `SqffRingCtx`, `FactorSet` | `factor/ctx` | sqff 因子链上下文（crate-internal） |
| `PolyFactorTower`, `CoeffRing` | `factor/tower` | FAC-G2；`factor_sqff_chain` = aux-lift → good_eval → sparse_bi |
| `factor_bivariate_flat` | `factor/sparse` | ℚ[main,aux] 二元分解，无嵌套 `factor_multivariate_rec` |
| `HenselPair` | `factor/ctx` | 二元 Hensel @ aux=0（crate-internal） |
| `EmbedFactorDraft` | `nested` | sparse_bi 重建 IR（crate-internal） |
| `FlatUni` | `nested` | **K[var]** flat；`C: PolyCoeff`（L1→`FieldCoeff`）；经典 Euclidean `div_rem` |
| `MultivariatePoly` | `nested` | 多元展示环边界；显式 leading-monomial `div_rem` |
| `PrimitivePart` | `nested` | `primitive_part_wrt` 结果 tagged |
| `DilationMap` | `nested` | sparse_bi dilation `{ aux_a, aux_b }` + apply/undo（crate-internal） |

**禁止：** 嵌套环热路径用 `Poly::div_rem` 验整除；见 [GIAC-poly-nested-ring-types](issues/GIAC-poly-nested-ring-types.md) Phase 1。

---

## 3. 公开但未 `lib.rs` re-export 的 API

以下 `pub fn` 存在，供 crate 内或将来导出；**跨 crate 请优先 §2 符号**。

| 符号 | 模块 | 层级 |
|------|------|------|
| `factor_multivariate`, `factor_into_poly` | `factor/multivariate` | **Stable (bounded)** |
| `gcd_univariate` | `univariate` | **Stable** |
| `subresultant_gcd` | `subresultant` | **Stable** |
| `quo_exact_wrt`, `quo_exact_coeff`, `univariate_div_rem_wrt`, `div_exact_coeff` | `subresultant` | **Stable (crate-internal)** — 嵌套环除法 |
| `content_wrt_impl`, `primitive_part_wrt_impl` | `subresultant` | **Stable (crate-internal)** |
| `coeff_wrt_poly`, `content_wrt`, `primitive_part_wrt`, `substitute_poly`, `square_free_wrt`, `derivative_wrt`, `term_with_var` | `factor/poly_uni` | **Stable (bounded)** |
| `factor_sqff_over_coeff_ring` | `factor/poly_uni` | **Partial** — upstream `do_factor_hensel` 链 |
| `try_sparse_factor`, `try_sparse_factor_bi` | `factor/sparse` | **Partial** — FAC-G1 |
| `find_good_eval`, `looks_irreducible_by_good_eval` | `factor/eval` | **Partial** — 好点种子 / 不可约快检 |
| `factor_univariate_flat`, `factor_univariate_pairs` | `factor/univariate` | **Partial** |
| `try_zassenhaus_factor` | `factor/zassenhaus` | **Partial** |
| `try_hensel_lift_bivariate` | `factor/hensel` | **Partial** — FAC-G3 |
| `try_unitary_factor` | `factor/unitary` | **Partial** — FAC-G1；`unitaryfactor` / `pzadic` / P2a；失败 `None` |
| `unitary_factor_rev` | `factor/unitary` | **Partial** — 核心 peel 循环（`vars_rev`）；crate 内由 `try_unitary_factor` 调用 |
| `try_factor_patterns` | `factor/patterns` | **Partial** — cyclotomic/二项式模式 |
| `try_factor_xn_minus_one` 等 | `factor/cyclotomic` | **Partial** |
| `factor_fpx`, `degree` | `factor/fpx` | **Stable**（模域） |
| `normalize_univariate_factors`, `try_lift_factors_in_aux_var` | `factor/hensel` | **Pipeline private** `pub(crate)` |
| `find_rational_root` | `factor/univariate` | **Pipeline private** `pub(crate)` |
| `dense::poly1::*`, `Poly1RingCtx`, `Poly1Order` | `dense/poly1` | **Stable (crate-internal)** — 稠密 poly1 算术；见 [GIAC-dense-poly1-refactor](issues/GIAC-dense-poly1-refactor.md) D1–D3 |
| `dense::convert::*` | `dense/convert` | **Stable (crate-internal)** — HighFirst ↔ ascending ↔ sparse（D4） |
| `RatioRingCtx`, `RatioRingOps` | `dense/ratio_ring` | **Stable (crate-internal)** — ℚ 系数环 |
| `modpoly_to_poly` | `factor/modular.rs` | **Pipeline private** `pub(crate)` |

---

## 4. Pipeline private（模块内 `fn try_*` / 辅助）

| 函数 | 文件 | 职责 |
|------|------|------|
| `factor_multivariate_rec` | `multivariate.rs` | 多变量分解主递归 |
| `factor_wrt_main_var` | `multivariate.rs` | 按主变元分解 |
| `square_free_wrt_impl` | `poly_uni.rs` | Yun sqff 主循环 |
| `matching_embed_factor`, `reconstruct_factor_two_aux` | `sparse.rs` | sparse_bi 嵌入重建 |
| `try_hensel_lift_interp` | `hensel.rs` | Hensel 插值 fallback |
| `hensel_lift_at_zero` | `hensel.rs` | y=0 处 Hensel |
| `try_lift_and_peel` | `unitary.rs` | pzadic peel → P2a `lift_factor_multi_eval` |
| `lift_factor_multi_eval` | `unitary.rs` | P2a 局部窗 Lagrange 抬升 |
| `pzadic` (`PzadicLift`) | `unitary.rs` | base-`B` digit 抬升 |
| `unitarize` / `ununitarize` | `unitary.rs` | 非 monic 首项尾链 |
| `trunc1_drop_var` / `untrunc1_insert_var` | `unitary.rs` | 常数项尾部 `trunc1` |
| `reverse_var_order` | `unitary.rs` | sqff 块 `reverse()`（U5 边界待补） |
| `try_factor_biquadratic`, `try_factor_two_cubics` | `univariate.rs` | 低次模式 |
| `try_nth_root`, `try_binomial_square` | `power.rs` | 完美幂 |

**已退役（`#[cfg(test)]`，不得上热路径）：**

| 函数 | 文件 | 说明 |
|------|------|------|
| `try_kronecker_bivariate`, `try_factor_bivariate_eval`, `try_lift_bivariate_from_eval` | `poly_uni.rs` | 由 sparse→Hensel 覆盖 |

---

## 5. 上游缺失（非临时函数 — 真算法债）

| 缺口 ID | upstream (`gausspol.cc`) | giac-rs 状态 |
|---------|---------------------------|--------------|
| **FAC-G1** | `try_sparse_factor` + `try_sparse_factor_bi` + `unitaryfactor` | **Partial** — P0–P2b + **U-P2c tier 复审** ✅；line25 ✅；U5 `reverse` 边界待补。原理：[giac-poly-factor-unitary-principles](giac-poly-factor-unitary-principles.md)；缺口：[unitaryfactor-gaps](issues/GIAC-poly-unitaryfactor-gaps.md) |
| **FAC-G2** | 参系数 `poly_factor` 塔 | **Partial** — `try_lift_factors_in_aux_var` 覆盖 L20 |
| **FAC-G3** | 混合次数二元 Hensel + fallback | **Partial** — L22 ✅（`hensel_lift_two_at_zero`） |
| — | partfrac 重复二次 / 实二次分裂 | **Partial** — 线性/重根/实分裂已覆盖；高次仍缺 |

**退役目标:** FAC-G1 落地后，缩小 §4 中互斥的 `try_*` 形状链，统一经 `factor_multivariate_rec` + sparse fallback。

---

## 6. 测试锚点

| 测试 | 文件 | 状态 |
|------|------|------|
| `testfactor_line16/17/21/12/24` | `factor/tracer.rs` | enabled |
| `testfactor_line20` | `factor/tracer.rs` | enabled ✅ |
| `testfactor_line22` | `factor/tracer.rs` | enabled ✅ |
| `testfactor_line25` / `testfactor_line26` | `factor/tracer.rs` | enabled ✅（unitary + P2a gate） |
| `p2a_line25_*`, `unitary_factor_line25_*` | `factor/unitary.rs` | 单元测试 |
| `sparse_factor_*`, `hensel_*` | `factor/sparse.rs`, `hensel.rs` | 单元测试 |

---

## 7. 跨 crate 契约

| 调用方 | API | 要求 |
|--------|-----|------|
| `giac-simplify::factor` | `factor_into`, `factor_poly`, `vars_in`, `ratio_perfect_sqrt` | 失败时 rootof/sqrt 临时路径 |
| `giac-calculus::partfrac_integrate` | `partfrac_rational_terms`, `try_linear_power` | factor 失败则积分失败 |
| `giac-calculus::risch` | `tresultant_*`, `hermite` 用 `Poly` | |
| `giac-solve` | `sturm_*`, `roots`, `gcd_univariate` | |
| `giac-core` | `expr_to_poly` / `poly_to_expr` | 仅多项式子类 |

---

## 8. 维护

1. 新增 `pub fn` → 标注 **Stable|Partial** → 更新 §2/§3
2. 新 `try_*` fallback → 登记 §4；若对应 FAC-G*，更新 upstream-gaps issue
3. 禁止在 giac-simplify / giac-calculus 复制 factor 形状表
4. `python3 scripts/annotate_api_tiers.py --inventory` 刷新 Per-file 表

```bash
cd giac-rs
python3 scripts/annotate_api_tiers.py              # 为新 fn 补 tier 注释（幂等）
python3 scripts/annotate_api_tiers.py --inventory # 刷新本文 Per-file 表
```

---

## Per-file function inventory (generated)

Regenerate: `python3 scripts/annotate_api_tiers.py --inventory`

### `chinrem.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `chinrem` | **Stable** | Chinese remainder two residues |
| `chinrem_lists` | **Stable** | CRT fold over lists |
| `x` | **Pipeline private** | `x` |
| `chinrem_two_linear_moduli` | **Pipeline private** | `chinrem_two_linear_moduli` |

### `dense/convert.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `reverse_coeffs` | **Stable (crate-internal)** | reverse coefficient order (HighFirst ↔ Ascending). |
| `ascending_to_dense_high_first` | **Stable (crate-internal)** | ascending dense → giac `poly1` HighFirst. |
| `dense_high_first_to_ascending` | **Stable (crate-internal)** | giac `poly1` HighFirst → ascending dense. |
| `sparse_ascending_to_dense_high_first` | **Stable (crate-internal)** | sparse `Poly` univariate in `var` → dense HighFirst. |
| `dense_high_first_to_sparse` | **Stable (crate-internal)** | dense HighFirst → sparse univariate in `var`. |
| `poly_from_ascending_coeffs` | **Pipeline private** | `poly_from_ascending_coeffs` |

### `dense/poly1.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `zero` | **Stable** | Poly zero |
| `one` | **Stable** | Poly one |
| `is_zero` | **Stable** | Poly is zero |
| `add` | **Stable** | Poly addition |
| `sub` | **Stable** | Poly subtraction |
| `neg` | **Stable** | Poly negation |
| `mul` | **Stable** | Poly multiplication |
| `inv` | **Pipeline private** | `inv` |
| `is_one` | **Stable** | Poly is one |
| `div_coeff` | **Pipeline private** | `div_coeff` |
| `to_high_first` | **Pipeline private** | `to_high_first` |
| `from_high_first` | **Pipeline private** | `from_high_first` |
| `poly_degree` | **Stable (crate-internal)** | degree in the given order (zero poly has degree 0). |
| `trim_high_first` | **Pipeline private** | `trim_high_first` |
| `trim` | **Stable (crate-internal)** | drop redundant leading zeros; keep at least one coefficient. |
| `trim_high_first_collect` | **Pipeline private** | `trim_high_first_collect` |
| `add` | **Pipeline private** | `add` |
| `sub` | **Pipeline private** | `sub` |
| `mul` | **Pipeline private** | `mul` |
| `neg` | **Pipeline private** | `neg` |
| `scale` | **Pipeline private** | `scale` |
| `div_rem` | **Stable (crate-internal)** | polynomial division; quotient and remainder in `order`. |
| `reduce_mod_monic` | **Stable (crate-internal)** | reduce `p` modulo monic `m` (leading coeff of `m` is ±1). |
| `ext_gcd` | **Stable (crate-internal)** | extended GCD: `(g, s)` with `s*a + t*b = g` (`t` omitted). |
| `inv_mod` | **Stable (crate-internal)** | multiplicative inverse of `a` modulo monic `m`. |

### `dense/ratio_ring.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `zero` | **Stable** | Poly zero |
| `one` | **Stable** | Poly one |
| `is_zero` | **Stable** | Poly is zero |
| `add` | **Stable** | Poly addition |
| `sub` | **Stable** | Poly subtraction |
| `neg` | **Stable** | Poly negation |
| `mul` | **Stable** | Poly multiplication |
| `inv` | **Pipeline private** | `inv` |

### `dense/tests.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `q` | **Pipeline private** | `q` |
| `ctx` | **Pipeline private** | `ctx` |
| `high_first_degree_and_trim` | **Pipeline private** | `high_first_degree_and_trim` |
| `high_first_add_sub` | **Pipeline private** | `high_first_add_sub` |
| `high_first_mul_x_plus_1_times_x_plus_2` | **Pipeline private** | `high_first_mul_x_plus_1_times_x_plus_2` |
| `high_first_neg_and_scale` | **Pipeline private** | `high_first_neg_and_scale` |
| `high_first_div_rem` | **Pipeline private** | `high_first_div_rem` |
| `high_first_reduce_mod_x_squared_minus_2` | **Pipeline private** | `high_first_reduce_mod_x_squared_minus_2` |
| `high_first_reduce_skips_non_monic_modulus` | **Pipeline private** | `high_first_reduce_skips_non_monic_modulus` |
| `high_first_inv_mod_and_ext_gcd` | **Pipeline private** | `high_first_inv_mod_and_ext_gcd` |
| `ascending_order_matches_reversed_high_first_mul` | **Pipeline private** | `ascending_order_matches_reversed_high_first_mul` |
| `reverse_coeffs_involution` | **Pipeline private** | `reverse_coeffs_involution` |
| `sparse_dense_high_first_roundtrip` | **Pipeline private** | `sparse_dense_high_first_roundtrip` |
| `dense_div_rem_matches_univariate_div_rem_wrt` | **Pipeline private** | `dense_div_rem_matches_univariate_div_rem_wrt` |

### `exp.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `bigint_pow` | **Pipeline private** | BigInt pow with overflow check |

### `factor/ctx.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `new` | **Stable** | `new` |
| `as_slice` | **Stable** | `as_slice` |
| `len` | **Stable** | `len` |
| `is_empty` | **Stable** | `is_empty` |
| `refs` | **Stable** | `refs` |
| `new` | **Stable** | `new` |
| `main_degree` | **Stable** | `main_degree` |
| `with_main_and_others` | **Stable** | `with_main_and_others` |
| `irreducible` | **Stable** | `Poly::irreducible` |
| `from_polys` | **Stable** | `Poly::from_polys` |
| `into_polys` | **Stable** | `into_polys` |
| `product_equals` | **Stable** | `product_equals` |
| `verify_divides_chain` | **Stable** | `verify_divides_chain` |
| `new` | **Stable** | `new` |
| `first_value` | **Stable** | `Poly::first_value` |
| `try_new` | **Partial** | optional algorithm path `try_new` |
| `as_view` | **Stable** | `as_view` |
| `div_rem_wrt_aux_indep` | **Stable** | `div_rem_wrt_aux_indep` |
| `try_new` | **Partial** | optional algorithm path `try_new` |
| `factor_set_product_equals` | **Pipeline private** | `factor_set_product_equals` |
| `aux_indep_factor_rejects_aux_dep` | **Pipeline private** | `aux_indep_factor_rejects_aux_dep` |
| `hensel_pair_requires_aux_indep` | **Pipeline private** | `hensel_pair_requires_aux_indep` |

### `factor/cyclotomic.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `divisors_u64` | **Pipeline private** | Positive divisors of `n`, sorted ascending. |
| `cyclotomic_poly` | **Stable** | n-th cyclotomic polynomial Φ_n(x) over ℚ. |
| `factor_xn_minus_one` | **Partial** | `x^n - 1 = ∏_{d\|n} Φ_d(x)`. |
| `factor_x2n_plus_xn_plus_1` | **Partial** | Factors of `x^(2n)+x^n+1 = (x^(3n)-1)/(x^n-1)` via cyclotomic selection. |
| `try_factor_x2n_plus_xn_plus_1` | **Partial** | detect and factor x^2n+x^n+1 |
| `is_x2n_plus_xn_plus_1_sparse` | **Pipeline private** | `is_x2n_plus_xn_plus_1_sparse` |
| `try_factor_xn_minus_one` | **Partial** | detect x^n-1 |
| `is_xn_minus_one_poly` | **Pipeline private** | shape test x^n-1 |
| `try_factor_xn_plus_one` | **Partial** | detect x^n+1 |
| `factor_xn_plus_one` | **Pipeline private** | `factor_xn_plus_one` |
| `cyclotomic_phi3` | **Pipeline private** | `cyclotomic_phi3` |
| `factor_x100_plus_x50_plus_1` | **Pipeline private** | `factor_x100_plus_x50_plus_1` |
| `factor_x10_minus_1` | **Pipeline private** | `factor_x10_minus_1` |

### `factor/eval.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `peval_at_main` | **Pipeline private** | substitute auxiliary vars with scalars; keep `main` univariate. |
| `find_good_eval` | **Pipeline private** | find evaluation preserving `main`-degree (upstream `find_good_eval`). |
| `eval_point_candidates` | **Pipeline private** | trial points: `start` first, then small integers / shifts. |
| `looks_irreducible_by_good_eval` | **Pipeline private** | upstream `do_factor_hensel`: two good evals, single factor → irreducible. |
| `find_good_eval_preserves_degree` | **Pipeline private** | `find_good_eval_preserves_degree` |
| `find_good_eval_skips_bad_start` | **Pipeline private** | `find_good_eval_skips_bad_start` |
| `irreducibility_probe_detects_x2_plus_y2_plus_1` | **Pipeline private** | `irreducibility_probe_detects_x2_plus_y2_plus_1` |
| `irreducibility_probe_rejects_reducible` | **Pipeline private** | `irreducibility_probe_rejects_reducible` |

### `factor/fpx.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `x_var` | **Pipeline private** | `x_var` |
| `mi` | **Pipeline private** | `mi` |
| `is_poly_one` | **Pipeline private** | `is_poly_one` |
| `degree` | **Stable** | total degree in `x` |
| `coeff` | **Pipeline private** | `coeff` |
| `from_coeffs` | **Pipeline private** | `from_coeffs` |
| `x_poly` | **Pipeline private** | `x_poly` |
| `one_poly` | **Pipeline private** | `one_poly` |
| `mod_poly` | **Pipeline private** | `mod_poly` |
| `div_exact` | **Stable** | exact division if remainder zero |
| `make_monic` | **Pipeline private** | `make_monic` |
| `derivative` | **Pipeline private** | `derivative` |
| `eval` | **Pipeline private** | `eval` |
| `powmod` | **Pipeline private** | `powmod` |
| `compose` | **Pipeline private** | `compose` |
| `subst_x_to_xp` | **Pipeline private** | `subst_x_to_xp` |
| `linear_factor` | **Pipeline private** | `linear_factor` |
| `new` | **Pipeline private** | `new` |
| `next_u64` | **Stable** | `Poly::next_u64` |
| `next_i64` | **Stable** | `Poly::next_i64` |
| `random_poly` | **Pipeline private** | `random_poly` |
| `square_free_yun` | **Pipeline private** | `square_free_yun` |
| `distinct_degree_factorization` | **Pipeline private** | `distinct_degree_factorization` |
| `extract_linear_factors` | **Pipeline private** | `extract_linear_factors` |
| `cantor_zassenhaus_block` | **Pipeline private** | `cantor_zassenhaus_block` |
| `factor_square_free` | **Pipeline private** | `factor_square_free` |
| `factor_fpx` | **Stable** | Full factorization in F_p[x] into monic irreducible factors (with repetition). |
| `poly` | **Pipeline private** | `poly` |
| `assert_product` | **Pipeline private** | `assert_product` |
| `fpx_x4_plus_1_mod_5` | **Pipeline private** | `fpx_x4_plus_1_mod_5` |
| `fpx_x6_minus_1_mod_7` | **Pipeline private** | `fpx_x6_minus_1_mod_7` |
| `fpx_x2_plus_1_mod_5` | **Pipeline private** | `fpx_x2_plus_1_mod_5` |
| `fpx_cubic_irreducible_mod_7` | **Pipeline private** | `fpx_cubic_irreducible_mod_7` |

### `factor/fpx_uni.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `monomial_pow` | **Pipeline private** | `monomial_pow` |
| `mod_int` | **Pipeline private** | `mod_int` |
| `univariate_degree` | **Pipeline private** | `univariate_degree` |
| `coeff_at` | **Pipeline private** | `coeff_at` |
| `set_coeff` | **Pipeline private** | `set_coeff` |
| `from_modint_coeffs` | **Pipeline private** | `from_modint_coeffs` |
| `from_bigint_coeffs` | **Pipeline private** | `from_bigint_coeffs` |
| `var_poly` | **Pipeline private** | `var_poly` |
| `make_monic` | **Pipeline private** | `make_monic` |
| `derivative` | **Pipeline private** | `derivative` |
| `at_modulus` | **Pipeline private** | `at_modulus` |
| `to_centered_int_poly` | **Pipeline private** | `to_centered_int_poly` |

### `factor/hensel.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `trim_rat` | **Pipeline private** | `trim_rat` |
| `rat_mul` | **Pipeline private** | `rat_mul` |
| `rat_div_rem` | **Pipeline private** | `rat_div_rem` |
| `rat_egcd` | **Pipeline private** | `rat_egcd` |
| `rat_add` | **Pipeline private** | `rat_add` |
| `rat_sub` | **Pipeline private** | `rat_sub` |
| `poly_univariate_rat` | **Pipeline private** | `poly_univariate_rat` |
| `poly_from_rat` | **Pipeline private** | `poly_from_rat` |
| `egcd_factor_list` | **Pipeline private** | `egcd_factor_list` |
| `truncate_y` | **Pipeline private** | `truncate_y` |
| `div_rem_x_over_qy` | **Pipeline private** | `div_rem_x_over_qy` via nested-ring API (divisor independent of aux) |
| `leading_coeff_x` | **Pipeline private** | `leading_coeff_x` |
| `scale_univariate_x` | **Pipeline private** | `scale_univariate_x` |
| `hensel_lift_two_at_zero` | **Pipeline private** | legacy y-degree truncation Hensel (2 factors); fallback when total-degree lift fails |
| `normalize_univariate_factors` | **Pipeline private** | `normalize_univariate_factors` |
| `truncate_total_degree` | **Pipeline private** | `truncate_total_degree` |
| `total_degree_poly` | **Pipeline private** | `total_degree_poly` |
| `scalar_at_y_zero` | **Pipeline private** | scalar value of `p(y=0)` when constant |
| `lcp_depends_on_aux` | **Pipeline private** | true when `lcp` depends on auxiliary var `aux` |
| `lcm_rat` | **Pipeline private** | lcm of two rationals |
| `lcm_poly_denoms` | **Pipeline private** | lcm of denominators of all coeffs in `p` |
| `scale_f0_factors` | **Pipeline private** | scale `f0[i]` by `lcoeff(y=0)/lc(f0[i])` (upstream `mulmodpoly` on `F0fact`) |
| `build_hensel_lift_seeds` | **Pipeline private** | `build_hensel_lift_seeds` |
| `egcd_factor_list_normalized` | **Pipeline private** | `egcd_factor_list` + lcm denominator `D` (upstream `lcmdeno` / `D`) |
| `hensel_lift_factor_loop` | **Pipeline private** | total-degree Hensel iteration (upstream `EZGCD_DEGONLY`, `b=0`) |
| `try_hensel_lift_factor` | **Pipeline private** | `try_hensel_lift_factor` |
| `normalize_hensel_lift_result` | **Pipeline private** | normalize lifted factors whose product is `p` or `p_adj` |
| `hensel_lift_at_zero` | **Pipeline private** | `hensel_lift_at_zero` |
| `as_rational_constant` | **Pipeline private** | `as_rational_constant` |
| `linear_root` | **Pipeline private** | `linear_root` |
| `factor_match_key_at` | **Pipeline private** | `factor_match_key_at` |
| `sort_factors_by_match_key` | **Pipeline private** | `sort_factors_by_match_key` |
| `factor_match_key` | **Pipeline private** | `factor_match_key` |
| `lagrange_interpolate_y` | **Pipeline private** | `lagrange_interpolate_y` |
| `lift_factor_from_evals` | **Pipeline private** | `lift_factor_from_evals` |
| `lift_factor_from_aux_evals` | **Pipeline private** | Hensel lift factor from auxiliary eval tracks |
| `lift_factor_from_aux_evals_with_rest` | **Pipeline private** | `lift_factor_from_aux_evals_with_rest` |
| `try_lift_factors_in_aux_var` | **Pipeline private** | lift bivariate factors via aux variable |
| `try_hensel_lift_interp` | **Pipeline private** | optional fallback `try_hensel_lift_interp` |
| `try_hensel_lift_bivariate` | **Partial** | Factor `p(x,y)` in ℚ[y][x]: Hensel lift at `aux=0`, both main orders, then interpolation fallback. |
| `hensel_three_linear_shifted` | **Pipeline private** | `hensel_three_linear_shifted` |
| `hensel_two_bilinear_factors` | **Pipeline private** | `hensel_two_bilinear_factors` |
| `hensel_egcd_two_factors` | **Pipeline private** | `hensel_egcd_two_factors` |
| `hensel_egcd_factor_list` | **Pipeline private** | `hensel_egcd_factor_list` |
| `div_rem_xy_by_x_minus_one` | **Pipeline private** | `div_rem_xy_by_x_minus_one` |
| `hensel_two_bilinear_at_zero` | **Pipeline private** | `hensel_two_bilinear_at_zero` |
| `factor_non_monic_product_at_zero` | **Pipeline private** | `factor_non_monic_product_at_zero` |
| `hensel_two_simple_non_monic` | **Pipeline private** | `hensel_two_simple_non_monic` |
| `aux_lift_line21` | **Pipeline private** | `aux_lift_line21` |
| `try_hensel_lift_factor_line22_debug_steps` | **Pipeline private** | `try_hensel_lift_factor_line22_debug_steps` |
| `try_hensel_lift_factor_line22` | **Pipeline private** | `try_hensel_lift_factor_line22` |
| `two_factor_hensel_at_zero_line22` | **Pipeline private** | `two_factor_hensel_at_zero_line22` |
| `interp_line22_mixed_bivariate` | **Pipeline private** | `interp_line22_mixed_bivariate` |
| `hensel_line22_mixed_bivariate` | **Pipeline private** | `hensel_line22_mixed_bivariate` |
| `hensel_uses_lift_at_zero_for_linears` | **Pipeline private** | `hensel_uses_lift_at_zero_for_linears` |

### `factor/mod.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `factor_into` | **Stable (bounded)** | Factor into irreducible polynomial factors over ℚ when possible. |
| `factor_poly` | **Stable (bounded)** | Integer-style factorization display (legacy `factor_poly`). |
| `factor_into_by_rational_roots` | **Stable (bounded)** | univariate via rational roots |
| `factor_poly_mod` | **Stable** | factor mod p display |
| `factor_mod_irreducibles` | **Stable** | Irreducible factors over F_p (monic, with repetition). |
| `as_perfect_power_quadratic_squared` | **Pipeline private** | `as_perfect_power_quadratic_squared` |
| `factor_multivariate_returns_none_without_hanging` | **Pipeline private** | `factor_multivariate_returns_none_without_hanging` |
| `factor_x_fourth_minus_one` | **Pipeline private** | `factor_x_fourth_minus_one` |
| `factor_x_cubed_plus_one` | **Pipeline private** | `factor_x_cubed_plus_one` |
| `factor_x_times_x_squared_plus_one` | **Pipeline private** | `factor_x_times_x_squared_plus_one` |
| `factor_x6_minus_y6` | **Pipeline private** | `factor_x6_minus_y6` |
| `factor_x100_plus_x50_plus_1` | **Pipeline private** | `factor_x100_plus_x50_plus_1` |
| `try_linear_power_detects_square` | **Pipeline private** | `try_linear_power_detects_square` |

### `factor/modular.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `factor_poly_mod` | **Stable** | Factor over ℤ/pℤ then lift display (giac `mod_factor` subset). |
| `modpoly_to_poly` | **Pipeline private** | PolyMod → Poly over ℤ/pℤ for display |
| `x` | **Pipeline private** | `x` |
| `factor_x4_plus_1_mod_5` | **Pipeline private** | `factor_x4_plus_1_mod_5` |
| `factor_x6_minus_1_mod_7` | **Pipeline private** | `factor_x6_minus_1_mod_7` |
| `factor_x2_plus_1_mod_5` | **Pipeline private** | `factor_x2_plus_1_mod_5` |
| `factor_x4_minus_1_mod_2` | **Pipeline private** | `factor_x4_minus_1_mod_2` |
| `factor_has_correct_product_mod_11` | **Pipeline private** | `factor_has_correct_product_mod_11` |

### `factor/multivariate.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `factor_into_poly` | **Stable (bounded)** | factor_multivariate ok→Some |
| `factor_multivariate` | **Stable (bounded)** | multivariate factorization |
| `factor_multivariate_rec` | **Pipeline private** | multivariate factor recursion (patterns→uni→main var) |
| `factor_wrt_main_var` | **Pipeline private** | `factor_wrt_main_var` |
| `factor_multivariate_rec_sqff` | **Pipeline private** | [`SqffFactorRecFn`] adapter for multivariate recursion |
| `factor_bivariate_product` | **Pipeline private** | `factor_bivariate_product` |
| `factor_bivariate_mixed` | **Pipeline private** | `factor_bivariate_mixed` |
| `factor_three_shifted_linears` | **Pipeline private** | `factor_three_shifted_linears` |
| `factor_var_power_times_linear` | **Pipeline private** | `factor_var_power_times_linear` |
| `factor_repeated_linear_pairs` | **Pipeline private** | `factor_repeated_linear_pairs` |

### `factor/patterns.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `try_factor_patterns` | **Partial** | cyclotomic/binomial pattern table |
| `factor_xn_minus_yn` | **Pipeline private** | `factor_xn_minus_yn` |
| `is_binomial_diff_power` | **Pipeline private** | `is_binomial_diff_power` |
| `factor_xn_minus_yn_explicit` | **Pipeline private** | `factor_xn_minus_yn_explicit` |
| `factor_xn_minus_one_display` | **Partial** | Legacy wrapper used by `factor_poly`. |

### `factor/poly_uni.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `coeff_wrt_poly` | **Stable** | Coefficient of `var^exp` as a polynomial in the remaining variables. |
| `content_wrt` | **Stable** | Content of `p` w.r.t. `var`: gcd of all x-coefficients in ℚ[others]. |
| `primitive_part_wrt` | **Stable** | `p / content_wrt(p, var)` in ℚ[others][var]. |
| `term_with_var` | **Stable** | coeff * var^exp as Poly |
| `derivative_wrt` | **Stable** | ∂p/∂var treating coefficients in ℚ[others]. |
| `square_free_wrt` | **Stable** | Square-free factorization w.r.t. `var` over ℚ[others] (Yun-style via gcd). |
| `square_free_wrt_impl` | **Pipeline private** | `square_free_wrt_impl` |
| `substitute_poly` | **Stable** | Substitute `sub_var -> sub_poly` in `p`. |
| `factor_sqff_over_coeff_ring` | **Partial** | Factor square-free `g` in ℚ[others][var] recursively. |
| `factor_sqff_over_coeff_ring_ctx` | **Pipeline private** | typed sqff factor chain |
| `as_constant` | **Pipeline private** | `as_constant` |
| `as_constant` | **Stable** | `Poly::as_constant` |
| `sqff_rec` | **Pipeline private** | `sqff_rec` |
| `sqff_chain_l20_via_factor_sqff_over_coeff_ring` | **Pipeline private** | `sqff_chain_l20_via_factor_sqff_over_coeff_ring` |
| `sqff_chain_l21_ternary_linears` | **Pipeline private** | `sqff_chain_l21_ternary_linears` |
| `content_wrt_xy_plus_y_squared` | **Pipeline private** | `content_wrt_xy_plus_y_squared` |
| `rational_content_vs_wrt_content` | **Pipeline private** | `rational_content_vs_wrt_content` |

### `factor/power.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `as_perfect_power` | **Stable** | If `p` is a perfect power, return `(base, exponent)`. |
| `try_nth_root` | **Pipeline private** | optional fallback `try_nth_root` |
| `try_binomial_square` | **Pipeline private** | optional fallback `try_binomial_square` |
| `try_linear_power` | **Partial** | detect (linear)^n |

### `factor/sparse.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `try_sparse_factor` | **Partial** | Sparse reconstruction from univariate factors at `other = 0` (FAC-G1). |
| `try_sparse_factor_at` | **Partial** | sparse factor with optional auxiliary evaluation point (upstream `b0`). |
| `try_sparse_factor_impl` | **Pipeline private** | single `(main, other)` attempt via typed sparse stages |
| `prepare` | **Pipeline private** | `prepare` |
| `into_system` | **Pipeline private** | `into_system` |
| `solve` | **Pipeline private** | `solve` |
| `leading_coeff_main` | **Pipeline private** | leading coefficient of `p` in `main` (poly in remaining vars) |
| `factor_unknown_count` | **Pipeline private** | non-leading term count in univariate factor `f` |
| `build_factor_templates` | **Pipeline private** | assign global `la` indices to non-leading terms |
| `scale_by_lcp_power` | **Pipeline private** | `lcp^(pow) * p` |
| `zero` | **Stable** | Poly zero |
| `is_zero` | **Stable** | Poly is zero |
| `add` | **Stable** | Poly addition |
| `mul` | **Stable** | Poly multiplication |
| `is_zero` | **Stable** | Poly is zero |
| `from_parts` | **Stable** | `Poly::from_parts` |
| `is_linear` | **Stable** | `Poly::is_linear` |
| `substitute` | **Pipeline private** | `substitute` |
| `merge_linear` | **Pipeline private** | `merge_linear` |
| `merge_bilinear` | **Pipeline private** | `merge_bilinear` |
| `poly_to_sparse_map` | **Pipeline private** | known `(main, other)` exponent map |
| `template_to_sparse_map` | **Pipeline private** | one factor template as sparse map with `la` unknowns |
| `sparse_map_mul` | **Pipeline private** | multiply sparse maps, tracking bilinear per monomial |
| `sparse_map_sub` | **Pipeline private** | `sparse_map_sub` |
| `build_sparse_equations` | **Pipeline private** | `build_sparse_equations` |
| `c_from_target_x0` | **Pipeline private** | `x^0` coefficients of `target` by `other`-degree |
| `lcp_coeffs_in_other` | **Pipeline private** | `lcp` as univariate poly in `other`: `[lcp_0, lcp_1, …]` |
| `x1_target_coeffs` | **Pipeline private** | `x^1` coefficients of `target` by `other`-degree |
| `solve_lcp_convolution_x1` | **Pipeline private** | `solve_lcp_convolution_x1` |
| `solve_sparse_two_factor` | **Pipeline private** | 2-factor path: `b_k = s_k - a_k`, bilinear in `a_k` |
| `a0_quadratic_roots` | **Pipeline private** | `a0_quadratic_roots` |
| `complete_two_factor_a` | **Pipeline private** | `complete_two_factor_a` |
| `bilinear_at_m` | **Pipeline private** | `bilinear_at_m` |
| `solve_a_j_at_m` | **Pipeline private** | `solve_a_j_at_m` |
| `rational_sqrt` | **Pipeline private** | `rational_sqrt` |
| `integer_sqrt_bigint` | **Pipeline private** | `integer_sqrt_bigint` |
| `solve_sparse_system` | **Pipeline private** | upstream iterative linear extraction + bilinear finish |
| `solve_linear_equations` | **Pipeline private** | extract one linear equation and solve via Gaussian elimination |
| `solve_bilinear_remaining` | **Pipeline private** | handle remaining bilinear equations (2-factor `A*B` terms) |
| `gaussian_elimination` | **Pipeline private** | Gaussian elimination over ℚ |
| `la_poly_from_sol` | **Pipeline private** | `la_poly_from_sol` |
| `build_factors` | **Pipeline private** | `build_factors` |
| `divide_poly_coeffs_by` | **Pipeline private** | divide each `main`-coefficient by `den` when exact |
| `try_adjust_sparse_scale` | **Pipeline private** | `try_adjust_sparse_scale` |
| `verify_sparse_factors` | **Pipeline private** | accept only `∏ f_i = p` (upstream `divbylgcd` applied above) |
| `try_sparse_factor_bi` | **Partial** | sparse factor via bivariate `eval_tn` embedding (FAC-G1, 2+ aux). |
| `embed_sorted_monomials` | **Pipeline private** | build `Poly` from sorted embed monomials + aux exponents |
| `bivariate_x_degrees_ok` | **Pipeline private** | distinct `main`-degrees with pairwise distinct coeffs (upstream `x_degrees`). |
| `factor_bivariate_flat` | **Stable** | `factor_bivariate_flat` |
| `select_bivariate_factor` | **Pipeline private** | pick sparsest bivariate factor candidate. |
| `matching_embed_factor` | **Pipeline private** | factor of `eval_tn(p)` with matching `main`-degree pattern. |
| `reconstruct_sparse_bi_monomials` | **Pipeline private** | `reconstruct_sparse_bi_monomials` |
| `reconstruct_factor_dual_embed` | **Pipeline private** | `reconstruct_factor_dual_embed` |
| `reconstruct_sparse_bi_factor` | **Pipeline private** | `reconstruct_sparse_bi_factor` |
| `exact_quo_wrt` | **Pipeline private** | exact quotient `p / factor` in ℚ[others][main]. |
| `dilate_poly` | **Pipeline private** | substitute `aux -> factor * aux` |
| `undilate_poly` | **Pipeline private** | undo dilation: `factor*aux -> aux` |
| `try_dilation_sparse_bi` | **Pipeline private** | upstream random dilation fallback (deterministic seeds) |
| `try_sparse_factor_bi_two_aux` | **Pipeline private** | two auxiliary variables; upstream `n` sweep then dilation fallback. |
| `try_sparse_factor_bi_two_aux_inner` | **Pipeline private** | optional fallback `try_sparse_factor_bi_two_aux_inner` |
| `try_sparse_factor_bi_single_n` | **Pipeline private** | one embed exponent vector `n` |
| `try_heuristic_factor_bivariate` | **Partial** | Heuristic factorization via large eval + `pzadic` lift (FAC-G1/G3). |
| `l22_poly` | **Pipeline private** | `l22_poly` |
| `lcp_convolution_deconv_line22` | **Pipeline private** | `lcp_convolution_deconv_line22` |
| `factor_bivariate_flat_bilinear` | **Pipeline private** | `factor_bivariate_flat_bilinear` |
| `factor_bivariate_flat_univariate_fallback` | **Pipeline private** | `factor_bivariate_flat_univariate_fallback` |
| `sparse_factor_tri_var` | **Pipeline private** | `sparse_factor_tri_var` |
| `sparse_factor_tri_var_sum_coeff` | **Pipeline private** | `sparse_factor_tri_var_sum_coeff` |
| `sparse_factor_bilinear` | **Pipeline private** | `sparse_factor_bilinear` |
| `sparse_factor_non_monic_quadratic` | **Pipeline private** | `sparse_factor_non_monic_quadratic` |
| `heuristic_factors_line22` | **Pipeline private** | `heuristic_factors_line22` |

### `factor/sqrt.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `quadratic_sqrt_factor_exprs` | **Partial** | Display factors `(x - (-b ± sqrt(disc))/(2a))` for a univariate quadratic. |
| `format_ratio` | **Pipeline private** | `format_ratio` |
| `format_sqrt_ratio` | **Pipeline private** | `format_sqrt_ratio` |
| `sqrt_factor_x2_minus_2` | **Pipeline private** | `sqrt_factor_x2_minus_2` |

### `factor/tower.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `new` | **Stable** | `new` |
| `from_slice` | **Stable** | `from_slice` |
| `inner_dim` | **Stable** | `inner_dim` |
| `as_slice` | **Stable** | `as_slice` |
| `new` | **Stable** | `new` |
| `from_sqff_ctx` | **Stable** | `from_sqff_ctx` |
| `main_degree` | **Stable** | `main_degree` |
| `flat_dim` | **Stable** | `flat_dim` |
| `all_vars` | **Stable** | `all_vars` |
| `is_parametric_tower` | **Stable** | `is_parametric_tower` |
| `unsplit_to_flat` | **Stable** | `unsplit_to_flat` |
| `split_factor` | **Stable** | `split_factor` |
| `try_factor_aux_lift` | **Partial** | optional algorithm path `try_factor_aux_lift` |
| `factor_sqff_chain` | **Stable** | `factor_sqff_chain` |
| `try_sparse_bi` | **Partial** | optional algorithm path `try_sparse_bi` |
| `try_factor` | **Partial** | optional algorithm path `try_factor` |
| `as_poly_factor_tower` | **Stable** | `as_poly_factor_tower` |
| `l20_poly` | **Pipeline private** | `l20_poly` |
| `l20_others` | **Pipeline private** | `l20_others` |
| `l21_poly` | **Pipeline private** | `l21_poly` |
| `tower_factor_sqff_chain_l20` | **Pipeline private** | `tower_factor_sqff_chain_l20` |
| `tower_factor_sqff_chain_l21` | **Pipeline private** | `tower_factor_sqff_chain_l21` |
| `tower_try_factor_skips_good_eval` | **Pipeline private** | `tower_try_factor_skips_good_eval` |
| `l24_poly` | **Pipeline private** | `l24_poly` |
| `tower_l24_irreducible_skips_sparse_bi` | **Pipeline private** | `tower_l24_irreducible_skips_sparse_bi` |
| `tower_l23_cubic_factors_are_irreducible` | **Pipeline private** | `tower_l23_cubic_factors_are_irreducible` |
| `tower_unsplit_is_flat_poly` | **Pipeline private** | `tower_unsplit_is_flat_poly` |
| `tower_aux_lift_factors_l20` | **Pipeline private** | `tower_aux_lift_factors_l20` |
| `tower_from_sqff_ctx` | **Pipeline private** | `tower_from_sqff_ctx` |
| `split_factor_tags_main` | **Pipeline private** | `split_factor_tags_main` |

### `factor/tracer.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `x` | **Pipeline private** | `x` |
| `y` | **Pipeline private** | `y` |
| `z` | **Pipeline private** | `z` |
| `b_var` | **Pipeline private** | `b_var` |
| `c_var` | **Pipeline private** | `c_var` |
| `assert_factors` | **Pipeline private** | `assert_factors` |
| `testfactor_line16_three_shifted_linears` | **Pipeline private** | `testfactor_line16_three_shifted_linears` |
| `testfactor_line17_two_quadratics` | **Pipeline private** | `testfactor_line17_two_quadratics` |
| `testfactor_line20_parametric_cubic_factor` | **Pipeline private** | `testfactor_line20_parametric_cubic_factor` |
| `testfactor_line21_three_linear_ternary` | **Pipeline private** | `testfactor_line21_three_linear_ternary` |
| `testfactor_line22_bivariate_mixed_degree` | **Pipeline private** | `testfactor_line22_bivariate_mixed_degree` |
| `testfactor_line12_two_cubics` | **Pipeline private** | `testfactor_line12_two_cubics` |
| `testfactor_line24_x6_minus_y6` | **Pipeline private** | `testfactor_line24_x6_minus_y6` |
| `testfactor_unitary_bilinear_gate` | **Pipeline private** | `testfactor_unitary_bilinear_gate` |
| `testfactor_line25_unitaryfactor_gate` | **Pipeline private** | `testfactor_line25_unitaryfactor_gate` |
| `testfactor_line26_y3_x2y_gate` | **Pipeline private** | `testfactor_line26_y3_x2y_gate` |

### `factor/unitary.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `initial` | **Pipeline private** | GCDHEU eval base 2·‖p‖∞+2 |
| `base` | **Pipeline private** | eval base accessor |
| `set_base` | **Pipeline private** | set eval base |
| `bump_sqff` | **Pipeline private** | sqff micro-bump base += 1 |
| `advance` | **Pipeline private** | upstream eval step ⌊base·73794/27011⌋+1 |
| `as_ratio` | **Pipeline private** | eval base as Ratio<BigInt> |
| `new` | **Pipeline private** | EvalBaseStream from poly norm |
| `current` | **Pipeline private** | current EvalBaseStream point |
| `current_mut` | **Pipeline private** | mutable current EvalBaseStream point |
| `next` | **Pipeline private** | advance outer eval base or stop when bits > 256 |
| `upstream_bases` | **Pipeline private** | first N bases on EvalBaseStream (tests) |
| `reverse_var_order` | **Pipeline private** | upstream tensor reverse on variable indices |
| `trunc1_drop_var` | **Pipeline private** | drop eval_var tail (upstream trunc1) |
| `untrunc1_insert_var` | **Pipeline private** | reinsert eval_var with zero exp (upstream untrunc1) |
| `eval_coeff_groups` | **Pipeline private** | group terms by eval_var exponent |
| `unitarize` | **Pipeline private** | scale to unitary leading coeff w.r.t. eval_var |
| `ununitarize` | **Pipeline private** | undo unitarize scaling factor |
| `pow_poly` | **Pipeline private** | integer exponentiation in Poly ring |
| `new` | **Pipeline private** | PzadicLift builder |
| `draft_from` | **Pipeline private** | build PzadicDraft from eval factor |
| `lift_candidates` | **Pipeline private** | pzadic lift candidates from draft |
| `pzadic` | **Pipeline private** | faithful base-B digit lift (dim+1 via eval_var) |
| `factor_sort_key` | **Pipeline private** | sort key (deg, lc) for eval factors |
| `sort_eval_factors` | **Pipeline private** | stable sort eval factor slots |
| `lagrange_interp_coeff` | **Pipeline private** | P2a Lagrange coeff in eval_var |
| `monic_wrt_main` | **Pipeline private** | normalize factor monic w.r.t. main |
| `lift_factor_multi_eval` | **Pipeline private** | P2a local-window multi-point coeff lift |
| `try_lift_and_peel` | **Pipeline private** | pzadic peel then P2a fallback + divides check |
| `sym_mod_digit` | **Pipeline private** | symmetric mod digit for pzadic expansion |
| `unreverse_factors` | **Pipeline private** | inverse of upstream `tensor::reverse()` on factors |
| `try_unitary_factor` | **Partial** | FAC-G1 last-resort; sparse/Hensel fallback; bounded GCDHEU eval stream |
| `unitary_factor_rev` | **Partial** | core unitaryfactor loop on vars_rev; pzadic peel + P2a fallback |
| `factor_constant_tail_into` | **Pipeline private** | recurse constant tail into factor list |
| `factor_constant_tail` | **Pipeline private** | factor tail when main degree → 0 |
| `try_peel_all_at_eval` | **Pipeline private** | batch peel all slots at one eval base |
| `factor_at_eval` | **Pipeline private** | factor eval image w.r.t. main |
| `verified_product` | **Pipeline private** | check factor product equals orig |
| `is_sqff_wrt_main` | **Pipeline private** | sqff test w.r.t. main var |
| `linfnorm` | **Pipeline private** | L∞ norm of Poly coefficients |
| `l22_y3_product` | **Pipeline private** | `l22_y3_product` |
| `l22_y3_x2y_product` | **Pipeline private** | `l22_y3_x2y_product` |
| `pzadic_lift_linear_via_digits` | **Pipeline private** | `pzadic_lift_linear_via_digits` |
| `reverse_roundtrip` | **Pipeline private** | `reverse_roundtrip` |
| `eval_base_stream_upstream_only` | **Pipeline private** | `eval_base_stream_upstream_only` |
| `p2a_line25_upstream_trajectory` | **Pipeline private** | `p2a_line25_upstream_trajectory` |
| `p2a_sample_window_anchors_below_base0` | **Pipeline private** | `p2a_sample_window_anchors_below_base0` |
| `unitary_bilinear_via_single_entry` | **Pipeline private** | `unitary_bilinear_via_single_entry` |
| `unitarize_extracts_leading_coeff` | **Pipeline private** | `unitarize_extracts_leading_coeff` |
| `unitary_factor_line25_l22_y3` | **Pipeline private** | `unitary_factor_line25_l22_y3` |
| `unitary_factor_line26_l22_y3_x2y` | **Pipeline private** | `unitary_factor_line26_l22_y3_x2y` |

### `factor/univariate.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `factor_univariate_flat` | **Partial** | Flat irreducible (or fully split) factor list. |
| `factor_univariate_pairs` | **Partial** | pairs with multiplicity |
| `factor_square_free` | **Pipeline private** | `factor_square_free` |
| `factor_by_rational_roots` | **Pipeline private** | `factor_by_rational_roots` |
| `factor_power_pairs` | **Stable (bounded)** | factors with multiplicities |
| `factor_power_pairs_core` | **Pipeline private** | `factor_power_pairs_core` |
| `find_rational_root` | **Pipeline private** | rational root via rational root theorem |
| `factor_quadratic` | **Pipeline private** | `factor_quadratic` |
| `try_factor_biquadratic` | **Pipeline private** | optional fallback `try_factor_biquadratic` |
| `monic_cubic_poly` | **Pipeline private** | `monic_cubic_poly` |
| `try_factor_two_cubics` | **Pipeline private** | optional fallback `try_factor_two_cubics` |
| `term_with_var` | **Stable** | coeff * var^exp as Poly |

### `factor/util.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `vars_in` | **Stable** | All variables appearing in `p`, lexicographically sorted. |
| `is_univariate_in` | **Stable** | `is_univariate_in` |
| `main_var` | **Stable** | Variable of minimum degree (giac `factor_multivar` main var heuristic). |
| `min_var_exponents` | **Stable** | Minimum exponent of each variable across all terms (missing var counts as 0). |
| `extract_var_power_factors` | **Stable** | Split `p = (∏ v^{e_v}) * rest` where `e_v` is the minimum exponent of `v` in `p`. |
| `coeff_gcd` | **Stable** | Integer gcd of rational coefficients. |
| `primitive_part` | **Stable** | divide out content |
| `linear_poly` | **Stable** | `linear_poly` |
| `monic_quadratic_poly` | **Stable** | `monic_quadratic_poly` |
| `integer_divisors` | **Stable** | `integer_divisors` |
| `integer_nth_root` | **Stable** | `integer_nth_root` |
| `rational_nth_root` | **Stable** | `rational_nth_root` |
| `ratio_perfect_sqrt` | **Stable** | detect perfect square Ratio |
| `rational_factor_pairs` | **Stable** | `rational_factor_pairs` |
| `is_monic_univariate` | **Stable** | `is_monic_univariate` |

### `factor/zassenhaus.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `is_int_poly` | **Pipeline private** | `is_int_poly` |
| `integer_coeffs` | **Pipeline private** | `integer_coeffs` |
| `mignotte_bound` | **Pipeline private** | `mignotte_bound` |
| `poly_mod_from_poly` | **Pipeline private** | `poly_mod_from_poly` |
| `is_square_free_mod` | **Pipeline private** | `is_square_free_mod` |
| `extgcd_mod` | **Pipeline private** | `extgcd_mod` |
| `egcd_factor_list_mod` | **Pipeline private** | `egcd_factor_list_mod` |
| `divides_exact_quotient` | **Pipeline private** | `divides_exact_quotient` |
| `coeff_div_mod` | **Pipeline private** | `coeff_div_mod` |
| `lift_correction` | **Pipeline private** | `lift_correction` |
| `hensel_lift_two` | **Pipeline private** | `hensel_lift_two` |
| `hensel_lift_n` | **Pipeline private** | `hensel_lift_n` |
| `polymod_product` | **Pipeline private** | `polymod_product` |
| `recover_factors_from_lifted` | **Pipeline private** | `recover_factors_from_lifted` |
| `extract_factors_via_combine` | **Pipeline private** | `extract_factors_via_combine` |
| `extract_combine_rec` | **Pipeline private** | `extract_combine_rec` |
| `verify_product` | **Pipeline private** | `verify_product` |
| `try_zassenhaus_factor` | **Partial** | Zassenhaus+Hensel lift |
| `cubic1` | **Pipeline private** | `cubic1` |
| `cubic2` | **Pipeline private** | `cubic2` |
| `zassenhaus_two_cubics` | **Pipeline private** | `zassenhaus_two_cubics` |
| `zassenhaus_sophie_germain_quartic` | **Pipeline private** | `zassenhaus_sophie_germain_quartic` |
| `zassenhaus_quadratic_times_cubic` | **Pipeline private** | `zassenhaus_quadratic_times_cubic` |

### `modint.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `new` | **Stable** | `new` |
| `from_i64` | **Stable** | `Poly::from_i64` |
| `add` | **Stable** | Poly addition |
| `sub` | **Stable** | Poly subtraction |
| `mul` | **Stable** | Poly multiplication |
| `inv` | **Stable** | `inv` |
| `is_zero` | **Stable** | Poly is zero |
| `is_one` | **Stable** | Poly is one |
| `smod` | **Stable** | symmetric mod for i64 |
| `irem` | **Stable** | integer remainder |
| `smod_positive` | **Pipeline private** | `smod_positive` |

### `modular.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `zero` | **Stable** | Poly zero |
| `one` | **Stable** | Poly one |
| `from_poly` | **Stable** | `from_poly` |
| `is_zero` | **Stable** | Poly is zero |
| `leading_term` | **Stable** | leading term by total degree |
| `add` | **Stable** | Poly addition |
| `sub` | **Stable** | Poly subtraction |
| `mul` | **Stable** | Poly multiplication |
| `div_rem` | **Stable** | multivariate division with remainder |
| `gcd` | **Stable** | Poly gcd via subresultant |
| `modp` | **Stable** | Poly → PolyMod mod p |

### `monomial.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `one` | **Stable** | Poly one |
| `var` | **Stable** | Poly univariate generator |
| `degree` | **Stable** | total degree |
| `exp_of` | **Stable** | `Poly::exp_of` |
| `is_const` | **Stable** | `is_const` |
| `iter` | **Stable** | `iter` |
| `mul` | **Stable** | Poly multiplication |
| `div_exact` | **Stable** | exact division if remainder zero |
| `is_dividing` | **Stable** | `is_dividing` |
| `cmp_lex` | **Stable** | `cmp_lex` |
| `divides` | **Stable** | `divides` |
| `lcm` | **Stable** | Poly lcm |

### `nested.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `new` | **Stable** | wrap a variable name as main indeterminate. |
| `as_var` | **Stable** | underlying [`Var`]. |
| `from` | **Stable** | `Poly::from` |
| `new` | **Stable** | view `p` as a coefficient-ring element. |
| `gcd` | **Stable** | gcd in ℚ[others] via subresultant (typed coefficient-ring path). |
| `exact_quo` | **Stable** | exact quotient in ℚ[others]; `None` if `den` does not divide `self`. |
| `new` | **Stable** | view `poly` as univariate in `main`. |
| `main_var` | **Stable** | main variable. |
| `degree` | **Stable** | degree w.r.t. `main` (coefficient-ring agnostic). |
| `coeff_at` | **Stable** | coefficient of `main^exp` as ℚ[others]. |
| `leading_coeff` | **Stable** | leading coefficient w.r.t. `main` (element of ℚ[others]). |
| `divides` | **Stable** | whether `self` divides `p` in ℚ[others][main]. |
| `exact_quo_dividing` | **Stable** | exact quotient `p / self` in ℚ[others][main]. |
| `div_rem_wrt_aux_indep` | **Stable** | `div_rem_wrt_aux_indep` |
| `eval_aux` | **Stable** | substitute `aux ↦ value` in ℚ[others][main]; `main` unchanged. |
| `new` | **Stable** | construct owned nested-ring element. |
| `as_view` | **Stable** | borrowed view. |
| `new` | **Stable** | view `poly` as univariate in `var`. |
| `as_poly` | **Stable** | underlying polynomial. |
| `into_poly` | **Stable** | consume and return inner [`Poly`]. |
| `var` | **Stable** | main variable. |
| `degree` | **Stable** | degree w.r.t. main variable (coefficient-ring agnostic). |
| `try_new` | **Stable** | construct when `poly` is genuinely univariate in `var` over ℚ. |
| `div_rem` | **Stable** | Euclidean `(q, r)` in ℚ[var] via [`univariate_div_rem_wrt`]. |
| `divides` | **Stable** | whether `divisor` divides `self` in ℚ[var]. |
| `exact_quo` | **Stable** | exact quotient `self / divisor` when remainder is zero. |
| `new` | **Stable** | wrap a general sparse polynomial. |
| `as_poly` | **Stable** | underlying [`Poly`] (explicit downgrade). |
| `into_inner` | **Stable** | consume and return inner [`Poly`]. |
| `div_rem` | **Stable** | multivariate leading-monomial division (display / heuristic contexts only). |
| `wrt` | **Stable** | `p / content_wrt(p, var)` in ℚ[others][var]. |
| `as_poly` | **Stable** | primitive part polynomial. |
| `main_var` | **Stable** | variable the part is primitive w.r.t. |
| `into_inner` | **Stable** | consume into `(poly, main)`. |
| `pair` | **Stable** | `Poly::pair` |
| `apply` | **Stable** | `Poly::apply` |
| `undo` | **Stable** | `undo` |
| `is_flat_univariate_in` | **Pipeline private** | `is_flat_univariate_in` |
| `substitute_wrt` | **Pipeline private** | substitute `sub_var ↦ sub_poly` (coefficients may depend on other vars). |
| `dilate_one` | **Pipeline private** | `dilate_one` |
| `undilate_one` | **Pipeline private** | `undilate_one` |
| `dilate_aux` | **Stable** | `dilate_aux` |
| `undilate_aux` | **Stable** | `undilate_aux` |
| `unit` | **Stable** | `n = [1, 1]` embed with fresh `t` variable. |
| `with_n` | **Stable** | copy with new exponent vector. |
| `embed` | **Stable** | map `p ∈ ℚ[main, aux…]` into ℚ[main, t] (typed embed result). |
| `view_main` | **Stable** | view embedded poly as univariate in `main` (coeffs in ℚ[t]). |
| `view_t` | **Stable** | view embedded poly as univariate in `t`. |
| `new` | **Stable** | construct from polynomial and embed configuration. |
| `as_poly` | **Stable** | underlying polynomial (explicit downgrade). |
| `into_poly` | **Stable** | consume and return the inner [`Poly`]. |
| `config` | **Stable** | embed configuration. |
| `view_main` | **Stable** | view as univariate in `main` (coeffs in ℚ[t]). |
| `view_t` | **Stable** | view as univariate in `t`. |
| `to_recon_draft` | **Stable** | `to_recon_draft` |
| `from_embedded` | **Stable** | `Poly::from_embedded` |
| `from_poly` | **Stable** | `Poly::from_poly` |
| `materialize` | **Stable** | `materialize` |
| `from_eval_factor` | **Stable** | `Poly::from_eval_factor` |
| `new` | **Stable** | `new` |
| `as_univariate_in` | **Stable** | `Poly::as_univariate_in` |
| `is_independent_of_var` | **Stable** | `is_independent_of_var` |
| `div_rem_wrt_aux_indep` | **Stable** | `div_rem_wrt_aux_indep` |
| `sorted_from_poly` | **Stable** | `Poly::sorted_from_poly` |
| `embed_monomials_to_poly` | **Pipeline private** | build `Poly` from sorted embed monomials + aux exponents |
| `term_with_var` | **Pipeline private** | shared with `poly_uni::term_with_var` |
| `coeff_wrt_impl` | **Pipeline private** | `coeff_wrt_impl` |
| `eval_tn_impl` | **Pipeline private** | `eval_tn` |
| `flat_uni_div_rem_vs_multivariate_div_rem` | **Pipeline private** | `flat_uni_div_rem_vs_multivariate_div_rem` |
| `primitive_part_wrt_tags_main` | **Pipeline private** | `primitive_part_wrt_tags_main` |
| `dilation_map_apply_undo` | **Pipeline private** | `dilation_map_apply_undo` |
| `multivariate_poly_div_rem_extract_powers` | **Pipeline private** | `multivariate_poly_div_rem_extract_powers` |
| `coeff_ring_gcd_divides_both` | **Pipeline private** | `coeff_ring_gcd_divides_both` |
| `univariate_in_divides_vs_div_rem` | **Pipeline private** | `univariate_in_divides_vs_div_rem` |
| `lifted_factor_peel_vs_div_rem` | **Pipeline private** | `lifted_factor_peel_vs_div_rem` |
| `univariate_in_generic_degree_matches_q` | **Pipeline private** | `univariate_in_generic_degree_matches_q` |
| `flat_uni_generic_degree` | **Pipeline private** | `flat_uni_generic_degree` |
| `coeff_ring_exact_quo` | **Pipeline private** | `coeff_ring_exact_quo` |
| `eval_aux_preserves_main` | **Pipeline private** | `eval_aux_preserves_main` |
| `tn_embed_maps_aux_to_t` | **Pipeline private** | `tn_embed_maps_aux_to_t` |
| `embed_factor_draft_materialize_with_aux_exps` | **Pipeline private** | `embed_factor_draft_materialize_with_aux_exps` |
| `div_rem_wrt_aux_indep_matches_univariate` | **Pipeline private** | `div_rem_wrt_aux_indep_matches_univariate` |

### `ops.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `gauss` | **Stable** | Gauss elimination on Poly rows |
| `content` | **Stable** | integer content of Poly |

### `partfrac.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `partfrac_terms` | **Stable (bounded)** | partial fraction terms |
| `partfrac_rational_terms` | **Stable (bounded)** | partfrac with poly part |
| `drop_zero_numerators` | **Pipeline private** | `drop_zero_numerators` |
| `partfrac_by_square_free` | **Pipeline private** | `partfrac_by_square_free` |
| `partfrac_affine_power_system` | **Pipeline private** | `partfrac_affine_power_system` |
| `denominator_power_factors` | **Pipeline private** | `denominator_power_factors` |
| `expand_sqff_factors` | **Pipeline private** | `expand_sqff_factors` |
| `partfrac_square_free_affine_numerators` | **Pipeline private** | `partfrac_square_free_affine_numerators` |
| `partfrac_one_quadratic` | **Pipeline private** | `partfrac_one_quadratic` |
| `solve_linear_system` | **Pipeline private** | `solve_linear_system` |
| `partfrac_mixed_affine` | **Pipeline private** | `partfrac_mixed_affine` |
| `linear_root` | **Pipeline private** | `linear_root` |
| `x` | **Pipeline private** | `x` |
| `find_rational_root_on_x_plus_one_times_x_fourth_minus_one` | **Pipeline private** | `find_rational_root_on_x_plus_one_times_x_fourth_minus_one` |
| `find_rational_root_on_x_plus_one_sq_times_x_sq_plus_one` | **Pipeline private** | `find_rational_root_on_x_plus_one_sq_times_x_sq_plus_one` |
| `factor_by_roots_x_plus_one_times_x_fourth_minus_one` | **Pipeline private** | `factor_by_roots_x_plus_one_times_x_fourth_minus_one` |
| `denominator_factors_x_plus_one_times_x_fourth_minus_one` | **Pipeline private** | `denominator_factors_x_plus_one_times_x_fourth_minus_one` |
| `partfrac_cubic_irreducible_denominator` | **Pipeline private** | `partfrac_cubic_irreducible_denominator` |
| `partfrac_three_quarters_over_x_fourth_minus_one` | **Pipeline private** | `partfrac_three_quarters_over_x_fourth_minus_one` |
| `partfrac_by_sqff_x_over_x_plus_one_times_x_fourth_minus_one` | **Pipeline private** | `partfrac_by_sqff_x_over_x_plus_one_times_x_fourth_minus_one` |
| `partfrac_x_over_x_plus_one_times_x_fourth_minus_one` | **Pipeline private** | `partfrac_x_over_x_plus_one_times_x_fourth_minus_one` |
| `partfrac_one_over_x_squared_minus_one` | **Pipeline private** | `partfrac_one_over_x_squared_minus_one` |
| `partfrac_x_over_repeated_linear` | **Pipeline private** | `partfrac_x_over_repeated_linear` |
| `partfrac_one_over_x_squared_minus_one_squared` | **Pipeline private** | `partfrac_one_over_x_squared_minus_one_squared` |
| `partfrac_biquadratic_half_angle_denominator` | **Pipeline private** | `partfrac_biquadratic_half_angle_denominator` |
| `partfrac_one_over_x_times_x_squared_plus_one` | **Pipeline private** | `partfrac_one_over_x_times_x_squared_plus_one` |
| `partfrac_xplus1_over_x_squared_minus_one` | **Pipeline private** | `partfrac_xplus1_over_x_squared_minus_one` |
| `partfrac_ck_int_05_denominator` | **Pipeline private** | `partfrac_ck_int_05_denominator` |
| `partfrac_one_over_x_minus_one_squared` | **Pipeline private** | `partfrac_one_over_x_minus_one_squared` |
| `partfrac_x_over_x_minus_one_squared` | **Pipeline private** | `partfrac_x_over_x_minus_one_squared` |
| `partfrac_one_over_x_squared_minus_four` | **Pipeline private** | `partfrac_one_over_x_squared_minus_four` |
| `partfrac_one_over_one_minus_x_squared` | **Pipeline private** | `partfrac_one_over_one_minus_x_squared` |

### `poly.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `ring_zero` | **Stable** | zero polynomial in ring `C`. |
| `ring_one` | **Stable** | unit polynomial in ring `C`. |
| `ring_constant` | **Stable** | constant polynomial in ring `C`. |
| `ring_var` | **Stable** | univariate generator in ring `C`. |
| `is_zero` | **Stable** | Poly is zero |
| `is_one` | **Stable** | Poly is one |
| `leading_term` | **Stable** | leading term by total degree |
| `leading_term_lex` | **Stable** | leading term with variable order |
| `term` | **Stable** | monomial × coefficient |
| `degree` | **Stable** | total degree |
| `degree_wrt` | **Stable** | univariate degree in `var` (independent of coefficient ring `C`). |
| `try_add` | **Stable** | fallible addition (required for non-ℚ coefficients). |
| `try_sub` | **Stable** | fallible subtraction |
| `try_neg` | **Stable** | fallible negation |
| `try_mul` | **Stable** | fallible multiplication |
| `try_pow` | **Stable** | fallible integer power |
| `zero` | **Stable** | Poly zero (ℚ) |
| `one` | **Stable** | Poly one (ℚ) |
| `constant` | **Stable** | Poly scalar constant (ℚ) |
| `var` | **Stable** | Poly univariate generator (ℚ) |
| `add` | **Stable** | Poly addition (ℚ) |
| `sub` | **Stable** | Poly subtraction |
| `neg` | **Stable** | Poly negation |
| `mul` | **Stable** | Poly multiplication |
| `mul_scalar` | **Stable** | scale Poly by rational |
| `pow` | **Stable** | Poly integer power |
| `content` | **Stable** | integer content of Poly |
| `primitive_part` | **Stable** | divide out content |
| `monic` | **Stable** | divide by leading coeff |
| `div_rem` | **Stable** | multivariate division with remainder |
| `div_exact` | **Stable** | exact division if remainder zero |
| `gcd` | **Stable** | Poly gcd via subresultant |
| `lcm` | **Stable** | Poly lcm |
| `horner` | **Stable** | Horner eval at rational point |
| `integer_content_gcd` | **Pipeline private** | gcd of rational coeffs as Ratio |
| `quo` | **Stable** | exact quotient Poly/ Poly |
| `rem` | **Stable** | remainder Poly/ Poly |
| `egcd` | **Stable** | extended gcd (s,t,g) |
| `simp2` | **Stable** | reduce fraction pair by gcd |
| `abcuv` | **Stable** | Bezout coeffs for au+bv=c |
| `x` | **Pipeline private** | `x` |
| `gcd_x3_x2` | **Pipeline private** | `gcd_x3_x2` |
| `quo_rem` | **Pipeline private** | `quo_rem` |
| `abcuv_linear_one` | **Pipeline private** | `abcuv_linear_one` |
| `generic_poly_try_add_matches_q` | **Pipeline private** | `generic_poly_try_add_matches_q` |

### `poly_coeff.rs`

| Function / trait | Tier | Ring | Description |
|------------------|------|------|-------------|
| `PolyCoeff` | **Stable** | 环 | sparse 系数环 |
| `FieldCoeff` | **Stable** | **域 K** | flat Euclidean / gcd；impl: `Ratio<BigInt>`, `AlgExtCPolyCoeff` |
| `coeff_inv`, `field_div` | **Stable** | 域 | default via `coeff_div` |
| `coeff_zero` | **Pipeline private** | — | `coeff_zero` |
| `coeff_one` | **Pipeline private** | — | `coeff_one` |
| `coeff_is_zero` | **Pipeline private** | — | `coeff_is_zero` |
| `coeff_is_one` | **Pipeline private** | — | `coeff_is_one` |
| `coeff_add` | **Pipeline private** | — | `coeff_add` |
| `coeff_sub` | **Pipeline private** | — | `coeff_sub` |
| `coeff_neg` | **Pipeline private** | — | `coeff_neg` |
| `coeff_mul` | **Pipeline private** | — | `coeff_mul` |
| `coeff_div` | **Pipeline private** | — | `coeff_div` |
| `coeff_zero` | **Stable** | — | `Poly::coeff_zero` |
| `coeff_one` | **Stable** | — | `Poly::coeff_one` |
| `coeff_is_zero` | **Stable** | — | `Poly::coeff_is_zero` |
| `coeff_is_one` | **Stable** | — | `Poly::coeff_is_one` |
| `ratio_coeff_ring` | **Pipeline private** | — | `ratio_coeff_ring` |

### `univ_wrt.rs`

| Function | Tier | Ring | Description |
|----------|------|------|-------------|
| `is_univariate_in` | **Stable** | — | only powers of `var` |
| `scalar_coeff_wrt` | **Stable** | K | coefficient of `var^exp` |
| `univariate_div_rem_wrt` | **Stable** | **K[var]** | Euclidean `(q,r)`; `d=0`/`lc=0` → Err |
| `quo_exact_wrt` | **Stable** | **K[var]** | exact quotient |
| `derivative_wrt` | **Stable** | **K[var]** | formal derivative |
| `egcd_wrt` | **Stable** | **K[var]** | extended gcd (monic) |
| `gcd_wrt` | **Stable** | **K[var]** | gcd (monic) |
| `content_scalars` | **Stable** | K | scalar content |
| `content_wrt` | **Stable** | **K[var]** | content w.r.t. var |
| `primitive_part_wrt` | **Stable** | **K[var]** | primitive part |
| `square_free_part_wrt` | **Stable** | **K[var]** | square-free part (Yun) |
| `quadratic_coeffs_wrt` | **Stable** | **K[var]** | `(a,b,c)` for deg-2 |
| `monic_wrt` | **Stable** | **K[var]** | monic normalize |
| `debug_assert_euclidean_post` | **Pipeline private** | — | debug-only post check |

### `resultant.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `resultant` | **Stable** | univariate resultant |
| `univariate_coeffs_ascending` | **Stable** | ascending univariate coefficient vector in ℚ[var]. |
| `univariate_coefficients` | **Pipeline private** | `univariate_coefficients` |
| `sylvester_det` | **Pipeline private** | `sylvester_det` |
| `det_rational` | **Pipeline private** | `det_rational` |
| `sylvester_det2` | **Pipeline private** | `sylvester_det2` |
| `coeff_at` | **Stable** | univariate coefficient at exponent |
| `univariate_degree` | **Stable** | degree w.r.t. var |
| `univariate_leading_coeff` | **Pipeline private** | `univariate_leading_coeff` |
| `roots` | **Stable (bounded)** | low-degree exact roots as Poly factors |
| `quadratic_abc` | **Pipeline private** | `quadratic_abc` |
| `quadratic_coeffs` | **Pipeline private** | `quadratic_coeffs` |
| `ratio_is_perfect_square` | **Pipeline private** | `ratio_is_perfect_square` |
| `int_isqrt` | **Pipeline private** | `int_isqrt` |
| `quadratic_roots` | **Pipeline private** | `quadratic_roots` |
| `is_xn_minus_one` | **Pipeline private** | `is_xn_minus_one` |
| `x` | **Pipeline private** | `x` |
| `resultant_shared_factor_is_zero` | **Pipeline private** | `resultant_shared_factor_is_zero` |
| `resultant_two_linear_polys` | **Pipeline private** | `resultant_two_linear_polys` |
| `resultant_constant_times_linear` | **Pipeline private** | `resultant_constant_times_linear` |
| `resultant_linear_times_constant` | **Pipeline private** | `resultant_linear_times_constant` |
| `resultant_quadratic` | **Pipeline private** | `resultant_quadratic` |
| `quadratic_abc_x2_minus_2` | **Pipeline private** | `quadratic_abc_x2_minus_2` |
| `roots_linear` | **Pipeline private** | `roots_linear` |
| `roots_zero_polynomial` | **Pipeline private** | `roots_zero_polynomial` |
| `roots_constant_nonzero_errors` | **Pipeline private** | `roots_constant_nonzero_errors` |
| `roots_x3_minus_one_real_root` | **Pipeline private** | `roots_x3_minus_one_real_root` |
| `roots_quadratic_perfect_square` | **Pipeline private** | `roots_quadratic_perfect_square` |
| `roots_quadratic_irrational_discriminant` | **Pipeline private** | `roots_quadratic_irrational_discriminant` |

### `subresultant.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `coeff_wrt` | **Stable** | coefficient Poly w.r.t. var^exp |
| `monomial_pow` | **Pipeline private** | `monomial_pow` |
| `term_with_var` | **Stable** | coeff * var^exp as Poly |
| `vars_in` | **Stable** | sorted variables in Poly |
| `vars_union` | **Pipeline private** | `vars_union` |
| `is_univariate_in` | **Pipeline private** | `is_univariate_in` |
| `main_var_for_gcd` | **Pipeline private** | `main_var_for_gcd` |
| `rational_primitive` | **Pipeline private** | `rational_primitive` |
| `div_exact_coeff` | **Stable** | `div_exact_coeff` |
| `quo_exact_coeff` | **Stable** | `quo_exact_coeff` |
| `univariate_div_rem_wrt` | **Stable** | `univariate_div_rem_wrt` |
| `quo_exact_wrt` | **Stable** | `quo_exact_wrt` |
| `pseudo_rem_wrt` | **Pipeline private** | `pseudo_rem_wrt` |
| `content_wrt` | **Stable** | content w.r.t. main var |
| `primitive_part_wrt` | **Stable** | primitive part w.r.t. var |
| `content_wrt_impl` | **Stable** | `content_wrt_impl` |
| `primitive_part_wrt_impl` | **Stable** | `primitive_part_wrt_impl` |
| `gcd_constant_wrt` | **Pipeline private** | `gcd_constant_wrt` |
| `subresultant_gcd_wrt` | **Pipeline private** | `subresultant_gcd_wrt` |
| `subresultant_gcd` | **Stable** | multivariate gcd subresultant |
| `scale_gcd_by_content` | **Pipeline private** | `scale_gcd_by_content` |
| `gcd_xy_and_y` | **Pipeline private** | `gcd_xy_and_y` |
| `gcd_bivariate_linear_and_quadratic` | **Pipeline private** | `gcd_bivariate_linear_and_quadratic` |
| `gcd_univariate_matches_subresultant` | **Pipeline private** | `gcd_univariate_matches_subresultant` |
| `content_wrt_y_of_xy` | **Pipeline private** | `content_wrt_y_of_xy` |

### `tresultant.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `num_minus_t_derivative` | **Stable** | RT numerator derivative |
| `embed_univariate_x` | **Pipeline private** | `embed_univariate_x` |
| `tresultant_eliminate_x` | **Stable** | eliminate x via t-resultant |
| `lagrange_poly` | **Pipeline private** | `lagrange_poly` |
| `eval_param_poly` | **Stable** | substitute parameter in Poly |
| `rational_roots_in_t` | **Stable** | rational roots in parameter t |
| `eval_univariate_at_t` | **Pipeline private** | `eval_univariate_at_t` |
| `rational_root_candidates` | **Pipeline private** | `rational_root_candidates` |
| `push_divisors` | **Pipeline private** | `push_divisors` |
| `biquadratic_res_conjugate_pairs` | **Partial** | RT biquadratic resolvent |
| `biquartic_conjugate_pairs` | **Partial** | RT biquartic resolvent |
| `pairs_from_symmetric_res_factors` | **Pipeline private** | `pairs_from_symmetric_res_factors` |
| `sqrt_rational_coeff_radicand` | **Pipeline private** | `sqrt_rational_coeff_radicand` |
| `extract_sqrt_factor` | **Pipeline private** | `extract_sqrt_factor` |
| `pure_biquartic_res_pairs` | **Pipeline private** | `pure_biquartic_res_pairs` |
| `ratio_perfect_sqrt` | **Stable** | detect perfect square Ratio |
| `integer_perfect_sqrt` | **Pipeline private** | `integer_perfect_sqrt` |
| `x` | **Pipeline private** | `x` |
| `t` | **Pipeline private** | `t` |
| `biquartic_pairs_one_over_x_fourth_plus_one` | **Pipeline private** | `biquartic_pairs_one_over_x_fourth_plus_one` |
| `biquadratic_pairs_x_fourth_plus_four_res` | **Pipeline private** | `biquadratic_pairs_x_fourth_plus_four_res` |
| `biquadratic_pairs_x_fourth_plus_x_squared_plus_one_res` | **Pipeline private** | `biquadratic_pairs_x_fourth_plus_x_squared_plus_one_res` |
| `tresultant_one_over_x_squared_plus_one` | **Pipeline private** | `tresultant_one_over_x_squared_plus_one` |
| `tresultant_one_over_x_fourth_plus_one` | **Pipeline private** | `tresultant_one_over_x_fourth_plus_one` |

### `univariate.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `univariate_coeffs` | **Pipeline private** | `univariate_coeffs` |
| `trim_coeffs` | **Pipeline private** | `trim_coeffs` |
| `poly_from_coeffs` | **Pipeline private** | `poly_from_coeffs` |
| `univariate_div_rem` | **Pipeline private** | `univariate_div_rem` |
| `coeffs_to_integer_primitive` | **Pipeline private** | `coeffs_to_integer_primitive` |
| `is_zero_int` | **Pipeline private** | `is_zero_int` |
| `int_exact_div_rem` | **Pipeline private** | `int_exact_div_rem` |
| `pseudo_remainder` | **Pipeline private** | `pseudo_remainder` |
| `trim_int` | **Pipeline private** | `trim_int` |
| `poly_from_int_coeffs` | **Pipeline private** | `poly_from_int_coeffs` |
| `monic_univariate` | **Pipeline private** | `monic_univariate` |
| `univariate_derivative` | **Stable** | derivative w.r.t. var |
| `square_free_factorization` | **Stable** | Yun square-free factors |
| `square_free_part` | **Stable** | product of square-free factors |
| `substitute_univariate` | **Stable** | substitute var → Poly |
| `odd_multiplicity_part` | **Stable** | odd multiplicity factor |
| `odd_part_core` | **Pipeline private** | `odd_part_core` |
| `gcd_univariate` | **Stable** | univariate gcd |
| `univariate_gcd` | **Pipeline private** | `univariate_gcd` |
| `univariate_div_exact` | **Pipeline private** | `univariate_div_exact` |
| `gcd_reduce` | **Pipeline private** | `gcd_reduce` |
| `univariate_rem` | **Pipeline private** | `univariate_rem` |
| `sturm_sequence` | **Stable** | Sturm chain |
| `eval_univariate_at` | **Stable** | Horner eval |
| `sign_of_ratio` | **Pipeline private** | `sign_of_ratio` |
| `sign_variations` | **Stable** | sign change count in sequence |
| `sturm_sign_variations_at` | **Stable** | Sturm sign count at point |
| `sturmab_count` | **Stable** | root count in (a,b) |
| `x` | **Pipeline private** | `x` |
| `sturm_x_cubed_plus_one_has_three_polys` | **Pipeline private** | `sturm_x_cubed_plus_one_has_three_polys` |
| `square_free_x_squared_times_cubic` | **Pipeline private** | `square_free_x_squared_times_cubic` |
| `sturmab_x_squared_times_cubic` | **Pipeline private** | `sturmab_x_squared_times_cubic` |
