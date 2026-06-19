# giac-poly API 稳定性分层

**规范来源:** [algorithm-expr-api.md](algorithm-expr-api.md)  
**上游缺口:** [issues/GIAC-simplify-poly-upstream-gaps.md](issues/GIAC-simplify-poly-upstream-gaps.md) §2  
**Expr ↔ Poly 边界:** 仅 **`giac-core::expr_to_poly` / `poly_to_expr`** 做跨表示转换；本 crate 全程 `Poly`。

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

| API | 语义 |
|-----|------|
| `Poly::div_rem` | 多元 leading-monomial 除法 |
| `subresultant::univariate_div_rem_wrt` | ℚ[others][var] 上一元除法 |
| `subresultant::quo_exact_wrt` | 上一元除法，余式非零 → `Err` |

**命名:** 公开 `try_*`（如 `try_hensel_lift_bivariate`）表示 **可选算法路径**，非 [algorithm-expr-api](algorithm-expr-api.md) 意义的临时 `shim_*`；失败时静默 `None`，调用方须处理。

**提交前复审（测试全绿后）：** 见 [algorithm-expr-api.md §6.2](algorithm-expr-api.md#62-测试通过后提交--合入前复审) — 检查 factor 等模块临时 fallback 是否净减少、新增 `fn` tier 是否已更新本文 Per-file 表。

---

## 2. Crate 公开 API（`lib.rs` re-export）

### 2.1 核心类型与环运算 — **Stable**

| 符号 | 模块 | 说明 |
|------|------|------|
| `Poly`, `Monomial`, `Var` | `poly`, `monomial` | 多项式表示 |
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
| `try_factor_patterns` | `factor/patterns` | **Partial** — cyclotomic/二项式模式 |
| `try_factor_xn_minus_one` 等 | `factor/cyclotomic` | **Partial** |
| `factor_fpx`, `degree` | `factor/fpx` | **Stable**（模域） |
| `normalize_univariate_factors`, `try_lift_factors_in_aux_var` | `factor/hensel` | **Pipeline private** `pub(crate)` |
| `find_rational_root` | `factor/univariate` | **Pipeline private** `pub(crate)` |
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
| **FAC-G1** | `try_sparse_factor` + `try_sparse_factor_bi` | **Partial** — 好点种子 + 2-aux MVP；sum-coeff / dilation / pzadic 待补 |
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

### `exp.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `bigint_pow` | **Pipeline private** | BigInt pow with overflow check |

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

### `factor/fpx.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `x_var` | **Pipeline private** | `x_var` |
| `mi` | **Pipeline private** | `mi` |
| `is_poly_one` | **Pipeline private** | `is_poly_one` |
| `degree` | **Stable** | total degree |
| `coeff` | **Pipeline private** | `coeff` |
| `set_coeff` | **Pipeline private** | `set_coeff` |
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
| `is_independent_of_y` | **Pipeline private** | `is_independent_of_y` |
| `leading_coeff_x` | **Pipeline private** | `leading_coeff_x` |
| `scale_univariate_x` | **Pipeline private** | `scale_univariate_x` |
| `div_rem_x_over_qy` | **Pipeline private** | `div_rem_x_over_qy` |
| `hensel_lift_two_at_zero` | **Pipeline private** | `hensel_lift_two_at_zero` |
| `normalize_univariate_factors` | **Pipeline private** | `normalize_univariate_factors` |
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
| `try_hensel_lift_bivariate` | **Partial** | Factor `p(x,y)` in ℚ[y][x]: Hensel lift at `y=0`, then interpolation fallback. |
| `hensel_three_linear_shifted` | **Pipeline private** | `hensel_three_linear_shifted` |
| `hensel_two_bilinear_factors` | **Pipeline private** | `hensel_two_bilinear_factors` |
| `hensel_egcd_two_factors` | **Pipeline private** | `hensel_egcd_two_factors` |
| `hensel_egcd_factor_list` | **Pipeline private** | `hensel_egcd_factor_list` |
| `div_rem_xy_by_x_minus_one` | **Pipeline private** | `div_rem_xy_by_x_minus_one` |
| `hensel_two_bilinear_at_zero` | **Pipeline private** | `hensel_two_bilinear_at_zero` |
| `factor_non_monic_product_at_zero` | **Pipeline private** | `factor_non_monic_product_at_zero` |
| `hensel_two_simple_non_monic` | **Pipeline private** | `hensel_two_simple_non_monic` |
| `aux_lift_line21` | **Pipeline private** | `aux_lift_line21` |
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
| `square_free_wrt_impl` | **Pipeline private** | Yun sqff loop; uses `quo_exact_wrt` |
| `substitute_poly` | **Stable** | Substitute `sub_var -> sub_poly` in `p`. |
| `factor_sqff_over_coeff_ring` | **Partial** | Factor square-free `g` in ℚ[others][var] recursively. |
| `try_factor_bivariate_eval` | **Temporary (retired)** | `#[cfg(test)]` eval+interp lift |
| `try_lift_bivariate_from_eval` | **Temporary (retired)** | `#[cfg(test)]` |
| `lift_univariate_factor` | **Temporary (retired)** | `#[cfg(test)]` |
| `try_kronecker_bivariate` | **Temporary (retired)** | `#[cfg(test)]` Kronecker embed |
| `kronecker_lift` | **Pipeline private** | Kronecker coeff decode (test-only caller) |
| `as_constant` | **Pipeline private** | `as_constant` |
| `as_constant` | **Stable** | `Poly::as_constant` |
| `content_wrt_xy_plus_y_squared` | **Pipeline private** | `content_wrt_xy_plus_y_squared` |
| `rational_content_vs_wrt_content` | **Pipeline private** | `rational_content_vs_wrt_content` |

### `factor/power.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `as_perfect_power` | **Stable** | If `p` is a perfect power, return `(base, exponent)`. |
| `try_nth_root` | **Pipeline private** | optional fallback `try_nth_root` |
| `try_binomial_square` | **Pipeline private** | optional fallback `try_binomial_square` |
| `try_linear_power` | **Partial** | detect (linear)^n |

### `factor/sqrt.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `quadratic_sqrt_factor_exprs` | **Partial** | Display factors `(x - (-b ± sqrt(disc))/(2a))` for a univariate quadratic. |
| `format_ratio` | **Pipeline private** | `format_ratio` |
| `format_sqrt_ratio` | **Pipeline private** | `format_sqrt_ratio` |
| `sqrt_factor_x2_minus_2` | **Pipeline private** | `sqrt_factor_x2_minus_2` |

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
| `coeff_wrt` | **Stable** | Coefficient of `var^exp` (quotient by `var^exp` on each matching term). |
| `monomial_pow` | **Pipeline private** | `monomial_pow` |
| `is_monic_univariate` | **Stable** | `is_monic_univariate` |

### `factor/zassenhaus.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `x_var` | **Pipeline private** | `x_var` |
| `monomial_x_pow` | **Pipeline private** | `monomial_x_pow` |
| `is_int_poly` | **Pipeline private** | `is_int_poly` |
| `integer_coeffs` | **Pipeline private** | `integer_coeffs` |
| `mignotte_bound` | **Pipeline private** | `mignotte_bound` |
| `coeff_mod` | **Pipeline private** | `coeff_mod` |
| `poly_mod_from_coeffs` | **Pipeline private** | `poly_mod_from_coeffs` |
| `poly_mod_from_poly` | **Pipeline private** | `poly_mod_from_poly` |
| `make_monic_mod` | **Pipeline private** | `make_monic_mod` |
| `at_modulus` | **Pipeline private** | `at_modulus` |
| `derivative_mod` | **Pipeline private** | `derivative_mod` |
| `is_square_free_mod` | **Pipeline private** | `is_square_free_mod` |
| `extgcd_mod` | **Pipeline private** | `extgcd_mod` |
| `egcd_factor_list_mod` | **Pipeline private** | `egcd_factor_list_mod` |
| `polymod_to_int_poly` | **Pipeline private** | `polymod_to_int_poly` |
| `divides_exact` | **Pipeline private** | `divides_exact` |
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
| `partfrac_one_over_one_minus_x_squared` | **Pipeline private** | `partfrac_one_over_one_minus_x_squared` |

### `poly.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `zero` | **Stable** | Poly zero |
| `one` | **Stable** | Poly one |
| `constant` | **Stable** | Poly scalar constant |
| `var` | **Stable** | Poly univariate generator |
| `is_zero` | **Stable** | Poly is zero |
| `is_one` | **Stable** | Poly is one |
| `leading_term` | **Stable** | leading term by total degree |
| `leading_term_lex` | **Stable** | leading term with variable order |
| `term` | **Stable** | monomial × coefficient |
| `degree` | **Stable** | total degree |
| `add` | **Stable** | Poly addition |
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
| `integer_content_gcd` | **Pipeline private** | `integer_content_gcd` |
| `quo` | **Stable** | exact quotient Poly/ Poly |
| `rem` | **Stable** | remainder Poly/ Poly |
| `egcd` | **Stable** | extended gcd (s,t,g) |
| `simp2` | **Stable** | reduce fraction pair by gcd |
| `abcuv` | **Stable** | Bezout coeffs for au+bv=c |
| `x` | **Pipeline private** | `x` |
| `gcd_x3_x2` | **Pipeline private** | `gcd_x3_x2` |
| `quo_rem` | **Pipeline private** | `quo_rem` |
| `abcuv_linear_one` | **Pipeline private** | `abcuv_linear_one` |

### `resultant.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `resultant` | **Stable** | univariate resultant |
| `univariate_coefficients` | **Pipeline private** | `univariate_coefficients` |
| `sylvester_det` | **Pipeline private** | `sylvester_det` |
| `det_rational` | **Pipeline private** | `det_rational` |
| `sylvester_det2` | **Pipeline private** | `sylvester_det2` |
| `coeff_at` | **Stable** | univariate coefficient at exponent |
| `univariate_degree` | **Stable** | degree w.r.t. var |
| `univariate_leading_coeff` | **Pipeline private** | `univariate_leading_coeff` |
| `roots` | **Stable (bounded)** | low-degree exact roots as Poly factors |
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
| `div_exact_coeff` | **Stable (crate-internal)** | exact quotient in coefficient ring when `b \| a` |
| `quo_exact_coeff` | **Stable (crate-internal)** | `div_exact_coeff` → `PolyResult` |
| `univariate_div_rem_wrt` | **Stable (crate-internal)** | division in ℚ[others][var] |
| `quo_exact_wrt` | **Stable (crate-internal)** | exact quotient w.r.t. `var`; use instead of `Poly::div_rem` for nested rings |
| `pseudo_rem_wrt` | **Pipeline private** | pseudo-remainder w.r.t. `var` |
| `content_wrt` | **Stable** | content w.r.t. main var |
| `primitive_part_wrt` | **Stable** | primitive part w.r.t. var |
| `content_wrt_impl` | **Stable (crate-internal)** | gcd of coefficient polys w.r.t. `var` |
| `primitive_part_wrt_impl` | **Stable (crate-internal)** | primitive part implementation |
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
