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

**提交前复审（测试全绿后）：** 见 [algorithm-expr-api.md §7.2](algorithm-expr-api.md#72-测试通过后提交--合入前复审) — 临时匹配净减少、新增 `fn` tier 已更新本文 §2–§5。

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

**批量标注 / 刷新清单：**

```bash
cd giac-rs
python3 scripts/annotate_api_tiers.py              # 为新 fn 补 tier 注释（幂等）
python3 scripts/annotate_api_tiers.py --inventory # 刷新本文 Per-file 表
```

---

## Per-file function inventory (generated)

Regenerate: `python3 scripts/annotate_api_tiers.py --inventory`

### `diff.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `diff` | **Stable** | symbolic differentiation (GIAC-113 / `giac-calculus`). |
| `diff_mul` | **Pipeline private** | product rule for `Mul`. |
| `diff_pow` | **Pipeline private** | power rule (integer exponent cases). |
| `diff_quotient` | **Pipeline private** | quotient rule via `Frac`. |
| `diff_sin` | **Pipeline private** | chain rule for `sin`. |
| `diff_cos` | **Pipeline private** | chain rule for `cos`. |
| `diff_ln` | **Pipeline private** | chain rule for `ln`. |
| `diff_exp` | **Pipeline private** | chain rule for `exp`. |
| `diff_atan` | **Pipeline private** | chain rule for `atan`. |
| `diff_tan` | **Pipeline private** | chain rule for `tan`. |
| `is_const_wrt` | **Pipeline private** | syntactic constness w.r.t. `var` (local copy). |
| `x` | **Pipeline private** | `x` |
| `xcas` | **Pipeline private** | `xcas` |
| `diff_simplified` | **Pipeline private** | `diff_simplified` |
| `diff_matches` | **Pipeline private** | `diff_matches` |
| `eval_diff_matches` | **Pipeline private** | `eval_diff_matches` |
| `diff_x_squared` | **Pipeline private** | `diff_x_squared` |
| `diff_sin_x_squared` | **Pipeline private** | `diff_sin_x_squared` |
| `diff_ln_times_x_squared` | **Pipeline private** | `diff_ln_times_x_squared` |
| `diff_constant_is_zero` | **Pipeline private** | `diff_constant_is_zero` |
| `diff_atan_x` | **Pipeline private** | `diff_atan_x` |
| `diff_frac_one_over_x` | **Pipeline private** | `diff_frac_one_over_x` |

### `eval_diff.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_diff` | **Stable** | evaluate `diff(expr, var)` or multivariate `derive(expr, [vars…])`. |
| `eval_diff_x_squared` | **Pipeline private** | `eval_diff_x_squared` |
| `eval_derive_multivariate` | **Pipeline private** | `eval_derive_multivariate` |
| `eval_diff_too_few_args` | **Pipeline private** | `eval_diff_too_few_args` |
| `eval_diff_list_multivariate` | **Pipeline private** | `eval_diff_list_multivariate` |
| `eval_diff_bad_variable` | **Pipeline private** | `eval_diff_bad_variable` |

### `eval_integrate.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_integrate` | **Stable** | evaluate `integrate(f, x)` or definite `integrate(f, x, a, b)`. |
| `eval_integrate_definite_one` | **Pipeline private** | `eval_integrate_definite_one` |
| `eval_integrate_inv_x` | **Pipeline private** | `eval_integrate_inv_x` |
| `eval_integrate_wrong_arity` | **Pipeline private** | `eval_integrate_wrong_arity` |
| `eval_integrate_definite_x_squared` | **Pipeline private** | `eval_integrate_definite_x_squared` |

### `expr_util.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `depends_on_var` | **Stable** | whether `e` syntactically depends on `var`. |
| `is_const_wrt` | **Stable** | whether `e` is constant with respect to `var` (syntactic). |
| `depends_on_var_symbol` | **Pipeline private** | `depends_on_var_symbol` |
| `depends_on_var_nested` | **Pipeline private** | `depends_on_var_nested` |
| `depends_on_var_pow_frac_func` | **Pipeline private** | `depends_on_var_pow_frac_func` |

### `integrate.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `integrate` | **Stable** | symbolic integration |
| `integrate_fueled` | **Pipeline private** | fueled recursive core. Public [`integrate`] mints a |
| `integrate_rewrite` | **Pipeline private** | rewrite re-entry: charge one fuel unit then recurse. |
| `try_as_rational` | **Stable** | normalize `Expr` to `(num, den)` rational form. |
| `integrate_frac` | **Stable** | integrate rational `num/den` w.r.t. `var`. |
| `integrate_reciprocal` | **Pipeline private** | integrate reciprocal. |
| `integrate_func` | **Pipeline private** | integrate func. |
| `integrate_sin_squared` | **Pipeline private** | integrate sin squared. |
| `integrate_cos_squared` | **Pipeline private** | integrate cos squared. |
| `integrate_sin_cos_product` | **Pipeline private** | integrate sin cos product. |
| `integrate_x_ln` | **Pipeline private** | integrate x ln. |
| `integrate_mul` | **Pipeline private** | integrate mul. |
| `count_var_factors` | **Pipeline private** | count var factors. |
| `integrate_pow` | **Pipeline private** | integrate pow. |
| `ln_abs` | **Pipeline private** | ln abs. |
| `xcas` | **Pipeline private** | `xcas` |
| `assert_integrate_matches` | **Pipeline private** | `assert_integrate_matches` |
| `eval_integrate_matches` | **Pipeline private** | `eval_integrate_matches` |
| `giac223_tanh_exp_frac` | **Pipeline private** | `giac223_tanh_exp_frac` |
| `giac223_exp_over_linear` | **Pipeline private** | `giac223_exp_over_linear` |
| `integrate_reciprocal` | **Pipeline private** | `integrate_reciprocal` |
| `integrate_one_minus_x_fourth` | **Pipeline private** | `integrate_one_minus_x_fourth` |
| `integrate_constant_and_sum` | **Pipeline private** | `integrate_constant_and_sum` |
| `integrate_const_times_x` | **Pipeline private** | `integrate_const_times_x` |
| `integrate_x_and_x_squared` | **Pipeline private** | `integrate_x_and_x_squared` |
| `integrate_one` | **Pipeline private** | `integrate_one` |
| `integrate_one_over_one_plus_x_squared` | **Pipeline private** | `integrate_one_over_one_plus_x_squared` |
| `integrate_frac_with_constant_numerator` | **Pipeline private** | `integrate_frac_with_constant_numerator` |
| `integrate_one_over_x_squared_plus_one_term_order` | **Pipeline private** | `integrate_one_over_x_squared_plus_one_term_order` |
| `integrate_x_over_x_squared_plus_one` | **Pipeline private** | `integrate_x_over_x_squared_plus_one` |
| `integrate_unsupported_returns_not_implemented` | **Pipeline private** | `integrate_unsupported_returns_not_implemented` |
| `integrate_const_over_quadratic` | **Pipeline private** | `integrate_const_over_quadratic` |
| `integrate_const_times_reciprocal` | **Pipeline private** | `integrate_const_times_reciprocal` |
| `integrate_product_of_consts_only` | **Pipeline private** | `integrate_product_of_consts_only` |
| `integrate_definite_bounds` | **Pipeline private** | `integrate_definite_bounds` |
| `integrate_product_one_plus_x_squared` | **Pipeline private** | `integrate_product_one_plus_x_squared` |
| `integrate_x_over_x_squared_plus_one_squared` | **Pipeline private** | `integrate_x_over_x_squared_plus_one_squared` |
| `integrate_x_squared_reciprocal` | **Pipeline private** | `integrate_x_squared_reciprocal` |
| `integrate_sin_times_cos` | **Pipeline private** | `integrate_sin_times_cos` |
| `integrate_one_over_cos_squared` | **Pipeline private** | `integrate_one_over_cos_squared` |
| `integrate_not_implemented_messages` | **Pipeline private** | `integrate_not_implemented_messages` |
| `integrate_cos_and_sin` | **Pipeline private** | `integrate_cos_and_sin` |
| `integrate_x_cubed` | **Pipeline private** | `integrate_x_cubed` |
| `integrate_exhausted_budget_returns_err` | **Pipeline private** | `integrate_exhausted_budget_returns_err` |
| `integrate_wide_polynomial_sum_does_not_exhaust` | **Pipeline private** | `integrate_wide_polynomial_sum_does_not_exhaust` |
| `integrate_nested_rewrite_terminates_within_budget` | **Pipeline private** | `integrate_nested_rewrite_terminates_within_budget` |

### `integrate_helpers.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `linear_coefficient` | **Pipeline private** | linear coefficient. |
| `affine_var_coeff` | **Pipeline private** | affine var coeff. |
| `is_exp_of_var` | **Stable** | detect `exp(var)` form. |
| `is_x_squared` | **Pipeline private** | is x squared. |
| `is_x_squared_plus_const` | **Pipeline private** | is x squared plus const. |
| `var_coefficient` | **Pipeline private** | var coefficient. |
| `ln_abs_expr` | **Stable** | build `ln(abs(arg))` expression. |
| `is_one` | **Pipeline private** | is one. |

### `integrate_heuristics.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `try_integrate_heuristic` | **Pipeline** | top-level sqrt / trig-fraction hooks before generic `integrate` dispatch. |
| `try_integrate_heuristic_fueled` | **Pipeline private** | fueled core; [`try_integrate_heuristic`] mints a fresh |
| `try_integrate_sin_kx_over_sin_x` | **Partial** | ∫ sin(k·x)/sin(x) dx via Chebyshev U_{k-1}(cos x) (GIAC-225). **退役：** Risch / partfrac. |
| `try_integrate_trig_power_product` | **Partial** | ∫ sin^m(x)·cos^n(x) dx by power reduction when m,n ≥ 1 (GIAC-225). **退役：** Risch / partfrac. |
| `trig_power_exponents` | **Pipeline private** | trig power exponents. |
| `trig_power_factor` | **Pipeline private** | `(is_sin, exponent)` for `sin(x)^n` / `cos(x)^n`. |
| `integrate_sin_sq_cos_4th` | **Pipeline private** | integrate sin sq cos 4th. |
| `expand_trig_power_product` | **Pipeline private** | expand trig power product. |
| `sin_squared_half_angle` | **Pipeline private** | sin squared half angle. |
| `cos_squared_half_angle` | **Pipeline private** | cos squared half angle. |
| `sin_multiple_of_var` | **Pipeline private** | sin multiple of var. |
| `linear_coefficient_int` | **Pipeline private** | linear coefficient int. |
| `chebyshev_u_cos_expr` | **Pipeline private** | Chebyshev U_n(cos x): sin((n+1)x)/sin(x). |
| `try_integrate_trig_deriv_ratio` | **Partial** | ∫ (a·sin + b·cos)'/(a·sin + b·cos) dx = ln\|a·sin + b·cos\| when numerator is d(den)/dx. **退役：** Risch / partfrac. |
| `try_integrate_var_over_sqrt_xsq_plus_c` | **Partial** | ∫ k·x/√(x²+c) dx = √(x²+c) when k=2. **退役：** Risch / partfrac. |
| `try_integrate_x_over_sqrt_affine` | **Partial** | ∫ x/√(x+a) dx = (2/3)(x+a)^(3/2) - 2√(x+a). **退役：** Risch / partfrac. |
| `try_integrate_x_times_sqrt_quadratic` | **Partial** | ∫ x·√(a+x²) dx = (a+x²)^(3/2)/3. **退役：** Risch / partfrac. |
| `try_integrate_sin2x_affine_over_cos2x` | **Partial** | ∫ (k·sin(2x)+c)/cos(2x) dx = −c/(2k)·ln\|c−k·sin(2x)\| (GIAC-normalized). **退役：** Risch / partfrac. |
| `sin_double_angle_expr` | **Pipeline private** | sin double angle expr. |
| `is_cos_sin_double_angle_expr` | **Pipeline private** | is cos sin double angle expr. |
| `sin_double_angle_affine_coeffs` | **Pipeline private** | sin double angle affine coeffs. |
| `expr_contains_sin_or_cos` | **Pipeline private** | detect `sin`/`cos` subexpressions. |
| `try_integrate_trig_rational_half_angle` | **Partial** | heuristic `try_integrate_trig_rational_half_angle`; **退役：** Risch / partfrac. |
| `weierstrass_sin` | **Pipeline private** | weierstrass sin. |
| `weierstrass_cos` | **Pipeline private** | weierstrass cos. |
| `weierstrass_dx_dt` | **Pipeline private** | weierstrass dx dt. |
| `sin_kx_in_t` | **Pipeline private** | sin kx in t. |
| `cos_kx_in_t` | **Pipeline private** | cos kx in t. |
| `replace_trig_with_t` | **Pipeline private** | replace trig with t. |
| `replace_symbol` | **Pipeline private** | replace symbol. |
| `is_trig_rational_in_x` | **Pipeline private** | is trig rational in x. |
| `sin_cos_coeffs` | **Pipeline private** | sin cos coeffs. |
| `trig_affine_parts` | **Pipeline private** | trig affine parts. |
| `sqrt_radicand` | **Pipeline private** | sqrt radicand. |
| `int_const_term` | **Pipeline private** | int const term. |
| `const_term_of_xsq_plus_const` | **Pipeline private** | const term of xsq plus const. |
| `const_term_of_x_plus_const` | **Pipeline private** | const term of x plus const. |
| `is_x_squared` | **Pipeline private** | is x squared. |
| `is_sin_double_angle` | **Pipeline private** | is sin double angle. |
| `is_cos_double_angle` | **Pipeline private** | is cos double angle. |
| `var_coefficient_int` | **Pipeline private** | var coefficient int. |
| `integrate_ck_int_11` | **Pipeline private** | `integrate_ck_int_11` |
| `sin_multiple_of_var_detects_sin_3x` | **Pipeline private** | `sin_multiple_of_var_detects_sin_3x` |
| `integrate_ck_int_43_sin_3x_over_sin_x` | **Pipeline private** | `integrate_ck_int_43_sin_3x_over_sin_x` |
| `integrate_ck_int_14_sin_sq_cos_4th` | **Pipeline private** | `integrate_ck_int_14_sin_sq_cos_4th` |
| `integrate_x_over_sqrt_xsq_plus_one` | **Pipeline private** | `integrate_x_over_sqrt_xsq_plus_one` |
| `integrate_x_times_sqrt_xsq_plus_one` | **Pipeline private** | `integrate_x_times_sqrt_xsq_plus_one` |
| `integrate_x_over_sqrt_x_plus_one` | **Pipeline private** | `integrate_x_over_sqrt_x_plus_one` |

### `integrate_try_rules.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `is_tan_of_var` | **Pipeline private** | is tan of var. |
| `try_integrate_tan_plus_tan_cubed` | **Partial** | heuristic `try_integrate_tan_plus_tan_cubed`; **退役：** Risch / partfrac. |
| `try_integrate_exp_over_linear_exp` | **Partial** | heuristic `try_integrate_exp_over_linear_exp`; **退役：** Risch / partfrac. |
| `try_integrate_exp_over_one_plus_exp2` | **Partial** | heuristic `try_integrate_exp_over_one_plus_exp2`; **退役：** Risch / partfrac. |
| `is_exp_of_double_x` | **Pipeline private** | is exp of double x. |
| `is_exp_x_squared` | **Pipeline private** | is exp x squared. |
| `try_integrate_one_over_cos_squared` | **Partial** | ∫ 1/cos(x)² dx = tan(x). **退役：** Risch / partfrac. |
| `try_integrate_sin_over_cos_sq_frac` | **Partial** | heuristic `try_integrate_sin_over_cos_sq_frac`; **退役：** Risch / partfrac. |
| `const_plus_const_times_exp` | **Pipeline private** | const plus const times exp. |
| `is_sin_double_angle` | **Pipeline private** | is sin double angle. |
| `try_integrate_sin2x_cos` | **Partial** | heuristic `try_integrate_sin2x_cos`; **退役：** Risch / partfrac. |
| `try_integrate_var_shifted_sqrt` | **Partial** | heuristic `try_integrate_var_shifted_sqrt`; **退役：** Risch / partfrac. |
| `var_times_x_squared_plus_const` | **Pipeline private** | var times x squared plus const. |
| `shifted_sqrt_base` | **Pipeline private** | shifted sqrt base. |
| `try_integrate_sin_over_cos_squared` | **Partial** | heuristic `try_integrate_sin_over_cos_squared`; **退役：** Risch / partfrac. |
| `try_integrate_tanh_exp_form` | **Partial** | heuristic `try_integrate_tanh_exp_form`; **退役：** Risch / partfrac. |
| `is_exp_minus_exp_neg` | **Pipeline private** | is exp minus exp neg. |
| `exp_term_sign` | **Pipeline private** | `Some(true)` for `+exp(x)`, `Some(false)` for `-exp(x)` / `exp(-x)` terms. |
| `is_exp_neg_of_var` | **Pipeline private** | is exp neg of var. |
| `is_exp_plus_exp_neg` | **Pipeline private** | is exp plus exp neg. |
| `try_integrate_exp_trig` | **Partial** | heuristic `try_integrate_exp_trig`; **退役：** Risch / partfrac. |
| `integrate_exp_sin` | **Pipeline private** | integrate exp sin. |
| `integrate_exp_cos` | **Pipeline private** | integrate exp cos. |
| `try_integrate_var_over_quadratic_squared` | **Partial** | heuristic `try_integrate_var_over_quadratic_squared`; **退役：** Risch / partfrac. |

### `limit.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_limit` | **Stable** | `limit(expr, var, point)` (GIAC-215). |
| `classify_limit_point` | **Pipeline private** | classify limit point as finite / ±∞. |
| `limit_expr` | **Pipeline private** | dispatch to known limits or `limit_engine` algebraic paths. |
| `try_known_limit` | **Pipeline private** | table lookup for limits needing evaluated point. |
| `try_known_limit_pointless` | **Pipeline private** | table lookup for point-independent classic limits. |
| `limit_finite` | **Pipeline private** | direct substitution at finite point (legacy path). |
| `limit_plus_infinity` | **Pipeline private** | elementary +∞ limits (legacy path). |
| `is_sin_over_x` | **Pipeline private** | detect `sin(x)/x` at 0. |
| `is_var_or_inverse` | **Pipeline private** | `x` or `x^{-1}`. |
| `is_one_minus_cos_over_x_squared` | **Pipeline private** | detect `(1-cos(x))/x^2` at 0. |
| `is_one_minus_cos_frac` | **Pipeline private** | frac form of `(1-cos(x))/x^2`. |
| `is_one_minus_cos_expr` | **Pipeline private** | detect `1 - cos(x)`. |
| `is_var_or_inverse_squared` | **Pipeline private** | `x^2` or `x^{-2}`. |
| `is_one_plus_one_over_x_power_x` | **Pipeline private** | detect `(1+1/x)^x` at +∞. |
| `is_one_plus_reciprocal_var` | **Pipeline private** | detect `1 + 1/x` base. |
| `is_x_squared` | **Pipeline private** | detect `var^2`. |
| `is_zero` | **Pipeline private** | zero test. |
| `is_one` | **Pipeline private** | one test. |
| `is_sin_squared` | **Pipeline private** | detect `sin(var)^2`. |
| `is_x_cubed` | **Pipeline private** | detect `var^3`. |
| `is_ln_one_plus_var` | **Pipeline private** | detect `ln(1+var)`. |
| `is_one_minus_cos_sin2_over_x3_ln1_plus_x` | **Pipeline private** | CK-INT-55 composite limit shape at 0. |
| `is_one_minus_2x_over_quadratic_pole` | **Pipeline private** | CK-INT-59 pole limit shape at 1. |
| `is_one_minus_kx` | **Pipeline private** | detect `1 - k*var`. |
| `int_coeff_times_var` | **Pipeline private** | extract integer coefficient of `var` in a product. |
| `is_x_squared_plus_x_minus_two` | **Pipeline private** | detect `var^2 + var - 2`. |
| `eval_limit_too_few_args` | **Pipeline private** | `eval_limit_too_few_args` |
| `limit_sin_over_x` | **Pipeline private** | `limit_sin_over_x` |
| `limit_ck_int_55` | **Pipeline private** | `limit_ck_int_55` |
| `limit_ck_int_59` | **Pipeline private** | `limit_ck_int_59` |
| `eval_parsed_limit` | **Pipeline private** | `eval_parsed_limit` |
| `limit_ck_int_58_parsed` | **Pipeline private** | `limit_ck_int_58_parsed` |
| `limit_ck_int_60_parsed` | **Pipeline private** | `limit_ck_int_60_parsed` |
| `limit_ck_int_61_parsed` | **Pipeline private** | `limit_ck_int_61_parsed` |
| `assert_limit` | **Pipeline private** | `assert_limit` |
| `rtest_limit_sin_over_x` | **Pipeline private** | `rtest_limit_sin_over_x` |
| `rtest_limit_one_minus_cos_over_x2` | **Pipeline private** | `rtest_limit_one_minus_cos_over_x2` |
| `rtest_limit_wester_one_plus_one_over_n_power_n` | **Pipeline private** | `rtest_limit_wester_one_plus_one_over_n_power_n` |
| `rtest_limit_a_over_n` | **Pipeline private** | `rtest_limit_a_over_n` |
| `rtest_limit_seven_pow_n_over_eight_pow_n` | **Pipeline private** | `rtest_limit_seven_pow_n_over_eight_pow_n` |
| `rtest_limit_four_pow_n_over_two_pow_2n` | **Pipeline private** | `rtest_limit_four_pow_n_over_two_pow_2n` |
| `rtest_limit_x_sqrt_conjugate` | **Pipeline private** | `rtest_limit_x_sqrt_conjugate` |
| `rtest_limit_x_over_x_pow_ln_x` | **Pipeline private** | `rtest_limit_x_over_x_pow_ln_x` |
| `rtest_limit_one_plus_one_over_x_sqrt` | **Pipeline private** | `rtest_limit_one_plus_one_over_x_sqrt` |
| `gruntz_exp_times_exp_diff_minus_one` | **Pipeline private** | `gruntz_exp_times_exp_diff_minus_one` |
| `gruntz_factored_exp_growth_unit` | **Pipeline private** | `gruntz_factored_exp_growth_unit` |
| `gruntz_three_x_five_x_root` | **Pipeline private** | `gruntz_three_x_five_x_root` |
| `gruntz_ck_int_60_ratio` | **Pipeline private** | `gruntz_ck_int_60_ratio` |
| `gruntz_exp_nested_diff` | **Pipeline private** | `gruntz_exp_nested_diff` |
| `rtest_limit_x_atan_x_over_x_plus_1` | **Pipeline private** | `rtest_limit_x_atan_x_over_x_plus_1` |
| `rtest_limit_minus_infinity_inv_x` | **Pipeline private** | `rtest_limit_minus_infinity_inv_x` |
| `rtest_limit_finite_cancel` | **Pipeline private** | `rtest_limit_finite_cancel` |

### `limit_engine/asymptotic.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `limit_at_plus_infinity` | **Pipeline** | `+∞` 极限主编排 |
| `limit_preprocessed_at_plus_infinity` | **Pipeline private** | limit preprocessed at plus infinity |
| `limit_exp_of_vanishing_frac_argument` | **Pipeline private** | limit exp of vanishing frac argument |
| `limit_exp_times_frac_quotient_at_plus_infinity` | **Pipeline private** | limit exp times frac quotient at plus infinity |
| `parse_signed_exp_frac_product` | **Partial** | 解析带符号 exp/分式积（快路径）; 退役: GIAC-limit-mrv-followup |
| `cancel_var_power_in_mul` | **Pipeline private** | cancel var power in mul |
| `extract_mul_sign` | **Pipeline private** | extract mul sign |
| `flatten_top_mul` | **Pipeline private** | flatten top mul |
| `is_limit_var` | **Pipeline private** | is limit var |
| `dominant_ln_exponent_at_plus_infinity` | **Partial** | 主导 ln 指数（快路径）; 退役: GIAC-limit-mrv-followup |
| `limit_const_rational_at_plus_infinity` | **Pipeline private** | limit const rational at plus infinity |
| `try_add_rational_to_frac` | **Pipeline private** | try add rational to frac |
| `split_mul_leading_coeff` | **Pipeline private** | split mul leading coeff |
| `as_frac_form` | **Pipeline private** | as frac form |
| `linear_var_term` | **Pipeline private** | linear var term |
| `float_to_expr` | **Pipeline private** | float to expr |
| `limit_add_at_plus_infinity` | **Pipeline private** | limit add at plus infinity |
| `limit_term_at_plus_infinity` | **Pipeline private** | limit term at plus infinity |
| `limit_at_plus_infinity_fallback` | **Partial** | MRV 不可用时的 `+∞` fallback; 退役: GIAC-limit-mrv-followup |
| `limit_at_zero_fallback` | **Pipeline** | `0` 极限 fallback（倒数换元 + 级数） |
| `limit_at_zero_rational_lead` | **Pipeline private** | limit at zero rational lead |
| `collapse_unit_powers` | **Pipeline private** | collapse unit powers |
| `limit_at_zero_rational_lead_escalating` | **Pipeline private** | limit at zero rational lead escalating |
| `limit_at_zero_from_series` | **Pipeline private** | limit at zero from series |
| `limit_at_zero_from_series_escalating` | **Pipeline private** | limit at zero from series escalating |
| `series_ordre_escalation` | **Pipeline private** | series ordre escalation |
| `asymptotic_series_at_infinity` | **Pipeline** | `+∞` 渐近级数 |
| `limit_rational_leading_at_infinity` | **Pipeline private** | limit rational leading at infinity |
| `is_sqrt` | **Pipeline private** | is sqrt |
| `sqrt_arg` | **Pipeline private** | sqrt arg |
| `limit_rational_over_sqrt_quotient_at_infinity` | **Pipeline private** | limit rational over sqrt quotient at infinity |
| `limit_sqrt_sum_quotient_at_infinity` | **Pipeline private** | limit sqrt sum quotient at infinity |
| `sqrt_sum_leading_linear_coeff` | **Pipeline private** | sqrt sum leading linear coeff |
| `is_monic_quadratic_leading` | **Pipeline private** | is monic quadratic leading |
| `is_indeterminate` | **Pipeline private** | is indeterminate |
| `peel_shared_u_inv_in_frac` | **Pipeline** | 分式中剥离共享 `u^-1` |
| `cancel_u_factors` | **Pipeline private** | cancel u factors |
| `flatten_mul` | **Pipeline private** | flatten mul |
| `normalize_limit_result` | **Pipeline private** | normalize limit result |
| `expr_to_ratio` | **Pipeline private** | expr to ratio |
| `u_exponent_bound` | **Pipeline private** | u exponent bound |
| `limit_from_fractional_u_valuation` | **Pipeline private** | limit from fractional u valuation |
| `is_inv_var` | **Pipeline private** | is inv var |
| `limit_var_over_x_pow_ln` | **Pipeline private** | limit var over x pow ln |
| `limit_poly_over_sqrt_at_infinity` | **Pipeline private** | limit poly over sqrt at infinity |
| `sqrt_term_in_expr` | **Pipeline private** | sqrt term in expr |
| `limit_exp_sum_nth_root` | **Pipeline private** | limit exp sum nth root |
| `rewrite_u_inv_sums` | **Pipeline private** | rewrite u inv sums |
| `peel_u_inv_factor` | **Pipeline private** | peel u inv factor |
| `is_u_inv` | **Pipeline private** | is u inv |
| `is_u_var` | **Pipeline private** | is u var |
| `is_one_plus_u_inv_sq` | **Pipeline private** | is one plus u inv sq |
| `u_negative_power_degree` | **Pipeline private** | u negative power degree |
| `simplify_reciprocal_sqrt` | **Pipeline private** | simplify reciprocal sqrt |
| `is_usable_limit` | **Pipeline private** | is usable limit |
| `contains_asym_var` | **Pipeline private** | contains asym var |
| `contains_zero_negative_power` | **Pipeline private** | contains zero negative power |
| `reciprocal_subst` | **Pipeline private** | reciprocal subst |
| `subst_map` | **Pipeline private** | subst map |
| `limit_from_rational_laurent` | **Pipeline private** | limit from rational laurent |
| `valuation_at_zero` | **Pipeline private** | valuation at zero |
| `limit_from_laurent_exponent` | **Pipeline private** | limit from laurent exponent |
| `sign_infinity` | **Pipeline private** | sign infinity |
| `limit_from_scaled_finite` | **Pipeline private** | limit from scaled finite |
| `sign_infinity_from_value` | **Pipeline private** | sign infinity from value |
| `laurent_terms_at_zero` | **Pipeline private** | laurent terms at zero |
| `laurent_terms_rational` | **Pipeline private** | laurent terms rational |
| `taylor_terms_as_laurent` | **Pipeline private** | taylor terms as laurent |
| `eval_at` | **Pipeline private** | eval at |
| `is_zero` | **Pipeline private** | is zero |
| `is_plus_infinity` | **Pipeline private** | is plus infinity |
| `is_minus_infinity` | **Pipeline private** | is minus infinity |
| `asymptotic_series_exp_at_infinity` | **Pipeline private** | `asymptotic_series_exp_at_infinity` |
| `asymptotic_maxima_sqrt_conjugate_shape` | **Pipeline private** | `asymptotic_maxima_sqrt_conjugate_shape` |
| `asymptotic_ck_int_56` | **Pipeline private** | `asymptotic_ck_int_56` |
| `asymptotic_ck_int_61` | **Pipeline private** | `asymptotic_ck_int_61` |
| `asymptotic_rational_infinity` | **Pipeline private** | `asymptotic_rational_infinity` |
| `assert_fallback` | **Pipeline private** | assert fallback |
| `rational_one` | **Pipeline private** | `rational_one` |
| `sqrt_conjugate_x` | **Pipeline private** | `sqrt_conjugate_x` |
| `ck_int_58` | **Pipeline private** | `ck_int_58` |
| `atan_over_x_plus_one` | **Pipeline private** | `atan_over_x_plus_one` |
| `one_plus_one_over_x_sqrt` | **Pipeline private** | `one_plus_one_over_x_sqrt` |

### `limit_engine/bounds.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `expr_nodes` | **Pipeline** | 表达式节点计数 |
| `expr_depth` | **Pipeline** | 表达式深度 |
| `too_heavy_for_expand` | **Pipeline** | expand 过重门禁 |
| `mrv_series_eligible` | **Pipeline** | MRV 级数资格检查 |
| `mrv_limit_eligible` | **Pipeline** | MRV 极限资格检查 |
| `mrv_rewrite_bounded` | **Pipeline** | MRV 换元后节点上限检查 |
| `expr_contains_exp` | **Pipeline** | 子树含 exp |
| `expr_contains_nested_exp` | **Pipeline** | 子树含嵌套 exp |
| `walk` | **Pipeline private** | walk |
| `expr_nodes_and_depth` | **Pipeline private** | `expr_nodes_and_depth` |
| `mrv_eligibility_flags` | **Pipeline private** | `mrv_eligibility_flags` |
| `nested_exp_detected` | **Pipeline private** | `nested_exp_detected` |
| `is_var_helper` | **Pipeline private** | `is_var_helper` |
| `expr_contains_exp_paths` | **Pipeline private** | `expr_contains_exp_paths` |
| `expr_nodes_relation_branch` | **Pipeline private** | `expr_nodes_relation_branch` |

### `limit_engine/ck_int_gruntz_fixture.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `inner` | **Pipeline private** | `inner` |
| `exp_inner` | **Pipeline private** | `exp_inner` |
| `ratio` | **Pipeline private** | `ratio` |
| `ck_int_60` | **Pipeline private** | `ck_int_60` |
| `ck_int_61` | **Pipeline private** | `ck_int_61` |

### `limit_engine/exp_diff.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `canonical_exp_diff` | **Stable** | x-layer exp-difference canonical form;唯一漂移收敛入口 |
| `exp_scale_times_exp_minus_one` | **Stable** | 规范构造器 `scale·(exp(ε)-1)`（Add 形） |
| `emit_exp_difference` | **Pipeline private** | emit exp difference |
| `detect_exp_difference` | **Pipeline private** | detect exp difference |
| `flatten_mul_factors` | **Pipeline private** | flatten mul factors |
| `detect_exp_difference_mul` | **Pipeline private** | detect exp difference mul |
| `flatten_mul_from_vec` | **Pipeline private** | flatten mul from vec |
| `sum_exp_logs` | **Pipeline private** | sum exp logs |
| `detect_exp_difference_add` | **Pipeline private** | detect exp difference add |
| `difference_add_exprs` | **Pipeline private** | difference add exprs |
| `canonical_add_term` | **Pipeline private** | canonical add term |
| `signed_add_terms` | **Pipeline private** | signed add terms |
| `peel_unit_negative_factor` | **Pipeline private** | peel unit negative factor |
| `rebuild_mul` | **Pipeline private** | rebuild mul |
| `flatten_mul_expr` | **Pipeline private** | flatten mul expr |
| `flatten_add_mul_terms` | **Pipeline private** | flatten add mul terms |
| `distribute_linear_mul_in_add` | **Pipeline private** | distribute linear mul in add |
| `distribute_linear_mul_term` | **Pipeline private** | distribute linear mul term |
| `simplify_balanced_frac` | **Pipeline private** | simplify balanced frac |
| `rebuild_signed_add` | **Pipeline private** | rebuild signed add |
| `exp_signed_pair` | **Pipeline private** | exp signed pair |
| `rewrite_exp_minus_scale_inv` | **Pipeline** | 特定 `exp(-s/x)` 倒数改写 |
| `simplify_add_sum` | **Pipeline** | 合并同类 Add 项 |
| `simplify_exp_argument_adds` | **Pipeline** | 合并 exp 参数中的加法 |
| `balance_exp_arguments_frac_var` | **Pipeline** | 平衡分式指数与 `-var` |
| `balance_frac_minus_var` | **Pipeline private** | balance frac minus var |
| `try_balance_frac_minus_var` | **Pipeline private** | try balance frac minus var |
| `linear_coeff_of_var` | **Pipeline private** | linear coeff of var |
| `classify_exp_at_plus_infinity` | **Partial** | +∞ 无符号 exp 增长分类（快路径）; 退役: GIAC-limit-mrv-followup |
| `unwrap_signed_frac` | **Stable** | 提取 `(sign, inner)` 供 +∞ 分类 |
| `classify_signed_exp_at_plus_infinity` | **Partial** | +∞ 带符号分式积 exp 增长分类（快路径）; 退役: GIAC-limit-mrv-followup |
| `exp_argument_tends_to_negative_infinity` | **Pipeline private** | exp argument tends to negative infinity |
| `is_positive_var_power_ge2` | **Pipeline private** | is positive var power ge2 |
| `is_neg_exp_of_positive_growth` | **Pipeline private** | is neg exp of positive growth |
| `exp_arg_grows_at_plus_infinity` | **Pipeline private** | exp arg grows at plus infinity |
| `algebraize_exp_vanishing_products` | **Partial** | `exp(f)·L` 代数化; 退役: GIAC-limit-exp-difference-unification |
| `algebraize_exp_mul_factors` | **Pipeline private** | algebraize exp mul factors |
| `algebraize_exp_log_pair` | **Pipeline private** | algebraize exp log pair |
| `flatten_mul_factors_slice` | **Pipeline private** | flatten mul factors slice |
| `first_order_exp_vanishing_epsilon` | **Partial** | Gruntz ε 展开（破坏 MRV peel 所需 `exp(·)-1` 形）; 退役: GIAC-limit-mrv-followup |
| `match_exp_times_exp_minus_one` | **Stable** | 识别 `exp(scale)·(exp(ε)-1)` 分解 |
| `unwrap_exp_arg` | **Stable** | 提取 `exp(arg)` / `-exp(arg)` 的参数 |
| `exp_func_arg` | **Pipeline private** | exp func arg |
| `exp_minus_one_epsilon` | **Stable** | 若因子为 `exp(ε)-1` 则返回 ε |
| `is_exp_minus_one_factor` | **Stable** | 谓词：是否为 `exp(ε)-1` 因子（只认 Add 形） |
| `exp_minus_one_inner_arg` | **Pipeline private** | exp minus one inner arg |
| `remove_add_term` | **Pipeline private** | remove add term |
| `vanishes_at_plus_infinity` | **Stable** | 谓词：`e→0` 当 `var→+∞` |
| `epsilon_vanishes_at_plus_infinity` | **Stable** | 谓词：ε 级小量（Gruntz 预处理） |
| `balance_epsilon_expr` | **Pipeline private** | balance epsilon expr |
| `exp_inner_vanishes_at_plus_infinity` | **Pipeline private** | exp inner vanishes at plus infinity |
| `exp_of_neg_exp_growth` | **Pipeline private** | exp of neg exp growth |
| `arg_grows` | **Pipeline private** | arg grows |
| `is_negated_var_power_mul` | **Pipeline private** | is negated var power mul |
| `is_neg_var_exp` | **Pipeline private** | is neg var exp |
| `exp_inner_arg` | **Pipeline private** | `exp_inner_arg` |
| `expr_mentions_symbol` | **Pipeline private** | `expr_mentions_symbol` |
| `is_exp_var_times_x_inv` | **Pipeline private** | `is_exp_var_times_x_inv` |
| `tree_contains_exp_of_sym` | **Pipeline private** | `tree_contains_exp_of_sym` |
| `fold_mul_shifted_difference_gruntz` | **Pipeline private** | `fold_mul_shifted_difference_gruntz` |
| `rewrite_exp_minus_w_inv_matches_remove_lnexp` | **Pipeline private** | `rewrite_exp_minus_w_inv_matches_remove_lnexp` |
| `fold_add_exp_difference` | **Pipeline private** | `fold_add_exp_difference` |
| `balance_frac_minus_var_rewrites_inner_minus_x` | **Pipeline private** | `balance_frac_minus_var_rewrites_inner_minus_x` |
| `first_order_vanishing_epsilon_rewrite` | **Pipeline private** | `first_order_vanishing_epsilon_rewrite` |
| `limit_preprocessed_gruntz_minus_one` | **Pipeline private** | `limit_preprocessed_gruntz_minus_one` |
| `ratio_preprocess_mrv_cancels_opposing_exp_mul` | **Pipeline private** | `ratio_preprocess_mrv_cancels_opposing_exp_mul` |
| `limit_exp_over_exp_ck_int_ratio` | **Pipeline private** | `limit_exp_over_exp_ck_int_ratio` |
| `fold_ck_int_61_shape` | **Pipeline private** | `fold_ck_int_61_shape` |
| `fold_nested_gruntz_add_difference` | **Pipeline private** | `fold_nested_gruntz_add_difference` |
| `fold_nested_gruntz_full_expr` | **Pipeline private** | `fold_nested_gruntz_full_expr` |
| `parse_neg_x_squared_in_exp` | **Pipeline private** | `parse_neg_x_squared_in_exp` |
| `limit_nested_gruntz_preprocessed` | **Pipeline private** | `limit_nested_gruntz_preprocessed` |
| `exp_scale_times_exp_minus_one_shape` | **Pipeline private** | `exp_scale_times_exp_minus_one_shape` |

### `limit_engine/mod.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `expr_has_nested_exp` | **Pipeline** | 嵌套 exp 检测（委托 bounds） |
| `limit_finite_algebraic` | **Pipeline** | 有限点代数极限 |
| `try_as_quotient` | **Pipeline private** | try as quotient |
| `is_negative_var_power` | **Pipeline private** | is negative var power |
| `positive_var_power` | **Pipeline private** | positive var power |
| `limit_quotient_plus_infinity` | **Pipeline private** | limit quotient plus infinity |
| `limit_plus_infinity_quick` | **Pipeline private** | limit plus infinity quick |
| `try_as_quotient_add_shared_power` | **Pipeline private** | try as quotient add shared power |
| `limit_quotient_finite` | **Pipeline private** | limit quotient finite |
| `limit_plus_infinity_algebraic` | **Pipeline** | `+∞` 代数极限 |
| `limit_minus_infinity_algebraic` | **Pipeline** | `-∞` 代数极限 |
| `limit_rational_finite` | **Pipeline private** | limit rational finite |
| `cancel_rational_pole` | **Pipeline private** | cancel rational pole |
| `limit_rational_infinity` | **Pipeline private** | limit rational infinity |
| `leading_ratio` | **Pipeline private** | leading ratio |
| `sign_infinity` | **Pipeline private** | sign infinity |
| `pole_infinity` | **Pipeline private** | pole infinity |
| `limit_via_reciprocal` | **Pipeline private** | limit via reciprocal |
| `subst_eval` | **Pipeline private** | subst eval |
| `expr_contains_exp` | **Pipeline private** | expr contains exp |
| `subst_map` | **Pipeline private** | subst map |
| `is_zero` | **Pipeline private** | is zero |
| `is_one` | **Pipeline private** | is one |
| `is_plus_infinity` | **Pipeline private** | is plus infinity |
| `is_minus_infinity` | **Pipeline private** | is minus infinity |
| `is_infinity` | **Pipeline private** | is infinity |
| `contains_zero_negative_power` | **Pipeline private** | contains zero negative power |
| `is_indeterminate` | **Pipeline private** | is indeterminate |
| `limit_line` | **Pipeline private** | limit line |
| `classify_test_point` | **Pipeline private** | classify test point |
| `engine_ck_int_56` | **Pipeline private** | `engine_ck_int_56` |
| `engine_ck_int_58` | **Pipeline private** | `engine_ck_int_58` |
| `engine_ck_int_59_via_finite` | **Pipeline private** | `engine_ck_int_59_via_finite` |
| `engine_ck_int_gruntz_ratio` | **Pipeline private** | `engine_ck_int_gruntz_ratio` |
| `engine_ck_int_60` | **Pipeline private** | `engine_ck_int_60` |
| `engine_ck_int_61` | **Pipeline private** | `engine_ck_int_61` |
| `engine_limit_minus_infinity_exp` | **Pipeline private** | `engine_limit_minus_infinity_exp` |
| `engine_limit_finite_rational_pole_cancel` | **Pipeline private** | `engine_limit_finite_rational_pole_cancel` |
| `engine_limit_indeterminate_at_one` | **Pipeline private** | `engine_limit_indeterminate_at_one` |
| `engine_limit_minus_infinity_poly` | **Pipeline private** | `engine_limit_minus_infinity_poly` |

### `limit_engine/mrv.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `is_empty` | **Pipeline private** | MRV 集合是否为空 |
| `merge` | **Pipeline private** | 合并两个 MRV 集合 |
| `mrv_at_plus_infinity` | **Stable** | 构建 `+∞` MRV 集合 |
| `collect_mrv` | **Pipeline private** | collect mrv |
| `ratio_tends_to_zero_at_plus_infinity` | **Stable** | 谓词：`n/d→0` at `+∞` |
| `vanishes_faster_than_at_plus_infinity` | **Stable** | 谓词：`a/b→0` at `+∞` |
| `mrv_compare` | **Stable** | `+∞` 增长比较（上游 `mrv_compare` 子集） |
| `growth_rank_detailed` | **Pipeline private** | growth rank detailed |
| `limit_const_pow_quotient_at_plus_infinity` | **Stable** | 常数底幂商在 `+∞` 的极限 |
| `asymptotic_linear_coeff_at_plus_infinity` | **Pipeline private** | asymptotic linear coeff at plus infinity |
| `leading_coeff_in_var` | **Pipeline private** | leading coeff in var |
| `limit_rational_const_at_plus_infinity` | **Stable** | 有理/常数 lead 在 `+∞` |
| `add_degree_in_var` | **Pipeline private** | add degree in var |
| `merge_mrv_pair` | **Pipeline private** | merge mrv pair |
| `growth_rank` | **Pipeline private** | growth rank |
| `linear_coeff_in_var` | **Stable** | 提取 `var` 的线性系数（含 Frac 指数） |
| `is_negative_const_expr` | **Stable** | 符号负常数谓词 |
| `try_const_f64` | **Stable** | 尝试将表达式读为 `f64` 常数 |
| `choose_mrv_w` | **Stable** | 从 MRV 集选取换元 `w` |
| `is_neg_linear` | **Pipeline private** | is neg linear |
| `is_neg_scaled_linear` | **Pipeline private** | is neg scaled linear |
| `linear_coeff` | **Pipeline private** | linear coeff |
| `expr_size` | **Pipeline private** | expr size |
| `mrv_nested_exp_inner` | **Pipeline private** | `mrv_nested_exp_inner` |
| `choose_mrv_seven_over_eight_pow_n` | **Pipeline private** | `choose_mrv_seven_over_eight_pow_n` |
| `mrv_exp_neg_x` | **Pipeline private** | `mrv_exp_neg_x` |
| `mrv_exp_neg_x_substitutes_ln_w` | **Pipeline private** | `mrv_exp_neg_x_substitutes_ln_w` |

### `limit_engine/mrv_lead_term.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `mrv_lead_exp_vanishing_plus_infinity` | **Pipeline private** | mrv lead exp vanishing plus infinity |
| `mrv_lead_term_plus_infinity` | **Stable** | MRV lead 管线入口（`+∞`） |
| `limit_from_mrv_lead_term` | **Stable** | 由 `MrvLeadTerm` 求极限 |
| `limit_unidirectional_plus_infinity` | **Stable** | 单向 `+∞` 极限（递归 lead） |
| `subst_symbol` | **Pipeline private** | subst symbol |
| `upscale_while_var_in_mrv` | **Pipeline private** | upscale while var in mrv |
| `omega_neg_linear_coeff` | **Pipeline private** | omega neg linear coeff |
| `neg_ln_w_from_omega` | **Pipeline private** | neg ln w from omega |
| `linear_inner_coeff_only` | **Pipeline private** | linear inner coeff only |
| `rewrite_in_mrv_w` | **Pipeline private** | rewrite in mrv w |
| `expr_eq` | **Pipeline private** | expr eq |
| `is_exp_neg_var` | **Pipeline private** | is exp neg var |
| `is_exp_pos_var` | **Pipeline private** | is exp pos var |
| `is_neg_var` | **Pipeline private** | is neg var |
| `rewrite_ln_w` | **Pipeline private** | rewrite ln w |
| `g_from_omega` | **Pipeline private** | g from omega |
| `omega_tends_to_zero_at_plus_inf` | **Pipeline private** | omega tends to zero at plus inf |
| `is_neg_var_linear` | **Pipeline private** | is neg var linear |
| `subst_w_inv` | **Pipeline private** | subst w inv |
| `is_ln_w` | **Pipeline private** | is ln w |
| `subst_ln_w_expr` | **Pipeline private** | subst ln w expr |
| `collect_ln_exprs` | **Pipeline private** | collect ln exprs |
| `rewrite_ln_exp_in_f` | **Pipeline private** | rewrite ln exp in f |
| `subst_expr_once` | **Pipeline private** | subst expr once |
| `peel_neg_ln_w_inv` | **Pipeline private** | peel neg ln w inv |
| `combine_lead_with_ln_inv` | **Pipeline private** | combine lead with ln inv |
| `is_series_coeff_undef` | **Pipeline private** | is series coeff undef |
| `normalize_series_coeff` | **Pipeline private** | normalize series coeff |
| `pnormal_series` | **Pipeline private** | pnormal series |
| `series_lead_at_zero` | **Stable** | `w=0` 级数 lead |
| `mrv_series_lead_loop` | **Pipeline private** | mrv series lead loop |
| `lead_from_peeled_core` | **Pipeline private** | `lead_from_peeled_core` |
| `mrv_series_lead_loop_inner` | **Pipeline private** | mrv series lead loop inner |
| `lead_coeff_ready` | **Pipeline private** | lead coeff ready |
| `depends_on_w` | **Pipeline private** | depends on w |
| `sign_infinity` | **Pipeline private** | sign infinity |
| `mrv_series_expansion_order_cap` | **Pipeline private** | `mrv_series_expansion_order_cap` |
| `limit_seven_pow_n_over_eight` | **Pipeline private** | `limit_seven_pow_n_over_eight` |
| `rewrite_exp_x_to_w_inv` | **Pipeline private** | `rewrite_exp_x_to_w_inv` |
| `rewrite_x_powers_in_mrv_w` | **Pipeline private** | `rewrite_x_powers_in_mrv_w` |
| `rewrite_ck61_normed_in_mrv_w` | **Pipeline private** | `rewrite_ck61_normed_in_mrv_w` |
| `ratio_mrv_lead_exp_vanishing_after_preprocess` | **Pipeline private** | `ratio_mrv_lead_exp_vanishing_after_preprocess` |
| `mrv_lead_ck_int_60_ratio` | **Pipeline private** | `mrv_lead_ck_int_60_ratio` |
| `mrv_lead_ck_int_61` | **Pipeline private** | `mrv_lead_ck_int_61` |

### `limit_engine/mrv_series_lead.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `try_as_quotient` | **Pipeline** | 提取 `(numerator, denominator)` |
| `normalize_expr_quotients` | **Pipeline** | 商式 Laurent 形态统一 |
| `unify_top_quotient` | **Pipeline** | 顶层 `Frac` → `Mul·den^-1` |
| `canonicalize_limit_entry` | **Pipeline** | 极限入口规范形 |
| `normalize_inverse_sums` | **Pipeline** | `a*(b+c)^-1` → `a/(b+c)` |
| `normalize_inverse_products` | **Pipeline private** | normalize inverse products |
| `normalize_inverse_with` | **Pipeline private** | normalize inverse with |

### `limit_engine/mrv_w.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `is_mrv_w_var` | **Stable** | 谓词：是否为 `_mrv_w` 符号 |
| `mrv_w_expr` | **Stable** | 规范原子 `Symbol(_mrv_w)` |
| `mrv_ln_w_expr` | **Stable** | 规范原子 `Ln(_mrv_w)` |
| `neg_ln_w_expr` | **Stable** | 规范原子 `-ln(w)` |
| `neg_ln_w_inv_expr` | **Stable** | 规范原子 `(-ln(w))^-1` |
| `is_expr_one` | **Stable** | 谓词：表达式为常数 1 |
| `is_neg_w_inv` | **Stable** | 谓词：Laurent 因子 `w^-1`（≠ `(-ln w)^-1`） |
| `is_negative_unit_exp` | **Pipeline private** | is negative unit exp |
| `is_negative_one_like` | **Pipeline private** | is negative one like |
| `expr_contains_w_var` | **Stable** | 谓词：子树含 `w` |
| `expr_contains_ln_w` | **Stable** | 谓词：子树含 `ln(w)` |
| `is_neg_ln_first_power` | **Stable** | 谓词：规范后等于 `(-ln(w))^1` |
| `canonical_mrv_coeff` | **Stable** | MRV 系数规范化；唯一漂移收敛入口 |
| `pending_for_series` | **Stable** | 系数仍含未处理 `ln(w)` / `(-ln(w))^k` |
| `has_neg_ln_inv` | **Stable** | 含 `(-ln(w))^-1`（需 peel + padd） |
| `decompose_mrv_coeff` | **Stable** | 语义分解 `(-ln(w))^k · ln(w) 加法部分 · rest` |
| `decompose_ln_w_coeff` | **Stable** | 仅 `k·ln(w)` 加法基分解 |
| `decompose_neg_ln_power` | **Pipeline private** | decompose neg ln power |
| `is_neg_ln_atom` | **Pipeline private** | is neg ln atom |
| `is_ln_w` | **Pipeline private** | is ln w |
| `ln_w_mul_coeff` | **Pipeline private** | ln w mul coeff |
| `rest_expr` | **Pipeline private** | rest expr |
| `drift_fold_ln_atoms` | **Temporary** | `drift_fold_ln_atoms` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A |
| `drift_is_neg_ln_shape` | **Temporary** | `drift_is_neg_ln_shape` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A |
| `drift_is_neg_ln_inv_shape` | **Temporary** | `drift_is_neg_ln_inv_shape` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A |
| `drift_is_neg_ln_expanded` | **Temporary** | `drift_is_neg_ln_expanded` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A |
| `drift_is_neg_ln_inv_expanded` | **Temporary** | `drift_is_neg_ln_inv_expanded` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A |
| `drift_is_negative_one_like` | **Temporary** | `drift_is_negative_one_like` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A |
| `drift_is_ln_w_symbol` | **Temporary** | `drift_is_ln_w_symbol` 漂移形识别; 退役: GIAC-limit-mrv-followup Phase 3A |
| `decompose_mrv_coeff_neg_ln_inv_cancels_in_ratio` | **Pipeline private** | `decompose_mrv_coeff_neg_ln_inv_cancels_in_ratio` |
| `canonical_mrv_coeff_preserves_neg_ln_to_first_power` | **Pipeline private** | `canonical_mrv_coeff_preserves_neg_ln_to_first_power` |
| `canonical_mrv_coeff_folds_drift_neg_ln` | **Pipeline private** | `canonical_mrv_coeff_folds_drift_neg_ln` |
| `is_neg_ln_first_power_atom_and_drift` | **Pipeline private** | `is_neg_ln_first_power_atom_and_drift` |
| `mrv_coeff_parts_pending_flags` | **Pipeline private** | `mrv_coeff_parts_pending_flags` |

### `limit_engine/preprocess.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `limit_preprocess_mrv` | **Pipeline** | MRV 路径预处理（无 Gruntz ε 改写） |
| `limit_preprocess_struct` | **Pipeline** | 结构预处理（含 ε 展开；勿用于 MRV peel 前） |
| `limit_preprocess_plus_infinity` | **Pipeline** | `+∞` 预处理编排 |
| `series_preprocess` | **Pipeline** | 级数路径预处理 |
| `surd2pow` | **Pipeline** | `sqrt` → 有理指数 |
| `sqrt_operand` | **Pipeline private** | sqrt operand |
| `negated_inner` | **Pipeline private** | negated inner |
| `try_sqrt_difference_frac` | **Pipeline private** | try sqrt difference frac |
| `try_sqrt_minus_x_frac` | **Pipeline private** | try sqrt minus x frac |
| `combine_mul_with_frac` | **Pipeline private** | combine mul with frac |
| `normalize_sqrt_conjugates` | **Pipeline** | 共轭有理化 `sqrt` 差分 |
| `fold_exp_zero_linear` | **Pipeline** | `exp(0·x+…)` 线性折叠 |
| `factor_exp_shifted_difference` | **Stable** | 转发 `canonical_exp_diff` |
| `merge_exp_quotients` | **Pipeline** | 合并 exp 商式 |
| `merge_exp_quotient_frac` | **Pipeline private** | `merge_exp_quotient_frac` |
| `preprocess_seven_pow_n_over_eight` | **Pipeline private** | `preprocess_seven_pow_n_over_eight` |
| `preprocess_surd2pow_sqrt` | **Pipeline private** | `preprocess_surd2pow_sqrt` |
| `preprocess_sqrt_conjugate_minus_var` | **Pipeline private** | `preprocess_sqrt_conjugate_minus_var` |
| `preprocess_nested_gruntz_exp_diff_factor` | **Pipeline private** | `preprocess_nested_gruntz_exp_diff_factor` |
| `preprocess_gruntz_exp_diff_factor` | **Pipeline private** | `preprocess_gruntz_exp_diff_factor` |

### `limit_engine/remove_lnexp.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `remove_lnexp` | **Stable** | 上游 `remove_lnexp`：`ln_expand` + `exp_series` |
| `try_collapse_w_inv_exp_plus_ln` | **Pipeline private** | try collapse w inv exp plus ln |
| `try_collapse_w_inv_exp_shift` | **Pipeline private** | try collapse w inv exp shift |
| `is_mrv_w_var_symbol` | **Pipeline private** | is mrv w var symbol |
| `try_rewrite_w_inv_times_exp_minus_one` | **Pipeline private** | try rewrite w inv times exp minus one |
| `try_rewrite_exp_minus_w_inv` | **Pipeline private** | try rewrite exp minus w inv |
| `lead_after_exp_ln_cancel` | **Pipeline private** | lead after exp ln cancel |
| `second_term_inner_plus_ln_expr` | **Pipeline private** | second term inner plus ln expr |
| `divide_lead_coeffs` | **Stable** | padd 商；同幂 `ln(w)`/`(-ln(w))^k` 相消 |
| `expr_contains_exp_or_ln` | **Stable** | 粗谓词：子树含 exp 或 ln |
| `fold_children` | **Pipeline private** | fold children |
| `ln_expand0` | **Pipeline private** | ln expand0 |
| `exp_series_expand` | **Pipeline private** | exp series expand |
| `exp_series_ln_w` | **Pipeline private** | exp series ln w |
| `unique_ln_subexpr` | **Pipeline private** | unique ln subexpr |
| `collect_ln` | **Pipeline private** | collect ln |
| `linear_decompose_wrt` | **Pipeline private** | linear decompose wrt |
| `is_integer_like` | **Pipeline private** | is integer like |
| `remove_lnexp_exp_ln_w` | **Pipeline private** | `remove_lnexp_exp_ln_w` |
| `remove_lnexp_exp_difference` | **Pipeline private** | `remove_lnexp_exp_difference` |
| `divide_lead_coeffs_neg_ln_w_inv_cancels` | **Pipeline private** | `divide_lead_coeffs_neg_ln_w_inv_cancels` |
| `divide_lead_coeffs_ln_w_cancels` | **Pipeline private** | `divide_lead_coeffs_ln_w_cancels` |
| `expr_contains_exp_or_ln_detects` | **Pipeline private** | `expr_contains_exp_or_ln_detects` |

### `limit_engine/simplify_util.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `simplify_limit_expr` | **Pipeline** | `ratnormal` 后 `normal`；失败保留输入 |
| `is_negative_const_expr` | **Stable** | 符号负常数谓词 |
| `int_pow_growth_sub_rank` | **Pipeline** | 整数底 `a^var` 的精确 sub-rank |

### `limit_engine/sparse_series.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `constant` | **Pipeline** | constant |
| `is_zero` | **Pipeline** | is zero |
| `lead` | **Pipeline** | lead |
| `term_count` | **Pipeline** | term count |
| `iter_terms` | **Pipeline** | iter terms |
| `to_expr` | **Pipeline** | to expr |
| `add` | **Pipeline** | add |
| `mul_with_cap` | **Pipeline** | mul with cap |
| `mul` | **Pipeline** | mul |
| `truncate_to_order` | **Pipeline private** | truncate to order |
| `truncate_to_order_cap` | **Pipeline private** | truncate to order cap |
| `truncate_display_order` | **Pipeline** | truncate display order |
| `map_coeffs` | **Pipeline** | map coeffs |
| `from_rational_laurent` | **Pipeline** | from rational laurent |
| `series_at_center` | **Pipeline** | 通用中心级数展开 |
| `series_sparse_to_expr` | **Pipeline private** | series sparse to expr |
| `series_at_zero` | **Stable** | `w=0` 稀疏级数（MRV 系数域） |
| `series_at_zero_order` | **Stable** | 指定阶数的 `w=0` 级数 |
| `series_spdiv_one` | **Stable** | 级数升阶除法 `spdiv(·,1)` |
| `series_at_zero_depth` | **Pipeline private** | series at zero depth |
| `series_sin` | **Pipeline private** | series sin |
| `series_cos` | **Pipeline private** | series cos |
| `series_sqrt` | **Pipeline private** | series sqrt |
| `series_atan` | **Pipeline private** | series atan |
| `series_atan_of_inv` | **Pipeline private** | series atan of inv |
| `is_inv_series_var` | **Pipeline private** | is inv series var |
| `series_exp_mrv` | **Pipeline private** | series exp mrv |
| `shift_series_exponents` | **Pipeline private** | shift series exponents |
| `series_exp_mrv_positive` | **Pipeline private** | series exp mrv positive |
| `series_ln_mrv` | **Pipeline private** | series ln mrv |
| `series_exp` | **Pipeline private** | series exp |
| `series_pow_int` | **Pipeline private** | series pow int |
| `series_div` | **Pipeline private** | series div |
| `series_inv` | **Pipeline private** | series inv |
| `merge_term` | **Pipeline private** | merge term |
| `simplify_series_coeff` | **Pipeline private** | simplify series coeff |
| `normalize_map` | **Pipeline private** | normalize map |
| `valuation_at_zero` | **Pipeline private** | valuation at zero |
| `bigint_to_i64` | **Pipeline private** | bigint to i64 |
| `factorial` | **Pipeline private** | factorial |
| `is_series_var` | **Pipeline private** | is series var |
| `arg_has_symbolic_ln_w` | **Pipeline private** | arg has symbolic ln w |
| `is_mrv_symbolic_pow` | **Pipeline private** | is mrv symbolic pow |
| `series_mrv_w_exp_minus_w_inv_cancel_reveals_sublead` | **Pipeline private** | `series_mrv_w_exp_minus_w_inv_cancel_reveals_sublead` |
| `sparse_series_rational_at_zero` | **Pipeline private** | `sparse_series_rational_at_zero` |
| `sparse_series_sin_at_zero` | **Pipeline private** | `sparse_series_sin_at_zero` |
| `sparse_series_at_center_shift` | **Pipeline private** | `sparse_series_at_center_shift` |
| `series_at_zero_order_escalates_beyond_default_cap` | **Pipeline private** | `series_at_zero_order_escalates_beyond_default_cap` |
| `sparse_series_atan_inv_over_one_plus_u` | **Pipeline private** | `sparse_series_atan_inv_over_one_plus_u` |
| `sparse_series_exp_at_zero` | **Pipeline private** | `sparse_series_exp_at_zero` |
| `sparse_series_cos_at_zero` | **Pipeline private** | `sparse_series_cos_at_zero` |
| `sparse_series_one_over_one_plus_x` | **Pipeline private** | `sparse_series_one_over_one_plus_x` |

### `limit_engine/util.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `is_expr_zero` | **Pipeline private** | `is_expr_zero` |
| `is_half_exponent` | **Pipeline private** | `is_half_exponent` |

### `partfrac_integrate.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `integrate_one_over_quadratic` | **Stable** | ∫ 1/(ax²+bx+c) dx for constant-coefficient denominator (degree 1 or 2). |
| `integrate_rational_partfrac` | **Pipeline private** | integrate rational via partial fractions (ℚ, then K fallback). |
| `integrate_q_partfrac_terms` | **Pipeline private** | integrate ℚ partfrac terms. |
| `integrate_k_partfrac` | **Pipeline private** | integrate K partfrac terms (`partfrac_rational_terms_over_k`). |
| `integrate_algext_linear_term` | **Pipeline private** | ∫ c/(a·x+b) dx with c, a, b ∈ K (AlgExtC). |
| `integrate_const_over_rational` | **Stable** | ∫ num/den dx for rational expressions (Hermite, Rothstein–Trager, partfrac). |
| `den_perfect_power_expr` | **Pipeline private** | detect `base^exp` denominator with `exp >= 2`. |
| `hermite_factor_sign` | **Pipeline private** | sign correction for Hermite quadratic factors. |
| `integrate_with_hermite` | **Pipeline private** | integrate via Hermite reduction on repeated quadratics. |
| `integrate_hermite_term` | **Pipeline private** | integrate one Hermite reduction term. |
| `integrate_rational_term` | **Pipeline private** | integrate one partial-fraction term. |
| `integrate_poly_over_linear` | **Pipeline private** | ∫ P(x)/(cx+d) dx with linear denominator. |
| `integrate_const_over_power` | **Pipeline private** | ∫ c / g^n dx for linear `g` and n >= 2. |
| `integrate_polynomial_over_linear_power` | **Pipeline private** | ∫ P(x)/g^n dx for `g = c·(ax+b)^n`. |
| `ratio_pow` | **Pipeline private** | raise a rational to an integer power. |
| `integrate_over_quadratic` | **Pipeline private** | ∫ rational over irreducible or repeated quadratic. |
| `integrate_over_quadratic_real_roots` | **Pipeline private** | ∫ (B·t+C)/(a·t²+b·t+c) dt when the quadratic has real roots (disc > 0). |
| `ratio_sqrt` | **Pipeline private** | exact square root of a perfect-square rational. |
| `sqrt_ratio_expr` | **Pipeline private** | build `sqrt(r)` as `Expr` (exact or nested `sqrt`). |
| `integrate_poly_over_linear_term` | **Pipeline private** | `integrate_poly_over_linear_term` |
| `partfrac_integrate_x_over_repeated_linear` | **Pipeline private** | `partfrac_integrate_x_over_repeated_linear` |
| `integrate_one_over_x_fourth_plus_one` | **Pipeline private** | `integrate_one_over_x_fourth_plus_one` |
| `integrate_one_over_x_fourth_plus_one_squared` | **Pipeline private** | `integrate_one_over_x_fourth_plus_one_squared` |
| `integrate_x_over_x_squared_plus_one_squared_expr` | **Pipeline private** | `integrate_x_over_x_squared_plus_one_squared_expr` |
| `hermite_integrate_x_over_x_squared_plus_one_squared` | **Pipeline private** | `hermite_integrate_x_over_x_squared_plus_one_squared` |
| `integrate_one_over_x_fourth_minus_one_squared` | **Pipeline private** | `integrate_one_over_x_fourth_minus_one_squared` |
| `integrate_x_over_x_plus_one_times_x_fourth_minus_one` | **Pipeline private** | `integrate_x_over_x_plus_one_times_x_fourth_minus_one` |
| `integrate_one_over_x_fourth_plus_one_fourth_power` | **Pipeline private** | `integrate_one_over_x_fourth_plus_one_fourth_power` |
| `partfrac_integrate_half_angle_rational_in_t` | **Pipeline private** | `partfrac_integrate_half_angle_rational_in_t` |
| `integrate_one_over_quadratic_one_minus_x_squared` | **Pipeline private** | `integrate_one_over_quadratic_one_minus_x_squared` |
| `integrate_one_over_x_squared_minus_two_via_k_partfrac` | **Pipeline private** | `integrate_one_over_x_squared_minus_two_via_k_partfrac` |
| `integrate_x_over_x_squared_minus_two_via_k_partfrac` | **Pipeline private** | `integrate_x_over_x_squared_minus_two_via_k_partfrac` |
| `integrate_ck_int_05_reciprocal` | **Pipeline private** | `integrate_ck_int_05_reciprocal` |

### `plugin.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_integrate` | **Stable** | delegate to [`eval_integrate`]. |
| `eval_diff` | **Stable** | delegate to [`eval_diff`]. |
| `eval_limit` | **Stable** | delegate to [`eval_limit`]. |
| `eval_series` | **Stable** | delegate to [`eval_series`]. |
| `eval_risch` | **Partial** | delegate to [`eval_risch`] (narrow Risch subset). |
| `install_calculus` | **Stable** | install the default calculus plugin on `ctx`. |
| `xcas_default` | **Stable** | full CAS context: linear algebra, solving, simplification, and calculus. |
| `eval_integrate_via_plugin` | **Pipeline private** | `eval_integrate_via_plugin` |
| `eval_limit_via_plugin` | **Pipeline private** | `eval_limit_via_plugin` |
| `eval_series_via_plugin` | **Pipeline private** | `eval_series_via_plugin` |

### `risch/algebraic_rt.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `integrate_monic_x4_plus_one` | **Pipeline private** | `integrate_monic_x4_plus_one` |
| `try_algebraic_rt_log_part` | **Partial** | RT log part when `Res_t` has algebraic conjugate pairs on an even monic quartic. **退役：** unified RT resultant handler. |
| `is_monic_x4_plus_one` | **Partial** | monic quartic `x^4+1` shape predicate. |
| `is_monic_even_quartic` | **Partial** | monic even quartic shape predicate (odd coefficients zero). |
| `even_quartic_quadratic_factors` | **Pipeline private** | factor monic even quartic into two quadratics. |
| `linear_sqrt_coeff_expr` | **Pipeline private** | linear coefficient from `√(u²)` rational or surd. |
| `sqrt_rational_coeff_radicand` | **Pipeline private** | split rational into outer coeff and inner radicand. |
| `extract_sqrt_factor_bigint` | **Pipeline private** | extract perfect-square factor from integer. |
| `negate_linear_expr` | **Pipeline private** | negate linear expression coefficient. |
| `quadratic_expr` | **Pipeline private** | build `x² + b·x + c` expression. |
| `sort_pairs_and_factors` | **Pipeline private** | align conjugate pairs with quadratic factors. |
| `b_lin_sign_key` | **Pipeline private** | sort key from linear `b` coefficient sign. |
| `alpha_re_key` | **Pipeline private** | sort key from conjugate root real part. |
| `integrate_from_pairs` | **Pipeline private** | sum log/atan contributions from conjugate pairs. |
| `default_x4_plus_one_pairs` | **Pipeline private** | hard-coded conjugate pairs for `x^4+1`. |
| `conjugate_pair_log_contribution` | **Pipeline private** | `Re(α)·ln\|G\| + 2·Im(α)·atan((2x+b)/√(4c-b²))`. |
| `alpha_re_expr` | **Pipeline private** | real part of algebraic root as expression. |
| `alpha_im_expr` | **Pipeline private** | imaginary part of algebraic root as expression. |
| `sqrt_ratio_expr` | **Pipeline private** | `√r` as expression (perfect square or surd ratio). |
| `ratio_perfect_sqrt` | **Pipeline private** | perfect rational square root if exists. |
| `integer_perfect_sqrt` | **Pipeline private** | integer perfect square root via binary search. |
| `integrand_one_over` | **Pipeline private** | `integrand_one_over` |
| `algebraic_rt_one_over_x4_plus_one_smoke` | **Pipeline private** | `algebraic_rt_one_over_x4_plus_one_smoke` |
| `algebraic_rt_one_over_x4_plus_one_deriv` | **Pipeline private** | `algebraic_rt_one_over_x4_plus_one_deriv` |
| `algebraic_rt_one_over_x4_plus_four_via_res_smoke` | **Pipeline private** | `algebraic_rt_one_over_x4_plus_four_via_res_smoke` |
| `algebraic_rt_one_over_x4_plus_four_via_res_deriv` | **Pipeline private** | `algebraic_rt_one_over_x4_plus_four_via_res_deriv` |
| `algebraic_rt_one_over_x4_plus_x2_plus_one_via_res_smoke` | **Pipeline private** | `algebraic_rt_one_over_x4_plus_x2_plus_one_via_res_smoke` |
| `algebraic_rt_one_over_x4_plus_x2_plus_one_via_res_deriv` | **Pipeline private** | `algebraic_rt_one_over_x4_plus_x2_plus_one_via_res_deriv` |

### `risch/hermite.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `hermite_reduce` | **Stable** | Hermite-reduce `numer / factor^mult` w.r.t. `var`. |
| `x_var` | **Pipeline private** | `x_var` |
| `x` | **Pipeline private** | `x` |
| `hermite_x_over_x_squared_plus_one_squared_term_shape` | **Pipeline private** | `hermite_x_over_x_squared_plus_one_squared_term_shape` |
| `hermite_x_over_x_squared_plus_one_squared` | **Pipeline private** | `hermite_x_over_x_squared_plus_one_squared` |
| `hermite_one_over_x_fourth_plus_one_squared_terms` | **Pipeline private** | `hermite_one_over_x_fourth_plus_one_squared_terms` |
| `hermite_one_over_x_fourth_plus_one_squared` | **Pipeline private** | `hermite_one_over_x_fourth_plus_one_squared` |
| `hermite_one_over_x_fourth_minus_one_squared_rem` | **Pipeline private** | `hermite_one_over_x_fourth_minus_one_squared_rem` |
| `hermite_one_over_x_fourth_plus_one_fourth_power_rem` | **Pipeline private** | `hermite_one_over_x_fourth_plus_one_fourth_power_rem` |
| `hermite_mult_one_unchanged` | **Pipeline private** | `hermite_mult_one_unchanged` |
| `hermite_two_step_cubic_power` | **Pipeline private** | `hermite_two_step_cubic_power` |

### `risch/mod.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_risch` | **Partial** | `risch(f,x)` minimal subset; delegates to `integrate` (GIAC-217). **退役：** full Risch decision procedure. |
| `eval_risch_delegates_to_integrate` | **Pipeline private** | `eval_risch_delegates_to_integrate` |

### `risch/pow2expln.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `pow2expln` | **Stable** | `pow2expln(e, x)`; subset of GIAC `subst.cc::pow2expln(e, x)`. |
| `x` | **Pipeline private** | `x` |
| `pow2expln_x_to_x` | **Pipeline private** | `pow2expln_x_to_x` |
| `pow2expln_integer_power_unchanged` | **Pipeline private** | `pow2expln_integer_power_unchanged` |
| `pow2expln_x_to_const_exp` | **Pipeline private** | `pow2expln_x_to_const_exp` |
| `pow2expln_nested_in_mul` | **Pipeline private** | `pow2expln_nested_in_mul` |

### `risch/rothstein_trager.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `try_algebraic_rt_even_quartic` | **Partial** | algebraic RT for monic even quartics with constant numerator. **退役：** merge into full RT pipeline. |
| `try_integrate_x4_plus_one` | **Partial** | `∫ k/(x^4+1) dx` via algebraic RT conjugate pairing. **退役：** general even-quartic RT. |
| `rothstein_trager_integrate` | **Partial** | integrate `numer / factor` when `factor` is square-free and partfrac failed. **退役：** full Rothstein–Trager with algebraic extensions. |
| `x_var` | **Pipeline private** | `x_var` |
| `x_id` | **Pipeline private** | `x_id` |
| `integrand_one_over_x4_plus_one` | **Pipeline private** | `integrand_one_over_x4_plus_one` |
| `rothstein_one_over_x_fourth_plus_one_smoke` | **Pipeline private** | `rothstein_one_over_x_fourth_plus_one_smoke` |
| `rothstein_deriv_equals_integrand` | **Pipeline private** | `rothstein_deriv_equals_integrand` |
| `rothstein_one_over_x_squared_plus_one` | **Pipeline private** | `rothstein_one_over_x_squared_plus_one` |

### `risch/tower.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `rlvarx` | **Stable** | logarithmic / exponential extension variables in `expr` that depend on `var`. |
| `risch_tower` | **Stable** | returns the tower (most complex extension first) when `expr` is elementary over `var`. |
| `collect_rlvarx` | **Pipeline private** | collect `exp`/`ln` extension atoms depending on `var`. |
| `is_exp_or_ln` | **Pipeline private** | `exp` or `ln` top-level form. |
| `contains_non_elementary_transcendental` | **Pipeline private** | detect non-elementary transcendentals (e.g. trig). |
| `push_unique` | **Pipeline private** | append extension atom if not already present. |
| `extension_rank` | **Pipeline private** | nesting depth for tower ordering. |
| `x` | **Pipeline private** | `x` |
| `rlvarx_exp_and_ln` | **Pipeline private** | `rlvarx_exp_and_ln` |
| `risch_tower_rejects_trig` | **Pipeline private** | `risch_tower_rejects_trig` |
| `risch_tower_polynomial_is_empty` | **Pipeline private** | `risch_tower_polynomial_is_empty` |
| `rlvarx_nested_exp` | **Pipeline private** | `rlvarx_nested_exp` |
| `risch_tower_pow2expln_exp_x` | **Pipeline private** | `risch_tower_pow2expln_exp_x` |
| `risch_tower_rejects_nested_trig` | **Pipeline private** | `risch_tower_rejects_nested_trig` |
| `rlvarx_independent_of_var` | **Pipeline private** | `rlvarx_independent_of_var` |

### `series.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_series` | **Stable** | `series(f,var,center,order)` / `taylor(f,var,center,order)` (GIAC-216). |
| `is_plus_infinity` | **Pipeline private** | detect `+infinity` center symbol. |
| `parse_series_location` | **Pipeline private** | parse `(var, center, order)` from 2- or 3-arg tails. |
| `parse_var_center` | **Pipeline private** | `x=0` or `x` alone defaults center to 0. |
| `series_order_arg` | **Pipeline private** | evaluate and coerce series order to `usize`. |
| `taylor_series` | **Pipeline private** | Taylor at finite center (diff path then `series_at_center`). |
| `taylor_series_diff` | **Pipeline private** | Taylor via repeated differentiation and substitution. |
| `series_term` | **Pipeline private** | single Taylor term `(coeff/k!) * (x - center)^k`. |
| `eval_at` | **Pipeline private** | substitute `var ↦ center` and fold elementary values. |
| `is_zero_arg` | **Pipeline private** | zero-test for substitution argument. |
| `fold_elementary` | **Pipeline private** | fold `sin(0)`, `cos(0)`, `exp(0)` at series coefficients. |
| `is_zero` | **Pipeline private** | zero-test for series coefficient. |
| `factorial` | **Pipeline private** | integer factorial for Taylor denominators. |
| `series_preprocess_pow2expln` | **Pipeline private** | `series_preprocess_pow2expln` |
| `series_exp_at_zero` | **Pipeline private** | `series_exp_at_zero` |
| `taylor_sin_x_equals_zero` | **Pipeline private** | `taylor_sin_x_equals_zero` |
| `series_at_plus_infinity_exp` | **Pipeline private** | `series_at_plus_infinity_exp` |
| `taylor_at_one` | **Pipeline private** | `taylor_at_one` |
| `series_too_few_args` | **Pipeline private** | `series_too_few_args` |
| `series_order_zero` | **Pipeline private** | `series_order_zero` |
| `series_two_arg_var_symbol_defaults_center_zero` | **Pipeline private** | `series_two_arg_var_symbol_defaults_center_zero` |
| `series_ln_one_plus_x` | **Pipeline private** | `series_ln_one_plus_x` |
| `series_ck_int_65_cos_exp` | **Pipeline private** | `series_ck_int_65_cos_exp` |
| `series_cos_taylor_diff_fallback` | **Pipeline private** | `series_cos_taylor_diff_fallback` |
| `series_taylor_at_center_two` | **Pipeline private** | `series_taylor_at_center_two` |
| `eval_series_direct_bad_order` | **Pipeline private** | `eval_series_direct_bad_order` |

### `test_verify.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `assert_deriv_equals_integrand` | **Pipeline private** | `assert_deriv_equals_integrand` |
| `assert_series_equiv_at` | **Pipeline private** | `assert_series_equiv_at` |
| `assert_taylor_equiv_at` | **Pipeline private** | `assert_taylor_equiv_at` |
| `assert_deriv_equals_integrand_smoke` | **Pipeline private** | `assert_deriv_equals_integrand_smoke` |
