# GIAC-expr-api 单测分层审计表

**规范：** [test-writing-spec.md](../test-writing-spec.md)（A / A′ / B / C）  
**父 issue：** [GIAC-expr-api-tech-debt.md](GIAC-expr-api-tech-debt.md) §3A / §3B / §3D  
**状态图例：** ✅ 符合 · ⚠️ 部分符合 / 待改进 · ❌ 不符合 · ➖ 未审

**最后更新：** 2026-06-23（4B + 3B diff/integrate）

---

## 1. 汇总

| Crate / 文件 | 测试数 | A | A′ | B | C | ⚠️/❌ |
|--------------|--------|---|---|---|---|-------|
| `giac-solve` | 16 | 10 | 0 | 4 | 1 | 2 |
| `giac-calculus` / `exp_diff.rs` | 14 | 4 | 4 | 8 | 1 | 2 |
| `giac-calculus` / `preprocess.rs` | 6 | 2 | 2 | 2 | 1 | 1 |
| `giac-calculus` / `integrate.rs` | 24 | 6 | 0 | 14 | 4 | 4 |
| `giac-calculus` / `diff.rs` | 6 | 4 | 0 | 6 | 2 | 0 |
| `giac-calculus` / `eval_diff.rs` | 5 | 3 | 0 | 2 | 0 | 0 |
| `giac-core` / `alg_ext.rs` | 19 | 1 | 0 | 17 | 1 | 2 |
| `giac-core` / `poly_roots.rs` | 24 | 0 | 0 | 24 | 0 | 0 |
| **合计** | **114** | **30** | **6** | **81** | **11** | **11** |

> **3B 下一批：** `series.rs` / `risch/` / `integrate_heuristics.rs` 仍含语义 `contains`。

---

## 2. giac-solve

| 测试 | 层 | 被测契约 | 断言手段 | 状态 | 备注 |
|------|-----|----------|----------|------|------|
| `plugin::eval_solve_via_plugin` | A | `eval(Solve)` + plugin | display `[1]` | ✅ | 二次方程重根 |
| `solve::solve_quadratic_double_root` | A | `eval(Solve)` | `assert_equiv` → 1 | ✅ | |
| `solve::solve_linear_system_delegates_to_linsolve` | A | `eval(Solve)` → linsolve | `eval_const_expr` + 算术 | ⚠️ | 通路正确；验证用 **测试内 ad-hoc parser**，应迁入 `test_verify` 或 stability doc 约定 `linsolve` 有理输出 |
| `rootof::quadratic_rootof_has_two_branches` | B | `quadratic_rootof_roots` | `assert_roots_zero_poly` | ✅ | |
| `rootof::biquadratic_t_fourth_minus_two_has_four_roots` | B | `biquadratic_rootof_roots` | 个数 | ✅ | |
| `rootof::rational_quadratic_still_uses_roots` | B | `roots` 有理路径 | 个数 | ✅ | |
| `rootof::solve_t_squared_minus_two_uses_rootof` | A | `eval(Solve)` + rootof | `assert_equation_solutions` | ✅ | |
| `rootof::solve_t_fourth_minus_two_uses_rootof` | A | 同上 | 代入 + `contains_algext` | ✅ | |
| `froot::froot_linear_factor` | A | `eval(froot)` | `assert_froot_has_root` | ✅ | |
| `froot::froots_flanex_line195` | A | `eval(froots)` | `assert_froot_has_root` | ✅ | |
| `realroot::realroot_x_squared_minus_two` | A | `eval(realroot)` | `assert_is_algext_or_rootof` | ✅ | |
| `realroot::realroot_x_fourth_minus_one` | A | `eval(realroot)` | `assert_realroot_has` | ✅ | |
| `sturm::sturm_x_cubed_plus_one_squared` | A | `eval(sturm)` | 仅 `eval` 成功 | ⚠️ | 缺语义断言（序列长度/系数） |
| `sturm::sturm_x_cubed_plus_one` | A | `eval(sturm)` | 同上 | ⚠️ | 同上 |
| `sturm::sturmab_counts_root_in_interval` | A | `eval(sturmab)` | display `"1"` | ✅ | 区间根数；C 成分可接受 |
| `fsolve::fsolve_x_squared_minus_two` | A | `eval(fsolve)` | 近似/结构 | ➖ | 待按 stability doc 复审 |

**giac-solve 待办**

1. `test_verify::eval_const_expr`（或 `linsolve` 输出改 `Int`/`Rat`）——消 ⚠️ linear system  
2. `sturm` 两条补 **A** 语义（根数 / Sturm 序列项数）  
3. `giac-solve-api-stability.md` 写清 `linsolve` 输出形态

---

## 3. giac-calculus / `limit_engine/exp_diff.rs`

| 测试 | 层 | 被测契约 | 断言手段 | 状态 | 备注 |
|------|-----|----------|----------|------|------|
| `fold_mul_shifted_difference_gruntz` | B | `canonical_exp_diff` | `match_*` + `assert_equiv` ε | ✅ | |
| `rewrite_exp_minus_w_inv_matches_remove_lnexp` | B | `exp_scale_times_exp_minus_one` | 结构 + `assert_equiv` ε | ✅ | 不用 `match_exp_*`（scale 非 exp） |
| `fold_add_exp_difference` | B | `canonical_exp_diff` | 规范构造器 + `assert_equiv` | ✅ | |
| `balance_frac_minus_var_rewrites_inner_minus_x` | C | `balance_exp_arguments_frac_var` | display golden | ⚠️ | 仅 **C**；Pipeline 步骤，宜补 B（`assert_equiv` 等价 Frac）或标 `_display` |
| `first_order_vanishing_epsilon_rewrite` | B | `first_order_exp_vanishing_epsilon` | 否定 `match_*` + 结构 | ✅ | Partial tier |
| `limit_preprocessed_gruntz_minus_one` | A′ | preprocess + limit | `limit → -1` | ✅ | |
| `ratio_preprocess_mrv_cancels_opposing_exp_mul` | A′ | MRV preprocess + limit | `limit → 1` | ✅ | 原 contains 负断言已改为结果 |
| `limit_exp_over_exp_ck_int_ratio` | A′ | `limit_at_plus_infinity(ratio)` | `"1"` | ✅ | |
| `fold_ck_int_61_shape` | B + ⚠️ | `canonical_exp_diff` + preprocess | B: `match_*`；pre: `tree_contains_exp_of_sym` | ⚠️ | preprocess 部分非契约化；宜改 A′ `limit(pre)==1` 或删弱断言 |
| `fold_nested_gruntz_add_difference` | B | `detect_exp_difference_add`（Pipeline private） | `emit` + `match_*` | ⚠️ | 测 private；Stable 覆盖靠 `canonical_exp_diff` 其它用例 |
| `fold_nested_gruntz_full_expr` | B | `canonical_exp_diff` | `match_*` | ✅ | |
| `parse_neg_x_squared_in_exp` | B | parse → AST | 结构 `matches!` | ✅ | 解析契约 |
| `limit_nested_gruntz_preprocessed` | A′ | struct preprocess + limit | `"1"` | ✅ | |
| `exp_scale_times_exp_minus_one_shape` | B | `exp_scale_times_exp_minus_one` | 结构 + `assert_equiv` | ✅ | |

---

## 4. giac-calculus / `limit_engine/preprocess.rs`

| 测试 | 层 | 被测契约 | 断言手段 | 状态 | 备注 |
|------|-----|----------|----------|------|------|
| `merge_exp_quotient_frac` | C | `merge_exp_quotients` | 多备选 display | ⚠️ | 宜补 B：`assert_equiv` → `exp(-n)` |
| `preprocess_seven_pow_n_over_eight` | A′ | `limit_preprocess_plus_infinity` + limit | `limit → 0` | ✅ | |
| `preprocess_surd2pow_sqrt` | B | `surd2pow` | `assert_equiv` | ✅ | |
| `preprocess_sqrt_conjugate_minus_var` | B | `normalize_sqrt_conjugates` | `Expr::Frac` | ✅ | |
| `preprocess_nested_gruntz_exp_diff_factor` | A′ | 多步 preprocess + limit | `limit → 1` | ✅ | 内含 B 片段（`match_*`） |
| `preprocess_gruntz_exp_diff_factor` | A′ | factor + struct pre + limit | `limit → -1` | ✅ | |

---

## 5. giac-calculus / `integrate.rs`

| 测试 | 层 | 被测契约 | 断言手段 | 状态 | 备注 |
|------|-----|----------|----------|------|------|
| `integrate_reciprocal` 等规则类 | **B** + **A** | `integrate` + `eval(Integrate)` | `assert_equiv` | ✅ | |
| `integrate_one_minus_x_fourth` | **B** | partfrac smoke | `is_ok()` | ✅ | 不 pin 形态 |
| `integrate_const_over_quadratic` 等 | **B** + **C** | display golden | `assert_eq` format | ✅ | |
| `integrate_unsupported_*` | **B** | 错误路径 | `NotImplemented` | ✅ | |
| `giac223_*` 冒烟 | **B** | 不 crash | `is_ok()` | ⚠️ | |

---

## 5b. giac-calculus / `diff.rs` + `eval_diff.rs`

| 测试 | 层 | 被测契约 | 断言手段 | 状态 | 备注 |
|------|-----|----------|----------|------|------|
| `diff_x_squared` / `diff_sin_x_squared` / `diff_constant_*` | **B** + **A** | `diff` + `eval(Diff)` | `assert_equiv` | ✅ | |
| `diff_ln_times_x_squared` / `diff_frac_one_over_x` | **B** + **A** + **C** | `diff` + `eval` | display golden + `eval`≡`diff` | ✅ | quotient 形未 `assert_equiv` 化简 |
| `diff_atan_x` | **B** | 通路 | `is_ok` | ✅ | |
| `eval_diff_x_squared` / `eval_derive_multivariate` | **A** | `eval(Diff/Derive)` | `assert_equiv` | ✅ | |
| `eval_diff_too_few_args` 等 | **B** | 错误路径 | `is_err` | ✅ | |

---

## 6. giac-core / `algebra/alg_ext.rs`（3D 部分）

| 测试 | 层 | 被测契约 | 断言手段 | 状态 | 备注 |
|------|-----|----------|----------|------|------|
| `algext_mul_squares_to_two` 等运算 | B | `AlgExt` 环运算 | `assert` / 等价 | ✅ | |
| `algext_sqrt_of_neg_sqrt2_is_complex` | B | 分支 + `eq_mod` | `assert_equiv` 分支² | ✅ | 3D 已修 |
| `algext_to_rootof_roundtrip_display` | B | `try_as_algext_data` | `eq_mod` | ✅ | |
| `algext_frac_via_eval` | A | `eval` 含 AlgExt | 通路 | ✅ | |
| 其余 `fold_algext_*` / `embed_*` | B | tower / fold API | 坐标 / 等价 | ✅ | |

---

## 7. giac-core / `algebra/poly_roots.rs`（3D）

| 测试类 | 层 | 被测契约 | 断言手段 | 状态 | 备注 |
|--------|-----|----------|----------|------|------|
| 二次 / 三次（`roots_quadratic_*`、`roots_cubic_*`、`roots_x3_minus_x_plus_1_vanish`） | B | `poly_algext_roots` | `verify_root` / `roots_all_vanish` | ✅ | |
| resolvent（`resolvent_*`、`f2_resolvent_split_*`） | B | resolvent 阶段 | `verify_root` + 个数 / `eq_mod` | ✅ | |
| 四次 `t⁴+t+1`（`roots_quartic_t4_plus_t_plus_1`、`euler_*`、`depressed_*`） | B | quartic + FieldSession | `verify_root` / `eq_mod` / 维数 ≤24 | ✅ | 无 `#[ignore]` |
| `adjoin_sqrt_*` / `field_session_dimension_bound_*` | B | `FieldSession` 契约 | `eq_mod` / `dim(L)` 断言 | ✅ | G4 门禁 |
| 双二次 `roots_biquadratic_t4_minus_2` | B | deg-4 biquadratic | `verify_root` | ✅ | |

**poly_roots 待办**

1. 可选补 **A**：`eval` 经 `poly_alg_from_expr` → `poly_algext_roots_for_ctx` 端到端（当前全为库函数 **B**）
2. 逐测行展开（上表为分类汇总；新增 `#[test]` 须拆行登记）

---

## 8. 基础设施阻塞项（审计结论）

| ID | 问题 | 影响测试 | 建议 |
|----|------|----------|------|
| **T1** | `giac-solve::xcas_default()` 无 simplify plugin | A 层 `assert_equiv` 对 `2*1^-1` 失败 | 测试 Context helper 或 `linsolve` 输出规范化 |
| **T2** | `match_exp_times_exp_minus_one` 仅认 `exp(L)*(exp(ε)-1)` | 误用于 `scale*(exp(ε)-1)` | 文档已澄清；用 `exp_minus_one_epsilon` |
| **T3** | `diff(integrate(f))` 未覆盖 `ln(abs)` 等 | integrate 无法用 A 验收 | 登记 calculus gap；暂 C / `assert_equiv` |
| **T4** | 缺统一 `test_verify`（calculus/core） | ad-hoc helper 扩散 | 按 solve 模式抽取；`poly_roots` 已用 `verify_root` 内聚 |

---

## 9. 维护

- 新增/修改 `#[test]`：**必须**更新本表对应行。  
- 3C（conformance 层）单独审计，见 tech-debt **GIAC-expr-api-3C**。  
- 全 crate 推广：复制 §2 表头，按模块追加小节。
