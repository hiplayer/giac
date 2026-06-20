# GIAC-expr-api 单测分层审计表

**规范：** [test-writing-spec.md](../test-writing-spec.md)（A / A′ / B / C）  
**父 issue：** [GIAC-expr-api-tech-debt.md](GIAC-expr-api-tech-debt.md) §3A / §3B / §3D  
**状态图例：** ✅ 符合 · ⚠️ 部分符合 / 待改进 · ❌ 不符合 · ➖ 未审

**最后更新：** 2026-06-20（3A/3B 首批落点）

---

## 1. 汇总

| Crate / 文件 | 测试数 | A | A′ | B | C | ⚠️/❌ |
|--------------|--------|---|---|---|---|-------|
| `giac-solve` | 16 | 10 | 0 | 4 | 1 | 2 |
| `giac-calculus` / `exp_diff.rs` | 14 | 4 | 4 | 8 | 1 | 2 |
| `giac-calculus` / `preprocess.rs` | 6 | 2 | 2 | 2 | 1 | 1 |
| `giac-calculus` / `integrate.rs` | 24 | 0 | 0 | 8 | 10 | 14 |
| `giac-core` / `alg_ext.rs` | 19 | 1 | 0 | 17 | 1 | 2 |
| **合计** | **79** | **17** | **6** | **39** | **14** | **21** |

> **integrate.rs** 多数为直接调 `integrate()`（**B**），几乎无 `eval(Integrate)`（**A**）——与规范 §3.1 差距最大，见 §5 阻塞项。

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
| `integrate_reciprocal` | B | `integrate` | `assert_equiv` vs `ln_abs` | ⚠️ | 语义 ✅；缺 **A** `eval(Integrate)` |
| `integrate_one` / `integrate_x_*` / `integrate_const_*` | B + C | `integrate` | display golden | ⚠️ | 宜 `assert_equiv` 或标 C |
| `integrate_constant_and_sum` | B + C | `integrate` | 常数 display + 和 display | ⚠️ | 和项 `ln(abs(x))+1*x` 为 C |
| `integrate_one_over_one_plus_x_squared` 等 exact golden | B + C | `integrate` | display | ⚠️ | 已知形态；可保留 C |
| `integrate_one_minus_x_fourth` 等 | C | `integrate` | `contains` | ❌ | 待 `assert_equiv` 或 A′ + gap |
| `integrate_x_over_x_squared_plus_one` 等 | C | `integrate` | `contains` | ❌ | `diff` 回代未就绪 |
| `integrate_unsupported_returns_not_implemented` | B | 错误路径 | `NotImplemented` | ✅ | |
| `integrate_not_implemented_messages` | B | 错误消息 | 子串 | ⚠️ | 错误消息 C 可接受 |
| `integrate_definite_bounds` | B | 定积分 | 结构/值 | ➖ | 待复审 |
| `giac223_tanh_exp_frac` / `giac223_exp_over_linear` | B | 不 crash | `is_ok()` | ⚠️ | 冒烟；缺语义 |
| **（缺）** | **A** | `eval(Integrate)` via plugin | — | ❌ | **全文件无 A**；见 `plugin::eval_integrate_via_plugin` 仅 1 条 |

**integrate 待办（3B 最大缺口）**

1. 每个规则类用例 duplicate：**A**（`eval(Integrate)` + `assert_equiv`）+ 保留 **B** 或直接 `integrate`  
2. `contains` 用例 → `assert_integrate_matches` 或登记 [GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md) + `diff` gap  
3. 冒烟测试改为最小 **A** 或移到 heuristics 模块

---

## 6. giac-core / `algebra/alg_ext.rs`（3D 部分）

| 测试 | 层 | 被测契约 | 断言手段 | 状态 | 备注 |
|------|-----|----------|----------|------|------|
| `algext_mul_squares_to_two` 等运算 | B | `AlgExt` 环运算 | `assert` / 等价 | ✅ | |
| `algext_sqrt_of_neg_sqrt2_is_complex` | B | 分支 + `eq_mod` | `assert_equiv` 分支² | ✅ | 3D 已修 |
| `algext_to_rootof_roundtrip_display` | B | `try_as_algext_data` | `eq_mod` | ✅ | |
| `algext_frac_via_eval` | A | `eval` 含 AlgExt | 通路 | ✅ | |
| 其余 `fold_algext_*` / `embed_*` | B | tower / fold API | 坐标 / 等价 | ✅ | Blocker：FieldSession（1A）未覆盖 `poly_roots` |

`poly_roots.rs` 单测：**➖** 待 1A 后单独增表。

---

## 7. 基础设施阻塞项（审计结论）

| ID | 问题 | 影响测试 | 建议 |
|----|------|----------|------|
| **T1** | `giac-solve::xcas_default()` 无 simplify plugin | A 层 `assert_equiv` 对 `2*1^-1` 失败 | 测试 Context helper 或 `linsolve` 输出规范化 |
| **T2** | `match_exp_times_exp_minus_one` 仅认 `exp(L)*(exp(ε)-1)` | 误用于 `scale*(exp(ε)-1)` | 文档已澄清；用 `exp_minus_one_epsilon` |
| **T3** | `diff(integrate(f))` 未覆盖 `ln(abs)` 等 | integrate 无法用 A 验收 | 登记 calculus gap；暂 C / `assert_equiv` |
| **T4** | 缺统一 `test_verify`（calculus/core） | ad-hoc helper 扩散 | 按 solve 模式抽取 |

---

## 8. 维护

- 新增/修改 `#[test]`：**必须**更新本表对应行。  
- 3C（conformance 层）单独审计，见 tech-debt **GIAC-expr-api-3C**。  
- 全 crate 推广：复制 §2 表头，按模块追加小节。
