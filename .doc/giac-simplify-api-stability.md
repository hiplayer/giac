# giac-simplify API 稳定性分层

**规范来源:** [algorithm-expr-api.md](algorithm-expr-api.md)  
**上游缺口:** [issues/GIAC-simplify-poly-upstream-gaps.md](issues/GIAC-simplify-poly-upstream-gaps.md) §1  
**代码注释约定:** 见下文 §1；`lib.rs` `//!` 含本文件子表摘要。

---

## 1. 源码注释格式

| 标记 | 用于 | 可见性 |
|------|------|--------|
| `/// **Stable** — …` | 有 I/O 契约、可跨模块调用、须有单测 | `pub` |
| `/// **Stable (bounded)** — …` | 稳定但范围受 giac-poly / 输入类限制 | `pub` |
| `/// **Partial** — …` | 规则表子集；注释写 **退役** 或扩展计划 | `pub` |
| `// **Temporary** — …` | `shim_*` / AlgExt stub；须写退役 issue | `fn` 私有 |
| `// **Pipeline private** — …` | 模块内递归/形状检测 | `fn` 私有 |

**未标注的私有函数:** 不应出现 — 运行 `giac-rs/scripts/annotate_api_tiers.py` 为每个 `fn` 补 `// **Pipeline private**`；Per-file 全表见本文末尾 **Per-file function inventory**。

**提交前复审（测试全绿后）：** 见 [algorithm-expr-api.md §7.2](algorithm-expr-api.md#72-测试通过后提交--合入前复审) — 检查临时匹配是否净减少、新增 `fn` tier 是否已写入本文 §2–§4。

---

## 2. Crate 公开 API（`lib.rs` re-export）

| 函数 / 类型 | 层级 | 模块 | 退役 / 边界 |
|-------------|------|------|-------------|
| `normal` | **Stable** | `expand` | 超越叶子 `expr_to_poly` 失败时原样返回 |
| `expand` | **Stable (legacy default)** | `expand` | 同 `expand_with_policy(..., Full)` |
| `expand_with_policy` | **Stable** | `expand` | `NoExpDistribute` 供 limit 管线 |
| `expand_polynomial` | **Stable** | `expand` | 不 distribute 过含 exp/ln 的 Add |
| `ratnormal` | **Stable** | `ratnormal` | 非有理子树 → `TypeError` |
| `factor` | **Stable (bounded)** | `factor` | 多变量 Hensel 缺口在 giac-poly FAC-G* |
| `ifactor` | **Stable** | `ifactor` | 大整数可能慢 |
| `assert_equiv` | **Stable** | `equiv` | 经 `normal`；radical 形态经 `canonical_radical_expr`（§3） |
| `is_zero` | **Stable** | `equiv` | |
| `sub` | **Stable** | `equiv` | |
| `texpand` | **Partial** | `trig` | 负整数倍角 → `NotImplemented` |
| `lin` | **Partial** | `trig` | 有界 exp 幂；非完整 `lin.cc` |
| `halftan` | **Partial** | `trig` | 单形状 `sin(2x)/(1+cos(2x))` |
| `install_simplify` | **Stable** | `plugin` | |
| `xcas_default` | **Stable** | `plugin` | |
| `DefaultAlgebraPlugin` | **Stable** | `plugin` | |

### 2.1 未实现的上游 API（非本 crate 导出）

规划见 `builtin-api-map.md`，**无 Rust 实现:**  
`tlin`, `trig2exp`, `trigsin`, `trigcos`, `tan2sincos`, `tcollect`, `tsimplify`, `lncollect`, 三角互化族, `reorder`, `truncate`, `epsilon2zero`, `non_recursive_normal`.

`simplify` 在 **`giac-core`**（AST flatten），语义 ≠ upstream `subst.cc`。

---

## 3. 子模块 Pipeline / 内部 API

| 模块 | Pipeline private（代表性） | 说明 |
|------|---------------------------|------|
| `expand.rs` | `expand_mul_pair`, `repeated_mul`, `expand_binomial`, `expr_contains_exp_ln` | 展开递归 |
| `ratnormal.rs` | `rational_parts`, `reduce_fraction`, `rational_add`, `rational_mul` | 有理式收集 |
| `factor.rs` | `factor_expr`, `factor_poly_form`, `flatten_mul`, `rational_num_den` | Expr→Poly→factor |
| `trig.rs` | `texpand_rec`, `expand_sin_arg`, `expand_cos_arg`, `integer_multiple`, `lin_rec`, `detect_halftan_tan` | 三角/指数规则表 |
| `ifactor.rs` | 内部试除循环 | |
| `equiv.rs` | `canonical_radical_expr`, `inv_sqrt_to_mul` (Stable crate-internal), `sub`/`is_zero`/`assert_equiv` | radical 预归一 + `normal` 等价判定 |

---

## 4. 临时 API（Temporary）

| 函数 | 文件 | 类型 | 退役条件 |
|------|------|------|----------|
| `ratnormal_algext` | `ratnormal.rs` | **shim** | [GIAC-algext-adoption](issues/GIAC-algext-adoption.md) A-03：`ext_reduce` 有理化 |
| `try_factor_quadratic_sqrt` | `factor.rs` | **Partial 内部** | `ctx.with_sqrt` 下二次 sqrt 分解 |

**已退役（2026-06）：** `try_factor_via_algext` / `try_factor_quadratic_rootof` 不得再接回 `factor()` 主路径；扩域线性分解见 `solve` / `factor_into_algext`（非通用 `factor` 展示）。

**已升格（2026-06-27，GIAC-expr-api-2C）：** `canonical_radical_expr` / `inv_sqrt_to_mul` 从 Temporary drift 升为 **Stable (crate-internal)**；I/O 契约见 `equiv.rs` 源码 doc。HITL 决策：选 A（升格），理由：命名已符 `canonical_*` 规约、行为不变、不波及 `normal`/golden、工作量最小；选项 B（并入 `normal`）会改 `normal` 行为且与 `ratnormal` 域职责重叠。

**禁止:** 在 `giac-calculus` / `limit_engine` 复制上述 drift 逻辑；应调用 `assert_equiv` 或扩展 owning 模块稳定 API。

---

## 5. 跨 crate 契约

| 调用方 | 使用的 Stable API | 注意 |
|--------|-------------------|------|
| `giac-calculus::limit_engine` | `ratnormal`, `normal`, `expand`（间接） | limit 前勿 `expand(Full)` 破坏 exp 差分形 |
| `giac-calculus::integrate` | `expand` | |
| `giac-conformance` | `assert_equiv` | L1 属性见 [conformance-testing.md §3.5](conformance-testing.md#35-命令-io-契约l1-normative) |
| `giac-simplify::factor` | `giac-poly::factor_into` | ℚ 上不可约因子保持多项式形；见 §5.1 |

### 5.1 `factor(expr)` — I/O 契约（normative）

| 字段 | 内容 |
|------|------|
| **Context** | `xcas_default()`（含 simplify plugin） |
| **输入** | 有理系数多项式 `Expr`；允许未展开 `Mul`/`Pow`；无 `AlgExt` 系数（有则 `factor_algext_form`） |
| **输出（数学）** | `expand(out) ≡ expand(in)`（与 [sympy_verify `factor`](../../giac-rs/tests/conformance/scripts/sympy_verify.py) 一致） |
| **输出（形态）** | 因子为 **ℚ[x₁,…,xₙ]** 中多项式；ℚ 上不可约则 **保留不可约因子** |
| **禁止** | 通用 `factor` 默认路径输出 `rootof` 线性因子；栈溢出 / UB |
| **L1 门禁** | `cargo nextest run --release -p giac-conformance --test giac_check_factor`（**不可**为绿而弱化 `sympy_verify.py`） |
| **L2 参考** | `giac-2.0.0/check/testfactor` + `factor.out`（字面参考；与 L1 冲突时 L1 + [known-divergences.md](known-divergences.md)） |
| **上游 C++** | `ezgcd.cc` / `gausspol.cc`（ℚ 分解；非 `ext_factor` 展示路径） |
| **边界** | `rootof` / 扩域分裂 → `solve`、`roots`、`factor_into_algext`（单独契约） |

**动 `factor` 实现：** PR 须引用本节；L1 失败见 [conformance-testing.md §3.6](conformance-testing.md#36-l1-失败处理须人工确认)。

---

## 6. 维护

新增 `pub` 函数时：

1. 源码加 `/// **Stable|Partial** — …`（或运行下方脚本自动补 **Pipeline private**）
2. 更新本文件 §2 表
3. 运行 `python3 scripts/annotate_api_tiers.py --inventory` 刷新 §Per-file inventory
4. 临时 shim/drift 写退役 issue

**批量标注 / 刷新清单：**

```bash
cd giac-rs
python3 scripts/annotate_api_tiers.py              # 为新 fn 补 tier 注释（幂等）
python3 scripts/annotate_api_tiers.py --inventory # 刷新本文 Per-file 表
```

---

## Per-file function inventory (generated)

Regenerate: `python3 scripts/annotate_api_tiers.py --inventory`

### `equiv.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `sub` | **Stable** | construct `a - b` as an expression tree. |
| `is_zero` | **Stable** | true when `normal(e)` simplifies to zero. |
| `assert_equiv` | **Stable** | true when `a` and `b` are mathematically equivalent under `normal`. |
| `canonical_radical_expr` | **Stable** | `canonical_radical_expr`: 1/sqrt(n) → (1/n)*sqrt(n) pre-normalizer for assert_equiv |
| `inv_sqrt_to_mul` | **Stable** | `inv_sqrt_to_mul`: Pow(Sqrt(n),-1) → Mul[Rat(1,n),Sqrt(n)] helper of canonical_radical_expr |
| `equiv_commutative_add` | **Pipeline private** | `equiv_commutative_add` |
| `equiv_expanded_square` | **Pipeline private** | `equiv_expanded_square` |
| `equiv_sqrt_half_forms` | **Pipeline private** | `equiv_sqrt_half_forms` |
| `canonical_radical_expr_inv_sqrt_form` | **Pipeline private** | `canonical_radical_expr_inv_sqrt_form` |
| `equiv_rational_reduced` | **Pipeline private** | `equiv_rational_reduced` |
| `is_zero_complex` | **Pipeline private** | `is_zero_complex` |
| `sub_fraction_difference` | **Pipeline private** | `sub_fraction_difference` |
| `equiv_mod_difference` | **Pipeline private** | `equiv_mod_difference` |
| `not_equiv_different` | **Pipeline private** | `not_equiv_different` |
| `sqrt2_minpoly` | **Pipeline private** | `sqrt2_minpoly` |
| `equiv_algext_same_field` | **Pipeline private** | `equiv_algext_same_field` |
| `equiv_algext_square_minus_two` | **Pipeline private** | `equiv_algext_square_minus_two` |

### `expand.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `expand` | **Stable** | distribute products over sums and expand powers of sums (`ExpandPolicy::Full`). |
| `expand_with_policy` | **Stable** | [`expand`] with an explicit policy. |
| `expand_polynomial` | **Stable** | polynomial-oriented expand: no distribution through `exp`/`ln` subtrees. |
| `expr_contains_exp_ln` | **Pipeline private** | used by `ExpandPolicy::NoExpDistribute`. |
| `expand_mul_pair` | **Pipeline private** | distribute one mul factor over add |
| `expand_pow` | **Pipeline private** | expand integer powers and mod-poly powers |
| `repeated_mul` | **Pipeline private** | repeated multiply for small integer power |
| `expand_binomial` | **Pipeline private** | binomial power via poly or Expr::pow |
| `normal` | **Stable** | expand then collect into polynomial normal form. |
| `modulus_from_expr` | **Pipeline private** | coerce Expr modulus to i64 |
| `normal_mod_power_displays_giac_style` | **Pipeline private** | `normal_mod_power_displays_giac_style` |
| `expand_square_of_sum` | **Pipeline private** | `expand_square_of_sum` |
| `normal_cancels_distributed_polynomial_terms` | **Pipeline private** | `normal_cancels_distributed_polynomial_terms` |
| `expand_distribute_mul_over_add` | **Pipeline private** | `expand_distribute_mul_over_add` |
| `expand_polynomial_distribute_rational_product` | **Pipeline private** | `expand_polynomial_distribute_rational_product` |
| `expand_polynomial_skips_exp_subtree` | **Pipeline private** | `expand_polynomial_skips_exp_subtree` |
| `expand_polynomial_skips_when_other_factor_has_exp` | **Pipeline private** | `expand_polynomial_skips_when_other_factor_has_exp` |
| `expand_pow_single_base` | **Pipeline private** | `expand_pow_single_base` |
| `expand_pow_cube_of_symbol` | **Pipeline private** | `expand_pow_cube_of_symbol` |
| `expand_complex_and_non_poly_normal` | **Pipeline private** | `expand_complex_and_non_poly_normal` |
| `expand_rhs_add_and_zero_power` | **Pipeline private** | `expand_rhs_add_and_zero_power` |
| `normal_fast_path_already_polynomial` | **Pipeline private** | `normal_fast_path_already_polynomial` |
| `expand_binomial_fallback_for_non_poly` | **Pipeline private** | `expand_binomial_fallback_for_non_poly` |
| `expand_binomial_fallback_semantic` | **Pipeline private** | `expand_binomial_fallback_semantic` |
| `expand_recurses_into_frac` | **Pipeline private** | `expand_recurses_into_frac` |

### `factor.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `factor` | **Stable (bounded)** | structural (`Mul`/`Pow`/`Frac`) then polynomial factorization. |
| `factor_expr` | **Pipeline private** | recursive factor on Mul/Pow/Frac |
| `flatten_mul` | **Pipeline private** | flatten Mul to factor vec |
| `factor_poly_form` | **Pipeline private** | normal→poly→factor_into chain |
| `factor_algext_form` | **Pipeline private** | path B: factor in K[var] when expr has algebraic coefficients. |

### `ifactor.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `ifactor` | **Stable** | factor `n` into primes: `p1^e1 * p2^e2 * ...`. |
| `is_probable_prime_u64` | **Pipeline private** | Miller-Rabin for u64 |
| `mod_mul` | **Pipeline private** | u64 modular multiply |
| `mod_pow` | **Pipeline private** | u64 modular exponentiation |
| `pollard_rho` | **Pipeline private** | Pollard rho on BigInt |
| `pollard_rho_u64` | **Pipeline private** | Pollard rho on u64 |
| `ifactor_small` | **Pipeline private** | `ifactor_small` |

### `plugin.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `normal` | **Stable** | expand then polynomial collect |
| `ratnormal` | **Stable** | single fraction in lowest terms |
| `expand` | **Stable** | expand with Full policy |
| `factor` | **Stable (bounded)** | structural then giac-poly factor |
| `ifactor` | **Stable** | integer prime factorization display |
| `texpand` | **Partial** | trig/exp/ln arg expand then algebraic expand |
| `halftan` | **Partial** | half-angle tan on narrow pattern |
| `lin` | **Partial** | linearize exp products and (exp+1)^2 |
| `install_simplify` | **Stable** | install [`DefaultAlgebraPlugin`] on `ctx`. |
| `xcas_default` | **Stable** | `giac-core::Context::xcas_default()` with simplification enabled. |
| `eval_normal_via_plugin` | **Pipeline private** | `eval_normal_via_plugin` |
| `eval_expand_binomial_via_plugin` | **Pipeline private** | `eval_expand_binomial_via_plugin` |
| `eval_factor_via_plugin` | **Pipeline private** | `eval_factor_via_plugin` |
| `eval_factor_x_squared_minus_two_irreducible` | **Pipeline private** | `eval_factor_x_squared_minus_two_irreducible` |
| `eval_factor_x_squared_minus_two_factorization` | **Pipeline private** | `eval_factor_x_squared_minus_two_factorization` |
| `eval_factor_x_squared_minus_two_stays_irreducible` | **Pipeline private** | `eval_factor_x_squared_minus_two_stays_irreducible` |

### `ratnormal.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `ratnormal` | **Stable** | normalize a rational expression to a single fraction in lowest terms. |
| `ratnormal_algext` | **Temporary** | shim: `AlgExt` ratnormal delegates to `eval` (GIAC-algext-adoption A-03). |
| `rational_parts` | **Pipeline private** | Expr to (num Poly, den Poly) |
| `rational_add` | **Pipeline private** | add rationals with common denominator |
| `reduce_fraction` | **Pipeline private** | gcd-reduce num/den Poly pair |
| `ratnormal_frac_form` | **Pipeline private** | `ratnormal_frac_form` |
| `ratnormal_negative_power` | **Pipeline private** | `ratnormal_negative_power` |
| `ratnormal_product_of_fractions` | **Pipeline private** | `ratnormal_product_of_fractions` |
| `ratnormal_polynomial_only` | **Pipeline private** | `ratnormal_polynomial_only` |
| `ratnormal_sum_of_fractions` | **Pipeline private** | `ratnormal_sum_of_fractions` |
| `ratnormal_division_by_zero` | **Pipeline private** | `ratnormal_division_by_zero` |
| `ratnormal_zero_numerator` | **Pipeline private** | `ratnormal_zero_numerator` |
| `ratnormal_empty_sum` | **Pipeline private** | `ratnormal_empty_sum` |
| `ratnormal_not_rational_error` | **Pipeline private** | `ratnormal_not_rational_error` |
| `ratnormal_non_integer_power_error` | **Pipeline private** | `ratnormal_non_integer_power_error` |
| `ratnormal_reduces_common_factor` | **Pipeline private** | `ratnormal_reduces_common_factor` |
| `ratnormal_reduces_common_factor_semantic` | **Pipeline private** | `ratnormal_reduces_common_factor_semantic` |
| `ratnormal_algext_square_minus_two` | **Pipeline private** | `ratnormal_algext_square_minus_two` |

### `test_verify.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `assert_factorization` | **Pipeline private** | `assert_factorization` |
| `assert_factorization_eval` | **Pipeline private** | `assert_factorization_eval` |

### `trig.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `texpand` | **Partial** | expand `sin`/`cos`/`exp`/`ln` arguments, then algebraically expand. |
| `texpand_rec` | **Pipeline private** | recursive texpand on Expr tree |
| `expand_sin_arg` | **Pipeline private** | sin angle-sum and n*x rules |
| `expand_cos_arg` | **Pipeline private** | cos angle-sum and n*x rules |
| `expand_exp_arg` | **Pipeline private** | exp of sum → product of exp |
| `expand_ln_arg` | **Pipeline private** | ln of product → sum of ln |
| `expand_sin_nx` | **Pipeline private** | sin(nx) for small integer n |
| `expand_cos_nx` | **Pipeline private** | cos(nx) for small integer n |
| `halftan` | **Partial** | half-angle tangent substitution on a narrow rational-trig pattern. |
| `halftan_half_angle_rational` | **Pipeline private** | Weierstrass tan(v/2) form |
| `tan_half` | **Pipeline private** | tan(v/2) Expr builder |
| `lin` | **Partial** | linearize exponentials: `(exp(x)+1)^2`, `exp(a)*exp(b)`, bounded `exp`-base powers. |
| `lin_rec` | **Pipeline private** | recursive lin on Expr tree |
| `contains_exp` | **Pipeline private** | subtree contains exp |
| `expand_integer_pow` | **Pipeline private** | expand exp-base integer power |
| `sin_expr` | **Pipeline private** | build sin Expr |
| `cos_expr` | **Pipeline private** | build cos Expr |
| `integer_multiple` | **Pipeline private** | detect n*x integer multiple |
| `detect_halftan_tan` | **Pipeline private** | detect sin(2x)/(1+cos(2x)) |
| `as_frac` | **Pipeline private** | view Expr as Frac pair |
| `unwrap_unit_mul_owned` | **Pipeline private** | peel unit coefficient from mul |
| `sin_double_angle` | **Pipeline private** | detect sin(2k*x) in halftan |
| `cos_double_angle` | **Pipeline private** | detect cos(2k*x) in halftan |
| `lin_exp_plus_one_pow` | **Pipeline private** | expand (exp+1)^2 only |
| `ctx` | **Pipeline private** | `ctx` |
| `texpand_cos_sum` | **Pipeline private** | `texpand_cos_sum` |
| `texpand_cos_sum_semantic` | **Pipeline private** | `texpand_cos_sum_semantic` |
| `texpand_cos_triple_angle` | **Pipeline private** | `texpand_cos_triple_angle` |
| `detect_halftan_direct` | **Pipeline private** | `detect_halftan_direct` |
| `halftan_sin_over_one_plus_cos` | **Pipeline private** | `halftan_sin_over_one_plus_cos` |
| `lin_exp_square` | **Pipeline private** | `lin_exp_square` |
