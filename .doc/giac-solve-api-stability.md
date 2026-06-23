# giac-solve API 稳定性分层

**规范来源:** [algorithm-expr-api.md](algorithm-expr-api.md)  
**Expr ↔ Poly:** [expr-poly-conversion.md](expr-poly-conversion.md)  
**算法缺口:** [issues/GIAC-algorithm-gaps-open.md](issues/GIAC-algorithm-gaps-open.md) SOL-G*  
**代码:** `giac-rs/crates/giac-solve/src/*`

---

## 1. 源码注释格式

| 标记 | 用于 | 可见性 |
|------|------|--------|
| `/// **Stable** — …` | builtin 求值入口 | `pub` |
| `/// **Stable (bounded)** — …` | 稳定但方程类/次数受限 | `pub` |
| `/// **Partial** — …` | 数值/超越 stub | `pub` |
| `// **Pipeline private** — …` | 方程→多项式、形状检测 | `fn` 私有 |

**未标注私有 fn:** `python3 scripts/annotate_api_tiers.py`；Per-file 表见本文末尾。

---

## 2. Crate 公开 API（`lib.rs` re-export）

| 函数 | 层级 | 模块 | 边界 |
|------|------|------|------|
| `eval_solve` | **Stable (bounded)** | `solve` | 一元多项式 + 线性系统 + 窄超越 |
| `eval_froot` | **Stable (bounded)** | `froot` | 有理根 |
| `eval_realroot` | **Stable (bounded)** | `realroot` | Sturm 实根隔离 |
| `eval_sturm` | **Stable** | `sturm` | Sturm 序列 |
| `eval_sturmab` | **Stable** | `sturm` | 区间内根个数 |
| `eval_fsolve` | **Partial** | `fsolve` | Newton 数值 stub |
| `eval_linsolve` | **Stable** | `giac_linalg` re-export | 线性系统 |
| `quadratic_rootof_roots` | **Stable (bounded)** | `rootof` | 二次 `rootof` 两支 |
| `biquadratic_rootof_roots` | **Partial** | `rootof` | 仅双二次；一般四次 NotImplemented |
| `install_solve` / `xcas_default` | **Stable** | `plugin` | |

---

## 3. I/O 契约

### `eval_solve`

| 字段 | 说明 |
|------|------|
| **输入** | `solve(equation, var)` 或 `solve([eqs], [vars])` |
| **输出** | `List<Expr>` 解 |
| **表示** | 低次：`Poly::roots` → `poly_to_expr`；代数扩域：`poly_algext_roots_for_ctx`；二次/双二次 fallback：`rootof` |
| **前提** | 方程经 `eval` 后可 `expr_to_poly`（ℚ 多项式）；`ctx` 传入 session cache（R5） |
| **禁止** | 对含超越式未走 MRV/级数路径而强行 poly |

### `poly_roots_as_exprs` / `try_algext_or_rootof_roots`（pipeline）

| 字段 | 说明 |
|------|------|
| **输入** | `&Poly`, `&Var`, `&Context` |
| **路径** | `giac_poly::roots` → 失败时 `poly_algext_roots_for_ctx(p, var, ctx)` → `quadratic_rootof_roots` / `biquadratic_rootof_roots` |
| **缓存** | `ctx.session()` extension cache（与 eval 同线程 `Context` 共享） |

### `quadratic_rootof_roots` / `biquadratic_rootof_roots`

| 字段 | 说明 |
|------|------|
| **输入** | `&Poly`（ℚ）, `&Var` |
| **输出** | `Vec<ExprArc>`（`rootof` / `AlgExt` 形） |
| **边界** | 双二次要求奇次项为零；一般四次 → `NotImplemented` |

---

## 4. Pipeline private（代表性）

| 模块 | 函数 | 说明 |
|------|------|------|
| `solve` | `equation_to_poly`, `try_transcendental_solve` | lhs-rhs → poly；`sin(x)=0` 等 |
| `solve` | `poly_roots_as_exprs`, `try_algext_or_rootof_roots` | ℚ roots → `poly_algext_roots_for_ctx(ctx)` |
| `froot` | 有理根搜索链 | 依赖 `giac-poly` |
| `rootof` | `rootof_expr` | Expr 构造 |

---

## 5. 演化 issue

| 缺口 | 说明 |
|------|------|
| 一般四次 | 经 `poly_algext_roots_for_ctx`（[GIAC-poly-p3-6](issues/GIAC-poly-p3-6-quartic-roots-gaps.md)）；e2e perf 见 `#[ignore]` |
| `froots` / assume | [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md) SOL-G1/G5 |
| 测试契约 | [GIAC-expr-api-tech-debt](issues/GIAC-expr-api-tech-debt.md) 3A |

---

## 6. 维护

```bash
cd giac-rs
python3 scripts/annotate_api_tiers.py
python3 scripts/annotate_api_tiers.py --inventory
```

---

## Per-file function inventory (generated)

Regenerate: `python3 scripts/annotate_api_tiers.py --inventory`

### `froot.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_froot` | **Stable (bounded)** | rational roots of univariate poly |
| `parse_froot_args` | **Pipeline private** | `parse_froot_args` |
| `append_factor_roots` | **Pipeline private** | `append_factor_roots` |
| `solve_factor_roots` | **Pipeline private** | S3: deg≤4 via S0 kernel; deg≥5 rootof branch |
| `froots_flanex_line195` | **Pipeline private** | `froots_flanex_line195` |
| `froot_quartic_t4_plus_t_plus_1` | **Pipeline private** | `froot_quartic_t4_plus_t_plus_1` |
| `froot_linear_factor` | **Pipeline private** | `froot_linear_factor` |

### `fsolve.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_fsolve` | **Partial** | Newton numeric solve stub |
| `equation_to_expr` | **Pipeline private** | `equation_to_expr` |
| `eval_at` | **Pipeline private** | `eval_at` |
| `expr_to_f64` | **Pipeline private** | `expr_to_f64` |
| `newton` | **Pipeline private** | `newton` |
| `fsolve_x_squared_minus_two` | **Pipeline private** | `fsolve_x_squared_minus_two` |

### `plugin.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_solve` | **Stable (bounded)** | solve via poly roots, rootof, or linsolve |
| `eval_linsolve` | **Stable** | `Poly::eval_linsolve` |
| `eval_fsolve` | **Partial** | Newton numeric solve stub |
| `eval_sturm` | **Stable** | Sturm sequence for univariate poly |
| `eval_sturmab` | **Stable** | root count in (a,b) via Sturm |
| `eval_realroot` | **Stable (bounded)** | real roots via Sturm isolation |
| `eval_froot` | **Stable (bounded)** | rational roots of univariate poly |
| `install_solve` | **Stable** | register DefaultSolvePlugin |
| `xcas_default` | **Stable** | Context with simplify plugin |
| `eval_solve_via_plugin` | **Pipeline private** | `eval_solve_via_plugin` |

### `realroot.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_realroot` | **Stable (bounded)** | rational roots + filtered exact algebraic (deg≤4 via S0) |
| `real_roots_with_multiplicity` | **Pipeline private** | merge rational + exact real algebraic roots |
| `rational_real_roots` | **Pipeline private** | `rational_real_roots` |
| `real_algebraic_roots` | **Pipeline private** | S6: exact roots filtered by Sturm count / casus cubic |
| `monic_cubic_discriminant` | **Pipeline private** | `monic_cubic_discriminant` |
| `merge_rational` | **Pipeline private** | `merge_rational` |
| `merge_root` | **Pipeline private** | `merge_root` |
| `roots_eq` | **Pipeline private** | `roots_eq` |
| `format_root_key` | **Pipeline private** | `format_root_key` |
| `expr_is_real` | **Pipeline private** | `expr_is_real` |
| `expr_is_zero` | **Pipeline private** | `expr_is_zero` |
| `realroot_x_squared_minus_two` | **Pipeline private** | `realroot_x_squared_minus_two` |
| `realroot_x_fourth_minus_one` | **Pipeline private** | `realroot_x_fourth_minus_one` |
| `realroot_x_cubed_minus_x_minus_1` | **Pipeline private** | `realroot_x_cubed_minus_x_minus_1` |

### `rootof.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `quadratic_rootof_roots` | **Stable (bounded)** | two rootof branches for quadratic |
| `biquadratic_rootof_roots` | **Partial** | biquadratic rootof; general quartic NotImplemented |
| `rootof_expr` | **Pipeline private** | `rootof_expr` |
| `x` | **Pipeline private** | `x` |
| `t_sq_minus` | **Pipeline private** | `t_sq_minus` |
| `quadratic_rootof_has_two_branches` | **Pipeline private** | `quadratic_rootof_has_two_branches` |
| `solve_t_squared_minus_two_uses_rootof` | **Pipeline private** | `solve_t_squared_minus_two_uses_rootof` |
| `solve_t_fourth_minus_two_uses_rootof` | **Pipeline private** | `solve_t_fourth_minus_two_uses_rootof` |
| `biquadratic_t_fourth_minus_two_has_four_roots` | **Pipeline private** | `biquadratic_t_fourth_minus_two_has_four_roots` |
| `rational_quadratic_still_uses_roots` | **Pipeline private** | `rational_quadratic_still_uses_roots` |

### `solve.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_solve` | **Stable (bounded)** | solve via poly roots, rootof, or linsolve |
| `equation_to_poly` | **Pipeline private** | `equation_to_poly` |
| `try_transcendental_solve` | **Pipeline private** | optional fallback `try_transcendental_solve` |
| `is_zero` | **Stable** | Poly is zero |
| `eval_const_expr` | **Pipeline private** | `eval_const_expr` |
| `solve_quadratic_double_root` | **Pipeline private** | `solve_quadratic_double_root` |
| `solve_linear_system_delegates_to_linsolve` | **Pipeline private** | `solve_linear_system_delegates_to_linsolve` |
| `solve_biquadratic_t4_minus_2` | **Pipeline private** | `solve_biquadratic_t4_minus_2` |
| `solve_quartic_t4_plus_t_plus_1` | **Pipeline private** | `solve_quartic_t4_plus_t_plus_1` |

### `solve_poly.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `solve_univariate_over_q` | **Pipeline private** | sqff × factor → per-factor roots |
| `solve_irreducible_factor` | **Pipeline private** | deg≤4 algext; deg≥5 single rootof branch |
| `irreducible_rootof_branch` | **Pipeline private** | S0 deg≥5 rootof (not biquadratic fallback) |
| `factor_poly_for_solve` | **Pipeline private** | sqff then factor_into; factor failure → irreducible piece |
| `rootof_expr` | **Pipeline private** | rootof([num], minpoly) |
| `dedup_expr_roots` | **Pipeline private** | collapse repeated roots (solve lists unique roots) |
| `poly_x5_minus_x_plus_1` | **Pipeline private** | `poly_x5_minus_x_plus_1` |
| `solve_irreducible_deg5_one_branch` | **Pipeline private** | `solve_irreducible_deg5_one_branch` |
| `solve_x2_plus_1_has_two_roots` | **Pipeline private** | `solve_x2_plus_1_has_two_roots` |
| `solve_x3_minus_x_plus_1_has_three_roots` | **Pipeline private** | `solve_x3_minus_x_plus_1_has_three_roots` |
| `solve_reducible_deg5_product` | **Pipeline private** | `solve_reducible_deg5_product` |
| `solve_reducible_deg5_product_zeros_poly` | **Pipeline private** | `solve_reducible_deg5_product_zeros_poly` |
| `solve_x4_minus_1_factor_descent` | **Pipeline private** | `solve_x4_minus_1_factor_descent` |
| `eval_solve_irreducible_deg5` | **Pipeline private** | `eval_solve_irreducible_deg5` |

### `sturm.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_sturm` | **Stable** | Sturm sequence for univariate poly |
| `eval_sturmab` | **Stable** | root count in (a,b) via Sturm |
| `poly_and_var` | **Pipeline private** | `poly_and_var` |
| `eval_to_rational` | **Pipeline private** | `eval_to_rational` |
| `sturm_x_cubed_plus_one_squared` | **Pipeline private** | `sturm_x_cubed_plus_one_squared` |
| `sturm_x_cubed_plus_one` | **Pipeline private** | `sturm_x_cubed_plus_one` |
| `sturmab_counts_root_in_interval` | **Pipeline private** | `sturmab_counts_root_in_interval` |

### `test_verify.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `list_items` | **Pipeline private** | `list_items` |
| `eval_at` | **Pipeline private** | `eval_at` |
| `assert_roots_zero_poly` | **Pipeline private** | `assert_roots_zero_poly` |
| `assert_is_algext_or_rootof` | **Pipeline private** | `assert_is_algext_or_rootof` |
| `assert_equation_solutions` | **Pipeline private** | `assert_equation_solutions` |
| `assert_list_has_equiv` | **Pipeline private** | `assert_list_has_equiv` |
| `froot_pairs` | **Pipeline private** | `froot_pairs` |
| `assert_froot_has_root` | **Pipeline private** | `assert_froot_has_root` |
| `realroot_entries` | **Pipeline private** | `realroot_entries` |
| `assert_realroot_has` | **Pipeline private** | `assert_realroot_has` |
