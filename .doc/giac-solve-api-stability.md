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
| **表示** | 低次：`Poly::roots` → `poly_to_expr`；二次/双二次：`rootof`；否则 Err |
| **前提** | 方程经 `eval` 后可 `expr_to_poly`（ℚ 多项式） |
| **禁止** | 对含超越式未走 MRV/级数路径而强行 poly |

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
| `froot` | 有理根搜索链 | 依赖 `giac-poly` |
| `rootof` | `rootof_expr` | Expr 构造 |

---

## 5. 演化 issue

| 缺口 | 说明 |
|------|------|
| 一般四次 | 应收敛到 `poly_algext_roots`（[GIAC-poly-p3-6](issues/GIAC-poly-p3-6-quartic-roots-gaps.md)） |
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
| `rational_num_den` | **Pipeline private** | Expr leaf to (num,den) Poly |
| `append_factor_roots` | **Pipeline private** | `append_factor_roots` |
| `solve_factor_roots` | **Pipeline private** | `solve_factor_roots` |
| `froots_flanex_line195` | **Pipeline private** | `froots_flanex_line195` |
| `froot_linear_factor` | **Pipeline private** | `froot_linear_factor` |

### `fsolve.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_fsolve` | **Partial** | Newton numeric solve stub |
| `equation_to_expr` | **Pipeline private** | `equation_to_expr` |
| `ident_from_expr` | **Pipeline private** | `ident_from_expr` |
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
| `eval_realroot` | **Stable (bounded)** | real roots via Sturm isolation |
| `algebraic_real_roots` | **Pipeline private** | `algebraic_real_roots` |
| `rational_real_roots` | **Pipeline private** | `rational_real_roots` |
| `realroot_x_squared_minus_two` | **Pipeline private** | `realroot_x_squared_minus_two` |
| `realroot_x_fourth_minus_one` | **Pipeline private** | `realroot_x_fourth_minus_one` |

### `rootof.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `quadratic_rootof_roots` | **Stable (bounded)** | two rootof branches for quadratic |
| `biquadratic_rootof_roots` | **Partial** | biquadratic rootof; general quartic NotImplemented |
| `rootof_expr` | **Pipeline private** | `rootof_expr` |
| `x` | **Pipeline private** | `x` |
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
| `ident_from_expr` | **Pipeline private** | `ident_from_expr` |
| `try_transcendental_solve` | **Pipeline private** | optional fallback `try_transcendental_solve` |
| `is_zero` | **Stable** | Poly is zero |
| `is_sin_of_var` | **Pipeline private** | `is_sin_of_var` |
| `solve_quadratic_double_root` | **Pipeline private** | `solve_quadratic_double_root` |
| `solve_linear_system_delegates_to_linsolve` | **Pipeline private** | `solve_linear_system_delegates_to_linsolve` |

### `sturm.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_sturm` | **Stable** | Sturm sequence for univariate poly |
| `eval_sturmab` | **Stable** | root count in (a,b) via Sturm |
| `poly_and_var` | **Pipeline private** | `poly_and_var` |
| `ident_from_expr` | **Pipeline private** | `ident_from_expr` |
| `eval_to_rational` | **Pipeline private** | `eval_to_rational` |
| `sturm_x_cubed_plus_one_squared` | **Pipeline private** | `sturm_x_cubed_plus_one_squared` |
| `sturm_x_cubed_plus_one` | **Pipeline private** | `sturm_x_cubed_plus_one` |
| `sturmab_counts_root_in_interval` | **Pipeline private** | `sturmab_counts_root_in_interval` |
