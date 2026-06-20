# giac-ode API 稳定性分层

**规范来源:** [algorithm-expr-api.md](algorithm-expr-api.md)  
**算法缺口:** [issues/GIAC-algorithm-gaps-open.md](issues/GIAC-algorithm-gaps-open.md)  
**代码:** `giac-rs/crates/giac-ode/src/*`

---

## 1. 源码注释格式

| 标记 | 用于 | 可见性 |
|------|------|--------|
| `/// **Stable** — …` | plugin 注册 | `pub` |
| `/// **Stable (bounded)** — …` | `desolve` 子集 | `pub` |
| `// **Pipeline private** — …` | ODE 解析 / 系数提取 | `fn` 私有 |

---

## 2. Crate 公开 API

| 函数 | 层级 | 模块 | 边界 |
|------|------|------|------|
| `eval_desolve` | **Stable (bounded)** | `desolve` | 常系数线性 ODE；`a2*y''+a1*y'+a0*y=f` |
| `install_ode` / `xcas_default` | **Stable** | `plugin` | |

---

## 3. I/O 契约 — `eval_desolve`

| 字段 | 说明 |
|------|------|
| **输入** | `desolve(equation, y(x))` |
| **输出** | `Relation(Eq, y(x), solution)` |
| **表示** | 全程 `ExprArc`；含积分常数 `c0`, `c1`, … |
| **边界** | 非常系数 / 非线性 → `TypeError` 或 `NotImplemented` |
| **禁止** | 用 `format_expr.contains("c0")` 作语义断言（见 tech-debt 3A） |

---

## 4. Pipeline private

| 函数 | 说明 |
|------|------|
| `parse_dep_fn`, `parse_linear_ode`, `solve_linear_ode` | 方程形状识别与求解 |
| `collect_derivatives`, `coeff_of_derivative` | 导数项收集 |

---

## 5. 维护

```bash
cd giac-rs
python3 scripts/annotate_api_tiers.py
python3 scripts/annotate_api_tiers.py --inventory
```

---

## Per-file function inventory (generated)

Regenerate: `python3 scripts/annotate_api_tiers.py --inventory`

### `desolve.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_desolve` | **Stable (bounded)** | linear constant-coefficient ODE subset |
| `parse_dep_fn` | **Pipeline private** | `parse_dep_fn` |
| `dep_fn_expr` | **Pipeline private** | `dep_fn_expr` |
| `parse_linear_ode` | **Pipeline private** | `parse_linear_ode` |
| `flatten_add` | **Pipeline private** | `flatten_add` |
| `add_expr` | **Pipeline private** | `add_expr` |
| `derivative_order` | **Pipeline private** | `derivative_order` |
| `is_dep` | **Pipeline private** | `is_dep` |
| `strip_derivative_factor` | **Pipeline private** | `strip_derivative_factor` |
| `solve_linear_ode` | **Pipeline private** | `solve_linear_ode` |
| `solve_first_order` | **Pipeline private** | `solve_first_order` |
| `is_one_expr` | **Pipeline private** | `is_one_expr` |
| `is_neg_var` | **Pipeline private** | `is_neg_var` |
| `is_var_expr` | **Pipeline private** | `is_var_expr` |
| `solve_second_order` | **Pipeline private** | `solve_second_order` |
| `expr_to_ratio` | **Pipeline private** | `expr_to_ratio` |
| `homogeneous_second_order` | **Pipeline private** | `homogeneous_second_order` |
| `particular_sin` | **Pipeline private** | `particular_sin` |
| `exp_rat_times_x` | **Pipeline private** | `exp_rat_times_x` |
| `trig_rat_times_x` | **Pipeline private** | `trig_rat_times_x` |
| `const_sym` | **Pipeline private** | `const_sym` |
| `ratio_sqrt` | **Pipeline private** | `ratio_sqrt` |
| `integer_sqrt` | **Pipeline private** | `integer_sqrt` |
| `is_zero_expr` | **Pipeline private** | `is_zero_expr` |
| `is_sin_of_var` | **Pipeline private** | `is_sin_of_var` |
| `desolve_harmonic` | **Pipeline private** | `desolve_harmonic` |

### `plugin.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `eval_desolve` | **Stable (bounded)** | linear constant-coefficient ODE subset |
| `install_ode` | **Stable** | register DefaultOdePlugin |
| `xcas_default` | **Stable** | Context with simplify plugin |
| `desolve_via_plugin` | **Pipeline private** | `desolve_via_plugin` |
