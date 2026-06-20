# Expr ↔ Poly 转换边界（B-01 / P0-C）

**状态:** normative（P0-C 已落地 giac-core）  
**实现:** `giac-rs/crates/giac-core/src/algebra/{poly.rs,poly_conv.rs}`；P1 全貌见 [giac-poly-p1-representation.md](giac-poly-p1-representation.md)  
**相关:** [GIAC-algext-adoption.md](issues/GIAC-algext-adoption.md) B-01、[GIAC-poly-algext-backlog.md](issues/GIAC-poly-algext-backlog.md) P0-C/P1-4、[algorithm-expr-api.md](algorithm-expr-api.md) §1

---

## 1. 两条路径

| 路径 | API | 目标类型 | 系数环 | 状态 |
|------|-----|----------|--------|------|
| **A — 有理多项式** | `expr_to_poly` | `giac_poly::Poly` | **ℚ** | ✅ 稳定 |
| **B — 代数系数多项式** | `poly_alg_from_expr` | `Poly<AlgExtCPolyCoeff>` | **K** / `AlgExtC` | ✅ 稳定 |

**原则：** 含 `rootof` / `AlgExt` / `AlgExtC` 的式子 **不得** 静默落入路径 A；须显式 `TypeError` 或走路径 B。

---

## 2. 路径 A：`expr_to_poly`（ℚ）

### 输入契约

- 式子为 **多项式形状**：`+`、`*`、整数/有理数幂、`Symbol` 变元。
- **系数 ∈ ℚ**（`Int` / `Rat`）。
- **禁止** 系数位置出现：
  - `Expr::AlgExt`
  - `Expr::AlgExtC`
  - `Func(RootOf, …)`（含嵌套于 `Add`/`Mul` 内）

### 输出

- `Ok(Poly)`：`giac-poly::Poly`，`terms: BTreeMap<Monomial, Ratio<BigInt>>`。
- `Err(EvalError::TypeError(msg))`：见 §4 错误消息表。

### 调用方

| Crate | 典型用途 |
|-------|----------|
| `giac-core` | `eval_poly`（gcd/resultant/factor/partfrac 等 builtin） |
| `giac-simplify` | `ratnormal`、`factor`（ℚ 因子）、`expand` |
| `giac-solve` | 多项式方程、`sturm`、`froot` 有理根 |
| `giac-calculus` | `partfrac_integrate`（ℚ 分式）、`limit_engine` 有理分支 |

**含代数系数的式子**（如 `x + rootof([1,0], poly1[1,0,-2])`）在此路径 **必须失败**；调用方应改走路径 B 或保持 `Expr` 层 `AlgExt` 运算。

---

## 3. 路径 B：`poly_alg_from_expr`（K，规划中）

### 目标（P1-4）

```text
poly_alg_from_expr(expr) → Poly<AlgExtC>
```

- 变元仍为 `x,y,…`（`Monomial` / `Var` 不变）。
- 系数为 `AlgExtC`（`giac-core` 塔顶域 + 可选虚部）。
- 逆变换：`algext_poly_to_expr` → `Expr::AlgExtC` / `rootof` 显示。

### 当前行为（P1-4 ✅）

```rust
poly_alg_from_expr(expr) → Result<PolyAlgExt, EvalError>   // Poly<AlgExtCPolyCoeff>
algext_poly_to_expr(poly) → Result<ExprArc, EvalError>
```

若式子 **不含** 代数系数 → `TypeError(ERR_POLY_ALG_NO_ALG_COEFF)`，提示改用 `expr_to_poly`。

**构造 API（非 ℚ）：** `Poly::<C>::ring_zero` / `ring_one` / `ring_constant` / `ring_var`；ℚ 路径仍用 `Poly::zero` / `one` / `constant` / `var`。

### 解锁算法（P2+）

| 算法 | 需要路径 B |
|------|------------|
| `partfrac(1/(x²-2), x)` disc>0 分裂 | ✅ |
| `factor(x²-2)` over ℚ 不可约 | ✅ |
| `Poly<AlgExtC>::roots` | ✅ |

---

## 4. 稳定 API 与错误消息

| 函数 | Tier | 说明 |
|------|------|------|
| `expr_to_poly` | **Stable** | ℚ 多项式；拒绝代数系数 |
| `poly_to_expr` | **Stable** | `Poly` → `Expr`（系数 ∈ ℚ） |
| `expr_contains_alg_coeff` | **Stable** | 谓词：是否含 `AlgExt`/`AlgExtC`/`rootof` |
| `poly_alg_from_expr` | **Stable** | `PolyAlgExt` 提升；含代数系数 |
| `algext_poly_to_expr` | **Stable** | `PolyAlgExt` → `Expr` |

| 常量 | `EvalError` | 场景 |
|------|-------------|------|
| `ERR_ALG_EXT_COEFF` | `TypeError` | `Expr::AlgExt` 作系数 |
| `ERR_ALG_EXT_C_COEFF` | `TypeError` | `Expr::AlgExtC` 作系数 |
| `ERR_ROOTOF_COEFF` | `TypeError` | `rootof(...)` 子树 |
| `ERR_POLY_ALG_NO_ALG_COEFF` | `TypeError` | `poly_alg_from_expr` 但式子纯 ℚ |

---

## 5. 调用方决策

```text
expr_contains_alg_coeff(e)?
  ├─ false → expr_to_poly(e)           // 现有 ℚ 管线
  └─ true  → poly_alg_from_expr(e)     // Poly<AlgExtCPolyCoeff>
             或保持 Expr 层 AlgExt 运算（eval / fold_algext_*）
```

**禁止：**

- 在 `giac-simplify` / `partfrac.rs` 对 `Poly<ℚ>` 手写 `rootof` 分裂补丁。
- 对含 `AlgExt` 的式子调用 `expr_to_poly` 并忽略 `TypeError`。

---

## 6. 与上游 giac 的对应

| upstream | giac-rs |
|----------|---------|
| `e2r`（有理式） | `expr_to_poly` |
| `e2r` + `_EXT` 系数 | `poly_alg_from_expr`（P1） |
| `r2e` | `poly_to_expr` / `algext_poly_to_expr` |

---

## 7. 验证

```bash
cd giac-rs
cargo test -p giac-core expr_to_poly poly_conv poly_alg
```
