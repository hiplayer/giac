# Expr 表示层并列架构 — upstream `gen` 全表

**状态:** normative（架构说明）  
**快照:** 2026-06-24（§3 `_GROB` 勘误、§6.6 Groebner 与 `_POLY` 对齐）  
**上游基线:** `giac/giac-2.0.0/src/{dispatch.h,gen.h,modpoly.h}`  
**相关:** [algorithm-expr-api.md](algorithm-expr-api.md) §1、[expr-poly-conversion.md](expr-poly-conversion.md)、[rust-migration-plan.md](rust-migration-plan.md) §3.3–3.4、[GIAC-algext-adoption.md](issues/GIAC-algext-adoption.md)、[GIAC-poly-flat-field-division-layering.md](issues/GIAC-poly-flat-field-division-layering.md)、[GIAC-poly-nested-ring-types.md](issues/GIAC-poly-nested-ring-types.md)、[giac-poly-api-stability.md](giac-poly-api-stability.md)

---

## 1. 问题陈述

upstream giac 用 **`gen`** 承载一切计算对象：每个值带 **`type` tag**（有时还有 **`subtype`**），多种表示 **并列**，**不是** OO 继承。

giac-rs 刻意 **拆分**：

- **用户层** — `Expr`（符号 AST，`Arc` 共享）
- **算法快路径** — `Poly` / `PolyMod` / `SparseSeries` / `AlgExt` 等 **不进** `Expr` 枚举（对标 `_POLY` / `_SPOL1` 等不默认暴露在用户树里）

本文列出 **upstream 存在的全部主要表示 tag**，映射到 giac-rs，并说明 **何时可转换、不可互换**。

---

## 2. 总览：Expr 在上，专用表示并列在下

```text
                         用户 / eval / solve / display
                                    │
                                    ▼
                         ┌─────────────────────┐
                         │   Expr（符号 AST）   │
                         └──────────┬──────────┘
                                    │ 窄转换（契约门控）
     ┌──────────┬──────────┬────────┼────────┬──────────┬──────────┐
     ▼          ▼          ▼        ▼        ▼          ▼          ▼
  Poly<ℚ>  Poly<AlgExt> PolyMod  SparseSeries  AlgExt   Matrix/List  evalf
  _POLY    _POLY+_EXT   modpoly  _SPOL1        _EXT     _VECT        _DOUBLE_
     │          │          │        │            │         │
     └─nested   └─flat     └─fpx/   └─limit/     └─标量    └─linalg/
       FAC        K[x]      Hensel    MRV          塔        groebner
```

**原则：**

1. 默认 **留在 Expr**；仅当算法 **明确要求** 某环/某局部结构时才投影下去。
2. 含 `sin`/`cos`/`exp` 的式子 **不得** 为认形状而强行 `expr_to_poly`。
3. **并列表示不可随意互换**（见 §5）；`series` 命令是 Func，`SparseSeries` 是结果载体。

---

## 3. Upstream `gen` 主类型全表（`dispatch.h`）

`enum gen_unary_types`（giac-2.0.0）— 每个 `gen` 值的 **顶层 tag**：

| Tag | C++ 载体 | 数学 / 用途 | giac-rs 对标 | 状态 |
|-----|----------|-------------|--------------|------|
| **`_INT_`** | 内联 `int val` | 小整数 | `Expr::Int` | ✅ |
| **`_DOUBLE_`** | 内联 `double` | 机器浮点 | `evalf` 结果 / 无独立 Expr 变体 | 🟡 数值层 |
| **`_FLOAT_`** | 内联扩展浮点 | 高精度浮点（平台相关） | 未单独建模 | ☐ |
| **`_ZINT`** | `mpz_t` | 大整数 | `Expr::Int`（`BigInt`） | ✅ |
| **`_REAL`** | `mpf_t` | 任意精度实数 | `evalf` | 🟡 |
| **`_FRAC`** | `fraction<gen>` | 有理数 \(a/b\)（分子分母可为 `gen`） | `Expr::Rat` / `Expr::Frac` | ✅ |
| **`_CPLX`** | `gen[2]` | 复数 \(re + im\)（分量是 `gen`） | `Expr::Complex` | ✅ |
| **`_IDNT`** | `identificateur` | 符号变元 `x` | `Expr::Symbol` | ✅ |
| **`_SYMB`** | `symbolic` | `sin(x)`、`exp(x)`、`f(args)`… | `Expr::Func` / `Add`/`Mul`/`Pow` | ✅ |
| **`_FUNC`** | `unary_function_ptr` | 未应用的内置函数名 | `FuncKind`（元数据） | ✅ |
| **`_POLY`** | `polynome`（稀疏张量） | 多元稀疏多项式 | `Poly<C>` | ✅ |
| **`_EXT`** | `gen[2]`（coords + minpoly） | 代数扩域元 / `rootof` | `Expr::AlgExt` / `AlgExtC` | ✅ |
| **`_ROOT`** | `real_complex_rootof` | 实/复根隔离表示（旧/特化） | 多归 `AlgExt` / `rootof` | 🟡 |
| **`_SPOL1`** | `sparse_poly1` | 一元稀疏 **截断** 级数 | `SparseSeries`（limit） | 🟡 |
| **`_MOD`** | `gen[2]` | \(a \bmod m\)（\(a,m\) 均为 `gen`） | `Expr::Mod` | ✅ |
| **`_VECT`** | `vecteur` | 列表/矩阵/程序…（见 §3.2） | `List` / `Matrix` / `Seq`… | 🟡 |
| **`_MAP`** | `gen_map` | 稀疏矩阵、字典 | 部分 `Matrix` subtype | ☐ |
| **`_STRNG`** | `string` | 字符串 | `Expr::Str` | ✅ |
| **`_USER`** | `gen_user` | 用户自定义对象 | 未移植 | ☐ |
| **`_GROB`** | `grob`（`draw`/`handle` 回调 + `void*`） | **XCAS 图形对象**（`DIMGROB`、`SUBGROB`…） | 不移植（与 `_POINTER_` 同类） | ☐ |
| **`_EQW`** | `eqwdata` | **排版**树（非计算语义） | 未移植 | ☐ |
| **`_POINTER_`** | `void*` | 外部裸指针（图形上下文等） | 不移植 | ☐ |

说明：

- **`_SYMB`** 与 **`_FUNC`**：`sin`  applied 是 `_SYMB`；单独 `sin` 作为函数对象可以是 `_FUNC`。
- **`_ROOT`** 与 **`_EXT`**：upstream 为两个 tag，数学上同属代数数；giac-rs **合并**为 `AlgExt` / `AlgExtC`（`root_index` 承担选枝），无并列 `_ROOT` 变体。
- **`_GROB` ≠ Groebner**：名字易混；Groebner 基全程用 **`_POLY`**（`vectpoly = vector<polynome>`），见 §3.1、§6.6。
- **`_EQW`**：仅 display / `latex` / `mathml`；**不参与** gcd/factor/integrate。

---

## 3.1 算法专用载体（**不是** 独立 `gen` tag）

这些类型在 upstream 算法内部使用，常由 `_POLY` / `_VECT` **转换**而来：

| 载体 | 定义 | 数学环 | giac-rs | 典型算法 |
|------|------|--------|---------|----------|
| **`modpoly`** | `typedef vecteur modpoly` | \(\mathbb{F}_p[x]\) 或 \(\mathbb{Z}/p\mathbb{Z}[x]\)（稠密） | **`PolyMod`** | `modfactor`、Hensel、`chinrem` |
| **`dense_POLY1`** | `_VECT` + subtype `_POLY1__VECT` | 一元稠密系数向量 | `giac-poly::dense`（D1–D4） | `horner`、模 gcd |
| **`polynome`** | `_POLY` | 稀疏多元多项式 | `Poly<C>` | gcd、factor、Groebner |
| **`vectpoly`** | `vector<polynome>` | Gröbner 基 / 理想生成元列表 | `Vec<Poly>` | `gbasis`、`greduce` |
| **`sparse_poly1`** | `_SPOL1` | 局部 Lauent 级数 | `SparseSeries` | `limit`、`series` |

`modularize(polynome) → modpoly`（`modpoly.h`）；`modp(Poly) → PolyMod`（giac-rs）。

---

## 3.2 `_VECT` 子类型（`comp_subtypes`）

同一 `_VECT` tag，**`subtype` 区分语义**（节选与 CAS 相关）：

| Subtype | 名称 | 用途 | giac-rs |
|---------|------|------|---------|
| `_SEQ__VECT` | 序列 | `(a,b,c)`、程序参数 | `Expr::Seq` |
| `_SET__VECT` | 集合 | `{1,2,3}` | 部分为 `List` |
| `_LIST__VECT` | 列表 | `[1,2,3]` | `Expr::List` |
| `_MATRIX__VECT` | 矩阵 | `[[1,2],[3,4]]` | `Expr::Matrix` / `GiacMatrix` |
| `_POLY1__VECT` | 一元稠密 | `poly1[1,0,-2]` | `Poly` / `dense::poly1` |
| `_VECTOR__VECT` | 几何向量 | 图形 | ☐ |
| `_ASSUME__VECT` | 假设 | `assume` | `Context` 部分 |
| `_PRG__VECT` | 程序 | `giac-prog` | ☐ |
| `_INTERVAL__VECT` | 区间 | 实分析 | ☐ |

**纪律：** 矩阵运算在 `giac-linalg`；符号矩阵元素仍是 `Expr`，数值快路径可用 `f64`（不进 `Expr`）。

---

## 4. giac-rs `Expr` 枚举与 upstream 对照

实现：`giac-rs/crates/giac-core/src/expr.rs`。

| `Expr` 变体 | upstream 来源 | 备注 |
|-------------|---------------|------|
| `Int` / `Rat` | `_INT_` / `_ZINT` / `_FRAC` | |
| `Frac` | `_FRAC` 或 `a/b` 未合并 | |
| `Complex` | `_CPLX` | `im` 可为 `Symbol(i)` 或 `AlgExt` |
| `Symbol` | `_IDNT` | |
| `Add` / `Mul` / `Pow` | `_SYMB`（`at_plus` / `at_prod` / `at_pow`） | 也可保留为 `Func` |
| `Func` | `_SYMB` + `_FUNC` | 超越函数、builtins |
| `AlgExt` / `AlgExtC` | `_EXT`（`_ROOT` 合并，见 §3 说明） | |
| `Mod` | `_MOD` | `a mod p` |
| `List` / `Matrix` / `Seq` | `_VECT` + subtype | |
| `Relation` | `_SYMB` `at_equal` 等 | |
| `Str` | `_STRNG` | |
| `Undefined` | 错误/未绑定 | |

**刻意不在 `Expr` 内：** `Poly`、`PolyMod`、`SparseSeries`（[rust-migration-plan.md](rust-migration-plan.md) §3.4）。

---

## 5. 并列表示：何时可转换、不可互换

| 从 → 到 | 是否可换 | 条件 / API |
|---------|----------|------------|
| Expr → Poly\<ℚ\> | 门控 | `expr_to_poly`；无超越、无 AlgExt 系数 |
| Expr → Poly\<AlgExt\> | 门控 | `poly_alg_from_expr` |
| Poly → Expr | ✅ | `poly_to_expr` / `algext_poly_to_expr` |
| Expr → PolyMod | 门控 | `modp`；系数整数 |
| PolyMod → Expr | ✅ | `poly_mod_to_expr` → `Expr::Mod` |
| Expr → SparseSeries | 局部 | `series_at_zero` 等；**有阶、有中心** |
| SparseSeries → Expr | 近似 | `to_expr`；**非**与 `sin` 语义相等 |
| Func(Sin) ↔ SparseSeries | ❌ | 仅 `series(sin(x),…)` **近似**；不可默认替换 |
| Poly ↔ PolyMod | 算法内 | `modp` / CRT lift；Hensel |
| Poly ↔ modpoly | upstream 内 | `modularize`；Rust 用 `PolyMod` |
| `_SPOL1` ↔ `_POLY` | 算法内 | `sparse_poly12gen` 等 |

---

## 6. 按数学对象的管线分区

### 6.1 多项式环（`_POLY` / `modpoly`）

| 环 | 表示 | Crate | 入口 |
|----|------|-------|------|
| ℚ[x] 稀疏 | `Poly` | `giac-poly` | `expr_to_poly` |
| K[x] flat | `Poly<AlgExt>` | `giac-core` | `poly_alg_from_expr` |
| ℚ[others][main] nested | `Poly` + `UnivariateIn` | `giac-poly::factor` | FAC 内部 |
| ℤ/pℤ[x] | `PolyMod` | `giac-poly` | `modp` |
| \(\mathbb{F}_p[x]\) 分解 | `PolyMod` | `factor/fpx` | `factor_poly_mod` |

Flat vs Nested 除法语义见 [GIAC-poly-flat-field-division-layering.md](issues/GIAC-poly-flat-field-division-layering.md)、[GIAC-poly-nested-ring-types.md](issues/GIAC-poly-nested-ring-types.md)。

### 6.2 超越函数（`_SYMB`）

`sin`/`cos`/`exp`/`ln` → **`Expr::Func`**；**不**进 `Poly` 系数。

专用重写（仍 Expr）：`texpand`、`trig2exp`（upstream `subst.cc`）；giac-rs 见 §8。

### 6.3 级数（`_SPOL1`）

- **命令** `series(...)` → `Expr::Func(Series, …)`（与 `sin` 并列的 **算子**）。
- **结果** `SparseSeries` → 与 `Func` **不同族**；局部截断，不可与 `sin(x)` 全局互换。

### 6.4 代数数（`_EXT`）

`rootof` / `AlgExt`；可进 `Poly` 系数（路径 B）。**不是**超越数。

### 6.5 模算术（`_MOD` + `modpoly`）

- `Expr::Mod(a, m)` — 符号层 \(a \bmod m\)。
- `PolyMod` — \(\mathbb{Z}/p\mathbb{Z}[x]\)；`factor mod p`、Hensel、模 Groebner。

upstream：`modpoly.cc`、`modfactor.cc`；`ezgcd.cc` Hensel 链。

### 6.6 Groebner（`_POLY` / `vectpoly`，**非** `_GROB`）

upstream `gbasis` / `greduce` 在 **`vectpoly`**（`vector<polynome>`，每个元素 tag 为 `_POLY`）上运算；**没有** 独立的 Gröbner tag。

giac-rs **`giac-groebner`** 对标同一模型：

| 步骤 | upstream | giac-rs | 状态 |
|------|----------|---------|------|
| 理想约化 | `greduce(p, G, order)` | `greduce(&Poly, &[Poly], vars)` | ✅ |
| 模约化 | 模 `p` 版 `reduce` | `greduce_mod` | ✅ |
| 基构造 | `gbasis`（F4 / 模重构） | 规划 | ☐ |

详见 [giac-groebner-api-stability.md](giac-groebner-api-stability.md)。

### 6.7 XCAS 图形（`_GROB` / `_POINTER_`）

`gen.h` 中 `struct grob` 为 **graphic object**（绘图回调 + 不透明数据），用于 XCAS 屏幕 API（`rpn.cc`：`DIMGROB`、`SUBGROB`、`GROBW`…），**不参与** CAS 代数。

giac-rs 与 [rust-migration-plan.md](rust-migration-plan.md) §3.1 一致：**不移植**。若将来做图形层，归 **`giac-geo` / emit**，不进 `Expr` 或 `Poly`。

### 6.8 数值（`_DOUBLE_` / `_REAL`）

`evalf`：Expr → 浮点；与符号层分离。

### 6.9 排版（`_EQW`）

`eqwdata` → `latex` / `mathml`；**不**参与代数运算。giac-rs 未移植。

---

## 7. 算法选路（简图）

```text
输入 Expr
    │
    ├─ factor/gcd/sturm/partfrac → Poly 或 Poly<AlgExt>（§6.1）
    ├─ factor mod p → PolyMod（§6.5）
    ├─ greduce / gbasis → Poly 或 PolyMod（§6.6）
    ├─ normal/ratnormal/texpand → 全程 Expr（§6.2）
    ├─ limit → SparseSeries / MRV（§6.3）
    ├─ integrate/diff → Expr（+ 局部 Poly 系数块）
    ├─ solve 超越 → Expr 重写链（§8）
    └─ evalf → 数值（§6.8）
```

---

## 8. 三角 / 指数重写（Expr 层）

| | **`texpand`** | **`trig2exp`** |
|---|---------------|----------------|
| 目的 | 角展开（恒等式） | Euler：\(\sin,\cos \to e^{ix}\) |
| 结果形态 | 仍为 `sin`/`cos` | `exp` + `i` |
| `trig2exp` 门控 | — | **`angle_radian`** |

触发：`trig2exp`/`tsimplify` 用户命令；`solve::rationalize`；Risch 前 `lin(trig2exp(...))` 等。  
**不**在 `factor`/`e2r` 默认路径。

giac-rs 缺口：[GIAC-simplify-poly-upstream-gaps.md](issues/GIAC-simplify-poly-upstream-gaps.md)。

---

## 9. 工程纪律

1. 新算法先查 §3 表：**upstream 用哪个 tag？** Rust 用哪条 API？
2. 含超越子树 **禁止** 静默 `expr_to_poly`。
3. Poly：**Flat / Nested** 分清（§6.1）。
4. **Mod**：`Expr::Mod` ≠ `PolyMod`；后者是多项式环投影。
5. 转换走稳定 API：[expr-poly-conversion.md](expr-poly-conversion.md)、[giac-poly-api-stability.md](giac-poly-api-stability.md)。

---

## 10. 验证与文档索引

```bash
cd giac-rs
cargo test -p giac-core algebra::poly expr_to_poly poly_alg_from_expr
cargo test -p giac-poly modular modp
cargo test -p giac-calculus limit_engine sparse_series
cargo test -p giac-simplify trig
```

| 主题 | 文档 |
|------|------|
| Expr ↔ Poly | [expr-poly-conversion.md](expr-poly-conversion.md) |
| 通用规范 | [algorithm-expr-api.md](algorithm-expr-api.md) |
| 迁移总图 | [rust-migration-plan.md](rust-migration-plan.md) §3 |
| AlgExt | [GIAC-algext-adoption.md](issues/GIAC-algext-adoption.md) |
| PolyMod / modp | [giac-poly-api-stability.md](giac-poly-api-stability.md) |
| Flat / Nested | [GIAC-poly-flat-field-division-layering.md](issues/GIAC-poly-flat-field-division-layering.md)、[GIAC-poly-nested-ring-types.md](issues/GIAC-poly-nested-ring-types.md) |
| Limit / 级数 | [limit-engine-expr-api.md](limit-engine-expr-api.md) |
| Groebner | [giac-groebner-api-stability.md](giac-groebner-api-stability.md) |
| 模块分工 | [module-division.md](module-division.md) |
