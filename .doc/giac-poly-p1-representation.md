# giac-poly P1 表示层（`Poly<C>` / AlgExt 桥接）

**状态:** normative（P1 已闭合）  
**实现:** `giac-rs/crates/giac-poly/src/{poly.rs,poly_coeff.rs,monomial.rs,nested.rs}`；`giac-core/src/algebra/{poly.rs,poly_alg_coeff.rs}`  
**相关:** [GIAC-poly-algext-backlog.md](issues/GIAC-poly-algext-backlog.md) §4、[expr-poly-conversion.md](expr-poly-conversion.md)、[GIAC-poly-nested-ring-types.md](issues/GIAC-poly-nested-ring-types.md)

---

## 1. 目标

在 **不破坏** 现有 `Poly`（ℚ）算法管线的前提下，引入：

1. 系数环抽象 `PolyCoeff`
2. 泛型稀疏多项式 `Poly<C>`
3. Expr 路径 B（`Poly<AlgExtC>`）
4. 单变量视图的泛型骨架（嵌套环算法仍 ℚ-only）
5. 变元/单项式与系数环解耦（P1-5）

---

## 2. 类型分层

```text
Monomial / Var              指数向量（与 C 无关）          P1-5 ✅
    ↓
PolyCoeff                   zero/one/add/mul/div           P1-1 ✅
    ├─ Ratio<BigInt>        giac-poly/poly_coeff.rs
    └─ AlgExtCPolyCoeff     giac-core/poly_alg_coeff.rs
    ↓
Poly<C: PolyCoeff>          BTreeMap<Monomial, C>          P1-2 ✅
    ├─ Poly (= Poly<ℚ>)     gcd/factor/partfrac 全管线
    └─ PolyAlgExt           giac-core 路径 B
    ↓
单变量视图（骨架）
    ├─ UnivariateIn<'a, C> / UnivariatePoly<C>   degree ✅；divides 仅 ℚ
    └─ FlatUni<C>                                  degree ✅；div_rem 仅 ℚ
```

---

## 3. P1 子项验收

| ID | 交付 | 位置 |
|----|------|------|
| **P1-1** | `PolyCoeff` trait + `Ratio<BigInt>` impl | `giac-poly/src/poly_coeff.rs` |
| **P1-2** | `Poly<C>` + `PolyQ`；`try_add/mul/…`；ℚ 用 `add/mul/…` | `giac-poly/src/poly.rs` |
| **P1-3** | `UnivariateIn<C>` / `UnivariatePoly<C>` / `FlatUni<C>` 骨架；**嵌套环整除仍 ℚ** | `giac-poly/src/nested.rs` |
| **P1-4** | `poly_alg_from_expr` / `algext_poly_to_expr` | `giac-core/src/algebra/poly.rs` |
| **P1-5** | `Var`/`Monomial` 不变；`Poly::degree_wrt` 与 C 无关 | `monomial.rs` + `poly.rs` |

---

## 4. API 速查

### 4.1 构造（系数环相关）

| 环 | 零/一/常数/变元 |
|----|----------------|
| **ℚ** | `Poly::zero()` / `one()` / `constant(r)` / `var(name)` |
| **任意 C** | `Poly::<C>::ring_zero()` / `ring_one()` / `ring_constant(c)` / `ring_var(name)` |

### 4.2 Expr 转换

| 路径 | API | 类型 |
|------|-----|------|
| A（ℚ） | `expr_to_poly` | `Poly` |
| B（K） | `poly_alg_from_expr` | `PolyAlgExt` = `Poly<AlgExtCPolyCoeff>` |
| B 逆 | `algext_poly_to_expr` | `Expr`（`AlgExt` / `rootof` / `AlgExtC`） |

门控：`expr_contains_alg_coeff` — 见 [expr-poly-conversion.md](expr-poly-conversion.md)。

### 4.3 单变量视图

| 类型 | 泛型 | ℚ 专用算法 |
|------|------|------------|
| `UnivariateIn<'a, C>` | `new`, `degree` | `coeff_at`, `divides`, `exact_quo_dividing`, `eval_aux` |
| `UnivariatePoly<C>` | `new`, `as_view` | （同上，经 `UnivariateIn<'a>`） |
| `FlatUni<C>` | `new`, `as_poly`, `degree` | `try_new`, `div_rem`, `divides`, `exact_quo` |

**纪律：** `Poly<AlgExtC>` 上 **禁止** 调用 `UnivariateIn::divides` / `FlatUni::div_rem`（无 impl）；P3-5 再扩展嵌套环 `K[others][main]`。

---

## 5. P1-5：单项式与变元

- [`Var`](giac-rs/crates/giac-poly/src/monomial.rs) = `Arc<str>`，排序/相等与系数无关。
- [`Monomial`](giac-rs/crates/giac-poly/src/monomial.rs) = `BTreeMap<Var, u64>`，全序由指数向量决定。
- `Poly<C>` 的 `terms` 键仍为 `Monomial`；**仅** 值类型 `Ratio<BigInt>` → `C` 变化。
- `MainVar` 标签主元变元，不编码系数环（嵌套环 issue 见 [GIAC-poly-nested-ring-types.md](issues/GIAC-poly-nested-ring-types.md)）。

---

## 6. 与上游 giac 的对应

| upstream (`gausspol`) | giac-rs P1 |
|-----------------------|------------|
| `_POLY` 系数 `gen` | `Poly<C>` 稀疏项 |
| `_EXT` 系数 | `PolyAlgExt` / `AlgExtCPolyCoeff` |
| `e2r` / `e2r`+`_EXT` | `expr_to_poly` / `poly_alg_from_expr` |
| `r2e` | `poly_to_expr` / `algext_poly_to_expr` |

P1 **不** 移植 `_EXT` 上的 gcd/factor/Hensel（→ P2/P3）。

---

## 7. 验证

```bash
cd giac-rs
cargo test -p giac-poly --lib
cargo test -p giac-core algebra::poly poly_alg_coeff
```

---

## 8. 下一步（P2）

- `PolyAlgExt::roots` 二次（P2-1）
- partfrac disc>0 分裂（P2-3）
- 嵌套环 `Poly<AlgExtC>[main]`（P3-5，长期）
