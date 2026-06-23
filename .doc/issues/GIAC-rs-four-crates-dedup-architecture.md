# GIAC-rs — 四 crate 去重架构确认（core / poly / calculus / solve）

**状态:** open（架构已确认，待竖切实施）  
**类型:** 架构决策 / 执行计划  
**来源:** `giac-core`、`giac-poly`、`giac-calculus`、`giac-solve` 重复代码审计（2026-06-23）  
**父项:** [GIAC-rs-crate-dedup-plan](GIAC-rs-crate-dedup-plan.md)（D0–D5 执行索引）  
**相关:** [GIAC-poly-algext-gcd-factor-priority](GIAC-poly-algext-gcd-factor-priority.md) T3-3、[expr-poly-conversion.md](../expr-poly-conversion.md)、[algorithm-expr-api.md](../algorithm-expr-api.md)  
**快照:** 2026-06-23

---

## 1. 问题陈述

`giac-core`、`giac-poly`、`giac-calculus`、`giac-solve` 在快速移植期形成 **平行算法路径** 与 **工具函数拷贝**，相似度 >50% 的块包括：

| 类别 | 典型位置 | 风险 |
|------|----------|------|
| 工具 fn 拷贝 | `ratio_to_expr`×6、`integer_sqrt`×5+、`rootof_expr`×2 | 修一处漏多处 |
| 二次多项式 | `quadratic_abc` / `quadratic_roots_formula` / `quadratic_rootof_roots` | 判别式逻辑三轨 |
| 双二次 / 偶次四次 | `biquadratic_roots` / `biquadratic_rootof_roots` / `try_factor_biquadratic` | 代换+分支重复 |
| giac-poly 内部 | `coeff_wrt`×3、`vars_in`×2、Yun sqff×3 | crate 内维护成本 |
| 求解编排 | `factor_poly_for_solve` ≈ `factor_univariate_pairs` | 薄包装重复 |

**已确认原则：**

1. **行为不变优先** — 合并不改数学语义；等价形走 `assert_equiv`。  
2. **单一权威落点** — 每层只有一个 owner（见 §3）。  
3. **进管线，不旁路** — 退役平行路径时扩展主路径，禁止第三个形状特例。  
4. **与 dedup-plan 对齐** — 本 issue 补 **四 crate 架构视图**；子项 ID **C\*** 映射到 [D0–D5](GIAC-rs-crate-dedup-plan.md) 已有竖切。

---

## 2. 分层架构（已确认）

```text
┌─────────────────────────────────────────────────────────────┐
│  giac-calculus / giac-solve  — Expr 编排（diff/limit/solve）   │
│  禁止：本地 ratio_to_expr、二次/四次求根、sqff×factor 骨架     │
└───────────────────────────┬─────────────────────────────────┘
                            │ expr_to_poly / poly_to_expr
                            │ poly_algext_roots_for_ctx
┌───────────────────────────▼─────────────────────────────────┐
│  giac-core  — Expr 枢纽 + AlgExt 塔 + deg≤4 代数求根          │
│  poly.rs（Expr↔Poly）· poly_roots.rs · field_session         │
└───────────────────────────┬─────────────────────────────────┘
                            │ Poly / PolyCoeff / dense 委托
┌───────────────────────────▼─────────────────────────────────┐
│  giac-poly  — 稀疏多项式代数 + ℚ 分解/有界根 + 泛型 dense     │
│  quadratic_abc · factor/* · univariate · subresultant        │
└─────────────────────────────────────────────────────────────┘

层 0 工具（giac-core::num_util）：ratio_to_expr、integer_nth_root
```

**明确不合并（保持分层）：**

| 项 | 理由 |
|----|------|
| `PolyCoeff` vs `AlgExtCPolyCoeff` | 不同系数环，有意分层 |
| `giac-core::field_arith` vs `giac-poly::dense` | 已委托 dense；field_arith 保留 ParentBlockRing 上下文 |
| `SparseSeries`（calculus） | limit 专用；giac-poly 无对应类型 |
| 各 crate `plugin.rs` | 组合模式，非重复 |

---

## 3. 重复块清单（相似度 >50%）

### 3.1 跨 crate 工具（~90–100%）

| 符号 | 权威落点 | 待删副本 |
|------|----------|----------|
| `ratio_to_expr` | `giac-core/algebra/poly.rs` | `field_arith::ratio_to_expr_arc`、`eval.rs` 私有版；calculus `limit_engine/{mod,asymptotic,sparse_series}.rs`、`partfrac_integrate.rs` |
| `integer_nth_root` / `integer_sqrt` | `giac-core/num_util`（新建或扩展） | `eval.rs`；`giac-poly/factor/sparse.rs`；`partfrac_integrate.rs` |
| `rootof_expr` | `giac-core/algebra/alg_ext.rs`（新建 `rootof_from_minpoly`） | `giac-solve/solve_poly.rs`、`rootof.rs` |
| `is_negative_rational` | `giac-core/algebra/poly_roots.rs`（唯一） | `field_session.rs` 副本 |

### 3.2 二次多项式（~60–80%）

| 位置 | 函数 | 输出 |
|------|------|------|
| `giac-poly/resultant.rs` | `quadratic_abc` → `quadratic_roots` | ℚ 有理根 `Poly` |
| `giac-poly/factor/univariate.rs` | `factor_quadratic` | 线性因子 |
| `giac-core/algebra/poly_roots.rs` | `quadratic_roots_formula` | `AlgExtCPolyCoeff` |
| `giac-solve/rootof.rs` | `quadratic_rootof_roots` | `rootof` Expr |

### 3.3 双二次 / 偶次四次（~55–70%）

| 位置 | 函数 | 目标 |
|------|------|------|
| `giac-core/poly_roots.rs` | `biquadratic_roots` | AlgExt 四根 |
| `giac-solve/rootof.rs` | `biquadratic_rootof_roots` | rootof 四根（生产退役） |
| `giac-poly/factor/univariate.rs` | `try_factor_biquadratic` | ℚ 因子分解 |
| `giac-poly/tresultant.rs` | `biquadratic_res_conjugate_pairs` | RT 共轭配对 |
| `giac-calculus/risch/algebraic_rt.rs` | `even_quartic_quadratic_factors` | Risch log 部分 |

### 3.4 giac-poly 内部（~70–95%）

| 符号 | 副本 | 收敛目标 |
|------|------|----------|
| `coeff_wrt` | `subresultant.rs`、`nested.rs::coeff_wrt_impl`、`factor/poly_uni.rs::coeff_wrt_poly` | `subresultant::coeff_wrt` 唯一 |
| `vars_in` | `factor/util.rs`、`subresultant.rs` | `factor/util::vars_in` 唯一 |
| Yun sqff | `univariate.rs`、`factor/poly_uni.rs`、`factor/fpx.rs` | 泛型 `square_free_yun<R>` |
| `integer_sqrt` | `factor/util.rs`、`factor/sparse.rs` | `integer_nth_root(n, 2)` |

### 3.5 求解编排（~50–60%）

| 本地 | 等价权威 |
|------|----------|
| `solve_poly::factor_poly_for_solve` | `giac-poly::factor_univariate_pairs` |
| `realroot::rational_real_roots` | `giac-poly::roots` + 重数追踪 |

---

## 4. 共享 API 设计（签名 + 约束）

### 4.1 层 0 — `giac-core::num_util`

```rust
pub fn ratio_to_expr(r: &Ratio<BigInt>) -> ExprArc;
pub fn integer_nth_root(n: &BigInt, k: u64) -> Option<BigInt>;
pub fn integer_sqrt(n: &BigInt) -> Option<BigInt>; // k = 2
```

- **Tier:** Stable（`giac-core` 公开或 `pub` re-export 自 `algebra::poly::ratio_to_expr`）
- **C\*** 子项:** C0-1、C0-2

### 4.2 层 1 — `giac-poly` 一元视图

```rust
pub fn coeff_wrt(p: &Poly, var: &Var, exp: u64) -> Poly;
pub fn vars_in(p: &Poly) -> Vec<Var>;
pub fn quadratic_abc(p: &Poly, var: &Var)
    -> Option<(Ratio<BigInt>, Ratio<BigInt>, Ratio<BigInt>)>;
pub fn is_biquadratic_monic(p: &Poly, var: &Var)
    -> Option<(Ratio<BigInt>, Ratio<BigInt>)>; // (b, c) for monic t⁴+bt²+c
```

- **Tier:** Stable（`quadratic_abc` 已存在；其余收敛后标 Stable）
- **C\*** 子项:** C1-1

### 4.3 层 2 — 泛型 Yun 无平方因子

```rust
pub trait SquareFreeRing: Clone {
    type Poly;
    fn zero(&self) -> Self::Poly;
    fn one(&self) -> Self::Poly;
    fn derivative(&self, p: &Self::Poly, var: &Var) -> Self::Poly;
    fn gcd(&self, a: &Self::Poly, b: &Self::Poly) -> Self::Poly;
    fn div_exact(&self, a: &Self::Poly, b: &Self::Poly) -> PolyResult<Self::Poly>;
    fn pow(&self, p: &Self::Poly, k: usize) -> Self::Poly;
}

pub fn square_free_yun<R: SquareFreeRing>(
    ring: &R,
    p: &R::Poly,
    var: &Var,
) -> PolyResult<Vec<(R::Poly, usize)>>;
```

- **impl:** `RatioRing`（ℚ）、`PolyModRing`（ℤ/pℤ）
- **Tier:** Stable（泛型核）+ Pipeline private（各 ring impl）
- **C\*** 子项:** C1-2

### 4.4 层 3 — 二次求解（分输出类型）

```rust
pub struct QuadraticCoeffs<C> {
    pub a: C,
    pub b: C,
    pub c: C,
}

pub trait QuadraticRing {
    type Coeff: Clone;
    fn half(&self) -> Self::Coeff;
    fn neg(&self, x: &Self::Coeff) -> Self::Coeff;
    fn add(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff>;
    fn mul(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff>;
    fn div(&self, a: &Self::Coeff, b: &Self::Coeff) -> PolyResult<Self::Coeff>;
    fn discriminant(&self, q: &QuadraticCoeffs<Self::Coeff>) -> PolyResult<Self::Coeff>;
    fn sqrt(&self, d: &Self::Coeff) -> PolyResult<Self::Coeff>;
}

// giac-poly — ℚ 有理根
pub fn quadratic_rational_roots(
    ring: &RatioRing,
    q: &QuadraticCoeffs<Ratio<BigInt>>,
) -> PolyResult<Vec<Ratio<BigInt>>>;

// giac-core — AlgExt + FieldSession
pub fn quadratic_algext_roots(
    session: &FieldSession,
    q: &QuadraticCoeffs<AlgExtCPolyCoeff>,
) -> Result<Vec<AlgExtCPolyCoeff>, EvalError>;
```

- **不强行合并** `quadratic_roots`（ℚ）与 `quadratic_roots_formula`（AlgExt）；共享 `QuadraticCoeffs` + `quadratic_abc` 提取。
- **C\*** 子项:** C2-1（映射 D1-1、D3-2）

### 4.5 层 4 — Expr rootof 构造（giac-core）

```rust
pub fn rootof_from_minpoly(num: &[i64], minpoly: &ExprArc) -> Result<ExprArc, EvalError>;
pub fn quadratic_rootof_branches(poly: &Poly, var: &Var) -> Result<Vec<ExprArc>, EvalError>;
```

- 内部：`quadratic_abc` 校验 → `univariate_poly_to_poly1_expr` → `AlgExtData::from_rootof`
- **C\*** 子项:** C0-3（映射 D3-2）

### 4.6 层 5 — 求解编排

- `giac-solve` 直接调用 `giac_poly::factor_univariate_pairs`（或 thin pub alias），删除 `factor_poly_for_solve`。
- **C\*** 子项:** C3-1

---

## 5. 各 crate 收缩职责

### giac-core

| 删除 | 保留 |
|------|------|
| `eval.rs` 私有 `ratio_to_expr`、`integer_sqrt` | `algebra/poly.rs` — Expr↔Poly 唯一权威 |
| `field_arith::ratio_to_expr_arc` | `algebra/poly_roots.rs` — `poly_algext_roots` deg 1–4 |
| `field_session::is_negative_rational` | `algebra/{ext_tower,field_session,alg_ext*}` |
| `split_depressed_quartic`（D4-3） | `eval_poly.rs` — giac-poly eval 薄包装 |
| | Plugin trait 定义；`num_util` 层 0 工具 |

### giac-poly

| 删除 | 保留 |
|------|------|
| `nested::coeff_wrt_impl`、`poly_uni::coeff_wrt_poly` | `Poly` / `nested` 环视图 |
| `subresultant::vars_in`（私有副本） | `subresultant_gcd`、`factor/*` |
| `sparse::integer_sqrt_bigint` | `quadratic_abc`、`roots`（ℚ deg≤2） |
| 独立 Yun 三份实现 | `partfrac`、`sturm`、`tresultant`、`dense/poly1` |

### giac-solve

| 删除 | 保留 |
|------|------|
| `biquadratic_rootof_roots` 生产路径 | `solve_univariate_over_q` 主编排 |
| 双份 `rootof_expr` | `poly_algext_roots_for_ctx` 回退（D3-1 ✅） |
| `quadratic_rootof_roots` 系数循环 | `sturm`、`realroot`、`fsolve` |
| `factor_poly_for_solve` | `rootof.rs` 仅测试/reference（最终可删） |
| `realroot::rational_real_roots` 线性扫描 | |

### giac-calculus

| 删除 | 保留 |
|------|------|
| 4 份 `ratio_to_expr` | `diff`、`integrate*`、`limit_engine/*` |
| `partfrac_integrate::integer_sqrt` | `partfrac_integrate` Hermite/RT **积分**层 |
| `algebraic_rt` 重复 `sqrt_rational_coeff_radicand`（若与 tresultant 等价） | `risch/*`、`SparseSeries` |

---

## 6. 子任务表

| ID | 标题 | 优先级 | 落点 | 映射 dedup-plan | 状态 |
|:--:|------|:--:|------|-----------------|:--:|
| C0-1 | `ratio_to_expr` 收拢到 `giac-core` | P0 | `algebra/poly.rs` + calculus 删 4 副本 | 扩 D0 | **done** (PR-A) |
| C0-2 | `integer_nth_root` / `integer_sqrt` 收拢 | P0 | `giac-core/num_util` | 扩 D0 | **done** (PR-B) |
| C0-3 | `rootof_from_minpoly` + solve 删双份 | P1 | `giac-core/algebra/alg_ext.rs` | D3-2 | open |
| C1-1 | `coeff_wrt` / `vars_in` giac-poly 唯一化 | P1 | `subresultant.rs`、`factor/util.rs` | — | **done** (PR-C) |
| C1-2 | `square_free_yun<R: SquareFreeRing>` | P2 | `giac-poly/univariate.rs` | 邻 D5 | open |
| C2-1 | `QuadraticCoeffs` + `quadratic_rational_roots` | P1 | `giac-poly/quadratic.rs`（新建） | D1-1 | partial |
| C2-2 | `quadratic_algext_roots` 抽自 `poly_roots` | P2 | `giac-core/algebra/poly_roots.rs` | D1-3 | open |
| C3-1 | `factor_poly_for_solve` → `factor_univariate_pairs` | P1 | `giac-solve/solve_poly.rs` | — | open |
| C3-2 | 退役 `biquadratic_rootof_roots` 生产路径 | P1 | `giac-solve/rootof.rs` | D3-2, [T3-3](GIAC-poly-algext-gcd-factor-priority.md) | open |
| C3-3 | `realroot` 改调 `giac_poly::roots` | P2 | `giac-solve/realroot.rs` | — | open |
| C4-1 | calculus 删工具拷贝 | P0 | `limit_engine/*`、`partfrac_integrate.rs` | C0-1/C0-2 | **done** (PR-A/B) |
| C4-3 | 删 `split_depressed_quartic` | P2 | `giac-core/poly_roots.rs` | D4-3 | open |
| C4-4 | 删 `field_session::is_negative_rational` 副本 | P1 | `giac-core` | — | open |

---

## 7. 推荐 PR 切片

| PR | 子项 | 预估 | 行为变化 |
|:--:|------|:--:|:--:|
| PR-A | C0-1 + C4-1（ratio_to_expr） | 小 | 无 | **done** |
| PR-B | C0-2 + C4-1（integer_nth_root） | 小 | 无 | **done** |
| PR-C | C1-1 | 中 | 无 | **done** |
| PR-D | C2-1 + C0-3 + C3-1 | 中 | 无 |
| PR-E | C3-2 + C4-3 + C4-4 | 中 | 无（能力已由 D3-1 覆盖） |
| PR-F | C1-2 | 大 | 无 |

**门禁（每 PR）：**

```bash
cd giac-rs && cargo test-timeout && cargo ci-clippy
```

---

## 8. 纪律

1. **禁止** 在 calculus/solve 新增第三个二次/四次形状分支。  
2. **禁止** 扩大 `giac-simplify` 对 `giac-solve` 的依赖；rootof 构造归 `giac-core`。  
3. **必须** Expr↔Poly 仅经 `giac-core::algebra::poly`。  
4. 新增/改动 fn 标 tier；合入前按 [algorithm-expr-api.md](../algorithm-expr-api.md) §7.2 复审。

---

## 9. 参考

- [GIAC-rs-crate-dedup-plan](GIAC-rs-crate-dedup-plan.md) — D0–D5 母计划  
- [GIAC-poly-algext-gcd-factor-priority](GIAC-poly-algext-gcd-factor-priority.md) — K 上 factor 后 T3-3 跨 crate 清理  
- [GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md) — API 分层母索引  
- [giac-poly-api-stability.md](../giac-poly-api-stability.md) — `quadratic_abc`、factor FAC-G*
