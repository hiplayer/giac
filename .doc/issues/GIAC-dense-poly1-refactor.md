# GIAC — 稠密 poly1 算术抽象（`field_arith` ↔ `giac-poly`）

**状态:** open  
**类型:** 重构计划 / AFK  
**相关:** [GIAC-lazy-common-tower-plan](GIAC-lazy-common-tower-plan.md) §12.9、[GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md)、[giac-poly-api-stability.md](../giac-poly-api-stability.md)、[expr-poly-conversion.md](../expr-poly-conversion.md)、[giac-tower-common-math.md](../giac-tower-common-math.md)  
**Rust 落点:** 新建 `giac-poly::dense::poly1`；收敛 `giac-core::algebra::field_arith`  
**快照:** 2026-06-20

---

## 1. 问题陈述

`giac-core::algebra::field_arith` 与 `giac-poly::univariate` **各自维护一套一元多项式算术**，算法同构、表示不同：

| 能力 | `field_arith` | `giac-poly` |
|------|---------------|-------------|
| 系数向量顺序 | **高次在前**（giac `poly1`） | **升幂**（`univariate_coeffs_ascending`） |
| 载体 | `CoordsQ = Vec<Ratio<BigInt>>` | sparse `Poly` 或内部 `Vec` |
| ±× / div/rem | `poly_add` … `poly_divrem` | `univariate_div_rem` 等 |
| mod / 域逆 | `poly_reduce`、`poly_inv_mod` | factor 管线用 remainder，无直接等价 |
| extgcd | `poly_ext_gcd` | `egcd`（sparse）、`pseudo_remainder` |

此外 `field_arith` **内部双份**：ℚ 常系数路径与 T3 `ParentCoeffRing` / `*_with_coeffs_in_field` / `*_blocks` 几乎平行实现（~600 行）。

**影响：**

- `ext_tower` / `AlgExtData` 热路径依赖稠密 `CoordsQ` + `poly1` 约定；与 sparse `Poly` 坐标不可混用（plan §12.9）。
- **T3+**（`adjoin(K, u²−α)`）与 **P3-6**（`Poly<AlgExtC>::roots`）若在现有双份代码上扩展，会继续复制 div/rem/reduce。
- 环运算按 [module-division.md](../module-division.md) 应对标 `ezgcd.cc` / `modpoly.cc`，归属 **`giac-poly`**，而非在 core 无限膨胀。

---

## 2. 目标与边界

### 2.1 目标

1. 新建 **`giac-poly::dense::poly1`**：稠密 `Vec<C>` 一元算术 + 显式 **poly1 顺序**。
2. 引入 **`Poly1RingOps`** 泛型环，统一 ℚ 路径与 T3 parent 系数路径（删除 `*_blocks` 副本）。
3. 提供 **高次在前 ↔ 升幂 ↔ sparse** 的可测转换 API，落实 §12.9 坐标审计。
4. `field_arith` 收敛为 **薄 facade + Expr/poly1 工具 + 扩域专用矩阵**（`char_poly_matrix` 等）。

### 2.2 非目标

- 不用 sparse `Poly` 替换 `ext_tower` 热路径上的 `CoordsQ`。
- 不把 `AlgExtCPolyCoeff` 迁入 `giac-poly`（环依赖；见 `poly_coeff.rs` 注释）。
- 不合并 `Expr`↔`Poly` 桥接（留在 `giac-core::algebra::poly` / `poly_conv`）。
- 不在本 issue 统一工作区 `ratio_to_expr` 私有副本（可另开）。
- 不重构 `char_poly_matrix` / `mult_matrix_of_element` 等矩阵块（与 poly 重复度低，另评估 `giac-linalg`）。

---

## 3. 架构

```text
┌─────────────────────────────────────────────────────────┐
│ giac-core::algebra                                      │
│  ext_tower / alg_ext  — 扩域语义、嵌入、registry         │
│  poly / poly_conv     — Expr ↔ Poly 边界（Stable）       │
│  field_arith          — 薄 facade + Expr/poly1 + 矩阵     │
└───────────────────────────┬─────────────────────────────┘
                            │ CoordsQ, ParentCoeffRing
┌───────────────────────────▼─────────────────────────────┐
│ giac-poly::dense::poly1  （新建）                        │
│  Poly1RingOps + 泛型 dense 算术                          │
│  Poly1Order: HighFirst | Ascending                       │
│  convert: ↔ sparse ascending                             │
└───────────────────────────┬─────────────────────────────┘
                            │ Ratio<BigInt>
┌───────────────────────────▼─────────────────────────────┐
│ giac-poly 现有层                                         │
│  sparse Poly, FlatUni, univariate, factor, roots…        │
└─────────────────────────────────────────────────────────┘
```

**职责**

| 层 | 内容 |
|----|------|
| `giac-poly::dense::poly1` | 稠密一元 ±×、div/rem、mod reduce、extgcd、inv_mod |
| `giac-core::field_arith` | `CoordsQ` 别名、`ParentCoeffRing`、`RatioRingOps` / `ParentBlockRing` adapter、Expr/`poly1` 解析、Newton 特征多项式、嵌入矩阵 |
| `giac-core::ext_tower` | 仅经 `field_arith` 公开面调用，不直接依赖 `dense` |

---

## 4. 核心 API（草案）

### 4.1 `Poly1RingCtx`（实现名；原草案称 `Poly1RingOps`）

静态 `zero()`/`one()` 无法携带 T3 `ParentCoeffRing` 上下文，故采用 **实例 trait** `Poly1RingCtx`（D3）。

```rust
// giac-poly/src/dense/poly1.rs — Stable (crate-internal)

pub trait Poly1RingCtx {
    type Coeff: Clone;
    fn zero(&self) -> Self::Coeff;
    fn one(&self) -> Self::Coeff;
    fn is_zero(&self, c: &Self::Coeff) -> bool;
    fn add / sub / neg / mul / inv(&self, ...) -> PolyResult<...>;
}

pub enum Poly1Order { HighFirst, Ascending }

pub fn trim / add / sub / mul / ...<R: Poly1RingCtx>(ctx: &R, ..., order: Poly1Order);
```

**实现（Phase D2–D3 ✅）** — **依赖方向：** `giac-poly` 不得依赖 `giac-core`；环 adapter 留在 core。

| 位置 | Impl | 元素类型 | 用途 |
|------|------|----------|------|
| `giac-poly::dense::poly1` | 泛型 `add` / `div_rem` / `reduce_mod_monic` / … | `R: Poly1RingCtx` | 算法唯一实现 |
| `giac-poly::dense::ratio_ring` | `RatioRingCtx` | `Ratio<BigInt>` | ℚ 系数 ctx |
| `giac-core::field_arith` | `ParentBlockRing<'a>` | `CoordsQ`（块 = parent 坐标） | 包装 `ParentCoeffRing` → 调 dense 泛型 |
| `giac-core::field_arith` | `ParentCoeffRing<'a>` | 闭包表 | **保留在 core**（指向 `ext_tower` `element_*`） |

**与 `PolyCoeff` 关系：** 不合并。`Poly1RingCtx` 服务稠密 `Vec`；`PolyCoeff` 服务 sparse 项系数。P3-6 若需 dense 层，再评估 `AlgExtCPolyCoeff` adapter（Phase D5，仍在 core）。

### 4.2 约定转换（Phase D4 ✅）

```rust
// giac-poly/src/dense/convert.rs — Stable (crate-internal)

pub fn reverse_coeffs(c: &[Ratio<BigInt>]) -> Vec<Ratio<BigInt>>;
pub fn ascending_to_dense_high_first / dense_high_first_to_ascending(...);
pub fn sparse_ascending_to_dense_high_first(p: &Poly, var: &Var) -> Vec<Ratio<BigInt>>;
pub fn dense_high_first_to_sparse(p: &[Ratio<BigInt>], var: &Var) -> Poly;
```

文档同步：[GIAC-lazy-common-tower-plan](GIAC-lazy-common-tower-plan.md) **§12.9** 坐标基表 + [giac-tower-common-math.md](../giac-tower-common-math.md) **§4**（可审计验证）Rust API 索引。

### 4.3 `field_arith` 目标形态

```text
field_arith.rs (~250–320 行；矩阵块仍留 core)
├── pub type CoordsQ = Vec<Ratio<BigInt>>
├── RatioRingOps + ParentBlockRing adapter → giac_poly::dense::poly1
├── pub fn poly_*  → dense 泛型 (RatioRingOps, HighFirst)
├── pub fn poly_*_with_coeffs_in_field → dense 泛型 (ParentBlockRing, HighFirst)
├── ParentCoeffRing<'a> 结构体（core）
├── Expr: poly1_coeffs, coords_to_expr, rationalize_poly1, canonical_poly1_expr, …
└── 扩域专用: char_poly_matrix, mult_matrix_of_element, kron_left, embed_in_square_extension, …
```

**删除：** `trim_leading_zero` + `trim_leading_zero_blocks`；`poly_divrem` + `poly_divrem_blocks` 等双份实现。

### 4.4 `poly_reduce` 语义约束（必须保留）

- **ℚ 路径**（`poly_reduce` / `RatioRingCtx`）：模多项式首项 **monic 或 ±1**。
- **T3 路径**（`poly_reduce_with_coeffs_in_field` / `ParentBlockRing`）：层 minpoly 经 `layer_minpoly_parent_coeffs` 嵌入 parent 后 **首项 block 恒为 `ring.one`**；generic `reduce_mod_monic` 虽接受 leading **−1**，T3 实际上不会触发该分支。
- 不泛化为任意 leading coeff 的通用 reduce，除非单独 issue + 全量回归。

---

## 5. 与路线图衔接

| 里程碑 | 关系 |
|--------|------|
| **T3+** | **登记** `adjoin(K, u²−α)` 可与 D1 并行；**层上 minpoly 约化 / `poly_*_with_coeffs_in_field`** 硬依赖 **D3** |
| **P3-6** | D5 评估 `FlatUni<AlgExtCPolyCoeff>` ↔ dense；roots 主路径仍 sparse + `PolyCoeff` |
| **塔 plan DoD** | 不改变 lazy common / `align_elements` 行为 |
| **expr-poly 边界** | `expr_to_poly` / `poly_alg_from_expr` 不变 |

**P3-6 分层策略**

| 场景 | 路径 |
|------|------|
| ext_tower 域内 × / inv | `CoordsQ` + dense high-first |
| `Poly<AlgExtC>::roots` resolvent | `FlatUni<C>` + sparse |
| 互操作 | 仅经 §4.2 转换 API；禁止 flatten θ 坐标直接挂塔句柄 |

---

## 6. 分阶段实施

| Phase | ID | 内容 | 验收 | PR |
|-------|-----|------|------|-----|
| **D0** | — | 本文档 + API tier 登记 | 索引可发现 | — |
| **D1** | **D1** | 新建 `giac-poly::dense::poly1` + `RatioRingCtx` + `HighFirst`；**新建** dense 单测 | `cargo test -p giac-poly dense::*` | PR1 ✅ |
| **D2** | **D2** | `field_arith` 改调 dense；删 ℚ 重复实现 | `cargo test -p giac-core` 全绿 | PR2 ✅ |
| **D3** | **D3** | `ParentBlockRing` 统一 T3；删 `*_blocks` | `t1b_k2_*`、`element_mul` 塔测 | PR3 ✅ |
| **D4** | **D4** | 约定转换 API + tower plan §12.9 / math §4 文档 + roundtrip 测 | `fold_algext_sum_rat_on_k2_*` + `dense::*` roundtrip | PR4 ✅ |
| **D5** | **D5** | P3-6 前：`FlatUni` ↔ dense 桥接评估 | 随 P3-6 立项 | 可选 |

**推荐 PR 顺序：** D1 → D2 → D3 → D4；（D5 不挡 T3+）。

**AFK 竖切：** 仅 **D1 + D2** 即可消除 ℚ 路径最大重复；**D3 须在 T3+ 层算术验收前完成**（勿在 `*_blocks` 副本上扩展 T3+）。

**T3+ 并行约定：** adjoin **登记** PR 可与 D1 并行；T3+ **验收**（层 minpoly 上 `element_*` / reduce）合并顺序：**D3 先于或同 PR 于** T3+ 算术验收。

---

## 7. 测试策略

1. **新建 dense 单测（D1）：** `giac-poly/src/dense/tests.rs` — 从 `field_arith` **行为**提取（非移植，因 core 侧无 mod tests）；固定小多项式锁死 **HighFirst** 下标语义（模块头注释引用 adoption `poly1` 定义）。
2. **约定 roundtrip（D4）：** 同一多项式 high-first ↔ ascending ↔ sparse 系数逐项相等。
3. **集成（D2+）：** 保留 `ext_tower` 全量；PR2 过渡期可 feature-gate **双跑**（旧/新 `poly_reduce` `eq_mod`），**D3 合并后删除**双跑 feature。
4. **属性（可选）：** 随机小次数 ℚ 多项式，dense `div_rem` 与 `univariate_div_rem` **经 §4.2 转换后**余式为零（两套 div/rem 不可字节级互换）。

---

## 8. API 稳定性

| 符号 | Tier | 说明 |
|------|------|------|
| `dense::poly1::*` | **Stable (crate-internal)** | 初版仅 `giac-core` 使用 |
| 约定转换 fn | **Stable (crate-internal)** | D4 后文档化 |
| `field_arith` 公开签名 | **不变** | facade 保持 `ext_tower` 稳定 |

合入前更新 [giac-poly-api-stability.md](../giac-poly-api-stability.md) Per-file 表（[algorithm-expr-api.md §7.2](../algorithm-expr-api.md) 复审清单）。

---

## 9. 风险与缓解

| 风险 | 严重度 | 缓解 |
|------|--------|------|
| 高/低次顺序 bug | 高 | D4 roundtrip 测；PR2 短期双跑 |
| `poly_reduce` monic 假设破坏 | 中 | 原注释 + `debug_assert`；不随意泛化 |
| `Poly1RingOps` vs `PolyCoeff` 两套 trait | 中 | 文档写清；D5 再评估 adapter |
| giac-poly 编译体积 | 低 | `mod dense` 独立；无新依赖 |

---

## 10. 明确不采用的方案

| 方案 | 原因 |
|------|------|
| dense 留在 `giac-core::dense_poly1` | 与 module-division「环运算在 poly」不一致 |
| macro 生成 ℚ / blocks 双份 | T3+ 仍难维护 |
| ext_tower 全面改 sparse `Poly` | 热路径分配 + §12.9 约定冲突 |
| 强行合并 `univariate_div_rem_wrt` 与 dense | 系数环为 ℚ[others]，语义不同 |

---

## 11. 完成定义（DoD）

- [x] `giac-poly/src/dense/poly1.rs` 存在且 `cargo test -p giac-poly` 含 dense 测
- [x] `field_arith` 无独立 `poly_add`/`poly_divrem`/`poly_reduce` 算法体（仅 facade）
- [x] T3 `poly_*_with_coeffs_in_field` 经单一 `Poly1RingCtx` 路径（D3）
- [x] §4.2 转换 API + roundtrip 测 + [GIAC-lazy-common-tower-plan §12.9](GIAC-lazy-common-tower-plan.md) / [giac-tower-common-math §4](../giac-tower-common-math.md) 索引更新（D4）
- [x] `giac-poly-api-stability.md` Per-file 登记 `dense`
- [x] PR2 双跑 feature（若有）在 **D3 合并后删除**（未引入双跑 feature）
- [x] `cargo test -p giac-core` 全绿；塔相关测无回归

---

## 12. Review 结论（2026-06-20）

**建议立项 D1–D4**；**D5 挂 P3-6**。

**优先级：** D1 → D2 → **D3（T3+ 层算术验收前）** → D4；避免在 `*_blocks` 副本上继续堆 adjoin 层逻辑。

**Review 修订（2026-06-20）：** `ParentBlockRing` / `ParentCoeffRing` 留 **giac-core**；dense 单测 **新建**非移植；文档索引指向 **tower plan §12.9**。

**已确认分层：** Expr 桥接留 core；环运算下沉 poly；扩域 registry/嵌入留 core — 与 [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) Phase 0–1 不冲突。
