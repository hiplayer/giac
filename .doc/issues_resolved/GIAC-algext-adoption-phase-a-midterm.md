# GIAC-AlgExt — Phase A + 中期基座（已归档）

**状态:** **resolved**（2026-06-23 代码快照）  
**续篇（仍 open）:** [GIAC-algext-adoption.md](../issues/GIAC-algext-adoption.md)  
**细项 backlog:** [GIAC-poly-algext-backlog.md](../issues/GIAC-poly-algext-backlog.md)、[GIAC-lazy-common-tower-plan.md](../issues/GIAC-lazy-common-tower-plan.md)  
**Rust 落点:** `giac-core::algebra::{alg_ext, alg_ext_c, ext_tower, field_session, poly_roots, poly_alg_coeff}`

---

## 1. 归档范围

本文件记录 **GIAC-algext-adoption** 中已验收交付：实代数数 MVP、塔式 `ExtensionField`、跨域 `common`、复代数标量 `AlgExtC` 内核、以及 `Poly<AlgExtCPolyCoeff>` 上 deg≤4 求根管线。未纳入：通用 `gcd`/`factor` over K、partfrac disc>0、`evalf`、参数 `PolyCoeff`、assume 引擎、删除全部过渡态特判。

---

## 2. 已交付数据模型（2026-06-23）

```rust
// giac-core — 已落地，非过渡态 min_poly-only
pub struct AlgExtData {
    pub field: Arc<ExtensionField>,
    pub coords: Vec<ExprArc>,
    pub root_index: Option<u32>,
}
// min_poly() 由 field.layer_min_poly_exprs() 导出

pub struct AlgExtCData {
    pub field: Arc<ExtensionField>,
    pub re: Vec<ExprArc>,
    pub im: Vec<ExprArc>,
    pub root_index: Option<u32>,
}
// Expr::AlgExt | Expr::AlgExtC

pub type PolyAlgExt = Poly<AlgExtCPolyCoeff>;  // giac-core::algebra::poly
```

**塔（阶段 2，S0–T4b）：** `ExtensionTower::{Base, Adj}`、`ExtensionField::common_over_q`（含 cache）、`align_elements` / `embedding_for`、`adjoin_irreducible`。详见 [GIAC-lazy-common-tower-plan.md](../issues/GIAC-lazy-common-tower-plan.md)（T3+ adjoin \(u^2-\alpha\) 仍 open）。

---

## 3. 功能清单（已验收）

| 类别 | 项 | 落点 | 验收 |
|------|----|------|------|
| **P0** | A-01 `eval` 同域/跨域 `AlgExt` 四则 | `alg_ext.rs` + `fold_algext_*` | `fold_algext_sum` lazy common(√2,√3)；`add` 经 `align_pair` |
| **P0** | A-02 `eval_frac` / `inv` | `eval.rs`, `AlgExtData::inv` | 分母 `AlgExt` → `inv`；零除 `DivisionByZero` |
| **P0** | A-03 `ratnormal` | giac-simplify | `AlgExt` 原子常数 |
| **P0** | A-04 `assert_equiv` | giac-simplify | 同域 `eq_mod` |
| **P0** | A-05 `display` / `to_rootof_expr` | `display.rs` | `AlgExt` + `AlgExtC` golden |
| **P1** | B-01 `expr_to_poly` 边界 | `poly_conv.rs` | 含代数系数 → `TypeError`；提升走 `poly_alg_from_expr` |
| **P1** | B-05 `ExtensionTower::common` | `ext_tower.rs` | `(√2)+(∛2)` deg 6；T4a 塔默认 |
| **P1** | B-06 **内核** `AlgExtC` | `alg_ext_c.rs`, `expr.rs` | `add/mul/inv/eq_mod`；`canonicalize_to_algext_c` |
| **P1** | B-02/03 **子集** | `poly_roots.rs`, `solve_poly.rs` | `poly_algext_roots` deg 1–4；solve deg≤4 走此路径 |
| **P1** | B-04 **子集** | `realroot.rs` | ℚ 上 sturm + `canonicalize_to_algext_c` 端点 |
| **表示** | `PolyCoeff` + `Poly<C>` | `giac-poly/poly_coeff.rs`, `nested.rs` | M1 ✅ |
| **桥接** | `poly_alg_from_expr` / `algext_poly_to_expr` | `giac-core::algebra::poly` | P1-4 ✅ |
| **过渡** | `fold_complex_algext_*` | `alg_ext.rs` | `Complex`+`AlgExt` 混合运算（待 2b eval 接线后删） |
| **过渡** | `algext_sqrt_branches` / `biquadratic_rootof_roots` | `alg_ext.rs`, `giac-solve/rootof.rs` | solve 主路径已迁 `poly_algext_roots`；文件仍保留测试/兼容 |
| **过渡** | `try_factor_quadratic_rootof` | giac-simplify | `factor(x²−2)` 钩子 |

---

## 4. 求根管线（P2-1/2/6 + P3-6 子集）

`poly_algext_roots` / `poly_algext_roots_for_ctx`（`poly_roots.rs`）：

| deg | 策略 | 测试代表 |
|-----|------|----------|
| 1 | 线性 | — |
| 2 | 二次 + √Δ | `roots_quadratic_x2_minus_2` |
| 3 | Cardano / pure cubic | `roots_cubic_t3_minus_2` |
| 4 | 双二次或 resolvent cubic | `roots_biquadratic_t4_minus_2`、quartic adjoin/deflate |

`solve_univariate_over_q`（`giac-solve/solve_poly.rs`）：sqff × factor → 每因子 `poly_algext_roots_for_ctx`（deg≤4）；deg≥5 不可约 → 单支 `rootof`。

**`PolyInK` + `FieldSession`：** 规范化与 `verify_root` 共用 prepare 路径（expr-api 1B）。

---

## 5. 阶段对照（§8.7 已关闭项）

| 阶段 | 状态 | 说明 |
|------|------|------|
| **A** | ✅ | `AlgExtData` + 同域/跨域运算 + eval 提升 |
| **2** | ✅ | `ExtensionField` + 塔 + `common` 缓存（T3+ 除外） |
| **2b** | **部分** | 类型与 `canonicalize` ✅；`eval` 未对 `AlgExtC` 做 `+−×÷` fold；`i` 未默认进塔 |
| **3–5** | open | 见主 issue |

---

## 6. 仍属过渡态（未在本归档关闭）

- `Expr::Complex(0, AlgExt)` 与 `AlgExt` / `AlgExtC` 并存
- `eval` 对 `AlgExtC` 仅原子传递，未走 `AlgExtCData::add/mul`
- `eval_frac` 未处理 `AlgExtC` 分母
- `giac-solve/rootof.rs` 二次/双二次特判仍存在
- partfrac disc>0 → `NotImplemented`
- 无 `AlgExtC::evalf`

---

## 7. 测试入口

```bash
cargo test -p giac-core alg_ext
cargo test -p giac-core ext_tower
cargo test -p giac-core poly_roots
cargo test -p giac-core alg_ext_c
cargo test -p giac-solve solve_poly
cargo test -p giac-solve rootof
```

conformance：`batch2_giac205_solve_quadratic_rootof`；`solve(t^4-2=0,t)` 四根。
