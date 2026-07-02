# P2 代数数论深化 — 做完整后可暴露的 API/功能

**状态:** open（跟踪 / 评估）
**类型:** 扩展功能规划（非 upstream 对齐）
**上游基线:** **`giac/giac-2.0.0` 无对标**（Pari/GP 才有完整类群 / 单位群 / LLL）
**相关:** [GIAC-core-upstream-gaps](GIAC-core-upstream-gaps.md) §4 Phase E（C-9..C-12）、[GIAC-poly-f5-fglm-hasse-lean4-verification](GIAC-poly-f5-fglm-hasse-lean4-verification.md)（Hasse 数学正确性线）
**Rust 落点:** `giac-rs/crates/giac-core/src/algebra/{class_group,unit_group,lattice,archimedean,number_field_arith,galois_automorphism,padic}.rs`
**快照:** 2026-07-02

---

## 问题陈述

[GIAC-core-upstream-gaps](GIAC-core-upstream-gaps.md) §4 Phase E（C-9..C-12）当前全为 `pub(crate)`，**仅服务四次求根的 `hasse_sqrt`**（`poly_roots.rs:2089` 唯一外部调用），**无用户级 API 暴露**。

本 issue 跟踪：**P2 代数数论基础设施做完整化后，可暴露为 giac 命令级函数 / 功能的清单**，以及对应的工作量评估与优先级。这些是 **giac-rs 比 upstream giac 多的差异化扩展**，对标 Pari/GP 命令，非 conformance 硬阻塞。

**定位：** upstream giac-2.0.0 **没有** `class_number` / `class_group` / `lll` / `is_principal` 等命令（Pari/GP 才有）。因此 P2 做完整不是「补齐 upstream 缺口」，而是「giac-rs 在代数数论上的扩展线」。

---

## 现状（snapshot 2026-07-02）

| 模块 | 行数 | 现状 pub API | 唯一外部调用 |
|------|------|-------------|-------------|
| `class_group.rs` | 666 | `pub(crate) ideal_is_principal` / `hasse_is_square` / `hasse_sqrt` | `poly_roots.rs:2089 hasse_sqrt` |
| `unit_group.rs` | 1046 | `pub(crate) Torsion` / `torsion` / `fundamental_units` | 内部 |
| `lattice.rs` | 449 | `pub(crate) lll` / `lll_with_transform` / `lll_reduce_overcomplete` / `lll_incremental` | 内部 |
| `archimedean.rs` | 803 | `pub(crate) Signature` / `field_signature` / `is_totally_positive` / `Embeddings` / `embeddings` / `eval_at_embedding` | 内部 |
| `number_field_arith.rs` | 1294 | `pub(crate) PrimeIdealRec` / 理想赋值扫描 | 内部 |
| `galois_automorphism.rs` | 303 | `pub(crate) conjugate_map_sends` / `try_galois_sqrt_second` | 内部 |
| `padic.rs` | 359 | `pub(crate) mod_pk` / `vp_bigint` / `rat_poly_mod_pk` / `poly_*mod_pk` | 内部 |

**`ponytail:` 现状边界：** `IDEAL_GEN_COORD_BOUND`（per degree）+ `N_J_MAX = 10⁷`；超过 → `None`（sound，**不区分「非主」vs「界太小」**）。

---

## 做完整后可暴露的函数 / 功能

| ID | 做完整后暴露的 giac 命令级函数 | 功能 | 对标 Pari/GP | 现状阻塞 |
|----|------------------------------|------|-------------|---------|
| **C-9** | `class_number(P)` | 数域 `ℚ(α)`（`P = α` 的 minpoly）的类数 `h_K` | `quadclassunit` / `bnfclassunit` | 有界生成元搜索 + sound-skip；超界 `None` 不区分「非主」vs「界太小」 |
| | `class_group(P)` | 类群结构 `Cl(K) ≅ C_{n₁} × … × C_{n_g}` + 生成元理想 | `bnfclassunit` | 同上 |
| | `is_principal(ideal)` | 理想是否主理想（含证书） | `bnfisprincipal` | `ideal_is_principal` 已 `pub(crate)`，可直接升 `pub` |
| **C-10** | `class_number_high_degree(P)` | deg≥3 数域类数（实二次 `r=1, d≥3` 现 `None`） | `bnfclassunit` | 完整 Step 7b Buchmann 亚指数 |
| | `certify_non_principal(P, ideal)` | 非主性证书（高次域） | — | 仅 deg-2 虚二次有穷举证书 |
| **C-11** | `lll(matrix)` | LLL 短向量约化基 | `qflll` / `qflllgram` | **✅ 已落地**（`lattice.rs::lll` 经 `eval_lll` 接 eval：`lll(matrix)` 命令，矩阵行→f64→LLL→整数/decimal `Rat` 还原；测试 `lll_command_reduces_2d_basis` `[[1,2],[3,4]]→[[1,0],[0,±2]]`）；未接类群/单位群生成元搜索（C-9 scope） |
| | `short_vector(lattice)` | 格最短向量 | `qfminim` | LLL 接入 `class_group` / `unit_group` 搜索（C-9 scope） |
| **C-12** | `evalf(rootof(α,P))` | 代数数 → 浮点（`rootof` / `AlgExt` / `AlgExtC` 求值） | `Mod` / `bestappr` | **✅ 已落地**（`archimedean::algext_evalf`/`algextc_evalf` + `eval_evalf` dispatch：`evalf(AlgExt/AlgExtC)` 经 embeddings[root_index.unwrap_or(0)] Horner 求值→decimal `Rat`/`Complex`；`evalf(sqrt(numeric))`→f64 sqrt；Add/Mul 折叠；`ponytail:` 无 `Float` variant→floats 为 `Rat` denom=10^digits，exact 有理数保持精确；embedding-0=最小实根（sqrt-like 域可能为负根），升级路径=rootof 带 index/evalf 取 index 参）；解锁 `evalf` 命令对代数数；测试 `evalf_algext_sqrt2_embedding0_is_negative_root` / `evalf_sqrt_of_two_is_positive` / `evalf_add_of_sqrt2_and_one` |

### eval dispatch 接线草图

P2 做完整后，`giac-core` eval 主路径新增分支（接 `Context` plugin trait）：

```rust
// crates/giac-core/src/eval.rs
FuncKind::ClassNumber => eval_class_number(args, ctx),  // class_number(poly) → BigInt
FuncKind::ClassGroup  => eval_class_group(args, ctx),   // class_group(poly)  → [n1,...,ng]
FuncKind::IsPrincipal => eval_is_principal(args, ctx),  // is_principal(ideal) → bool
FuncKind::Lll         => eval_lll(args, ctx),           // lll(matrix) → 约化基
FuncKind::Evalf       => evalf 已有，对 AlgExt/RootOf 分支补 C-12
```

---

## 价值与定位

- **C-12 `evalf` 最实用**：解锁 `evalf(rootof)`、数值验证代数数运算、`approx_rootof`，是 conformance / 用户交互高频需求。**独立可做**（不依赖 C-9/10/11）。
- **C-11 LLL 最易暴露**：`lll` 已实现且 `pub(crate)`，升 `pub` + 接 eval 即成 `lll(matrix)` 命令，**工作量最小**。
- **C-9/10/11 类群 / LLL 搜索**：数论研究级功能，giac 用户面窄，但补齐后 giac-rs 在代数数论上**超过 upstream giac**（对标 Pari）。现状已为四次 √-判定做了 sound-skip 版，做完整是「有界 → 完整证书」的升级。

---

## 建议优先级

```text
1. C-12 AlgExtC::evalf          独立、高频、解锁 conformance 代数数数值验证
2. C-11 LLL pub 暴露 + eval 接线  已实现，升 pub + eval 接线，工作量最小
3. C-9 完整 Buchmann 类群        数论深化，需亚指数算法，大工程
4. C-10 高次非主认证             依赖 C-9
```

P2 整体**非 conformance 硬阻塞**（upstream giac 也没这些），是 giac-rs 的差异化扩展线。按 Hasse lean4 验证线（[GIAC-poly-f5-fglm-hasse-lean4-verification](GIAC-poly-f5-fglm-hasse-lean4-verification.md)）的数学正确性约束推进。

---

## 分阶段验收

### 阶段 1 — C-12 `AlgExtC::evalf`（独立可做）

- [x] `archimedean::algext_evalf` / `algextc_evalf`：embeddings[root_index.unwrap_or(0)] Horner 求值 → `nalgebra::Complex<f64>`
- [x] `eval_evalf` dispatch + `FuncKind::Evalf` + parser/display `evalf`
- [x] `evalf(AlgExt/AlgExtC)` → decimal `Rat` / `Complex`；`evalf(sqrt(numeric))` → f64 sqrt；Add/Mul 折叠
- [x] `ponytail:` 无 `Float` variant → floats 为 `Rat` denom=10^digits（exact 有理数保持精确）；embedding-0=最小实根（升级=rootof 带 index / evalf 取 index 参）
- [x] 单元测试：`evalf_algext_sqrt2_embedding0_is_negative_root`、`evalf_sqrt_of_two_is_positive`、`evalf_add_of_sqrt2_and_one`、`evalf_algext_constant_is_exact`、`evalf_rational_stays_exact`

### 阶段 2 — C-11 LLL pub 暴露

- [x] `eval_lll` dispatch + `FuncKind::Lll` + parser/display `lll`（`lattice::lll` 保持 `pub(crate)`，eval.rs 同 crate 调用）
- [x] `lll(matrix)` giac 命令：矩阵行→f64→LLL→整数/decimal `Rat` 还原
- [x] 单元测试：`lll_command_reduces_2d_basis`（`[[1,2],[3,4]]→[[1,0],[0,±2]]`）、`lll_command_rejects_non_matrix`

### 阶段 3 — C-9 完整 Buchmann 类群

- [ ] Buchmann 亚指数类群算法
- [ ] `class_number(P)` / `class_group(P)` / `is_principal(ideal)` 暴露
- [ ] 砍掉 `IDEAL_GEN_COORD_BOUND` / `N_J_MAX` sound-skip 边界
- [ ] 单元测试：`class_number(x^2+5)` = 2（虚二次 ℚ(√−5)）

### 阶段 4 — C-10 高次非主认证

- [ ] 完整 Step 7b（依赖 C-9）
- [ ] `certify_non_principal(P, ideal)`
- [ ] 单元测试：高次域非主理想证书

---

## 验证命令

```bash
cd giac-rs
cargo nextest run --release -p giac-core class_group -- --include-ignored   # 类群（sound-skip 现状）
cargo nextest run --release -p giac-core lattice                            # LLL
cargo nextest run --release -p giac-core alg_ext                            # AlgExt / AlgExtC
cargo nextest run --release -p giac-conformance giac_check_cas              # 含 rootof conformance
```

---

## 参考

- 主缺口索引：[GIAC-core-upstream-gaps](GIAC-core-upstream-gaps.md) §4 Phase E
- Hasse 数学正确性线：[GIAC-poly-f5-fglm-hasse-lean4-verification](GIAC-poly-f5-fglm-hasse-lean4-verification.md)
- API 分层：[giac-core-algebra-api-stability.md](../giac-core-algebra-api-stability.md)
- Pari/GP 对照：`bnfclassunit` / `bnfisprincipal` / `qflll` / `qfminim`
