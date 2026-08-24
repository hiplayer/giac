# GIAC `galoisconj4_main` 完整移植 — issue 跟踪

**状态:** open（P0–P2 ✅；P3c S₄/F₃₆ ✅；P3a ✅; P3b ✅; P5/P6/P9-upstream-shadow ✅; P7 ◐；**P9 G6 探针域 ✅**）
**类型:** AFK（除 P9 golden 脚本可 HITL 审 Pari 基线）  
**父项:** [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md) **R39r**（Galois `galoisconj` 矩阵，◐）  
**上游基线:** Pari `pari/src/basemath/galconj.c`（`galoisconj4_main` L2988）、`Zp.c`、`FpX.c`、`bibli2.c`、`base2.c`、`nffactor.c`  
**Rust 落点:** `giac-rs/crates/giac-core/src/algebra/galois_conj.rs` + 子模块 `galoisconj4/`（`analysis` `borne` `lift` `frobenius` `testlift` `fixed_field` `gen` `perm` `trace` `types`）、`padic/`（`zpx` `fpx_factor` `fpx_vandermonde` `types`）、`archimedean.rs`（`ArchBudget` / 复根单路径）  
**快照:** 2026-08-24 · giac-rs **`2884573`**

**说明：** `galoisconj4` 是 **p-adic Frobenius 提升 + 模 Vandermonde**，**不是** `buch2.c` 的 LLL。`galoisconj_in_field` 已对齐 Pari `galoisconj_monic`（g4→g1）；`EmbAutPerms` 在 `galoisinit` 失败时仍保留 n≤8 arch 启发式（BNF 专用，非 conjugates 主路径）。

---

## 总览（依赖 DAG）

```text
P0a ──┬──► P1a ──► P2 ──► P3a ──► P3c ──► P4 ──► P5 ──► P6
P0b ──┤         ▲              ▲
P0c ──┘         │              │
P0d ────────────┴──► P1b ──────┘
P7（nfroots）可与 P3 尾期并行，P5 前需完成
P8（pr_orbit）依赖 P5 接线
P9（golden）依赖 P5；P6/P7 完成后扩金值矩阵
```

| ID | 标题 | 状态 | 阻塞 | 估时 | Pari 锚点 |
|----|------|------|------|------|-----------|
| **P0a** | `ZpX_liftroot` + `ZpX_roots` | **✅** | — | 3–4d | `Zp.c` L755–829 |
| **P0b** | `FpX_roots` + 按度因子计数 | **✅** | — | 2d | `FpX_factor.c` |
| **P0c** | `FpV_invVandermonde` | **✅** | P0b | 2d | `FpX.c` L1865 |
| **P0d** | `ZX_disc_all` + `indexpartial` + 复根精度 | ✅ | — | 2–3d | `polarit3.c` / `base2.c` / `bibli2.c` |
| **P1a** | `galoisanalysis` + `numberofconjugates` | ✅ | — | 3d | `galconj.c` L1084–1235, L3061 |
| **P1b** | `galoisborne` + `initgaloisborne` | ✅ | — | 2d | `galconj.c` L245–281 |
| **P2** | Frobenius 提升链 | **✅** | — | 8–10d | `galconj.c` L290–612, L1967–2110 |
| **P3a** | `galoisgen` 循环 + 固定域 | **◐** | — | 5d | `galconj.c` L2196–2222, L2781+ |
| **P3b** | `galoisgenlift` / 幂零扩张 | **✅** | P3a | 4d | `galconj.c` L2252+, L2703+ |
| **P3c** | A₄ / S₄ / F₃₆ 快路 | **✅** | P2 | 2d | `galconj.c` L2794–2822 |
| **P4** | `permtopol` + `galoisvecpermtopol` | **✅** | P5 | 3d | `galconj.c` `vectopol` 族 |
| **P5** | `galoisconj4_main` 编排 + `GaloisConjugates` 接线 | **✅** | P3a,P3b,P3c,P4 | 3d | `galconj.c` L2988–3057 |
| **P6** | `GaloisInit` 缓存 + `GaloisAutPerms` 同源 | **✅** | P5 | 2d | `galoisinit` |
| **P7** | `galoisconj1` / `nfroots` 回退 | **◐** | P0a,P1a | 5–6d | `galconj.c` L37–63; `nffactor.c` |
| **P8** | `pr_orbit_fill` → `be_honest`（`pr_orbit_fill` ✅） | open | P5 | 1d | `buch2.c` `be_honest`（`pr_orbit_fill` ✅） |
| **P9** | Pari golden harness + 文档 | **◐** | P5 | 2d | conformance · **upstream-shadow ✅** |

**合计（串行上界）：** ~10–14 人周；P0/P7 可并行减日历时间。

---

## galoisconj 测试验收规格（P4–P9 共用）

**目的：** 避免把「管线阶段单测绿」或「MonicZx 中间层计数」误当成 `galoisconj` 已验收。本文与 [test-writing-spec.md](../test-writing-spec.md) §1（A/B/C）、[rust-migration-plan.md](../rust-migration-plan.md) §6.4（算法偏离）互补；**giac 53 golden 不覆盖** `giac-core` 内 `galoisconj`（见 P9）。

### 权威层级

```text
数学真值（ℚ-自同构、共轭为根、群复合）
        │
Pari 登记基线（R39r：nfgaloisconj / galoisinit / numberofconjugates）
        │
giac-core field 层谓词（try_insert_conjugate 等，见下）
        │
管线阶段 B 测（perm_to_pol、Frobenius 阶、borne 数量级…）
        │
✗ 不可作唯一依据：PowerBasisElement 个数、perm 候选数、emb 启发式命中
```

**绑定：** 数学语义 + 本表谓词 + Pari 金值（P9）；**不绑定** giac C++ 行级实现，也不绑定 `galconj.c` 内部临时结构。

### 返回值语义（写测试前必读）

| 概念 | 规格 |
|------|------|
| `galoisconj(nf)` / `GaloisConjugates` | `K=ℚ(α)` 中 **σ(α)** 的列表（`HighFirstQ`，幂基 high-first） |
| 列表长度 | `numberofconjugates(T)`（**非**一般等于 `deg T`；非 Galois 时 `< deg`） |
| 恒等自同构 | `σ(α)=α` → `field.generator_coords()`，**排在最后** |
| **不是** 平凡 Galois | 域元素整数 `1`（`PowerBasisElement` low-first `[1,0,…]` ≠ σ(α)） |
| MonicZx 层 `galoisconj1` / `perm_to_pol` | **中间表示**；验收必须在 field 层经 `power_basis_to_high_first` + `try_insert_conjugate` |

坐标：`PowerBasisElement` = Pari low-first `{1,α,…}`；`ExtensionField::generator_coords()` = high-first，α 在 `v[n-2]`。见 `field_arith::generator_coords` 与 `ExtensionField::generator_coords` 注释，**勿混用**。

### Field 层验收谓词（A′ 集成 — 唯一关单依据）

**入口：** `galoisconj_in_field` / `GaloisConjugates::compute`（不是裸 `galoisconj1(t)`）。

对每个域 `K=ℚ(α)`、首一不可约 `T`、嵌入 `emb`（`n≤8` easy 路径需 `emb.reliable`；`n>8` 回退路径不强制）：

| ID | 谓词 | 实现锚点 |
|----|------|----------|
| **G1** | `len(conjs) == expected_conjugate_count(T)` | `number_of_conjugates` |
| **G2** | 每个 `σ`：`m(σ(α))=0`（在 `K` 内） | `sigma_preserves_minpoly` / `try_insert_conjugate` |
| **G3** | 每个 `σ`：`σ` 为 `K/ℚ` 自同构（在生成元上） | `is_field_automorphism_on_generator` |
| **G4** | `conjs.identity() == field.generator_coords()` | `reorder_conjugates_identity_last` |
| **G5** | （Galois 且 `len>1`）`σ ∘ τ` 与列表中某元素一致 | `galoisapply` / 矩阵复合（P6+ 推荐补测） |
| **G6** | （P9）与 Pari `nfgaloisconj(nf)` **逐项**幂基坐标相等 | golden harness |

**关单最低集：** G1–G4 全过；G5 对 cyclic / WSS 表内域；G6 由 P9 覆盖。

### 分层单测：什么能证什么

| 层级 | 测谁 | 可断言 | 不能替代 |
|------|------|--------|----------|
| **A′** | `galoisconj_in_field` | G1–G5 | — |
| **B** | `number_of_conjugates`、`perm_to_pol`、`galois_analysis` | 与 Pari debug 一致的数量级/阶/素数 | G1–G4 |
| **B** | `galoisconj1(t)` 仅 MonicZx | `nf_roots` 分裂、Frobenius 置换阶 | **禁止** 只断言 `galoisconj1.len() >= c` |
| **C** | trace / 快照 | 诊断 | 语义 |

### 探针域矩阵（crate 内单测最低）

| 域 | deg | `numberofconjugates` | 主路径 | 必验 |
|----|-----|----------------------|--------|------|
| ℚ(i) | 2 | 2 | easy / g4 | G1–G4 |
| x³−3x+1 | 3 | 3 | g4 | G1–G4，G5 阶 3 |
| ℚ(∛11) | 3 | **1** | `expected==1` 快路径 | G1,G4（仅 α） |
| Φ₁₁ | 10 | 10 | g4 | G1–G4 |
| x¹⁰−10x+1 | 10 | **1** | `expected==1`（n>8） | G1,G4；证明 P7 不挡路 |
| x⁴+1 等 WSS | 4 | 4 | g4 / nilp | G1–G5 |

`c>1` 且 `n>8` 的非 Galois 域（`1 < numberofconjugates < deg`）在 P9 金值表登记前：**`#[ignore]` 语义测 + issue**，禁止用弱计数让 CI 绿（[test-writing-spec.md](../test-writing-spec.md) §6 双轨）。

### 写测试前检查单（AFK / PR）

```text
□ 返回值是 σ(α) 还是别的？恒等在 field 层长什么样？
□ 个数用 numberofconjugates 还是 Gal 阶？
□ 断言在 galoisconj_in_field 还是 MonicZx 中间层？
□ Pari 对照多项式是否已登记？无则 ignore + 本 issue 或 P9 表
□ 是否覆盖 G2–G3（不只 G1 个数）？
□ 性能：split-prime / Frobenius 扫描是否有预算上界？
```

### 反模式（P7 教训登记）

| ❌ | ✅ |
|----|-----|
| `PowerBasisElement [1,0,…]` 表示平凡 Gal | `expected==1` → `generator_coords()` |
| `galoisconj1` 候选个数 ≥ c 即过关 | `galoisconj1_to_field_conjugates` 后 G1–G3 |
| 入口强制 `emb.reliable` 挡掉 n>8 回退 | easy 需 emb；g4/g1 不依赖 |
| 扫 5000 素数无早停当实现细节 | `find_totally_split_prime` 预算 + 失败 `None` 可登记 |
| 管线全绿 = galoisconj 验收 | P9 / G6 或探针域 G1–G4 |

**偏离：** 与 Pari 行为不同但数学正确 → [known-divergences.md](../known-divergences.md)；暂不可证 → sound-skip + `#[ignore]`。

---

## P0a — `ZpX_liftroot` + `ZpX_roots` ✅（2026-07-10）

**模块：** `giac-core/src/algebra/padic/zpx.rs`（`padic.rs` → `padic/mod.rs`）

**落地：**
- `zpx_liftroot` / `zpx_liftroot_typed` → `ZpRootLift`
- `zpx_roots` — `FpPolynomial::split_part` + `FpRoots` + Hensel
- `FpPolynomial::split_part` — BigInt Frobenius `x^p mod f`

**余量（分阶段，非 ponytail）：**
- 无 `ZpX_liftfact` / `ZpX_liftroots_full`（重因子、部分分裂需 MultiLift）
- `x³−3x+1` 验收用 **p=17**（mod 7 不可分，与 Pari 一致返回空）

**验收：**

- [x] `(x²-2) mod 7^e`，`e=1..5`，根金值一致
- [x] 分裂三次 `x³−3x+1 mod 17` 三根提升（`e=1..3`）
- [x] `cargo test -p giac-core zpx` 5 绿；`hensel_lift_factor` 回归绿

---

## P0b — `FpX_roots` + 按度因子计数 ✅（2026-07-10）

**模块：** `padic/fpx_factor.rs`

**落地：**
- `fp_x_factor` → `FpFactorization`；`fp_x_nbfact_by_degree` → `DistinctDegreeCounts`
- `FpIrreducibleFactor::roots_in_fp` — 线性 + `FpX_quad_root`（`fp_sqrt` / Legendre）
- `zpx_roots` 经 `fp_x_roots_typed`

**余量：** `p: i64`（`giac_poly` 边界）；未走 Pari `ddf_Shoup` 快路

**验收：**

- [x] `x²−2 mod 7`：`D[1]=2`；`x³−3x+1 mod 7` 不可约 `D[3]=1`；mod 17 分裂 `D[1]=3`
- [x] `fp_x_roots` 与 P0a 金值一致

---

## P0c — `FpV_invVandermonde` ✅（2026-07-10）

**模块：** `padic/fpx_vandermonde.rs`

**落地：** `fpv_inv_vandermonde_typed` → `InvVandermonde`（`FpProductTree` 内部）

**验收：**

- [x] `n=3,5,10`：`V·M ≡ I (mod p)`
- [x] 可选 `den` 标量乘子

---

## P0d — `ZX_disc_all` + `indexpartial` + 复根精度 ✅（2026-07-10）

**模块：** `padic/zx.rs`；`galoisconj4/borne.rs`；`archimedean.rs`（`ArchBudget` / `arch_budget`）

**做什么：** 非 monogenic 域的 `den`（`indexpartial`）；`galoisborne` 所需 `embed_roots` / `vandermondeinverse` 精度界。复根嵌入走 **单一路径**：Sturm 实根 → `polish_root_rat_until` → Laguerre+deflation（Pari `roots_aux` / `clean_roots` 对齐的 `ArchBudget`）；已移除 companion 特征值 / 线程超时双路径。

**验收：**

- [x] `ℚ(√5)` 非极大：`indexpartial` → `den=2`（与 Pari `nfdisc` 一致）
- [x] `galoisborne` 对 `x³-3x+1` 输出 `valabs`/`ladicabs` 与 Pari debug 同级数量级（f64 ±1）
- [x] `archimedean` embeddings 单测 23 绿（含 `ℚ(√2)`、三次一实二复、A₄/S₄ arch 嵌入）

**Blocked by:** 无

---

## P1a — `galoisanalysis` + `numberofconjugates` ✅（2026-07-10）

**模块：** `galoisconj4/analysis.rs`

**做什么：** 扫素数得 Frobenius 阶、totally split `l`、WSS 判定；失败返回 `None`（非 Galois / 非 WSS）。`numberofconjugates` 供 `galoisconj1` 预算。

**验收：**

- [x] `ℚ(∛11)` → `None`（非 Galois）
- [x] `Φ₁₁` → 通过，`l` 与 Pari 一致
- [x] `ℚ(i)`、`x³-3x+1` 通过

**诊断：** `GIAC_GA_TRACE=1` 或测试内 `set_ga_trace_for_test(true)` 打印 Frobenius / `calcul_l` 环进度。

**Blocked by:** P0b

---

## P1b — `galoisborne` + `initgaloisborne` ✅（2026-07-10）

**模块：** `galoisconj4/borne.rs` + `types.rs`（`GaloisBorne`）

**做什么：** 复根嵌入、`matrixnorm`、`logint` 得 `valsol`/`valabs`/`ladicsol`/`ladicabs`/`bornesol`。

**验收：**

- [x] 结构体字段与 Pari `struct galois_borne` 语义一一对应
- [x] `Φ₅`、`Φ₁₁` 精度界单测（与 `gp` 打印对照或快照）

**Blocked by:** P0d

---

## 类型栈（2026-07-10 复审 · `b4c6c42`）

`padic/types.rs` + `galoisconj4/types.rs` — 算法中间结果均有命名类型，非裸 `Vec<BigInt>`：

| 类型 | 数学对象 | 用于 |
|------|----------|------|
| `FpModulus` | 𝔽_p | 所有 mod-p 运算 |
| `FpPolynomial` | monic f ∈ 𝔽_p[x] | split_part、因子、Vandermonde 树 |
| `FpIrreducibleFactor` | 不可约因子 | `.roots_in_fp()` |
| `FpFactorization` | 平方因子分解 | `fp_x_factor` |
| `DistinctDegreeCounts` | D[d] 按度计数 | `galoisanalysis`（P1a） |
| `FpRoots` | f 在 𝔽_p 中的根集 | `fp_x_roots_typed` |
| `PadicPrecision` / `ZpRootLift` | 根 mod p^e | `zpx_liftroot_typed` |
| `FpProductTree` | 乘积树 | Vandermonde 内部 |
| `InvVandermonde` | V⁻¹ 的 Lagrange 列 | `galoisconj4` permtopol |
| `ZqPolynomial` | f ∈ (ℤ/Qℤ)[x] | `GaloisLift.t_mod_q`、Bezout |
| `ZqQuotientElement` | (ℤ/Qℤ)[x]/(T) 元素 | `pauto` 幂、`frobeniusliftall` |
| `BezoutLiftFactors` | Bezout 提升 cofactor | `init_test_lift` |
| `AutomorphismPowers` | x, aut, aut∘aut, … (`FpXQ_autpowers`) | `fpxq_autpowers` |
| `ComboProductCache` | `C` / `Cd` 组合积缓存 | `frobeniusliftall` |
| `PermTestMatrix` | Vandermonde 测试行 | `galois_test_perm` |
| `GaloisAnalysis` | Frobenius 扫描 / WSS 判定 | `galois_analysis` |
| `GaloisBorne` | p-adic + archimedean 系数界 | `galois_borne` |
| `GaloisLift` / `GaloisTestLift` / `GaloisPermTest` | Pari lift 状态机 | P2 全链 |
| `NewtonSumMatrix` / `SymmetricPolynomial` / `SympolOrbitValues` | 固定域 sympol | `fixed_field_sympol` |
| `OrbitImageValues` | sympol 轨道像 `PL`（**非**共轭根） | `galois_gen_fixed_field0` |
| `FactorImageIndex` | `get_image` → `gf->psi[g]`（**1-based**） | `galois_gen_lift` |
| `PsiCofactorDegrees` | Pari `gf->psi`（slot `0` 空） | `galois_frobenius_lift` |
| `PermTestPvOrderTable` | `td->order[n]` → `Vmatrix` 行 | `testpermutation` |
| `FixedFieldOrbits` / `FixedFieldPrep` | 轨道 + 固定多项式 prep | `galois_gen_fixed_field0` |
| `ComplexEmbeddings` / `VandermondePrep` | 复根 + T′(α_i) 积 | `init_galois_borne` |
| `ArchBudget` | Sturm 深度 / 残差 / Laguerre 容差 | `archimedean` 复根 |

**已消除的 ponytail（galoisconj4 / padic / archimedean）：**

- ~~`fp_x_roots_brute`（p < 10⁵ 穷举）~~ → `FpX_quad_root`（Legendre + Tonelli–Shanks `fp_sqrt`）
- ~~`p.to_u64()` Frobenius~~ → `poly_x_pow_mod_f(&BigInt)` 全精度
- ~~`fp_x_normalize` 静默 `lc_inv=1`~~ → `FpPolynomial::from_zx_monic` 仅在可逆时缩放
- ~~`poly_divmod_p` 首项为 0 时死循环~~ → 每步 `trim` + 零首项退出（`p=11` ramified 触发）
- ~~`GaloisTestLift` 裸 `Vec<Vec<BigInt>>`~~ → `BezoutLiftFactors` / `AutomorphismPowers` / `ComboProductCache` / `PermTestMatrix`
- ~~`fixed_field_sympol` 简化单射判定~~ → Pari `sympol_is1to1` + `vecsmall_is1to1`（`perm.c`）
- ~~`archimedean` companion 特征值 + 250ms 线程超时双路径~~ → Sturm→Laguerre+deflation 单路径 + `ArchBudget`

**登记余量（非 ponytail，为 Pari 管线分阶段）：**

| 余量 | 说明 | 关闭于 |
|------|------|--------|
| `FpFactorization` 经 `giac_poly::factor_mod_irreducibles` | 要求 `p: i64`；非 Pari `ddf_Shoup` 快路，但因子/计数数学正确 | 可选 P1 前 BigInt 模因子分解 |
| `zpx_roots` 无 `ZpX_liftfact` | 重因子/部分分裂时与 Pari 路径不等价 | P0a 余量或 P5 |
| `legendre_symbol` 非二次剩余 → 空根集 | 正确语义，非捷径 |

---

## P2 — Frobenius 提升链 ✅（2026-07-10 · `b4c6c42`）

**模块：** `galoisconj4/lift.rs`、`frobenius.rs`、`testlift.rs`；`padic/zpx.rs` 增 `zpxq_lift_monomorphism`

**落地：**
- `GaloisLift` / `PadicRootEmbedding` / `RootPermutation` / `FrobeniusLift` / `FrobeniusFind`
- `init_lift`, `galois_do_lift`, `galois_do_lift_n`, `zpxq_lift_monomorphism`（线性 Hensel）
- `init_test_lift`（`BezoutLiftFactors` + `AutomorphismPowers` + `ComboProductCache`）
- `frobenius_lift_all`, `pol_to_perm_test`, `galois_frobenius_test`, `galois_test_perm`
- `galois_frobenius_lift_nilp`, `galois_find_frobenius` / `galois_find_frobenius_full`
- `GIAC_G4_TRACE=1` / `GIAC_GA_TRACE=1` 阶段计时（`trace.rs`）

**余量：** `monoratlift` 早停；WSS 非循环大域端到端仍依赖 P3b

**验收：**

- [x] `x³-3x+1`：Frobenius 置换阶 3
- [x] `Φ₁₁`：非平凡 Frobenius 提升，置换长度 10
- [x] `galois_test_perm` 对已知 S₃ 子群 membership 正确（`x³−3x+1`：A₃ 通过、对换拒绝；Frobenius 置换通过）
- [x] `cargo test -p giac-core --release --lib galoisconj4` **91/91**

**Blocked by:** — · **可启动 P3b / P5**

---

## P3a — `galoisgen` 循环 + 固定域 ✅（2026-07-10）

**模块：** `galoisconj4/gen.rs`、`fixed_field.rs`、`lift_gen.rs`

**落地：**
- `galois_gen` 主入口：Frobenius 阶 = n → `galois_gen_cyclic`；否则 → `galois_gen_fixed_field`
- `fixed_field_orbits` / `fixed_field_sympol`（Pari `sympol_is1to1` + `vecsmall_is1to1` + `fixedfieldsurmer`）
- `NewtonSumMatrix` / `SympolOrbitValues` / `sympol_eval` / `galois_gen_fixed_field0`（sympol + `sigma` + `t_mod_factors`）
- `get_image` / `galois_gen_fixed_field_rec` / `galois_gen_lift` + `galois_gen_lift_auto`（`testpermutation` intheadlong 主路径，对齐 Pari `galconj.c`）
- `flxq_minpoly`（`fixed_field_fact_mod`）；`fixed_poly_mod_p`（`get_image` 的 `P mod p`）
- 循环群：`cyclic_group_elts` + `perm_to_pol` 生成共轭列

**余量：** `galoisgenlift_nilp` / `is_central`；`galoisconj4_main` 接线；golden harness

**已移除 ponytail（2026-07-13）：** deg-2 `get_image` 兜底、`brute_orbit_perm_lift`、psi 多索引重试

**验收：**

- [x] 循环三次域：群阶 = 3，生成元置换正确（`galois_gen_cubic_cyclic_order_3`）
- [x] `Φ₁₁`：ℤ/10 循环群生成（`galois_gen_phi11_cyclic_order_10`）
- [x] `x⁴+1`：`galois_gen_fixed_field0` sympol + 固定多项式 prep ~2s（`galois_gen_fixed_field0_x4_plus_1_prep`）
- [x] 非循环 WSS 域 `x⁴+1` 完整 `galois_gen` 群阶 4（`galois_gen_x4_plus_1_order_4`）
- [x] `cargo test -p giac-core --release --lib galoisconj4` **91/91**

**Blocked by:** — · **可启动 P5**

---

## P3b — `galoisgenlift_nilp` / 中心扩张 ✅（2026-07-13）

**模块：** `galoisconj4/nilp.rs`、`lift_gen.rs`、`gen.rs`、`fixed_field.rs`

**落地：**
- `PcPresentation` / `pcgrp_lift` / `pcgrp_insert` / `pc_to_perm` / `genorbit` / `nilp_froblift`
- `PariVec1`（1-based `t_VECSMALL`）/ `PcWord` 固定 Pari 索引语义
- `galois_gen_lift_nilp_full` 主循环；`GaloisPermCache`（`permtoaut` + `pc_evalcache`）
- `FixedFieldGaloisStep` 扩展 PG2–PG5；度 2 固定域快捷路径始终带完整 nilp 元数据
- `is_central` 时 `bad ∨ dis` 传播；轨道搜索预算 `NILP_ORBIT_MAX`；失败回退 `galoisgenliftauto`
- Frobenius easy 分支去重（避免 `perm_inv(frob)==frob` 误命中）

**验收：**

- [x] 度 4 非循环 WSS（`x⁴+1`）群阶 4（P3a：`testpermutation` 对齐 Pari）
- [x] `galois_gen_lift_nilp_x4_plus_1_order_4` — 强制 nilp 路径 + `galoisgenliftauto` 回退
- [x] `is_central_extension` 谓词单测；`cargo test -p giac-core --release --lib galoisconj4` **91/91**
- [x] Pari 语义对齐：`PariVec1`/`PcWord`、`pcgrp_insert`、`brl_add`、`pc_exp`、`id_factor`、`search_pc.br` 临时更新
- [x] 度 >104 真实中心扩张多项式（`ga_easy` 关闭；polcyclo(107) 度 106，`#[ignore]` ~33s；analysis 探针 `galoisconj_analysis_degree_106_gt_104` 绿）
- [x] `genorbit`：`k=1`+H 分裂 → `None`（Pari 同条件）；`k≥2` 单测；见 `known-divergences.md` DIV-103

**Blocked by:** — · **P3b 主线已关闭**（nilp 模块完整实现；`#[ignore]` 仅性能相关）

---

## P3c — A₄ / S₄ / F₃₆ 快路 ✅（2026-07-14 · F₃₆ 2026-08-03）

**模块：** `galoisconj4/specials.rs`；`perm.rs` 增 `vec_perm_orbits`；`padic/zx.rs` `indexpartial`

**落地：**
- `a4_galois_gen` + `try_special_galois_gen` — Pari `a4galoisgen`（deg 12, `ga.ord==3`, `!p4`）
- Phase1/2 因子枚举退出条件对齐 Pari（`a≠0` 时停）
- `indexpartial` / `ZpX_reduced_resultant_fast` 对齐 → A₄ den 与 Pari 一致
- `galois_borne`：f64 Vandermonde 范数比 Pari REAL 略紧时 `valsol+=1`（见下表余量）
- **`galoisanalysis` `norm_o`**：按 Pari `Fpe=p^e` 累乘（勿从 `o` 起乘裸素数）→ deg-24 与 Pari `p=73,ord=6,deg=3`
- **`padic_root_embedding`**：`L` 用 `valabs`，`Lden` 仍 mod `ladicsol`（对齐 Pari `makeLden`）
- **`roots_to_monic_poly`**：`(x−r)·f`（曾误用 `(1−rx)` 导致系数反转）
- **`sympol_aut_evalmod`**：`f∘σ`（Pari `FpX_FpXQ_eval`）；`Sp` 禁止 `from_zx_monic`（会破坏域元素）→ `FpPolynomial::from_zx`
- **`vectopol` / `PowerBasisElement`**：Pari `gdiv(·,den)` → `Ratio` 系数；`get_image` / Fp 路径用 `to_fp`（`RgX_to_FpX`）
- **嵌套 fixed-field：** `Pgb.l = gb->l`；PL 升/降精度；autos 回映到 PL 序（`gens_on_pl_roots`）→ `x⁴+1` get_image / nilp 已绿
- **`testpermutation` 不变量**（见下节 + `lift_gen.rs` / `testlift.rs` 注释）→ deg-24 WSS golden 绿
- **`FpX_ffisom` / `FpXQ_ffisom_inv` / `FpXV_ffisom`**（`padic/fpx_ffisom.rs`）：`deg(P)|deg(Q)` 嵌入/同构
- **`FpX_ffintersect`**（`padic/fpx_ffintersect.rs`）：Allombert 全路径 — special / cyclo / Hilbert-90
- **`FpXQ_sqrtn`**（`padic/fpxq_sqrtn.rs`）：`gen_Shanks_sqrtn` on `𝔽_q^*` — cyclo 分支无 root fallback
- **`FpXV_chinese`**（`padic/fpx_chinese.rs`）+ **`mkliftpow`**（`lift.rs`）：S₄ `liftp` 列 = CRT(`trans(misom)`) → `automorphismlift`
- **`s4releveauto` / `lincomb` / `s4makelift` / `s4test`** + **`FqC_FqV_mul`**（`specials.rs`）
- **`s4_galois_gen` / `try_s4`** — prep + σ→τ→φ 搜索接线（`galconj.c` L1542–1665）；验收 `s4_galois_gen_orders_24_p4` + **`galoisconj_golden_s4_degree_24`**（Pari `s4galoisgen` 探针，orders `[2,2,3,2]`，~15s release）
- **`bezout_lift_fact`**：按因子 `hensel_lift_factor` + CRT 幂等元（非裸 cofactor）；**`FpXQ_autpowers`**：合成幂 `[x,σ,σ∘σ,…]`（曾误为乘法幂 `[1,σ,σ·σ,…]`）
- **`f36releveauto2/4` / `f36_galois_gen` / `try_f36`** — prep + σ→τ→ρ 搜索（`galconj.c` L1699–1876）；验收 `f36_galois_gen_orders_36_p4` + **`galoisconj_golden_f36_degree_36`**（Pari `f36galoisgen` 探针，orders `[3,3,4]`，~80s release）
- 余量：`FpX_nbroots` 已对齐 Pari `Flx_nbroots`（Frobenius gcd，非全因子分解）；`calcul_l` 与 upstream 同语义
- 余量：`bezout_lift_fact` 无 Pari `MultiLift` 产品树（逐因子 Hensel；S₄ n=6 可接受，大 g 可换树）

### S₄ `s4galoisgen` 数学结构（2026-08-03）

探针 `s4_galois_probe_poly()`：`Gal= S₄`，`identify=[24,12]`，PC orders `[2,2,3,2]`。

| 阶段 | 数学对象 | 代码 |
|------|----------|------|
| 门控 | `ord=3`（8 个三次 Frobenius）且 `p4≠0`（6 个四次因子 mod p4）⇒ S₄ 非 A₄ | `S4GaloisCandidate` |
| Prep | `T mod p4` → 6 个四次不可约 ⇒ `FpXV_ffisom` → `mkliftpow` + bezout 幂等元 | `s4_galois_gen` 前半 |
| σ | 三对因子 `(a,b)` 的 `s4releveauto` + `(ℤ/4ℤ)³` 上 `lincomb` → 3-阶置换 | σ 搜索 |
| τ, φ | `rot3` 因子标号 + 类似搜索 | τ、φ |
| 输出 | `r1=στ`, `r2=φστφ`, `r3=φσ`, `r4=σ` | `special_galois_gen_result(...,[2,2,3,2])` |

**Rust 探针见证（`galois_analysis_s4_probe_matches_pari`）：** `plift=97`, `p4=109`（扫素数 `p0=53`，首个 ord-4），Frobenius 首次 split `l=83`，`calcul_l` 后 `l=643`。Pari 文献常写 `@29/@31` 为更小素数上的同型分解，非 `galoisanalysis` 扫序首个见证。

### F₃₆ `f36galoisgen` 数学结构（2026-08-03）

探针 `f36_galois_probe_poly()`：`Gal= F₃₆`（`SmallGroup(36,9)`），`identify=[36,9]`，PC orders `[3,3,4]`。

| 阶段 | 数学对象 | 代码 |
|------|----------|------|
| 门控 | `ord=3`（12 个三次 Frobenius）且 `p4≠0`（9 个四次因子 mod p4）⇒ F₃₆ 非 A₄/S₄ | `F36GaloisCandidate` |
| Prep | `T mod p4` → 9 个四次不可约 ⇒ `FpXV_ffisom` → `mkliftpow` + bezout 幂等元 | `f36_galois_gen` 前半 |
| σ | 四对因子 `(a,b)` 的 `f36releveauto2` + `(ℤ/4ℤ)⁴` 上 `lincomb` → 3-阶置换 | σ 搜索 |
| τ | `f36releveauto4` + `rot3` 因子标号 + `(ℤ/4ℤ)²` 搜索 | τ |
| ρ | 固定 `w[4][6]` 表 + `rot4(sp[3],sp[5],sp[8],sp[7])` 重标号 + `(ℤ/4ℤ)³` 搜索 | ρ |
| 输出 | `r1[τ[i]]=ρ[i]`, `r2[i]=τ⁻¹[ρ[i]]`, `r3=τ` | `special_galois_gen_result(...,[3,3,4])` |

**Rust 探针见证（`f36_galois_gen_prep_diagnostic`）：** `plift=181`, `p4=73`, `l=701`, `ord=deg=3`；9 个四次因子 @ p4。ρ 阶段 `rot4` 须对齐 Pari 1-based `(3,5,8,7)` → 0-based `(2,4,7,6)`（非 `(2,4,7,5)`）。

### `testpermutation` / `galoisgenliftauto` 移植不变量（2026-07-24）

对标 Pari `galconj.c` `testpermutation` / `Vmatrix`。破坏后典型症状：`nn` 穷尽、`headlong_ok≈1`、`galois_gen_lift` 失败（V₄/`x⁴+1` 可能仍绿）。

| # | 不变量 | 说明 |
|---|--------|------|
| 1 | **`ar` 全量重算** | 每轮 `ar[a+1]=0` 后重算 `ar[a]…ar[1]`。**有意不跟** Pari 增量后缀（`ar_from`）；数学等价，去掉循环携带状态。勿再引入半截增量除非对照 Pari 中间量回归 |
| 2 | **`G[cx]=F[orbit_id]`** | `gel(G,cx)=gel(F,coeff(B,i,j))`；`B` 是 Frobenius **轨道编号**。禁止 `find(cycle.contains(root))` |
| 3 | **`umael(W,a,b)=L[b]·y[a]`** | Pari 列主序；禁止转置。对合可能掩盖错误 |
| 4 | **验收** | deg-24 **S₄ `s4galoisgen`**（`galoisconj_golden_s4_degree_24`，orders `[2,2,3,2]`）+ WSS（`galoisconj_golden_s4_wss_degree_24`） |

**代码锚点：** `lift_gen.rs`（`test_permutation`、`f_cycle_for_coeff`）；`testlift.rs`（`vmatrix_headlong`）；`types.rs`（`HeadlongPvMatrix`）。

**验收：**
- [x] `vec_perm_orbits` 单测
- [x] A₄ 快路代码 + 接线（`try_special_galois_gen`）
- [x] `a4_galois_gen_order_12` 端到端
- [x] `galoisconj_golden_a4_degree_12`（`galoisconj4_main`）
- [x] `x⁴+1` get_image / `galois_gen` / nilp
- [x] deg-24 WSS e2e（`galoisconj_golden_s4_wss_degree_24`）
- [x] deg-24 S₄ `s4galoisgen` e2e（`galoisconj_golden_s4_degree_24`，PC orders `[2,2,3,2]`）
- [x] `FpXV_ffisom` + S₄/F₃₆ 门控接线
- [x] `FpX_ffintersect` 真嵌入（`deg(P)|deg(Q)`；special + cyclo + Hilbert-90）
- [x] `FpXV_chinese` + `mkliftpow`（`mkliftpow_x4_plus_1_mod_5`）
- [x] `s4releveauto` / `lincomb` / `s4makelift` / `s4test`（`s4_make_lift_and_test_frobenius_cubic`）
- [x] `try_s4` σ→τ→φ（`s4_galois_gen_orders_24_p4`，Pari [24,12] 探针）
- [x] F₃₆ 快路端到端（`f36_galois_gen` / `try_f36`，`f36_galois_gen_orders_36_p4` + `galoisconj_golden_f36_degree_36`）

**登记余量：**

| 余量 | 说明 | 关闭于 |
|------|------|--------|
| `valsol += 1`（f64） | Pari 用 REAL/`ceil_safe`；`x⁴+1` 在精确 den=4 时 f64 少 1 个 `l`-digit | 多精度 arch 范数 |
| `zpx_roots` 排序 | ≠ Pari `galoisinit` 根序 → Pari sigma 单测 ignore | 可选根序对齐 |
| S₄ / F₃₆ 快路 | `try_s4` / `try_f36` ✅ | — |
| `bezout_lift_fact` | 逐因子 Hensel，无 `MultiLift` 树 | 大 `g` 时换产品树 |
| `galois_analysis` `plift`/`p4` | 见证素数依赖扫序起点（`p0≈2n`）；与 Pari 同循环但首个 `@p4` 可差（如 109 vs 31） | 可选：对齐 `improves()` 序 |

**Arch（2026-07-14，对标 Pari `QX_complex_roots` / `fujiwara_bound`）：**
- Sturm 隔离界：`min(Cauchy, 2·Fujiwara)`（宽 Cauchy 会丢大根）
- verify：相对残差 `|P|/Σ|a||x|ⁱ`；`|x|>1` 在倒数多项式上算
- 偶多项式：`y=x²` 代换；负实 `y` 走复平方根（勿 clamp 到 0）

**Blocked by:** —（`FpXV_ffisom` 已用于 S₄/F₃₆）

---

## P9 — Pari golden harness + 文档 ◐（2026-07-14）

**落地：**
- `scripts/galoisconj_golden.sh` — Rust `galoisconj_golden_counts` + 可选 `PARI_GOLDEN=1` gp 对照
- `main.rs` `golden`：ℚ(i)、x³−3x+1、Φ₁₁、x⁴+1、∛11、**A₄ deg-12**、**S₄ `s4galoisgen` deg-24**、**WSS deg-24**、**F₃₆ `f36galoisgen` deg-36**、**deg-106 中心扩张 analysis**）

**验收：**
- [x] `./scripts/galoisconj_golden.sh` 绿（探针域）
- [x] deg 12 golden 绿
- [x] `embeddings_s4_degree_24_all_real`（arch）
- [x] deg 24 WSS golden 绿（`galois_gen_lift` / `testpermutation` 对齐 Pari）
- [x] deg 36 golden 绿（`galoisconj_golden_f36_degree_36`，orders `[3,3,4]`；**default nextest `#[ignore]` ~80s**）
- [x] G6 逐项坐标 = Pari `nfgaloisconj`（探针域：ℚ(i)、x³−3x+1、∛11、x⁴+1、Φ₁₁；`g6_pari_nfelt` + `permtopol_e2e`）

**落地（2026-08-03）：**
- `galoisconj4/pari_golden.rs` — `assert_g6_nfelt_multiset`（low-first 整数 nfelt，multiset）
- `galois_conj.rs` `g6_pari_nfelt` — field 层 G6；**x³−3x+1 high-first = `[1,0,-3,1]`**（勿与 MonicZx low-first `[1,-3,0,1]` 混用）
- `main.rs` `permtopol_e2e` — MonicZx 层 Q(i)/三次/Φ₁₁ nfelt multiset（尾部零归一化）
- `scripts/galoisconj_golden.sh` — 追加 `g6_pari_nfelt`、`upstream_shadow`

**Blocked by:** P5（field 层 G1–G4 已绿；G6 探针域已绿；非 WSS / deg>8 扩金值待 P6/P7）

## P4 — `permtopol` + `galoisvecpermtopol` ✅（2026-07-10 · Pari nfelt 2026-08-03）

**模块：** `galoisconj4/perm.rs`

**落地：** `vec_permute` / `vec_to_pol` / `perm_to_pol` / `perm_cycles` / `cyclic_group_elts` / `root_perms`；`galois_gen_cyclic` 已用 `perm_to_pol` 生成共轭；**`galois_vec_perm_to_pol`** + `PermToPolPrep::galois_vec_perm_to_pol` 批量 API。

**余量：** ~~arch easy 主路径~~ 已移除（P5-upstream）；conjugates 仅 g4/g1

**验收：**

- [x] `perm_to_pol` 三次非平凡 Frobenius 单测（mod 17 Vandermonde）
- [x] `perm_cycles` / `cyclic_group_elts` 阶 3 单测
- [x] `galois_vec_perm_to_pol` 批量 = 逐 `perm_to_pol`；= `galois_init.conjugates`（Q(i)/Φ₁₁/三次）
- [x] `ℚ(i)`：共轭 `α ↦ −α`（`galoisconj4_main` + Pari nfelt multiset）
- [x] `Φ₁₁`：10 个共轭 `m(σ(α))=0` 全过 + Pari nfelt multiset
- [x] 与 Pari `galoisconj(nf)` 逐项 nfelt 等价（Q(i)、Φ₁₁ 登记基线；G6 扩金值见 P9）

**Blocked by:** —

---

## P5-upstream — `galoisconj_monic` 路由对齐 ✅（2026-08-03）

**做什么：** 删除 n≤8 arch easy 主路径；`galoisconj_in_field` ≡ Pari `galoisconj_monic`（deg 快捷 → g4 → g1）。

**落地：**
- `galoisconj_monic_g4_g1` — 显式 g4→g1 核心
- 删除 `galoisconj_easy` / `sigma_alpha_from_arch_slot_perm` 等 conjugates 用 arch 启发式
- `galoisconj_upstream_routing_shadow` — 探针域 shadow 测试
- `EmbAutPerms::compute_from_arch_heuristic` 保留（BNF 专用）

**验收：**
- [x] 探针域 G1–G4 + G6 仍绿
- [x] shadow 测试

**Blocked by:** — · Phase B = P6 GaloisInit 缓存 ✅

---

## P9-upstream-shadow — CI 回归门禁 + 文档 ✅（2026-08-03）

**做什么：** Phase C — 防 arch easy 回退；entry / snapshot / public API 同源 shadow；文档与 R39r 状态更新。

**落地：**
- `galoisconj_upstream_only` — Pari `galoisconj_monic` 参考路径（deg 快捷 → g4 → g1）
- `upstream_shadow` 模块（4 测）：routing / snapshot / public API / g4_g1 核心
- `arch_emb_heuristic` — `galois_root_perms_large_n_*` 迁出并标注 BNF 专用
- `scripts/galoisconj_golden.sh` — 追加 `upstream_shadow` 必跑
- `.doc/giac-galoisconj4-pari-port.md` — Pari↔Rust 函数对照
- `known-divergences.md` DIV-105 — arch easy 删除登记

**验收：**
- [x] `./scripts/galoisconj_golden.sh` 含 shadow 绿
- [x] 5 探针域 entry ≡ upstream-only ≡ snapshot ≡ `GaloisConjugates`
- [ ] 非 WSS G6 扩表（P7 后）；`upstream-shadow` 4 测已绿
- [ ] R39r 父文档矩阵全 ✅（g4 主线已绿，P8 余量仍 open；父文档已标注「余量：P8」）

**Blocked by:** —

---

## P5 — `galoisconj4_main` 编排 + `GaloisConjugates` 接线 ✅（2026-08-03）

**模块：** `galoisconj4/mod.rs`；改 `galois_conj.rs`

**做什么：** 完整 `galoisconj4_main` 流水线；`GaloisConjugates::compute` 对齐 Pari `galoisconj_monic`（g4→g1，无 arch easy）。

**落地（2026-08-03）：**
- `galoisconj_in_field_inner` 对齐 Pari `galoisconj_monic`：deg-1/2 快捷 → `galois_init`/`galoisconj4_main` → `galoisconj1`；**G1 严格 `len == numberofconjugates`**
- **P5-upstream：** 删除 n≤8 arch easy；`galoisconj_monic_g4_g1` + shadow 测试
- field 层探针矩阵单测 `galoisconj_probe_*_g1_g4`：ℚ(i)、x³−3x+1、ℚ(∛11)、Φ₁₁、x⁴+1

**验收：**

- [x] 探针域 G1–G4 + G6
- [x] R39a grow 回归

**Blocked by:** —

---

## P6 — `GaloisInit` 缓存 + `GaloisAutPerms` 同源 ✅（2026-08-03）

**模块：** `galois_conj.rs` `FieldGaloisSnapshot`；`class_group.rs` `GrhRelCache`

**落地：**
- `FieldGaloisSnapshot::compute` — 至多一次 `galois_init`；`init()` + `conjugates()` 同源
- `GaloisConjugates` / `AutomorphismMatrices` / `GaloisAutPerms::compute` 改读 snapshot
- `GrhRelCache.galois_snapshot` + `GaloisAutPerms::from_snapshot`

**验收：** [x] `field_galois_snapshot_shared_init_and_conjugates`（cubic，`init` Some）；[x] `field_galois_snapshot_deg2_has_no_init`（Q(i) 短路无 `galoisinit`）；[x] grow/R39r 回归绿

**余量：** g1-only 域（无 `galoisinit`）`EmbAutPerms` 仍走 arch 启发式

---

## P7 — `galoisconj1` / `nfroots` 回退 ◐（2026-07-13）

**模块：** `galoisconj4/galoisconj1.rs`；`analysis::find_totally_split_prime`；`galois_conj.rs` 路由

**做什么：** `galoisconj4_main` 失败 → `galoisconj1`（`numberofconjugates` + `nfroots`）；`nfsqff` ROOTS 模式。

**落地（◐）：**

- `nf_roots` / `galoisconj1`；Frobenius 置换 + `perm_to_pol`；`galoisconj1_to_field_conjugates`
- 路由：`g4 → g1`（P5-upstream ✅）；`expected==1` → `[α]` 快路径（不调用 `nfroots`）
- `n>8` 不强制 `emb.reliable`

**验收（见上文 [galoisconj 测试验收规格](#galoisconj-测试验收规格p4p9-共用)）：**

- [x] ℚ(∛11)：`len=1`，G4（`galoisconj_q_cbrt11`）
- [x] x¹⁰−10x+1：`n>8`，`len=1`（`galoisconj_degree10_non_galois_single_conjugate`）
- [x] x³−3x+1：`galoisconj1` 三次 Galois `len=3`（MonicZx B 测）
- [ ] 非 WSS、`1 < c < n`：探针域 G1–G4（待 P9 多项式）；`galoisconj1` 模块已就绪，`nf_roots` + `perm_to_pol_prep` 单测绿）
- [ ] 与 Pari `nfgaloisconj` 失败回退一致（sound-skip → `known-divergences.md`）

**余量：** `find_totally_split_prime` 扫 `GA_PRIME_SCAN_CAP` 无早停；`galoisconj1` 按候选 pbe 计数非 field 验证后计数；`c>1,n>8` 非 Galois 未端到端证明。

**Blocked by:** P0a, P1a（可与 P3 并行）· **P9 扩金值后标 ✅**

---

## P8 — `pr_orbit_fill` → `be_honest`（2026-08-24）

**模块：** `class_group.rs` / `bnf.rs`

**做什么：** grow `be_honest` 用已有 `pr_orbit_fill` 跳过非正规素理想轨道（Pari `buch2.c`）。


**落地：**
- `pr_orbit_fill` 已实现并单测通过（`galois_conj.rs:61`）
- `be_honest` 高层函数尚未实现

**验收：**

- [x] `pr_orbit_fill_q_i_split_p5_marks_both_ideals` 绿
- [ ] `be_honest` 实现 + grow 接线
- [ ] grow 路径在 ℚ(i) 分裂素上不误拒关系（集成测或登记探针域）

**Blocked by:** P6 ✅（Galois 缓存已落地）

---

## P9 — Pari golden harness + 文档

**做什么：** `scripts/galoisconj_golden.sh`（或 Rust integration）批跑 Pari vs giac；更新 [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md) R39r → ✅；**[giac-galoisconj4-pari-port.md](../giac-galoisconj4-pari-port.md)** 函数对照表；修正全文「galoisconj4=LLL」表述。

**状态（2026-07-13）：** 见上文 **P9 ◐** 节（脚本 + 探针域 G1/G4 已绿；deg 12/24/36 `#[ignore]`）。

**金值矩阵（2026-08-24 验证）：**

| 域 | deg | Gal |
|----|-----|-----|
| ℚ(i) | 2 | ℤ/2 |
| x³−3x+1 | 3 | ℤ/3 |
| ℚ(∛11) | 3 | 平凡 |
| Φ₅ | 4 | ℤ/4 |
| x⁴−17 | 4 | ℤ/2 |
| Φ₁₁ | 10 | ℤ/10 |
| deg-12 A₄ | 12 | A₄ |
| deg-24 S₄ | 24 | S₄ |

**验收：**

- [x] 上表探针域 golden 绿 + **G6**（ℚ(i)、x³−3x+1、∛11、x⁴+1、Φ₁₁）；deg 12/24/36 计数已绿；F₃₆ `#[ignore]` ~80s
- [x] **upstream-shadow** CI 门禁（4 测：`routing_matches_upstream_only`、`g4_g1_core_matches_upstream_only`、`snapshot_matches_entry_and_upstream`、`public_api_matches_snapshot`）
- [x] [giac-galoisconj4-pari-port.md](../giac-galoisconj4-pari-port.md) 函数对照
- [x] 与 [galoisconj 测试验收规格](#galoisconj-测试验收规格p4p9-共用) G1–G4 在探针域上一致
- [ ] 非 WSS / `1<c<n` 扩 G6 金值（P7 后）
- [ ] R39r 父文档状态更新（父文档已标注「余量：P8」；R39r 主线 ✅）

**Blocked by:** P6/P7（扩金值矩阵）

---

## 建议 PR / 分支命名

```text
r39r-g4-p0a-zpx-roots
r39r-g4-p0c-fpv-vandermonde
r39r-g4-p1-analysis
r39r-g4-p2-frobenius
r39r-g4-p3-galoisgen
r39r-g4-p5-main-wire
r39r-g4-p6-galois-init-cache
r39r-g4-p9-upstream-shadow
r39r-g4-p7-nfroots
r39r-g4-p9-golden
```

---

## 验证命令（全程）

```bash
cd giac-rs
cargo test -p giac-core --release --lib galoisconj4   # 91 绿 + 7 ignored（P2–P3c ✅；P3b ✅；P7 ◐）
cargo test -p giac-core --release --lib galois_conj # 28 绿（G1–G4 探针 + G6 golden + upstream-shadow）
cargo test -p giac-core --release --lib archimedean::tests  # 23 绿（P0d 复根 + A₄/S₄ arch 嵌入）
cargo test -p giac-core --lib zpx      # P0a 起
cargo test -p giac-core --lib galois   # P5 起全量
cargo test -p giac-core --lib

# 诊断 trace（固定域 prep / Frobenius / analysis）
GIAC_G4_TRACE=1 GIAC_GA_TRACE=1 cargo test -p giac-core --release --lib galois_gen_fixed_field0_x4_plus_1_prep -- --nocapture
```
