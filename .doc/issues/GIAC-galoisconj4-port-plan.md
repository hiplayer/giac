# GIAC `galoisconj4_main` 完整移植 — issue 跟踪

**状态:** open（P0–P2 ✅；P3a/P4 ◐；**下一步 P3b `galoisgenlift`**）
**类型:** AFK（除 P9 golden 脚本可 HITL 审 Pari 基线）  
**父项:** [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md) **R39r**（Galois `galoisconj` 矩阵，◐）  
**上游基线:** Pari `pari/src/basemath/galconj.c`（`galoisconj4_main` L2988）、`Zp.c`、`FpX.c`、`bibli2.c`、`base2.c`、`nffactor.c`  
**Rust 落点:** `giac-rs/crates/giac-core/src/algebra/galois_conj.rs` + 子模块 `galoisconj4/`（`analysis` `borne` `lift` `frobenius` `testlift` `fixed_field` `gen` `perm` `trace` `types`）、`padic/`（`zpx` `fpx_factor` `fpx_vandermonde` `types`）、`archimedean.rs`（`ArchBudget` / 复根单路径）  
**快照:** 2026-07-10 · giac-rs **`b4c6c42`**

**说明：** `galoisconj4` 是 **p-adic Frobenius 提升 + 模 Vandermonde**，**不是** `buch2.c` 的 LLL。现有 `galois_conj.rs` arch 槽 / 根置换枚举保留为 `galoisconj_easy`（deg≤8 快路），完整移植后由本计划逐项替换余量。

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
| **P3b** | `galoisgenlift` / 幂零扩张 | open | P3a | 4d | `galconj.c` L2252+, L2703+ |
| **P3c** | A₄ / S₄ / F₃₆ 快路 | open | P2 | 2d | `galconj.c` L2794–2822 |
| **P4** | `permtopol` + `galoisvecpermtopol` | **◐** | P5 | 3d | `galconj.c` `vectopol` 族 |
| **P5** | `galoisconj4_main` 编排 + `GaloisConjugates` 接线 | open | P3a,P3b,P3c,P4 | 3d | `galconj.c` L2988–3057 |
| **P6** | `GaloisInit`（flag=1）+ `GaloisAutPerms` 改读 | open | P5 | 2d | `galoisinit` |
| **P7** | `galoisconj1` / `nfroots` 回退 | open | P0a,P1a | 5–6d | `galconj.c` L37–63; `nffactor.c` |
| **P8** | `pr_orbit_fill` → `be_honest` | open | P5 | 1d | `buch2.c` `be_honest` |
| **P9** | Pari golden harness + 文档 | open | P5 | 2d | conformance |

**合计（串行上界）：** ~10–14 人周；P0/P7 可并行减日历时间。

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
- [x] `archimedean` embeddings 单测 15 绿（含 `ℚ(√2)`、三次一实二复）

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
| `AutomorphismPowers` | 1, aut, aut², … | `fpxq_autpowers` |
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
- [x] `cargo test -p giac-core --release --lib galoisconj4` **27/27**

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
- [x] `cargo test -p giac-core --release --lib galoisconj4` **33/33**

**Blocked by:** — · **可启动 P3c / P5**

---

## P3b — `galoisgenlift_nilp` / 中心扩张 ◐

**模块：** `galoisconj4/lift_gen.rs`、`gen.rs`

**做什么：** Pari `galoisgenlift_nilp`（`is_central` 分支）、`fixedfieldinclusion`；大域组合搜索。

**验收：**

- [x] 度 4 非循环 WSS（`x⁴+1`）群阶 4（P3a：`testpermutation` 对齐 Pari）
- [ ] 幂零/中心扩张多项式
- [ ] 失败路径返回 `None` 不 panic（大 n 组合爆炸）

**Blocked by:** —

---

## P3c — A₄ / S₄ / F₃₆ 快路

**模块：** `galoisconj4/specials.rs`

**做什么：** `a4galoisgen`、`s4galoisgen`、`f36galoisgen`（deg 12 / 24 / 36）。

**验收：**

- [ ] 各度至少 1 个 Pari 登记多项式：`nfgaloisconj` 个数 = 群阶
- [ ] 与通用 `galoisgen` 结果一致（交叉验）

**Blocked by:** P2（可与 P3a 并行，P5 前需齐）

---

## P4 — `permtopol` + `galoisvecpermtopol` ◐（2026-07-10 · `b4c6c42`）

**模块：** `galoisconj4/perm.rs`

**落地：** `vec_permute` / `vec_to_pol` / `perm_to_pol` / `perm_cycles` / `cyclic_group_elts` / `root_perms`；`galois_gen_cyclic` 已用 `perm_to_pol` 生成共轭。

**余量：** 未替换 `galois_conj.rs` 中 `sigma_alpha_from_arch_slot_perm` 主路径；无 `galoisvecpermtopol` 批量 API

**验收：**

- [x] `perm_to_pol` 三次非平凡 Frobenius 单测（mod 17 Vandermonde）
- [x] `perm_cycles` / `cyclic_group_elts` 阶 3 单测
- [ ] `ℚ(i)`：共轭 `α ↦ -α`（经 `galoisconj4_main` 端到端）
- [ ] `Φ₁₁`：10 个共轭，`m(σ(α))=0` 全过
- [ ] 与 Pari `galoisconj(nf)` 逐项 `nfelt` 等价

**Blocked by:** P5（模块内 API 已就绪）

---

## P5 — `galoisconj4_main` 编排 + `GaloisConjugates` 接线

**模块：** `galoisconj4/mod.rs`；改 `galois_conj.rs`

**做什么：** 完整 `galoisconj4_main` 流水线；`GaloisConjugates::compute`：先 `galoisconj_easy`（n≤8），不足则 `galoisconj4_main`；删除 `GALOIS_CONJ_ENUM_MAX_N` 硬顶对大 n 的退化。

**验收：**

- [ ] `cargo test -p giac-core --lib galois` 全绿且 case 数增加
- [ ] `Φ₁₁`：10 共轭经 **galoisconj4** 路径（非仅 arch 启发式）
- [ ] `x³-11` 平凡 Gal 仍 1 非平凡共轭
- [ ] R39a grow 回归：`add_relation_galois_orbit` ℚ(i) / ∛11 仍绿

**Blocked by:** P3a, P3b, P3c, P4

---

## P6 — `GaloisInit`（flag=1）+ `GaloisAutPerms` 改读

**模块：** `galoisconj4/types.rs`；`galois_conj.rs`；`class_group.rs`

**做什么：** 组装 Pari `galoisinit` 对象（`pol/T/L/M/den/group/cyc`）；`GrhRelCache.galois` 一次 `galoisconj4_main(flag=1)`；`FbAutPerms`/`EmbAutPerms` 从 `L/M/ladic` 派生，去掉重复矩阵重建。

**验收：**

- [ ] `GaloisAutPerms::compute` 与 P5 共轭列一致
- [ ] `ideal_perm_under_galois_q_i` 等 R39r 单测无回归
- [ ] `EmbAutPerms` 与 `FbAutPerms` 同源（文档登记若仍双路径）

**Blocked by:** P5

---

## P7 — `galoisconj1` / `nfroots` 回退

**模块：** `nf_roots.rs` 或扩 `number_field_arith.rs`

**做什么：** `galoisconj4_main` 失败 → `galoisconj1`（`numberofconjugates` + `nfroots`）；`nfsqff` ROOTS 模式。

**验收：**

- [ ] 非 WSS 域：`galoisanalysis` 失败后走 `nfroots`，个数 = `numberofconjugates`
- [ ] 与 Pari `nfgaloisconj` 失败回退行为一致（sound-skip 登记见 `known-divergences.md`）

**Blocked by:** P0a, P1a（可与 P3 并行）

---

## P8 — `pr_orbit_fill` → `be_honest`

**模块：** `class_group.rs` / `bnf.rs`

**做什么：** grow `be_honest` 用已有 `pr_orbit_fill` 跳过非正规素理想轨道（Pari `buch2.c`）。

**验收：**

- [ ] 单测已有 `pr_orbit_fill_q_i_split_p5_marks_both_ideals` 仍绿
- [ ] grow 路径在 ℚ(i) 分裂素上不误拒关系（集成测或登记探针域）

**Blocked by:** P5（Galois 缓存稳定后接线）

---

## P9 — Pari golden harness + 文档

**做什么：** `scripts/galoisconj_golden.sh`（或 Rust integration）批跑 Pari vs giac；更新 [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md) R39r → ✅；新增 `.doc/giac-galoisconj4-pari-port.md` 函数对照表；修正全文「galoisconj4=LLL」表述。

**金值矩阵（最低）：**

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

- [ ] 上表全部 golden 绿（或 sound-skip 登记）
- [ ] `cargo test -p giac-core --lib` 全绿
- [ ] R39r 父文档状态更新

**Blocked by:** P5（P6/P7 完成后扩表）

---

## 建议 PR / 分支命名

```text
r39r-g4-p0a-zpx-roots
r39r-g4-p0c-fpv-vandermonde
r39r-g4-p1-analysis
r39r-g4-p2-frobenius
r39r-g4-p3-galoisgen
r39r-g4-p5-main-wire
r39r-g4-p7-nfroots
r39r-g4-p9-golden
```

---

## 验证命令（全程）

```bash
cd giac-rs
cargo test -p giac-core --release --lib galoisconj4   # 27 绿（P2–P3a）
cargo test -p giac-core --release --lib archimedean::tests  # 15 绿（P0d 复根）
cargo test -p giac-core --lib zpx      # P0a 起
cargo test -p giac-core --lib galois   # P5 起全量
cargo test -p giac-core --lib

# 诊断 trace（固定域 prep / Frobenius / analysis）
GIAC_G4_TRACE=1 GIAC_GA_TRACE=1 cargo test -p giac-core --release --lib galois_gen_fixed_field0_x4_plus_1_prep -- --nocapture
```
