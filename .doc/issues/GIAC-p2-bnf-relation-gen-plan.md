# Buchmann 关系生成对齐计划（giac-rs vs Pari upstream）

**状态:** open（方案 / 跟踪）
**类型:** #7 内禀 GRH 路径子计划
**父文档:** [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md)（R13 砖 7e ✅ 后）
**上游基线:** Pari `pari/src/basemath/buch2.c` — `small_norm` / `rnd_rel` / `RELAT` / `FBgen(LIMC)`
**Rust 落点:** `giac-rs/crates/giac-core/src/algebra/class_group.rs`（主）；`grh.rs`（7f）
**快照:** 2026-07-06

---

## 目标

在保持 **测试即规格 + sound-skip** 前提下，把关系生成与外层控制循环逐步对齐 Pari `buch2.c`，优先解锁 **deg≥3 h>1**（如 `x⁴−17` h=2）。**不追求行级移植**；绑数学语义，算法可自选并登记偏离。

---

## 现状快照（砖 7e 后）

```text
GRH deg≥3 主路径:
  factor_base_norm_bounds [M_K, LIMC]
  → ramification_relations
  → enumerate_relations_lli（确定性 ∏𝔭^a + LLL）
  → SNF → regulator_from_relation_arch → certify_hr_product
  → Bnf::try_from_grh_relations

仍缺:
  small_norm（deg≥3）
  rnd_rel（随机补关系）
  增量 RELAT（need=1，非整批重跑）
  7f GRHchk / LIMC 倍增 / goto START
  7h 复签名 arch
```

**7e 本质：** 外层对 `(因子基界 × exp_bound)` 做笛卡尔积重试，仍是 **批处理**，不是 Pari 的 **边加关系边验**。

---

## 优劣势分析

### 我们的优势

#### 数学与正确性

| 优势 | 说明 |
|------|------|
| **关系双证** | 每条关系要求 `γ ∈ J` 且 `|N(γ)| = N(J)`，比 Pari「先猜后 HNF 修」更严；错关系更难混入 SNF |
| **联合证书谓词** | `certify_hr_product(h', R')` 一次证关系格 + 单位格；与 Pari `compute_R` 数学同构（阈值 0.5–1.5，略宽于 Pari ~0.75–1.3） |
| **sound-skip 契约** | `None` = 不可证，非错答；符合项目「giac 有 bug 不复现」原则 |
| **deg-2 独立交叉证** | imag reduced forms / real 解析类数，不依赖 Buchmann 完备性也能给出 `h` |
| **SNF 幺模修复** | DIV-101 后 `UnimodMat` 类型化，关系格操作不易再破坏 gcd |

#### 工程与可维护

| 优势 | 说明 |
|------|------|
| **模块边界** | `grh.rs`（解析）、`bnf.rs`（arch/getfu）、`class_group.rs`（因子基+关系+SNF） |
| **可测性** | 每砖有锚定单测；改策略可回归（cubic/quartic h=1、`x⁴−17` probe） |
| **关系携带 γ** | `(Vec<BigInt>, γ)` 直接服务 `bnfisprincipal` 生成元重组 |
| **LLL 基础设施** | `ideal_from_prime_factors` + `lll_principal_generator` 已是 deg-agnostic 的 rnd_rel 核心 |
| **deg=2 模板** | `enumerate_relations_deg2` 可泛化为 deg≥3 `small_norm` |

#### 策略层面

| 优势 | 说明 |
|------|------|
| **不必绑死 buch2.c** | 只补数学缺口（增量 need + 第二关系源），不必搬 PRECI/HNF 一体状态机 |
| **确定性子集可保留** | `enumerate_relations_lli` 作可复现的「确定性 rnd_rel」，测试友好 |
| **7c fallback** | `getfu_analytic_fallback` 在 arch lift 失败时仍有 sound 退路（登记 known-divergence） |

### 我们的劣势

#### 关系生成（核心短板）

| 劣势 | 影响 |
|------|------|
| **仅一条源（deg≥3）** | 只有 `∏ 𝔭^a` + LLL；缺 `small_norm`（小系数 γ 扫范数） |
| **确定性枚举** | 固定指数网格；h>1 域上关键关系可能在网格外 |
| **批处理** | 一次枚举完再 SNF；失败则整批作废、外层加大 `exp_bound` |
| **无 rnd_rel** | 组合空间大时随机采样覆盖更好 |
| **LLL 天花板** | 只查归约基 `n` 个向量；生成元为非平凡 ℤ-组合时漏关系（sound 但不完备） |
| **`N_J_MAX = 10⁷`** | 大范数理想直接 skip |
| **非负指数 only** | 数学上 SNF 够用，但枚举集比 Pari 小 |

#### 外层控制循环

| 劣势 | 影响 |
|------|------|
| **无 GRHchk / LIMC2** | 因子基规模可能不足；大域卡死或 sound-skip |
| **无 goto START** | stuck 时不能动态倍增 LIMC |
| **无 PRECI** | arch/regulator f64 精度不足时无法像 Pari 提精度重算 |
| **7e 重试是笛卡尔积** | `[M_K, LIMC] × [4,6,8]` 有限档，非自适应 need |

#### 域类型与 API

| 劣势 | 影响 |
|------|------|
| **r₂>0 无 arch** | 7h 未做；含复嵌入的 GRH 全线 skip |
| **非极大序** | `power_order_is_maximal` 门控 |
| **`bnf_for_field` 无缓存** | 每次重算 `certified_class_data_with_gens` |
| **`bnfisprincipal` deg-2 only** | deg≥3 用户面未接 |
| **getfu 非精确** | f64 solve + fallback，非 `RgM_solve` |

### 对照总表

| 维度 | giac-rs | Pari | 对齐后目标 |
|------|---------|------|------------|
| 关系正确性 | 双证，严 | HNF 管线内修 | **保持双证** |
| 关系完备性 | 批处理+有限网格 | need 循环+双源 | 增量+双源 |
| 类群证书 | SNF + GRH 联合 | `compute_R` | 已有 ✅ |
| 单位来源 | 关系 arch + getfu | 同 | 已有 ✅（全实） |
| 因子基 | M_K+LIMC 两档 | GRHchk 动态 | 7f |
| 可复现测试 | 强（确定性枚举） | 弱（rnd） | seed 固定 rnd |
| 失败语义 | sound-skip | 重试至成功 | 保留 skip，减少 None |

---

## 对齐策略：不替换，叠加

**原则：** 保留 `enumerate_relations_lli` 作为 **确定性 rnd_rel**；叠加 **small_norm**、**rnd_rel**、**增量 RELAT**；因子基仍走 **7f**。

```text
关系源（概念层，可不抽 trait）:
  1. ramification_relations          — 已有
  2. enumerate_relations_small_norm    — 新（deg≥3，泛化 deg2）
  3. enumerate_relations_lli           — 保留
  4. rnd_rel_batch                     — 新

控制器（改 buchmann_grh_certified_data）:
  for norm_bound in factor_base_bounds:
    ideals = prime_ideals_below_norm_bound(...)
    relations = [ramification]
    loop need:
      relations += next_batch(source_i)
      if try_certify_relations(relations): return Some
      if budget_exhausted: break
    // 7f: goto START with larger LIMC
```

**不必为对齐而：**

- 删掉 LLL 路径改回纯坐标穷举（大域不可行）
- 强行支持负指数（数学上非必须，见 alignment 文档 R8 sound 边界）
- 逐行移植 `buch2.c` 状态机

---

## 分阶段计划

### 砖 7e-i：增量 RELAT + rnd_rel（优先） ✅

**状态：** giac-rs 已落地（`buchmann_grh_grow_relations` / `rnd_rel_one` / `try_certify_grh_relations`）。

**估时：** ~1–2 天 · **落点：** `class_group.rs`

**任务：**

1. 拆 `certified_class_data_grh_attempt`：
   - `try_certify_relations(field, ideals, relations) -> Option<CertClassData>`
   - `grow_relations_loop` 增量 append
2. `rnd_rel_batch(field, ideals, exp_bound, n_trials, seed)`：
   - 随机 `a ∈ [0, exp_bound]^k` → `ideal_from_prime_factors` → `lll_principal_generator`
   - 固定 seed 保证 CI 可复现
3. 控制逻辑：
   ```text
   relations = ramification
   for exp in BUCH_LLL_EXP_BOUNDS:
     relations += enumerate_relations_lli(..., exp)   // 可分批
     if try_certify: return
     relations += rnd_rel_batch(..., RND_TRIALS, seed)
     if try_certify: return
   ```

**验收：**

- `cargo test -p giac-core --lib` 全绿；clippy 绿
- `x⁴−17` h=2 probe：尽量 `Some`；仍 `None` 时 sound-skip 可接受直至 7e-ii
- 新测：`rnd_rel_*`（小域 + 固定 seed）

**风险：** 随机性需 seed；勿引入无界循环

---

### 砖 7e-ii：deg≥3 small_norm ✅

**状态：** giac-rs 已落地（`enumerate_relations_small_norm` + `relation_from_gamma_basis_factorization`；接入 `buchmann_grh_grow_relations`）。

**估时：** ~2–3 天 · **落点：** `class_group.rs`

**任务：**

1. `enumerate_relations_small_norm(field, ideals, coord_bound)`：
   - 泛化 `enumerate_relations_deg2`：`norm_of` + `compute_ideal_valuations_full` + `all_prime_factors_in_set`
   - 按 `Σ|c_i|` 或 `|N(γ)|` 递增；`coord_bound` 与 `exp_bound` 同步升
2. 生成顺序：`ramification → small_norm(B) → lli → rnd_rel`

**验收：**

- cubic/quartic h=1 回归
- `x⁴−17` 硬断言 h=2（若 GRH 路径稳定）
- 性能：deg=4、B≤16 单测有 timeout 守卫

**风险：** `(2B+1)^n` 爆炸；deg≥3 valuation sound-skip

---

### 砖 7f：GRHchk + LIMC 倍增

**估时：** ~2–4 天 · **落点：** `grh.rs` + `class_group.rs`

**任务：**

1. `grhchk(field) -> (limc, limc2, …)`（移植 Pari `GRHchk`/`GRHok` 参数）
2. `factor_base_norm_bounds` 扩为递增序列直至 `GRHok` 或上限
3. stuck → `goto START` 换更大 LIMC

**验收：**

- 大 |D| 不 hang（上限/timeout）
- `factor_base_norm_bounds_*` 测试扩域

---

### 并行后续（7h / 精确化 / API）

| 砖 | 内容 |
|----|------|
| 7h | `arch_log_of_element` 支持 r₂>0 |
| getfu | `RgM_solve` 或有理精确解 |
| API | `bnf` session 缓存；`bnfisprincipal` deg≥3 |

---

## 成功标准（里程碑）

| 里程碑 | 标准 |
|--------|------|
| **M1** | 7e-i 落地；`x⁴−17` probe 有概率变 `Some` |
| **M2** | deg≥3 h>1 ≥2 个域（含 `x⁴−17`）GRH `Some` 稳定 |
| **M3** | 7f 后大域不 infinite loop |
| **M4** | eval `bnfunits`/`bnfisunit` 在 cached `Bnf` 下 r≥2 稳定 |
| **M5** | [known-divergences.md](../known-divergences.md) 登记 fallback/threshold 与 Pari 差 |

---

## 明确不做（YAGNI）

- 逐行移植 `buch2.c` PRECI 全链（除非 M2 仍卡死）
- 负指数分数理想枚举
- 为对齐去掉 sound-skip
- giac C++ 侧 bnf（项目无原生 bnf）

---

## 建议执行顺序

```text
7e-i (增量 + rnd_rel)  →  探 x⁴−17
        ↓ 若仍 skip
7e-ii (small_norm)     →  再探
        ↓ 若因子基不够
7f (GRHchk/LIMC)       →  大域
        ↓ 并行
7h / RgM_solve / bnf 缓存
```

**第一刀最小 diff：** 只改 `buchmann_grh_certified_data` 控制流 + `rnd_rel_batch`（~150 行），不动 SNF/arch/getfu。

---

## 验证命令

```bash
cd giac-rs
cargo test -p giac-core --lib class_number_general_cert_grh
cargo test -p giac-core --lib factor_base_norm_bounds
cargo test -p giac-core --lib enumerate_relations
cargo clippy -p giac-core -- -D warnings
```

---

## 参考

- [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md) — R8–R13 砖序
- [known-divergences.md](../known-divergences.md) — 有意偏离登记
- Pari：`pari/src/basemath/buch2.c`（`small_norm` / `rnd_rel` / `Buchall_param`）
