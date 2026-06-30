# GIAC — Hasse 平方判定：完整 Scholz 设计与实施路线

**状态:** 跟踪文档（Step 1-7a 已实现；Step 7a = `(C)` 直接主性（生成元搜索）+ Scholz 管线 `hasse_is_square`；剩余 = Step 5 复域 torsion/log-map + r≥2 index-1 证书 + Step 7b 完整类群（Buchmann，证非主性 `Some(false)`）+ `hasse_is_square` 接入 `sqrt_fmodule`）。完整 Scholz ≈ 重写 PARI `bnfinit`，按 7 步推进，每步独立可测、sound-false-pass-on-skip（与 (A) 同架构）。
**关联:** [GIAC-poly-f5-fglm-complete-square-decision](issues/GIAC-poly-f5-fglm-complete-square-decision.md)（总 issue / Phase C 优先级）、[giac-poly-f5-fglm-ideal-valuation-proof](giac-poly-f5-fglm-ideal-valuation-proof.md)（(A) 赋值公式证明）、`giac-rs/crates/giac-core/src/algebra/number_field_arith.rs` + `poly_roots.rs::sqrt_fmodule`
**快照:** 2026-06-30

---

## 0. 数学约束（完整 Scholz）

`u ∈ (K×)² ⟺` 以下四条**同时**成立：

| 条件 | 含义 | 状态 |
|---|---|---|
| **(A_fin)** | 所有有限素理想 `v_𝔭(u)` 为偶 | 部分（good/tame/wild 已落地 + dedekind 分支素数 bug 已修；缺大素因子 + `p∣index`(Montes) + tower） |
| **(A_inf)** | u **全正**：每个实嵌入 `σ(u) > 0` | ✅ Step 1（`archimedean::is_totally_positive`） |
| **(C)** | `J = ∏𝔭^{v_𝔭/2}` 在 `Cl(K)` 中主 | ✅ Step 7a（`ideal_is_principal` 直接主性：生成元搜索 `|N(γ)|=N(J)` + 赋值成员校验；小范数 `J` 全过，大范数/非主 `None` sound-skip。缺 Step 7b 完整类群 → 证非主性 `Some(false)` + 大范数 LLL 短向量） |
| **(B)** | `η = u/γ² ∈ (𝔬_K×)²`（γ 是 J 的生成元） | ✅ Step 6（`unit_is_square`：Dirichlet 分解 + mod-2 指数 + torsion 平方；certified 域类型 r=0 虚二次 / r=1 全实二次 全过；r≥2/复域 sound-false-pass） |

> **关键修正**：[总 issue](issues/GIAC-poly-f5-fglm-complete-square-decision.md) 原 P4 范围（"LLL + 基本单位 + mod-2"）**漏了 (C) 类群主性 + (A_inf) 全正**。完整 Scholz ≈ 重写 PARI `bnfinit`。本路线按 7 步推进，每步独立可测、sound-false-pass-on-skip（与 (A) 同架构：任一 skip → 续走 `sqrt_base_case`，ℚ 验证兜底，永不 false-reject）。

---

## 1. 接口审计（需标注的歧义 / 缺口）

| 接口 | 文件:行 | 问题 / 标注 | 状态 |
|---|---|---|---|
| `ideal_valuations_parity_ok → bool` | `number_field_arith.rs:98` | 丢弃 `(𝔭_i, v_𝔭_i(u))`、`g_i`、Hensel `G_i`、`u_int`、`d_u` —— (B)/(C) 全需要。**Step 2 重构为 `IdealValuationScan` 结构体** | 待 Step 2 |
| `generator_minpoly_low() -> Option<LowFirstQ>` | `ext_tower.rs:500` | `pub(crate)`，tower 返 `None`。**标注：仅单层 K=ℚ(α)；`None` = skip（sound false-pass），非 error；Step 2/3 key off 此 `None` bail tower** | ✅ 已标注 |
| `FieldEmbedding` | `ext_tower.rs:1044` | **命名冲突**：这是子域 ℚ-线性嵌入，**非** archimedean `σ_i:K→ℂ`。后者尚未实现（Step 3，`archimedean.rs`）；doc 注明勿混淆 | ✅ 已标注 |
| `field_arith::generator_coords(n)` vs `ExtensionField::generator_coords()` | `field_arith.rs:162` / `ext_tower.rs:466` | 同名不同物不同序：free fn 返 low-first `CoordsQ` 幂基单位向量（建乘阵用）；method 返 high-first `HighFirstQ` 字段生成元 α。两侧 doc 互指 + 显式转换路径 `LowFirstQ::from_high`/`HighFirstQ::new` | ✅ 已标注 |
| `krylov_minpoly_coords -> CoordsQ` | `poly_roots.rs:1103` | 输入 high-first `HighFirstQ`，输出 low-first monic `CoordsQ`；无标签 `CoordsQ` 隐藏方向翻转（foot-gun）。doc 注明 + 显式转换路径 | ✅ 已标注 |
| `clear_denoms_low` (private) | `number_field_arith.rs:371` | (B) 需 `u_int = d_u·u ∈ ℤ[α]`；私有，Step 2 一并暴露 | 待 Step 2 |
| **缺失** | — | 高精度 log / 类群 —— Step 5–7 缺（archimedean 嵌入 Step 3、LLL Step 4 已落地） | Step 5–7 |

---

## Step 1 —— signature (r₁, r₂) + 全正 (A_inf)【✅ 已实现】

**位置**：`number_field_arith.rs` 旁新子模块 `algebra/archimedean.rs`（注册 `pub(crate) mod archimedean;` 于 `algebra/mod.rs`）。纯 ℚ 算术，**不依赖** LLL/f64/嵌入/类群。

### 函数

**`field_signature(field: &Arc<ExtensionField>) -> Option<(usize, usize)>`**
- tower（`generator_minpoly_low() = None`）→ `None`（调用方 skip，sound）。
- `m_α` → `poly_from_low` → `giac_poly::Poly`。Cauchy 根界 `B = 1 + max|a_i|/|lc|`（新增 ~10 行）。
- `r1 = sturmab_count_rational_poly(m_α_poly, &x, &(-B), &B)?`（已存在，`poly_alg_sturm.rs:89`）。
- `r2 = (d_k - r1) / 2`。

**`is_totally_positive(field: &Arc<ExtensionField>, u_high: &HighFirstQ) -> bool`**
- tower → `true`（skip，sound）。
- `(r1, r2) = field_signature(field)`；`r1 == 0`（复域）→ `true`（复嵌入恒平方，无全正约束）。
- 否则：Sturm **二分隔离** m_α 的 r1 个实根到互不相交区间（~80 行，纯 `Ratio<BigInt>` 递归用 `sturmab_count_rational_poly` 计数二分）；每个区间细化到 `u_poly`（u 作 ℚ[t]，`LowFirstQ::from_high`）在两端点同号且非零（精确有理 Horner）；任一实根处 `sign(u(ρ)) ≤ 0` → 返 `false`。

### 接线（`poly_roots.rs:~2055`，紧接 `ideal_valuations_parity_ok` 之后、`if d_f == d_k` 之前）

```rust
if !super::super::number_field_arith::archimedean::is_totally_positive(field, u_coords) {
    return None;
}
```

**Soundness**：`u = δ²` 的每个实嵌入 `σ(δ²) = σ(δ)² ≥ 0`；若 `σ(δ) = 0` 则 `σ(u) = 0` → `N(u) = 0`，已被早先 `norm.is_zero` 拦截。故 `δ²` 全正 ⟹ **永不误拒真平方**。

### 测试（PARI `/home/kanli.hu/upstream/pari/gp` 交叉验证，serial `--test-threads=1`）

- **signature**：`ℚ(√2)` → (2,0)；`m=t³+t+1`（1 实根）→ (1,1)；`m=t⁴+t+1`（0 实根, 2 复对）→ (0,2)。PARI `nf.r1`/`nf.r2`。
- **全正拒**：`ℚ(√2)`，`u=-1` → N=1（norm 过）、(A_fin) 过（无素数整除 1），但 `σ₁(-1)=-1<0` → (A_inf) **拒**。PARI `nfeltissquare(bnfinit(t^2-2), -1)=0`。**严格强于 norm+(A_fin)**。
- **不误拒**：`u=(1+√2)²=3+2√2`（全正真平方）→ 过；`u=2`（`σ₁(2)=σ₂(2)=2>0`）→ 过。
- **复域**：`m=t³+t+1`，`u=-1` → 无实嵌入 → 全正返 `true`（不约束）。

### 验证
`cargo test -p giac-core --lib algebra::number_field_arith -- --test-threads=1` + clippy。

---

## 步骤 2–7（路线图）

### Step 2 —— 重构 (A) 为 `IdealValuationScan`【✅ 已实现】

- **做什么**：`ideal_valuations_parity_ok` 拆成 `compute_ideal_valuations(field, u_high, norm) -> IdealValuationScan`，返回结构体：
  ```
  IdealValuationScan {
      u_int_low: Vec<BigInt>, d_u: BigInt,
      signature: (usize, usize), totally_positive: bool,
      entries: Vec<IdealValEntry>,  // { p, g_i: PolyMod, e_i, f_i, v_p_u: i64, status: Processed|Skipped }
      parity_ok: bool,
  }
  ```
  `ideal_valuations_parity_ok` 保留为薄包装（向后兼容）。
- **约束**：skip 语义保留（sound false-pass）；Skipped 素数标记 —— (C)/(B) 仅在所有素数 Processed 时才有意义，否则 skip。
- **实现偏离**：`signature`/`totally_positive` **未并入** scan —— 留在 `archimedean`（Step 1）解耦，避免 `sqrt_fmodule` 热路径双重计算 (A_inf)。完整 pre-(B)/(C) bundle 在 Step 5+ 组装。`all_processed` 含义标注为"所有*被检*素数 Processed"（`v_p(N)=0` 素数未检 —— Step 7 (C) 完备性前的已知限制，已写入 struct doc）。
- **文件**：`number_field_arith.rs`。**工作量**：中。**依赖**：Step 1。**PARI**：`nfeltval` 逐素数。

### Step 3 —— 数值嵌入（f64）【✅ 已实现】

- **做什么**：`archimedean::embeddings(field) -> Option<Embeddings>`（r1 实根 + 2·r2 复共轭，共 d_k 个 `σ_i(α)` as `nalgebra::Complex<f64>`）。实根用 Step 1 的 Sturm 二分隔离 → `RealRootCtx::real_roots_f64` 精化到 f64 端点重合（自然 f64 精度极限）；**复根用 companion matrix + `nalgebra::complex_eigenvalues`**（`nalgebra` 升为 `giac-core` 直接依赖，~30 行）。`eval_at_embedding(u_high, σ_α) -> Complex<f64>`：u_low 的 Horner。
- **约束**：f64 ~15 位精度；d≤12 + 有界系数下根精度足够判符号 + 粗 log。**verify-and-skip**：`|m_α(σ_α)| > ROOT_RESID_TOL(1e-6)` 或 `real+complex ≠ d_k` → `Embeddings.reliable=false`，下游 skip（sound）。复根过滤 `|im|>COMPLEX_IM_TOL(1e-7)` —— 实根权威归 Sturm，companion 对实根的数值噪声被丢弃；genuinely-near-real 复根被误滤会触发 count 不匹配 → reliable=false → skip（sound，已写 `ponytail:` 上限+升级路径）。
- **接口形态**：`Embeddings { sig, roots, reliable }` —— 把"共轭 + 是否可信"打包，避免调用方各自重算 sig 或漏检 reliable。`embeddings` 对 d_k≤1/tower 返 `None`（ℚ 无生成元可嵌，归 base case）。`eval_at_embedding` 接 `HighFirstQ`（与 (A_inf)/(A_fin) 入口一致），内部 `from_high` 转 low-first Horner。
- **文件**：`archimedean.rs`（+ `giac-core/Cargo.toml` 增 `nalgebra`）。**工作量**：中。**依赖**：Step 1。**PARI**：`polroots(m_α)` 对照（√2 / t³+t+1 / t⁴+t+1 三例全过，tol 1e-7~1e-9）。
- **测试**：5 个新增 —— `embeddings_{q_sqrt2,cubic,quartic}_match_pari`（实/混合/全复三类签名）、`embeddings_eval_at_embedding_recovers_generator`（u=α ⟹ σ(u)=σ(α)）、`embeddings_eval_at_embedding_constant_is_rational`（u=5 ⟹ σ=5）。全 14 archimedean 测试通过；340 giac-core lib release 全过（11.11s）；#7 bar `quartic_a4_galois_dim_le_12` 0.11s 无回归；clippy 干净。

### Step 4 —— 朴素 LLL（f64，dim ≤ 5）【✅ 已实现】

- **做什么**：`lattice::lll(basis: Vec<Vec<f64>>) -> Option<Vec<Vec<f64>>>` —— 实格 LLL 规约（δ=3/4, η=1/2），dim ≤ `r1+r2-1 ≤ 5`。标准 LLL + f64 Gram-Schmidt（`gso` 给 μ/b*/||b*ᵢ||²），~90 行。每轮：找最深的 `|μ_{i,j}|>η` 做一次 size-reduction（`b_i -= round(μ)·b_j`）后重算 GSO；无 size-reduction 时扫 Lovász 找首个违反并 `swap(b_i,b_{i-1})`；两者皆无 → 规约完成。
- **接口偏离**：spec 写 `(Vec<Vec<f64>>, bool)`；实际用 `Option<Vec<Vec<f64>>>`（None=skip）—— 与 `field_signature`/`embeddings` 的 sound-skip `Option` 模式一致，Step 5 失败即 skip，bool 的"返原基"语义对消费者无额外价值（调用方持有原基）。
- **约束**：f64 LLL 对近退化基数值脆弱 —— `ponytail:` 注明上限（bigfloat/dashu）；**post-reduction sanity** `lovasz_ok`：所有 GSO 范数 finite-positive AND `|μ_{i,j}|≤η+tol` AND Lovász within tol，否则返 `None`（下游 skip，sound）。`LLL_MAX_ITERS=1000` 兜底（精确算术 LLL 必终止，触顶即数值失败信号）。ragged/NaN 输入直接 None。
- **文件**：新 `lattice.rs` 子模块（纯 f64 泛用格规约，与 `archimedean` 的 field→嵌入解耦）。**工作量**：中。**依赖**：无（纯 f64）。**PARI**：`qflll` 对照（2D `[[1,2],[3,4]]→[[1,0],[0,±2]]` 对齐；3D 用 |det| 保持 + Lovász 不变量，δ 不同不强求逐向量一致）。
- **测试**：8 个新增 —— 2D/3D PARI 对照、already-reduced、single/empty、degenerate→None、NaN→None、ragged→None。全 348 giac-core lib release 11.12s；#7 bar 0.11s 无回归；clippy 干净。

### Step 5 —— 基本单位系 + 挠群（f64+verify）【✅ 已实现（r=0 虚二次 + r=1 全实二次 + r≥2 全实内部）；剩余 sound-skip】

- **做什么**：`unit_group::fundamental_units(field) -> Option<Vec<HighFirstQ>>` + `unit_group::torsion(field) -> Option<Torsion>`。
- **已落地范围**：
  - **全实二次域 + 极大幂基**（`r₂=0`, `r=1`, `ℤ[α]=𝔬_K`）：
    - `torsion`：`μ(K)={±1}`（全实域无非实单位根，exact）。
    - `fundamental_units`：`ℤ[α]` 有界枚举（`UNIT_COORD_BOUND=256`），精确 `|N|=1` 滤子（二次闭式 `a²−c₁ab+c₀b²`），用 Step 3 嵌入取极小 regulator 的单位 ε。PARI `bnfinit(t²-d).fu` 对照（d=2,3,6,7 全过，regulator 误差 <1e-6）。
  - **虚二次域 + 极大幂基**（`d_k=2`, `r₂=1`, `r=0`）：`torsion` 按 `D_K=disc(m_α)` 分类（门限 `power_order_is_maximal`）：`D_K=-3→6 阶 {±1,±ω,±ω²}`、`D_K=-4→4 阶 {±1,±i}`、`else→2 阶 {±1}`，非实根作 `ℚ(α)` 元素返回；`fundamental_units` 返 `Some([])`（空基 = 平凡满自由基，无需认证）。
  - **全实 degree 3–4（`r≥2`）**：`fundamental_units_real_multi`（内部）— 坐标有界枚举 + f64 norm 预滤 + 精确 `|N|=1` + 按 log-norm 排序 + **incremental LLL**（`lattice::lll_incremental`）折入全部生成元得满 `r` 阶 ℤ-基 + 变换 `T`，`compose_unit` 提升为精确 ℚ(α) 单位并 `|N|=1` 验证。PARI regulator 对照（2 cubic + 1 quartic 全过）。**公开 `fundamental_units` 仍返 `None`**：缺 index-1 运行时证书（Friedman `R>g(1/|D|)` 对小判别式太弱），子格基喂 (B) 会 under-count mod 2 ⟹ unsound。
- **sound-skip（None）**：非实域 rank≥1（复 torsion via `x^k−1` 根 + 复 log-map 延后）；非极大虚二次序（非实根非 ℤ[α]-整，需 index）；非极大全实序；全实 degree≥5；大 regulator（坐标 > bound）；`c₀,c₁` 不入 i64。
- **新增基础设施（Step 5/7 共用）**：
  - `number_field_arith::poly_discriminant_low` —— 稠密 ℚ Sylvester 结式 + 判别式（`disc=(−1)^{n(n−1)/2}·Res(m,m′)/lc`）。
  - `number_field_arith::power_order_is_maximal(field) -> Option<bool>` —— 在每个 `p | disc(m_α)` 跑 Dedekind 判据认证 `[𝔬_K:ℤ[α]]=1`；`Some(true)`=极大，`Some(false)`=非极大/不可认证（下游 skip），`None`=tower/d≤1。Step 7 类群复用。
  - **`dedekind_index_ok` bug 修复**：原代码把 `M=∏g_i^{e_i}` 累加 **mod p**（丢了 p-adic 高位），污染 `T=(m−M)/p mod p`，对**分支素数**给错判据（全实二次 p=2：d=3 误判 index|、d=5 误判 index∤）。改为 **mod p²** 累加（规范 [0,p) 提升，gcd 结果与提升无关 — Cohen §6.2.6）。原 (A) 测试从未检分支素数故漏网；新增 `scan_ramified_prime_dedekind_fix_processes_p2_in_q_sqrt3` 回归守护（ℚ(√3), u=2, p=2 现判 Processed/v_𝔭=2）。
- **约束**：f64+verify 框架就位（嵌入 + 精确范数）；二次闭式范数已精确，无需 ℚ 重建。`ponytail:` 上限 `UNIT_COORD_BOUND=256`（大 regulator → None，sound）；升级路径 Buchmann 亚指数单位搜索。
- **文件**：`unit_group.rs` + `lattice.rs` + `number_field_arith.rs`。**依赖**：Step 1、3、4。**PARI**：`bnfinit(f).fu`/`.reg`/`.tu`、`poldisc`、`nf.disc`。
- **测试**：r≥2 多单位 3（cubic×2 + quartic×1 PARI regulator）+ 虚二次 torsion 5（Q(i) 4 阶 / Q(√-3) 6 阶 / Q(√-5) 2 阶 / 非极大 None / degree3 复 None）+ 全实二次 fund_units 4 PARI 对照。全 370 giac-core lib release 21.40s；clippy 干净。
- **剩余（后续回合）**：r≥2 全实域 index-1 运行时证书（Friedman/Mellin 或接类群）→ 公开 API 解锁；非实域 rank≥1 复 torsion（`x^k−1` 在 K 中的根，接 `poly_algext_roots`）+ 复 log-map（Step 3 已有复嵌入）。

### Step 6 —— 单位 mod-2 求解（B）【✅ 已实现】

- **做什么**：`unit_group::unit_is_square(field, η, units, torsion) -> bool`。Dirichlet 分解 `η = ζ·∏ε_i^{a_i}`：f64 Minkowski log-map `L(η) = Σ a_i L(ε_i)` 的 `r×r` 系统（Gauss-Jordan）+ 四舍五入 `a_i`；精确 `ζ = η·∏ε_i^{-a_i}`（`compose_unit`）；验证 `ζ^|μ(K)|=1`（exact ⟹ ζ∈μ(K)）；`η∈(𝔬_K×)² ⟺` 所有 `a_i` 偶 AND `ζ∈μ(K)²`（torsion 平方：枚举 `ζ'∈roots` 查 `ζ'²=ζ`，exact）。全正已由 (A_inf) 保证。
- **已落地范围**：certified 域类型全过 — r=0 虚二次（`units=[]`，纯 torsion 平方判定）+ r=1 全实二次（单单位 + `{±1}` torsion）。PARI `nfeltissquare(bnfinit(f), η)` 对照：ℚ(√2) ε²/ε⁴/1 yes、ε/-1 no；ℚ(√3) ε² yes、ε/ε³/-1 no；ℚ(i) -1=i² yes、±i no；ℚ(√-3) α/(-α-1) yes、(α+1)/(-α)/-1 no（μ₆²=μ₃）。
- **约束**：**verify-and-skip**：embeddings 不可信 / `r×r` 系统奇异 / 舍入 `a_i` 不恢复 μ(K) 元 / `|a_i|>1e15` → 返 `true`（sound false-pass，`sqrt_base_case` ℚ 证书兜底）。**前提**：调用方须传 *certified full* 基（子格基 under-count mod 2 ⟹ unsound）；r≥2 公开 API 因 Step 5 仍 `None` 故 (B) 暂不接 r≥2。
- **文件**：`unit_group.rs`。**工作量**：中。**依赖**：Step 3（嵌入）、Step 5（单位系 + torsion）。**PARI**：`nfeltissquare(B,η)`。
- **测试**：4 个新增（ℚ(√2) / ℚ(√3) / ℚ(i) / ℚ(√-3) 各一，覆盖 r=0/r=1 × torsion 阶 2/4/6 × square/non-square）。全 374 giac-core lib release 21.45s；clippy 干净。

### Step 7 —— 类群 2-挠 + J 主性（C）【✅ Step 7a 直接主性已实现；Step 7b 完整类群剩余】

- **Step 7a（已实现）—— 直接主性 via 生成元搜索**：`class_group::ideal_is_principal(field, &scan) -> Option<HighFirstQ>`。核心等价：`J = ∏𝔭^{v_𝔭/2}` 主 ⟺ ∃ `γ ∈ J` with `|N(γ)| = N(J)`（`(γ) ⊆ J` + 等范数 ⟹ `(γ)=J`）。实现：bounded `ℤ[α]` 枚举（`IDEAL_GEN_COORD_BOUND` 按度 2/3/4 = 256/32/12；`(2B+1)^d ≲ 3·10⁵`）+ f64 `|N|` 预滤 + 精确 Krylov `|N(γ)|=N(J)` + 成员校验 `γ ∈ J`（复用 (A_fin) `compute_ideal_valuations` 算 `v_𝔭(γ) ≥ e_𝔭`，按 `(p, g_i)` 匹配素理想）。门限 `power_order_is_maximal`（`𝔬_K=ℤ[α]` 使枚举覆盖 𝔬_K）+ `scan.parity_ok && all_processed` + `N(J) ≤ 10⁷`。`N(J)=1` ⟹ `J=𝔬_K` 平凡主。**sound-skip**：生成元不在搜索界内 / `J` 真非主 / 大范数 → `None`（不可区分，均 defer `sqrt_base_case`）。
- **Step 7a 管线 `class_group::hasse_is_square(field, u) -> Option<bool>`**：Scholz √-判定 `Some(true)`=认证平方 / `Some(false)`=认证非平方 / `None`=skip。链：`(A_fin)` parity → `(A_inf)` 全正 → `(C)` J 主+γ → `(B)` `η=u/γ²` 单位平方（Step 6 `unit_is_square` 改 `Option<bool>` 认证版）。各 `Some(false)` = 平方 obstruction（奇赋值 / 非全正 / η 非单位平方），均 exact 认证 → sound。
- **Step 7a 测试（PARI `nfeltissquare` 交叉验证，13 个）**：ℚ(√2) u=2/6+4√2 yes、u=2+2√2/-1 no（(A_inf)）；ℚ(√3) u=3 yes、u=6-3√3 no（(B) ε=2-√3 全正非平方）；ℚ(i) -1 yes、i no；ℚ(√-5) u=2 → `None`（类数 2，𝔭/2 非主，sound skip；PARI 0）；`ideal_is_principal` 隔离测试（(√2) 主 / 𝔭 非主 skip / 平凡 J）；0=0²。全 387 giac-core lib release 21.77s；clippy 干净。
- **Step 7b（剩余）—— 完整 Buchmann 类群**：`class_group_2torsion` + `bnfisprincipal`-式 J 主性。Minkowski 界 `B_M = (4/π)^{r2}·d!/d^d·√|disc|`；枚举 `norm ≤ B_M` 素理想 → 收集关系（小范理想 LLL 短向量找主生成元）→ HNF/SNF → 类群结构 → Sylow-2。**解锁**：证非主性 `Some(false)`（当前只能 `None` skip）；大范数 `J` 的 LLL 短向量生成元搜索（lift `N_J_MAX=10⁷` ceiling）；r≥2 index-1 单位证书（接类群）。需 `disc(m_α)` + (A) 完备（P3 大素因子 + P5 `p∣index` Montes）。~400-600 行，多轮。
- **文件**：`class_group.rs`（Step 7a 已建）。**依赖**：Step 1-6。**PARI**：`bnfinit(f).clgp` + `bnfisprincipal(B,J)` + `nfeltissquare(B,u)`。

### 终态接线（Step 7 完成后）

`sqrt_fmodule` 调用链：norm-square → **(A_fin)** → **(A_inf)**[Step 1] → **(C)** J 主性 + 取 γ[Step 7] → η=u/γ² → **(B)** 单位 mod-2[Step 6]。任一失败 `return None`；任一 skip → 续走 `sqrt_base_case`（ℚ 验证兜底，sound）。

```mermaid
flowchart TD
    S1["Step 1: signature r1,r2 + totally-positive (A_inf)<br/>Sturm + bisection, pure Q arithmetic"]
    S2["Step 2: refactor (A) to IdealValuationScan<br/>expose valuations, g_i, u_int, d_u"]
    S3["Step 3: numerical embeddings<br/>real Sturm-bisection + complex companion-Eigen, f64"]
    S4["Step 4: naive LLL (f64, dim <= 5)"]
    S5["Step 5: fundamental units + torsion<br/>LLL on Minkowski log-lattice, f64+verify"]
    S6["Step 6: unit mod-2 (B)<br/>express eta, check even exp + torsion square"]
    S7["Step 7: class group 2-torsion (C)<br/>Buchmann sub-exponential + principality of J"]
    S1-->S2-->S3-->S4-->S5-->S6
    S7-->S6
```

---

## 进度跟踪

| 步骤 | 状态 | 备注 |
|---|---|---|
| Step 1 signature + (A_inf) | ✅ 已实现 | `algebra/archimedean.rs`：`field_signature` + `is_totally_positive`（自洽稠密有理 Sturm + 有理 Horner 二分）；接线 `poly_roots.rs` (A_fin) 之后；PARI 交叉验证 9 测试全绿；#7 bar 0.63s |
| Step 2 IdealValuationScan | ✅ 已实现 | `compute_ideal_valuations -> IdealValuationScan`（暴露 entries (p,g_i,e_i,f_i,v_𝔭), u_int_low, d_u, parity_ok, all_processed）；`ideal_valuations_parity_ok` 薄包装向后兼容；4 scan 测试全绿；signature/totally_positive 解耦留 archimedean |
| Step 3 数值嵌入 | ✅ 已实现 | `archimedean::embeddings`（实 Sturm 二分 + 复 companion-Eigen，f64）+ `eval_at_embedding`；`Embeddings.reliable` verify-and-skip；PARI `polroots` 对照 5 测试 |
| Step 4 朴素 LLL | ✅ 已实现 | `lattice::lll`（δ=3/4 f64 LLL + `lovasz_ok` post-sanity + `LLL_MAX_ITERS`）；后续扩展 `lll_with_transform`/`lll_reduce_overcomplete`/`lll_incremental`（Step 5 多单位用）+ ill-conditioning 守卫 |
| Step 5 基本单位系 + 挠群 | ✅ 已实现（剩余 sound-skip） | r=0 虚二次 torsion（D_K 分类）+ r=1 全实二次 fund_units + r≥2 全实内部 `fundamental_units_real_multi`（incremental LLL，PARI regulator 对照）；公开 API r≥2 仍 None（缺 index-1 证书）；复域 rank≥1 + 非极大虚二次 None |
| Step 6 单位 mod-2 (B) | ✅ 已实现 | `unit_group::unit_is_square`（Dirichlet 分解 + mod-2 指数 + torsion 平方，sound false-pass）；certified 域类型 r=0/r=1 PARI `nfeltissquare` 对照 4 测试；r≥2 待 Step 5 解锁 |
| Step 7a (C) 直接主性 + Scholz 管线 | ✅ 已实现 | `class_group::ideal_is_principal`（生成元搜索 `|N(γ)|=N(J)`+赋值成员校验，sound-skip）+ `hasse_is_square`（A_fin→A_inf→C→B 管线，`Option<bool>` 认证）；PARI `nfeltissquare` 交叉验证 13 测试；ℚ(√-5) 非主 J sound-skip |
| Step 7b (C) 完整类群 (Buchmann) | 待实现 | `class_group_2torsion` + `bnfisprincipal` 式主性；证非主性 `Some(false)` + 大范数 LLL 短向量 + r≥2 index-1 证书；需 (A) 完备 |

## 参考

- Keith Conrad, *The Local-Global Principle* — https://kconrad.math.uconn.edu/blurbs/gradnumthy/localglobal.pdf
- Cohen, *A Course in Computational Algebraic Number Theory* §6（Dedekind / 素理想）、§4（范数与赋值）、§5/§6（Minkowski / 类群 / 单位 Dirichlet）
- Buchmann, *A subexponential algorithm for the determination of class groups and regulators of algebraic number fields* (1990) — Step 7 类群算法
- Pohst / Zassenhaus, *Algorithmic Algebraic Number Theory* — 单位搜索 / LLL 应用
- PARI/GP `bnfinit` / `bnfisprincipal` / `bnfissunit` / `nfeltissquare` — 完备 oracle（首选）
