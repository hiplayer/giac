# 实二次域解析类数证书 — 原理与流程

**状态:** reference（已落地，commit `498a929` / `766423d`）
**落点:** `giac-rs/crates/giac-core/src/algebra/class_group.rs`
**相关:** [GIAC-p2-bnf-pari-alignment](issues/GIAC-p2-bnf-pari-alignment.md) §2a-S2b-full-real、[giac-poly-f5-fglm-hasse-proof](giac-poly-f5-fglm-hasse-proof.md)（Hasse 数学正确性线）
**快照:** 2026-07-03

---

## 1. 证书证明什么（soundness target）

对实二次域 `K = ℚ(√d)`（基本判别式 `D = D_K > 0`），证书返回 `h(K) ∈ ℤ_{>0}`，且**保证等于真类数**，否则拒绝（`None`）。soundness 是硬约束 —— 永不返回错的 `h`。

### 为何实二次 `h > 1` 比已解锁情形难

| 情形 | 证书 | 为何 sound |
|------|------|-----------|
| `h = 1` | `all_basis_primes_principal`（每个基素理想都主） | Minkowski：每类有 `norm ≤ M_K` 代表 ⟹ 全主 ⟹ 平凡类群 |
| 虚二次 `h > 1` | `count_reduced_forms`（既约**定**型二元二次型计数） | 有限计数（定型 ⟹ 有限），是定理 |
| **实二次 `h > 1`** | **解析类数公式**（本文档） | 见 §2 |

实二次的根难处：**范数型是不定的**（`N(a+bα) = a² − c₁ab + c₀b²` 可正可负），生成元轨迹无界 ⟹ 既没有「定型型有限计数」，也没有「坐标界内枚举即穷尽主理想」。

Buchmann 关系格法给出的是

\[
h_{\text{buchmann}} = [\mathbb{Z}^k : \langle R\rangle] = h\cdot[L:\langle R\rangle] \;\ge\; h
\]

即**仅上界** —— 枚举到的关系 `\langle R\rangle` 可能只是真关系格 `L = \ker(\mathbb{Z}^k \to \mathrm{Cl}(K))` 的**真子格**（`[L:\langle R\rangle] > 1`）。`h > 1` 时无法判定枚举是否完整 ⟹ 旧代码 sound-skip。

---

## 2. 证书的数学内容

**Dirichlet 实二次类数公式**（`D = D_K > 0` 基本判别式，`χ_D` 为二次特征）：

\[
h\cdot R = \frac{\sqrt D}{2}\,L(1,\chi_D),\qquad
L(1,\chi_D) = -\frac{1}{\sqrt D}\sum_{a=1}^{D-1}\chi_D(a)\ln\left|\sin\frac{\pi a}{D}\right|
\]

（`L(1,χ_D)` 闭式来自偶特征的高斯和 `τ(χ_D)=√D` 与 `|1−ζ_D^a|=2|\sin(πa/D)|`；`\sum χ_D(a)=0` 消去 `\ln 2`。）消去 `\sqrt D` 得计算闭式：

\[
\boxed{\;h = -\frac{1}{2R}\sum_{a=1}^{D-1}\chi_D(a)\,\ln\left|\sin\frac{\pi a}{D}\right|\;}\qquad
R = \ln|\sigma_{\max}(\varepsilon_0)|,\quad \chi_D(a) = (D|a)\;\text{Kronecker}
\]

**这是闭式定理** —— 不做关系枚举，**没有「关系格完备性」问题**。证书的 soundness 落在三点：

1. **公式是定理**（Dirichlet）。
2. **`R` 精确**：`\varepsilon_0` 是真基本单位。`fundamental_units` 找不到时返 `None` ⟹ **绝不用错的 `R`**（用错单位 ⟹ `R` 错 ⟹ `h` 错 —— 关键防线）。
3. **舍入无歧义**：`h` 是整数，f64 算出实数 `h_{\text{real}}` 后取最近整数；若 `|h_{\text{real}} − \mathrm{round}| > 0.01`（距半整数太近）⟹ `None`。

### 数值验证（Pari 金值）

| 域 | `D` | `\varepsilon_0` | `R` | `h`（公式 / Pari） |
|----|-----|-----------------|-----|--------------------|
| ℚ(√2) | 8 | `1+√2` | `ln(1+√2)≈0.8814` | 1 |
| ℚ(√3) | 12 | `2+√3` | ≈1.3170 | 1 |
| ℚ(√5) | 5 | `(1+√5)/2` | ≈0.4812 | 1 |
| ℚ(√6) | 24 | `5+2√6` | ≈2.2924 | 1 |
| ℚ(√7) | 28 | `8+3√7` | ≈2.7686 | 1 |
| ℚ(√10) | 40 | `3+√10` | ≈1.8184 | **2** |
| ℚ(√15) | 60 | `4+√15` | ≈2.0634 | **2** |

`√10`、`√15` 即本轮解锁（原 sound-skip）。`\sqrt 2`、`\sqrt 3` 手算逐项核验闭式（`\chi_8, \chi_{12}` 对称 `[1,−1,−1,1]`，和 `= −2R` ⟹ `h=1`）。

### 精度分析

f64（~15–16 位）。和 `Σ χ_D(a)\ln|\sin(πa/D)|` 有 `D` 项，每项 `O(\ln D)`，部分和达 `O(D\ln D)`，最终和 `O(\sqrt D)`（因 `h∼\sqrt D` 量级）。抵消比 `∼\sqrt D\ln D`。对 `D = 10⁶`：绝对误差 `∼ 10⁵·2^{-52} ∼ 2\text{e−11}`，`h ∼ 10³` 量级 ⟹ 舍入远无歧义。`D ≤ 10⁶` 下 f64 充分。

---

## 3. 设计的流程（两层，镜像虚二次架构）

```
                         ┌─ Layer 1: 独立 oracle（解析，自身 sound）────────┐
class_number(P) ─────────┤                                                     │
   real-quad dispatcher ─▶ class_number_real_quad_analytic(field)              │
                          │   1. 门控：deg-2 极大序，D = disc(m_α) > 0          │
                          │   2. ε₀ = fundamental_units(field)?                 │
                          │      (None → 拒绝，无 R)                            │
                          │   3. R = ln max|σ(ε₀)|  (archimedean::embeddings)   │
                          │   4. Σ_{a<D} χ_D(a)·ln|sin(πa/D)|  (f64)            │
                          │   5. h = −sum/(2R)，舍入 + 歧义检查                  │
                          │   → Some(h) 或 None                                 │
                          └─────────────────────────────────────────────────────┘
                          ┌─ Layer 2: Buchmann 交叉校验（defense-in-depth）────┐
class_number_general ─────┤  1. ramification + enumerate_relations_deg2         │
   (测试 / 管线验证)      │  2. h_buchmann = class_number_from_relations        │
                          │     = ∏ SNF 不变因子  (恒为 h 的倍数, 上界)         │
                          │  3. cert_h = class_number_real_quad_analytic?       │
                          │  4. cert_h == h_buchmann ?                          │
                          │     ├─ 相等 → Some  (管线无 bug + 关系格完备)       │
                          │     └─ 不等 → None  (关系不全 或 管线 bug)          │
                          └─────────────────────────────────────────────────────┘
```

**关键架构决定**：用户级 `class_number` 经 **Layer 1 直接返回**（如同虚二次经 `count_reduced_forms` 直接返回）—— sound，**不依赖 Buchmann 管线**。Layer 2 的交叉校验在 `class_number_general` 里，是**管线一致性验证**（在测试中演练），不是用户结果的 soundness 闸门。

### 代码锚点

| 角色 | 位置 |
|------|------|
| 解析证书主体 | `class_group.rs::class_number_real_quad_analytic` |
| Kronecker 符号（i128 快路径） | `class_group.rs::kronecker_symbol` |
| Kronecker 符号（BigInt 兜底，复用原语） | `class_group.rs::kronecker_symbol_bigint` |
| 基本单位 | `unit_group.rs::fundamental_units` — 实二次全经连分数（**无坐标枚举**）：`√d`-CF ±1 Pell（`d≡2,3mod4`）/ `(1+√d)/2`-CF ±4 半整数（`d≡1mod4`），`cf_period_cap=4√d+8` 可证界（Lagrange 周期 `≤2√d`） |
| Regulator 计算 | 复用 `embeddings` + `σ_max` 模式（同 `eval.rs::eval_bnfregulator`） |
| Buchmann 交叉校验 | `class_group.rs::class_number_general`（real-quad 分支） |
| dispatcher | `class_group.rs::class_number` |

---

## 4. 「完备性」在此设计里的精确含义

问题的要害。有**两个**「完备性」概念，本设计用一个把另一个变免费：

**(a) 关系格完备性**（Buchmann 本义）：`\langle R\rangle = L`，即枚举关系生成 `\mathbb{Z}^k \to \mathrm{Cl}(K)` 的**完整核**。实二次 `h > 1` 暴力枚举固定坐标界只能得子格 ⟹ `h_{\text{buchmann}} = h\cdot m`（`m > 1`）—— 旧代码无法证明。

**(b) 解析证书的「完备性」**：闭式定理**无需枚举**，没有 (a) 的问题。它只有两种「不完备」：① 找不到基本单位（`R` 未知）② 精度歧义。两者都 → `None`。

**核心洞察**：用 (b) 的闭式公式**绕过** (a) 的关系格完备性问题 —— 不需要证明 `\langle R\rangle = L`，因为 `h` 来自恒等式而非搜索。然后 Layer 2 的交叉校验 `h_{\text{buchmann}} == h_{\text{analytic}}` **同时**证明两件事：

- **管线无 bug**（Hensel 赋值 / SNF / k×k 子式 gcd 与独立 oracle 吻合），
- **关系格完备** `\langle R\rangle = L`（因 `h_{\text{buchmann}} = h\cdot[L:\langle R\rangle]` 且 `= h` ⟹ `[L:\langle R\rangle] = 1`）。

即：独立 oracle 让「关系格完备性」从**需要直接证明的难题**变成**交叉校验免费送出的推论**。这正是虚二次 `forms` 证书的同一套路 —— 虚二次「闭式」是组合计数，实二次「闭式」是解析三角和 + regulator，两者都是定理的有限计算。

---

## 5. 与原方案（单位约化完备性证书）的区别

对齐文档原 `ponytail:` 设想的实二次证书是**「单位约化完备性证书」**（代数路径）：用基本单位 `\varepsilon_0` 做 fundamental domain 枚举 `norm ≤ M_K` 的主理想，导出坐标界

\[
B^* = \left\lceil\sqrt{M_K \cdot \sigma_1(\varepsilon_0)}\,\right\rceil
\]

证明在该界内枚举穷尽 ⟹ `\langle R\rangle = L`。这是**直接证 (a)**，但 `B^*` 对大 regulator 爆炸（`\varepsilon_0 \sim 10^{37}` ⟹ `B^* \sim 10^{37}`，枚举不可行）。

本轮改走**解析路径**（闭式公式证 `h`，交叉校验得 (a)）：

| 路径 | 证明对象 | 瓶颈 |
|------|---------|------|
| 代数（单位约化） | 直接证 (a) `\langle R\rangle = L` | `B^*` 界内穷尽枚举 ⟹ 大 regulator 爆炸 |
| **解析（本文档）** | 闭式得 `h`，交叉校验免费得 (a) | 只需 `R`（基本单位），与大 regulator 无关（除找单位本身） |

两条路都解锁实二次 `h > 1`，但解析路径的瓶颈只在「找基本单位」（`fundamental_units` 256-bound），不在「枚举完备性界」。

---

## 6. 边界（`ponytail:`，证书**不**覆盖的）

- **基本单位求解——全经连分数，无坐标枚举（`ponytail:` 架构升级）**：基本单位的坐标上界是 `exp(Θ(√D))`（大 regulator 的字面含义，如 ℚ(√163) `ε~10⁸`、ℚ(√d) `ε~10³⁵`），**任何多项式坐标界都覆盖不了**——这正是旧的 256-bound 暴力枚举的根本缺陷（要么太小漏大 regulator，要么取指数界不可行）。连分数把**工作量**与**输出大小**解耦：CF 周期 `ℓ ≤ 2√d`（Lagrange）⟹ `O(√d)` 渐近分数即找到单位，与单位大小无关。`unit_group.rs` 现按 `c1` 奇偶分派两条 CF：
  - `c1` 偶（`d≡2,3mod4`，`𝓞_K=ℤ[√d]`）：`fundamental_unit_via_cf_sqrt_d(d)` —— `√d`-CF，首 `p_n²−d·q_n²=±1` 渐近分数 = 极小 Pell 解 = 基本单位。
  - `c1` 奇（`d≡1mod4`，`𝓞_K=ℤ[(1+√d)/2]`）：`fundamental_unit_via_cf_half(d)` —— `(1+√d)/2`-CF，首 `(2p_n−q_n)²−d·q_n²=±4` 渐近分数 = 基本单位。**关键**：`(1+√d)/2`-CF 的 Legendre 判据是 `d>4`（覆盖**所有** squarefree `d≡1mod4`，`d≥5`，含 `d=5`——`√d`-CF±4 方案 `d>16` 门漏掉的唯一例外）；且所有 `𝓞_K` 单位在该 CF 下统一是 `±4`（整数单位 `(x,y)` 偶 ↔ `(2x,2y)` 仍 `±4`）⟹ 首 `±4` 渐近分数 = 基本单位（无 `±1`/`±4` 混合、无 min-|ε| 追踪、无 `d>16` 门、无 `d=5` 特例）。
  - `cf_period_cap = 4√d+8`（≥4× 余量，可证界，取代旧 `2·10⁶` 安全网）；BigInt 渐近分数（坐标可超 i128），CF 状态 `(m,dd,a)` 留 i128（界 `2√d`）。
  - **解锁** `class_number(ℚ(√163))=1`（`d≡3mod4`，R≈18.669）、`class_number(ℚ(√265))=2`（`d≡1mod4` 大 regulator，R≈9.4046）、`class_number(ℚ(√5))=1`（`d=5` 极大序，R≈0.4812，原 sound-skip）、`class_number(ℚ(√409))=1`（`d≡1mod4` 大 regulator，素数 disc 409 > 281，ε 坐标 ~10¹⁰，R≈26.1342，原认证侧 sound-skip）；`bnfregulator`/`bnfunits` 同步受益。
- **`power_order_is_maximal` 极大性认证——`disc` 完全分解（`ponytail:` 架构升级）**：原实现用 `small_primes(60)`（前 60 个素数，到 281）+ 「remainder>1 → `Some(false)`」，对 `disc` 含 `>281` 素因子的域（如 ℚ(√409)，disc=409 素数）在**认证侧** sound-skip——即便 Dedekind 判据本身对任意大素数成立、CF 也能找到单位。现 `number_field_arith.rs::distinct_prime_factors_i64(disc)` 把 `|disc|` 完全分解为素因子（试除到 `min(√|disc|, 10^7)` + Miller-Rabin `is_prime_i64` 收尾），对**每个**素因子跑 Dedekind。`ponytail:` 上界 `|disc| > 10^14` 且剩余 composite > 试除上限 → `None` → `Some(false)` sound-skip；升级路径 = Pollard-rho。覆盖范围内（`|disc| ≤ 10^14`，含所有在域二次 `disc`）任意大素数 `disc` 均可证。
- **`D > 10^6`**：求和可行性 cap（10⁶ 次迭代），**非 soundness 缺口**，纯算力。可升级为 Poisson 求和 / baby-step giant-step 求 `L(1,χ)`。
- **`deg ≥ 3`**：无解析闭式（高次 Dedekind L 函数）⟹ `None`，仍需完整 Buchmann 完备性证书（C 路径）。
- **f64 精度**：舍入歧义 ⟹ `None`（已验证 `D ≤ 10^6` 下精度远够；大 regulator 的 `R` 用 f64 表达 `ε ~ 10^8` 仍只有 `~10⁻¹⁶` 相对误差 ⟹ `h` 误差 `~10⁻¹⁴`，舍入无歧义）。

---

## 7. Kronecker 符号实现要点

- `kronecker_symbol(a: i128, n: i128) -> i64`：i128 全程，2-adic 分解 + Jacobi 互反律。覆盖 `D ≤ 10^6`（求和 cap）且 i128 本身覆盖到 `~10^{38}`，远超任何二次域判别式。**解析证书热路径直接用此版**。
- `kronecker_symbol_bigint(a: &BigInt, n: &BigInt) -> i64`：`pub(crate)` 复用原语，两操作数落 i128 时走 i128 快路径（零 BigInt 分配，热路径快 `~100×`），超 i128 才走 BigInt 算法（`.bit(0)` 判奇偶、`>>= 1` 除 2、`%` 取模，无需 `num_integer` import）。镜像本文件 `det_minor_i64 → det_bareiss` 的 dispatch 模式。供未来 `bnr*`（conductor / 射线类群特征）、`galois*`（分裂域二次剩余）大判别式消费。
- BigInt 路径 2-adic 分解用 `while !bit(0) { >>=1 }` 计数 —— 对 `n = 2^k` 大 `k` 是 `O(k)` 次 BigInt 缩位移位；可升级为 `to_u64_digits()` 首非零字 + `trailing_zeros` 一次取尽，当前无超 i128 的 `2^k` 型消费者，留待需要时。

---

## 参考

- 对齐跟踪：[GIAC-p2-bnf-pari-alignment](issues/GIAC-p2-bnf-pari-alignment.md) §2a-S2b-full-real
- 虚二次对照（组合闭式）：`class_group.rs::count_reduced_forms` / `class_number_imag_quad_forms`
- Buchmann 管线：`class_group.rs::enumerate_relations_deg2` / `ramification_relations` / `class_number_from_relations` / `snf_bigint`
- Hasse 数学正确性线：[giac-poly-f5-fglm-hasse-proof](giac-poly-f5-fglm-hasse-proof.md)
- 上游对照：Pari/GP `bnfinit(P).no` / `quadclassunit(D).no` / `bnfinit(P).reg`
