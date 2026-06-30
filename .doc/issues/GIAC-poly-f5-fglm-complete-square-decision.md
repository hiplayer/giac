# GIAC — F5 FGLM / F-module √-reduction 完备平方判定算法

**状态:** open（Phase A 已落地 2026-06-30；Phase B-D 待立项实现）
**类型:** 算法升级 / 完备性闭环
**来源:** [GIAC-poly-f5-fglm-debug-postmortem](GIAC-poly-f5-fglm-debug-postmortem.md) §3.3 follow-up 范畴之外、[GIAC-poly-f5-fglm-fullchain-fuel-audit](GIAC-poly-f5-fglm-fullchain-fuel-audit.md)（fuel 是封顶，本 issue 是替换算法）
**根因:** 当前 `sqrt_base_case` 是「fast path + 8 次 lift 截断」的启发式；soundness 由 ℚ 验证保证，completeness 由 fuel=8 / `nf>3` 跳过 / `target_bits=320` / `small_primes(60)` 四个工程封顶兜底，无完备性定理
**Rust 落点:** `giac-core::algebra::poly_roots`（`sqrt_fmodule` / `sqrt_base_case`），新模块 `giac-core::algebra::number_field_arith`（待建）
**相关:** [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)、[GIAC-poly-f5-fglm-math-background](GIAC-poly-f5-fglm-math-background.md)
**快照:** 2026-06-29

---

## 0. 背景与目标

### 0.1 当前实现是什么

`sqrt_base_case(field, u, m_gen, fuel: &Fuel)`（`poly_roots.rs:1460`）解决：u ∈ K=ℚ(gen) 是否平方，是则返回 δ 使 δ²=u。

策略：模小素数 p → factor m_gen mod p → 各分量 Tonelli 开 √ → 2ⁿᶠ 个符号组合 CRT 回 δ₁ → p-adic Newton 提升到 mod p^320 → 有理重构 → **ℚ 上验证 δ²=u**。`fuel=8` 封顶总 lift 次数。

### 0.2 数学保证现状

| 性质 | 现状 |
|---|---|
| **Soundness**（返回 Some ⟹ 真是 √u） | ✅ 可证（第 1564 行 ℚ 等式验证是 certificate） |
| **Completeness**（u 真是平方 ⟹ 必找到） | ❌ 无定理。四个工程封顶（fuel=8 / nf>3 skip / target_bits=320 / small_primes(60)）均可静默错过真平方 |

### 0.3 目标

把 completeness 也升级到可证：用 Hasse 局部-全局原理替代 fuel 截断。**保持 soundness 不变**（ℚ 验证作为最后关卡）。

### 0.4 不做的事

- 不改 `sqrt_fmodule` 的 F-module 递归结构（属于 [GIAC-poly-f5-fglm-fullchain-fuel-audit](GIAC-poly-f5-fglm-fullchain-fuel-audit.md)）
- 不替换 p-adic 提升作为「构造 δ」的工具（它本身可证，只是当前用法被 fuel 截断）
- 不形式化证明（Lean / Coq）——本 issue 只交付可读的纸面证明 sketch + 实现；形式化是后续独立 issue

---

## 1. 数学骨架

### 1.1 Hasse 局部-全局（数域版）

**定理**：u ∈ K× 是平方 ⟺ 对所有素理想 𝔭 of 𝔬_K，u 在 K_𝔭 是平方。

**判定等价于两条同时成立**：

| 条件 | 含义 | 判定方法 |
|---|---|---|
| **(A)** 主理想 (u) 的素理想分解里所有指数 v_𝔭(u) 为偶 | (u) = 𝔞²，𝔞 是 fractional ideal | 算 (u) 在每个 𝔭 的赋值 |
| **(B)** u 的单位部分在 𝔬_K×/(𝔬_K×)² 里平凡 | 局部单位部分 mod 2 = 0 | Dirichlet 单位群结构 + mod-2 线性方程 |

(A) 是比 norm 更强的必要条件。norm 平方 ⟸ (A) 成立（但反向不成立，故 (A) 严格更强）。(A) ∧ (B) ⟺ u 是平方（完备）。

### 1.2 素理想分解 (A) 的具体计算

对 K = ℚ(α)（m_α 给定）：

1. 对每个有理素数 p：
   - 若 p ∤ disc(m_α)（good prime）：m_α mod p = ∏ g_i（不可约因子分解，已有 `factor_mod_irreducibles`）。每个 g_i 对应一个素理想 𝔭_i = (p, g_i(α))，剩余次数 f_i = deg(g_i)，分裂指数 e_i = 1（Dedekind）。
   - 若 p ∣ disc(m_α)（bad prime）：需 Dedekind/Kummer 的 ramification 分析，单独处理。
2. 对每个 𝔭_i，算 v_{𝔭_i}(u)：
   - 把 u 视为 ℚ[t]/(m_α) 元素，模 g_i 约化得 ū_i ∈ 𝔽_{p^{f_i}}。
   - v_{𝔭_i}(u) = v_p(N_{𝔽_{p^{f_i}}/𝔽_p}(ū_i))（norm 形式赋值公式，对 good prime 成立）。
3. (A) 通过 ⟺ 所有 v_{𝔭_i}(u) 偶。

**关键观察**：只需对 p ≤ 某个上界（由 |N_{K/ℚ}(u)| 决定）扫描——超出此界的 p 必有 v_𝔭(u)=0。终止性可证。

### 1.3 单位群 (B) 的具体计算

Dirichlet 单位定理：𝔬_K× ≅ μ(K) × ℤ^{r₁+r₂−1}，r₁/r₂ 是 K 的实/复 embedding 数。

1. 算 K 的 signature (r₁, r₂)（由 m_α 的实根数 + 复对数）。
2. 求基本单位系 ε₁,…,ε_{r₁+r₂−1}：**难件**，一般用 Minkowski 嵌入 + LLL 求近似关系，再精确化。
3. 把 u 的单位部分 η = u / ∏_𝔭 𝔭^{v_𝔭(u)/2} 表达在基本单位系 + 挠群生成元下。
4. 列 mod-2 线性方程组，判定 η ≡ η'² mod 𝔬_K× 是否有解。

### 1.4 复杂度对比

| 步骤 | 当前 fast path | 完备算法 |
|---|---|---|
| Negative filter | norm 平方（已落地） | norm 平方 + (A) 全指数偶（更强） |
| Positive 构造 | p-adic lift + ℚ 验证 | p-adic lift + ℚ 验证（不变） |
| 完备性来源 | fuel 截断（不可证） | (A) ∧ (B) ⟺ 平方（可证） |

---

## 2. 当前缺口盘点

### 2.1 已有（无需新写）

| 子件 | 位置 |
|---|---|
| `factor_mod_irreducibles(p, m_α) → Vec<PolyMod>` | `giac-poly::factor::mod` |
| `field.discriminant()` | `field_arith.rs`（待确认 API 名） |
| `field.dimension()` / generator minpoly 访问 | `ExtensionField` |
| norm 早退 `is_q_square(norm)` | `poly_roots.rs:1955-1968`（P5/F1 已落地） |
| p-adic Newton lift | `padic_sqrt_lift` `poly_roots.rs:1583` |
| 有理重构 | `rational_reconstruct` `poly_roots.rs:1327` |
| ℚ 验证 δ²=u | `poly_roots.rs:1559-1567` |

### 2.2 缺件（需新写）

| 子件 | 用途 | 工程量估 | 依赖 |
|---|---|---|---|
| **C1** trial division 素因子分解（BigInt → Vec<i64>） | 算 disc 的素因子 | 小（~30 行） | 无 |
| **C2** 素理想分解（good prime 路径） | (A) 的输入 | 中（~150 行） | C1 + `factor_mod_irreducibles` |
| **C3** 素理想赋值 v_𝔭(u)（norm 形式公式） | (A) 判定 | 中（~80 行） | C2 |
| **C4** bad prime ramification（p ∣ disc） | 补全 (A) | 中大（~200 行，需 Kummer/Dedekind） | C1 |
| **C5** K 的 signature (r₁, r₂)（实根计数） | (B) 的预备 | 小（~50 行，已有实根隔离算法可复用） | 无 |
| **C6** LLL 实现（或引入 `nalgebra` 已有） | (B) 基本单位系 | 中大（~400 行 或 引入 crate） | `nalgebra`（已允许） |
| **C7** 基本单位系求取 | (B) 核心 | 大（~600 行） | C5 + C6 |
| **C8** 单位部分 mod-2 求解 | (B) 判定 | 小（~100 行） | C7 |
| **C9** 接线进 `sqrt_fmodule`（filter 前置） | 落地 | 小（~30 行） | C2-C4 + C8 |

### 2.3 LLL 依赖决策点

giac-rs 工程规范（`giac-rust-engineering.mdc`）允许 `nalgebra`，禁止 `faer`。LLL 实现：
- 选项 a：自己写朴素 LLL（~200 行，d ≤ 12 时足够）
- 选项 b：引入 `oorandom` / `vecmath` 之外的 LLL crate（需审 WASM 兼容）
- 选项 c：复用 giac-poly 已有的 reduction primitive（待查）

**推荐 a**——d ≤ `POLY_ROOTS_DIM_HARD` 时朴素 LLL 完全够用，无新依赖、WASM 安全。

---

## 3. 实现优先级（四阶段）

### Phase A — 修当前 fast path 的静默截断（先做，立竿见影）

| ID | 内容 | 工作量 | 价值 |
|---|---|---|---|
| **A1** | `target_bits` 动态化：从 m_gen/u 系数上界推导 `2·d·coeff_bits + log2(d!) + 16`，取 max(320, …) | 小（~10 行） | 消除「真平方坐标 >160 bit 时 RRecon 静默失败」 |
| **A2** | 放宽 `nf > 3` 跳过为「combos ≤ fuel.remaining()+1 才进」，配合 fuel 上调到 64 | 小（~5 行） | 覆盖完全分裂素数（Chebotarev density 1/12） |

> **A1 上界出处**（2026-06-30 订正）：`2·d·coeff_bits + ⌈log2(d!)⌉ + 16` **不是任何论文里的定理**，是 self-derived 工程启发式上界，由三条经典结论拼装而成——
> - **(i) 多项式乘积系数增长**：δ² 未规约系数 ≤ d·(2C-bit 乘积) ⇒ `2·coeff_bits` 项 [von zur Gathen & Gerhard, *Modern Computer Algebra* 3e, §3.1]；
> - **(ii) Hadamard 行列式不等式**：规约 mod m_α 的 d×d trace/Sylvester 系统行列式因子 ≤ d!·(max coeff)^d ⇒ `·d` + `⌈log2(d!)⌉` 项 [Hadamard 1893, Bull. Sci. Math.；Cohen, *A Course in Computational Algebraic Number Theory*, Springer 1993, §2.2.4/p.50 —— giac 自己的 `algo.tex` 也引此条]；
> - **(iii) 绝对乘性高度 H(α²)=H(α)²**（精确）：⇒ H(δ)=√H(u)，是「δ 坐标只需 ~C bit 而非 ~2C bit」的结构性理由 [Silverman, *The Arithmetic of Elliptic Curves* 2e, §VIII.5；Bombieri–Gubler, *Heights in Diophantine Geometry*, §1.5]。
>
> **为何不直接用 naive 严谨界**：在 Vandermonde `V_{ki}=σ_k(α)^i` 上跑 Cramer+Hadamard 得 `|a_i| ≤ d^(d/2)·2^(R·d(d-1))·H(δ)/|disc|^(1/2)` ≈ `2^(d²·C)`，即 `target_bits ~ d²·C`——可证但比启发式慢 ~30×（padic_lift 多 2 次 pk 平方、pk 大 ~4×）。**又紧又严谨**的界需显式算 `disc(m_α)` 取回 `−(1/2)·bit(disc)` 项，归入 **Phase B**（number-field arith 模块）。当前启发式：比 naive 严谨界紧、floor 320 保证小用例不回归、由 `sqrt_base_case_large_coord_above_160_bits`（321-bit 坐标）实测覆盖。

**DoD**：
- `quartic_a4_galois_dim_le_12` 仍绿
- 新增大系数 A₄ 测（构造 √Δ 坐标 >160 bit）绿
- 非平方 u 单测耗时 ≤ 4s（fuel=64 worst case）

**不阻塞**：完备算法（Phase B-D）。

### Phase B — (A) 素理想分解 + 指数全偶判定（核心 negative filter）

| ID | 内容 | 工作量 | 价值 |
|---|---|---|---|
| **B1** | 新建 crate 模块 `giac-core::algebra::number_field_arith` | 小 | 模块边界 |
| **B2** | C1 trial division 素因子 | 小 | 基础 |
| **B3** | C2 素理想分解（good prime） | 中 | (A) 主件 |
| **B4** | C3 素理想赋值 + 指数全偶判定 | 中 | (A) 输出 |
| **B5** | C9 接线：`sqrt_fmodule` 在 norm filter 后、`sqrt_base_case` 前调 (A)；不通过直接 `return None` | 小 | 落地 |
| **B6** | 测试：构造非平方但 norm 平方的 u（如 K=ℚ(√2,√3) 中 u=2·3·√6 的某个非平方倍），验证 (A) 抓住 | 中 | 防回归 |

**DoD**：
- (A) 不通过 ⟹ `sqrt_base_case` 不被调用（用计数器或 `#[cfg(test)]` hook 验证）
- 真平方测（A₄ dim-6 + 大系数）仍绿
- paper proof sketch：good prime 路径的 v_𝔭 公式正确性（写到 `.doc/giac-poly-f5-fglm-ideal-valuation-proof.md`）

**不阻塞**：(B) 单位群。Phase B 完成时已有完备的 negative filter（比 norm 更强），消除大部分 false-None；只是无 positive 完备性保证。

### Phase C — (B) 单位群 mod-2 求解（完备性闭环）

| ID | 内容 | 工作量 | 价值 |
|---|---|---|---|
| **C1** | C5 K 的 signature | 小 | (B) 预备 |
| **C2** | 朴素 LLL（d ≤ POLY_ROOTS_DIM_HARD） | 中大 | (B) 工具 |
| **C3** | C7 基本单位系求取 | 大 | (B) 核心 |
| **C4** | C8 单位部分 mod-2 线性求解 | 小 | (B) 判定 |
| **C5** | C9 接线：(A) ∧ (B) 联合判定；通过才进 p-adic 构造 | 小 | 落地 |
| **C6** | C4 bad prime ramification 补全 (A) | 中大 | (A) 完整 |
| **C7** | 完备性测试：随机生成 K（d ≤ 8），随机 u ∈ K，对比完备算法 vs brute-force（Gröbner 求解 x²=u in K） | 中 | 完备性证据 |

**DoD**：
- 完备性测试 N=1000 随机用例 100% 一致
- paper proof sketch：(A) ∧ (B) ⟺ u 平方的 Hasse 论证写到 `.doc/giac-poly-f5-fglm-hasse-proof.md`
- 任意非平方 u 进入 `sqrt_fmodule` 一定在 (A)+(B) 阶段 None，不进 fuel 搜索（计数器验证）

### Phase D — 形式化（可选，独立排期）

Lean / Mathlib 形式化：
- (A) 公式正确性（good prime 路径）
- Hasse 局部-全局在 `IsSquare` 上的刻画
- 接到 Mathlib 现有的 `NumberField` / `Padic` 发展中

**不在本 issue 范围**。立项时建 `GIAC-poly-f5-fglm-lean-formalization.md`。

---

## 4. 落地建议

### 4.1 推荐执行顺序

```
Phase A（先做，~1 天）
  ├─ A1 target_bits 动态化        ← 立刻消除静默 RRecon 失败
  └─ A2 nf>3 放宽 + fuel=64       ← 覆盖完全分裂素数

Phase B（核心 negative filter，~1 周）
  ├─ B1-B4 素理想分解 + 指数全偶
  └─ B5-B6 接线 + 测试

Phase C（完备性闭环，~3-4 周）
  ├─ C1-C5 (B) 单位群路径
  ├─ C6 bad prime 补全
  └─ C7 完备性测试

Phase D（形式化，独立）
```

### 4.2 与现有 issue 的关系

- **不取代** [GIAC-poly-f5-fglm-fullchain-fuel-audit](GIAC-poly-f5-fglm-fullchain-fuel-audit.md)：fuel 审计是「封顶失败路径代价」，本 issue 是「替换算法使失败路径不发生」。两者可并行——fuel 在 (A)+(B) 漏判时（实现 bug）兜底。
- **取代** [postmortem](GIAC-poly-f5-fglm-debug-postmortem.md) §3.3 的「local Fuel」精神：fuel 是临时工程封顶，本 issue 是它的最终替代品。
- **依赖** [GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md) 的 P2 双边 newtype（已落地）：素理想分解需要严格的 HighFirst/LowFirst 区分，否则 norm/赋值算错（postmortem #5/#7 同款坑）。

### 4.3 何时可以删 `fuel`

- Phase B 完成：可保留 fuel 作为「实现 bug 兜底」，但 (A) 不通过的 u 不会进 fuel。
- Phase C 完成：可降 fuel 上限到 4（仅作 implementation safety net），或彻底删 fuel 改回完整 60 素数 × 2ⁿᶠ 枚举（此时 (A)∧(B) 已保证只有真平方进搜索，fuel 截断不再需要）。

### 4.4 风险登记

| 风险 | 缓解 |
|---|---|
| 素理想赋值公式在 bad prime 错 | Phase B 先只走 good prime，bad prime 退回当前 fast path；Phase C C6 补全 |
| 基本单位系求取数值不稳 | 用精确 BigInt 算术 + LLL；d ≤ POLY_ROOTS_DIM_HARD 上界可控 |
| 完备性测试 oracle 缺 | 用 giac C++ 上游输出对比（golden）+ Mathematica/PARI 交叉验证 |
| 工程量超预算 | Phase A+B 已是「major win」（消除静默 bug + 更强 negative filter）；Phase C 可推迟 |

---

## 5. 参考

### 5.1 理论参考

| 主题 | 出处 |
|---|---|
| Hasse 局部-全局原理（一般 global field） | Keith Conrad, *The Local-Global Principle* — https://kconrad.math.uconn.edu/blurbs/gradnumthy/localglobal.pdf （最易读入门） |
| Hasse 原理的数域版（每完成处局部 ⟺ 全局） | Neukirch, *Algebraic Number Theory* III.6 |
| 素理想分解 / Dedekind 算法 | Cohen, *A Course in Computational Algebraic Number Theory* §6.2 |
| 单位群 + 基本单位系 + LLL | Cohen §6.5；de Weger, *Algorithms for Diophantine Equations* |
| 类群 / bnfinit 算法背景 | PARI tutorial *Computing class groups and class fields* — https://pari.math.u-bordeaux.fr/Events/PARI2026/talks/bnfinit.pdf |
| giac-rs 现有数学背景 | [GIAC-poly-f5-fglm-math-background](GIAC-poly-f5-fglm-math-background.md) |
| postmortem 上下文 | [GIAC-poly-f5-fglm-debug-postmortem](GIAC-poly-f5-fglm-debug-postmortem.md) §3.3 |

### 5.2 已知实现（参考 / 交叉验证 oracle）

| 实现 | API | 算法路径 | 用途 |
|---|---|---|---|
| **PARI/GP** | `nfeltissquare(nf, x, &y)` | Hasse 完备：`bnfinit` 算类群 + 基本单位系，`idealprimedec` + `nfeltval` 做素理想赋值，单位部分 `bnfissunit` mod-2 求解 | **首选 oracle**，与 Phase C 完备算法对齐；含 bad prime（ramified）路径 |
| **PARI/GP** | `nfislocalpower(K, P, x, n)` | 单素理想局部判定（绕过 bnfinit 重算） | Phase B 单素理想赋值单元测试 oracle |
| **Sage** | `NumberFieldElement.is_square(root=True)` → `(bool, sqrt)` | 包装 PARI；`K.ideal(x).factor()` + `K.unit_group()` | Python 接口跑大规模随机测试（Phase C C7） |
| **Magma** | `IsSquare(x)` for `RngOrdElt` | Hasse + 理想论（同 PARI 路径） | 第三方交叉验证（学术机构许可） |
| **SymPy** | `AlgebraicNumber` **无 `is_square`** | 退路：在 `QQ[θ]` 上分解 `x²−u`（多项式因子分解路径，非 Hasse） | 仅作小 d 对照；不靠它做完备性 oracle（issue #25087 明确未实现） |
| **giac C++ 上游** | `solve.cc` / `usual.cc` 中 `is_square` 等价物 | **未实现完备算法**，与 giac-rs 当前 fast path 同源 | golden 等价对照（项目原则：先达 giac 等价再扩展） |

**关键结论**：完备算法的成熟实现是 **PARI/GP（`bnfinit` + `nfeltissquare`）**。giac-rs Phase C 实现完成后，应以 PARI 作为完备性 oracle（`gp -q` 跑 `nfeltissquare(bnfinit(f), u)` 对比 Rust 输出）。

#### 5.2.1 PARI/GP 性能 benchmark（实测，2026-06-29，`/home/kanli.hu/upstream/pari/gp` v2.18.1 dev 30935）

| 操作 | d | 用时 | 备注 |
|---|---|---|---|
| `nfinit(f)`（基础 nf 结构） | 4 | ~0 ms | 不含类群 / 单位群 |
| `bnfinit(f)`（完备 bnf 结构） | 8 | ~40 ms | 含类群 + 基本单位系 |
| `nfeltissquare(nf, x)`（基础路径） | 4 | ~0 ms | 仅 norm + 局部 |
| `nfeltissquare(bnf, x)`（完备路径） | 8 | ~1 ms | 完整 Hasse 局部-全局 |

**对比 giac-rs 现状**：真平方 fast path 命中 ~120 ms，非平方 fuel=8 截断 ~2 s。PARI 完备判定 **比 giac-rs fast path 快 2 个数量级**——完备算法在 d ≤ 8 范围内**比启发式 fast path 更快**，因为不需要 fuel 截断 + 重复 p-adic 提升尝试，直接走 Hasse 一次性判定。

**giac-rs Phase C 目标**：d ≤ `POLY_ROOTS_DIM_HARD`（=12）单次完备判定 ≤ 50 ms，比 PARI 慢 ≤ 50× 可接受（Rust 无 bnfinit 缓存 + 自写朴素 LLL）。

**API 正确性验证**（全部通过）：

| 域 | 元素 | 预期 | 实测 |
|---|---|---|---|
| ℚ(√2) | 2, 8, 9 | 1, 1, 1 | ✅ |
| ℚ(√2) | 3 | 0 | ✅ |
| ℚ(∛2) | α, α² | 0, 1 | ✅（√(α²)=α） |
| ℚ(√2,√3) | 2, 3, 6 | 1, 1, 1 | ✅（√6=√2·√3） |
| ℚ(√2,√3) | 5 | 0 | ✅ |
| ℚ(√2,√3,√5) d=8 | 30 | 1 | ✅（√30=√2√3√5） |
| ℚ(√2,√3,√5) d=8 | 7 | 0 | ✅ |

**已知坑**：`nffactor(KL, x^2-D)` 报「incorrect priority」——PARI 要求外部多项式用与域生成元**不同**的变量名（`subst(x^2-D, x, y)` 规避）。不影响 `nfeltissquare` 本身。

### 5.3 形式化现状（Phase D 准备）

| 系统 | 状态 |
|---|---|
| **Lean 4 / Mathlib** | 有 `Mathlib.NumberTheory.NumberField.Basic`（数域 + 整数环 + integral basis）、`Mathlib.Tactic.NormNum.IsSquare`（ℕ/ℤ/ℚ 的 `IsSquare`）、adele 环局部紧性已形式化（*Formalising the local compactness of the adele ring*, AFM 2024） |
| **Lean 4 / Mathlib 缺件** | **「数域元素的 Hasse 局部-全局 `IsSquare` 刻画」尚未形式化**；素理想赋值的可计算版本、单位群 Dirichlet 结构、基本单位系算法均未进 Mathlib |
| **Coq / Isabelle** | 无现成数域平方判定形式化（Mathematical Components 库有代数数论基础设施但无此具体定理） |

**Phase D 风险**：不是「已有 Mathlib 定理 + 接线」的工作量，而是「形式化新定理 + 算法」的论文级工作。立项时需明确这点，预期 1-3 个月独立项目，本 issue 的 Phase A-C 不阻塞 Phase D。

### 5.4 giac-rs 内部参考

#### 5.4.1 代码位置索引

| 子件 | 位置 |
|---|---|
| 当前 fast path 实现 | `giac-rs/crates/giac-core/src/algebra/poly_roots.rs:1460` (`sqrt_base_case`) |
| norm 早退（已落地） | `poly_roots.rs:1955-1968` |
| p-adic Newton lift | `poly_roots.rs:1583` (`padic_sqrt_lift`) |
| ℚ 验证 δ²=u | `poly_roots.rs:1559-1567` |
| `factor_mod_irreducibles`（素理想分解可复用） | `giac-poly::factor::mod` |
| `ExtensionField` API | `giac-rs/crates/giac-core/src/algebra/field_arith.rs`、`ext_tower.rs` |

#### 5.4.2 PARI/GP oracle 实测命令（Phase C C7 交叉验证用）

PARI/GP 二进制：`/home/kanli.hu/upstream/pari/gp`（v2.18.1 dev 30935，2026-06-16 编译）

**最小 smoke 测试**：

```bash
/home/kanli.hu/upstream/pari/gp -q <<'EOF'
K = nfinit(x^2-2);
print(nfeltissquare(K, 2));   \\ 1
print(nfeltissquare(K, 3));   \\ 0
{y} = nfeltissquare(K, 2, &y); print(y);   \\ sqrt(2) 的坐标
EOF
```

**完备路径 oracle（Phase C C7 大规模随机测试）**：

```bash
/home/kanli.hu/upstream/pari/gp -q <<'EOF'
\\ 对一个 (K, u) 对判定 + 返回 sqrt
f = x^4 - 10*x^2 + 1;          \\ K = Q(sqrt2, sqrt3)
B = bnfinit(f);
u = 6;                          \\ 6 = 2*3, sqrt(6) = sqrt2*sqrt3 in K
ok = nfeltissquare(B, u, &r);
print(ok, "  sqrt = ", r);     \\ 1  sqrt = ...
EOF
```

**批量 oracle 调用模板（Rust 测试集对照）**：

```bash
\\ 输入：每行 (minpoly, u_coords) 对
\\ 输出：每行 (ok, sqrt_coords) 对，与 Rust `sqrt_fmodule` 输出 diff
/home/kanli.hu/upstream/pari/gp -q <<'EOF'
pairs = [
  [x^4 - 10*x^2 + 1, 6],
  [x^4 - 10*x^2 + 1, 5],       \\ expect 0
  [x^3 - 2, 4],                \\ expect 1, sqrt=2
  [x^3 - 2, 2],                \\ expect 0
  \\ ... Phase C C7 扩展到随机生成 K + 随机 u
];
for (i = 1, #pairs, \
  f = pairs[i][1]; u = pairs[i][2]; \
  B = bnfinit(f); \
  ok = nfeltissquare(B, u, &r); \
  print("f=", f, "  u=", u, "  ->  ok=", ok, "  sqrt=", r); \
);
EOF
```

**性能 benchmark 复现**：

```bash
/home/kanli.hu/upstream/pari/gp -q <<'EOF'
T = gettime(); B = bnfinit(x^8 - 12*x^6 + 23*x^4 - 12*x^2 + 1);
print("bnfinit d=8: ", gettime()-T, " ms");
T = gettime(); print("is_square(30): ", nfeltissquare(B, 30), "  ", gettime()-T, " ms");
EOF
```

**PARI API 速查**（Phase C 实现参照）：

| PARI 函数 | 数学含义 | giac-rs 对应 Phase C 子件 |
|---|---|---|
| `nfinit(f)` | 基础数域结构 | 已有 `ExtensionField` |
| `bnfinit(f)` | +类群 + 基本单位系 | Phase C C2-C7 |
| `nf.disc` / `nf.r1` / `nf.r2` | 判别式 / 实嵌入数 / 复嵌入对数 | Phase C C1 (signature) |
| `idealprimedec(K, p)` | p 在 K 的素理想分解 | Phase B B3 |
| `nfeltval(K, x, P)` | v_P(x) 素理想赋值 | Phase B B4 |
| `bnfissunit(B, x)` | x 在基本单位系 + 挠群下的对数 | Phase C C4 |
| `nfeltissquare(B, x, &r)` | 完备判定 + 构造 sqrt | Phase C 总成 |
| `nfislocalpower(K, P, x, n)` | 单素理想局部判定 | Phase B 单元测试 oracle |
