# GIAC — Phase B (A) 素理想赋值公式证明 sketch

**状态:** Phase B 落地配套（good-prime 完整路径：inert + split-prime Hensel）+ Phase C C4（tame bad-prime Dedekind-Kummer，`p ∤ index` 且 `e_i < p`）；wild / `p | index` / 单位群待 Phase C+
**关联:** [GIAC-poly-f5-fglm-complete-square-decision](GIAC-poly-f5-fglm-complete-square-decision.md) §1.2（条件 (A)）、`giac-rs/crates/giac-core/src/algebra/number_field_arith.rs`
**快照:** 2026-06-30

---

## 0. 目的

为 Phase B 落地的 `ideal_valuations_parity_ok` 提供纸面证明 sketch：**good-prime 路径**（inert + split）的 `v_𝔭(u)` 公式正确性，以及该路径作为 *negative filter* 的 soundness（返回 `false` ⟹ `u` 确非平方）。完备性（`(A) ∧ (B) ⟺ u` 平方）是 Phase C 的事，不在本文。

## 1. 记号与前置

- `K = ℚ(α)`，`m_α ∈ ℚ[t]` 为 `α` 在 ℚ 上的首一极小多项式，`d = deg m_α = [K:ℚ]`。
- `u ∈ K×`，`u = u(t) ∈ ℚ[t]`（在基 `{1,α,…,α^{d−1}}` 下的坐标多项式）。
- 有理素数 `p`。**good prime**：`p ∤ disc(m_α)` 且 `p ∤ lc(m_α)`（`m_α` mod `p` 平方自由、次数不降）。**inert prime**：good 且 `m_α mod p` 在 `𝔽_p[t]` 中不可约（单一素理想 `𝔭` 在 `p` 上，剩余次数 `f_𝔭 = d`，分歧指数 `e_𝔭 = 1`）。
- `v_𝔭(u)`：`u` 在 `K` 的素理想 `𝔭` 处的离散赋值。
- `N_{K/ℚ}(u) = ∏_{σ:K↪ℚ̄} σ(u) ∈ ℚ×`：全局范数。

### 1.1 Dedekind 定理（good prime）

对 good prime `p`，`m_α mod p = ∏_i g_i`（`𝔽_p[t]` 中互异首一不可约），且 `𝔬_K ⊗ ℤ_p = ℤ_p[α]`（`p` 不整除指数 `[𝔬_K:ℤ[α]]`，因为 good prime 不整除 `disc(m_α)`，而 `disc(K)·[𝔬_K:ℤ[α]]² = disc(m_α)`）。于是 `p` 在 `𝔬_K` 中的素理想分解

```
(p) = ∏_i 𝔭_i^{e_i},   𝔭_i = (p, g_i(α)),   e_i = 1,   f_i = deg g_i
```

是 *unramified* 的（所有 `e_i = 1`）。inert 即 `r = 1`（单一因子，`f_1 = d`）。

### 1.2 范数与赋值的关系

对任意 `u ∈ K×`，

```
v_p(N_{K/ℚ}(u)) = Σ_{𝔭 | p} f_𝔭 · v_𝔭(u)                                  …(★)
```

这是局部-全局范数公式的标准推论：`N_{K/ℚ}(u) = ∏_{𝔭|p} N_{K_𝔭/ℚ_p}(u)`（按 `p` 分量），而
`v_p(N_{K_𝔭/ℚ_p}(u)) = f_𝔭 · v_𝔭(u)`（unramified 局部域 `K_𝔭/ℚ_p`，`e=1`，`v_𝔭(p)=1`，
`N_{K_𝔭/ℚ_p}(π_𝔭) = p^{f_𝔭}`，故 `v_ℚ_p(局部范数) = f_𝔭·v_𝔭`）。
[Neukirch, *Algebraic Number Theory* III.6；Cohen, *A Course in Computational ANT* §4.2 / §6.2]

## 2. inert 路径赋值公式（核心）

**命题 (A-inert).** 设 `p` 为 `K` 的 good **inert** prime，`u ∈ K×`。则 `𝔭 = (p, m_α mod p)` 是 `p` 上唯一素理想，`f_𝔭 = d`，且

```
v_𝔭(u) = v_p(N_{K/ℚ}(u)) / d.                                              …(†)
```

**证明.** inert ⟹ `r = 1` ⟹ (★) 右端只有一项：`v_p(N_{K/ℚ}(u)) = f_𝔭·v_𝔭(u) = d·v_𝔭(u)`。解出 `v_𝔭(u) = v_p(N_{K/ℚ}(u))/d`。整数性自动成立：`v_𝔭(u) ∈ ℤ`，故 `d | v_p(N_{K/ℚ}(u))`。 ∎

**实现对应**（`number_field_arith.rs::ideal_valuations_parity_ok`）：对每个小素数 `p | N_{K/ℚ}(u)`，
- 跳过 `p | lc(m_α)`（degree drop = bad）；
- `factor_mod_irreducibles(m_α, p)`：若因子次数和 `≠ d` 或有重复因子 ⟹ bad，跳过；
- `facs.len() == 1`（inert）⟹ `v_𝔭(u) := vp_norm / d`，其中 `vp_norm = v_p(N_{K/ℚ}(u))` 由 `N(u)` 的分子/分母 trial-division 得到（(†) 不需要 `u` 的坐标，只需要全局范数 `N(u)`，因此**无需 Hensel / 局部范数 / resultant**）；
- `v_𝔭(u)` 奇 ⟹ 返回 `false`。

**为什么不需要 `u_int` / 去分母**：(†) 用的是 `u` 的 *全局* 范数 `N_{K/ℚ}(u) ∈ ℚ×`，它已包含 `u` 的分母信息（`v_p(N(u))` 对 `p | denom(u)` 取负值，正确反映 `v_𝔭(u) < 0` 的情形）。所以 `u` 的有理坐标直接经 `krylov_minpoly_coords` → `N(u) = (−1)^d·c0^{d/d_f}` 即可，无需额外整数化。这是 inert 路径比 split 路径轻得多的原因。

## 2a. split 路径赋值公式（Hensel + 局部范数）

**命题 (A-split).** 设 `p` 为 `K` 的 good **split** prime，`m_α mod p = ∏_{i=1}^r g_i`（`r > 1`，`g_i` 互异首一不可约，`f_i = deg g_i`，`Σ f_i = d`）。对每个 `i`，记 `u_int = d_u·u ∈ ℤ[α]`（`d_u = lcm(u 的分母)`）。则

```
v_{𝔭_i}(u) = v_p( N_{K_{𝔭_i}/ℚ_p}(u_int) ) / f_i  −  v_p(d_u),          …(‡)
```

其中 `N_{K_{𝔭_i}/ℚ_p}(u_int) = det( 乘 u_int 在 (ℤ/p^N)[t]/(G_i) )`，`G_i` 是 `g_i` 的 Hensel 提升到 `mod p^N` 的首一因子（`G_i | m_α mod p^N`）。精度 `N = v_p(N_{K/ℚ}(u_int)) + 1` 足够（见下）。

**证明.**
1. **Dedekind**（§1.1）：`𝔭_i = (p, g_i(α))`，`e_i = 1`，`f_i = deg g_i`。`u_int ∈ ℤ[α] ⊂ 𝔬_K`（整数），故 `v_{𝔭_i}(u_int) ≥ 0`。
2. **局部范数 = resultant**：`K_{𝔭_i} ≅ ℚ_p[t]/(G_i)`（`G_i` 是 `g_i` 的 Hensel 提升，`m_α` 在 `ℚ_p` 上分裂为 `∏ G_i`，`G_i` 不可约）。`N_{K_{𝔭_i}/ℚ_p}(u_int) = ∏_{G_i(ρ)=0} u_int(ρ) = Res(G_i, u_int)`（`G_i` 首一 ⟹ `lc(G_i)^{deg u_int} = 1`）。实现用 `det(乘 u_int 的 f_i×f_i 矩阵)`（= resultant，Bareiss over ℤ 再 `mod p^N`）。
3. **赋值公式**：unramified 局部域 `K_{𝔭_i}/ℚ_p`，`v_{𝔭_i}(u_int) = v_p(N_{K_{𝔭_i}/ℚ_p}(u_int))/f_i`（同 (★) 的局部版，`e=1`）。
4. **去分母**：`u = u_int/d_u` ⟹ `v_{𝔭_i}(u) = v_{𝔭_i}(u_int) − v_{𝔭_i}(d_u) = v_{𝔭_i}(u_int) − v_p(d_u)`（`d_u ∈ ℤ`，`v_{𝔭_i}(d_u) = e_i·v_p(d_u) = v_p(d_u)`）。合并 3、4 得 (‡)。
5. **精度足够**：`v_p(N_{K_{𝔭_i}/ℚ_p}(u_int)) = f_i·v_{𝔭_i}(u_int) ≤ Σ_j f_j·v_{𝔭_j}(u_int) = v_p(N_{K/ℚ}(u_int)) =: V`（最后一步由 (★) 用于 `u_int`）。故取 `N = V+1` 时 `v_p(local_norm) < N`，`mod p^N` 下 `v_p` 精确（不溢出）。 ∎

**实现对应**（`number_field_arith.rs::ideal_valuations_parity_ok` split 分支）：
- `u_int_low, d_u = clear_denoms_low(u_low)`；`vp_nu_int = v_p(N(u)) + d·v_p(d_u)`（即 `V`）。
- 若 `vp_nu_int > PREC_CAP(60)` ⟹ 跳过（Hensel 代价过高，sound false-pass）；`N = vp_nu_int + 1`。
- 对每个 `g_i`：`G_i = hensel_lift_factor(m_α, g_i, p, N)`（Zassenhaus 二次 Hensel，`(G,H,A,B)` 四元组，**reduced corrections** `δg = (b·e) rem g`、`δh = (a·e) rem h` 保持 `G,H` 首一；见下注）；`local_norm = local_norm_mod_pk(u_int, G_i, p^N)`（`det(乘 u_int 的 f_i×f_i 矩阵)`，Bareiss over ℤ，`mod p^N`）；`v_{𝔭_i}(u) = v_p(local_norm)/f_i − v_p(d_u)`，奇 ⟹ `false`。

> **Hensel reduced-correction 注**：标准二次 Hensel 步 `g' = g + m·(t·e mod m)` 不对 `t·e` 取 `mod g`，会使 `g'` 的首项系数偏离 1（`t·e` 的 `deg ≥ f` 项污染），强制 monic 则破坏 `g'·h' ≡ m_α mod m²` 不变量。正确做法是 **reduced** Hensel：`δg = (t·e) rem g`（除以首一 `g`，余式 `deg < f`），`g' = g + m·δg` 自动保持首一；witness 同理 `δa = (t·q) rem g`。这是 `hensel_lift_sqrt2_mod_17_cubed` 单测（`√2 mod 17³ = 4290`，PARI `factorpadic` 交叉验证）捕捉到的关键修复。

**PARI/GP oracle 交叉验证**（`/home/kanli.hu/upstream/pari/gp`）：
- `u = 187+102√2 = 17·(3+√2)²`，`N(u) = 119²`（norm 平方，norm filter 放行）。split `p=17`（`2` 是 QR mod 17：`6²=2`，`t²−2=(t−6)(t+6)`）。PARI `nfeltval(K,u,pr[1])=1`、`nfeltval(K,u,pr[2])=1`（均奇），`nfeltissquare=0`。split 路径抓住两个奇赋值 ⟹ (A) `false` ⟹ `sqrt_fmodule` bail，`sqrt_base_case` 不调用（`sqrt_fmodule_b6_split_non_square_skips_base_case` 计数器 DoD）。
- split 真平方 `u = (7(1+√2))² = 147+98√2`，`N(u) = 7⁴`，split `p=7`（`2` QR mod 7）。`v_{𝔭_i}(u) = 2·v_{𝔭_i}(7(1+√2)) = 2·1 = 2`（偶）⟹ (A) 放行（**不误拒真平方**，`ideal_valuations_parity_split_true_square_even_val`）。

## 3. Soundness（negative filter 正确性）

**命题 (sound).** 若 `ideal_valuations_parity_ok` 返回 `false`，则 `u` 在 `K` 中**不是**平方。

**证明.** 返回 `false` 当且仅当存在某 good inert prime `p | N(u)` 使 `v_𝔭(u)` 奇（由 (†)）。若 `u = δ²`（`δ ∈ K×`），则 `v_𝔭(u) = v_𝔭(δ²) = 2·v_𝔭(δ)` 恒为偶——矛盾。故 `u` 非平方。 ∎

**关键**：`v_𝔭(δ²) = 2·v_𝔭(δ)` 是赋值的 *加性* 与 *整性* 的直接推论，对任意 `𝔭`（含 ramified）都成立。因此 (A) 的「全偶」是 `u` 为平方的**必要条件**；其否定（存在奇赋值）就是 `u` 非平方的**充分**证书。这是 (A) 作为 negative filter 永不 false-reject 的根源。

## 4. 比 norm filter 更强的具体刻画（even-degree 场景）

norm filter 检查 `N_{K/ℚ}(u)` 是否 ℚ-平方，即 `v_p(N(u))` 对所有 `p` 为偶。对 `d` 偶的 inert prime `p`：

```
v_p(N(u)) 偶  ⇏  v_𝔭(u) 偶
```

反例：`v_p(N(u)) = d·v_𝔭(u)`。`d` 偶时 `v_p(N(u))` 偶 ⟺ `v_𝔭(u)·d` 偶 ⟺ `v_𝔭(u)` 任意（因 `d` 已贡献因子 2）。具体地 `v_p(N(u)) ≡ d (mod 2d)`（即 `v_𝔭(u) = 1`）时 `v_p(N(u))` 偶但 `v_𝔭(u)` 奇——norm filter 放行，(A) 抓住。

**B6 测试**就是此情形：`K = ℚ(√2)`（`d=2`），`u = 3`，`N(u) = 9 = 3²`，inert good prime `p = 3`（`2` 是 `mod 3` 非剩余 ⟹ `t²−2` 在 `𝔽_3` 不可约）。`v_3(9) = 2`，`v_𝔭(3) = 2/2 = 1` 奇 ⟹ (A) 返回 `false`，`sqrt_fmodule` 在 `sqrt_base_case` 之前 bail（DoD：计数器验证 `sqrt_base_case` 调用数不变）。norm filter 因 `9` 是 ℚ-平方而放行——正是 (A) 比 norm 更强的体现。`√3 ∉ ℚ(√2)`（`3 = a²+2b²+2ab√2` 无 ℚ 解），故 `3` 确为非平方，soundness 保真。

## 5. 不完备性（为何 Phase B 只是 negative filter）

Phase B 的 (A) **不完备**，有二处 sound-but-not-complete 的跳过：

1. **wild prime**（`(m_α mod p)' = 0`，即 `m_α mod p` 是 p-次幂；或更一般地某 `e_i ≥ p` wild-ish ramification）：`factor_mod_irreducibles` 的 Yun squarefree 分解在 char-p 不能做 p-次根抽取，给错因子（会 false-reject 真平方）。`factor_with_multiplicities_tame` 检测 `(m_α mod p)' = 0` 或恢复出的 `e_i ≥ p` 即返回 `None` → 跳过（sound false-pass）。对 `f=1` wild prime，全局 norm filter 已能抓 `v_𝔭` 奇（`v_p(N) = v_𝔭`），跳过不损完备性；`f` 偶的 wild prime 需 char-p p-次根技巧 + Montes/Round 4，延后至 Phase C+。
2. **`p | [𝔬_K : ℤ[α]]`（index）的 bad prime**：Dedekind 单步分解不适用（`g_i` 不直接给素理想结构），`dedekind_index_ok` 判据检测 `p | index` 即跳过（sound false-pass；需 Montes/Round 4，延后）。
3. **大素数 / 高赋值素数**：只扫 `small_primes(60)` 中整除 `N(u)` 的素数；且对 `v_p(N(u_int)) > PREC_CAP(60)` 的素数跳过（Hensel 代价过高）。`N(u)` 的 *大* 素因子不被 trial-division 发现，或 Hensel 精度超 cap ⟹ 跳过。对大系数 `u`（如 A1 测试 `N(u) ~ 2^{13000}`），大素因子无法分解或超 cap，(A) 退化为 vacuous-pass，回退到现有 `sqrt_base_case` fast path（由 ℚ 验证保 soundness）。

> **tame bad-prime 已补全（Phase C C4，2026-06-30）**：`p ∤ index` 且所有 `e_i < p`（tame ramification）的 bad prime 现由 Dedekind-Kummer 路径处理 —— `factor_with_multiplicities_tame` 用 `b = f/gcd(f,f')` 取平方自由部分（避免对非平方自由输入调用 `factor_mod_irreducibles` 而挂死），重复除法恢复 `e_i`；`dedekind_index_ok` 判 `p ∤ index`；ramified Hensel 提升 `g_i^{e_i} → G_i`（首一，次数 `e_i·f_i`，`ℚ_p`-不可约），局部范数 `det(乘 u_int mod G_i)`，`v_{𝔭_i}(u) = v_p(local_norm)/f_i − v_p(d_u)`。r=1 ramified-single 用 fast path `v_𝔭 = v_p(N)/f_1`（对任意 `e_1` 成立）。PARI 交叉验证 `m=t⁴+2t²+4` at `p=3`（`e=2, f=2`）：`v_p3(α²+1)=1` 奇、`N=9` ℚ-平方 ⟹ (A) 拒；`−3=(α²+1)²` 真平方 ⟹ (A) 放行。

二处跳过都是 *negative filter 的安全失败*：跳过一个可能失败的判定只能导致「本可 bail 却进入搜索」（false-pass，浪费 fuel），不会导致「误拒真平方」（false-reject，破坏 soundness）。最终的 ℚ 等式验证 `δ²=u`（`poly_roots.rs:1559`）是兜底 certificate，保证整体 soundness 不被任何 (A) 缺口破坏。

> **split-prime 已补全**（2026-06-30）：原列于 Phase C C6 的 split-prime Hensel 路径已在 Phase B 落地（§2a），不再是缺口。

## 6. 真平方不被误拒的验证

- inert 真平方 `ideal_valuations_parity_true_square_inert_even`：`K=ℚ(√2)`，`δ = 9+6√2`（`N(δ)=9`），`u = δ² = 153+108√2`，`N(u) = 81 = 3⁴`。inert `p=3`：`v_𝔭(u) = v_3(81)/2 = 2` 偶 ⟹ (A) 通过（与 `v_𝔭(u)=2·v_𝔭(δ)=2·1` 一致）。
- split 真平方 `ideal_valuations_parity_split_true_square_even_val`：`u=(7(1+√2))²=147+98√2`，split `p=7`，`v_{𝔭_i}=2` 偶 ⟹ (A) 放行。
- 单位 `u=(1+√2)²`（`N=1`）：无素数整除 ⟹ vacuous-pass。
- A₄ dim-12 tower 真平方测（`quartic_a4_galois_dim_le_12`）因 `generator_minpoly_low() = None`（tower）而整段跳过 (A)，保持原 fast path 绿。

## 7. 待办（Phase C）

- **C6 bad-prime Kummer**：**tame 部分（`p ∤ index`，所有 `e_i < p`）已落地 2026-06-30**（见 §5 注）。剩余：wild prime（`(m_α mod p)'=0` p-次幂，需 char-p p-次根技巧）、wild-ish（`e_i ≥ p`）、`p | index`（需 Montes/Round 4）。落地 wild 前需先修复 `giac-poly::factor_mod_irreducibles` 的 Yun 在 tame perfect square 上的挂死（当前用 `b = f/gcd(f,f')` 平方自由预处理绕过；wild 需 p-次根，Yun 本身不支持）。ramified-multi（`r>1` 且某 `e_i>1`）路径已实现但缺直接单测（r=1 ramified + split r>1 已覆盖）。
- **大素因子 / 高精度 cap**：trial-division 上界外的 `N(u)` 素因子（需 `N(u)` 的完整素分解或 ECM）；`v_p(N(u_int)) > 60` 的素数（提高 cap 或改用迭代赋值法降低 Hensel 精度需求）。
- **(B) 单位群**：即使 (A) 全完备，单位部分 `η = u/∏𝔭^{v_𝔭(u)/2}` 是否在 `𝔬_K×/(𝔬_K×)²` 平凡仍需 Dirichlet 单位群 + LLL（C5-C8）。
- 完备性测试 `N=1000` 随机用例 vs PARI `nfeltissquare(bnfinit(f), u)`（C7）。

## 8. 参考

- Keith Conrad, *The Local-Global Principle* — https://kconrad.math.uconn.edu/blurbs/gradnumthy/localglobal.pdf
- Neukirch, *Algebraic Number Theory* III.6（范数-赋值公式、Dedekind 分解）
- Cohen, *A Course in Computational Algebraic Number Theory* §6.2（Dedekind / 素理想分解算法）、§4.2（范数与赋值）
- von zur Gathen & Gerhard, *Modern Computer Algebra* 3e §15（Hensel lifting，Phase C 用）
