# GIAC-poly — Hasse-Scholz 平方判定管线的 Lean 4 验证方案（未来）

**状态:** draft / 长期
**类型:** 形式化 / 架构
**前置:** [giac-poly-f5-fglm-hasse-proof.md](../giac-poly-f5-fglm-hasse-proof.md) Step 1-7a + 7b(部分) + 端到端接线 已落地（认证域 r=0/r=1 规格冻结后再证「实现≈规格」）
**相关:** [GIAC-poly-quartic-lean4-verification.md](GIAC-poly-quartic-lean4-verification.md)（同套三层证书 / fixture registry 约定）、[GIAC-poly-f5-fglm-ideal-valuation-proof](giac-poly-f5-fglm-ideal-valuation-proof.md)（(A) 赋值公式证明）、[GIAC-poly-f5-fglm-complete-square-decision](GIAC-poly-f5-fglm-complete-square-decision.md)
**工具:** Lean 4 + Mathlib4（主）；Rust CI **不阻塞** on `lake build`（独立 job）

---

## 1. 目标与非目标

### 1.1 要证什么

在 **代数数论 / Dedekind 域** 层证明下列命题（名称与 Rust 模块 `class_group.rs::hasse_sqrt` 的出口对齐，便于交叉引用）：

| 层级 | 命题（示意） | 对应 Rust 出口 |
|------|----------------|-----------|
| **S0** | **完备性（Hasse 原理）**：`u ∈ (K×)² ⟺ (A_fin)∧(A_inf)∧(C)∧(B)` | 管线数学根基（`hasse_sqrt` 的存在依据） |
| **S1-true** | `Some((true, δ))` ⟹ `δ² = u`（精确） | 出口 TRUE ①（norm=0）/ ②（末尾 `δ=γ·√η`） |
| **S1-false** | `Some((false, _))` ⟹ `¬IsSquare u` | 出口 FALSE ①②③④（4 个 obstruction） |
| **O1** | ① 奇 `v_𝔭(u)` ⟹ 非平方（`v_𝔭(δ²)=2·v_𝔭(δ)` 偶） | `(A_fin)` `!parity_ok` |
| **O2** | ② 某实嵌入 `σ(u)<0` ⟹ 非平方（`σ(δ²)=σ(δ)²≥0`） | `(A_inf)` `!is_totally_positive` |
| **O3** | ③ `J` 非主 ⟹ 非平方（`u=δ² ⟹ J=(δ)` 主） | `(C)` `!is_principal` |
| **O4** | ④ `η=u/γ²` 非单位平方 ⟹ 非平方（`u=δ²∧(γ)=J ⟹ δ/γ∈𝔬_K× ∧ (δ/γ)²=η`） | `(B)` `!eta_is_sq` |
| **C-suff** | `ideal_is_principal` 找到 `γ` ⟹ `(γ)=J`（`|N(γ)|=N(J)∧γ∈J ⟹ 等理想`） | `(C)` `Some((true,Some(γ)))` |
| **C-exh** | 虚二次定域：正定范型特征值界 `√(N(J)/λ_min)≤界` ⟹ 盒内无 `γ` = 全空间无 `γ` | `imag_quad_search_exhaustive` → `Some((false,None))` |
| **B-recon** | `unit_sqrt` 重建 `√η=ζ'·∏ε_i^{a_i/2}` ⟹ `(√η)²=η` | `(B)` `Some((true,Some(√η)))` |
| **N-norm** | `norm_of`（Krylov）= `N_{K/ℚ}(x)` | `norm_of` |
| **N-val** | 好素数 `v_𝔭(u) = v_p(N(u))/f_𝔭` | `compute_ideal_valuations` |

**soundness 总声明（验证的真正目标）：**

```lean
-- TRUE 出口: 返回 δ ⟹ δ² = u  (永不错)
theorem hasse_sqrt_true_sound :
  hasse_sqrt K u = some (true, some δ) → δ * δ = u

-- FALSE 出口: 返回 false ⟹ u 不是平方  (永不错)
theorem hasse_sqrt_false_sound :
  hasse_sqrt K u = some (false, none) → ¬ IsSquare u
```

`None` 出口**无需验证**——不声称任何事，soundness 平凡成立。Lean 只验**有声称的出口**。

### 1.2 不证什么（第一版）

| 不做 | 原因 | 替代 |
|------|------|------|
| `recover_free_exponents` 的 **f64 log-map** | 数值启发式，无 exact 语义 | 当 **oracle**；靠 `ζ^|μ(K)|=1` exact 检查当证书 |
| `ideal_is_principal` 的 **盒内枚举搜索** | bounded search，不完备 | 当 **oracle**；找到 `γ` 后靠 `|N(γ)|=N(J)∧γ∈J` exact 检查当证书 |
| `fundamental_units` **r≥2 全实 index-1 证书** | 运行时证不出（Friedman/regulator 下界太弱） | **缺口**——Lean 证 r=0/r=1 域 sound，r≥2 Rust 返 `None` 是诚实反映 |
| 复域 rank≥1 torsion / 复 log-map | 未实现 | 范围外（Step 5 剩余） |
| `(A_fin)` 坏素数 / `p∣index` / 大素数（P3）/ `p∣denom` 分支 | 这些分支**跳过**（`!all_processed → None`） | Lean 只验 good/tame/wild 精确分支；跳过分支不声称 |
| Rust `ext_tower` / `ExtensionField` **逐行提取** | 成本极高 | shallow spec + 证书检查器（§2） |
| 通用「算法终止 / 复杂度」 | 用 Rust release 测 + `N_J_MAX=10⁷` 工程界 | 同 quartic doc §1.2 |

**原则：** Mathlib 证 **数学规格 + 证书检查器**；f64/搜索当 oracle 找证书；Rust 测证 **实现满足规格**；文档用 `Theorem giac.Hasse.O3` 链两者。

---

## 2. 架构：三层证书 + Scholz 证书机制

### 2.1 三层（同 quartic doc §2）

```text
┌─────────────────────────────────────────────────────────┐
│  L-数学层 (Lean 4 / Mathlib)                             │
│  定理 S0, S1-true/false, O1-O4, C-suff, C-exh, B-recon  │
│  + 证书检查器 cert_checker : ScholzCert → Prop           │
└───────────────────────────┬─────────────────────────────┘
                            │ 文档引用定理名 + 假设列表
┌───────────────────────────▼─────────────────────────────┐
│  S-规格层 (Lean 结构体 ScholzCert / Rust doc 镜像)       │
│  ScholzCert：TRUE 带见证 (γ,a,ζ',δ) / FALSE 带 obstruction│
└───────────────────────────┬─────────────────────────────┘
                            │ cert_checker 测试 + δ²=u 自检
┌───────────────────────────▼─────────────────────────────┐
│  I-实现层 (giac-core class_group / unit_group / poly_roots)│
│  hasse_sqrt / ideal_is_principal / unit_sqrt (f64+搜索)  │
│  f64/搜索 = 找证书的 oracle; exact 检查 = 签字           │
└─────────────────────────────────────────────────────────┘
```

### 2.2 证书机制（与代码的核心连接点）

**关键设计：** `hasse_sqrt` 的每个 `Some((true/false))` 出口附带一个 **`ScholzCert` 证书对象**，携带 exact 见证。证书可被一个**独立的、小而稳的 checker** 重验——checker 只做精确 ℚ(α) 算术（无 f64、无搜索），可整体嵌入 Lean 验证。这是 Pocklington 素性证书的传统套路：**搜索找证书，checker 签字**。

```rust
// Rust 侧（未来加；当前 hasse_sqrt 只返 Option<(bool,Option<δ>)>）
pub(crate) enum ScholzCert {
    Square {
        gamma: HighFirstQ,       // (C) 生成元: (γ)=J
        a_vec: Vec<i64>,         // (B) 自由指数 (全偶)
        zeta_prime: HighFirstQ,  // (B) torsion 平方根: ζ'²=ζ
        delta: HighFirstQ,       // 重建 √u = γ·ζ'·∏ε_i^{a_i/2}
    },
    NonSquare(Obstruction),
}

pub(crate) enum Obstruction {
    OddValuation { p: i64, g_i: PolyMod, v_p_u: i64 },            // ①
    NotTotallyPositive { real_emb_index: usize },                 // ② (符号证据)
    NonPrincipalIdeal { j_factors: Vec<JFactor>, bound_evidence },// ③ (虚二次特征值界)
    NonUnitSquare { a_odd_index: Option<usize>, zeta: HighFirstQ },// ④
}
```

**Lean checker（要证 sound）：**

```lean
-- TRUE 证书: 检查 γ∈J ∧ |N(γ)|=N(J) ∧ all_even a ∧ ζ'²=ζ ∧ δ=γ·ζ'·∏ε_i^(a/2) ⟹ δ²=u
theorem cert_square_sound (c : ScholzCert.square) :
  cert_square_check K u c = true → δ_of c * δ_of c = u

-- FALSE 证书: 每种 obstruction 检查通过 ⟹ ¬IsSquare
theorem cert_nonsquare_sound (c : Obstruction) :
  cert_nonsquare_check K u c = true → ¬ IsSquare u
```

**为什么这样设计连接代码：**
- f64 `recover_free_exponents` 找的 `a_vec` 可能因舍入错 —— 但 checker 验 `η·∏ε_i^(-a_i) = ζ ∧ ζ^|μ(K)|=1`，这等价于"`a_i` 是真指数"（前提：单位基 full）。错舍入 ⟹ checker 拒。**f64 不进 Lean，证书进 Lean。**
- 盒搜索找不到 `γ` —— 在虚二次域，checker 验"特征值界 ≤ 枚举界 ∧ 盒内无 `γ`"⟹ `C-exh` 定理 ⟹ 非主。在不定域/高次，**不发 FALSE 证书**（Rust 返 `None`，无证书，无声称）。
- `δ²=u` 自检（已在端到端测试里）就是 `cert_square_sound` 的运行时实例。

---

## 3. 仓库与目录（建议）

复用 [GIAC-poly-quartic-lean4-verification.md](GIAC-poly-quartic-lean4-verification.md) §3 的 `giac-proofs/` Lake 项目，新增 `Hasse/` 子树：

```text
giac-proofs/
  Giac/
    Root/            # 已有：四次路径 (Q0-Q4, T1-T2)
    Tower/           # 已有
    Hasse/
      Basic.lean          # 数域 K=ℚ(α)、范数 N_{K/ℚ}、素理想赋值 v_𝔭
      Valuations.lean     # N-val: 好素数 v_𝔭 = v_p(N)/f_𝔭
      Obstructions.lean   # O1-O4: 4 个平方 obstruction 引理
      Principal.lean      # C-suff: |N(γ)|=N(J)∧γ∈J ⟹ (γ)=J
      EigenBound.lean     # C-exh: 正定二次型特征值界 (虚二次非主认证)
      UnitSquare.lean     # B-recon: √η=ζ'·∏ε_i^(a/2) 重建
      ScholzIff.lean      # S0: Hasse 原理 (A_fin)∧(A_inf)∧(C)∧(B) ⟺ IsSquare
      Cert.lean           # ScholzCert 结构 + cert_square/nonsquare_check + sound
    Examples/
      QSqrt2.lean         # ℚ(√2) u=2 / u=6+4√2 (TRUE 证书实例)
      QI.lean             # ℚ(i) u=-1 (TRUE，torsion 平方)
      QSqrtNeg5.lean      # ℚ(√-5) u=2 (FALSE ③ 非主) / u=4 (TRUE γ=±2)
      QSqrtNeg6.lean      # ℚ(√-6) u=2 (FALSE ③)
      QSqrt3.lean         # ℚ(√3) u=3 / u=6-3√3 (FALSE ④ η=ε 非单位平方)
  README.md               # 定理索引 (含 Hasse 定理表)
```

**依赖：** `require mathlib from git`；不 require Rust。复用 quartic 项目的 `lean-toolchain` pin。

---

## 4. 数学形式化路线（分阶段）

### Phase H0 — 数域基础设施 + 范数（2-3 周）

**Mathlib 已有：** `NumberField`, `IsDedekindDomain`, `Ideal`, `norm`, `MinimalPolynomial`, `Polynomial.roots`.

**新建（`Giac.Hasse.Basic` / `Valuations`）：**

```lean
-- K = ℚ(α), α 的 minpoly = m_α
variables (K : Type*) [Field K] [NumberField K] [FiniteDimensional ℚ K]
variable (α : K) (hα : IsIntegral ℚ α) (hmin : minpoly ℚ α = m_α)

-- N-norm: Krylov minpoly 常数项 → 范数
lemma norm_of_eq (x : K) : norm_of x = N_{K/ℚ} x
-- N-val: 好素数 (p∤disc, p∤index) 赋值公式
lemma val_good_prime (p : ℕ) (hp : GoodPrime p m_α) (𝔭 : Ideal 𝔬_K) (h𝔭 : 𝔭 ∣ p) :
  v_𝔭 x = (v_p (N_{K/ℚ} x)) / f_𝔭
```

**与 Rust：** `norm_of` (class_group.rs) / `compute_ideal_valuations` (number_field_arith.rs) good-prime 分支。

### Phase H1 — 4 个 obstruction 引理（O1-O4）（3-4 周，核心数学）

**这是整个验证的数学核心——每条都是"平方的必要条件 + 逆否"，短而 Mathlib 友好：**

```lean
-- O1: 奇赋值
theorem not_square_of_odd_val (𝔭 : Ideal 𝔬_K) (h : Odd (v_𝔭 u)) : ¬ IsSquare u := by
  -- δ² ⟹ v_𝔭(δ²) = 2·v_𝔭(δ) is Even; 逆否

-- O2: 非全正
theorem not_square_of_not_tot_pos (σ : K →ₐ[ℚ] ℝ) (h : σ u < 0) : ¬ IsSquare u := by
  -- δ² ⟹ σ(δ²) = σ(δ)² ≥ 0; 逆否

-- O3: J 非主
theorem not_square_of_non_principal (h : ¬ IsPrincipal J) : ¬ IsSquare u := by
  -- δ² ⟹ J = ∏𝔭^(v_𝔭(δ)) = (δ) principal; 逆否

-- O4: η 非单位平方
theorem not_square_of_eta_non_square (γ : K×) (hγ : (γ) = J) (h : ¬ IsSquare η 𝔬_K×) :
    ¬ IsSquare u := by
  -- δ² ∧ (γ)=J ⟹ u/γ² = (δ/γ)² 且 δ/γ ∈ 𝔬_K× (同理想) ⟹ η 单位平方; 逆否
```

**Mathlib 支撑：** `IsSquare`, `Ideal.IsPrincipal`, `Algebra.norm`, `Ideal.valuation`（若 Mathlib 版本足够；否则自写 Dedekind 赋值）。O3 需"平方的理想是主理想"引理（`Ideal.principal_of_square`）。

### Phase H2 — (C) 主性 + 虚二次非主认证（C-suff, C-exh）（4-6 周）

```lean
-- C-suff: 找到 γ ⟹ J 主
theorem principal_of_gen_membership (γ : 𝔬_K) (hN : |N γ| = N J) (hmem : γ ∈ J) :
    IsPrincipal J := by
  -- (γ) ⊆ J (成员) ∧ [(γ):J] = |N γ|/N J = 1 (Dedekind: 指数=范数比) ⟹ (γ)=J

-- C-exh: 虚二次正定范型特征值界
theorem pdqf_bound (Q : QuadraticForm ℤ) (hQ : Q.PositiveDefinite) (N_J : ℕ) :
    Q (a,b) = N_J → a^2 + b^2 ≤ N_J / λ_min Q
-- 推论: 盒内 (|a|,|b| ≤ B) 无 γ ∧ √(N_J/λ_min) ≤ B ⟹ 全空间无 γ ⟹ J 非主
theorem non_principal_of_empty_box (B : ℤ) (hB : √(N_J/λ_min) ≤ B)
    (hempty : ∀ a b ∈ [-B,B], Q(a,b)=N_J → ¬ member J) : ¬ IsPrincipal J
```

**Mathlib 支撑：** `QuadraticForm`, `PositiveDefinite`, 矩阵特征值（`Matrix.eigenvalues` 实对称）。C-exh 是 `imag_quad_search_exhaustive`（class_group.rs）的数学根基。

### Phase H3 — (B) 单位平方重建（B-recon）（3-4 周）

```lean
-- Dirichlet 单位分解 η = ζ·∏ε_i^a_i (假设 fundamental_units 是 full 基 — index-1 证书)
variables (ε : Fin r → 𝔬_K×) (hfull : UnitBasisFull ε)  -- ← 缺口 (§1.2, §10)
variable (ζ : μ K) (a : Fin r → ℤ)

-- B-recon: √η = ζ'·∏ε_i^(a_i/2) 重建
theorem unit_sqrt_recon (hζ' : ζ'^2 = ζ) (heven : ∀ i, Even (a i))
    (hη : η = ζ * ∏ i, ε i ^ a i) :
    (√η)^2 = η := by
  -- √η = ζ'·∏ε_i^(a/2), (√η)² = ζ'²·∏ε_i^a = ζ·∏ε_i^a = η

-- torsion 平方: ζ'² = ζ 的枚举检查
theorem torsion_sqrt_correct (ζ' : μ K) : torsion_sqrt ζ = some ζ' → ζ'^2 = ζ
```

**关键：** `hfull : UnitBasisFull ε` 是**前提**——Lean 证的是"**若**单位基 full，**则**重建 sound"。Rust r=0/r=1 域满足（torsion-only / 单单位 + index-1 由 maximality 保证）；r≥2 缺此前提 ⟹ Rust 返 `None` ⟹ 不发证书 ⟹ 无声称 ⟹ 无需证。

### Phase H4 — Hasse 原理 + 证书检查器（S0, Cert）（4-6 周）

```lean
-- S0: Hasse 原理 (核心定理, 综合 O1-O4 + C-suff + B-recon)
theorem scholz_iff : IsSquare u ↔
    (∀ 𝔭, Even (v_𝔭 u)) ∧ TotallyPositive u ∧
    IsPrincipal (∏ 𝔭 ^ (v_𝔭 u / 2)) ∧
    IsSquare (u / γ^2) 𝔬_K×

-- Cert: 证书检查器 sound
theorem cert_square_sound (c : ScholzCert.square) :
    cert_square_check K u c = true → (δ_of c)^2 = u
theorem cert_nonsquare_sound (o : Obstruction) :
    cert_nonsquare_check K u o = true → ¬ IsSquare u
```

**S0 文献锚点：** Conrad, *The Local-Global Principle*; Cohen, *CANT* §5.3. Mathlib 可能无现成"数域平方 Hasse"条目 → 自写，综合 local square (`IsSquare` at `𝔭`) + archimedean + class group + units。这是**研究级 Mathlib 推进**，可与 Mathlib 上游合作。

---

## 5. 与 Step 1-7 时间对齐

```text
giac-rs 工程                              Lean 4（giac-proofs/Hasse）
────────────────────────────────────────────────────────────────
Step 1 (A_inf) 全正            ────────►  H1 O2 (可并行启动)
Step 2-3 (A_fin) 赋值 + 嵌入    ────────►  H0 N-val (好素数分支)
Step 4-5 单位系 + torsion       ────────►  H3 B-recon (r=0/r=1 only)
Step 6 (B) unit_is_square       ────────►  H3 B-recon + torsion_sqrt
Step 7a (C) ideal_is_principal  ────────►  H2 C-suff
Step 7b(部分) 虚二次非主认证    ────────►  H2 C-exh (特征值界)
端到端 hasse_sqrt → sqrt_fmodule ──────►  H4 S0 + Cert (证书检查器)
Step 7b(余) 完整 Buchmann       ────────►  远期 (§10): 类群 SNF + 不定域非主
```

**建议：** **端到端接线稳定后再冻结** `ScholzCert` 结构（`hasse_sqrt` 返类型 + 证书字段），避免证了旧出口语义。

---

## 6. 证明战术与 Mathlib 复用

| 目标 | 推荐战术 / 库 |
|------|----------------|
| 范数恒等 / 赋值公式 | `ring`, `field_simp`, `linear_combination` |
| `ℚ` 上具体 minpoly (`t²-2`, `t²+5`) | `native_decide` / `Poly` + `decide` |
| 理想包含 / 主性 | `Ideal.IsPrincipal`, `Ideal.span`, `DedekindDomain` |
| 赋值 | `Ideal.valuation` (若可用) 或自写 `v_𝔭` |
| 正定二次型特征值 | `QuadraticForm`, `Matrix.eigenvalues` (实对称), `PositiveDefinite` |
| 单位 / Dirichlet | `IsUnit`, `unitGroup`, `rank` (Mathlib 单位群部分发展中) |
| 扩域 / 嵌入 | `IntermediateField`, `Algebra.norm`, `Polynomial.Splits` |
| obstruction 逆否 | `mt`, `by_contra` + 上述库 |

**文献锚点（证明草图用）：**
- Keith Conrad, *The Local-Global Principle* — https://kconrad.math.uconn.edu/blurbs/gradnumthy/localglobal.pdf
- Cohen, *A Course in Computational Algebraic Number Theory* §4-6 (赋值 / Minkowski / 类群 / 单位)
- Neukirch, *Algebraic Number Theory* (Dedekind 域 / 类群)
- PARI/GP `bnfinit` / `nfeltissquare` — oracle 对照（非 Lean 依赖）

---

## 7. Rust ↔ Lean 对照表（维护用）

### 7.0 核心对照

| Lean 定理 | Rust 函数 | Rust 出口/测试 |
|-----------|-----------|-----------|
| `Hasse.Obstructions.O1` | `hasse_sqrt` `(A_fin)` `!parity_ok` | `hasse_q_sqrt2_neg1_not_square_via_ainf` 等 |
| `Hasse.Obstructions.O2` | `hasse_sqrt` `(A_inf)` `!is_totally_positive` | `hasse_q_sqrt2_u_2_plus_2sqrt2_not_square` |
| `Hasse.Obstructions.O3` | `hasse_sqrt` `(C)` `!is_principal` | `hasse_q_sqrtneg5_u2_nonprincipal_certified_false` |
| `Hasse.Obstructions.O4` | `hasse_sqrt` `(B)` `!eta_is_sq` | `hasse_q_sqrt3_u_6_minus_3sqrt3_not_square` |
| `Hasse.Principal.C_suff` | `ideal_is_principal` `Some((true,γ))` | `ideal_is_principal_q_sqrt2_j_sqrt2_finds_generator` |
| `Hasse.EigenBound.C_exh` | `imag_quad_search_exhaustive` | `ideal_is_principal_q_sqrtneg5_j_nonprincipal_certified_false` |
| `Hasse.UnitSquare.B_recon` | `unit_sqrt` `Some((true,√η))` | `unit_is_square_q_sqrt2_matches_pari` 等 |
| `Hasse.ScholzIff.S0` | `hasse_sqrt` 全链 | `hasse_*` 端到端 |
| `Hasse.Cert.square_sound` | `δ²=u` 自检 | `sqrt_fmodule_top_recovers_*` (6 个端到端) |
| `Hasse.Norm.norm_of` | `norm_of` (class_group) | `norm_of` 各测试间接 |
| `Hasse.Valuations.val_good_prime` | `compute_ideal_valuations` good 分支 | `scan_*` 测试 |

### 7.1 与 Rust 测试的对应关系（分三层，同 quartic doc §7.1）

```text
         Lean                              Rust
  ─────────────────────            ─────────────────────────
  theorem S (∀ 假设 H, P)    ←→    #[test] + 构造满足 H 的输入
  #check QSqrt2.U2            ←→    ℚ(√2) u=2 同一系数
  Examples.lean 无 sorry      ←→    该测 enabled 且绿
```

#### 层 A — 系数 / 定义对齐（1:1，最易）

| Fixture ID | Lean | Rust 测试 | 断言内容 |
|------------|------|-----------|----------|
| `FIX-QSQRT2` | `Examples.QSqrt2.minpoly` | `hasse_q_sqrt2_*` | `m_α = t²-2` (high-first `[1,0,-2]`) |
| `FIX-QSQRT5NEG` | `Examples.QSqrtNeg5.minpoly` | `hasse_q_sqrtneg5_*` | `m_α = t²+5` (`[1,0,5]`) |
| `FIX-QI` | `Examples.QI.minpoly` | `hasse_q_i_*` | `m_α = t²+1` (`[1,0,1]`) |

**落地：** `giac-proofs/Giac/Fixtures.lean` 与 `class_group::tests` 顶部注释写 `// FIX-QSQRT2`；CI 不要求互跑，人工 diff minpoly 系数。

#### 层 B — 数学性质对齐（定理 ↔ 属性测）

Lean 证 **∀… P**；Rust 用 exact 算术对**具体 u** 检查同一 P。

| Lean 定理 | Rust 测 | 对应方式 |
|-----------|---------|----------|
| `O3` (J 非主 ⟹ 非平方) | `hasse_q_sqrtneg5_u2_nonprincipal_certified_false` | 定理实例化：ℚ(√-5) u=2，J=𝔭 非主 ⟹ `Some(false)` |
| `C_exh` (特征值界) | `ideal_is_principal_q_sqrtneg5_j_nonprincipal_certified_false` | Rust 算 `λ_min=1`，界 `√2+1=2≤256` ⟹ 穷尽 ⟹ 非主 |
| `C_suff` (找到 γ ⟹ 主) | `ideal_is_principal_q_sqrtneg5_j_p2_principal_finds_2` | γ=±2, `|N|=4=N(J)`, `v_𝔭(2)=2=e` ⟹ 主 |
| `B_recon` (`√η²=η`) | `sqrt_fmodule_top_recovers_2_plus_sqrt2_in_q_sqrt2` | δ=2+√2, δ²=6+4√2=u (exact 自检) |
| `S1_true` (`δ²=u`) | `sqrt_fmodule_top_recovers_*` (6 个) | 每个 assert `element_mul(δ,δ)==u` |
| `S1_false` (`¬IsSquare`) | `hasse_q_sqrtneg5_u2_*` | `Some(false)` ⟹ PARI `nfeltissquare=0` |

**关键：** Rust 不证「对所有域成立」，只证 **Lean 定理在 fixture 上的 witness 成立**。定理变 → 测必须仍绿。

#### 层 C — 实现契约（仅 Rust；Lean 给界/证书）

| Rust 测 | Lean 覆盖 | 说明 |
|---------|----------------|------|
| `sqrt_fmodule_non_square_bails_within_4s` (A₄ tower) | 不覆盖 (tower) | Lean `hasse_sqrt` 对 tower 返 `None`（无单生成元）；Rust 测 F-module 行为 |
| `sqrt_fmodule_b6_*_skips_base_case` | `O1` 间接 | (A_fin) 短路 base-case；Lean 验 obstruction，不验 fuel |
| `imag_quad_search_exhaustive` 的 `λ_min` f64 计算 | `C_exh` 数学界 + Rust 算术 | Lean 证界公式；Rust 算 λ_min 数值（f64 不进 Lean，但 `λ_min` 可在 ℚ 上精确算 — 二次型特征值有闭式） |
| `fundamental_units` r≥2 返 `None` | **缺口** | Lean `UnitBasisFull` 前提无法满足；Rust 诚 skip（§10） |

#### 不能 1:1 的部分

| 缺口 | 原因 | 替代 |
|------|------|------|
| `recover_free_exponents` f64 log-map | 数值启发式 | 当 oracle；`ζ^|μ(K)|=1` exact 检查当证书（cert 里带 `a_vec` + `ζ`） |
| `ideal_is_principal` 盒内枚举 | bounded search | 当 oracle；找到 γ 后 `|N(γ)|=N(J)∧γ∈J` 当证书 |
| `NORM_PREFILTER_TOL` / `IDEAL_GEN_COORD_BOUND=256` 等常数 | 工程界 | Lean 用 `∃ B, √(N/λ_min) ≤ B` 存在量词；Rust 用 256 实例化 |
| f64 嵌入 `embeddings` | 数值 | `eval_at_embedding` 不进 Lean；用 `Algebra.norm` + 精确根符号代替 |

### 7.2 建议：Hasse 测试注册表（同 quartic §7.2）

```yaml
# hasse-test-registry.yaml（示意）
- id: FIX-QSQRT2
  field: "Q[sqrt(2)]"
  minpoly_high: [1, 0, -2]
  lean: Giac.Examples.QSqrt2
  rust:
    - hasse_q_sqrt2_u2_is_square                  # TRUE ② δ=α
    - hasse_q_sqrt2_u_6_plus_4sqrt2_is_square     # TRUE ② δ=2+√2 (B-recon)
    - hasse_q_sqrt2_u_2_plus_2sqrt2_not_square    # FALSE ② (A_inf)
    - hasse_q_sqrt2_neg1_not_square_via_ainf      # FALSE ②
    - sqrt_fmodule_top_recovers_sqrt2_in_q_sqrt2
    - sqrt_fmodule_top_recovers_2_plus_sqrt2_in_q_sqrt2
  lean_theorems:
    - Giac.Hasse.Obstructions.O2
    - Giac.Hasse.UnitSquare.B_recon
    - Giac.Hasse.Cert.square_sound
  status: partial  # 待 Lean 侧 H1/H3/H4 落地

- id: FIX-QSQRT5NEG
  field: "Q[sqrt(-5)]"
  minpoly_high: [1, 0, 5]
  lean: Giac.Examples.QSqrtNeg5
  rust:
    - hasse_q_sqrtneg5_u2_nonprincipal_certified_false   # FALSE ③ (C-exh)
    - hasse_q_sqrtneg5_u4_principal_certified_true       # TRUE ② γ=±2
    - ideal_is_principal_q_sqrtneg5_j_nonprincipal_certified_false
    - ideal_is_principal_q_sqrtneg5_j_p2_principal_finds_2
    - sqrt_fmodule_top_recovers_2_in_q_sqrtneg5
    - sqrt_fmodule_top_q_sqrtneg5_u2_non_square_none
  lean_theorems:
    - Giac.Hasse.Obstructions.O3
    - Giac.Hasse.EigenBound.C_exh
    - Giac.Hasse.Principal.C_suff
    - Giac.Hasse.Cert.nonsquare_sound
  status: partial

- id: FIX-QSQRT3
  field: "Q[sqrt(3)]"
  minpoly_high: [1, 0, -3]
  rust:
    - hasse_q_sqrt3_u3_is_square                     # TRUE ②
    - hasse_q_sqrt3_u_6_minus_3sqrt3_not_square      # FALSE ④ (B: η=ε 非单位平方)
    - sqrt_fmodule_top_recovers_sqrt3_in_q_sqrt3
  lean_theorems:
    - Giac.Hasse.Obstructions.O4
    - Giac.Hasse.UnitSquare.B_recon
  status: partial

- id: FIX-QI
  field: "Q[i]"
  minpoly_high: [1, 0, 1]
  rust:
    - hasse_q_i_neg1_is_square                       # TRUE ② (torsion 平方)
    - hasse_q_i_i_not_square                         # FALSE ④ (i 非 μ₄ 平方)
    - sqrt_fmodule_top_recovers_i_in_q_i
  lean_theorems:
    - Giac.Hasse.UnitSquare.torsion_sqrt
    - Giac.Hasse.Obstructions.O4
  status: partial

- id: FIX-QSQRT6NEG
  field: "Q[sqrt(-6)]"
  minpoly_high: [1, 0, 6]
  rust:
    - hasse_q_sqrtneg6_u2_nonprincipal_certified_false
    - hasse_q_sqrtneg6_u4_principal_certified_true
  lean_theorems:
    - Giac.Hasse.EigenBound.C_exh
  status: partial
```

**工作流：** 同 quartic — 改 `hasse_sqrt` 出口语义 → 同步改 Lean `ScholzCert` + registry；PR 模板勾选。

### 7.3 Lean checker ↔ Rust `δ²=u` 自检的强连接

**这是与代码最直接的连接点：** 每个端到端测试里的 `assert!(element_mul(&delta,&delta) == u)` 就是 `cert_square_sound` 的**运行时实例**。Lean 证 checker 对**任意**证书 sound；Rust 测证 checker 在**具体** fixture 上接受正确的 δ。

```rust
// Rust 测试 (poly_roots.rs::sqrt_fmodule_top_recovers_*) 已有的自检:
let sq = field.element_mul(&delta, &delta).expect("δ² computable");
let ok = sq.iter().zip(u.iter()).all(|(a, b)| a == b);
assert!(ok, "{label}: δ² must equal u (end-to-end soundness)");
```

```lean
-- Lean (Giac.Hasse.Cert):
theorem cert_square_sound (c : ScholzCert.square) :
    cert_square_check K u c = true → (δ_of c)^2 = u
-- Rust 的 `ok` = `cert_square_check` 在具体 c 上的运行结果
```

未来可加一个 `#[test] fn cert_checker_matches_rust_selfcheck` 对同一 `(K,u,δ)` 同时跑 Rust 自检和（若 checker 提取到 Rust）证书检查器，断言两者一致。

---

## 8. CI、版本与维护

同 quartic doc §8：Pin `lean-toolchain` + Mathlib commit；`sorry` 政策（H0 + 至少一个 `Examples` 无 sorry 进 main）；不在 Rust pre-commit 跑 Lean；Hasse 定理变更 = 破坏性变更，同步改 hasse-proof.md 与本表。

---

## 9. 人力与里程碑（粗估）

| 里程碑 | 内容 | 人周 |
|--------|------|------|
| M0 | `giac-proofs/Hasse/` 骨架 + H0 范数/赋值 | 2-3 |
| M1 | H1 O1-O4 (4 个 obstruction，核心数学) | 3-4 |
| M2 | H2 C-suff + C-exh (主性 + 特征值界) | 4-6 |
| M3 | H3 B-recon + torsion_sqrt (r=0/r=1 域) | 3-4 |
| M4 | H4 S0 (Hasse iff) + Cert 检查器 + Examples 无 sorry | 4-6 |

**总计：** 约 **4-6 人月**（熟悉 Mathlib 代数数论方向）；可与 quartic 项目 **并行**（共享 `giac-proofs` 基础设施）；M1 (obstruction) 可最先落地，因短而独立。

---

## 10. 若将来要加强

| 路径 | 说明 | 难度 |
|------|------|------|
| **A. r≥2 全实域 index-1 单位证书** | Lean 证 Friedman 型 regulator 下界 ⟹ 单位基 full ⟹ 解锁 r≥2 的 (B) 证书 | 研究级（Mathlib 缺口，可上游合作） |
| **B. 完整 Buchmann 类群 (Step 7b 余)** | Lean 形式化 Minkowski 界 + 关系格 + SNF ⟹ 类群结构 ⟹ 不定域/高次非主 `Some(false)` 证书 | 研究级（类群算法形式化几乎无先例） |
| **C. (A_fin) 坏素数 / Montes** | Lean 形式化 Montes 带分 ramification 理论 ⟹ 解锁 `p∣index` 分支 | 研究级（Montes 理论形式化空白） |
| **D. 证书 checker 提取到 Rust** | Lean → Rust 提取 `cert_square_check`/`cert_nonsquare_check`，Rust 运行时验证书 | 中（需提取工具成熟） |
| **E. 深嵌入 hasse_sqrt** | Aeneas / Lean4Lean 提取整个 `hasse_sqrt` ⟹ Rust refines Lean spec | 高（CAS ℚ(α) 表示 vs Mathlib 差太大，同 quartic §10） |

**r=0/r=1 域（当前已实现 + 测试覆盖）是 Lean 验证的天然第一版范围**；r≥2 / 不定域非主 / 坏素数是远期，依赖 Mathlib 代数数论推进。

---

## 11. 索引

| 主题 | 文档 |
|------|------|
| **Scholz 管线实施（Step 1-7 + 端到端）** | [giac-poly-f5-fglm-hasse-proof.md](../giac-poly-f5-fglm-hasse-proof.md) |
| (A) 赋值公式证明 | [giac-poly-f5-fglm-ideal-valuation-proof.md](../giac-poly-f5-fglm-ideal-valuation-proof.md) |
| 总 issue / Phase C 优先级 | [GIAC-poly-f5-fglm-complete-square-decision.md](GIAC-poly-f5-fglm-complete-square-decision.md) |
| **同套 quartic Lean 4 验证方案** | [GIAC-poly-quartic-lean4-verification.md](GIAC-poly-quartic-lean4-verification.md) |
| 塔 L0-L7 分层 | [giac-tower-common-math.md](../giac-tower-common-math.md) §3 |
| 测试即规格 | [conformance-testing.md](../conformance-testing.md) |
| 符号偏离 | [known-divergences.md](../known-divergences.md) |
