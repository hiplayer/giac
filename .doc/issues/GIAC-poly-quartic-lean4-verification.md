# GIAC-poly — 四次开方算法的 Lean 4 验证方案（未来）

**状态:** draft / 长期  
**类型:** 形式化 / 架构  
**前置:** [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) **F3 绿**（算法规格冻结后再证「实现≈规格」）  
**相关:** [giac-tower-common-math.md](../giac-tower-common-math.md) §3、[GIAC-algext-adoption](GIAC-algext-adoption.md) §8  
**工具:** Lean 4 + Mathlib4（主）；Rust CI **不阻塞** on `lake build`（独立 job）

---

## 1. 目标与非目标

### 1.1 要证什么

在 **抽象域论 / 多项式** 层证明下列命题（名称与 Rust 模块对齐，便于交叉引用）：

| 层级 | 命题（示意） | 对应 Rust |
|------|----------------|-----------|
| **Q0** | 抑郁化 \(t \mapsto t-a/4\) 保持根集（仿射） | `depress_quartic` |
| **Q1** | Ferrari resolvent \(R\) 的三根 \(\alpha,\beta,\gamma\) 与 depressed \(x^4+px^2+qx+r\) 的系数关系 | `build_resolvent_cubic` |
| **Q2** | Euler 公式：在含 \(\sqrt\alpha,\sqrt\beta,\sqrt\gamma\) 的扩域中，四元组 \(\frac12(\varepsilon_1\sqrt\alpha+\varepsilon_2\sqrt\beta+\varepsilon_3\sqrt\gamma)\)（\(\varepsilon_i\in\{\pm1\}\)，\(\varepsilon_1\varepsilon_2\varepsilon_3=1\)）为零化子 | `euler_four_roots_from_triple` |
| **Q3** | 当 \(q\neq0\) 时 \(\sqrt\gamma = -q/(\sqrt\alpha\sqrt\beta)\)（在约定平方根分支下） | `euler_derived_sqrt_gamma` |
| **Q4** | 四根 **集合** 等于原四次在分裂域中的根集（在可解 / 根式可表情形） | `quartic_roots` e2e |
| **T1** | 若 \(u=\varepsilon^2\) 且 \(\varepsilon\in L\)，则 blind adjoin \(u^2-\alpha\) 次数严格冗余（存在同态意义下更小扩域） | F1 `sqrt_in_field` |
| **T2** | Resolvent 三根可落在 \(L_3\) 的数学陈述（\(d_L\le 6\) resolvent；全流程 \(\le 24\) 硬顶；**12 仅 Gal≅A₄ 已证时**） | F2 / F4 |

### 1.2 不证什么（第一版）

| 不做 | 原因 |
|------|------|
| Rust `ext_tower` / `FieldSession` **逐行提取**（Aeneas / Lean export） | 成本极高；CAS 多项式表示与 Mathlib `Polynomial` 差太大 |
| `coords` 张量基与 giac `poly1` 顺序 **字节级一致** | 属表示层；用 `eq_mod` / 同构陈述代替 |
| giac C++ golden **字面一致** | 见 [known-divergences.md](../known-divergences.md) DIV-080 |
| 通用「算法终止 / 复杂度」 | 用 Rust release 测 + `dim ≤ B` 工程界 |
| 五次及以上 Abel–Ruffini | 范围外 |

**原则：** Mathlib 证 **数学规格**；Rust 用测试证 **实现满足规格**；文档用 `Theorem giac.Q2` 链两者。

---

## 2. 架构：三层证书

```text
┌─────────────────────────────────────────────────────────┐
│  L-数学层 (Lean 4 / Mathlib)                             │
│  定理 Q0–Q4, T1–T2：存在性 + 公式正确 + 次数界          │
└───────────────────────────┬─────────────────────────────┘
                            │ 文档引用定理名 + 假设列表
┌───────────────────────────▼─────────────────────────────┐
│  S-规格层 (Lean 结构体 / Rust doc 镜像)                  │
│  QuarticRootsSpec：输入 (p,q,r) → 输出四元组 + 不变量    │
└───────────────────────────┬─────────────────────────────┘
                            │ root_vanishes / eq_mod 测试
┌───────────────────────────▼─────────────────────────────┐
│  I-实现层 (giac-core poly_roots / field_session)         │
│  F1–F3 代码；不 machine-check 到 Lean                    │
└─────────────────────────────────────────────────────────┘
```

**可审计闭环：** 每条 Lean 定理在 `.doc` 与 `poly_roots.rs` 模块注释写 `/** Lean: Giac.Quartic.Q2 */`；Rust 测 `verify_root` 覆盖定理假设的 **实例化**（如 \(t^4+t+1\)）。

---

## 3. 仓库与目录（建议）

独立 Lake 项目，不污染 `giac-rs` 编译：

```text
giac-proofs/                    # 新 git 子目录或 sibling repo
  lakefile.lean
  lean-toolchain                # 与 Mathlib 兼容 pin（如 v4.16.0）
  Giac/
    Root/
      Basic.lean                # 抑郁化、resolvent 多项式定义
      ResolventCubic.lean       # Q1
      Euler.lean                # Q2, Q3
      QuarticRoots.lean         # Q4 汇总
    Tower/
      Adjoin.lean               # 对齐 ext_tower L2
      SqrtInField.lean          # T1
      DegreeBound.lean          # T2, F4
    Examples/
      T4PlusTPlus1.lean         # t^4+t+1 实例（#eval 可选）
  README.md                     # lake build 说明；定理索引
```

**依赖：** `require mathlib from git`；不 require Rust。

**CI：** 可选 workflow `lean-proofs.yml`：`lake exe cache get && lake build`；失败不挡 `cargo test`。

---

## 4. 数学形式化路线（分阶段）

### Phase P0 — 多项式与抑郁化（2–3 周）

**Mathlib 已有：** `Polynomial`, `Polynomial.eval`, `Polynomial.roots`, `Field`, `Splits`.

**新建定义（`Giac.Root.Basic`）：**

```lean
-- 示意；正式版用 Mathlib 命名规范
noncomputable def depress (P : ℚ[X]) (a : ℚ) : ℚ[X] := ...
lemma depress_roots (P a) :
  (depress P a).roots = (P.map (C · - a)).roots  -- 或正确的仿射换元陈述
```

**验收：** `depress (X^4 + X + 1) 0` 系数为 `(0,0,1,1)` 用 `native_decide` 或 `ring` 战术。

**与 F1–F5：** 与 F3 前可并行；不依赖塔实现。

---

### Phase P1 — Resolvent 三次（Q1）（3–4 周）

对 monic depressed \(f = X^4 + pX^2 + qX + r\)：

\[
R(X) = X^3 - pX^2 - 4rX + (4pr - q^2)
\]

**要证的引理：**

1. `resolvent_poly p q r` 定义与 Rust `build_resolvent_cubic` 系数一致（`ℚ` 上 `ext` 或 `decide`）。
2. **经典事实：** 若 \(\alpha\) 为 \(R\) 的根，则 \(X^4+pX^2+qX+r\) 在 \(\mathbb{Q}(\sqrt\alpha)\) 上可分解为两个二次（或给出分裂域次数整除 24）。  
   - Mathlib 可能无现成 Ferrari 条目 → 自写引理，参考 *Dummit & Foote* §13.4 / Lang §IV.6。
3. **实例：** \(p=0,q=1,r=1\) 时 \(R = X^3 - 4X - 1\)（与 `resolvent_golden_t4_plus_t_plus_1` 一致）。

**难点：** 在 **一般系数** \(p,q,r\) 下证不可约性 / 分裂域结构需假设；建议定理带 `Hypothesis`：

```lean
variable (p q r : ℚ)
variable (hR : Irreducible (resolvent_poly p q r))  -- 或 weaker: separable

theorem resolvent_root_parametrizes_quartic_splitting ...
```

**Rust 链接：** `resolvent_golden_t4_plus_t_plus_1` → Lean `Examples.T4PlusTPlus1.resolvent_eq`.

---

### Phase P2 — Euler 四根公式（Q2, Q3）（4–6 周）

在域 \(F\) 上，设 \(\alpha,\beta,\gamma\) 为 \(R\) 的三根；设 \(s_\alpha,s_\beta,s_\gamma\in F\) 满足 \(s_\alpha^2=\alpha\) 等。

**定义（与 Rust 一致）：**

```lean
def eulerRoot (sα sβ sγ : F) (ε₁ ε₂ ε₃ : ℤ) : F :=
  (ε₁ • sα + ε₂ • sβ + ε₃ • sγ) / 2

def eulerSignTable : List (ℤ × ℤ × ℤ) :=
  [(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)]
```

**核心定理 `euler_root_vanishes`：**

```lean
theorem euler_root_vanishes (hγ : sγ^2 = γ) (hε : ε₁ * ε₂ * ε₃ = 1) ... :
  eval (depressed_poly p q r) (eulerRoot sα sβ sγ ε₁ ε₂ ε₃) = 0
```

**定理 `euler_gamma_from_q`（\(q\neq0\)）：**

```lean
theorem sqrt_gamma_relation (hq : q ≠ 0) (hα hβ : sα^2 = α ∧ sβ^2 = β) :
  sγ^2 = γ → sα * sβ * sγ + q = 0  -- 在约定分支下；或陈述为 γ = q²/(αβ)
```

**分支约定：** Lean 中 **不** 形式化 `approx_real_sign`；用存在量词「∃ 一组平方根满足关系」+ Rust F3 **确定性** 分支作为 witness 构造的实例测。

**与 DIV-077：** 定理对 **任意** 满足 \(s^2=u\) 的平方根成立；canonical 分支是 witness 选取，不是数学必要性。

---

### Phase P3 — 塔与 `sqrt_in_field`（T1, T2）（与 F1/F2 同步，4–8 周）

复用 [giac-tower-common-math.md](../giac-tower-common-math.md) L2–L5 路线。

**T1（F1 数学内核）：**

```lean
theorem sqrt_exists_in_extension (L : IntermediateField ℚ K) (u : L)
    (hu : IsSquare u) :
  ∃ ε : L, ε^2 = u

theorem blind_adjoin_redundant ... :
  -- 若 u = ε², ε ∈ L，则 adjoin X^2 - u 给出 [L(√u):L]=1（在合适不可约假设下）
```

**T2（F2/F4 次数界）：**

对固定规格（Cardano 一条 casus 根 + 同场另两根，不二次 adjoin）：

```lean
theorem quartic_resolvent_tower_finrank_le :
  finrank ℚ L₃ ≤ 24   -- D_HARD；finrank ≤ 12 仅 Gal≅A₄ 已证时的推论
```

**实现：** 用 `IntermediateField.adjoin` 塔、`FiniteDimensional.finrank_mul` 组合；**不** 建模 Rust `Arc<ExtensionField>`。

**Rust：** F4 `field_session_dimension_bound_quartic` 是 T2 的 **工程实例**；Lean 证上界，Rust 测 tight 例。

---

### Phase P4 — 汇总定理 Q4（2 周）

```lean
theorem quartic_roots_exist (f : ℚ[X]) (hf : f.natDegree = 4) (hsep : Separable f) :
  ∃ (L : Type*) [Field L] [FiniteDimensional ℚ L],
    ∀ x, f.eval x = 0 → x ∈ algebraMap ℚ L '' (univ : Set ℚ) ∨
    x ∈ rootsInL f L   -- 四根在 L 中
```

**弱版（第一版足够）：** 在 **可解** 四次（resolvent 分裂）假设下，Euler 四元组 **覆盖** \(f\) 在分裂域中的全部根（集合相等）。

**强版（后期）：** 与 Rust `dedup_roots` 对应的 **重数 / 共轭配对**。

---

## 5. 与 F1→F5 的时间对齐

```text
giac-rs 工程                         Lean 4（giac-proofs）
────────────────────────────────────────────────────────────
F1 sqrt_in_field          ────────►  P3 T1（可并行启动）
F2 resolvent 三根 L₃      ────────►  P3 T2 + P1 实例化
F3 Euler 规格冻结         ────────►  P2 Q2,Q3 冻结签名；写 Examples
F4 维数 bound             ────────►  P3 T2 对齐常数 B
F5 ext_tower API          ────────►  可选：refactor 不改定理陈述
```

**建议：** **F3 合并后再冻结** Lean 中 `eulerSignTable` 与 `sqrt_gamma` 公式，避免证了旧穷举路径。

---

## 6. 证明战术与 Mathlib 复用

| 目标 | 推荐战术 / 库 |
|------|----------------|
| 系数恒等 | `ring`, `field_simp`, `linear_combination` |
| \( \mathbb{Q} \) 上具体多项式 | `native_decide`（小度）、`MvPolynomial` + `decide` |
| 扩域次数 | `IntermediateField.adjoin`, `finrank_mul`, `FiniteDimensional` |
| 平方根存在 | `IsSquare`, `exists_sq_eq` |
| 分裂域 | `Polynomial.Splits`, `splittingField`（若 Mathlib 版本足够新） |
| 不可约实例 | `irreducible_X_pow_sub_C` 等 + 组合；\(t^4+t+1\) 可单独 `sorry`→逐步填 |

**文献锚点（证明草图用）：**

- Dummit & Foote, *Abstract Algebra*, §13.4（四次公式）
- Garling, *Galois Theory*（resolvent cubic）
- Cohen, *Computational Algebraic Number Theory* §3（工程 adjoin，非 Lean 必需）

---

## 7. Rust ↔ Lean 对照表（维护用）

| Lean 定理 | Rust 函数 | Rust 测试 |
|-----------|-----------|-----------|
| `Giac.Root.Basic.depress_coeffs` | `depress_quartic` | `depressed_t4_plus_t_plus_1_coeffs` |
| `Giac.Root.ResolventCubic.coeffs` | `build_resolvent_cubic` | `resolvent_golden_t4_plus_t_plus_1` |
| `Giac.Root.Euler.vanishes` | `euler_four_roots_from_triple` | `euler_four_roots_vanish` |
| `Giac.Root.Euler.gamma_relation` | `euler_derived_sqrt_gamma` | `euler_gamma_relation_holds`（待加） |
| `Giac.Tower.SqrtInField.exists` | `try_sqrt_in_field` | `adjoin_sqrt_squares_*` |
| `Giac.Tower.DegreeBound.le24` | `FieldSession::working` | `field_session_dimension_bound_quartic` |
| `Giac.Tower.DegreeBound.le12_a4` | （可选）Gal≅A₄ 分测 | 须附 Gal 证明，非 `t⁴+t+1` 默认 |
| `Giac.Root.QuarticRoots.complete` | `quartic_roots` | `roots_quartic_t4_plus_t_plus_1` |

### 7.1 与 Rust 测试的对应关系（能，分三层）

**结论：** Lean 证 **全称命题**；Rust 测 **同一命题在具名实例上的见证（witness）**。二者通过 **共享 fixture ID** 对齐，不是同一个可执行文件跑两遍。

```text
         Lean                              Rust
  ─────────────────────            ─────────────────────────
  theorem T (∀ 假设 H, P)    ←→    #[test] + 构造满足 H 的输入
  #check T4PlusTPlus1          ←→    t^4+t+1 同一系数向量
  Examples.lean 无 sorry       ←→    该测 enabled 且绿
```

#### 层 A — 系数 / 定义对齐（1:1，最易）

Lean 与 Rust 断言 **同一 ℚ 多项式系数**；无需扩域实现。

| Fixture ID | Lean | Rust 测试 | 断言内容 |
|------------|------|-----------|----------|
| `FIX-T4P1` | `Examples.T4PlusTPlus1.depressed` | `depressed_t4_plus_t_plus_1_coeffs` | \((p,q,r)=(0,1,1)\) |
| `FIX-T4P1` | `Examples.T4PlusTPlus1.resolvent` | `resolvent_golden_t4_plus_t_plus_1` | \(R=z^3-4z-1\) |
| `FIX-Z3M4Z1` | （可选）`Examples.ResolventZ3` | `resolvent_one_cubic_root_z3_minus_4z_minus_1` | 一根零化 resolvent |

**落地：** 在 `giac-proofs/Giac/Fixtures.lean` 与 `poly_roots::tests` 顶部注释写 `// FIX-T4P1`；CI 不要求互跑，人工 diff 系数即可。

#### 层 B — 数学性质对齐（定理 ↔ 属性测）

Lean 证 **∀… P(x)=0**；Rust 用 `verify_root` / `eq_mod` 对 **具体根** 检查同一 P。

| Lean 定理 | Rust 测 | 对应方式 |
|-----------|---------|----------|
| `Euler.vanishes` | `euler_four_roots_vanish` | 定理实例化：F3 固定分支下四根 `verify_root(dep)` |
| `Euler.gamma_relation` | `euler_gamma_relation_holds`（F3 待加） | Rust 构造 \((s_\alpha,s_\beta,s_\gamma)\) 验 \(s_\alpha s_\beta s_\gamma + q = 0\) |
| `ResolventCubic.root_vanishes` | `resolvent_cubic_all_roots_vanish` | 三根零化 R（F2 新路径） |
| `SqrtInField.squares` | `adjoin_sqrt_squares_one_resolvent_root` | Rust 验 \(\varepsilon^2 \equiv u\)（`eq_mod`）≈ Lean 的 \(ε^2=u\) |
| `QuarticRoots.complete` | `roots_quartic_t4_plus_t_plus_1` | 四根零化原式 + `len==4` |

**关键：** Rust 不证明「对所有域成立」，只证明 **Lean 定理在 FIX-T4P1 上的 witness 成立**。定理变 → 测必须仍绿。

#### 层 C — 实现契约（仅 Rust；Lean 给上界）

| Rust 测 | Lean 是否覆盖 | 说明 |
|---------|----------------|------|
| `field_session_dimension_bound_quartic` | T2 给 **上界** \(B\) | Lean 证 \(\mathrm{finrank}\le B\)；Rust 测实际 `dimension()` |
| `adjoin_sqrt_after_quadratic_roots_formula`（F1 待加） | T1 存在性 | Lean 不跑 adjoin 代码；Rust 验 ε²=u |
| `roots_biquadratic_t4_minus_2` | 非四次主链 | 可另开 `Biquadratic` 引理或标「仅 Rust」 |
| `approx_real_sign` / 6×flip 穷举 | **不对应** | F3 删除后 Lean 不形式化 |

#### 不能 1:1 的部分

| 缺口 | 原因 | 替代 |
|------|------|------|
| `coords` / `rootof` 打印 | Lean 用抽象域元素 | `eq_mod` + DIV-083 |
| giac golden 四根字面 | DIV-080 | 集合零化 |
| 全 `cargo test` 矩阵 | 测太多、非四次 | 只登记 **QuarticRegistry** 子集 |

### 7.2 建议：Quartic 测试注册表（单一来源）

在 `.doc/issues/GIAC-poly-quartic-roots-F1-F5.md` 或 `giac-proofs/README.md` 维护：

```yaml
# quartic-test-registry.yaml（示意）
- id: FIX-T4P1
  poly: "t^4+t+1"
  lean: Giac.Examples.T4PlusTPlus1
  rust:
    - depressed_t4_plus_t_plus_1_coeffs
    - resolvent_golden_t4_plus_t_plus_1
    - euler_four_roots_vanish      # 需 F3
    - roots_quartic_t4_plus_t_plus_1
  lean_theorems:
    - Giac.Root.Basic.depress_coeffs
    - Giac.Root.ResolventCubic.coeffs
    - Giac.Root.Euler.vanishes
    - Giac.Root.QuarticRoots.complete
  status: partial  # euler/quartic e2e ignore 中

- id: FIX-SQRT-RESOLVENT
  lean: Giac.Tower.SqrtInField
  rust:
    - adjoin_sqrt_squares_one_resolvent_root
    - adjoin_sqrt_after_resolvent_deflate_only
  status: green
```

**工作流：**

1. 改 `build_resolvent_cubic` → 同步改 Lean `resolvent_poly` + `FIX-T4P1` 两侧。
2. PR 模板勾选：「已更新 registry / Lean Examples」。
3. 可选 CI：`scripts/check-quartic-registry.sh` 解析 Rust `// FIX-*` 注释与 yaml 一致。

### 7.3 CI 矩阵（解耦但可对账）

```text
job rust-poly-roots:     cargo test -p giac-core poly_roots::tests --release
job lean-giac-proofs:    cd giac-proofs && lake build
job registry-audit:      可选；yaml ↔ 注释一致性
```

两 job **独立**；对账靠 registry + 同一 fixture ID，不靠「Lean 调用 Rust」。

---

## 8. CI、版本与维护

1. **Pin** `lean-toolchain` + Mathlib commit；每月可 Dependabot 式升级（单独 PR）。
2. **`sorry` 政策：** 允许在 P1–P2 用 `sorry` 搭骨架；进 `main` 前 P0 + 至少一个 `Examples` 无 `sorry`。
3. **不** 在 Rust pre-commit 跑 Lean。
4. 定理变更 = **破坏性变更**：同步改 F1–F5 文档与本表。

---

## 9. 人力与里程碑（粗估）

| 里程碑 | 内容 | 人周 |
|--------|------|------|
| M0 | `giac-proofs` 骨架 + P0 `depress` | 1–2 |
| M1 | P1 resolvent 定义 + \(t^4+t+1\) 实例无 `sorry` | 2–3 |
| M2 | P2 Euler `vanishes` 对一般 \(p,q,r\)（带假设） | 4–6 |
| M3 | P3 T1 + T2 与 F1/F4 常数对齐 | 4–8 |
| M4 | Q4 汇总 + README 定理索引 | 2 |

**总计：** 约 **3–6 人月**（熟悉 Mathlib 的域论方向）；可与 giac-rs F1–F3 **并行** 由不同人维护。

---

## 10. 若将来要加强到「实现证书」

仅作远期选项，**不** 纳入第一版：

| 路径 | 说明 |
|------|------|
| **A. 精简模型** | Lean 中定义 `CoordsQ` + `element_mul` 的 **规范函数**，Rust 测与 Lean `#eval` 小域一致 |
| **B. 反例搜索** | Lean `slim_check` / 外部 QuickCheck 生成 \((p,q,r)\) 反例 |
| **C. 提取** | Aeneas / Lean4Lean 提取 `element_mul` — 仅当塔 API 冻结数年 |

行业参照：HACL*、CryptoMiniSat 证核心引理；**开放 CAS 几乎不 full proof kernel**。

---

## 11. 索引

| 主题 | 文档 |
|------|------|
| **giac-proofs 说明（depression / Euler 定理 / lake build）** | [giac-proofs/.doc/README.md](../../giac-proofs/.doc/README.md) |
| 工程实施 F1–F5 | [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) |
| 塔 L0–L7 分层 | [giac-tower-common-math.md](../giac-tower-common-math.md) §3 |
| 符号偏离 | [known-divergences.md](../known-divergences.md) DIV-076–083 |
| 测试即规格 | [conformance-testing.md](../conformance-testing.md) |
