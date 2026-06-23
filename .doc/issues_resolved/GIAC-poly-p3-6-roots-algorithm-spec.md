# P3-6 — `poly_algext_roots` 算法规格（normative）

**状态:** **resolved**（2026-06-23 验收；deg 1–4 normative 管线已落地）  
**续篇（plan）:** [GIAC-poly-p3-6-roots-algorithm-spec.md](../issues/GIAC-poly-p3-6-roots-algorithm-spec.md) §12（F5 / F4′）  
**父项:** [GIAC-poly-p3-6-quartic-roots-gaps.md](../issues/GIAC-poly-p3-6-quartic-roots-gaps.md)  
**实现:** [GIAC-poly-quartic-roots-F1-F5](../issues/GIAC-poly-quartic-roots-F1-F5.md)、[GIAC-poly-roots-field-session-plan](../issues/GIAC-poly-roots-field-session-plan.md)  
**Rust 落点:** `giac-core::algebra::{field_session, poly_roots}`  
**约定:** [algorithm-expr-api.md §6.3](../algorithm-expr-api.md)、[giac-poly-nested-ring.mdc](../../.cursor/rules/giac-poly-nested-ring.mdc)  
**验收:** `cargo test -p giac-core poly_roots::tests --release` — **25 passed, 0 ignored**（~11s）

---

## 0. 设计原则（收紧）

1. **一条三次分裂 API** — 独立三次、resolvent 三次、（间接）四次内层 **共用** `split_monic_cubic_roots_in_session`；**禁止**求三根时走 Cardano `u+v` 叠层（会污染 session，见 C1 `dim=18`）。
2. **F2 A′ 分工** — **ω 共轭仅用于纯三次** \(t^3+a_0\)（及 depressed \(p_{\mathrm dep}=0\)）；**\(p_{\mathrm dep}\neq 0\)**（含 casus resolvent `z³−4z−1`、`x³−x+1`）→ **不可约 adjoin 一根 + deflate 二次 + F1 √Δ**。**禁止**对 \(p_{\mathrm dep}\neq 0\) 用 \(\omega^k z_0\)（数学上不成立，见 §3.2 注）。
3. **开方唯一入口** — 生产路径仅 `FieldSession::sqrt_in_field` / `sqrt_disc`（后者内部仅调 `sqrt_in_field`）；**禁止** `algext_square_roots` 出现在 `poly_roots` 热路径；`adjoin_sqrt_new` **仅**在 `sqrt_in_field` blind fallback 内。
4. **维数硬顶** — 常量 `POLY_ROOTS_DIM_HARD=24`、`POLY_ROOTS_DIM_RESOLVENT=6`（`field_session.rs`）；超过硬顶立即 `Err`。**禁止**在无 Galois 群 / 分裂域次数证据时，把中间值（如 12）写成全体四次 normative 目标（见 §5.4）。
5. **根输出契约** — 返回根列表经 `dedup_roots` + `canonical_sort_roots`；单测/调试对每个根 `verify_root(p, var, r)`。

---

## 1. 全局参数（收紧）

符号：\(K = \texttt{session.ambient()}\)，\(L = \texttt{session.working()}\)，\(d_K = \dim(K/\mathbb{Q})\)，\(d_L = \dim(L/\mathbb{Q})\)。

| 参数 | Rust 常量 | 值 | 含义 |
|------|-----------|-----|------|
| **`D_PURE_CUBIC`** | — | **3** | 纯三次 \(t^3+a_0\)：\(\mathbb{Q}(\sqrt[3]{|a_0|},\omega)\) |
| **`D_CUBIC_SPLIT`** | — | **6** | 一般三次 / S3 分裂（\(d_K=1\)）：不可约 \(z_0\) + deflate 二次一次 adjoin |
| **`D_RESOLVENT`** | `POLY_ROOTS_DIM_RESOLVENT` | **6** | resolvent 三次结束（\(d_K=1\)）；**normative** |
| **`D_QUARTIC_A4_REF`** | `POLY_ROOTS_DIM_QUARTIC_OUT` | **12** | **参考值**：仅当 \(\mathrm{Gal}(f)\cong A_4\)（故 \([K_f:\mathbb{Q}]=12\)）**已证**时的分裂域次数；**非**全体四次出口断言 |
| **`D_HARD`** | `POLY_ROOTS_DIM_HARD` | **24** | 全管线硬失败（\(\mathrm{Gal}(f)\le S_4\) 时 \([K_f:\mathbb{Q}]\mid 24\)）；**normative** |
| **`T_ROOTS`** | — | **10s** | `cargo test -p giac-core poly_roots::tests --release` |
| **`TRY_SQRT_PAIRWISE_DIM`** | `TRY_SQRT_PAIRWISE_DIM_CEILING` | **12** | F1 1c 二元 ±1 枚举（ponytail interim） |
| **`TRY_SQRT_TRIPLE_DIM`** | `TRY_SQRT_TRIPLE_DIM_CEILING` | **9** | F1 1d 三元 ±1 枚举（ponytail interim） |
| **`TRY_SQRT_QUADRATIC_COMBO_DIM`** | `TRY_SQRT_QUADRATIC_COMBO_DIM_CEILING` | **6** | F4 eᵢ±k·g 小组合（ponytail interim） |

**interim 说明：** 上表三行是 **性能旋钮**，漏检 → `adjoin`；**非终态**。终态见 §12（F5 + Galois 替代枚举）。

**相对上界（\(d_K>1\)）：** 保守写 `d_L <= d_K * D_CUBIC_SPLIT`（或阶段表值 × \(d_K\)）。

**已废弃参数（初版规格）：** `D_CUBIC_CASUS` 作「\(\mathbb{Q}(z_0,\omega)\)」仅适用于 **纯三次**；`MAX_ADJOIN_CARDANO` 不再用于 **三根分裂** 路径（单根 API `one_cubic_root` 仍可用 Cardano）。

---

## 2. 入口管线（不变）

```text
poly_algext_roots_for_ctx(p, var, ctx)
  K  = infer_field(p)
  session = ctx.session().fork_ambient(K)
  p' = monic_univariate(normalize_coeffs(p, session), var, session)
  roots_dispatch(session, p', var)

roots_dispatch(session, p, var):
  d = deg_wrt(p, var)
  d=1 → linear_root
  d=2 → quadratic_roots_formula
  d=3 → split_monic_cubic_roots_in_session      // §3
  d=4 → quartic_roots                          // §5
  d>4 → Err(NotImplemented)
  → canonical_sort_roots(dedup_roots(...))
```

---

## 3. 核心：`split_monic_cubic_roots_in_session`（C1/C2 + F2）

**Tier:** Pipeline private（唯一三次分裂入口）  
**薄包装:** `cubic_roots`、`resolvent_cubic_roots_in_session` → 仅调本节。

### 3.1 前置条件

| 字段 | 要求 |
|------|------|
| **输入** | `p` monic，`deg(p,var)=3`，系数 ∈ `session.ambient()` |
| **输出** | 3 根（dedup 后 ≤3），`canonical_sort_roots` |
| **后验** | 每根 `verify_root`；`d_L ≤ D_HARD`；典型 \(d_K=1\) 时 `d_L ≤ D_CUBIC_SPLIT` |

### 3.2 算法（normative，与 `poly_roots.rs` 一致）

```text
split_monic_cubic_roots_in_session(session, p, var):

  // --- 分支 A：仿射纯三次 x³ + a₁x + a₀ 且 a₂=a₁=0 ---
  if coeff(p,var,2)=0 and coeff(p,var,1)=0:
    return finish(f2_pure_cubic_roots_in_session(session, coeff(p,var,0)))   // §3.3 A′

  (p_dep, q_dep, shift, _) ← cubic_depressed_parts(session, p, var)

  // --- 分支 B：depressed t³ + q_dep（p_dep=0）---
  if p_dep = 0:
    rs ← f2_pure_cubic_roots_in_session(session, q_dep)
    return finish(rs 每根 + shift)

  // --- 分支 C：一般 / casus / resolvent（p_dep ≠ 0）---
  return f2_split_cubic_via_deflate_in_session(session, p, var, p_dep, q_dep, shift)   // §3.4

finish(rs):
  if d_L > D_HARD: Err(NotImplemented("poly roots dim bound"))
  return canonical_sort_roots(dedup_roots(rs))
```

**数学注（F2 A′ 修正，2026-06-23）：** 对 depressed \(t^3+p t+q\) 且 \(p\neq 0\)，若 \(z_0\) 是多项式的一根，则 **另两根一般不是** \(\omega z_0\)、\(\omega^2 z_0\)（反例：resolvent \(z^3-4z-1=0\)）。Cardano casus 的三实根是 \(u\omega^k+v\omega^{-k}\) 型，**不是** \(\omega^k z_0\)。故 **A′ ω 分支仅限** \(p_{\mathrm dep}=0\)（已在分支 A/B 覆盖）。casus resolvent 与 `x³−x+1` 均走分支 C。

### 3.3 `f2_pure_cubic_roots_in_session`（F2 A′ — ω）

```text
f2_pure_cubic_roots_in_session(session, a0):   // monic t³ + a0
  β ← session.adjoin_cbrt(neg(a0))
  ω ← session.adjoin_primitive_cube_root_of_unity()
  return dedup [β, ω·β, ω²·β]
  ASSERT d_L <= d_K * D_PURE_CUBIC   // 典型 3
```

### 3.4 `f2_split_cubic_via_deflate_in_session`（F2 + F1）

**适用：** \(p_{\mathrm dep}\neq 0\) — 含 Cardano \(\Delta<0\) 的 resolvent、S3 不可约三次（如 `x³−x+1`，\(\Delta>0\) 但 Gal=S3）。

```text
f2_split_cubic_via_deflate_in_session(session, p, var, p_dep, q_dep, shift):
  z0 ← adjoin_one_cubic_root_depressed(session, p_dep, q_dep, shift)
       // = casus_adjoin_cubic_root：不可约 t³+p_dep·t+q 叠一层，禁止 Cardano u+v
  quad ← deflate_monic(session, p, var, z0)
  (r1, r2) ← quadratic_roots_formula(session, quad, var)   // §4，√Δ 仅 sqrt_in_field
  return finish [lift(z0), r1, r2]
  ASSERT d_L <= D_CUBIC_SPLIT * scale   // x³−x+1、resolvent 实测 ≤6
```

**禁止：** 三根分裂路径调用 `one_cubic_root`（Cardano）；禁止对 \(p_{\mathrm dep}\neq 0\) 使用 `adjoin_primitive_cube_root_of_unity` 凑第三根。

### 3.5 `one_cubic_root`（单根 API，非分裂入口）

**Tier:** Pipeline private — **仅**单根测例 / 四次前探；**不**用于 `split_monic_cubic_roots_in_session`。

```text
one_cubic_root(session, p, var):
  // 可保留 Cardano / casus_adjoin 单支逻辑（求一根）
  // 禁止在此求第二、第三根
```

### 3.6 与 resolvent 的关系

```text
resolvent_cubic_roots_in_session(session, R, z_var):
  return split_monic_cubic_roots_in_session(session, R, z_var)

quartic_roots 在 resolvent 后:
  debug_assert!(d_L <= POLY_ROOTS_DIM_RESOLVENT)   // d_K=1 → 6
```

---

## 4. 二次：`quadratic_roots_formula`（F1 收紧）

**前置:** monic \(x^2+bx+c\)，系数 ∈ `session.working()`。

```text
quadratic_roots_formula(session, p, var):
  Δ ← b² - 4c
  sqrt_Δ ← sqrt_disc(session, Δ)
  return [(-b + sqrt_Δ)/2, (-b - sqrt_Δ)/2]
```

**`sqrt_disc`:**

```text
sqrt_disc(session, u):
  if is_negative_rational(u): return mul_formal_i(sqrt_in_field(abs(u)))
  return sqrt_in_field(u)
```

---

## 5. 四次：`quartic_roots`（F2 resolvent + F3 Euler）

### 5.1 子情形

| 子情形 | 条件 | 算法 | \(d_L\) 上界（\(d_K=1\)） |
|--------|------|------|---------------------------|
| 双二次 | \(a_3=a_1=0\) | `biquadratic_roots` | **4** |
| 一般 | else | resolvent §3 + Euler §5.2 | **\(\le 24\)**（`D_HARD`）；若 \(\mathrm{Gal}\cong A_4\) 则 \([K_f:\mathbb{Q}]=12\) 为**数学推论**，非未证时的实现目标 |

### 5.2 一般四次

```text
quartic_roots(session, p, var):
  dep ← depress_quartic(...)
  R ← build_resolvent_cubic(...)
  z_roots ← split_monic_cubic_roots_in_session(session, R, _z)   // §3
  ASSERT d_L <= D_RESOLVENT
  rs ← euler_depressed_quartic_roots(session, dep, var, z_roots)
  ASSERT d_L <= D_HARD                    // 出口门禁（§5.4）
  return rs   // 4 根，verify 每个
```

### 5.3 Euler（F3 确定性）

```text
(α, β, γ) ← canonical_sort 后 z_roots[0..3]
√α ← sqrt_in_field(α);  √β ← sqrt_in_field(β)
√γ ← sqrt_in_field(γ) 或 euler_derived_sqrt_gamma(q≠0)

固定四行 ε ∈ {(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)}
x ← (ε₁√α + ε₂√β + ε₃√γ)/2 + shift
```

**禁止:** `sqrt_principal` 生产路径（已改为 `sqrt_in_field`）；`permutations` / `approx_real_sign` / 搜索环。

### 5.4 维数门禁与数学含义（F4）

**Rust 落点:** `quartic_roots` 末尾、`sqrt_in_field_euler_second`、`sqrt_in_field`（`field_session.rs`）。

| 常量 | 值 | 阶段 | 违反时 |
|------|-----|------|--------|
| `D_RESOLVENT` | 6 | resolvent 三根输出后 | `debug_assert` / F2 回归 |
| `D_HARD` | 24 | **出口门禁**：`quartic_roots` 返回前；`sqrt_in_field*` blind adjoin 前 | `NotImplemented("poly roots dim bound")` |
| `D_QUARTIC_A4_REF` | 12 | **无门禁**；仅作 \(\mathrm{Gal}\cong A_4\) 已证时的参考 / 分测 | — |

**出口门禁（normative）：**

```text
if dim(session.working()) > D_HARD:   // D_HARD = 24
  Err(NotImplemented("poly roots dim bound"))
```

**数学含义（系数域 \(K=\mathbb{Q}\)，\(d_K=1\)）：**

1. **\(d_L = \dim_{\mathbb{Q}} L\)** 是当前 Session 工作域 **L** 作为 \(\mathbb{Q}\)-向量空间的维数。单调扩域下 \(d_L\) 等于已 adjoin 的代数层次数之积（在实现采用的塔表示下）。

2. **不可约四次 \(f\) 的分裂域** \(K_f/\mathbb{Q}\) 满足 \([K_f:\mathbb{Q}] \mid 24\)，因为 \(\mathrm{Gal}(f)\) 是 \(S_4\) 的子群。上界 **24** 对应 \(\mathrm{Gal}(f)\cong S_4\)（“全对称”四次）。

3. **\(d_L > 24\) 拒绝输出** 的语义：在 **F2 resolvent + F3 Euler + F1 sqrt_in_field** 管线下，若 Session 已超过 \(S_4\) 分裂域的最大次数，则视为 **错误叠塔**（多余 blind adjoin、未识别的已有平方根等），**不是** “数学上需要更大域” 的合法情形——合法四次实代数根应能在 \([L:\mathbb{Q}]\le 24\) 内表达。

4. **\(d_L \le 24\) 允许输出** 并不保证 \(L = K_f\)（分裂域）；只保证维数未超过四次绝对上界。实现仍要求四根 `verify_root`。

5. **`D_QUARTIC_A4_REF = 12`（非通用目标）：** 当且仅当 \(\mathrm{Gal}(f)\cong A_4\) 时，分裂域满足 \([K_f:\mathbb{Q}]=12\)。**不得**在未证明 \(\mathrm{Gal}(f)\) 时，把 12 写成全体四次（含 `t⁴+t+1`）的 normative 目标。`POLY_ROOTS_DIM_QUARTIC_OUT` 保留为 Rust 常量名，语义同本参考值。

6. **`t⁴+t+1` 实测 \(d_L=24\)：** 与 \(\mathrm{Gal}(f)\) 为 \(S_4\) 阶子群、\([K_f:\mathbb{Q}]=24\) 的**工作假设**一致（仓库内**未**形式证明 Gal 群）。当前 Euler 路径在 resolvent 后 \(d_L=6\)，两次二次 blind adjoin 后 \(d_L=24\)；四根 `verify_root` ✅，在 `D_HARD` 内。**不**视为未达「12 目标」的回归失败。

7. **与 `sqrt_in_field_euler_second` 的关系：** 第二 Euler 开方先 `try_sqrt(v)`、`try_sqrt(v/u)`、F4′ Galois；仅在 \(d_L < D_HARD\) 时允许对 \(v\) 再 blind adjoin。

**测例绑定：** `field_session_dimension_bound_quartic` / `_tight` 断言 `t⁴+t+1` 全流程 `d_L ≤ 24` 且四根 verify。

---

## 6. `sqrt_in_field` / `try_sqrt_in_field`（F1）

**落点:** `field_session.rs` + `ext_tower::try_square_root_in_field`。

```text
sqrt_in_field(session, u):
  u ← lift(u)
  if let Some(ε) = try_sqrt_in_field(u): return ε
  if dim(L) >= D_HARD: Err(NotImplemented("poly roots dim bound"))
  adjoin_sqrt_new(u); assert ε² eq_mod u

try_sqrt_in_field:
  1a  塔层 deg-2 生成元 scan（含 k·g, k≤8；dim>8 时 k≤4 ponytail）
  1b  operational 基 eᵢ，试 eᵢ² ≡ u
  1c  dim ≤ 12：{0,±1} 两元线性组合
  1d  dim ≤ 9：{0,±1} 三元线性组合
  1e  dim ≤ 6：quadratic-gen 小组合（F4）
  miss → None

Euler 第二开方（`sqrt_in_field_euler_second`）：仅 1a+1b（shallow），miss 则 `adjoin_sqrt_new`（跳过 1c–1e 与 `sqrt_in_field` 重复全搜）。**interim**；终态 §12.2。
```

---

## 7. API 映射（Rust 落点，2026-06-23）

| 规格函数 | Rust | 状态 |
|----------|------|------|
| `split_monic_cubic_roots_in_session` | `poly_roots.rs` | ✅ |
| `f2_pure_cubic_roots_in_session` | 同上 | ✅ |
| `f2_split_cubic_via_deflate_in_session` | 同上 | ✅ |
| `adjoin_one_cubic_root_depressed` | → `casus_adjoin_cubic_root` | ✅ |
| `resolvent_cubic_roots_in_session` | 薄包装 §3 | ✅ |
| `finish_cubic_roots` | `POLY_ROOTS_DIM_HARD` 门禁 | ✅ |
| `one_cubic_root` | 单根 only | ✅ |
| `is_cubic_casus` → ω·z₀ | — | **废弃**（§3.2 注） |
| `quadratic_roots_formula` / `sqrt_disc` | F1 路径 | ✅ |
| `euler_depressed_quartic_roots` | `sqrt_in_field` | ✅ |
| `POLY_ROOTS_DIM_*` | `field_session.rs` | ✅ |

---

## 8. 验收测例（与参数绑定）

| 测例 | 断言 | 状态 |
|------|------|------|
| `roots_cubic_t3_minus_2` | 3 根；A′ ω | ✅ |
| `roots_x3_minus_x_plus_1_vanish` | 3 根 verify；分支 C | ✅ |
| `field_session_dimension_bound_cubic_x3_minus_x_plus_1` | `d_L ≤ 6` | ✅ |
| `resolvent_cubic_all_roots_vanish` | 3 根；`d_L ≤ 6` | ✅ |
| `f2_resolvent_split_no_cardano_stack` | resolvent 无 Cardano 叠层 | ✅ |
| `resolvent_dim_bound_t4_plus_t_plus_1` | `d_L ≤ 6` | ✅ |
| `adjoin_sqrt_after_quadratic_roots_formula` | F1 ε²=u | ✅ |
| `adjoin_sqrt_after_sqrt_disc_on_resolvent` | F1 | ✅ |
| `roots_quartic_t4_plus_t_plus_1` / `euler_four_roots_vanish` | 4 根 | ✅ |
| `field_session_dimension_bound_quartic` | `d_L ≤ 24` | ✅ |
| `field_session_dimension_bound_quartic_tight` | `t⁴+t+1`：`d_L ≤ 24` + 四根 verify | ✅ |

**命令:** `cargo test -p giac-core poly_roots::tests --release` — **0 ignore**（2026-06-23：**25 passed**）。

---

## 9. 禁止清单（算法层）

| ❌ | ✅ |
|----|-----|
| 三根分裂走 Cardano `u+v` + 多次 cbrt | `adjoin_one_cubic_root_depressed` + deflate |
| \(p_{\mathrm dep}\neq 0\) 时 \(\omega^k z_0\) | 分支 C：deflate + `quadratic_roots_formula` |
| resolvent / cubic 两套分裂逻辑 | 统一 `split_monic_cubic_roots_in_session` |
| `algext_square_roots` in `poly_roots` 热路径 | `sqrt_in_field` |
| `adjoin_sqrt_new` 在 Cardano/Euler 多处 | 仅 `sqrt_in_field` fallback |
| `dim ≥ D_HARD` 仍 blind adjoin | `NotImplemented("poly roots dim bound")` |
| solve 层 workaround C1 | C1 已绿（S0 代入测例同步） |

---

## 10. 与 F1–F5 文档的差异

[GIAC-poly-quartic-roots-F1-F5.md](../issues/GIAC-poly-quartic-roots-F1-F5.md) §1.5 初稿写「casus resolvent → ω·z₀」。**本规格优先**（实现验证）：仅 **纯三次** 用 ω；**casus resolvent**（\(p\neq 0\)）走 §3.4。若需统一 F1–F5 正文，以本节为准修订该 issue。

---

## 12. 升级路径（plan）— 代数层替代 ponytail 枚举

**现状（interim，2026-06）：** F1 `try_sqrt_in_field` 用塔层扫描 + 有界 ±1 枚举（§6，dim 上界见 §1）；Euler 第二开方用 shallow + blind adjoin。根 **正确**（`verify_root`），但 S₄ 四次（如 `t⁴+t+1`）实测 `d_L=24`，枚举 **不是** 长期方案。

**终态（plan，不阻塞 P3-6 验收）：**

| 优先级 | 项 | 替代对象 | 落点 | 收益 |
|--------|-----|----------|------|------|
| **P1** | **F5 — 结构开方** | §6 的 1c–1e ±1 暴力枚举 | `ext_tower::try_square_root_in_field` | 按塔层/最小多项式/嵌入判定平方根；\(O(\text{层数}+\dim)\)，可删或极大缩小 enum ceiling |
| **P2** | **F4′ — Galois 第二 Euler 开方** | `sqrt_in_field_euler_second` 的 shallow+adjoin | `galois_automorphism.rs` + `poly_roots.rs` | 在 **已证或应存在** 于当前 L 的平方根时避免多余 adjoin；**不**承诺未证 Gal 时 `t⁴+t+1` 必 \(d_L=12\) |
| P3 | F5 收敛 | `FieldSession::try_sqrt_in_field` 薄封装 F5 | `field_session.rs` | 与 partfrac / 未来 gcd 共用 |

**依赖顺序：** F1 绿（✅）→ **F5 结构探测** → **F4′ Galois**（可并行设计，实现建议 F5 先行）→ 收紧/删除 §1 三行 `TRY_SQRT_*` interim 常量。

**验收（plan DoD）：**

- [ ] F4′：在 **Gal≅A₄ 已证** 的四次（或单独分测）上 \(d_L \le D_QUARTIC_A4_REF\) 且四根 verify；**或** 登记为何需要 \(d_L>12\)
- [ ] `t⁴+t+1`：保持 \(d_L \le D_HARD\) + 四根 verify（**现行 baseline**，非 12）
- [ ] `try_square_root` 无 dim³ 枚举热路径；`poly_roots::tests` 仍 ≤ `T_ROOTS`
- [ ] [GIAC-poly-quartic-roots-F1-F5](../issues/GIAC-poly-quartic-roots-F1-F5.md) §F5 / F4′ 与本节同步

**非 plan：** 继续压低 enum dim 上界换 CI 速度 — 仅 interim，不叠加。

---

## 11. 变更日志

| 日期 | 变更 |
|------|------|
| 2026-06-23 | 初版：§3 统一三次分裂；D_* 参数表 |
| 2026-06-23 | **§12 plan：** F5 结构开方 + F4′ Galois 替代 ponytail 枚举（interim 枚举/shallow 非终态） |
| 2026-06-23 | **§5.4 收紧：** 删除无证据的全体四次 `d_L=12` 目标；`D_QUARTIC_A4_REF` 仅 Gal≅A₄ 已证时；`t⁴+t+1` baseline \(d_L\le 24\) |
| 2026-06-23 | **归档：** P3-6 normative 管线验收通过 → `issues_resolved/` |
