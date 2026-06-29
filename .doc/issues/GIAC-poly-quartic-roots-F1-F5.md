# GIAC-poly — 四次 `roots` 塔修复（F1→F5）

**状态:** open  
**类型:** 实现 / AFK 可抓取  
**父项:** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) **P3-6**（通用四次 `Poly<AlgExtC>::roots`）  
**相关:** [GIAC-algext-adoption](GIAC-algext-adoption.md) §8.2–8.4、[giac-tower-common-math.md](../giac-tower-common-math.md)、[known-divergences.md](../known-divergences.md) DIV-076–083 / DIV-081  
**Rust 落点:** `giac-core::algebra::{field_session, poly_roots, ext_tower}`  
**快照:** 2026-06-23（M1/M2 绿；F3–F4 绿；F5 结构开方进行中）

---

## 0. 问题陈述

`poly_roots::quartic_roots`（Euler + Ferrari resolvent）在 **resolvent 分裂后** 再 `adjoin_sqrt` 时破坏 \(\varepsilon^2=u\)（塔嵌入与已有平方根未识别）。表现（**2026-06-23 已修复**）：

- `adjoin_sqrt_squares_one_resolvent_root`、`adjoin_sqrt_after_resolvent_deflate_only` ✅
- `sqrt_disc` / `quadratic_roots_formula` 后再 adjoin ✅（F1 `sqrt_in_field`）
- `euler_four_roots_vanish`、`roots_quartic_t4_plus_t_plus_1` ✅（F3 单路径）

**根因（摘要）：** 盲目 `adjoin_irreducible(u²−α)` 叠塔，未在 **当前 L** 内查找已有平方根；resolvent 二次分裂又额外扩域，维数爆炸 + 生成元语义错位。

**目标：** `solve(t^4+t+1=0,t)` 四根 `eq_mod` 零化；release 下 `poly_roots::tests` 无 ignore、无分钟级穷举。**✅ 已达成**（25/25 release ~3s）。

### 0.1 进度快照（2026-06-23）

| 阶段 | 完成度 | 代码现状 |
|------|--------|----------|
| **F1** | ✅ **M1** | `sqrt_in_field` / `try_sqrt_in_field` → `ext_tower::try_square_root_in_field`；4 项 DoD 测绿 |
| **F2** | ✅ **M2** | `resolvent_cubic_roots_in_session` = `split_monic_cubic_roots_in_session`；纯三次 ω；p≠0 irreducible adjoin + deflate + F1 √Δ |
| **F3** | ✅ **M3** | Euler 单路径；无 flip / `approx_real_sign`；e2e 无 ignore |
| **F4** | ✅ **M4** | `field_session` 契约 doc；维数 bound 测（硬顶 24；resolvent ≤6） |
| **F5** | 🟡 | S3–S6 结构开方已落地；F4′ Galois（避免多余 adjoin，**不**绑未证的 d_L=12） |

---

## 1. 依赖顺序（必须按序）

```text
F1  sqrt_in_field（含 ε²=u 回归）
  ↓
F2  resolvent 三根留在 L₃
  ↓
F3  Euler 去穷举 + e2e 绿
  ↓
F4  session 契约 + 维数上界测
  ↓
F5  ext_tower 通用开方 API（可选）
```

| Issue | 阻塞 | 解除 |
|-------|------|------|
| **F1** | F2、F3、F4 部分 | — |
| **F2** | F3 | F1 |
| **F3** | P3-6 验收、e2e | F1 + F2 |
| **F4** | 长期回归门禁 | F1（契约与 F2 维数断言可并行，建议 F1 后） |
| **F5** | 无（可选重构） | F1 语义稳定后 |

---

## 1.5 实现方案复审（2026-06-22）

> 复审结论：**依赖顺序不变**；F1 算法可分两阶段落地；F2 应走 **A′（casus + ω 共轭）** 而非裸 `quadratic_roots_formula`；F3 在 F1+F2 绿后主要是删代码；F4 可与 F3 同 PR；F5 推迟到 F1 稳定后。

### F1 — `sqrt_in_field` 算法（推荐）

**入口契约（normative）：**

```text
sqrt_in_field(session, u) -> Result<ε, Err>
  1. try_sqrt_in_field(u)  // Option
  2. 若 Some(ε): lift(ε); assert mul(ε,ε) eq_mod u; return ε
  3. 否则 adjoin_irreducible(γ²−embed(u)) via algext_square_roots; assert ε²=u; return ε
```

**`try_sqrt_in_field` 两阶段（先易后全）：**

| 阶段 | 做法 | 复杂度 | 覆盖 |
|------|------|--------|------|
| **1a 塔层扫描** | 沿 `ExtensionTower::Adj` 链：每层 deg-2 生成元 \(g\) 满足 \(g^2=\alpha\)；若 `u eq_mod α` 返回 \(\pm g\) | \(O(\text{层数})\) | 根式管线中绝大多数 adjoin |
| **1b 基元平方** | 对 operational 基 \(\{e_i\}\)：若 `mul(e_i,e_i) eq_mod u` 返回 \(e_i\) | \(O(d)\) | 复合坐标下的平方 |
| **1c 有界枚举**（ponytail 天花板） | `dim(L)≤B`（四次管线 **B=24**）时，对基线性组合系数 \(\in\{0,\pm1\}\) 的 \(v\) 试 \(v^2\) | \(O(3^d)\) 最坏；仅 fallback | 漏检时升级路径 |

**落点：** 探测核心先写 `ext_tower::try_square_root_in_field(field, coords)`（F5 前驱）；`FieldSession::try_sqrt_in_field` 薄封装 + `lift`。`algext_square_roots` **保留**为「必定叠层」底层，与 `try_*` 分工（同 F5 文档意图）。

**`sqrt_disc` 简化：** 删除「adjoin 后验失败再 adjoin \(-\Delta\)」；统一 `sqrt_in_field(Δ)`，虚部仅 `is_negative_rational` + `mul_formal_i`。

**反模式（禁止）：** 加大 `approx_real_sign` 搜索；在 F1 未绿前 unignore e2e；用 `cubic_roots` 代替 F2 新 API 凑绿。

### F2 — resolvent 三根（推荐路径 **A′**）

**数学修正：** \(t^4+t+1\) 的 resolvent \(R(z)=z^3-4z-1\) 为 **casus**（\(\Delta_R<0\)）。单根域 \(\mathbb{Q}(z_0)\) 次数 **3**；另两根 **不在** 该域，需 resolvent **分裂域** \(\mathbb{Q}(z_0,\omega)\)，次数 **6**。故「L₃」在规格上指 **resolvent 分裂域**，非「仅一条 Cardano 根」。

**推荐实现 `resolvent_cubic_roots_in_session`：**

```text
z0 = one_cubic_root(R)
若 R 为 casus（Δ<0 或 deflate 二次 Δ 不在 Q(z0) 内平方）:
  ω = adjoin_primitive_cube_root_of_unity()   // 已有 API
  z1 = ω·z0, z2 = ω²·z0   // 域内 mul，无新 adjoin
否则（Cardano 实根 + 分裂已在 L）:
  quad = deflate_monic(R, z0)
  z1,z2 = quadratic_roots_formula 但 √Δ 走 sqrt_in_field（F1）
返回 dedup [z0,z1,z2]；三根 verify_root(R)
```

- **不**在 resolvent 阶段调用全 `cubic_roots`（会 `quadratic_roots` adjoin 叠塔）。
- **不**对 casus 走 `quadratic_roots_formula`+blind √Δ（当前 `quartic_roots` 仍犯此错）。
- 与 `pure_cubic_roots` / `roots_cubic_t3_minus_2` 同模式，复用 ω 分支。

**维数上界（F2/F4，有证据者 normative）：**

| 阶段 | `t^4+t+1` 预期 `dim(working)` | 依据 |
|------|-------------------------------|------|
| resolvent 结束 | **≤6** | F2 规格 + 实测 |
| 全流程出口 | **≤24** | \(\mathrm{Gal}(f)\le S_4 \Rightarrow [K_f:\mathbb{Q}]\mid 24\)；`D_HARD` |
| Gal≅A₄ 时 | **12** | **数学推论** \([K_f:\mathbb{Q}]=12\)；须 **单独证明 Gal** 后才可作分测目标，**非** `t⁴+t+1` 默认 |

PR 中记录实测 `dimension()`；超过 `D_HARD` 必失败。

### F3 — Euler 确定性分支（F1+F2 后）

1. **Resolvent 排序：** `z_roots` 按固定键排序（建议：先 `eq_mod` 去重，再按 `flatten_min_poly_over_q` 字典序或 Vieta \((z_0+z_1+z_2, z_0 z_1 z_2)\) 二元组）；固定 \((\alpha,\beta,\gamma)=(z_0,z_1,z_2)\)。
2. **平方根：** `pick_sqrt_euler` → `sqrt_in_field`；仅当 \(u\in\mathbb{Q}\) 且 \(u<0\) 时 `i·sqrt(|u|)`；**删除** `approx_real_sign` 生产路径。
3. **符号表：** 固定四行 \((\varepsilon_1,\varepsilon_2,\varepsilon_3)\in\{(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)\}\)；\(q\neq0\) 时 \(\sqrt\gamma=-q/(\sqrt\alpha\sqrt\beta)\)（已有 `euler_derived_sqrt_gamma`）。
4. **删除** `permutations` / `flip_*` / `set_working` 搜索环；`euler_depressed_quartic_roots` 单路径调用。

### F4 / F5

- **F4：** 与 F3 同 PR；`field_session.rs` 模块 doc 补 §契约；`field_session_dimension_bound_quartic` 断言 **硬顶 24**（`t⁴+t+1` baseline）。
- **F5：** F1 合并后 refactor；**不阻塞** P3-6。

### 下游接线（F3 后，非本 issue 阻塞）

| 项 | 文档 | 说明 |
|----|------|------|
| `solve(t^4+t+1)` | [GIAC-rs-crate-dedup-plan](GIAC-rs-crate-dedup-plan.md) D3-1 | `eval_solve` 已 wired；e2e 随 F3 unignore |
| `general quartic rootof` | [P3-6 gaps](GIAC-poly-p3-6-quartic-roots-gaps.md) | `rootof.rs` 仍 `NotImplemented` |

---

## F1 — `sqrt_in_field` / `try_existing_square_root`

**状态:** ✅ **M1 完成**

### 问题

`FieldSession::adjoin_sqrt` 总是 adjoin 新元；当 \(u\) 已在 **L** 中为某元平方时，应返回 **嵌入 L 的** \(\varepsilon\) 且保证 \(\varepsilon^2 \equiv u\)，而非新层 \(\varepsilon'\) 满足 \((\varepsilon')^2\) 在错误子域/mod 下不等于 \(u\)。

### 范围

| 做 | 不做 |
|----|------|
| `try_sqrt_in_field(u) -> Option<AlgExtCPolyCoeff>`（或 `sqrt_in_field` 先 try 再 adjoin） | 一般 \(n\) 次根（仅平方） |
| `adjoin_sqrt` 改为：`try` → 成功则 `lift` + 验 \(\varepsilon^2=u\)；失败再 adjoin | 改 giac C++ `common_EXT` 行级行为 |
| `sqrt_principal` / `sqrt_disc` / `pick_sqrt_euler` 统一走该入口 | F2 resolvent 三根算法（留给 F2） |

### 任务

1. 在 `ext_tower` 实现 `try_square_root_in_field`（§1.5 阶段 1a→1b→1c）；`field_session::try_sqrt_in_field` / `sqrt_in_field` 薄封装。
2. `FieldSession::adjoin_sqrt`：`try_sqrt_in_field` → 验 `mul(ε,ε) eq_mod u` → 否则 `adjoin_irreducible(u²−α)`。
3. `sqrt_disc`：去掉「adjoin 后验失败再盲 adjoin \(-\Delta\)」的脆弱双路径，优先 `sqrt_in_field`。
4. 扩展回归：**在** `quadratic_roots_formula` / `sqrt_disc` **之后** adjoin 仍满足 ε²=u（当前失败场景）。

### 验收（DoD）

- [x] 新测 `adjoin_sqrt_after_quadratic_roots_formula`（或等价名）通过
- [x] 新测 `adjoin_sqrt_after_sqrt_disc_on_resolvent`（resolvent 一条根 + formula 二次后 adjoin）通过
- [x] 既有 `adjoin_sqrt_squares_one_resolvent_root`、`adjoin_sqrt_after_resolvent_deflate_only` 仍绿
- [x] `cargo test -p giac-core poly_roots::tests --release` 无新增 ignore

### 文件

- `giac-core/src/algebra/field_session.rs`
- `giac-core/src/algebra/ext_tower.rs`（若探测在塔层）
- `giac-core/src/algebra/poly_roots.rs`（`sqrt_disc`、测试）

### 偏离

- 闭合 [DIV-081](../known-divergences.md)「实现缺口」子项；[DIV-077](../known-divergences.md) 复审

---

## F2 — resolvent 三根留在 L₃（无二次分裂叠塔）

**状态:** ✅ **M2 完成**  
**依赖:** F1（`quadratic_roots_formula` 分支需 `sqrt_in_field`）

### 问题

当前 `quartic_roots`：`one_cubic_root` + `deflate_monic` + `quadratic_roots_formula` 得到 **三条** resolvent 根，但 formula 路径会 **再次 adjoin** √Δ，破坏 L₃ 上 Euler 所需坐标系。全 `cubic_roots` 亦 over-split。

### 范围

| 做 | 不做 |
|----|------|
| 在 **L₃**（一条 Cardano/casus 根 + 同场另外两根）内得到 resolvent 三零点 | 替换 Cardano 为完全不同三次算法（除非必要） |
| 优先：**Galois 共轭** / **Vieta 对称** / **deflate 后域内求根**（不 `quadratic_roots` adjoin） | F3 Euler 符号表（留给 F3） |
| `quartic_roots` 仅消费 `[z0,z1,z2] ⊂ L₃` | 全次数 `roots` 重构 |

### 任务

1. 设计 API：`resolvent_cubic_roots_in_session(session, R) -> [z0,z1,z2]`，**不** 使 `session.working()` 超过「一条 casus + 必要 ω」的约定维数上界（见 F4）。
2. **推荐路径 A′**（§1.5）：casus → `adjoin_primitive_cube_root_of_unity` + \(z_k=\omega^k z_0\)；非 casus → deflate + `quadratic_roots_formula` 且 √Δ 仅经 `sqrt_in_field`。PR 说明所选分支。
3. `quartic_roots` 改调 `resolvent_cubic_roots_in_session`；删除对裸 `quadratic_roots_formula` / `cubic_roots` 的 resolvent 依赖。
4. 单测：`resolvent_cubic_all_roots_vanish` 在 **不调用** `quadratic_roots_formula` 下仍绿；三根 `eq_mod` 互异且零化 \(R\)。

### 验收（DoD）

- [x] `resolvent_cubic_all_roots_vanish` 走 `resolvent_cubic_roots_in_session`
- [x] `session.working().dimension()` 在 `t^4+t+1` resolvent 阶段 ≤ **6**（`resolvent_dim_bound_t4_plus_t_plus_1`）
- [x] F1 回归全套仍绿
- [x] `x³−x+1` 三次：`dim≤6` + verify（`f2_one_cubic_root_for_deflate` 禁止 Cardano 叠塔）

### 文件

- `giac-core/src/algebra/poly_roots.rs`（`one_cubic_root`、`deflate_monic`、`quartic_roots`）
- `giac-core/src/algebra/field_session.rs`（必要时 `adjoin_primitive_cube_root_of_unity` 契约）

### 偏离

- [DIV-080](../known-divergences.md)、[DIV-081](../known-divergences.md)

---

## F3 — Euler 固定分支、去穷举、e2e 绿

**状态:** ✅ **M3 完成**  
**依赖:** F1 + F2

### 问题

`euler_depressed_quartic_roots`：**6 置换 × flip 搜索** + `pick_sqrt_euler` 数值 `approx_real_sign`。忽略测 `euler_four_roots_vanish`、`roots_quartic_t4_plus_t_plus_1`（当前 `--ignored` 秒失败 `NotImplemented("quartic euler")`，非慢搜）。

### 范围

| 做 | 不做 |
|----|------|
| 固定 **canonical** \((\alpha,\beta,\gamma)\) 赋值（resolvent 三根排序/Vieta） | 6×flip 穷举 |
| 固定 Euler 符号 \((\varepsilon_1,\varepsilon_2,\varepsilon_3)\) 表（\(\varepsilon_1\varepsilon_2\varepsilon_3=+1\)） | 改 Ferrari 为 Depressed+Biquadratic 特判表 |
| \(\sqrt\gamma = -q/(\sqrt\alpha\sqrt\beta)\) 在 L 内用 F1 | 追 giac 四根字面（见 DIV-080） |
| 去掉 `#[ignore]`；`poly_algext_roots(t^4+t+1)` e2e | `giac-solve` P4-6 全管线（可 follow-up） |

### 任务

1. `euler_depressed_quartic_roots`：删除 `permutations` / `flip_*` 循环；单路径 `euler_four_roots_from_triple`。
2. `pick_sqrt_euler`：在 F1+F2 下改为 **确定性** 分支（例如：先试实 `sqrt_in_field`；负则 `i·sqrt(|u|)`；禁止 deg≤6 数值搜索，或仅作 `debug_assert` fallback）。
3. 测 `euler_gamma_relation_holds`：\(\sqrt\alpha\sqrt\beta\sqrt\gamma + q \equiv 0\)。
4. 取消 ignore：`euler_four_roots_vanish`、`roots_quartic_t4_plus_t_plus_1`。
5. release 全 suite：`cargo test -p giac-core poly_roots::tests --release` **<5s**（无 hang）。

### 验收（DoD）

- [ ] `euler_four_roots_vanish` enabled 且绿
- [ ] `roots_quartic_t4_plus_t_plus_1` enabled 且绿（`poly_algext_roots` 四根 + `verify_root`）
- [ ] **P3-6** 核心：`solve(t^4+t+1=0,t)` 数学规格满足（crate 内 roots 层；solve 接线可另 PR）
- [ ] 无 6×flip 搜索代码路径（或 `#[cfg(test)]` 对照仅）

### 文件

- `giac-core/src/algebra/poly_roots.rs`
- `giac-core/src/algebra/field_session.rs`（`sqrt_principal`）

### 偏离

- [DIV-077](../known-divergences.md)、[DIV-080](../known-divergences.md) 标 **accepted**（Rust 自研规格）

---

## F4 — `FieldSession` 契约 + 维数上界测

**状态:** ✅ **M4 完成**  
**依赖:** F1（可与 F3 并行落地）

### 问题

`FieldSession` 的 **ambient / working / checkpoint** 语义分散；四次管线曾维数爆炸（3→6→12→24）。需可审计契约，防止回归再次引入 blind adjoin 或 session 泄漏。

### 范围

| 做 | 不做 |
|----|------|
| 文档化 + 测试 **Session 不变量** | 重写整个 `ext_tower` |
| `set_working` / `bump_to` / `lift` 契约（何时允许回退、checkpoint 合法性） | F5 公共 API 形态 |
| 维数上界断言：`t^4+t+1` 全流程 `dim(working) ≤ B` | 证明 assistant |

### 任务

1. 在 `field_session.rs` 模块 doc 或 `.doc/issues/` 短节写 **契约**：
   - `ambient` 不变；`working` 单调扩域（除 `set_working(checkpoint)`）；
   - `lift`/`align` 后两操作数同 `working`；
   - `adjoin_*` 必须经 F1 `sqrt_in_field` 或显式不可约 adjoin。
2. 测试 `field_session_dimension_bound_quartic`：`t^4+t+1` 从 ℚ 到四根输出，`working.dimension()` ≤ **24**（resolvent ≤6；见算法规格 §5.4）。
3. 测试 `set_working_restores_adjoin`：checkpoint 后 adjoin 不污染已放弃分支（Euler 删除穷举后可为 smoke）。
4. （可选）`FieldSession::checkpoint()` / `restore()` 类型安全包装，替代裸 `set_working`。

### 验收（DoD）

- [ ] 契约文字入 doc + `field_session` 模块注释
- [ ] 维数 bound 测在 CI 默认路径绿
- [ ] F3 e2e 绿时维数测仍绿

### 文件

- `giac-core/src/algebra/field_session.rs`
- `giac-core/src/algebra/poly_roots.rs`（维数测）
- 可选：[giac-tower-common-math.md](../giac-tower-common-math.md) §4.3 增一行 P3 优先级

### 偏离

- [DIV-081](../known-divergences.md) compositum/coords 字面仍 open；维数 **数学** 正确性见 [giac-tower-common-math.md](../giac-tower-common-math.md)

---

## F4′ — Euler 第二开方 via Galois（plan，依赖 F5 或最小自同构 API）

**状态:** 🟡 部分（`galois_automorphism.rs` + `sqrt_in_field_euler_second` 已接线；二次层 σ(κ) **open**）  
**依赖:** F3 绿（✅）；F5 或 `ExtensionField` 上 **Gal(L/K)** 作用

### 问题

Euler 第二开方在 shallow/ratio 失败后可能 **blind adjoin**，使 Session 维数超过分裂域所需（多余叠塔）。F4′ 在 **当前 L** 内用 resolvent 根的 Galois 共轭找已有平方根，避免能识别的重复 adjoin。

**不声称：** 对未证明 \(\mathrm{Gal}(f)\) 的 `t⁴+t+1`，全流程必 \(d_L=12\)。实测 \(d_L=24\) 与 \(S_4\) 型分裂域 \([K_f:\mathbb{Q}]=24\) 的工作假设一致，在 `D_HARD` 内合法（算法规格 §5.4）。

### 任务（plan）

1. `sqrt_in_field_euler_second`：shallow/ratio miss 后试 `galois_automorphism::try_galois_sqrt_second`，再 blind adjoin。
2. 完成二次扩张上 σ(κ) 延拓（当前缺口）。
3. **分测：** 仅在 **Gal≅A₄ 已证** 的四次上断言 \(d_L \le 12\)；`t⁴+t+1` 保持 \(d_L \le 24\) baseline。

### 验收（DoD）

- [ ] Gal≅A₄ **已证** 分测：\(d_L \le 12\) + 四根 verify（或文档说明例外）
- [ ] `t⁴+t+1`：\(d_L \le 24\) + 四根 verify（**现行** `_tight` 测已覆盖）
- [ ] 删除或降级无收益的 Euler 双 adjoin interim 路径（有证据时）

**规格：** [GIAC-poly-p3-6-roots-algorithm-spec.md](GIAC-poly-p3-6-roots-algorithm-spec.md) §12（plan）；normative 见 [issues_resolved](../issues_resolved/GIAC-poly-p3-6-roots-algorithm-spec.md)。

---

## F5 — `ext_tower` 结构开方 API（plan：替代 F1 枚举）

**状态:** 🟡 **M5 部分完成**（F5 S0–S6 + S4 全量 ✅；F4′ Galois σ 待办；`t⁴+t+1` 实测 dim=24）  
**依赖:** F1 语义稳定（✅）

### F5 reframe（2026-06-27）— √Δ 探测的真正缺口定位

> **结论摘要：** 原「fix √Δ_C probe in dim-4 ℚ(α) for A₄→12」前提数学上不可能 — A₄ dim 12 = 4·3，故 √Δ_C ∉ dim-4（已验证：`N(α³−α²+9)=81` 是平方但 ℚ(α) 内无根；dim-4 probe 正确返回 None）。**真正缺口**是 dim-12 塔域 ℚ(α,β) 的 √Δ_Q 探测（现返回 None，应为 Some）→ 盲 adjoin → 24。该塔探测归约为 flat dim-4 ℚ(α) 内 √(·) 子问题。

对 A₄ 四次 `x⁴+8x+12`（resolvent `z³−48z−64`，disc=576²）实测 `diag_a4_sqrt_probe_gap`：

| 探测 | 域 | `try_sqrt_in_field` | 数学期望 | 结论 |
|------|----|--------------------|----------|------|
| √Δ_C | dim-4 ℚ(α) | **None** | None（正确） | A₄ 分裂域 dim 12 = 4·3，**非** 4·2 ⇒ √Δ_C ∉ ℚ(α)（否则四根全在 ℚ(α) ⇒ dim≤4，矛盾）。验证：`N(α³−α²+9)=81`（平方，必要条件满足）但暴力搜 [−16,16] 无整数根 ⇒ 确不在 dim-4。**dim-4 探测正确返回 None，无需修复。** |
| √Δ_Q | dim-12 ℚ(α,β) | **None** | **Some**（应为分裂域本身） | **真正缺口**：adjoin-deflate 已达 dim 12，但 F5 启发式在 dim-12 **塔域**漏检 √Δ_Q ⇒ 盲 adjoin ⇒ 24。 |

**结论：** 原「补 F5 √Δ_C 在 dim-4 ℚ(α)」前提错误；A₄→12 的杠杆是 **dim-12 塔域 ℚ(α,β) 的 √Δ_Q 探测**。该塔探测的 x₀ 坐标递归归约为 **flat dim-4 ℚ(α) 内 √(·)** 子问题（3 个 quadric in 3 未知数 over ℚ(α)）。

**已落地（giac-groebner）：** `groebner_basis_lex`（Buchberger + Gebauer-Möller product/chain criteria），小例正确（1077 测绿）。但 **lex 序对 x²=u 的 4 坐标 quadric 中间度爆炸**（120s 未完成；rank-deficient 系统即使 generic Vandermonde 换坐标亦然）。flat-over-ℚ solver `proto_try_sqrt_flat_over_q` 正确但过慢，待 **FGLM**（grevlex GB → lex 线性代数转换）提速。

**下一步选项：** (A) FGLM：grevlex Buchberger + FGLM 转换（~250 行，稳）；(B) 3-var-over-ℚ(α) 定向消元（resultant/GCD，比 4-var 更可控）；(C) 其它。`proto_try_sqrt_flat_over_q` 及其 helper（`generic_vandermonde`/`substitute_linear`/`proto_subst`/`proto_rational_roots`）保留为 WIP（`#[allow(dead_code)]`）。

**已选定并立项：** 方案 (A) FGLM over 系数域 `F=ℚ(α)`。详细 plan 与优先级见 [GIAC-poly-f5-fglm-over-coefficient-field](../issues_resolved/GIAC-poly-f5-fglm-over-coefficient-field.md)（P0 groebner 泛型 `C: FieldCoeff` → P1 grevlex → P2 FGLM → P3 塔域 S7 探测 → P4 回归；~4–5d）。


### 问题

F1 当前 `try_square_root_in_field` 在 layer/basis 扫描之后依赖 **有界 ±1 线性枚举**（§6 1c–1e，dim 上界见算法规格 §1）。这是 **ponytail interim**：漏检 → blind adjoin，根仍对但可能多余扩域、且 dim³ 路径曾拖慢 CI。

**终态（plan）：** 用 **塔结构 / 子域嵌入 / 最小多项式** 判定 \(u\) 是否为平方，**删除或仅保留极小 fallback 枚举**。详见 [GIAC-poly-p3-6-roots-algorithm-spec.md](GIAC-poly-p3-6-roots-algorithm-spec.md) §12；normative 背景见 [issues_resolved](../issues_resolved/GIAC-poly-p3-6-roots-algorithm-spec.md) §12。

### F5 实施计划（细项）

#### 目标 API（多模块共用）

```text
ExtensionField::try_square_root_in_field(field, u: CoordsQ) -> Option<CoordsQ>
  语义：若 ∃ε∈L, ε²≡u (mod minpoly)，返回其一；否则 None。绝不 adjoin。

FieldSession::try_sqrt_in_field / sqrt_in_field
  薄封装：lift → coords → try_square_root_in_field → AlgExtCPolyCoeff
  sqrt_in_field：None 且 dim<D_HARD → adjoin_sqrt_new

algext_square_roots(u)   // 保留
  语义：必定叠新层；display / 旧路径；roots 热路径禁止直接调用
```

**共用方（现状 → 终态）：**

| 模块 | 今日 | F5 后 |
|------|------|-------|
| `poly_roots` | `session.sqrt_in_field` / `try_sqrt_*` | 不变入口，底层换结构探测 |
| `field_session` | 委托 `ext_tower` | 仅 lift/align，无枚举逻辑 |
| `partfrac` / 积分塔（未来） | 可能直调 `algext_square_roots` | 改 `sqrt_in_field` 或 `try_*` |
| `giac-simplify/factor` | `from_rootof` | 不在 F5 范围 |

#### 「结构开方」= 五类探测（替代 1c–1e 枚举）

按序执行；命中即返回。复杂度目标 **O(层数·dim + dim²)**，**无 dim³**。

| 级 | 名称 | 数学 | 实现 sketch | 覆盖示例 |
|----|------|------|-------------|----------|
| **S0** | 零元 | \(u=0\) | 已有 | 0 |
| **S1** | 二次层生成元 | 塔链 deg-2：±g, ±k·g | 已有 `try_sqrt_layer_generators` | √2, 2√2, √8 in ℚ(√2) |
| **S2** | 基元平方 | ±eᵢ | 已有 `try_sqrt_basis_squares` | u 为某一基坐标 |
| **S3** | **子域下降** | \(u\in F\subset L\) 且 \(F\) 为真子域 → 在 \(F\) 递归 try，再 embed | `try_subfield_embedding` + 按 parent 链递归 | u 仅含低层 coords |
| **S4** | **二次扩张闭式** | 顶层 \(L=F(\sqrt d)\)，\(u=a+b\omega\)；解 \((a+b\omega)^2=u\) 的 \(a,b\in F\) | 写 `try_sqrt_quadratic_layer(field, u)` | ratio 型、\(F(\sqrt\alpha)\) 内元素 |
| **S5** | **层生成元代数** | 有限集 \(\{g_i\}\)（各 deg-2 层）：试 \(g_i g_j\)、\(g_i/g_j\)、\(c_1 g_i+c_2 g_j\)（\(c\in\{\pm1,\pm2\}\)） | 推广现有 `try_sqrt_small_combo`，**不**扫全基 | resolvent 后组合根 |
| **S6** | ponytail fallback | dim≤6：二元 ±1；**禁止**三元枚举 | 仅 debug/兜底 | 漏网 |

**删除（F5 DoD）：** `TRY_SQRT_TRIPLE_*`、dim>6 的 quadratic combo、dim>12 的 pairwise 枚举（或改为 S6 仅 test-only）。

#### 分期（AFK 可抓取）

| 期 | 内容 | 估时 | 验收 |
|----|------|------|------|
| **F5.0** | API 契约入 `ext_tower` 模块 doc + `algorithm-expr-api`；列共用调用点 grep | 0.5d | 模块 doc ✅；调用表 / algorithm-expr-api ⬜ |
| **F5.1** | **S3 子域下降** | 1d | ✅ `try_sqrt_subfield_descent` + `try_square_root_sqrt2_in_q_sqrt2_sqrt3` |
| **F5.2** | **S4 二次层闭式** | 1–1.5d | ✅ `u=b²·d` + 一般 `(a+bω)²=u` |
| **F5.3** | **S5 层生成元代数**（替代 F4 enum 块） | 1d | ✅ `try_sqrt_layer_generator_algebra` + `try_square_root_sqrt6_in_q_sqrt2_sqrt3` |
| **F5.4** | 删 dim³ / 收紧 S6；`try_sqrt_in_field_shallow` 合并回 full try | 0.5d | ✅ 无 dim³；S6 dim≤6；`poly_roots::tests` ~3s 25/25 绿 |
| **F5.5** | `FieldSession` 薄化；登记 `giac-core-algebra-api-stability` | 0.5d | diff 仅 delegation ⬜ |

**合计 F5：** ~4–5d（与 F4′ 可并行设计，**实现建议 F5.1→F5.3 先于 F4′**）。

#### F4′（并行轨，依赖 S3/S4 或最小 `apply_automorphism`）

| 步 | 内容 |
|----|------|
| F4′.1 | `ExtensionField::automorphisms_fixing_subfield` 或 resolvent 根置换表（ponytail：三根固定 Vieta） |
| F4′.2 | `sqrt_in_field_euler_second`：shallow miss → 对 \(\sigma(\sqrt\alpha)\) 试 \(\sigma(\alpha)=\beta\) |
| F4′.3 | 二次层 σ(κ)；Gal≅A₄ **分测** \(d_L\le 12\)；`t⁴+t+1` baseline \(d_L\le 24\) |

#### 总 DoD（F5 + F4′）

- [x] 热路径无 dim³ 枚举
- [x] `cargo test -p giac-core poly_roots::tests --release` ≤10s，25/25 绿
- [ ] Gal≅A₄ 分测：\(d_L\le 12\)（F4′）；`t⁴+t+1`：\(d_L\le 24\) + 四根 verify ✅
- [x] §1 `TRY_SQRT_*` interim 常量删除或仅 S6 dim≤6（`TRY_SQRT_PAIRWISE_DIM_CEILING=6`）
- [x] [known-divergences.md](../known-divergences.md) 登记「枚举→结构」closed — 见 [GIAC-poly-p3-6-quartic-roots-gaps](GIAC-poly-p3-6-quartic-roots-gaps.md) G1/G5

### 范围

| 做 | 不做 |
|----|------|
| `ExtensionField::try_square_root(coords) -> Option<CoordsQ>` 或 `AlgExtData` 级 API | 阻塞 F3（可选） |
| `FieldSession::adjoin_sqrt` 薄封装调用塔 API | n 次根、任意 `rootof` 合并 |
| 与 F1 测试共享核心断言 ε²=u | 改 T4a common 默认 |

### 任务

1. 从 F1 抽出探测核心到 `ext_tower.rs`（或 `field_arith`）。
2. `#[doc]` 说明：与 `algext_square_roots`（**总是 adjoin**）的分工。
3. 塔层单测：ℚ(√2) 内 √8、√2 等；与 `align_elements` 组合。

### 验收（DoD）

- [x] F1 行为不变（仅重构或薄化）
- [x] `ext_tower` 新测 ≥2 例（`sqrt2_in_q_sqrt2_sqrt3`、`sqrt6_in_q_sqrt2_sqrt3`）
- [ ] [giac-core-algebra-api-stability.md](../giac-core-algebra-api-stability.md) 登记 tier（若存在该文档）

### 文件

- `giac-core/src/algebra/ext_tower.rs`
- `giac-core/src/algebra/field_session.rs`

---

## 2. 测试矩阵（跨 F1–F4）

| 测试名 | 归属 | 目标 | **现状** |
|--------|------|------|----------|
| `adjoin_sqrt_squares_one_resolvent_root` | F1 | L₃ 单根 adjoin | ✅ 绿 |
| `adjoin_sqrt_after_resolvent_deflate_only` | F1 | deflate 不破坏 | ✅ 绿 |
| `adjoin_sqrt_after_quadratic_roots_formula` | F1 | formula 后 adjoin | ✅ 绿 |
| `adjoin_sqrt_after_sqrt_disc_on_resolvent` | F1 | resolvent + √Δ 后 adjoin | ✅ 绿 |
| `resolvent_cubic_all_roots_vanish` | F2 | 三根 + dim≤6 | ✅ `resolvent_cubic_roots_in_session` |
| `f2_resolvent_split_no_cardano_stack` | F2 | resolvent dim≤6 | ✅ 绿 |
| `roots_x3_minus_x_plus_1_vanish` | F2/C1 | 一般三次 | ✅ dim≤6 |
| `euler_gamma_relation_holds` | F3 | \(\sqrt\alpha\sqrt\beta\sqrt\gamma+q=0\) | ✅ 绿 |
| `euler_four_roots_vanish` | F3 | depressed 四根 | ✅ 绿 |
| `roots_quartic_t4_plus_t_plus_1` | F3 | 公开 API | ✅ 绿 |
| `field_session_dimension_bound_quartic` | F4 | 维数 ≤24 | ✅ 绿 |
| `field_session_dimension_bound_quartic_tight` | F4 | 实测 dim=24 | ✅ 绿 |

---

## 3. 优先级与里程碑

### 3.1 优先级总表

| 优先级 | ID | 内容 | 估时 | 解除条件 |
|--------|-----|------|------|----------|
| **P0** | **F1** | `sqrt_in_field` + F1 回归测（含 2 个新测） | 2–4d | — |
| **P1** | **F2** | `resolvent_cubic_roots_in_session`（路径 A′）+ 改 `resolvent_cubic_all_roots_vanish` | 2–3d | F1 绿 |
| **P2** | **F3** | Euler 单路径 + unignore e2e + `euler_gamma_relation_holds` | 1–2d | F1+F2 绿 |
| **P2** | **F4** | Session 契约 + `field_session_dimension_bound_quartic` | 0.5–1d | F1（与 F3 同 PR） |
| **P3** | D3-1 验收 | `solve_quartic_t4_plus_t_plus_1` unignore | 0.5d | F3 绿 |
| **P4** | **F5** | `ext_tower::try_square_root` 抽出重构 | 1–2d | F1 稳定 |
| **P5** | rootof | `general quartic rootof` 删 NotImplemented | follow-up | F3 绿 |

### 3.2 里程碑（可勾选）

- [x] **M1 — F1 绿：** 4 项 F1 DoD + 既有 2 测仍绿
- [x] **M2 — F2 绿：** `resolvent_cubic_all_roots_vanish` 走新 API；`t^4+t+1` resolvent 阶段 `dim≤6`
- [x] **M3 — P3-6 核心：** `roots_quartic_t4_plus_t_plus_1` 绿；release `poly_roots::tests` <5s、0 ignore
- [x] **M4 — 门禁：** F4 维数测 + Session 契约入 `field_session` doc
- [ ] **M5 — 可选：** F5 ✅（S4 全量）；F4′ Galois **open**（σ(κ) 二次层；`t⁴+t+1` baseline dim=24，非未达 12 之失败）

### 3.3 AFK 抓取顺序

| 顺序 | Issue | 估时 | 备注 |
|------|-------|------|------|
| 1 | **F1** | 2–4d | 先实现 1a+1b；`adjoin_sqrt_after_quadratic_roots_formula` 作红灯测 |
| 2 | **F2** | 2–3d | 与 F1 可同 PR，**宜分 commit**；先 casus/ω 再 Cardano 分支 |
| 3 | **F3** | 1–2d | 删 `permutations`/`approx_real_sign`；验收 release <5s |
| 4 | **F4** | 0.5–1d | 可与 F3 同 PR |
| 5 | **F5** | 1–2d | 可选 |

**禁止：** 在 F1 未绿前 unignore e2e 并加大 `approx_real_sign` 搜索。

**建议首 PR 范围：** F1（完整）+ F1 红灯测；不夹 F3 unignore。

---

## 4. 索引

| 主题 | 文档 |
|------|------|
| 工程实施 F1–F5 | [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) |
| **Lean 4 验证（未来）** | [GIAC-poly-quartic-lean4-verification](GIAC-poly-quartic-lean4-verification.md) |
| 塔数学 | [giac-tower-common-math.md](../giac-tower-common-math.md) |
