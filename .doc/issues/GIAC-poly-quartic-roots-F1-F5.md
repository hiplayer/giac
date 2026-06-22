# GIAC-poly — 四次 `roots` 塔修复（F1→F5）

**状态:** open  
**类型:** 实现 / AFK 可抓取  
**父项:** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) **P3-6**（通用四次 `Poly<AlgExtC>::roots`）  
**相关:** [GIAC-algext-adoption](GIAC-algext-adoption.md) §8.2–8.4、[giac-tower-common-math.md](../giac-tower-common-math.md)、[known-divergences.md](../known-divergences.md) DIV-076–083 / DIV-081  
**Rust 落点:** `giac-core::algebra::{field_session, poly_roots, ext_tower}`  
**快照:** 2026-06-20

---

## 0. 问题陈述

`poly_roots::quartic_roots`（Euler + Ferrari resolvent）在 **resolvent 分裂后** 再 `adjoin_sqrt` 时破坏 \(\varepsilon^2=u\)（塔嵌入与已有平方根未识别）。表现：

- `adjoin_sqrt_squares_one_resolvent_root`、`adjoin_sqrt_after_resolvent_deflate_only` ✅
- `sqrt_disc` / `quadratic_roots_formula` / 全 `cubic_roots` 后再 adjoin ❌
- `euler_four_roots_vanish`、`roots_quartic_t4_plus_t_plus_1` **`#[ignore]`**

**根因（摘要）：** 盲目 `adjoin_irreducible(u²−α)` 叠塔，未在 **当前 L** 内查找已有平方根；resolvent 二次分裂又额外扩域，维数爆炸 + 生成元语义错位。

**目标：** `solve(t^4+t+1=0,t)` 四根 `eq_mod` 零化；release 下 `poly_roots::tests` 无 ignore、无分钟级穷举。

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

## F1 — `sqrt_in_field` / `try_existing_square_root`

### 问题

`FieldSession::adjoin_sqrt` 总是 adjoin 新元；当 \(u\) 已在 **L** 中为某元平方时，应返回 **嵌入 L 的** \(\varepsilon\) 且保证 \(\varepsilon^2 \equiv u\)，而非新层 \(\varepsilon'\) 满足 \((\varepsilon')^2\) 在错误子域/mod 下不等于 \(u\)。

### 范围

| 做 | 不做 |
|----|------|
| `try_sqrt_in_field(u) -> Option<AlgExtCPolyCoeff>`（或 `sqrt_in_field` 先 try 再 adjoin） | 一般 \(n\) 次根（仅平方） |
| `adjoin_sqrt` 改为：`try` → 成功则 `lift` + 验 \(\varepsilon^2=u\)；失败再 adjoin | 改 giac C++ `common_EXT` 行级行为 |
| `sqrt_principal` / `sqrt_disc` / `pick_sqrt_euler` 统一走该入口 | F2 resolvent 三根算法（留给 F2） |

### 任务

1. 在 `ext_tower` 或 `field_session` 实现 **L 内平方根探测**（候选：对 `dim(L)≤N` 枚举/矩阵法；或沿塔层检查 generator 幂次；文档写明算法与复杂度上界）。
2. `FieldSession::adjoin_sqrt`：`try_sqrt_in_field` → 验 `mul(ε,ε) eq_mod u` → 否则 `adjoin_irreducible(u²−α)`。
3. `sqrt_disc`：去掉「adjoin 后验失败再盲 adjoin \(-\Delta\)」的脆弱双路径，优先 `sqrt_in_field`。
4. 扩展回归：**在** `quadratic_roots_formula` / `sqrt_disc` **之后** adjoin 仍满足 ε²=u（当前失败场景）。

### 验收（DoD）

- [ ] 新测 `adjoin_sqrt_after_quadratic_roots_formula`（或等价名）通过
- [ ] 新测 `adjoin_sqrt_after_sqrt_disc_on_resolvent`（resolvent 一条根 + formula 二次后 adjoin）通过
- [ ] 既有 `adjoin_sqrt_squares_one_resolvent_root`、`adjoin_sqrt_after_resolvent_deflate_only` 仍绿
- [ ] `cargo test -p giac-core poly_roots::tests --release` 无新增 ignore

### 文件

- `giac-core/src/algebra/field_session.rs`
- `giac-core/src/algebra/ext_tower.rs`（若探测在塔层）
- `giac-core/src/algebra/poly_roots.rs`（`sqrt_disc`、测试）

### 偏离

- 闭合 [DIV-081](../known-divergences.md)「实现缺口」子项；[DIV-077](../known-divergences.md) 复审

---

## F2 — resolvent 三根留在 L₃（无二次分裂叠塔）

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
2. 实现路径（择一或组合，PR 说明）：
   - **A.** `z0 = one_cubic_root`；`z1,z2` 由 \(z^2-(z_0+z_1)z+z_0 z_1\) 系数在 L₃ 内用 `sqrt_in_field` + 有理运算；
   - **B.** 三次因子在 L₃ 上 split 为线性因子（域内根公式，无新 adjoin）；
   - **C.** 登记于 [known-divergences](../known-divergences.md) 的若与 giac 不同的 **数学等价** 分支。
3. 删除 `quartic_roots` 对 `quadratic_roots_formula` 的 resolvent 依赖；`cubic_roots` 不作为四次子调用。
4. 单测：`resolvent_cubic_all_roots_vanish` 在 **不调用** `quadratic_roots_formula` 下仍绿；三根 `eq_mod` 互异且零化 \(R\)。

### 验收（DoD）

- [ ] `resolvent_cubic_all_roots_vanish` 走新路径
- [ ] `session.working().dimension()` 在 `t^4+t+1` resolvent 阶段 ≤ F4 约定上界（暂定 **≤12** 或 PR 中论证的 tight bound）
- [ ] F1 回归全套仍绿

### 文件

- `giac-core/src/algebra/poly_roots.rs`（`one_cubic_root`、`deflate_monic`、`quartic_roots`）
- `giac-core/src/algebra/field_session.rs`（必要时 `adjoin_primitive_cube_root_of_unity` 契约）

### 偏离

- [DIV-080](../known-divergences.md)、[DIV-081](../known-divergences.md)

---

## F3 — Euler 固定分支、去穷举、e2e 绿

### 问题

`euler_depressed_quartic_roots`：**6 置换 × flip 搜索** + `pick_sqrt_euler` 数值 `approx_real_sign`；单次失败 quartic 运行 ~15–35s。忽略测 `euler_four_roots_vanish`、`roots_quartic_t4_plus_t_plus_1`。

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
2. 测试 `field_session_dimension_bound_quartic`：`t^4+t+1` 从 ℚ 到四根输出，`working.dimension()` ≤ **B**（与 F2 PR 共同敲定 B，建议 12–24 有文档理由）。
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

## F5 — `ext_tower` 通用开方 API（可选）

### 问题

F1 若在 `field_session` 内实现探测，逻辑可能重复；长期宜在 `ExtensionField` 暴露 **「在域内求平方根」** 供 `poly_roots`、partfrac、未来 `gcd` 共用。

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

- [ ] F1 行为不变（仅重构或薄化）
- [ ] `ext_tower` 新测 ≥2 例
- [ ] [giac-core-algebra-api-stability.md](../giac-core-algebra-api-stability.md) 登记 tier（若存在该文档）

### 文件

- `giac-core/src/algebra/ext_tower.rs`
- `giac-core/src/algebra/field_session.rs`

---

## 2. 测试矩阵（跨 F1–F4）

| 测试名 | F1 | F2 | F3 | 说明 |
|--------|----|----|-----|------|
| `adjoin_sqrt_squares_one_resolvent_root` | ✅ | — | — | L₃ 单根 adjoin |
| `adjoin_sqrt_after_resolvent_deflate_only` | ✅ | — | — | deflate 不破坏 |
| `adjoin_sqrt_after_quadratic_roots_formula` | ✅ 新增 | — | — | formula 后 adjoin |
| `resolvent_cubic_all_roots_vanish` | — | ✅ | — | 三根无 formula adjoin |
| `euler_four_roots_vanish` | 依赖 | 依赖 | ✅ unignore | depressed 四根 |
| `roots_quartic_t4_plus_t_plus_1` | 依赖 | 依赖 | ✅ unignore | 公开 API |
| `field_session_dimension_bound_quartic` | — | ✅ | ✅ | F4 |

---

## 3. AFK 抓取建议

| 顺序 | Issue | 估时 | 备注 |
|------|-------|------|------|
| 1 | **F1** | 2–4d | 阻塞一切；先读 `adjoin_sqrt_after_*` 与 `ext_tower::adjoin_irreducible` |
| 2 | **F2** | 2–3d | 与 F1 同 PR 可，但宜分 commit |
| 3 | **F3** | 1–2d | 删代码为主；验收 release 耗时 |
| 4 | **F4** | 0.5–1d | 可与 F3 同 PR |
| 5 | **F5** | 1–2d | 可选；F1 稳定后重构 |

**禁止：** 在 F1 未绿前 unignore e2e 并加大 `approx_real_sign` 搜索（已证无效且慢）。

---

## 4. 索引

| 主题 | 文档 |
|------|------|
| 工程实施 F1–F5 | [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) |
| **Lean 4 验证（未来）** | [GIAC-poly-quartic-lean4-verification](GIAC-poly-quartic-lean4-verification.md) |
| 塔数学 | [giac-tower-common-math.md](../giac-tower-common-math.md) |
