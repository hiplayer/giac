# GIAC — F5 FGLM postmortem 整改优先级

**状态:** P1–P5 全部落地 / 已解决（2026-06-29，P6 转出独立 issue，验收全绿，移入 `issues_resolved`）
**类型:** 整改计划 / 优先级排序
**来源:** [GIAC-poly-f5-fglm-debug-postmortem](../issues/GIAC-poly-f5-fglm-debug-postmortem.md)（根因分析 + 解决方案草案）
**相关:** [GIAC-dense-poly1-refactor](../issues/GIAC-dense-poly1-refactor.md)（P2 落点）、[GIAC-poly-f5-fglm-over-coefficient-field](GIAC-poly-f5-fglm-over-coefficient-field.md)（P3 主体，已 resolved）、[GIAC-poly-f5-fglm-fullchain-fuel-audit](../issues/GIAC-poly-f5-fglm-fullchain-fuel-audit.md)（P6 独立 issue，open）
**快照:** 2026-06-29

---

## 1. 问题清单（来自 postmortem §2）

| 根因 | 占调试时间 | 已爆发错误 | 解决方案（postmortem §） |
|---|---|---|---|
| **A. 约定不可见（intra-world）** | ~60% | #5, #7（踩两次） | §3.1 双边 order newtype + 过渡 guardrail |
| **B. 隐式双义参数** | ~15% | #7 norm 过滤放错位置 | §3.2 `FmoduleMode` enum |
| **C. 搜索函数无代价上限** | ~20% | #6 100s 超时 | §3.3 `sqrt_base_case` 接 `Fuel`（局部） |
| **D. 跨算术世界无桥接（inter-world）** | 0%（本轮未触发） | — | dense-poly1-refactor §4.2 转换 API（D4，另案） |
| **E. 机械 Rust 摩擦** | ~5% | #3, #4 | 不改（正常日常） |

> **D 不在本 issue 处理**：dense-poly1-refactor D4 是其根治防线，已在该计划排期。本 issue 只管 A/B/C 三类已爆发根因 + 两个 follow-up（F1 norm 复活、F2 全链 fuel 审计）。

---

## 2. 优先级表（按 ROI = 影响 / 工作量）

| 序 | 项 | 根因 | 影响 | 工作量 | 依赖 | ROI |
|---|---|---|---|---|---|---|
| **P1** ✅ | §3.1 第 0 步 过渡 guardrail（`// order:` 注释 + `lint-coordsq-order.sh` ratchet） | A | 封住**新增**约定错配（不破签名，lint 强制新 fn 必须声明 order） | 小（3 注释 + 1 lint 脚本 + 1 allowlist 88 条） | 无 | ★★★★★ 已落地 |
| **P2** ✅ | §3.1 第 1 步 双边 `LowFirstQ`/`HighFirstQ` newtype（const-generic `FieldCoords<ORD>`） | A（根治） | 消灭 #5/#7 复发，编译期抓所有 HighFirst↔LowFirst 错配（60% 调试时间） | 实际 ~150 站点（`element_*`/`align_pair`/`try_sqrt_*`/`sqrt_fmodule`/F-module 管线全迁），见 postmortem §2.5 inventory | P1 先行 | ★★★★ 已落地 |
| **P3** ✅ | §3.3 局部 `Fuel`（`sqrt_base_case` 接 `Fuel::new(N)`） | C | 封住 #6 已爆发超时，`LIFT_BUDGET=8` const 升级为参数 | 小（1 函数签名 + 调用点） | 无 | ★★★★ 已落地 |
| **P4** ✅ | §3.2 `FmoduleMode` enum（`sqrt_fmodule` + `sqrt_base_case` 签名重构） | B | 消灭 norm 过滤放错位置，`m_gen_low` 冗余不一致 | 小中（2 函数签名 + norm 过滤归位） | 无 | ★★★ 与 P4 milestone 顺手 |
| **P5** ✅ | F1 norm 预过滤复活 | —（纯 perf） | 早退非平方 u，省 `LIFT_BUDGET` 次 lift | 极小（恢复 #7 删除的过滤） | **P2**（约定错配消除后才安全） | ★★★ §3.1 落地后免费收益 |
| **P6** | F2 全链 fuel 审计（60 素数外循环 / Newton 步数 / `factor_mod_irreducibles`） | C（残余） | 封住 `sqrt_fmodule` 递归链上所有无界搜索 | 中大（改递归结构） | P3（P3 是其局部子集） | ★★ 建议开独立 issue |

---

## 3. 执行顺序与阶段

```text
Phase 0（立即，独立，封已爆发）
  ├─ P1  过渡 guardrail ✅      ← 已落地：lint + // order: 注释 + allowlist ratchet
  └─ P3  sqrt_base_case Fuel ✅  ← 已落地：release 6 quartic 测 0.30s 绿

Phase 1（dense-poly1-refactor 首切片，根治 A）
  └─ P2  双边 newtype          ← 影响面见 postmortem §2.5 inventory
                                  DoD 见 postmortem §3.1 末（复用 diag_fmodule_sqrt_recovery / quartic_a4_galois_dim_le_12）

Phase 2（P2 落地后，免费 perf）
  └─ P5  norm 预过滤复活 ✅      ← top-gated：仅 `m_gen_low.is_none()`（顶层）算 N_{K/ℚ}(u) 早退

顺手（与 P4 milestone 并行，改动局部）
  └─ P4  FmoduleMode enum ✅     ← sqrt_fmodule(mode: FmoduleMode) + sqrt_base_case 从 field 取 m_gen

独立 issue（超出 postmortem 修复范畴）
  └─ P6  全链 fuel 审计          ← 改 sqrt_fmodule 递归结构，单独排期
```

**关键路径：** P1 → P2 → P5（P2 是 P5 的硬依赖；P5 在 P2 未落地前**禁止**启用，否则 #7 误杀复发）。
**可并行：** P3、P4 与 Phase 0/1 互不依赖，可并行做。
**P6 独立：** 不阻塞本 issue 闭合，单独开 issue 跟踪。

---

## 4. 各项验收（DoD）

| 项 | DoD |
|---|---|
| **P1** ✅ | `scripts/lint-coordsq-order.sh` 存在并接入 `ci-clippy.sh` + `cargo ci-clippy` 绿；`field.element_mul` / `poly_from_low` / `krylov_minpoly_coords` 标注 `// order:` 注释；88 个未标注 fn 进 `scripts/coordsq-order-allowlist` ratchet。`debug_assert` 不可行（`CoordsQ` 无运行期 order 值），并入 P2 newtype。 |
| **P2** ✅ | postmortem §3.1 末 DoD：`diag_fmodule_sqrt_recovery` + `quartic_a4_galois_dim_le_12` 仍 green（`cargo test -p giac-core` 305 passed / 0 failed）**且** 反向误用编译失败已验证 —— `sqrt_fmodule` 的 `f_basis: Vec<HighFirstQ>` 拒收 raw `CoordsQ`（`poly_roots.rs:1933: error[E0308]` 探针确认）。newtype 用 const-generic `FieldCoords<{HIGH_FIRST}>`/`FieldCoords<{LOW_FIRST}>` 统一定义，`Deref→Vec<Ratio>` 保留借用侧零摩擦，owned 边界显式 `.0` 标记。 |
| **P3** ✅ | `sqrt_base_case` 签名含 `fuel: &Fuel`（crate `Fuel` 用 `&Fuel` 不可变引用 + `Cell` 内部可变，非 `&mut`）；`const LIFT_BUDGET` 删除，`8` 移至调用方 `Fuel::new(8)`；release 6 quartic 测 0.30s 绿（≤0.4s/测 DoD 满足），`diag_fmodule_sqrt_recovery` 绿 |
| **P4** ✅ | `sqrt_fmodule(field, u, mode: FmoduleMode)`，`enum FmoduleMode { Top, Recurse }`（无数据载体，因 `sqrt_base_case` 改从 `field.generator_minpoly_low()` 取 m_gen，mode 纯控制流判别）；`sqrt_base_case` 去掉 `m_gen_low` 参，新增 `ExtensionField::generator_minpoly_low() -> Option<CoordsQ>`（high-first `min_poly_q` 反转成 low-first，仅单层 base extension / `common_over_q` 命中，tower/block 返回 None bail）；norm 过滤（P5）挂 `FmoduleMode::Top` 分支（`matches!(mode, FmoduleMode::Top)`）；top 基例 bail（外部构建的 K 可能 u≠存储生成元，m_gen 不可安全取），仅 `Recurse` 进 base case 读 field。全 quartic 测绿（`cargo test -p giac-core --lib` 306 passed / 0 failed）。 |
| **P5** ✅ | `N_{K/ℚ}(u)=(−1)^{d_k}·c0^{d_k/d_f}` 早退过滤恢复（`c0`=low-first monic minpoly 常数项，`sqrt_fmodule` 内 `is_q_square` 谓词）。P2-followup 落地后**全层启用**（原 top-gated 限制解除）：递归边界改传**真 HighFirst** 字节（`HighFirstQ::from_low`，反转），krylov 的 `element_mul` 在每层都正确读 u′ → norm 在递归层也正确 → 过滤对 Top 与 Recurse 都安全。base case 端到端 `LowFirstQ`（`sqrt_base_case(field, u: &LowFirstQ, fuel) -> Option<LowFirstQ>`，`generator_minpoly_low() -> Option<LowFirstQ>`）。附带修了 d_k==1 分支 `BigInt::sqrt` 对负数 panic 的潜在 bug（mislabel 旧路径绕开，修正路径触发后加负数守卫）。DoD：A₄ 真平方仍命中（`quartic_a4_galois_dim_le_12` + `diag_fmodule_sqrt_recovery` 绿）+ biquadratic 两测绿（`x⁴+1`/`t⁴−2`）+ `is_q_square_predicate` 单元测 + `cargo test -p giac-core --lib` 306 passed / 0 failed。perf 收益：非平方 u 在任意递归深度跳过 `sqrt_base_case` 的 Fuel 搜索。 |
| **P6** | 独立 issue DoD（待立项）：`sqrt_fmodule` 递归链所有循环（素数外循环、Newton 步数、`factor_mod_irreducibles`）均经 `Fuel` 受控 |

---

## 5. 关闭条件（revisit trigger）

本 issue 在以下**全部**满足后关闭（与 postmortem §0.2 一致）：

- [x] P1 落地 → 过渡 guardrail（lint + `// order:` 注释 + allowlist ratchet）就位，封住新增约定错配
- [x] P1 + P2 落地 → 关闭根因 A
- [x] P3 落地 → 关闭根因 C（局部）
- [x] P4 落地 → 关闭根因 B
- [x] P5 落地 → 收割 perf 收益（全层启用，P2-followup 解除 top-gated 限制）
- [x] P6 立独立 issue → 根因 C 残余转出（[GIAC-poly-f5-fglm-fullchain-fuel-audit](../issues/GIAC-poly-f5-fglm-fullchain-fuel-audit.md)，postmortem §C 残余待该 issue 闭合）

**根因 D** 不在本 issue 关闭条件内（归 dense-poly1-refactor D4）。

---

## 6. 风险

| 风险 | 缓解 |
|---|---|
| P2 newtype churn 大、跨 crate 签名变更 | P1 guardrail 先行封新增；P2 作为 dense-poly1-refactor 首切片单独 PR，编译器驱动小步验证 |
| P5 在 P2 未落地前误启用 → #7 复发 | P5 硬依赖 P2，PR 顺序强制；P5 启用前必须跑 `quartic_a4_galois_dim_le_12` |
| P4 `FmoduleMode` 与 P2 newtype 同时改 `sqrt_fmodule` 签名 → 冲突 | 错峰：P4 在 Phase 1 之后或同 PR；两者改的是不同参数（mode vs `CoordsQ`→`LowFirstQ`），可合并不冲突 |
| P6 全链 fuel 改递归结构、回归面大 | 单独 issue + 独立 PR，不阻塞本 issue |
