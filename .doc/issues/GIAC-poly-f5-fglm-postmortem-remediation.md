# GIAC — F5 FGLM postmortem 整改优先级

**状态:** open
**类型:** 整改计划 / 优先级排序
**来源:** [GIAC-poly-f5-fglm-debug-postmortem](GIAC-poly-f5-fglm-debug-postmortem.md)（根因分析 + 解决方案草案）
**相关:** [GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md)（P2 落点）、[GIAC-poly-f5-fglm-over-coefficient-field](../issues_resolved/GIAC-poly-f5-fglm-over-coefficient-field.md)（P3 主体，已 resolved）
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
| **P4** | §3.2 `FmoduleMode` enum（`sqrt_fmodule` + `sqrt_base_case` 签名重构） | B | 消灭 norm 过滤放错位置，`m_gen_low` 冗余不一致 | 小中（2 函数签名 + norm 过滤归位） | 无 | ★★★ 与 P4 milestone 顺手 |
| **P5** | F1 norm 预过滤复活 | —（纯 perf） | 早退非平方 u，省 `LIFT_BUDGET` 次 lift | 极小（恢复 #7 删除的过滤） | **P2**（约定错配消除后才安全） | ★★★ §3.1 落地后免费收益 |
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
  └─ P5  norm 预过滤复活        ← 仅在 P2 验证 #5/#7 编译失败后启用

顺手（与 P4 milestone 并行，改动局部）
  └─ P4  FmoduleMode enum       ← 都在 sqrt_fmodule/sqrt_base_case

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
| **P4** | `sqrt_fmodule(field, u, mode: FmoduleMode)`；`sqrt_base_case` 不再收 `m_gen_low`（从 `field` 取）；norm 过滤挂在 `Mode::Top` 分支；全 quartic 测绿 |
| **P5** | `N_{K/ℚ}(u)` 早退过滤恢复；非平方 u 提前退出（diag 计数 lift 次数下降）；A₄ 真平方仍命中（`quartic_a4_galois_dim_le_12` 绿） |
| **P6** | 独立 issue DoD（待立项）：`sqrt_fmodule` 递归链所有循环（素数外循环、Newton 步数、`factor_mod_irreducibles`）均经 `Fuel` 受控 |

---

## 5. 关闭条件（revisit trigger）

本 issue 在以下**全部**满足后关闭（与 postmortem §0.2 一致）：

- [x] P1 落地 → 过渡 guardrail（lint + `// order:` 注释 + allowlist ratchet）就位，封住新增约定错配
- [x] P1 + P2 落地 → 关闭根因 A
- [ ] P3 落地 → 关闭根因 C（局部）
- [x] P3 落地 → 关闭根因 C（局部）
- [ ] P4 落地 → 关闭根因 B
- [ ] P5 落地 → 收割 perf 收益
- [ ] P6 立独立 issue → 根因 C 残余转出

**根因 D** 不在本 issue 关闭条件内（归 dense-poly1-refactor D4）。

---

## 6. 风险

| 风险 | 缓解 |
|---|---|
| P2 newtype churn 大、跨 crate 签名变更 | P1 guardrail 先行封新增；P2 作为 dense-poly1-refactor 首切片单独 PR，编译器驱动小步验证 |
| P5 在 P2 未落地前误启用 → #7 复发 | P5 硬依赖 P2，PR 顺序强制；P5 启用前必须跑 `quartic_a4_galois_dim_le_12` |
| P4 `FmoduleMode` 与 P2 newtype 同时改 `sqrt_fmodule` 签名 → 冲突 | 错峰：P4 在 Phase 1 之后或同 PR；两者改的是不同参数（mode vs `CoordsQ`→`LowFirstQ`），可合并不冲突 |
| P6 全链 fuel 改递归结构、回归面大 | 单独 issue + 独立 PR，不阻塞本 issue |
