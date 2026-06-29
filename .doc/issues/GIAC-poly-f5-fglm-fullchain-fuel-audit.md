# GIAC — F5 FGLM / F-module √-reduction 全链 fuel 审计

**状态:** open
**类型:** 实现 / 性能护栏（AFK 可抓取）
**来源:** [GIAC-poly-f5-fglm-debug-postmortem](GIAC-poly-f5-fglm-debug-postmortem.md) §3.3 follow-up **F2**、[GIAC-poly-f5-fglm-postmortem-remediation](../issues_resolved/GIAC-poly-f5-fglm-postmortem-remediation.md) **P6**
**根因:** C（搜索函数无代价上限）残余 —— postmortem §3.3 的「局部」只封了 `sqrt_base_case`，本 issue 处理「全链」
**Rust 落点:** `giac-core::algebra::poly_roots`（`sqrt_fmodule` / `sqrt_base_case` / `padic_sqrt_lift`）
**相关:** [GIAC-poly-f5-fglm-over-coefficient-field](../issues_resolved/GIAC-poly-f5-fglm-over-coefficient-field.md)（P3 主体，已 resolved）、crate `Fuel` 抽象（`EvalError` 管线、S7 `Fuel::new(8)`）
**快照:** 2026-06-29

---

## 0. 背景

postmortem §3.3 把 #6（100s 超时）的修复分为**局部**与**全链**：

- **局部（P3，✅ 已落地）**：`sqrt_base_case` 接 `Fuel::new(8)`，封住已爆发的 `(prime, combo)` p-adic 搜索。release 6 quartic 测 0.30s 绿。
- **全链（P6，本 issue）**：`sqrt_fmodule` 递归链上**还有三个无界搜索**，未受 `Fuel` 控制。改递归结构超出了「postmortem 修复」范畴，故独立排期。

本 issue 是 postmortem §C 残余的关闭条件（见 postmortem §0.2 revisit trigger 最末一行：全链 fuel 审计独立 issue **闭合** → 关闭 §C 残余）。

---

## 1. 范围：三个无界循环

均在 `giac-core/src/algebra/poly_roots.rs` 的 `sqrt_fmodule` 递归链上：

| # | 位点 | 现状 | 无界性 |
|---|------|------|--------|
| **1** | `sqrt_base_case` `let primes = small_primes(60);`（`:1488`）+ `for p in primes` 外循环 | magic number `60` 硬编码；非平方 u 会穷举全部 60 素数 | 素数个数无参数化、无 fuel 上限；`small_primes(n)`（`:1814`）本身可造任意多 |
| **2** | `padic_sqrt_lift` `while pk.bits() <= target_bits`（`:1610`），`target_bits = 320`（`:1489`） | Newton 倍步 `p^{2j}`，循环到 `p^k > 2^320` | `k`（步数）无硬上限 —— 对小素数 `p`，`k` 可很大；`target_bits` 是 magic number |
| **3** | `sqrt_base_case` `factor_mod_irreducibles(&m_gen_poly, p)`（`:1501`） | `giac-poly` 内部对 `m_gen mod p` 做无界因子搜索 | `factor_mod_irreducibles` 内部循环未经调用方 fuel 传递 |

> **注**：位点 1 的 `(prime, combo)` fan-out 已被 P3 的 `Fuel::new(8)` **部分**封住（`fuel.enter()` 在 combo 内循环 `:1557` 处计次）。但**素数外循环本身**（`for p in primes`）仍跑满 60，且 `fuel` 只在 combo 层计次、不在素数层计次 —— 一个素数若 yield 0 viable combo，外循环继续下一个素数，无总预算封顶。

---

## 2. DoD

- `sqrt_fmodule` 递归链上三个循环均经 `Fuel` 受控：
  - 位点 1：`small_primes(60)` 的 `60` 参数化（调用方传入或从 `Fuel` 派生）；素数外循环计入总 fuel 预算。
  - 位点 2：Newton lift 步数 `k` 有硬上限（`Fuel` 或显式 `max_lift_steps`）；`target_bits` 参数化。
  - 位点 3：`factor_mod_irreducibles` 接受调用方 fuel（或其内部搜索有界且该界经调用方控制）。
- 超限时走 `None` / `Err` 而非静默穷举（与 P3 一致：`miss ≠ wrong answer`，caller 回退到 blind adjoin）。
- release quartic 测仍绿且不超时（≤0.4s/测，沿用 P3 DoD）；`diag_fmodule_sqrt_recovery` + `quartic_a4_galois_dim_le_12` 绿。
- 新增一个「非平方输入超限早退」回归测（构造一个会让旧路径穷举到超时的非平方 u，断言新路径在 fuel 耗尽时快速 `None`）。

---

## 3. 依赖与风险

- **依赖**：P3（✅ 已落地，是本 issue 的局部子集 —— `Fuel` 抽象与 `sqrt_base_case` 的 `fuel: &Fuel` 签名已就位，本 issue 把 fuel 往上游递归链传）。
- **风险**：改 `sqrt_fmodule` 递归结构、`padic_sqrt_lift` 签名、可能触及 `giac-poly::factor_mod_irreducibles` 公开 API（跨 crate）→ 回归面大。**单独 PR**，不阻塞 postmortem/remediation issue 闭合。
- **关键路径无关**：P6 不阻塞任何已落地项（P1–P5）；本 issue 闭合仅关 postmortem §C 残余。

---

## 4. 执行建议

1. 先做位点 1（最小、最局部）：`60` 参数化 + 素数外循环计 fuel，封住「非平方 u 跑满 60 素数」。
2. 位点 2：`padic_sqrt_lift` 加 `max_lift_steps`（或从 `Fuel` 派生 `target_bits → max_steps`）。
3. 位点 3：评估 `factor_mod_irreducibles` 是否需要跨 crate API 改（如加 `&Fuel` 参数），或在调用方包一层超时探测。
4. 每步跑 release quartic 测 + 新增回归测，编译器/测驱动。

---

## 5. 关闭条件

- [ ] 位点 1/2/3 均经 `Fuel` 受控
- [ ] release quartic 测绿且不超时
- [ ] 非平方输入超限早退回归测落地
- [ ] 本 issue 闭合 → 关闭 postmortem §C 残余（postmortem §0.2 最末一行）→ postmortem 可移入 `issues_resolved`
