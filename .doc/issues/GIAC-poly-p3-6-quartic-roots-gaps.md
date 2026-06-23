# GIAC-poly — P3-6 `poly_algext_roots` 缺口索引

**状态:** open（P3-6 **DoD ✅**；F4′ / F5 plan 仍跟踪于 [F1-F5](GIAC-poly-quartic-roots-F1-F5.md)）  
**类型:** 缺口 / AFK 索引  
**父项:** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) **P3-6**  
**算法规格（已归档）:** [issues_resolved/GIAC-poly-p3-6-roots-algorithm-spec.md](../issues_resolved/GIAC-poly-p3-6-roots-algorithm-spec.md)  
**plan 续篇:** [GIAC-poly-p3-6-roots-algorithm-spec.md](GIAC-poly-p3-6-roots-algorithm-spec.md) §12  
**Rust 落点:** `giac-core::algebra::{poly_roots, field_session, ext_tower}`  
**快照:** 2026-06-23（C1–C4、F1–F3、S2/S3 落地后复审）

---

## 0. P3-6 范围

normative 能力（backlog §5.1）：

```text
poly_algext_roots / poly_algext_roots_for_ctx
  deg 1 → linear
  deg 2 → quadratic + sqrt(Δ)
  deg 3 → split_monic_cubic_roots_in_session（§3 算法规格）
  deg 4 → resolvent cubic + Euler（四根 eq_mod 零化）
  deg ≥ 5 → NotImplemented（solve S0 用 rootof 一支，非 P3-6）
```

**非 P3-6：** `giac-solve` 的 factor/sqff 降次（[GIAC-solve-polynomial-pipeline-plan](GIAC-solve-polynomial-pipeline-plan.md) S0–S4）。

---

## 1. 症状总表（2026-06-23 复审）

| 输入 | `poly_roots` `verify_root` | `solve` S0 | 阻塞 issue |
|------|---------------------------|------------|------------|
| `x²+1=0` | ✅ | ✅ | — |
| `x³−2=0` | ✅ | — | — |
| **`x³−x+1=0`** | ✅ `roots_x3_minus_x_plus_1_vanish` | ✅ 代入 | — |
| resolvent `z³−4z−1=0` | ✅ | — | — |
| **`t⁴+t+1=0`** | ✅ `roots_quartic_t4_plus_t_plus_1` | ✅ `solve_quartic_t4_plus_t_plus_1` | — |
| `t⁴−2=0` | ✅ | ✅ | — |
| `x⁵−x+1=0` | N/A | ✅ 单支 `rootof` | — |
| `(x²+1)(x³−x+1)=0` | ✅ | ✅ 5 根代入 | — |

**结论：** P3-6 **核心代数根**与 solve **S0/S2 代入**已绿；剩余主要是 **F4′ Galois（避免多余 adjoin）** 与 solve **S4–S6** 下游。

---

## 2. 四次缺口 G0–G5（→ F1–F5）

| ID | 标题 | 状态 | 验收 |
|----|------|------|------|
| **G0** | `FieldSession` + `poly_algext_roots_for_ctx` | ✅ | field-session-plan §9 |
| **G1** | `sqrt_in_field` / `try_sqrt_in_field` | ✅ | F1 回归测绿 |
| **G2** | resolvent 三根（纯三次 ω；p≠0 deflate+F1） | ✅ | `f2_resolvent_split_no_cardano_stack` |
| **G3** | Euler 单路径；`t⁴+t+1` 四根 | ✅ | `euler_four_roots_vanish` 无 ignore |
| **G4** | Session 契约 + 维数上界 | ✅ | 硬顶 **≤24** ✅；`t⁴+t+1` 实测 **dim=24**（baseline，非 Gal 未证之「12 目标」失败） |
| **G5** | 热路径无裸 `align_coeff` | ✅ | algorithm-expr-api §6.3 |

**详细任务：** [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) §F1–F5。

---

## 3. 三次缺口 C1–C4（已关闭，2026-06-23）

| ID | 标题 | 状态 | 证据 |
|----|------|------|------|
| **C1** | `x³−x+1` verify_root | ✅ | `roots_x3_minus_x_plus_1_vanish` |
| **C2** | 统一 `split_monic_cubic_roots_in_session` | ✅ | 算法规格 §3 |
| **C3** | 三次维数 bound | ✅ | `field_session_dimension_bound_cubic_x3_minus_x_plus_1` |
| **C4** | solve 代入 | ✅ | `solve_poly` 7/7 绿 |

---

## 4. Issue 依赖图（剩余）

```text
F4 closed (t⁴+t+1 dim≤24 baseline ✅)
  └─ field_session_dimension_bound_quartic_tight ✅

F4′ open (Galois σ(κ); A₄-only d_L≤12 须 Gal 已证)

solve 下游:
  S4 删 rootof 旁路 ✅（生产走 S0；rootof.rs 仅测试）
  S6 realroot ✅（Sturm 计数 + casus 三次过滤）
  S7 fsolve
  general quartic rootof (rootof.rs)
```

**AFK 建议顺序（plan）：** **F5** 结构开方 → **F4′** Galois σ(κ) → S7 → 文档同步。**不要** 把未证 Gal 的 `d_L=12` 当全体四次门禁。

---

## 5. 与 solve 管线的关系

| solve 计划 | 状态 |
|------------|------|
| S0 共享内核 | ✅ |
| S1 `eval_solve` → S0 | ✅ |
| S2 删 stale ignore | ✅（`solve_quartic_t4_plus_t_plus_1`） |
| S3 `froot` deg≤4 | ✅（接 `solve_irreducible_factor`） |
| S4 删 `rootof.rs` fallback | ✅（`solve_poly` S0；`rootof.rs` 测试/参考） |
| S6 `realroot` | ✅（Sturm 计数 + casus 三次；全 Sturm 隔离待扩展） |
| S7 `fsolve` | 🔴 |

---

## 6. 验收矩阵（P3-6 总 DoD）

| 检查 | 测例 | 状态 |
|------|------|------|
| 二次 / 纯三次 / 一般三次 | 见 §1 | ✅ |
| resolvent / 双二次 / 一般四次 | 见 §1 | ✅ |
| release 套件 | `cargo test -p giac-core poly_roots::tests --release` | **25 passed, 0 ignored** |
| solve | `cargo test -p giac-solve --release` | **27 passed** |
| F4 tight bound | `field_session_dimension_bound_quartic_tight` | ✅ **dim≤24**（`t⁴+t+1` S₄） |

---

## 7. 交叉引用

| 文档 | 关系 |
|------|------|
| [GIAC-poly-p3-6-roots-algorithm-spec.md](../issues_resolved/GIAC-poly-p3-6-roots-algorithm-spec.md) | normative 算法（**resolved**） |
| [GIAC-poly-p3-6-roots-algorithm-spec.md](GIAC-poly-p3-6-roots-algorithm-spec.md) | §12 plan（F5 / F4′，open） |
| [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) | F1–F5（§0.1 待同步） |
| [GIAC-solve-polynomial-pipeline-plan](GIAC-solve-polynomial-pipeline-plan.md) | S0–S7 |

---

## 8. 变更日志

| 日期 | 变更 |
|------|------|
| 2026-06-23 | 初版；C1–C4 索引 |
| 2026-06-23 | 链算法规格 §3 |
| 2026-06-23 | **F4/S4/S6：** Euler 第二开方门禁改 `D_HARD=24`；`t⁴+t+1` 实测 dim=24；S4/S6 绿 |
| 2026-06-23 | **算法规格归档：** normative §0–§11 → `issues_resolved/`；§12 plan 留 issues |
