# GIAC-open-algo — giac-rs 当前算法未实现缺口总览

**状态:** open  
**类型:** 索引 / AFK 跟踪  
**上游基线:** **`giac/giac-2.0.0`**（`check/` 黄金回归 + `src/*.cc` 算法实现；giac-rs conformance 默认 `GIAC_VERSION_DIR` 指向此目录）  
**相关:** [GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md)、[GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md)、[GIAC-limit-exp-difference-unification](GIAC-limit-exp-difference-unification.md)、[phase4-issues.md](../phase4-issues.md) §2.3  
**验收:** 本 issue 中 **P0** 项对应的 `#[ignore]` 清零或改为 enabled conformance

---

## 问题陈述

本 issue 汇总 **用户可见的 `NotImplemented` / 上游 giac-2.0.0 已有而 giac-rs 未移植的算法**，按域分类，并指向子 issue 或 phase4 编号。**禁止**在此清单外新增 per-case 形状函数（见 [GIAC-limit-exp-diff P3](GIAC-limit-exp-difference-unification.md#p3--禁止-per-gruntz-形状表)）。

**快照日期:** 2026-06-18（`cargo test --workspace` 全绿 + 23 ignored；其中 **4** 个为算法 `#[ignore]`，见 [§7](#7-非算法类-ignore23-中的-19-个)）

---

## 0. 与 upstream giac-2.0.0 对照

### 0.1 模块映射与规模（估算移植度）

| 域 | giac-2.0.0 源文件 | ~行数 | giac-rs 对应 | ~行数 | 对照结论 |
|----|-------------------|-------|--------------|-------|----------|
| **极限 / 级数** | `src/series.cc` | 3713 | `limit_engine/*` + `series.rs` | ~6200 | MRV / Gruntz **骨架在**；`mrv_lead_term` 完整 ordre 循环、双向极限、完整 `limit_symbolic_preprocess` **未齐** |
| **积分** | `src/intg.cc` | 7344 | `integrate.rs` + `partfrac_integrate` + 启发式 | ~1500 | 表驱动 + partfrac/Hermite/RT **子集**；换元 / 特殊函数 / 定积分假设 **大量未移植** |
| **Risch** | `src/risch.cc` | 1107 | `risch/*.rs` | ~800 | 塔/Hermite/RT **窄路径**；`in_risch` / `remains_to_integrate` 主流程 **无** |
| **ODE** | `src/desolve.cc` | 2539 | `giac-ode/desolve.rs` | ~300 | 常系数线性 ODE **子集** |
| **分解** | `src/gausspol.cc`（`do_factor_hensel` 等） | 8137 | `giac-poly/src/factor/*` | ~6000 | 一元/二元 Zassenhaus+Hensel **多数绿**；高次变元 Hensel + **启发式** `try_sparse_factor` **缺** |
| **求解** | `src/solve.cc` | 11712 | `giac-solve/*` | ~1200 | 多项式 + `rootof` + Newton **有**；超越/系统/假设 **弱** |
| **杂项** | `src/misc.cc`（`froot` 等） | 11175 | — | — | `froot`/`froots` **未注册** |

**结论：** giac-rs Phase 4 目标是 **headless CAS 子集**，不是 2.0.0 全库 1:1；缺口应以上游 **check 黄金行能否 eval** 为准，而非行数比例。

### 0.2 check 黄金回归对照（算法缺口锚点）

giac-rs conformance 与 triple 测试对齐 `giac-2.0.0/check/`（见 `tests/conformance/src/lib.rs` `GIAC_VERSION_DIR`）。

#### `check/testlimit`（26 条 `limit`）

| 行 | 上游表达式（摘要） | giac-2.0.0 | giac-rs | 缺口 ID |
|----|-------------------|------------|---------|---------|
| L2 | `exp(x)*(exp(1/x+…+exp(-x²))-exp(1/x-exp(-exp(x))))` @ inf | ✅ → `1` | ❌ **`#[ignore]`** | **LIM-G1** |
| L1 | `exp(x)*(exp(1/x-exp(-x))-exp(1/x))` @ inf | ✅ → `-1` | ✅（preprocess + 快路径） | — |
| L10,L21 | `(exp(inner)-exp(x))/x` @ inf（CK-INT-61） | ✅ → `-exp(2)` | ✅（MRV + w 层特例） | LIM-G2 通用路径仍缺 |
| L15 | `exp(inner)/exp(x)` @ inf（比值） | ✅ → `1` | ✅（`try_limit_exp_over_exp_via_quotient`） | 待 MRV 主路径（Phase C） |
| L11 | `(3^x+5^x)^(1/x)` | ✅ | ✅ | — |
| 其余 ~21 条 | 代数 / MRV / 倒数级数 | ✅ | ✅（Maxima 14 子集 12/14） | LIM-G5/G6 对更全表仍可能缺 |

上游 **`series.cc`** 能力 giac-rs **尚未完整移植**（216e 路线 A）：

- `mrv_lead_term` 内 **upscale** 循环（`ln(x)→x`, `x→exp(x)`）及 `exp(x)` 防递归守卫
- ordre **`ordre*1.5+1` 升阶**重算级数（非仅增大参数名）
- `remove_lnexp` → `ln_expand` / `exp_series` 系数 **`padd` 合并**
- `limit_symbolic_preprocess`（trig 归一、`factorial2gamma`、直接代入探测）
- `+∞` **unidirectional** 主路径（**不做** `x=1/u`）；有限点单侧才倒数换元

#### `check/testintegrate`（67 条，含 limit/series 行）

| 行 | 内容 | 上游 | giac-rs enabled | 缺口 ID |
|----|------|------|-----------------|---------|
| L5 | `integrate(1/(3*x*(x²+x+1)*(x-1)³),x)` | ✅ | ❌ eval 失败 | **CAL-G2** / GIAC-212b |
| L51–52 | `assume(t>2),integrate(…),purge(t)` | ✅ | parse only | **CAL-G4** / GIAC-204b |
| L55 | `assume(x>0),integrate(ln(…),t,0,inf),purge(x)` | ✅ | parse only | **CAL-G4** |
| L61–62 | CK-INT-60/61 极限（同 testlimit） | ✅ | ✅ conformance | — |
| L1–4,6–50 等 | 有理 / 三角 / exp 积分 | ✅ 多数 | ✅ 多数（~30/31 enabled eval） | CAL-G1 真 Risch 仍无 |

#### `check/testfactor`（31 条）

| 行 | 表达式 | 上游 | giac-rs | 缺口 ID |
|----|--------|------|---------|---------|
| L20 | `(x+b+c)*(x²-xb+b²-xc-bc+c²)` 参三元 | ✅ | ❌ **`#[ignore]`** | **FAC-G2** |
| L22 | `(3x-y²+y-5)*(xy+3x-y²-1)` 二元混合次数 | ✅ | ❌ **`#[ignore]`** | **FAC-G3** |
| L17–19,21 | 其它二元/三元 | ✅ | ✅ tracer 覆盖 | — |

上游 **`gausspol.cc` `do_factor_hensel`** 在 Hensel 失败时有 **`try_sparse_factor` 启发式**（注释 L6863–6864）；giac-rs **无等价 fallback**，故 FAC-G1/G3 在相同行上 fail。

#### 其它 check 文件（索引）

| check 文件 | 条数级 | giac-rs 接入 | 主要缺口 |
|------------|--------|--------------|----------|
| `testlimit` | 26 | CK-INT-55–66 部分 enabled | 全表 SymPy 门 **`#[ignore]`**（非算法） |
| `testintegrate` | 67 | **31** enabled eval | assume/purge、CK-INT-05 |
| `testpartfrac` | 1 | partfrac 单测 | 高阶因子 |
| `flanex` 等 | — | 未全接 | **SOL-G1** `froots` |

### 0.3 缺口 ↔ 上游函数速查

| 缺口 ID | giac-2.0.0 已有（入口） | check 锚点 |
|---------|-------------------------|------------|
| LIM-G1–G4 | `series.cc`: `unidirectional_limit`, `mrv_lead_term`, `remove_lnexp`, `series__SPOL1` | `testlimit` L2 |
| FAC-G1,G3 | `gausspol.cc`: `do_factor_hensel`, `try_hensel_lift_factor`, `try_sparse_factor` | `testfactor` L22 |
| FAC-G2 | `gausspol.cc`: 多变量 `factor` 塔 | `testfactor` L20 |
| CAL-G1 | `intg.cc`/`risch.cc`: `do_risch`, `in_risch` | 一般超越积分 |
| CAL-G2 | `intg.cc` + `gausspol` factor | `testintegrate` L5 |
| CAL-G4 | `usual.cc` `giac_assume`, `prog.cc` `_purge` | L51–52, L55 |
| SOL-G1 | `misc.cc` `_froot` | `flanex` |
| SOL-G5 | 同上 assume/purge | 同 CAL-G4 |

### 0.4 giac-rs 已对齐 2.0.0 的代表性能力（非缺口）

避免「全盘未实现」误解——以下在 **同一 check 行** 上与上游一致或 SymPy 验证通过：

- `limit(sin(x)/x,x,0)`、`(1-cos)/x²`、`(1+1/n)^n` @ inf
- CK-INT-55–59、62–66 limit/series 行（conformance enabled）
- CK-INT-60/61 Gruntz 极限（`testintegrate` L61–62）
- `testfactor` L17–19、L21 及多数一元分解
- Phase 3 线代（`test_linalg*` bin 脚本对标）

---

## 1. 极限 / 级数 / Gruntz（P0）

| ID | 能力 | 现状 | 阻塞测试 / 用例 | 跟踪 |
|----|------|------|-----------------|------|
| **LIM-G1** | 嵌套 `exp` 差分 @ `+∞` | `NotImplemented("limit")`；fold 可工作，MRV 主路径未通 | `maxima_rtest::gruntz_exp_nested_diff` **`#[ignore]`** | [GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md)、[GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md) G3/G4 |
| **LIM-G2** | MRV 换元后 `SparseSeries` 主项收敛 | `series_lead_at_zero` 复杂式回退 lead-only 特例 | CK-INT-61 当前绿（特例链）；通用路径未收敛 | [GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md) |
| **LIM-G3** | `(-ln(w))⁻¹` / `x⁻¹` 换元后的 padd 消去 | 主项系数残留 `ln(w)` | 嵌套 gruntz、部分 CK 类 | 216e §已知实现缺口 |
| **LIM-G4** | upstream ordre 递增 / upscale / spdiv | `try_order` 被 cap，升阶无效 | 慢路径 / 超时风险 | 216e 路线 A |
| **LIM-G5** | 一般 `-∞` / 双向 / 方向极限 | 部分 `NotImplemented("limit")` | Maxima 子集外 | [GIAC-limit-maxima-upstream-alignment](GIAC-limit-maxima-upstream-alignment.md) |
| **LIM-G6** | `mrv_compare` 完整比较 | 部分用 `f64` 估值 | `3^x` vs `5^x` 等已绿；更一般式待补 | limit-maxima §后续 |

**Maxima gruntz 14 条：12/14 ✅；仍 ignore 1 条（LIM-G1）。**

**目标管线（勿再加 `try_limit_gruntz_*`）：**

```text
preprocess → MRV → remove_lnexp → SparseSeries lead → limit
```

临时快路径（待 Phase C 退役）：`exp_diff::try_limit_*_preprocessed` — [GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md) G1。

---

## 2. 多项式分解（P1）

| ID | 能力 | 现状 | 阻塞测试 | 跟踪 |
|----|------|------|----------|------|
| **FAC-G1** | 高次变元 Hensel 提升 | 混合次数二元式失败 | `hensel::hensel_line22_mixed_bivariate` **`#[ignore]`** Issue 2.3 | `giac-poly/src/factor/hensel.rs` |
| **FAC-G2** | 参三元 `factor` tower | 无 `poly_factor` 塔 | `tracer::testfactor_line20_parametric_cubic_factor` **`#[ignore]`** Issue 2.4 | `giac-poly/src/factor/tracer.rs` |
| **FAC-G3** | 二元混合次数（Hensel 或启发式） | 同 FAC-G1 | `tracer::testfactor_line22_bivariate_mixed_degree` **`#[ignore]`** Issue 2.3 | 同上 |

**影响：** 高阶/多变量 `factor` → 积分 `partfrac` 管道（**GIAC-212b**）、一般有理式分解。

---

## 3. 微积分 — 积分 / 微分 / Risch（P1–P2）

| ID | 能力 | 现状 | 代表用例 | 跟踪 |
|----|------|------|----------|------|
| **CAL-G1** | 真 Risch 主流程 | `eval_risch` ≡ 启发式 `integrate` | — | **GIAC-229–231** |
| **CAL-G2** | 一般有理式 `integrate frac` | 依赖 `factor` + partfrac 缺口 | **CK-INT-05**、**INT-A06** | **GIAC-212b**、**GIAC-210** |
| **CAL-G3** | `derive(atan(...))` 等 | `NotImplemented: diff` | 扩展三角微分 | phase4 §2.3 |
| **CAL-G4** | 含参 / 反常定积分 | 无 `assume`/`purge` | CK-INT-50/54 | **GIAC-204b** |
| **CAL-G5** | 一般 `desolve` | 常系数 ODE 子集 ✅；非线性/变系数 ❌ | ODE-H03 部分 | **GIAC-218** |

---

## 4. 求解 / 根（P1–P2）

| ID | 能力 | 现状 | 代表用例 | 跟踪 |
|----|------|------|----------|------|
| **SOL-G1** | `froot` / `froots` | 未注册 | `check/flanex` L195 | **GIAC-232** |
| **SOL-G2** | `realroot` 无理根隔离 | 仅有理根 | 区间算术完整版 | **GIAC-207** |
| **SOL-G3** | 一般超越方程 `solve` | 最小子集（如 `sin(x)=0`） | 参数解、多分支 | **GIAC-208** |
| **SOL-G4** | `sturm` 重因子 / 平方因子 | `sturm((x^3+1)^2)` 等 | **GIAC-206** | `giac-solve/src/sturm.rs` |
| **SOL-G5** | 语句级 `assume` / `purge` | 未解析多语句程序 | CK-INT-50/54 | **GIAC-204b** |

---

## 5. 级数（非 MRV）（P2）

| ID | 能力 | 现状 | 说明 |
|----|------|------|------|
| **SER-G1** | 一般 `series` @ 非 0 / 含 `ln`/`exp` 深层 | 部分 `NotImplemented("series")` | 经典 Taylor 子集 ✅（**GIAC-216**） |
| **SER-G2** | `sparse_series` 完整 `in_series__SPOL1` 语义 |  bounded `MAX_SERIES_ORDER` | 与 LIM-G2 重叠 |

---

## 6. 线性代数 / 其它（P3）

| ID | 能力 | 现状 | 跟踪 |
|----|------|------|------|
| **LA-G1** | 符号特征值 / Jordan 完整 | 部分子集 | phase3 |
| **MISC-G1** | `gbasis`（Groebner 基） | 仅 `greduce` | phase5 |
| **MISC-G2** | WASM / 插件全覆盖 | 冒烟级 | **GIAC-222** |

---

## 7. 非算法类 ignore（23 中的 ~19 个）

**不计入本 issue 的算法债** — CI / 运维策略：

| 类别 | 数量 | 原因 |
|------|------|------|
| `giac_probe_*` | 12 | 手动 SymPy 积分探针 |
| `giac_check_integrate*` 全量 SymPy | 3 | SymPy 慢/易挂；CI 用逐条 `giac_check_integrate_ck_int_*` |
| `giac_check_integrate_inventory` | 3 | 仅 eval 清单，人工启用行 |
| `giac_check_limit` 全量报告 | 1 | 全表 SymPy |
| `phase4_integrate_table` SymPy 门禁 | 1 | 同上；eval 门启用 `giac_check_integrate_table_enabled` |

运行 ignored：`cargo test --workspace -- --include-ignored` 或 `cargo test-timeout -- --include-ignored`（需 nextest 传参）。

---

## 优先级与建议顺序

```text
P0  LIM-G1/G2/G3/G4  →  un-ignore gruntz_exp_nested_diff；退役 exp_diff try_limit_*
P1  FAC-G1/G2/G3      →  un-ignore 3 个 factor 测试；解锁 GIAC-212b
P1  CAL-G1/G2         →  Risch + 一般 partfrac 积分
P1  SOL-G1/G4/G5      →  flanex / sturm / assume
P2  其余
```

---

## 分阶段验收

### Phase A — 文档（本 issue）

- [x] 汇总算法缺口 vs 非算法 ignore
- [x] 与 upstream **giac-2.0.0** 对照（[§0](#0-与-upstream-giac-200-对照)）
- [x] 链接子 issue（216e、limit-pipeline、phase4）
- [ ] 新算法债 PR 须更新本表或子 issue，不得只加 `#[ignore]` 无条目

### Phase B — P0 极限

- [ ] [GIAC-216e](GIAC-216e-mrv-series-lead-convergence.md) 路线 A 落地
- [ ] `gruntz_exp_nested_diff` un-ignore
- [ ] [GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md) Phase 3：`try_limit_*` 退役

### Phase C — P1 factor + 积分

- [ ] FAC-G1/G3：Hensel 或启发式 fallback
- [ ] FAC-G2：parametric factor tower
- [ ] CAL-G2：CK-INT-05 / INT-A06 enabled

---

## 验证命令

```bash
cd giac-rs
cargo test-timeout                                    # 默认跳过 23 ignored
cargo test -p giac-calculus maxima_rtest -- --include-ignored
cargo test -p giac-poly hensel_line22 testfactor_line -- --include-ignored
cargo test -p giac-conformance giac_check_integrate_ck_int -- --include-ignored  # 无；均已 enabled
```

---

## 参考

- `NotImplemented` 分布：`giac-calculus/limit_engine/*`、`integrate.rs`、`giac-poly/factor/*`、`giac-solve/*`
                     - 上游基线：**`giac/giac-2.0.0`** — `src/series.cc`（limit/MRV）、`src/gausspol.cc`（factor/Hensel）、`src/intg.cc` + `src/risch.cc`（积分）、`src/usual.cc`/`src/prog.cc`（assume/purge）、`src/misc.cc`（`froot`）
- check 黄金：`giac/giac-2.0.0/check/testlimit`、`testintegrate`、`testfactor`
- 回归：`giac-rs/crates/giac-calculus/src/limit.rs` `maxima_rtest`；`tests/conformance/`（`GIAC_VERSION_DIR` → `giac-2.0.0`）
