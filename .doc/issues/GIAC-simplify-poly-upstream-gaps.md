# GIAC-simplify-poly — giac-simplify / giac-poly / giac-calculus 上游对齐缺口与 API 分层

**状态:** open  
**类型:** 索引 / AFK 跟踪  
**上游基线:** **`giac/giac-2.0.0`**（`check/` 黄金 + `usual.cc` / `sym2poly.cc` / `gausspol.cc` / `intg.cc` / `series.cc`）  
**相关:** [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)、[GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md)、[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)、[giac-simplify-api-stability.md](../giac-simplify-api-stability.md)、[giac-poly-api-stability.md](../giac-poly-api-stability.md)  
**验收:** FAC-G1–G3 tracer un-ignore ✅；simplify/poly 稳定 API 表与源码 `/// **Stable**` 一致（进行中）

**快照日期:** 2026-06-24

---

## 问题陈述

giac-rs Phase 4 目标是 **headless CAS 子集**，不是 giac-2.0.0 全库 1:1。本 issue 汇总 **giac-simplify、giac-poly、giac-calculus** 相对 upstream 的算法缺口，并指向各 crate 的 **稳定 / Partial / 临时 API** 分层文档。

缺口判定以 **check 黄金行能否 eval** 为准；禁止新增 per-case 形状函数（见 [GIAC-limit-exp-diff P3](GIAC-limit-exp-difference-unification.md#p3--禁止-per-gruntz-形状表)）。

---

## 0. 当前回归快照（2026-06-24）

| 域 | check / 单测 | giac-rs | 说明 |
|----|--------------|---------|------|
| **limit** | Maxima rtest 17 条 | **17/17 ✅** | 含 Gruntz（`gruntz_exp_nested_diff` 等） |
| **limit** | CK-INT-55–66 | **enabled ✅** | conformance eval 门禁 |
| **integrate** | `testintegrate` 表 | **35 enabled / 31 disabled** | enabled eval 门禁绿（含 CK-INT-05） |
| **factor** | `testfactor` tracer | **9/9 ✅**（`--include-ignored`） | L16–L17、L20–L22、L24–L26 全绿 |
| **normalize** | `testnormalize` 188 条 | 经 `normal` 间接覆盖 | `non_recursive_normal` 未注册 |

**Phase B 结论（FAC-G 门禁）：** testfactor L20/L22/L25/L26 均已 un-ignore 且绿；FAC-G1–G3 **tracer 验收已闭合**。一般参系数 `poly_factor` 塔、3+ aux `sparse_bi`、`factor_multivariate` 边界 `reverse()` 仍为 Partial（见 §2.2、§5）。

---

## 1. giac-simplify — 上游对齐

**上游:** `usual.cc`、`sym2poly.cc`、`lin.cc`、`subst.cc`（见 [module-division.md](../module-division.md) §2）

### 1.1 已对齐（Stable）

| 能力 | 上游 | Rust |
|------|------|------|
| `normal` / `expand` | `usual.cc` | `expand.rs` — 多项式路径；`ExpandPolicy::NoExpDistribute` 为 limit 保留 exp 形状 |
| `ratnormal` | 有理式规范化 | `ratnormal.rs` — 不含 `exp`/`sin` 等超越子树 |
| `factor` | `ezgcd.cc` 入口 | `factor.rs` → `giac-poly` |
| `ifactor` | `ifactor.cc` | `ifactor.rs` |
| 等价判定 | conformance 规格 | `equiv.rs` — `assert_equiv` / `is_zero` |

### 1.2 部分对齐（Partial）

| API | 缺口 |
|-----|------|
| `texpand` | 负整数倍角 → `NotImplemented`；无完整 `tlin` 链 |
| `lin` | 仅 `(exp(x)+1)^2`、`exp(a)*exp(b)`、有界 exp 幂 |
| `halftan` | 仅 `sin(2*x)/(1+cos(2*x))` |
| `factor` | 依赖 giac-poly 一般多元塔（FAC-G2 一般情形）；`try_factor_quadratic_sqrt` 仍为 Expr 侧临时路径 |

### 1.3 未移植

`builtin-api-map.md` 规划在 giac-simplify、crate 内 **无实现**：

`tlin`, `trig2exp`, `trigsin`, `trigcos`, `tan2sincos`, `tcollect`, `tsimplify`, `lncollect`, `acos2asin` 族, `reorder`, `truncate`, `epsilon2zero`, **`non_recursive_normal`**

**语义分叉:** `simplify` 在 `giac-core` 仅为 AST flatten，**不是** upstream `subst.cc` 化简链。

**技术债:** `ratnormal` 对 `AlgExt` 走 `eval` stub（[GIAC-algext-adoption](GIAC-algext-adoption.md) A-03）。

→ API 分层见 [giac-simplify-api-stability.md](../giac-simplify-api-stability.md)

---

## 2. giac-poly — 上游对齐

**上游:** `gausspol.cc`、`ezgcd.cc`、`modfactor.cc`、`sym2poly.cc`（partfrac）

### 2.1 已对齐（Stable 子集）

| 能力 | 说明 |
|------|------|
| `Poly` 环运算 | `quo`/`rem`/`egcd`/`content`/`gauss` |
| 一元分解 | Zassenhaus + 有理根；cyclotomic / 模式表 |
| 二元 Hensel | `try_hensel_lift_bivariate`；L16–L22 tracer 绿 |
| 模分解 | `factor_poly_mod` / `factor_fpx` |
| partfrac | `partfrac_rational_terms` — 线性因子链；deg≤3 不可约 ✅ |
| Sturm / 结果式 | `sturm_*`、`resultant`、`tresultant_*`（Rothstein–Trager 用） |
| nested-ring P0 | `BivariateEmbed`、`SqffRingCtx`、`FlatUni` — 见 [GIAC-poly-p0-backlog](GIAC-poly-p0-backlog.md) §1 ✅ |

### 2.2 算法缺口（FAC-G*）

| ID | upstream | giac-rs | check 锚点 |
|----|----------|---------|------------|
| **FAC-G1** | `try_sparse_factor` + `try_sparse_factor_bi` + `unitaryfactor` | **Partial** — P0–P2b ✅；line25/26 gate ✅；U5 `reverse` 边界、3+ aux `sparse_bi` 仍缺 | L25 unitary；L22 由 Hensel 主路径覆盖 |
| **FAC-G2** | 参三元 `poly_factor` 塔 | **Partial** — `try_lift_factors_in_aux_var` 覆盖 L20 ✅；一般参系数塔仍缺 | testfactor L20 ✅ |
| **FAC-G3** | 二元混合次数 Hensel 或 fallback | **Partial** — L22 ✅（`hensel_lift_two_at_zero` @ `main=y,aux=x`） | testfactor L22 ✅ |

**partfrac 连带缺口：**

| 情形 | 算法 | 管线 |
|------|------|------|
| sqff 因子 `deg > 3` | `NotImplemented` | giac-poly ℚ |
| 实二次 `disc > 0` | **giac-core K 切片 ✅** | integrate **未接** K 路径；giac-poly ℚ 仍 `NotImplemented` |
| 线性 / 重根 / `disc≤0` 二次 | ✅ | giac-poly ℚ |

→ API 分层见 [giac-poly-api-stability.md](../giac-poly-api-stability.md)  
→ P0 剩余项见 [GIAC-poly-p0-backlog](GIAC-poly-p0-backlog.md)

### 2.3 表示层（架构，非上游算法）

[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md) Phase 1–3 **已勾选**；剩余：

- **U5:** `factor_multivariate` 边界 `reverse_var_order`（`|vars|≥3` Partial；二元 reverse `p` 会破坏 line25/26）
- **sparse_bi 3+ aux:** 当前 aux 两两配对；三 aux 同 embed 仍缺
- **P2+:** `Poly<AlgExtC>` 上 gcd/factor/partfrac 闭环 — [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md)

---

## 3. giac-calculus — 上游对齐（摘要）

**上游:** `derive.cc`、`intg.cc`、`risch.cc`、`series.cc`

本 issue 主责 simplify/poly；calculus 细节仍见 [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)。摘要：

| 域 | 测试现状 | 算法深度缺口 |
|----|----------|--------------|
| **limit** | Maxima 17/17、CK limit 绿 | MRV SparseSeries ordre/upscale/padd 未完整（LIM-G2–G4）；CK-61 仍部分依赖形状特例；asymptotic 快路径待退役 |
| **integrate** | 35/66 enabled | 真 Risch（CAL-G1）；disc>0 有理积分 **算法在、管线未接**；31 disabled 换元/assume |
| **diff** | 基础超越 | 一般幂、asin/sqrt 等 |
| **series** | Taylor 子集 | 深层 ln/exp @ 非 0（SER-G1/G2） |
| **assume/purge** | 未实现 | CK-INT-50/54（CAL-G4，giac-prog/core） |

---

## 4. 跨 crate 依赖

```text
giac-simplify::factor ──► giac-poly::factor_into
                              │
                              ├─ integrate 未接 K partfrac ──► ∫1/(x²−k) 类失败（算法在 giac-core ✅）
                              ├─ FAC-G1 尾部（U5 / 3+ aux）──► 多元 factor 覆盖率
                              └─ FAC-G2 一般塔 ──► 参系数一般情形

giac-simplify::normal/ratnormal ──► giac-calculus::limit_engine
giac-simplify::expand ──► giac-calculus::integrate
giac-core::partfrac_rational_terms_over_k ──► （待接）giac-calculus::partfrac_integrate
giac-poly::partfrac (ℚ) ──► giac-calculus::integrate（Hermite / RT / 有理项）
```

---

## 5. 优先级（2026-06-24）

**原则：** P0 = **算法能力核对 → 管线接线 / 算法补全**；禁止先 enable conformance 或 un-ignore（见 [GIAC-algorithm-gaps-open §8](GIAC-algorithm-gaps-open.md#8-p0-判定算法--接线--测试)）。

```text
P0  CAL-G2 接线   ✅ integrate K 回落（2026-06-24）；残余：混合 sqff / 非常数分子

P0  LIM-G2/G3     MRV 主路径收敛
                  非：先扩 check 门禁

P1  FAC-G1 尾部   U5 reverse；sparse_bi 三 aux
P1  CAL-G1        Risch 真主流程
P1  API tier 同步 annotate_api_tiers ↔ stability doc

P2  partfrac deg>3 高次不可约；simplify tlin/trig2exp；integrate 31 disabled
P3  assume/purge；LIM-G4/G5
```

### 5.1 已闭合（勿再排 P0）

| 项 | 验收 |
|----|------|
| FAC-G2 L20 | `testfactor_line20` ✅ |
| FAC-G3 L22 | `testfactor_line22` ✅ |
| FAC-G1 line25/26 | `testfactor_line25` / `line26` ✅ |
| LIM-G1 嵌套 exp | `gruntz_exp_nested_diff` ✅ |
| CK-INT-05 | integrate 表 enabled eval ✅ |
| nested-ring P0 | [GIAC-poly-p0-backlog](GIAC-poly-p0-backlog.md) §1 ✅ |

---

## 6. 分阶段验收

### Phase A — 文档（本 issue）

- [x] 汇总 simplify / poly / calculus 上游缺口
- [x] [giac-simplify-api-stability.md](../giac-simplify-api-stability.md) — 稳定 / Partial / 临时 API
- [x] [giac-poly-api-stability.md](../giac-poly-api-stability.md) — 同上
- [x] `annotate_api_tiers.py` 覆盖 simplify / poly / calculus / core-algebra / solve / ode / groebner（2026-06-23）
- [ ] 源码 `/// **Stable|Partial|Temporary**` 与文档表增量同步
- [ ] PR 合入前按 [algorithm-expr-api.md §7.2](../algorithm-expr-api.md#72-测试通过后提交--合入前复审) 复审
- [ ] 更新 [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md) 过时条目（FAC-G L20/L22、LIM-G1、CK-INT-05 已绿）

### Phase B — FAC-G1–G3（tracer 门禁 ✅）

- [x] `try_sparse_factor`（@0 + `find_good_eval`）
- [x] `try_sparse_factor_bi`（2-aux sum-coeff dual-embed + dilation）
- [x] 参三元 factor tower（FAC-G2 L20 由 aux-lift 覆盖）
- [x] tracer L22 / L20 / L25 / L26 un-ignore
- [x] 退役 `try_factor_bivariate_eval` / `try_kronecker_bivariate` 热路径

#### Phase B 架构诊断（2026-06-24）

**L22 现状（FAC-G3 ✅）：**

| 路径 | L22 结果 | 说明 |
|------|----------|------|
| `try_sparse_factor` (x,y) / (y,x) | None | `lcp` 为 poly；模板解出 spurious 分支，verify 正确拒绝 |
| `try_hensel_lift_factor` (x,y) | None | `lcp=3y+9` 依赖 aux → upstream 同拒 |
| `hensel_lift_at_zero` (y,x) | **✅** | `hensel_lift_two_at_zero` 在 `main=y,aux=x` 下完成混合次数 lift |
| `try_factor_bivariate_eval` | 退役 | `#[cfg(test)]` |
| `try_kronecker_bivariate` | 退役 | 同上 |

**upstream 正确顺序**（`gausspol.cc` `do_factor_hensel` L6855–7048）：

```text
sqff → 随机 eval 不可约快检 → try_sparse_factor(v0) → try_sparse_factor_bi
     → try_hensel_lift_factor(p, F0, v0)  // 含 bivariate lc 先验
     → unitaryfactor / pzadic 启发式
```

**仍缺（降为 P1，非 tracer 阻塞）：**

- **U5:** `factor_multivariate` 边界 `reverse()` — 详见 [GIAC-poly-unitaryfactor-gaps](GIAC-poly-unitaryfactor-gaps.md) §U5
- **sparse_bi 3+ aux:** 三 aux 同 embed sum-coeff（当前仅两两配对）
- **FAC-G2 一般情形:** 完整参系数 `poly_factor` 塔（L20 特判已绿）

**临时算法可否删除：**

| 位置 | 函数 | 结论 |
|------|------|------|
| `giac-simplify` | `ratnormal_algext`, `canonical_radical`, `inv_sqrt_to_mul` | **暂留** — AlgExt / equiv 契约未替代 |
| `giac-simplify` | `try_factor_quadratic_sqrt` | **暂留** — giac-poly 未覆盖 Expr 侧 sqrt |
| `giac-poly/poly_uni` | `try_factor_bivariate_eval`, `try_kronecker_*` | **已退役**（`#[cfg(test)]`） |
| `giac-poly/poly_uni` | 静默 `Ok(vec![g])` | **语义正确**（不可约）— 链未含随机 eval / 全 sparse_bi / 边界 reverse，属能力缺口非 bug |

### Phase C — simplify 三角/化简链

- [ ] `tlin` / `trig2exp` 最小子集
- [ ] `non_recursive_normal` 注册

### Phase D — partfrac K 路径接线

- [x] **算法切片** `partfrac_rational_terms_over_k` + `factor_univariate_over_k`
- [x] **P4-3** `integrate_rational_partfrac`：ℚ 失败 → K 回落（2026-06-24）
- [x] **P4-2** `eval_partfrac`：ℚ 失败 → K 回落
- [ ] 扩展路由：非常数分子、sqff 链中单块 disc>0 二次
- [ ] conformance 评估（算法闭合后再 enable）

---

## 7. 验证命令

```bash
cd giac-rs
cargo test-timeout
cargo test -p giac-calculus maxima_rtest --lib
cargo test -p giac-conformance giac_check_integrate_table_enabled --test giac_check_integrate
cargo test -p giac-poly testfactor_line -- --include-ignored
cargo test -p giac-poly --lib unitary sparse_factor partfrac
```

---

## 8. 参考

- 模块映射：[module-division.md](../module-division.md)
- 表示层规范：[algorithm-expr-api.md](../algorithm-expr-api.md)
- calculus 分层：[giac-calculus-api-stability.md](../giac-calculus-api-stability.md)
- MRV 后续：[GIAC-limit-mrv-followup](GIAC-limit-mrv-followup.md)
- AlgExt 待办：[GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md)
