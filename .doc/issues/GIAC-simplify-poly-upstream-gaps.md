# GIAC-simplify-poly — giac-simplify / giac-poly / giac-calculus 上游对齐缺口与 API 分层

**状态:** open  
**类型:** 索引 / AFK 跟踪  
**上游基线:** **`giac/giac-2.0.0`**（`check/` 黄金 + `usual.cc` / `sym2poly.cc` / `gausspol.cc` / `intg.cc` / `series.cc`）  
**相关:** [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)、[GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md)、[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)、[giac-simplify-api-stability.md](../giac-simplify-api-stability.md)、[giac-poly-api-stability.md](../giac-poly-api-stability.md)  
**验收:** FAC-G1–G3 tracer un-ignore；simplify/poly 稳定 API 表与源码 `/// **Stable**` 一致

**快照日期:** 2026-06-19

---

## 问题陈述

giac-rs Phase 4 目标是 **headless CAS 子集**，不是 giac-2.0.0 全库 1:1。本 issue 汇总 **giac-simplify、giac-poly、giac-calculus** 相对 upstream 的算法缺口，并指向各 crate 的 **稳定 / Partial / 临时 API** 分层文档。

缺口判定以 **check 黄金行能否 eval** 为准；禁止新增 per-case 形状函数（见 [GIAC-limit-exp-diff P3](GIAC-limit-exp-difference-unification.md#p3--禁止-per-gruntz-形状表)）。

---

## 0. 当前回归快照（2026-06-19）

| 域 | check / 单测 | giac-rs | 说明 |
|----|--------------|---------|------|
| **limit** | Maxima rtest 17 条 | **17/17 ✅** | 含 Gruntz（`gruntz_exp_nested_diff` 等） |
| **limit** | CK-INT-55–66 | **enabled ✅** | conformance eval 门禁 |
| **integrate** | `testintegrate` 表 | **35 enabled / 31 disabled** | enabled eval 门禁绿（含 CK-INT-05） |
| **factor** | `testfactor` tracer | L20/L22 **enabled ✅** | — |
| **normalize** | `testnormalize` 188 条 | 经 `normal` 间接覆盖 | `non_recursive_normal` 未注册 |

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
| `factor` | 依赖 giac-poly FAC-G1–G3 |

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
| 二元 Hensel | `try_hensel_lift_bivariate` — 多数 check L16–L21 绿 |
| 模分解 | `factor_poly_mod` / `factor_fpx` |
| partfrac | `partfrac_rational_terms` — 线性因子链 |
| Sturm / 结果式 | `sturm_*`、`resultant`、`tresultant_*`（Rothstein–Trager 用） |

### 2.2 算法缺口（FAC-G*）

| ID | upstream | giac-rs | check 锚点 |
|----|----------|---------|------------|
| **FAC-G1** | `try_sparse_factor` + `try_sparse_factor_bi` 启发式 fallback | **Partial** — `@0`/`find_good_eval`；`sparse_bi` 2-aux sum-coeff ✅（`reconstruct_factor_dual_embed`）；`unitaryfactor` P0–P2a ✅（line25 gate） | L22 Hensel；line25 unitary |
| **FAC-G2** | 参三元 `poly_factor` 塔 | **Partial** — `try_lift_factors_in_aux_var` 覆盖 L20 | testfactor L20 ✅ |
| **FAC-G3** | 二元混合次数 Hensel 或 fallback | **Partial** — `try_hensel_lift_factor` + `(main,aux)` 双序 | testfactor L22 ✅ |

**partfrac 连带缺口:** 高次因子仍 `NotImplemented`；重复/实二次已部分覆盖（`disc=0` 重根、`disc>0` 有理分裂）。

→ API 分层见 [giac-poly-api-stability.md](../giac-poly-api-stability.md)

### 2.3 表示层缺口（架构，非上游算法）

`Poly` 擦除环上下文导致 `div_rem` / `quo_exact_wrt` 混用（FAC-G1 tri_var sum-coeff 等）。类型化路线图见 **[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)**：

- P0：`BivariateEmbed`、`EmbedFactorDraft`、`factor/*` 禁嵌套环 `div_rem`
- P1：`SqffRingCtx`、`FactorSet`、`SparseAtZero`、`HenselPair`
- P2：`FlatUni` / `MultivariatePoly` 边界分裂

---

## 3. giac-calculus — 上游对齐（摘要）

**上游:** `derive.cc`、`intg.cc`、`risch.cc`、`series.cc`

本 issue 主责 simplify/poly；calculus 细节仍见 [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)。摘要：

| 域 | 测试现状 | 算法深度缺口 |
|----|----------|--------------|
| **limit** | Maxima 17/17、CK limit 绿 | MRV SparseSeries ordre/upscale/padd 未完整；asymptotic 快路径待退役（LIM-G2–G4） |
| **integrate** | 35/66 enabled | 真 Risch（CAL-G1）；31 disabled（换元/反三角/多参/assume）；partfrac 已补 `∫P/(cx+d)`（`fdeg==1, ndeg>=1`） |
| **diff** | 基础超越 | 一般幂、asin/sqrt 等 |
| **series** | Taylor 子集 | 深层 ln/exp @ 非 0（SER-G1/G2） |
| **assume/purge** | 未实现 | CK-INT-50/54（CAL-G4，giac-prog/core） |

---

## 4. 跨 crate 依赖

```text
giac-simplify::factor ──► giac-poly::factor_into
                              │
                              ├─ FAC-G1/G3 失败 ──► partfrac 缺口 ──► integrate 失败
                              └─ FAC-G2 参系数塔 ──► testfactor L20

giac-simplify::normal/ratnormal ──► giac-calculus::limit_engine
giac-simplify::expand ──► giac-calculus::integrate
```

---

## 5. 优先级

```text
P0  giac-poly FAC-G1/G3     → un-ignore testfactor L22；解锁 partfrac/积分
P1  giac-poly FAC-G2        → un-ignore testfactor L20
P1  giac-calculus limit     → 216e SparseSeries 主路径；退役 asymptotic 快路径
P1  giac-calculus Risch     → CAL-G1 真主流程
P2  giac-simplify           → tlin/trig2exp/non_recursive_normal
P2  giac-calculus integrate → 换元表、反三角、ibp
P3  assume/purge            → CAL-G4（giac-prog）
```

---

## 6. 分阶段验收

### Phase A — 文档（本 issue）

- [x] 汇总 simplify / poly / calculus 上游缺口
- [x] [giac-simplify-api-stability.md](../giac-simplify-api-stability.md) — 稳定 / Partial / 临时 API
- [x] [giac-poly-api-stability.md](../giac-poly-api-stability.md) — 同上
- [ ] 源码 `/// **Stable|Partial|Temporary**` 与文档表同步（incremental PR）
- [ ] PR 合入前按 [algorithm-expr-api.md §6.2](../algorithm-expr-api.md#62-测试通过后提交--合入前复审) 复审：临时匹配净减少 + 新增 fn tier
- [ ] 更新 [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md) 过时条目（LIM-G1、CK-INT-05 已绿）

### Phase B — FAC-G1–G3

- [x] `try_sparse_factor`（@0 + `find_good_eval` 好点重试；L22 仍 None — 见架构注）
- [x] `try_sparse_factor_bi`（2-aux：`eval_tn` + sum-coeff dual-embed + dilation；`others.len()>=2` 接入）
- [x] 参三元 factor tower（FAC-G2 L20 由 aux-lift 覆盖）
- [x] tracer L22 un-ignore
- [x] tracer L20 un-ignore
- [x] 退役 `try_factor_bivariate_eval` / `try_kronecker_bivariate` 热路径

#### Phase B 架构诊断（2026-06-19，更新）

**L22 现状（FAC-G3 ✅）：**

| 路径 | L22 结果 | 说明 |
|------|----------|------|
| `try_sparse_factor` (x,y) / (y,x) | None | `lcp` 为 poly；模板解出 spurious `∏f=lcp·p` 分支，verify 正确拒绝 |
| `try_hensel_lift_factor` (x,y) | None | `lcp=3y+9` 依赖 aux → upstream 同拒 |
| `hensel_lift_at_zero` (y,x) | **✅** | `hensel_lift_two_at_zero` 在 `main=y,aux=x` 下完成混合次数 lift |
| `try_factor_bivariate_eval` | 退役 | 与 upstream 链重复，eval 匹配不可靠 |
| `try_kronecker_bivariate` | 退役 | L22 上 23s+ 仍失败 |

**upstream 正确顺序**（`gausspol.cc` `do_factor_hensel` L6855–7048）：

```text
sqff → 随机 eval 不可约快检 → try_sparse_factor(v0) → try_sparse_factor_bi
     → try_hensel_lift_factor(p, F0, v0)  // 含 bivariate lc 先验
     → unitaryfactor / pzadic 启发式
```

giac-rs 缺口：

- **已对齐：** 双主元 `(main,aux)`；`try_sparse_factor` + `find_good_eval`；`try_sparse_factor_bi`（2-aux sum-coeff + dilation）；`try_hensel_lift_factor` + `hensel_lift_two_at_zero`；`unitaryfactor` P0–P2a（line25 gate ✅）
- **仍缺：** `factor_multivariate` 边界 `reverse()`（U5）；3+ aux `sparse_bi`；完整参系数 `poly_factor` 塔（FAC-G2 一般情形）。详见 [GIAC-poly-unitaryfactor-gaps](GIAC-poly-unitaryfactor-gaps.md)

**临时算法可否删除：**

| 位置 | 函数 | 结论 |
|------|------|------|
| `giac-simplify` | `ratnormal_algext`, `canonical_radical`, `inv_sqrt_to_mul` | **暂留** — AlgExt / equiv 契约未替代 |
| `giac-simplify` | `try_factor_quadratic_*` | **暂留** — giac-poly 未覆盖 Expr 侧 sqrt/rootof |
| `giac-poly/poly_uni` | `try_factor_bivariate_eval`, `try_kronecker_*` | **已退役**（`#[cfg(test)]`）— FAC-G3 由 sparse→Hensel 覆盖 |
| `giac-poly/poly_uni` | 静默 `Ok(vec![g])` | **语义正确**（不可约）— 但链未含随机 eval / sparse_bi / pzadic，属能力缺口非 bug |

**FAC-G2 L20**（参系数 `b,c`）：`try_lift_factors_in_aux_var` 已覆盖 L20；完整 `poly_factor` 系数环塔仍缺（一般参系数情形）。

### Phase C — simplify 三角/化简链

- [ ] `tlin` / `trig2exp` 最小子集
- [ ] `non_recursive_normal` 注册

---

## 7. 验证命令

```bash
cd giac-rs
cargo test-timeout
cargo test -p giac-calculus maxima_rtest --lib
cargo test -p giac-conformance giac_check_integrate_table_enabled --test giac_check_integrate
cargo test -p giac-poly testfactor_line -- --include-ignored
```

---

## 8. 参考

- 模块映射：[module-division.md](../module-division.md)
- 表示层规范：[algorithm-expr-api.md](../algorithm-expr-api.md)
- calculus 分层：[giac-calculus-api-stability.md](../giac-calculus-api-stability.md)
