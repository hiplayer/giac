# GIAC-poly-nested-ring — 嵌套环表示层类型化与除法 API 分家

**状态:** open  
**类型:** 架构 / 表示层  
**上游基线:** **`giac/giac-2.0.0`**（`gausspol.cc` `do_factor_hensel` / `try_sparse_factor` / `try_sparse_factor_bi`）  
**相关:** [giac-poly-api-stability.md](../giac-poly-api-stability.md) §1、[GIAC-simplify-poly-upstream-gaps](GIAC-simplify-poly-upstream-gaps.md) §2、Cursor [giac-poly-nested-ring.mdc](../../.cursor/rules/giac-poly-nested-ring.mdc)  
**验收:** `factor/*` 嵌套环热路径无裸 `Poly::div_rem`；embed/sparse 重建不经 `Poly` 往返；`nested::tests` + `testfactor_line` / `sparse_factor` 全绿

**快照日期:** 2026-06-19

---

## 问题陈述

`giac-poly` 的 `Poly` 是 **intentionally erased** 的通用稀疏表示：同一类型承载 ℚ[x]、ℚ[others][main]、embed 像 ℚ[main,t] 等。环结构、主元、管线阶段默认不在类型里，须靠 API 契约与调用方纪律维持。

FAC-G1/G3 排错中反复暴露：**`Poly::div_rem`（多元 leading-monomial 商）≠ ℚ[others][main] 整除**（`quo_exact_wrt`）。根因是类型未编码环上下文，导致在 `factor/*` 用错除法语义（如 `rest.div_rem(f)` 验因子、`lcpt/lc` 用 `div_rem`）。

已落地 MVP（`giac-rs/crates/giac-poly/src/nested.rs`）：

| 数学对象 | Rust 类型 | 整除 / 商 |
|----------|-----------|-----------|
| ℚ[others][main] | `MainVar` + `UnivariateIn` / `UnivariatePoly` | `.divides()` / `.exact_quo_dividing()` |
| ℚ[others] 系数 | `CoeffRingPoly` | `.exact_quo()` |
| sparse_bi embed 配置 | `TnEmbed` | `.embed()` → 裸 `Poly`（**缺口**） |
| (main,t) 重建 IR | `EmbedMonomial`（crate 内） | monomial 逐项对比 |

本 issue 跟踪 **将嵌套环数学编码进 Rust 类型**、**除法 API 按环分家**、**各阶段专用 IR**，避免 per-case `if` 补丁。

---

## 背景：为何是类型问题

管线排错顺序（优先于改 bug）：

1. **在哪个环？** 多元商环 / ℚ[others][main] / ℚ[others] / embed 像
2. **主元是谁？** `main`（及 aux、t、embed 的 `n`）
3. **哪一阶段？** sqff / sparse / sparse_bi / Hensel / 启发式 fallback

只有三者清楚后，才问「`div_rem` 为何余式非零」。多数是**除法语义用错**，不是算术 bug。

回归锚点：`nested::tests::univariate_in_divides_vs_div_rem` — 证明对嵌套环对象 `div_rem` 与 `quo_exact_wrt` 结论相反。

---

## 目标（分阶段）

### Phase 1 — 类型 + 门禁（P0）

#### 1.1 `UnivariateOver<Var>` 命名对齐

- 现状：`UnivariateIn` + `MainVar`（语义即用户所称 `UnivariateOver<Var>`）
- 动作：可选 `type alias UnivariateOver<'a> = UnivariateIn<'a>` 或统一重命名；**禁止** `pub poly` 直接 `.div_rem` 于新代码

#### 1.2 `BivariateEmbed`（嵌入**结果**类型）

今天 `TnEmbed::embed() -> Poly` 丢失 `(main, t, n, aux)` 上下文。

```rust
// 目标形态（示意）
struct BivariateEmbed {
    poly: Poly,
    embed: TnEmbed,
}
// 仅暴露 .view_main() / .view_t() / .to_recon_draft()
// 禁止隐式 From<BivariateEmbed> for Poly
```

#### 1.3 `EmbedFactorDraft`（embed 阶段专用 IR）

```rust
struct EmbedFactorDraft {
    monos: Vec<EmbedMonomial>,
    aux_exps: Vec<(u64, u64)>,  // 与 monos 等长
    embed: TnEmbed,
    seldegs: Vec<u64>,
}
// fn materialize(self) -> Option<UnivariatePoly>
```

monomial 循环只改 draft，**禁止**中途 `Poly ↔ embed` 往返。

#### 1.4 除法 API 按环分家 + 门禁

| 环 | 允许 API | 禁止 |
|----|----------|------|
| ℚ[others][main] | `UnivariateIn::divides` / `quo_exact_wrt` | `Poly::div_rem` |
| ℚ[others] 系数 | `CoeffRingPoly::exact_quo` | `div_rem` 归一化 lc |
| ℚ[main] 平坦一元 | `FlatUni`（新）或专用 `div_rem_wrt` | 与 nested 混用 |
| 多元展示环 | `MultivariatePoly` 边界类型 | 渗入 sparse/Hensel 热路径 |

**门禁：**

- `factor/*` 对嵌套环对象禁止 `Poly::div_rem` — `clippy::disallowed_methods` 或 `NestedPoly` 私有 wrapper（不 `Deref` 到 `Poly`）
- 登记 [giac-poly-api-stability.md](../giac-poly-api-stability.md)

**待迁移 `div_rem` 落点（2026-06-19 审计）：**

| 文件 | 行/函数 | 目标 API |
|------|---------|----------|
| `sparse.rs` | ~1474 heuristic `rest.div_rem(f)` | `UnivariateIn::divides` |
| `hensel.rs` | `div_rem_x_over_qy` | `UnivariateIn::div_rem_wrt_aux_indep` |
| `univariate.rs` | 有理根 `rest.div_rem(&lin)` | `FlatUni` |
| `zassenhaus.rs` | `rat_div_rem` 链 | 保留在模/平坦一元上下文，标 `FlatUni` |
| `fpx.rs` / `power.rs` / `cyclotomic.rs` / `util.rs` | 各自上下文 | 标环边界或迁 `FlatUni` |

---

### Phase 2 — 管线上下文（P1）

#### 2.1 `SqffRingCtx` / `FactorTower`

封装 `factor_sqff_over_coeff_ring(g, var, others, …)` 的裸 `&Var` + `&[Var]`：

```text
{ poly, main: MainVar, others: AuxVars }
```

`FactorRecFn` 改为 `fn(SqffRingCtx) -> …`；`try_sparse_factor` / `try_hensel_lift_bivariate` / 不可约快检共用。

**FAC-G2 `PolyFactorTower`**（`factor/tower.rs`）：`factor_sqff_chain` = aux-lift → good_eval → sparse_bi（不可约快检在 sparse 前，避免 L24 类空洞）。

#### 2.2 `FactorSet` / `SqffComponent`

```rust
struct SqffComponent { factor: UnivariatePoly, multiplicity: usize }
// FactorSet::product_equals(&Poly) — 内部全 UnivariateIn::divides 链
```

替代混用 `Vec<Poly>`、`Vec<(Poly, usize)>`、`Option<Vec<Poly>>`。

#### 2.3 `SparseAtZero` + `SparseSystem`

FAC-G1 二元 sparse：`substitute_poly @ 0`、模板、`SparseEquation` 收成阶段对象，求解输出 `SparseSolution` → 单次 `UnivariatePoly`，不经 `Poly` 往返。

#### 2.4 `GoodEval` / `EvalPoint`

`find_good_eval` 返回包装「赋值后仍保持 main 次数」的 `GoodEval { values, preserved_main_degree }`。

#### 2.5 `HenselPair`

```text
{ p, main, aux, f0: UnivariateIn, g0: UnivariateIn }
```

`f0`/`g0` 构造时约束 aux 上常系数（`is_independent_of_y` 进类型）。

#### 2.6 `SubstMainAux`

`eval_aux(p: &UnivariateIn, aux, value: &CoeffRingPoly) -> UnivariateIn` — 保证 main 不变。

---

### Phase 3 — 边界环（P2）

| 类型 | 用途 |
|------|------|
| `FlatUni` | ℚ[main] 平坦一元；`factor_univariate_flat`、Zassenhaus 一元像；**唯一**允许经典 `div_rem` 的 factor 子路径之一 |
| `MultivariatePoly` | cyclotomic、展示；与 nested **无隐式 `From`** |
| `PrimitivePart<Wrt MainVar>` | `primitive_part_wrt` 结果 tagged |
| `DilationMap` | sparse_bi dilation `{ aux_a, aux_b }` + apply/undo |
| `DualEmbedSnapshot` | 双嵌入快路径；标明仅单项式系数 embed |

---

## 设计原则

1. **降级要显式** — `as_poly()` / `into_inner()`；不要 `Deref<Target=Poly>`
2. **一种数学环一种除法** — 类型上只暴露该环的商/余接口
3. **阶段 IR 不经过 `Poly` 往返** — embed / sparse 求解 / monomial 重建用 draft 类型
4. **新增类型要有误用测试** — 参考 `nested::tests::univariate_in_divides_vs_div_rem`
5. **禁止**在 `sparse`/`hensel` 调用方用 `div_rem` 打补丁绕过环语义（见 [algorithm-before-patch.mdc](../../.cursor/rules/algorithm-before-patch.mdc)）

---

## 分阶段验收清单

### Phase 1

- [x] `BivariateEmbed` + `EmbedFactorDraft` 落地；`sparse_bi` 重建改走 draft
- [x] `factor/*` 嵌套环路径 `div_rem` 清零或迁 `UnivariateIn` / `CoeffRingPoly`（`sparse` heuristic、`hensel` `div_rem_x_over_qy`）
- [x] `hensel.rs` `div_rem_x_over_qy` → `UnivariateIn::div_rem_wrt_aux_indep`
- [x] `factor/mod.rs` 嵌套环门禁注释；`lib.rs` 导出 `BivariateEmbed` / `UnivariateOver`
- [x] [giac-poly-api-stability.md](../giac-poly-api-stability.md) §2.6 + `nested.rs` 契约
- [x] `cargo test -p giac-poly nested sparse_factor testfactor_line` 全绿

### Phase 2

- [x] `SqffRingCtx` + `FactorSet` 接入 `poly_uni` 因子链
- [x] `SparseAtZero` / `SparseSystem` 替代 sparse 内裸 `Poly` 拼装
- [x] `GoodEval` + `HenselPair` + `UnivariateIn::eval_aux`
- [x] `PolyFactorTower` / `CoeffRing`（FAC-G2）；`factor_sqff_chain`；good_eval 门控进 tower
- [x] `CoeffRingPoly::gcd`；`content_wrt_impl` / `primitive_part_wrt_impl` 类型化
- [x] `factor_bivariate_flat` 替代 sparse_bi 内嵌 `factor_multivariate_rec`
- [x] `unitaryfactor` / `pzadic` 二元 MVP 迁入 typed 上下文（`factor/unitary.rs`；`poly_uni` 链接入）
- [ ] `unitaryfactor` 完整上游尾链（见 [GIAC-poly-unitaryfactor-gaps](GIAC-poly-unitaryfactor-gaps.md) U-P0–P1）
- [x] `cargo test -p giac-poly` Phase 2 回归全绿

### Phase 3

- [x] `FlatUni` vs `MultivariatePoly` 显式分裂
- [x] `DilationMap`、`PrimitivePart` 按需

---

## 验证命令

```bash
cd giac-rs
cargo test -p giac-poly nested:: --lib
cargo test -p giac-poly sparse_factor --lib
cargo test -p giac-poly testfactor_line -- --include-ignored
cargo clippy -p giac-poly -- -D warnings  # Phase 1 门禁合入后
```

---

## 参考

- 实现：`giac-rs/crates/giac-poly/src/nested.rs`
- 已接线：`factor/sparse.rs`（sparse_bi）、`factor/poly_uni.rs`（`coeff_wrt_poly`）
- Cursor 规则：[giac-poly-nested-ring.mdc](../../.cursor/rules/giac-poly-nested-ring.mdc)
- 上游链：[GIAC-simplify-poly-upstream-gaps](GIAC-simplify-poly-upstream-gaps.md) Phase B 架构诊断
- unitaryfactor 缺口：[GIAC-poly-unitaryfactor-gaps](GIAC-poly-unitaryfactor-gaps.md)
