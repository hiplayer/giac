# GIAC-core — `giac-core` 对 upstream 主要功能缺口总览

**状态:** open
**类型:** 索引 / AFK 跟踪
**上游基线:** **`giac/giac-2.0.0`**（`src/alg_ext.cc`、`gausspol.cc`、`usual.cc`、`subst.cc`、`misc.cc`、`sym2poly.cc`、`prog.cc`）；**§4 数论另对标 Pari/GP**（`pari/src/basemath/{buch2,hnf_snf,bnfunits}.c` 等，giac-2.0.0 无对标）
**相关:** [GIAC-algext-adoption](GIAC-algext-adoption.md)、[GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md)、[GIAC-poly-f5-fglm-hasse-lean4-verification](GIAC-poly-f5-fglm-hasse-lean4-verification.md)、[GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md)、[GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md)、[giac-core-algebra-api-stability](../giac-core-algebra-api-stability.md)
**Rust 落点:** `giac-rs/crates/giac-core/src/{eval,expr,simplify,context,stmt,display}.rs` + `algebra/*`（29 文件，~35k 行）
**快照:** 2026-07-09

---

## 问题陈述

本 issue 跟踪 **`giac-core` crate 相对 upstream giac-2.0.0 主要功能的未移植 / 未完成缺口**。`giac-core` 是中央 crate，管：

- **Expr / eval / display / context / stmt** — 表达式树、求值分发、显示、上下文、语句执行
- **代数数类型系统** — `AlgExt` / `AlgExtC` / `ExtensionTower` / `FieldSession`（对标 upstream `ref_algext` / `_EXT`）
- **K 上多项式算法** — `poly_alg_*`、`poly_roots`（系数为代数数时的 gcd / factor / partfrac / sturm / resultant / 求根）
- **代数数论基础设施** — `class_group` / `bnf` / `grh` / `hnf_spec` / `ideal` / `unit_group` / `padic` / `lattice` / `archimedean` / `number_field_arith` / `galois_automorphism`（Hasse √-判定 + **Pari `bnf*` 对标线**，详见 [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md)）
- **`simplify`** — 目前仅 AST flatten，**非** upstream `subst.cc` 化简链
- **plugin trait** — algebra / calculus / linalg / ode / solve（实际算法在其它 crate）

**不在本 issue 范围**（由其它 crate 跟踪）：

| 域 | 跟踪 issue |
|----|-----------|
| limit / series / integrate / diff / risch | [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md) §1/§3、[GIAC-limit-layered-pipeline](GIAC-limit-layered-pipeline.md) |
| factor / gcd / partfrac over ℚ | [GIAC-simplify-poly-upstream-gaps](GIAC-simplify-poly-upstream-gaps.md) §2 |
| solve / sturm / froot / realroot | [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md) §4 |
| desolve | [GIAC-algorithm-gaps-open](GIAC-algorithm-gaps-open.md) §3 CAL-G5 |
| linalg | phase3 |
| groebner | phase5 |

**判定准则：** 缺口以 **upstream `at_*` 能否在 giac-core eval 落地** + **能否砍掉过渡特判** 为准，非行数比例。

---

## 0. 与 upstream giac-2.0.0 对照

### 0.1 模块映射与规模

| 域 | giac-2.0.0 源 | ~行数 | giac-core 对应 | ~行数 | 对照结论 |
|----|--------------|-------|---------------|-------|----------|
| 代数扩域 | `alg_ext.cc` | 2124 | `algebra/{alg_ext,alg_ext_c,ext_tower,field_arith,field_session,common_minimal,compositum_session}` | ~7700 | 塔式 `_EXT` + `AlgExtC` **中期基座已落地**；eval fold / `i` 进塔 ✅；evalf ✅（C-12，`evalf(AlgExt/AlgExtC)`→decimal `Rat`/`Complex`） |
| K 上多项式 | `gausspol.cc`（`_EXT` 系数段） | — | `algebra/poly_alg_*`、`poly_roots` | ~6100 | deg≤4 roots ✅；gcd/factor T0–T3 ✅；deg≥5 solve 路径 `rootof` 一支 ✅；partfrac 重根+非线性 ✅；quartic Euler √-决策 ✅ |
| 化简 | `subst.cc::simplify`、`usual.cc` `tlin`/`halftan`/`lin` | — | `simplify.rs` | 299 | **仅 AST flatten**；`tlin` NotImplemented |
| assume/purge | `usual.cc::giac_assume`、`prog.cc::_purge` | — | `context.rs`、`stmt.rs` | ~710 | 语句级栈 ✅；关系假设 + `symbol_roles` ✅；`check_assume` solve 侧 ✅；e2r/factor 参数系数阻塞于 C-16 |
| 数值求根 | `misc.cc::proot` | — | `eval.rs:eval_proot` | — | `NotImplemented("proot")` |
| 代数数论 | upstream giac **无对应**；**Pari/GP 有完整 `nf*`/`ideal*`/`bnf*`/`bnr*`** | — | `class_group`（~6800）/ `hnf_spec`（~2530，移植 `hnf_snf.c`）/ `bnf` / `grh` / `ideal` / `unit_group` / `padic` / `lattice` / `archimedean` / `number_field_arith` / `galois_automorphism` | ~17k+ | **giac 无对标；Pari 对标 ◐**：Buchmann 主线 R8–R26 已落地（GRH `LIMC2` + `rnd_rel` + `HnfSpecState` + `grh_hr_check`）；deg-2 极大序类群 100% sound；deg≥3 h=1 + 锚域 ℚ(∛11) h=2 已出证；**用户命令** `class_number`/`class_group`/`ideal*`/`bnfisprincipal`/`bnfunits`/`bnfregulator`/`bnfisunit`/`nfdisc`/`nfbasis` 等已暴露；仍缺 `bnfinit` 一体对象、`bnr*` 整层、非极大序、大域 h>1 普遍覆盖（见 [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md)） |

**结论：** 相对 **giac-2.0.0**，对齐缺口集中在 **(1) simplify 真化简链**、**(2) proot**、**(3) C-4d/C-16 参数系数 `check_assume`**（AlgExtC eval / assume 语句级 / K 上多项式主体已 ✅）。**§4 数论**单独对标 **Pari**（非 giac upstream），按 [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md) 跟踪。

### 0.2 已对齐 upstream 的代表性能力（非缺口）

避免「全盘未实现」误解：

- `AlgExt` 同域 / 跨域 `+−×`、`inv_EXT`（Phase A ✅）
- `ExtensionTower::common` + lazy cache（B-05 ✅，对标 `common_EXT`）
- `AlgExtC` 类型 + `canonicalize_to_algext_c`（B-06 部分 ✅）
- `Poly<AlgExtC>` 表示 + Expr 桥（P1 ✅）
- K 上 gcd / div_rem / monic / egcd / sqff / split_quadratic（T0–T3 ✅）
- `poly_algext_roots` deg 1–4 + `solve_poly` deg≤4（P3-6 ✅）
- K 上 resultant / sturm 计数（P3-4 部分 ✅）
- partfrac disc>0 二次分裂 + integrate K 回落（P4-2/3 ✅）
- partfrac over K 重根 + 非线性（C-7 ✅）
- quartic Euler √-决策 / A₄ dim≤12（C-8 ✅）
- Pari 对标：`class_number`/`class_group`/`ideal*`/`bnfisprincipal` 等用户命令（§4.1 ✅）

---

## 1. AlgExtC eval 接线（P0）

对标 upstream `ext_add` / `ext_mul` / `inv_EXT` 在 `eval` 主路径的完整接入。

| ID | 能力 | 现状 | 上游 | 阻塞 |
|----|------|------|------|------|
| **C-1** | `eval` 对 `AlgExtC` fold `+−×` | ✅ `eval_add`/`eval_mul` 经 `any_complex_algext` 门控走 `alg_ext_c::fold_algextc_{sum,product}` | `ext_add`/`ext_mul` | `fold_complex_algext_*` 已删；`rootof.rs` 过渡特判 / `try_factor_quadratic_rootof` 已无 |
| **C-2** | `eval_frac` 分母为 `AlgExtC` | ✅ `eval_frac` `AlgExtC` 分母支走 `AlgExtC::inv` | `inv_EXT` | — |
| **C-3** | `i` 进塔（`adj(t²+1)`） | ✅ `canonicalize_to_algext_c` 接收 `sym("i")` → `AlgExtC{0,1,ℚ}`；`to_expr` 无递归（`coords_as_rational` 经 `embed_rational` 定位常数位） | `complex_mode` + `_EXT` | — |

**验收：** `eval((1+√2*i)+(3-√2*i))` 直接走 `AlgExtC::add`；`eval(1/(1+√2*i))` 走 `AlgExtC::inv`；`fold_complex_algext_*` 删除。

**跟踪：** [GIAC-algext-adoption](GIAC-algext-adoption.md) §8.7 阶段 2b 收尾。

---

## 2. assume / purge 语句级（P0）

对标 upstream `usual.cc::giac_assume`、`prog.cc::_purge`。

| ID | 能力 | 现状 | 阻塞 |
|----|------|------|------|
| **C-4a** | 语句级 `assume` / `purge` 栈 | ✅ `stmt.rs::exec_stmt`/`exec_script` 已落地（`eval.rs` 裸 `eval()` 调用保留 guard） | CK-INT-50/54、CAL-G4、SOL-G5（阻塞于 C-4d） |
| **C-4b** | 关系假设 `A≠0`、`A>B` | ✅ `ParsedRelation` + `Context::is_assumed_{nonzero,positive,negative}` | `solve(A*x+B=0,x)` 需 `assume(A≠0)`（阻塞于 C-4d） |
| **C-4c** | `symbol_roles`（Parameter vs Variable） | ✅ `Context::symbol_roles` + `assume(sym,"parameter"\|"variable")` + `role_of`/`is_parameter` | `diff(...,x)` 时 A,B 视为常数（阻塞于 C-4d） |
| **C-4d** | 算法侧 `check_assume`（e2r / factor / solve 前） | **PR1 ✅**（`abs(var)` 区间简化）；**PR1.5 ✅**（`abs(arg)` arg 符号推理：var/var^k/ln(u)，`abs(ln(var²))` from `|var|>1`）；**C-13 slice ✅**（`exp(c*ln(u))→u^c` + `sqrt(var²)→var` + `fold_ratio`）；**CK-INT-50 ✅ enabled**；**solve 侧 ✅**（`filter_roots_by_var_sign`：`solve(x²=4,x)`+`assume(x>0)`→`[2]`，按 `ctx.is_assumed_{positive,negative}(var)` 过滤常数根，代数/rootof 根未知符号则保留）；**残留**：e2r/factor 侧参数系数 `check_assume` 阻塞于 C-16（`PolyCoeff`，符号 A,B 当系数；`solve(A*x+B=0,x)+assume(A≠0)` 现 `expr_to_poly` 把 A,B 当变量→factor fail）、AlgExt/rootof 根符号过滤需 C-12 `evalf` | `sym2poly.cc::check_assume` 对标 |

**跟踪：** [GIAC-algext-adoption](GIAC-algext-adoption.md) §8.9、[parser-token-map.md](../parser-token-map.md) §6、GIAC-204b。

---

## 3. K 上多项式管线收尾（P1）

对标 upstream `gausspol.cc` `ext_factor` / `ext_factor_nodegck` 完整路径。

| ID | 能力 | 现状 | 上游 | 阻塞 |
|----|------|------|------|------|
| **C-5** | `Poly<AlgExtC>` deg≥5 factor | solve 路径已走 `irreducible_rootof_branch_algext`（`rootof(α,P)` 一支）；`poly_algext_roots` deg≥5 `NotImplemented` **按设计保留**（Abel–Ruffini） | `ext_factor` 高次 → `rootof(α,P)` 一支 | —（已落地） |
| **C-6** | T3+ `adjoin(K, u²−α)` parent 域系数层 | **主体 ✅**（T3+a `adjoin(K₁,u²−√2)` 登记 + `element_*` parent-coeff 算术 + flatten + `solve(t⁴+t+1)` conformance 全过；`min_poly_parent_blocks: Vec<CoordsQ>` 存 parent 域 op 坐标，对齐 upstream `_EXT` minpoly 系数=gen）；**残留 compositum-of-two-parent-block-towers ✅**（`is_embedded_rational` 蕴含方向修正——识别 parent-block 塔的嵌套有理常数位；`verify_embedded_generator` ParentBlocks 分支复合 `try_subfield_embedding(parent,operand) ∘ embedding_for(operand,pair)`——共享父非直接 operand 时；`build_adjoin_parent_coeffs` 形态校验 monic/block-len/degree；测试 `compositum_two_parent_block_towers_over_shared_parent`）；**类型诚实度 P2**（`min_poly_parent_blocks` 裸 `CoordsQ` 靠运行时校验，升 `ParentBlock`/`FieldElement` 绑定型见 [GIAC-field-element-bound-model](GIAC-field-element-bound-model.md)） | `common_EXT` 嵌套路径（minpoly 系数=gen，原生支持任意嵌套深度） | —（已落地；原 gap doc "多处 NotImplemented / 维度爆" 描述过时，2026-07 校正） |
| **C-7** | partfrac over K 重根 / 非线性 | **重根 ✅**（C-7a：`partfrac_affine_power_system_over_k` K-线性系统，cover-up 不处理 multiplicity>1）；**非线性 ✅**（C-7b：irreducible-over-K degree≥2 因子 → degree-(d−1) 多项式分子，同一 K-线性系统泛化，d≤3；与 upstream `gausspol.cc` `pf` deg≤2 闭式 / `Tpartfrac` Taylor 数学等价，按项目规则绑数学语义非行级实现） | `sym2poly.cc` partfrac + `ext_factor` | `integrate` K 路径因子链已通 |
| **C-8** | quartic Euler resolvent √-决策边界 | **已落地**（2026-07 校正）：`try_sqrt_in_field` 在 dim-12 塔域 ℚ(α,β) 找到 √Δ_Q（A₄ splitting field=dim 12，√Δ_Q ∈ 该域），`quartic_roots_by_adjoin_deflate` 不 adjoin、dim 稳定 12；`quartic_a4_galois_dim_le_12` ✅（非 ignore，dim≤`POLY_ROOTS_DIM_QUARTIC_OUT`）；Euler path `euler_four_roots_vanish`/`roots_quartic_t4_plus_t_plus_1` ✅。过时注释已清理。原 doc "多处 NotImplemented / A4 dim>12 / 求根失败" 描述过时（NotImplemented 为 fallback 错误路径，正常不触达） | `gausspol.cc` 不走此路径（upstream 用 `proot` 数值） | —（已落地） |

**跟踪：** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) P3-5 / F4′、[GIAC-poly-f5-fglm-hasse-lean4-verification](GIAC-poly-f5-fglm-hasse-lean4-verification.md)。

---

## 4. 代数数论基础设施深化（P2）

**注：** upstream giac-2.0.0 **无对应实现**；本节 **upstream 对标 = Pari/GP**（`nf*` / `ideal*` / `bnf*` / `bnr*`）。起源为 Hasse √-判定扩展，现按 [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md) 主线推进 Buchmann。**数学正确性**见 [GIAC-poly-f5-fglm-hasse-lean4-verification](GIAC-poly-f5-fglm-hasse-lean4-verification.md)；**已暴露命令清单**见 [GIAC-p2-algebraic-number-theory-api](GIAC-p2-algebraic-number-theory-api.md)。

### 4.1 已暴露 Pari 对标命令（eval，2026-07-09）

`nfdisc` / `nfrootsof1` / `nfsign` / `nfbasis` · `class_number` / `class_group` · `bnfunits` / `bnfregulator` / `bnfisunit` · `idealhnf` / `idealmul` / `idealnorm` / `idealred` / `idealinv` / `idealpow` · `bnfisprincipal` · `lll` / `evalf`

**仍无用户命令：** `bnfinit`（一体缓存对象）、`bnfcond`、`idealaddtoone` / `idealchinese`、整个 `bnr*` / `galois*` 族。

### 4.2 缺口表（对标 Pari）

| ID | 能力 | 现状 | upgrade path |
|----|------|------|--------------|
| **C-9** | 完整 Buchmann / `bnfinit` | **◐ R8–R38 已落地**：… + R38 FP 幂基 Minkowski 行缓存 | `bnfinit` 用户对象；非极大序；冷启动 perf（~5.7s vs Pari ~3ms） |
| **C-10** | 非主性认证 / `bnfisprincipal` 完备 | **◐** deg-2 极大 + M_K-smooth ✅；deg≥3 h=1 ✅；h>1 ℚ(∛11) `(2)→(None,[1])` / `(4)→γ` ✅；`q>M_K` sound-skip ✅ | 更多 deg≥3 h>1 锚域；大素理想纳入关系库 |
| **C-11** | LLL 短向量 | **✅ slice + 类群接入**：`eval_lll`（`lll(matrix)`）；`idealred` / `ideal_is_principal` LLL 回退（R7）；`enumerate_relations_lli`（R8） | Minkowski 嵌入 LLL（Pari `idealred` 风格）；`qfminim` |
| **C-12** | `AlgExtC::evalf` | **✅ 已落地**：`archimedean::algext_evalf`/`algextc_evalf` + `eval_evalf` dispatch | `horner_rootof` / `proot` 浮点逼近 |
| **C-18** | `bnfinit` 用户对象 | 内部 `Bnf` + `CertClassData` + `Context` session 缓存；**无** `FuncKind::Bnfinit` | 对标 Pari `buchall_end` 可查询对象；`bnfisunit`/`bnfregulator` 只读缓存 |
| **C-19** | `bnr*` 射线类群 + 类域论 | ❌ 完全无 | 新 `ray_class.rs`；依赖 `idealaddtoone`/`idealchinese` |
| **C-20** | Pari 工程残差 | ◐ R39c–p NF/FP + `idealHNF_mul` 整数 HNF ✅；grow 仍 ~13.6s | 见 [GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md) **R39+ backlog** |

**`ponytail:` 现状边界：** R7 已砍掉 deg≥5 / `N(J)>10⁷` 主性 sound-skip（LLL 回退）；现存 sound-skip 主因 = GRH 关系不足、非极大序、`power_order_is_maximal≠true`、无 `PREC` 环。deg-2 有界路径 `BUCH_COORD_BOUND=64` 仍约束大 regulator 实二次 h>1 **结构**（`class_number` 仍经解析公式返 h）。

---

## 5. simplify 真化简链 / 显示 / 参数化（P3）

对标 upstream `subst.cc::simplify`、`usual.cc` 三角化简链。

| ID | 能力 | 现状 | 上游 |
|----|------|------|------|
| **C-13** | `simplify` 真化简链 | `simplify.rs`（299 行）仅 flatten Add/Mul + 合并数字系数；**C-13 slice ✅**（C-4d 路径局部 `exp(c*ln(u))→u^c` + `sqrt(var²)→var` + `fold_ratio`，`eval_integrate` 入口用，非全局 simplify 链） | `subst.cc::simplify`（`reorder`/`canonical_form`/`evalf` 化简） |
| **C-14** | `tlin` | `eval.rs` `NotImplemented("tlin")` | `usual.cc::tlin` 三角线性化 |
| **C-15** | `proot` 数值求根 | `eval.rs` `NotImplemented("proot")` | `misc.cc::proot` |
| **C-16** | 参数系数 A,B（`PolyCoeff`） | `min_poly`/`coords` 不能含符号参数 | 阶段 3 / Phase C；与 C-4 协同 |
| **C-17** | `rootof` 显示稳定排序 | A-05 partial；golden 字符串顺序仍需稳定化 | `symb_rootof` |

**跟踪：** [GIAC-simplify-poly-upstream-gaps](GIAC-simplify-poly-upstream-gaps.md) §1.3、[cas-long-term-vision.md](../cas-long-term-vision.md) §5.2–5.3（参数化）。

---

## 优先级与建议顺序

```text
P0  C-1/C-2/C-3   AlgExtC eval fold + frac + i 进塔
                   → 砍 fold_complex_* / rootof.rs / try_factor_quadratic_rootof
P0  C-4           assume/purge 语句级 + 关系假设 + symbol_roles + check_assume
                   → 解锁 CK-INT-50/54、CAL-G4、SOL-G5
P1  C-5/C-6/C-7   K 上 deg≥5 factor / T3+ adjoin / partfrac 重根+非线性（C-5/C-6/C-7 ✅）
P1  C-8           quartic Euler √-决策边界（Hasse lean4 验证线）
P2  C-9..C-12,C-18..C-20  Pari Buchmann / bnfinit / bnr*（非 conformance 硬阻塞；细节见 GIAC-p2-bnf-pari-alignment R27+）
P3  C-13..C-17    simplify 真化简链 / tlin / proot / 参数化 / 显示稳定
```

---

## 分阶段验收

### Phase A — 文档（本 issue）

- [x] 汇总 `giac-core` vs upstream 主要功能缺口
- [x] 与 upstream giac-2.0.0 对照（[§0](#0-与-upstream-giac-200-对照)）
- [x] 链接子 issue（algext-adoption、poly-algext-backlog、hasse-lean4）
- [ ] 新 `giac-core` 算法债 PR 须更新本表或子 issue，不得只加 `NotImplemented` 无条目

### Phase B — P0 AlgExtC eval 接线

- [x] C-1 `eval` AlgExtC fold `+−×`
- [x] C-2 `eval_frac` AlgExtC 分母
- [x] C-3 `i` 进塔（`adj(t²+1)`）
- [x] 删除 `fold_complex_algext_*` / `rootof.rs` 过渡特判 / `try_factor_quadratic_rootof`
- [x] [GIAC-algext-adoption](GIAC-algext-adoption.md) 阶段 2b 收尾

### Phase C — P0 assume/purge

- [x] C-4a 语句级 `assume`/`purge` 栈（`stmt.rs` — `exec_stmt`/`exec_script` 已落地；`eval.rs` 裸调用 guard 保留）
- [x] C-4b 关系假设 `A≠0`、`A>B`（`ParsedRelation` + `Context::is_assumed_{nonzero,positive,negative}`）
- [x] C-4c `symbol_roles`（Parameter vs Variable；`assume(sym, "parameter"|"variable")`）
- [~] C-4d `check_assume` 接 e2r / factor / solve（跨 crate，分 PR）— **PR1 ✅** abs(var)；**PR1.5 ✅** abs(arg) arg 符号；**C-13 slice ✅** exp(c*ln(u))→u^c + sqrt(var²)→var + fold_ratio；**CK-INT-50 ✅ enabled**；**solve 侧 ✅** filter_roots_by_var_sign（assume(var sign) 过滤常数根）；e2r/factor 参数系数侧阻塞于 C-16、AlgExt/rootof 根符号过滤阻塞于 C-12，待续
- [~] CK-INT-50/54 enabled — **CK-INT-50 ✅**（conformance 绿）；**CK-INT-54 暂缓**（`integrate(ln(x²+t²)/(1+t²),t,0,∞)`=π·ln(1+x) 阻塞于 `integrate_frac` ln-numerator + 无穷限特殊积分，非 C-4d scope；通用 Risch transcendental 是大工程，特化 hack 不雅，待后续评估）

### Phase D — P1 K 上管线收尾

- [x] C-5 `Poly<AlgExtC>` deg≥5 factor → `rootof(α,P)` 一支（solve 路径 `irreducible_rootof_branch_algext` 已落地；`poly_algext_roots` deg≥5 按设计 NotImplemented）
- [x] C-6 T3+ `adjoin(K, u²−α)` parent 域系数层（主体已实现：登记 + `element_*` + flatten + `solve(t⁴+t+1)` 全过；残留 compositum-of-two-parent-block-towers 已补：`is_embedded_rational` 蕴含修正 + `verify_embedded_generator` 父嵌入复合 + `build_adjoin_parent_coeffs` 校验 + 测试 `compositum_two_parent_block_towers_over_shared_parent`；类型诚实度 P2 跟踪 [GIAC-field-element-bound-model](GIAC-field-element-bound-model.md)）
- [x] C-7a partfrac over K **重根**（`partfrac_affine_power_system_over_k` K-线性系统）
- [x] C-7b partfrac over K **非线性**（irreducible-over-K degree≥2 → degree-(d−1) 多项式分子；K-线性系统泛化 d≤3；测试 `partfrac_nonlinear_quadratic_over_k` `x/((x²−2)(x²−3)) → −x/(x²−2)+x/(x²−3)`）
- [x] C-8 quartic Euler √-决策边界（A4 Galois dim≤12）— **已落地**（2026-07 校正）：`try_sqrt_in_field` 在 dim-12 塔域 ℚ(α,β) 找到 √Δ_Q（`diag_a4_sqrt_probe_gap` 实测 dim-12 √Δ_Q probe = Some），`quartic_roots_by_adjoin_deflate` 不 adjoin、dim 稳定 12；`quartic_a4_galois_dim_le_12` ✅（非 ignore）；Euler path `euler_four_roots_vanish`/`roots_quartic_t4_plus_t_plus_1` ✅。过时注释（F4′.3b "probe MISSES"、reframe "currently None"）已清理。`try_galois_sqrt_second` 保留作 Euler fallback 子路径（adjoin-deflate 优先时不触达，但 Euler test 可达）。

### Phase E — P2 代数数论深化（对标 Pari）

- [~] C-9 完整 Buchmann / `bnfinit` — **◐ R8–R36**：… + R36 grow 热路径 embedding 去重；锚域 release ~5.7s
- [~] C-10 非主性认证 — **◐** `bnfisprincipal` deg-2 极大 + deg≥3 h=1 ✅；h>1 非主金值 / 大素理想待续
- [x] C-11 LLL — **已接入类群**（`eval_lll` + `ideal_is_principal` LLL 回退 + `enumerate_relations_lli` + `idealred`）
- [x] C-12 `AlgExtC::evalf` — **已落地**（详见 [GIAC-p2-algebraic-number-theory-api](GIAC-p2-algebraic-number-theory-api.md)）
- [ ] C-18 `bnfinit` 用户对象
- [ ] C-19 `bnr*` 射线类群 + 类域论
- [~] C-20 Pari 工程残差（`rnd_rel_par` 串行入库 ✅ / Buchmann perf R35+）

### Phase F — P3 simplify / 数值 / 参数化

- [~] C-13 `simplify` 真化简链 — **slice ✅**（C-4d 路径局部 exp/log/sqrt，`eval_integrate` 入口用）；完整 `subst.cc::simplify` 链（reorder/canonical_form/evalf 化简）仍 P3 未做
- [ ] C-14 `tlin`
- [ ] C-15 `proot`
- [ ] C-16 参数系数 A,B（`PolyCoeff` + Phase C）
- [ ] C-17 `rootof` 显示稳定排序

---

## 验证命令

```bash
cd giac-rs
cargo test-timeout                                # 默认跳过 ignored
cargo test -p giac-core alg_ext                   # AlgExt / AlgExtC 单元
cargo test -p giac-core ext_tower                 # 塔 / common
cargo test -p giac-core poly_roots                # K 上求根
cargo test -p giac-core class_group              # Buchmann / GRH / bnfisprincipal
cargo test -p giac-core grh                     # GRHchk / HR 证书
cargo test -p giac-core hnf_spec                # Pari hnfspec/hnfadd 移植
cargo test -p giac-core ideal                   # 理想 HNF 算术
cargo test -p giac-core class_number_general_cert_grh  # ℚ(∛11) h=2 等探针
cargo test -p giac-conformance giac_check_cas     # 含 rootof / AlgExt conformance
```

---

## 参考

- 上游基线（giac）：**`giac/giac-2.0.0`** — `src/alg_ext.cc`、`gausspol.cc`、`usual.cc`、`subst.cc`、`misc.cc`、`sym2poly.cc`、`prog.cc`
- 上游基线（§4 数论）：**Pari/GP** — `pari/src/basemath/buch2.c`、`hnf_snf.c`、`bnfunits.c`；对照 [giac-buchmann-classical-vs-pari](../giac-buchmann-classical-vs-pari.md)
- check 黄金：`giac/giac-2.0.0/check/testcas`（`rootof` 用例）、`testintegrate` L51–55（assume/purge）
- API 分层：[giac-core-algebra-api-stability.md](../giac-core-algebra-api-stability.md)
- 子 issue：[GIAC-algext-adoption](GIAC-algext-adoption.md) §8、[GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md)、[GIAC-poly-f5-fglm-hasse-lean4-verification](GIAC-poly-f5-fglm-hasse-lean4-verification.md)、**[GIAC-p2-bnf-pari-alignment](GIAC-p2-bnf-pari-alignment.md)**（Pari Buchmann 主线，R8–R27+）
- 已知偏离：[known-divergences.md](../known-divergences.md)
