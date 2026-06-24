# GIAC-poly — `AlgExt` / `Poly<AlgExtC>` 接入待办

**状态:** open  
**类型:** 索引 / AFK  
**相关:** [GIAC-algext-adoption](GIAC-algext-adoption.md) §8、[GIAC-poly-algext-gcd-factor-priority](../issues_resolved/GIAC-poly-algext-gcd-factor-priority.md)（T0–T3 ✅）、[GIAC-poly-p3-6-quartic-roots-gaps](GIAC-poly-p3-6-quartic-roots-gaps.md)、[GIAC-poly-p0-backlog](GIAC-poly-p0-backlog.md) §2 partfrac、[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)、[giac-poly-api-stability.md](../giac-poly-api-stability.md)  
**上游参考:** `giac/giac-1.5.0/src/alg_ext.cc`、`gausspol.cc` `algext_convert`、`_EXT` 系数多项式  
**Rust 落点:** `giac-core::algebra::{alg_ext, ext_tower, alg_ext_c}` → `giac-poly::Poly<AlgExtC>`  
**快照:** 2026-06-24（**复审** 2026-06-24：T0–T3、P3-6 DoD、solve S0 已落地）

---

## 1. 问题陈述

`giac-poly` 主路径仍为 **`Poly<ℚ>`**（`type PolyQ = Poly<Ratio<BigInt>>`）；**`Poly<AlgExtC>`** 表示层与 K 上 gcd/factor/roots 已接入（P1 + T0–T3）。下列能力仍无法在纯 ℚ 上完成，须扩域闭环：

| 缺口 | 例 | 为何需要扩域 |
|------|-----|--------------|
| partfrac disc>0 实二次分裂 | `partfrac(1/(x²-2), x)` | 分母 `x²-2` 在 ℚ 不可约，须 `A/(x−√2)+B/(x+√2)` |
| factor 二次无理 | `factor(x²-2)` | 根为 `rootof`，非 ℚ 系数 |
| roots 高次 / 双二次 | `roots(t⁴−2=0, t)` | 根落在 ℚ(α) |
| **solve 通用四次** | `solve(t⁴+t+1=0,t)` | 需 K 上 factor/roots + resolvent；非双二次特判 |
| gcd / sturm 代数系数 | `gcd(x²−2, x−√2)` | 系数、余式 ∈ K |

`rootof` 是用户/API 语法；内部一等类型为 **`AlgExt` / `AlgExtC`**（`giac-core`）。目标架构见 [GIAC-algext-adoption](GIAC-algext-adoption.md) §8.4：**`Poly<AlgExtC>`** 上 gcd / factor / roots / partfrac 闭环。

**原则：** 不在现有 `Poly<ℚ>` 上打 `rootof` 补丁；新能力走 **`Poly<AlgExtC>`**（或过渡态仅实代数 `AlgExt` 子集）。

---

## 2. 依赖总览

```text
giac-core 标量闭环（AlgExtC / ExtensionTower::common）
    ↓
giac-poly 表示层（PolyCoeff + Poly<C> + Expr 桥接）
    ↓
最小算法切片（roots → partfrac → factor/gcd）
    ↓
全管线（sturm / 嵌套环 / 跨 crate 接线）
```

**现状快照：**

| 组件 | 状态 |
|------|------|
| `AlgExtData` + 同域 `+−×` | ✅ giac-core Phase A |
| `AlgExtC` / `ExtensionTower::common` | ✅ P0-A/B（塔 T4a 默认） |
| `Poly<AlgExtC>` 表示 + Expr 桥 | ✅ P1-1…P1-5 |
| K 上 gcd / factor / sqff | ✅ P3-1…P3-3（[gcd-factor-priority](../issues_resolved/GIAC-poly-algext-gcd-factor-priority.md) T0–T3） |
| `poly_algext_roots` deg≤4 | ✅ P2-1/2/6、**P3-6 DoD**（[p3-6 gaps](GIAC-poly-p3-6-quartic-roots-gaps.md)） |
| solve / froot 统一管线 | ✅ P4-1/4/6/7（S0；`rootof.rs` 仅测试/参考） |
| partfrac deg≤3 不可约（ℚ） | ✅ |
| partfrac disc>0 二次分裂 | **算法 ✅** + **接线 ✅**（常数/非常数分子）；DIV-084 |
| `∫` 有理 partfrac K 回落 | **✅** P4-3（结构测；`diff(ln rootof)` → B-T3） |
| `field_arith` ↔ `giac-poly` 稠密 poly1 重复 | ☐ [GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md) D1–D4 |

### 2.1 任务优先级（AFK）

**已完成（本 issue 主体）：** P0-A…D、P1-1…5、P2-1…6、P3-1…3、P3-6（DoD）、P3-7、P4-1/4/6/7。详见各 Phase 表 ✅ 列。

| 优先级 | ID | 任务 | 依赖 | 验收 | 估时 |
|:--:|:---:|------|------|------|:--:|
| **P0** | **P4-2** | `eval_partfrac` K 路由：**非常数分子**（如 `x/(x²−2)`） | P2-3/4 ✅ | `expand(partfrac(f,x))≡f`；`partfrac_k_route` 扩例 | ✅ |
| **P0** | **P4-3** | `integrate` 消费上项 partfrac 项 | P4-2 | `∫x/(x²−2)dx` 结构绿；C-02 / `diff` 链待 B-T3 | ✅ |
| **P1** | **P4-5** | `known-divergences` + 同步 [p0-backlog](GIAC-poly-p0-backlog.md) / 本 issue | — | DIV-084/085；partfrac disc>0 追踪 | ✅ |
| **P2** | **F5** | 结构开方 API（`ext_tower` 通用开方收尾） | P3-6 ✅ | [F1-F5](GIAC-poly-quartic-roots-F1-F5.md) §F5 S0–S6 | ✅ S0–S6；F5.5 api-stability 薄化 ⬜ |
| **P2** | **F4′** | Galois σ(κ)：避免四次多余 adjoin | F5 语义稳定 | 维数/塔审计；**不**绑未证 `d_L=12` | 1–2w |
| **P3** | **D3** | 稠密 `poly1` 抽象（`field_arith` 收敛） | — | [dense-poly1](GIAC-dense-poly1-refactor.md) D3 门禁 | ✅ |
| **P4** | **P3-4** | `resultant` / `sturm` over `Poly<AlgExt>`；`realroot` 区间 | P3-1 ✅ | `realroot(x²−2)`；扩 S6 全 Sturm 隔离 | 🟡 resultant + sturm K + realroot 计数 ✅；全 Sturm 隔离 ☐ |
| **P5** | **S7** | `fsolve` 数值互补 | P4-6 ✅ | solve 管线；不替代 `rootof` | 按需 |
| **长期** | **P3-5** | 嵌套环 `Poly<AlgExtC>[main]` | P1-3、FAC | [nested-ring-types](GIAC-poly-nested-ring-types.md) | M4 |

**明确不做（保持登记，勿重开）：**

| 项 | 理由 |
|----|------|
| `Poly<AlgExtC>` 上 Hensel / sparse_bi / unitaryfactor | 上游 `_EXT` 亦不走（§6） |
| 通用五次及以上根式闭式 | Abel–Ruffini；deg≥5 不可约 → `rootof(α)` 一支 |
| `eval(factor(x²−2))` 在 ℚ 上分裂 | 有意偏离 `B-FACTOR-ALGEXT`；K 分裂经 `factor_into_algext` / solve |
| 删 `rootof.rs` 文件 | 非阻塞；生产已走 S0，保留测试/参考 |

**推荐抓取顺序：** P4-2 → P4-3 → P4-5 → F5 → F4′ → D3 → P3-4。

```text
P4-2/3  partfrac 收尾（M2 关门）
    ↓
P4-5    文档同步
    ↓
F5→F4′  四次塔优化（非阻塞主路径）
    ↓
D3      稠密 poly1（T3+ 前）
    ↓
P3-4    sturm/realroot over K（M4）
```

---

## 3. Phase 0 — 前置（giac-core，阻塞 giac-poly）

| ID | 事项 | 说明 | 验收 | 状态 |
|----|------|------|------|------|
| **P0-A** | `AlgExtC` + `canonicalize`（B-06） | 复代数数一等类型；`Complex(0,AlgExt)` / `AlgExt` 统一 | `(i·√2)² = −2`；`canonicalize` 往返 | ✅ |
| **P0-B** | `ExtensionTower::common` 重构（B-05） | `AlgExtData.field` + 塔式/子域 common + cache | `(√2)+(∛2)` 同 compositum | ✅ [T4a/T4b](GIAC-lazy-common-tower-plan.md) |
| **P0-C** | `expr_to_poly` / `poly_to_expr` 边界（B-01） | 含 `AlgExt`：**拒绝**；提升走 **`poly_alg_from_expr`** | 契约见 [expr-poly-conversion.md](../expr-poly-conversion.md) | ✅ |
| **P0-D** | `AlgExt` 同域 `inv` / 跨域 `common`（A-01/A-02 收尾） | partfrac 线性方程组要求域内除法 | `1/rootof(√2)` 不 panic | ✅ T0-1 |

---

## 4. Phase 1 — giac-poly 表示层

| ID | 事项 | 说明 | 验收 |
|----|------|------|------|
| **P1-1** | **系数环 trait** `PolyCoeff` | 抽象 `zero/one/add/mul/div/is_zero`；impl `Ratio<BigInt>`、`AlgExtC` | ✅ `poly_coeff.rs` + `poly_alg_coeff.rs` |
| **P1-2** | **`Poly<C: PolyCoeff>` 骨架** | 从 `Poly<Ratio<BigInt>>` 泛化；`type PolyQ = Poly<Ratio<BigInt>>` 别名 | ✅ `cargo test -p giac-poly` 全绿 |
| **P1-3** | **单变量视图** `UnivariatePoly<C>` / `FlatUni<C>` | 嵌套环 **先仅 ℚ 系数**；AlgExt 版延后 Phase 3 | ✅ 泛型骨架 + ℚ `divides`/`div_rem`；见 [giac-poly-p1-representation.md](../giac-poly-p1-representation.md) §4.3 |
| **P1-4** | **Expr ↔ Poly 桥** `poly_alg_from_expr` / `algext_poly_to_expr` | 新 API，不污染 `expr_to_poly`；系数为 `AlgExtC` | ✅ `giac-core::algebra::poly` |
| **P1-5** | **单项序 / `Var` / `Monomial`** | 变元仍是 `x,y,…`；系数环升级不影响指数向量 | ✅ `Poly::degree_wrt` + [giac-poly-p1-representation.md](../giac-poly-p1-representation.md) §5 |

---

## 5. Phase 2 — 最小算法切片（先解锁 conformance）

| ID | 事项 | 说明 | 验收 | 状态 |
|----|------|------|------|------|
| **P2-1** | **`Poly<AlgExtC>::roots` 二次** | `gausspol` + `alg_ext` | `roots(x²−2,x)` → `[±√2]` | ✅ `poly_algext_roots` |
| **P2-2** | **`Poly<AlgExtC>::roots` 双二次** | 同上 | `t⁴−2=0` 四根 | ✅ |
| **P2-3** | **partfrac：disc>0 实二次分裂** | `sym2poly` partfrac + `_EXT` | `partfrac(1/(x²−2),x)` 两项一次分母 | ✅ 算法；**接线** P4-2 Partial |
| **P2-4** | **partfrac 线性方程组 over K** | 域内 `solve_linear_system` | 系数、右端 ∈ `AlgExtC` | ✅ T2-4 |
| **P2-5** | **`factor` 二次无理（ℚ 上不可约）** | `gausspol` `algext_convert` | `factor(x²−2)` → `(x−α)(x+α)` | ✅ T1-3；删 `try_factor_quadratic_sqrt` |
| **P2-6** | **`Poly<AlgExtC>::roots` 三次** | resolvent / `gausspol` | 一般三次 + `t³−2` 等 | ✅ |

**原推荐顺序（已走完）：** P2-1 → P2-4 → P2-3 → P2-6 → P2-2 → P2-5。剩余见 §2.1 **P4-2/3**。

### partfrac 与 ℚ 路径的关系

| 情形 | 需要 AlgExt？ | 状态 |
|------|---------------|------|
| `x/(x³+2)` 三次不可约 | 否 | ✅ deg≤3 sqff |
| `1/(x²−2)` disc>0 | **是** | ✅ 常数分子；**`x/(x²−2)`** → P4-2 |
| 线性 / 重根 / disc≤0 二次 | 否 | ✅ |
| sqff 因子 deg>3 不可约 | 否（单项式）/ 是（若需分裂） | 单项式 ✅；高次分裂走 K factor（P3-3 ✅） |

---

## 6. Phase 3 — 核心多项式算法（B-02 终态）

**实施顺序（T0→T3、上游 `gcd_ext`/`ext_factor` 对齐）：** [GIAC-poly-algext-gcd-factor-priority](../issues_resolved/GIAC-poly-algext-gcd-factor-priority.md) ✅

| ID | 事项 | 说明 | 验收 | 状态 |
|----|------|------|------|------|
| **P3-1** | **`gcd` / `quo` / `rem` over `Poly<AlgExtC>`** | 子结果式或模 gcd；同扩域 | `gcd(x²−2, x−√2)` | ✅ T1-1 |
| **P3-2** | **`square_free` / `content` / `primitive_part`** | 一元 sqff 先于 factor | partfrac sqff 链在 K 上 | ✅ T1-2 |
| **P3-3** | **`factor` 一元 over K** | 一次/二次分裂/有理根/不可约 witness | `factor(x⁴−4)` | ✅ T2-1（Zassenhaus 子集 **不做**） |
| **P3-4** | **`resultant` / `sturm` over `Poly<AlgExt>`** | 实代数 Sturm（B-04） | `realroot(x²−2)` 区间形式 | 🟡 `resultant_wrt_algext` + `sturm_*` + `realroot` K 计数 ✅；全 Sturm 隔离 ☐ |
| **P3-5** | **嵌套环 `Poly<AlgExtC>[main]`** | `UnivariateIn<C>` 泛化；FAC 管线最后接 | 参系数 + 代数系数塔 | ☐ 长期 |
| **P3-6** | **`Poly<AlgExtC>::roots` 通用四次** | resolvent cubic + K 上二次 split | `solve(t⁴+t+1=0,t)` 四根 `eq_mod` | ✅ DoD；**F4′/F5** 优化见 §2.1 P2 |
| **P3-7** | **`factor` + 有理根降次（deg≥5 前置）** | `factor_into` / sqff 与 solve 共用 | `(x²+1)(x³−x+1)` 五根 | ✅ T2-3 |

**过渡期明确不做：** `Poly<AlgExtC>` 上的 Hensel / sparse_bi / unitaryfactor — 上游 `gausspol` 对 `_EXT` 系数也极受限。  
**明确不做（数学）：** Abel–Ruffini — **无**「通用五次根式闭式」；deg≥5 走 §5.1 兜底策略。

---

## 7. Phase 4 — 跨 crate 接线

| ID | 事项 | crate | 验收 | 状态 |
|----|------|-------|------|------|
| **P4-1** | `giac-solve::roots` 改调 `Poly<AlgExtC>::roots` | giac-solve | 生产走 `poly_algext_roots_for_ctx` | ✅ S0 |
| **P4-2** | `giac-core::eval_partfrac` 含 `AlgExt` 分支 | giac-core + giac-poly | `partfrac(1/(x²−2),x)`；`x/(x²−2)` eval 绿 | ✅ |
| **P4-3** | `giac-calculus::integrate` 消费 AlgExt partfrac 项 | giac-calculus | `∫1/(x²−2)dx`；`∫x/(x²−2)dx` 结构绿 | ✅ |
| **P4-4** | `giac-simplify::factor` 二次 `rootof` 迁入 `giac-poly` | giac-simplify | 删 `try_factor_quadratic_*` | ✅ |
| **P4-5** | 登记 `known-divergences`；更新本 issue 与 p0-backlog | 文档 | DIV-084/085 | ✅ |
| **P4-6** | **`solve` 统一管线：factor 降次 + 递归 `roots`** | giac-solve | deg≤4 → `poly_algext_roots`；deg≥5 → `rootof(α)` | ✅ S0 |
| **P4-7** | **`froot` / `froots` 与 solve 共用因子根** | giac-solve | deg≤4 接 `solve_irreducible_factor` | ✅ S3 |

---

## 5.1 Solve 管线：四次与 deg≥5（normative）

**原则：** 根的类型统一为 **`AlgExt` / `AlgExtC`**（经 `poly_algext_roots`）；生产 solve 已走 S0，**禁止**再增 `rootof.rs` 形状表或 `giac_poly::roots` 旁路。

### 分层策略

```text
solve(P)   P ∈ Poly<ℚ>, 变元 t, deg n
  │
  ├─ sqff / factor 降次（P3-7）
  │     ├─ 有理根 → ℚ
  │     └─ 因子 Pᵢ，deg dᵢ
  │           ├─ dᵢ ≤ 4 → Poly<AlgExtC>::roots（P2-*/P3-6）
  │           └─ dᵢ ≥ 5，Pᵢ 在 ℚ 不可约
  │                 ├─ （可选，后期）Galois 可解子类 → 根式特判
  │                 └─ 默认 → rootof(生成元, Pᵢ) 一支；要小数 → fsolve
  │
  └─ 输出 List<Expr>，元素为 AlgExt / AlgExtC / ℚ（同场根应对齐 field）
```

### 四次（必做）

| 子情形 | 依赖 | 说明 | 状态 |
|--------|------|------|------|
| 二次、双二次 | P2-1、P2-2 | 替代 `quadratic_rootof_roots` / `biquadratic_rootof_roots` | ✅ |
| 三次 | P2-6 | 四次 resolvent 内层；`solve(t³−2=0)` | ✅ |
| 一般四次（含奇次项） | P3-6 + **T3+** | K 上 resolvent + split | ✅ DoD；F4′ 优化 open |
| factor 先降次 | P3-7、P4-6 | `(t²+1)(t²+2)` 等 | ✅ |

### 五次及更高（策略，非根式公式）

| 立场 | 内容 |
|------|------|
| **不做** | 通用五次（及以上）根式闭式 — Abel–Ruffini |
| **要做** | 可约 → factor 递归；不可约 → **`rootof` 精确表示**（`min_poly=Pᵢ` 已在 adoption §2） |
| **可选后期** | 可根式解特判（`x⁵−a`、分圆、soluble Galois）；`sturm`/`realroot` + `rootof` 隔离实根 |
| **数值互补** | `fsolve` 与精确 `rootof` 并存，不替代 |

### 与塔计划（[GIAC-lazy-common-tower-plan](GIAC-lazy-common-tower-plan.md) §11）的阻塞

- **T1–T2：** solve 会话内子域嵌入；多根少 flatten common。✅
- **T3：** 塔顶 parent 系数 `element_*`（T1b 等）。✅
- **T3+：** 四次一般式 register \(u^2-\alpha\) adjoin。✅ 主路径；**F4′** 减少多余 adjoin 仍 open。
- **S0：** 多根列表 `align` / `eq_mod` 不依赖 `embed_a` 调用顺序。✅

### 四次基建完成后的自然增量

统一 `roots` + P4-6 已落地；**待接竖切**见 §2.1：**P4-2/3**（partfrac）、**P3-4**（sturm/realroot）、**S7**（fsolve）。**不**另开 `giac-solve` 补丁线。

---

## 8. 里程碑

| 里程碑 | 范围 | 交付 | 状态 |
|--------|------|------|------|
| **M1** | P0-A/C + P1-1…5 | `Poly<AlgExtC>` 骨架 + Expr 桥接 | ✅ |
| **M2** | P2-1/4/3 + P4-2/3 | `partfrac(1/(x²−2))`；`x/(x²−2)`；∫ K 回落 | ✅ |
| **M3** | P3-1/2/3/6 + P4-1/2/6 | factor/gcd；solve 统一；通用四次 | ✅ |
| **M4** | P3-5 + P3-4 + S7 | 嵌套环；sturm/realroot；fsolve | 🟡 P3-4 resultant + realroot K-Sturm |

```text
M1  Poly<AlgExtC> 骨架                          ✅
M2  partfrac disc>0 + ∫ 对接                    ✅
M3  factor/gcd + solve 统一 + 四次               ✅
M4  嵌套环 + sturm + fsolve                     open
```

---

## 9. 过渡态纪律（避免走弯路）

1. **禁止** 在 `partfrac.rs` 手写 `rootof` 分裂后仍用 `Poly<ℚ>` 解方程 — 数学上不成立。
2. **禁止** 在 `giac-simplify` 再增 `try_partfrac_quadratic_rootof` 形状表 — 归 `Poly<AlgExtC>` 主路径。
3. **可先** 做 **仅实代数** 子集（`AlgExt`，`im=0`）；复根 partfrac 延后到 `AlgExtC`。
4. **assume 不进系数** — `A>0` 选 `rootof` 分支在 `Context`，不在 `coords`（[GIAC-algext-adoption](GIAC-algext-adoption.md) §8.9）。
5. **禁止** 为通用五次及以上写根式闭式 — 见 §5.1；不可约因子用 `rootof(α)` + 可选 `fsolve`。

---

## 10. 与现有 backlog 的关系

| 现有项 | AlgExt 接入后 |
|--------|---------------|
| [GIAC-poly-p0-backlog](GIAC-poly-p0-backlog.md) partfrac disc>0 | P2-3 ✅；**收尾** P4-2/3、P4-5 |
| partfrac deg>3 不可约（ℚ 单项式） | ✅ ℚ 单项式；高次分裂走 K factor（P3-3 ✅） |
| [GIAC-poly-unitaryfactor-gaps](GIAC-poly-unitaryfactor-gaps.md) U5 / sparse_bi | 与 AlgExt **正交**；可并行 |
| e2r / r2e | 未实现；AlgExt 边界见 **P0-C** ✅ |
| [GIAC-poly-p3-6-quartic-roots-gaps](GIAC-poly-p3-6-quartic-roots-gaps.md) | P3-6 DoD ✅；F4′/F5/S7 跟踪 |

---

## 11. 验证

```bash
cd giac-rs

# Phase 0–1（giac-core / giac-poly 表示）
cargo test -p giac-core alg_ext
cargo test -p giac-core alg_ext_c
cargo test -p giac-poly --lib poly_coeff

# K 上 gcd / factor / roots（T0–T3 + P3-6）
cargo test -p giac-core poly_alg_ops
cargo test -p giac-core poly_alg_factor
cargo test -p giac-core poly_roots::tests --release

# partfrac K 路由（P4-2 门禁）
cargo test -p giac-core --test partfrac_k_route
cargo test -p giac-core partfrac

# 跨 crate（P4-*）
cargo test -p giac-solve --release
cargo test -p giac-calculus partfrac_integrate
cargo test -p giac-conformance --test giac_check_integrate  # C-02
```

**门禁：** `cargo test-timeout` 全绿；M2 关门增 `x/(x²−2)` partfrac + integrate `assert_equiv`。

---

## 12. 任务 ID 速查

| Phase | IDs | 未完成 |
|-------|-----|--------|
| 0 前置 | P0-A … P0-D | — |
| 1 表示 | P1-1 … P1-5 | — |
| 2 最小切片 | P2-1 … P2-6 | — |
| 3 核心算法 | P3-1 … P3-7 | **P3-4**、**P3-5**；P3-6 **F4′** |
| 4 接线 | P4-1 … P4-7 | — |

**AFK 抓取顺序：** 见 §2.1。

与 [GIAC-algext-adoption](GIAC-algext-adoption.md) 任务表对应：B-01→P0-C，B-02→P3-*，B-03→P2-*/P4-1，B-04→P3-4，B-05→P0-B，B-06→P0-A，C-02→P2-3/P4-3。

---

## 13. 变更日志

| 日期 | 变更 |
|------|------|
| 2026-06-24 | P4-2/3/5：非常数分子 K partfrac、`expr_to_rational_polys`；DIV-084/085；M2 ✅ |
