# Giac CAS 长远愿景

**状态:** normative（架构约束）  
**相关:** [rust-migration-plan.md](rust-migration-plan.md)（工程迁移路径）、[functional-coverage.md](functional-coverage.md)、[conformance-testing.md](conformance-testing.md)  
**定位:** 本文定义 CAS **数学边界、能力分层与参数化语义**；`rust-migration-plan.md` 定义 **Phase 0–5 工程交付**。二者并列，不互相替代。

---

## 1. 核心结论

**「任意初等函数、系数在 ℚ、任意插入 A/B/C/D、总能给出闭式结果」不能作为可实现的完备目标。**

限制来自数学不可判定性，而非工程能力上限。长远目标应表述为：

> **在 ℚ 系数、任意有限符号参数下，对每个运算给出 {精确闭式, 代数数表示, 级数/渐近, 数值, 或结构化的失败原因}。**

| 运算 | 初等函数 + ℚ 系数 | 现实 |
|------|-------------------|------|
| 加减乘除、有理式化简 | ✅ 可判定 | MVP / Phase A |
| 求导 | ✅ 初等函数的导数仍在初等闭包内 | Phase A（微积分） |
| 积分 | ❌ 不存在通用初等闭式算法 | Liouville + Risch 仅覆盖子集 |
| 化简到 canonical form | ❌ 一般不可判定 | Richardson 定理（含 exp/ln/sin） |
| 极限 | ⚠️ exp-log 子类可判定（Gruntz/MRV） | Phase B；见 [GIAC-limit-exp-difference-unification](issues/GIAC-limit-exp-difference-unification.md) |
| 解方程 | ⚠️ 多项式在 ℚ 可判定；超越方程一般不行 | Phase A–B |

---

## 2. 「在 ℚ 上」的语义分层

同一表述有三层含义，能力承诺不同：

```
L1  系数环 = ℚ
    Poly<Ratio<BigInt>>；normal / gcd / factor / 多项式 solve

L2  表达式系数 ∈ ℚ，含超越函数
    sin(x) + A*exp(B*x)；diff ✅；limit ⚠️；integrate ❌（一般）

L3  值域 / 结果 ∈ ℚ
    ∫ sin(x) dx = -cos(x)     ← 结果不在 ℚ
    ∫ exp(-x²) dx             ← 无初等闭式
```

| 层级 | 含义 | 能力承诺 |
|------|------|----------|
| **L1** | 系数在 ℚ | ✅ 一等能力：`Expr::Rat`、`Poly<Rat>`、`ModInt` |
| **L2** | 符号式系数在 ℚ，可含 exp/ln/trig | ✅ 有理式部分完备；超越部分按 §4 分层 |
| **L3** | 结果必须是 ℚ 中的数或式 | ❌ 对积分、极限、三角化简通常不成立；允许 `AlgExt`、`rootof`、`Unevaluated` |

**Normative：** 对外文档与 API 契约以 **L1 + L2** 为准；**不**承诺 L3。

---

## 3. 与 rust-migration-plan 的关系

```
rust-migration-plan.md          cas-long-term-vision.md（本文）
─────────────────────          ─────────────────────────────
Phase 0–5（~7 月）              Phase A（= 迁移计划 Phase 0–5）
53 golden / WASM headless       Phase B（初等函数引擎，+12–18 月）
Expr / crate 结构               Phase C（参数化 CAS，+6–12 月）
测试即规格                      Phase D（研究级超越，无 deadline）
```

| 文档 | 回答的问题 |
|------|------------|
| [rust-migration-plan.md](rust-migration-plan.md) | 何时交付、测什么、crate 怎么拆 |
| 本文 | 数学上能承诺什么、参数 A/B/C/D 语义、失败时返回什么 |

**Phase A 完成定义**（不变）：见 `rust-migration-plan.md` §10——53 项 conformance 全绿、不依赖 C++ giac、`cargo test` + `cargo ci-clippy` 持续全绿。

---

## 4. 初等函数闭包与运算承诺

### 4.1 闭包定义

可承诺的初等函数闭包 **E**：

```
E = 有理式( ℚ[x₁,…,xₙ, A, B, C, D, …] )
  ∪ { exp, ln, sin, cos, tan, … } 的有限复合
  ∪ rootof / AlgExt
```

### 4.2 对 f ∈ E 的分运算承诺

| 运算 | 承诺 | 失败语义 |
|------|------|----------|
| `+ − × ÷`, `expand`, `normal`（有理式） | ✅ 完备 | — |
| `gcd`, `factor`, `solve`（多项式） | ✅ 完备（显式变量列表） | `DivisionByZero` 等 |
| `diff(f, x)` | ✅ 结果在 E 中（或 Err） | 参数 A,B,… 视为常数 |
| `integrate(f, x)` | ⚠️ 在 E 中 **或** `NotElementary` | Risch 可判定子集；其余诚实失败 |
| `limit(f, x→a)` | ⚠️ exp-log：Gruntz/MRV；其余：级数/启发式 | 见 §6 极限架构 |
| `normal(f)`（含超越） | ⚠️ 启发式 + `assert_equiv` 验证 | 不承诺全局零判定 / 最简形 |

### 4.3 不可判定项（禁止对外承诺）

- 任意初等表达式的全局 canonical 化简
- 任意初等函数的初等原函数
- 任意初等表达式的恒等判定（零判定）
- 任意超越方程的闭式解

---

## 5. 长远阶段规划

### Phase A — 基线（= rust-migration-plan Phase 0–5）

**目标：** 可替换 headless giac CAS；WASM 同语义。

**范围：** 见 [rust-migration-plan.md](rust-migration-plan.md) §1.1（P0–P14）、§5（Phase 0–5）。

**核心类型（与迁移计划 §3.2 对齐）：**

- `Expr`：`Int` / `Rat` / `Frac` / `Symbol` / `Add` / `Mul` / `Pow` / `Func` / `AlgExt` / …
- `Poly<Rat>`、`AlgExtData`、`Context`、`ModInt`
- 符号层与数值层分离；`Arc<Expr>` 不可变 + 结构共享

**验收：** 53 项 golden；`assert_equiv` 处理等价不同形（见 [conformance-testing.md](conformance-testing.md)）。

---

### Phase B — 初等函数引擎（+12–18 月，估算）

**目标：** 在 exp-log 与三角子类上接近 Maxima/Sage 可靠度；**非**「任意初等完备」。

| 模块 | 目标 | 依赖 issue / 参考 |
|------|------|-------------------|
| `giac-simplify` | 统一 canonical 策略 + 规则 DAG | `assert_equiv` 为门禁 |
| `giac-calculus::integrate` | 规则表 → Hermite/Rothstein → **真 Risch 子集** | GIAC-229–231；[phase4-issues.md](phase4-issues.md) |
| `giac-calculus::limit` | Gruntz/MRV **单管线** | [GIAC-limit-exp-difference-unification](issues/GIAC-limit-exp-difference-unification.md)、[GIAC-216e](issues/GIAC-216e-mrv-series-lead-convergence.md) |
| `giac-solve` | 多项式 / 超越方程分离路由 | 超越部分 → `Err` 或 `RootOf` |
| `giac-core::alg_ext` | 代数数闭包：+ − × ÷、minpoly、`rootof` | `testcas` / geo |

**极限架构约束（normative，摘自 limit issue P1–P3）：**

1. **同一理论，两层规则：** `factor_exp_shifted_difference`（x 层）与 `remove_lnexp`（w 层）收敛为共享规则表，不得各维护独立特例。
2. **快路径有保质期：** 临时 `limit_factored_exp_growth` 类快路径在 MRV 主路径收敛后删除或仅留 `#[cfg(test)]`。
3. **禁止 per-gruntz 形状表：** 新 gruntz 用例驱动 `remove_lnexp` / preprocess，不得新增 `limit_*_at_infinity` 按用例命名函数。

**积分路径（目标管线）：**

```
integrate 请求
  → 规则表（代数函数、三角、有理式换元等）
  → partfrac / Hermite / Rothstein
  → risch() 兜底（真 Risch 子集）
  → NotElementary | Exact(expr)
```

**验收（Phase B 增量）：**

- Maxima gruntz rtest 扩展集（目标 ≥ 当前 giac-rs 对齐集全绿）
- `giac_check_integrate` 全绿或 documented 等价偏离
- 无新增 per-gruntz 极限形状函数（P3）

---

### Phase C — 参数化 CAS（+6–12 月，估算）

**目标：** **A、B、C、D 及任意 Ident 作为一等符号参数**，运算语义明确、可测试。

#### 5.1 符号参数模型

```rust
// 概念层 — 扩展 rust-migration-plan §3.2 Context

pub enum SymbolRole {
    Variable,   // x, t — 微积分 / 极限的活跃变量
    Parameter,  // A, B, C, D — 默认不参与 groebner 主元选择等
    Constant,   // π, e — 已知超越常数
}

pub struct Context {
    pub vars: HashMap<Ident, Arc<Expr>>,
    pub assumptions: Vec<Assumption>,       // assume(A, integer), assume(B>0)
    pub symbol_roles: HashMap<Ident, SymbolRole>,
    pub active_var: Option<Ident>,          // diff / integrate / limit 默认变量
}
```

#### 5.2 参数化能力矩阵

| 操作 | 对任意参数 A,B,C,D | 说明 |
|------|-------------------|------|
| 解析 / AST | ✅ | 无硬编码变量表；任意 `Ident` |
| 算术、`expand`、`normal`（有理式） | ✅ | 系数环 ℚ(A,B,…) |
| `gcd` / `factor` / `solve`（多项式） | ✅ | 变量列表显式传入 `lvar` / `e2r` |
| `diff` | ✅ | 参数视为常数 |
| `subst` / `:=` / 程序 | ✅ | Phase A `giac-prog` |
| `integrate` | ⚠️ | 参数保留在结果中；积分本身仍受 §4.2 限制 |
| `limit` | ⚠️ | 单分支：`assume` 选路；全分支：`piecewise`（§5.3）；缺信息 → `NeedsAssume` |
| `solve`（含参） | ⚠️ | 如 `solve(A*x=B,x)` → `x=B/A`，需 `A≠0` 或 `Failed(NeedsAssume)` 或 `piecewise` |

#### 5.3 assume、piecewise 与参数分情形

参数 A,B,C 与 `assume` / `piecewise` 解决**不同层次**的问题；与 [GIAC-algext-adoption.md](issues/GIAC-algext-adoption.md) §8.9 对齐。

| 机制 | 存放位置 | 回答的问题 | Phase |
|------|----------|------------|-------|
| **`SymbolRole::Parameter`** | `Context::symbol_roles` | A,B 是参数还是变量？ | C |
| **`assume`** | `Context::assumptions` | **当前会话**里，对参数**额外假定什么**？ | C |
| **`piecewise`** | `Expr` AST（结果式子内） | **所有参数区域**上的完整答案是什么？ | D |
| **`NeedsAssume`** | `CasResult::Failed` | 缺假设、无法安全选支时**停住** | C |

**A,B 本身不是假设**：未写 `assume` 时参数仍可在式中出现，算法须保守（见下）。

##### 5.3.1 `assume` — 会话级分支选择

`assume` 把**已知但未写进式子**的信息注入 `Context`，供各算法在**分叉处选一支**继续算。

**概念模型（扩展 `Context::assumptions`）：**

```rust
pub enum Assumption {
    Integer(Ident),
    Real(Ident),
    Positive(Ident),
    // Phase C 扩展（概念层）：
    // NonZero(Ident),
    // Relation(ExprArc),   // A>0, A≠B, A>B, …
}

// 语句级：assume(A>0); …; purge(A)  — 见 GIAC-204b
```

**Normative（Phase C）：**

- 无 `assume` 时：**保守**——不将 `sqrt(A^2)` 化为 `A`，不默认 `A≠0`，不在 MRV 增长比较中猜 `A>B`。
- 有 `assume(A, integer)` / `assume(B>0)` / `assume(A≠0)` 时：化简、积分、求解、极限可走对应分支。
- 需要假设但未给出 → `CasResult::Failed(NeedsAssume(...))`，**非静默错答**。
- **`purge(A)`** 清除与 A 相关的绑定与假设；同式子可能从「可算」变回 `NeedsAssume`（Context 状态变化，非数学矛盾）。
- **assume 不进 `AlgExt` 坐标 / 不进 `min_poly`**——假设只在 `Context`，不在代数扩域系数里编码「A>0」。

**典型用途：**

| 场景 | assume 作用 |
|------|-------------|
| `normal(sqrt(A^2))` | `assume(A>0)` → `A` |
| `solve(A*x+B=0, x)` | `assume(A≠0)` → `x = -B/A` |
| `integrate(abs(A*x), x)` | `assume(A>0)` → 去绝对值分支 |
| `limit(..., x, +inf)` 含 `exp(A*x)-exp(B*x)` | `assume(A>B)` → MRV 选 `exp(A*x)` 为 dominant，走单支 Gruntz/MRV |

**与 limit / MRV 的关系：** `exp_diff`、`mrv_w` 的代数 fold（`exp(A)-exp(B)`、`ln(w)` padd）**不读** `assumptions`；`assume` 仅在**下游决策点**介入——MRV 比较谁更快、选 `scale_log` / ε 是否 → 0、lead 系数约分。详见 [exp-diff-expr-api.md](exp-diff-expr-api.md)、[limit-engine-expr-api.md](limit-engine-expr-api.md)。

##### 5.3.2 `piecewise` — 定义与功能

**定义：** `piecewise` 是 **E 内的一等分段表达式**，用**显式条件列表**表示「参数或自变量落在不同区域时，值不同」。它是**结果的一部分**（写在 `Expr` 里），不是 `Context` 里的临时假设。

**语法（对标 upstream giac / Xcas）：**

```text
piecewise(c₁, e₁, c₂, e₂, …, cₙ, eₙ [, e_default])
when(c, e_then, e_else)          // 二元分支；嵌套 when ≡ 多支 piecewise
```

- 偶数个主体参数时：`c₁,e₁,c₂,e₂,…` — 依次检验条件，**第一个为真**的 `eᵢ` 为值。
- 奇数个且最后一项无条件：最后一项为 **default**（所有前面条件均为假时）。
- 条件 `cᵢ` 为逻辑式（等式、不等式、`and`/`or`）；Phase D 前可限制为参数上的多项式不等式 + 等式。

**概念 AST（Phase D 目标）：**

```rust
/// 分段值：按顺序匹配，首个条件为真者生效；optional default 兜底
pub struct Piecewise {
    pub branches: Vec<(ExprArc /* cond */, ExprArc /* value */)>,
    pub default: Option<ExprArc>,
}

// 纳入 Expr 枚举，或 FuncKind::Piecewise；打印 round-trip 为 piecewise(...)
```

**功能（为何需要 piecewise）：**

1. **完整回答含参问题** — 参数在分界面（如 `A=B`）两侧极限/解不同；`piecewise` 一次给出**全部分支**，而非只算 assume 下的一支。
2. **参数空间的「中断点」** — 如 \(\lim_{x\to+\infty}(e^{Ax}-e^{Bx})/e^{Bx}\) 在 `A=B` 与 `A≠B` 处行为不同；数学上应分情形，不应强行合并为一个无参式。
3. **替代静默错答** — 当算法能枚举有限分支但用户未 `assume` 时，优先 `Exact(piecewise(...))`，而非猜一支或返回错误数值。
4. **可组合** — `diff` / `integrate` / `limit` / `subst` 对 `piecewise` 按分支传播（各支独立运算，条件不变或按规则合并）；`eval` / 数值代入在条件可判定后落单支。

**算法保证：** 语法不隐含「分类完备」；分支生成的可判定性与承诺范围见 **§5.3.5**（L0–L3）。

**与 `assume` 对比：**

| | `assume` | `piecewise` |
|--|----------|-------------|
| 存哪 | `Context`（会话） | `Expr`（答案） |
| 语义 | 「**在此假设下**继续算」 | 「**在所有区域上**答案如此」 |
| 分支数 | 隐含选 **1** 支 | 显式列出 **多** 支 |
| 生命周期 | `purge` 可清 | 持久于结果式 |
| Phase | C（语句级 + 算法读 Context） | D（AST + 分支传播） |

##### 5.3.3 决策流：assume / NeedsAssume / piecewise

含参运算遇到**无法比较**或**多分支**时的 normative 路由：

```text
含参请求（limit / solve / normal / …）
  │
  ├─ Context 已有 assume，且足以唯一选支
  │     → 走单支算法 → CasResult::Exact(expr)     // 参数可仍留在式中
  │
  ├─ 无 assume，但条件可枚举为有限分支（Phase D，**L1**）
  │     → 各支分别计算 → CasResult::Exact(piecewise(...))
  │
  ├─ 无 assume，需用户指定才安全（Phase C）
  │     → CasResult::Failed(NeedsAssume { hint: "A≠0" | "compare A and B" | … })
  │
  └─ 已知不可判定 / 未实现
        → Failed(Undecidable | NotImplemented)
```

**极限示例（概念）：**

```text
limit((exp(A*x) - exp(B*x)) / exp(B*x), x, +infinity)

assume(A > B)  →  Exact(+infinity)           // 或经 MRV 的等价形
assume(A < B)  →  Exact(0)
assume(A = B)  →  Exact(0)                   // ε 小量路径

无 assume（Phase D 目标）→
  Exact(piecewise(
    A > B,  +infinity,
    A < B,  0,
    A = B,  0
  ))
```

Phase C 在 `piecewise` AST 未就绪时：无 `assume(A>B)` 等 → **`NeedsAssume`**，不静默选 `exp(A*x)` 或 `exp(B*x)` 为主导项。

##### 5.3.4 验收

**Phase C（assume）：**

- 参数化 golden：`integrate(A*sin(x),x)`、`solve(A*x+B=0,x)`、`normal((A*x)^2/A)` **分 assume 有无**
- 与 `test_subst`、`assume` 在 flanex / testintegrate 中的行为一致或 documented 偏离
- 语句级 `assume` / `purge`：**GIAC-204b**

**Phase D（piecewise，增量）：**

- `piecewise` / `when` 解析与打印 round-trip
- `diff` / `subst` 对简单 `piecewise` 分支传播
- 含参 `limit` / `solve` 在分界面用 `piecewise` 覆盖（替代仅 `NeedsAssume`）
- 与 upstream `prog.cc` `piecewise` / `when` 语义一致或登记 [known-divergences.md](known-divergences.md)
- 分支生成仅在其 **L 层级**（§5.3.5）承诺范围内声称「完备」；超越该层 → `NeedsAssume` / `Undecidable`

##### 5.3.5 分支生成的算法保证边界

**Normative：** `piecewise` 语法是**答案表示**；「分类讨论是否正确、是否完备」取决于**谁生成** `(cᵢ, eᵢ)` 以及问题落在哪一层。**禁止**对外暗示「启用 piecewise 即自动保证含参分类完备」。

**分层（L-piecewise）：**

| 层级 | 分支来源 | 典型问题 | 算法基础 | 保证（实数语义下） | giac-rs Phase |
|------|----------|----------|----------|-------------------|---------------|
| **L0** | 用户 / 程序**手写** `piecewise` | 任意（条件由用户负责） | 无自动生成；`eval`/`subst` 按首真条件选支 | **表示语义**正确：条件为真时取对应支；**不**保证用户列全分支 | D（AST） |
| **L1** | 算法**有限枚举**显式比较 | `A≶B`、`A=0`；含参 `limit` 中 MRV 主导项三分 | 各支调用原算法（Gruntz/MRV、solve、normal） | **各支在对应条件下**结果正确；**不**保证找全所有分界面 | C–D |
| **L2** | **实代数**胞腔分解 | 系数 ∈ ℚ(A,B,…)；条件为多项式等式/不等式；含参 `solve`、符号 `realroot` | Sturm / 结果ants / **CAD** / 实闭域 **QE**（子集） | 在**多项式假设**下，胞腔**互斥且覆盖**参数空间（给定变量序）；各胞内分支固定 | D+（Sturm 已有；CAD/QE 后置） |
| **L3** | 含 **exp/ln/trig** 的含参极限、积分、化简 | `limit(..., A,B)`、`integrate(abs(A*x),x)` 无 assume | Gruntz/MRV、启发式；**无**通用完备分类器 | **不承诺**自动 `piecewise` 完备；缺信息 → `NeedsAssume`；已知不可判定 → `Undecidable`（§4.3） | B–C |

**正确性含义（按层）：**

- **L0：** 契约在 `piecewise`/`when` 的 **AST 语义**与分支传播规则；与用户手写 `assume` 等价于只算一支。
- **L1：** `Exact(piecewise(...))` 正确 ⟺ 每个 `(cᵢ,eᵢ)` 满足「在 `cᵢ` 为真的参数区域内，`eᵢ` 等于该运算的数学结果」；**不要求** `⋃cᵢ` 覆盖整个参数空间。
- **L2：** 额外要求分支条件来自 **CAD/QE/Sturm** 输出，开胞内公式 sign-invariant；登记实现所支持的 **多项式次数/变量数** 上限。
- **L3：** 仅承诺 **不静默错答**（§5.3.3）；可输出**已知有限** L1 分支（如 `A>B|A<B|A=B`），但**不**声称分界面完备。

**与 limit / `exp_diff` / MRV：** L1/L3 的分叉点在 MRV 增长比较、主导 `exp` 选择、lead 约分；**不在** `canonical_exp_diff` / `canonical_mrv_coeff` fold 层（见 [exp-diff-expr-api.md](exp-diff-expr-api.md)）。

**upstream giac 对照（`to_piecewise`，非完备引擎）：**

giac `piecewise(expr [, x])` 在 `prog.cc` 中可对式子调用 `to_piecewise(e, x)`（`signalprocessing.cc`）。算法概要：

```text
abs/sign → Heaviside → 收集 Heaviside 线性变元的断点 zᵢ
  → 若断点经 evalf 均为 numeric：在 (zᵢ, zᵢ₊₁) 上 interval 化简 → 拼 piecewise(x < zᵢ, …)
  → 若仍含符号断点或 Heaviside：放弃，原式返回
```

| 项 | giac `to_piecewise` | giac-rs 目标（本文 L 层） |
|----|---------------------|---------------------------|
| 适用对象 | 单变量 `x` 上含 `abs`/Heaviside 的式子 | L0 任意；L1 参数比较；L2 多项式参数 |
| 断点 | **数值化**后区间分割 | L2：符号 CAD/QE；L1：显式 `A≶B` |
| 含参 `A,B` | 断点非 numeric 则**失败退回** | L2 专门处理；L3 不冒充 L2 |
| `piecewise` 求值 | `prog.cc`：运行时首真条件选支 | 同 L0 语义 |
| 含参 `limit` | 不自动生成 parametric piecewise | L1 手工枚举 + L3 `NeedsAssume` |

**giac-rs 对外表述（normative）：**

1. 文档与 API **必须标注**结果所属的 L 层（或 `Failed` 原因）。
2. **L2 未实现前**，含参 `solve`/`limit` 不得输出声称「参数空间完备」的 `piecewise`。
3. **L3** 默认路由：`assume` → 单支；否则 `NeedsAssume`；仅当算法**显式识别**有限 L1 分支集时才 `Exact(piecewise(...))`。
4. 引入 CAD/QE 时单独 milestone + 测试集，不与 Phase D 基础 `piecewise` AST 混为一谈。

**参考实现：** giac-2.0.0 `signalprocessing.cc`（`to_piecewise`、`flatten_piecewise`）；`prog.cc`（`_piecewise` 求值）。

---

### Phase D — 可选超越（研究级，不设 deadline）

不阻塞 Phase A 交付；独立 milestone + 独立测试集。

| 方向 | 说明 |
|------|------|
| 完整 Risch（超越情形） | 超出 Phase B 子集 |
| 更广初等极限 | Gruntz 全套件；SymPy/Maxima rtest 对标 |
| `gbasis` | 当前迁移计划后置；`greduce` 先行 |
| `piecewise` | §5.3.2–5.3.3：AST、分支传播、含参 limit/solve |
| 条件化简 | 依赖完整 assume 引擎（§5.3.1） |

---

## 6. 输出契约：`CasResult`

长远 API 不应只有 `Expr` 或 panic；**五种合法输出**：

```rust
pub enum CasResult {
    /// 精确符号结果（参数可保留；Phase D 可含 `piecewise` 节点，见 §5.3.2）
    Exact(Arc<Expr>),
    /// 代数数（rootof / AlgExt）
    Algebraic(Arc<AlgExtData>),
    /// 级数 / 渐近主项（极限、级数展开）
    Series { lead: Arc<Expr>, rest: Option<Arc<Expr>> },
    /// 数值（evalf / fsolve）
    Numeric(Numeric),
    /// 结构化失败
    Failed {
        kind: FailKind,
        partial: Option<Arc<Expr>>,
        hint: String,
    },
}

pub enum FailKind {
    NotElementary,    // 无初等闭式（积分等）
    NotImplemented,   // 算法未实现
    NeedsAssume,      // 缺少 assume（如 A≠0、A>B）；见 §5.3.3
    Undecidable,      // 已知不可判定
    Timeout,          // 资源上限
}
```

**Normative：**

- Phase A 内部可仍用 `Result<Arc<Expr>, EvalError>`；Phase B 起新 API 与 WASM 边界优先 `CasResult`。
- WASM 层（`eval_to_string`）对 `Failed` 须输出可解析、可测试的错误前缀，而非空串或挂起。

---

## 7. 内核架构（目标态）

```
                    ┌─────────────────────────────────┐
                    │  giac-parse → Stmt / Expr       │
                    └───────────────┬─────────────────┘
                                    │
                    ┌───────────────▼─────────────────┐
                    │  Context（vars, assume, roles）   │
                    └───────────────┬─────────────────┘
                                    │
         ┌──────────────────────────┼──────────────────────────┐
         │                          │                          │
         ▼                          ▼                          ▼
  ┌─────────────┐          ┌─────────────┐          ┌─────────────┐
  │ 多项式快路径  │          │ 符号 AST     │          │ 数值快路径   │
  │ Poly<Rat>   │◄─e2r/r2e─►│ Arc<Expr>   │─evalf───►│ Numeric     │
  │ PolyMod     │          │ AlgExt      │          │ Matrix<f64> │
  └─────────────┘          └──────┬──────┘          └─────────────┘
                                  │
              ┌───────────────────┼───────────────────┐
              ▼                   ▼                   ▼
        giac-simplify        giac-calculus        giac-solve
        normal/trig          diff/limit/          solve/sturm
                             integrate/series
                                  │
                                  ▼
                            CasResult
```

**设计要点（与迁移计划 §3.4 一致并延伸）：**

1. 符号层与数值层分离。
2. 多项式快路径不进用户-facing `Expr` 枚举；`to_poly()` 桥接。
3. 代数数 `AlgExt` 是一等公民，不可推迟。
4. 参数符号与活跃变量由 `Context` 显式管理，非全局隐式状态。
5. 算法可偏离 giac 实现；验证靠数学等价 + 测试 + [known-divergences.md](known-divergences.md)。

---

## 8. 直接 FAQ

### Q1：CAS 最终要写成什么样？

- **内核：** `Expr` + `Poly<Rat>` + `AlgExt` + 分层 `Numeric`
- **运算：** 按 {多项式可判定, 微积分子集, 启发式+验证, 数值} 四层路由
- **输出：** 精确符号优先；失败时 `CasResult::Failed`
- **部署：** Native CLI + WASM 同语义（见 [rust-migration-plan.md](rust-migration-plan.md) §WASM）
- **参数：** A/B/C/D 与 x 同为 `Ident`；`Context` 区分 `SymbolRole` 与 `Assumption`

### Q2：能否对任意初等函数在 ℚ 上给出结果，并支持任意 A,B,C,D？

| 子问题 | 答案 |
|--------|------|
| A,B,C,D 任意插入 | ✅ Phase C 明确目标 |
| 系数在 ℚ | ✅ L1/L2 一等能力 |
| 任意初等、总有闭式 | ❌ 数学不可能 |
| 可替代承诺 | ✅ 最大可判定子集 + 诚实失败（§4、§6） |

### Q3：与 giac C++ 的关系？

Phase A：**测试即规格**，golden 等价或 `assert_equiv` + documented 偏离。  
Phase B+：**数学语义** 为权威；giac 为参考实现，非算法行级标准（见 [rust-migration-plan.md](rust-migration-plan.md) §6.4）。

---

## 9. 里程碑摘要

| 里程碑 | 标志 | 与迁移计划关系 |
|--------|------|----------------|
| **M-A** | 53 conformance + WASM | = 迁移 M4 |
| **M-B1** | Gruntz/MRV 单管线；limit issue Phase C 验收 | Phase B |
| **M-B2** | 真 Risch 子集；integrate check 全绿或等价 | Phase B |
| **M-C** | 参数化 golden + `CasResult` WASM 边界 | Phase C |
| **M-D** | 独立 research milestone（可选） | Phase D |

---

## 10. 文档维护

| 变更类型 | 更新本文 | 更新 rust-migration-plan |
|----------|----------|--------------------------|
| 新增数学能力承诺 | ✅ | 若需新 crate/phase 则 ✅ |
| 工程 phase 工期调整 | 仅交叉引用 | ✅ |
| 已知不可判定边界澄清 | ✅ | — |
| golden / crate 结构 | 交叉引用 | ✅ |

---

## 参考

- [rust-migration-plan.md](rust-migration-plan.md) — 工程迁移主方案
- [conformance-testing.md](conformance-testing.md) — `assert_equiv` 规格
- [external-test-resources.md](external-test-resources.md) — 外部 fixture 与抽取脚本
- [phase4-issues.md](phase4-issues.md) — 微积分 / Risch 移植 issue 树
- [issues/GIAC-algext-adoption.md](issues/GIAC-algext-adoption.md) — 参数 A,B 与 assume（§8.9）
- [exp-diff-expr-api.md](exp-diff-expr-api.md)、[limit-engine-expr-api.md](limit-engine-expr-api.md) — limit MRV / exp 差分（assume 不介入 fold 层）
- [issues/GIAC-limit-exp-difference-unification.md](issues/GIAC-limit-exp-difference-unification.md) — 极限统一管线
- Richardson, *Some Undecidable Problems Involving Elementary Functions*
- Gruntz, *On Computing Limits in a Symbolic Manipulation World*
- Liouville / Risch — 初等积分可判定性
