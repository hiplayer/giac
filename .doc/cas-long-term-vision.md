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
| `limit` | ⚠️ | 可能需要 `assume` 或分情形（`piecewise` 后期） |
| `solve`（含参） | ⚠️ | 如 `solve(A*x=B,x)` → `x=B/A`，需 `A≠0` 或 `Failed(NeedsAssume)` |

#### 5.3 assume 与保守化简

**Normative：**

- 无 `assume` 时：**保守**——不将 `sqrt(A^2)` 化为 `A`，不假设 `A≠0`。
- 有 `assume(A, integer)` / `assume(B>0)` 时：化简、积分、求解可走对应分支。
- 需要假设但未给出 → `CasResult::Failed(NeedsAssume(...))`，非静默错答。

**验收（Phase C）：**

- 参数化 golden 集：`integrate(A*sin(x),x)`、`solve(A*x+B=0,x)`、`normal((A*x)^2/A)` 分 assume 有无
- 与 `test_subst`、`assume` 在 flanex / testintegrate 中的行为一致或 documented 偏离

---

### Phase D — 可选超越（研究级，不设 deadline）

不阻塞 Phase A 交付；独立 milestone + 独立测试集。

| 方向 | 说明 |
|------|------|
| 完整 Risch（超越情形） | 超出 Phase B 子集 |
| 更广初等极限 | Gruntz 全套件；SymPy/Maxima rtest 对标 |
| `gbasis` | 当前迁移计划后置；`greduce` 先行 |
| `piecewise` | 参数分情形、`limit` / `solve` 分支 |
| 条件化简 | 依赖完整 assume 引擎 |

---

## 6. 输出契约：`CasResult`

长远 API 不应只有 `Expr` 或 panic；**五种合法输出**：

```rust
pub enum CasResult {
    /// 精确符号结果（参数可保留在式中）
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
    NeedsAssume,      // 缺少 assume（如 A≠0）
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
- [issues/GIAC-limit-exp-difference-unification.md](issues/GIAC-limit-exp-difference-unification.md) — 极限统一管线
- Richardson, *Some Undecidable Problems Involving Elementary Functions*
- Gruntz, *On Computing Limits in a Symbolic Manipulation World*
- Liouville / Risch — 初等积分可判定性
