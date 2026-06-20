# giac-poly — `unitaryfactor` / `pzadic` / P2a 数学原理

**类型:** 算法原理（FAC-G1 尾部）  
**上游基线:** `giac/giac-2.0.0` `gausspol.cc` — `unitaryfactor` L6701、`pzadic` L3893、`unitarize` L6783  
**实现:** `giac-rs/crates/giac-poly/src/factor/unitary.rs`  
**缺口索引:** [GIAC-poly-unitaryfactor-gaps](issues/GIAC-poly-unitaryfactor-gaps.md)（实现状态 / 优先级）  
**环类型:** [GIAC-poly-nested-ring-types](issues/GIAC-poly-nested-ring-types.md)（`UnivariateIn::divides` 等）  
**快照日期:** 2026-06-19

---

## 1. 问题陈述

多元多项式 `p ∈ ℚ[x₁,…,xₙ]` 在 sparse / Hensel 等路径失败后，上游 `do_factor_hensel` 调用 **`unitaryfactor`** 作有界启发式兜底。

核心子问题：**赋值降维 + 系数抬升**。

1. 选赋值变量 `eval_var` 与大整数 `B`
2. 代入得 `p(B, ·) ∈ ℚ[main]`（其余变量已消去或进入递归）
3. 在 `main` 上分解赋值像，得一元因子 `f ∈ ℚ[main]`
4. 把 `f` **抬回** `ℚ[eval_var][main]` 中的真正因子 `f̂`
5. 用精确整除验证 `f̂ | p`，剥因子后递归

第 4 步是数学难点：赋值后系数是**常数**，真因子系数是 **`eval_var` 的多项式**。

---

## 2. 环与变量约定（giac-rs）

| 符号 | 含义 |
|------|------|
| `vars_rev` | 上游反序变量表；`eval_var = vars_rev[0]`，`main = vars_rev.last()` |
| `p` | 当前 sqff 块，`p ∈ ℚ[eval_var, …, main]`（嵌套环） |
| `ev` | `p` 在 `eval_var ↦ B` 后的像，`ev ∈ ℚ[main]`（系数为常数） |
| `fz` | `ev` 在 `ℚ[main]` 上的因子列表 |
| 抬升目标 | `f̂ ∈ ℚ[eval_var][main]`，使得 `f̂(B, x) = f(x)` 且 `f̂ | p` |

**整除语义：** 抬升后的 peel 必须用 `UnivariateIn::divides`（`ℚ[eval_var][main]` 上一元精确除），**禁止** `Poly::div_rem`（多元 leading-term 除法）。

---

## 3. 外层赋值轨迹（upstream）与 GCDHEU 背景

giac-rs **不做** 全局 `2..N` 小整数扫描。外层仅在 upstream 轨迹上换点；该轨迹与 giac 上游 `unitaryfactor`、`gcdheu` 共用同一套赋值点启发式，源自 **GCDHEU**（启发式多项式 GCD）。

### 3.1 公式（giac-rs / 测试对齐）

| 步骤 | 公式 | 实现 |
|------|------|------|
| 初值 | `x₀ = 2·‖p‖∞ + 2`（有理系数时再加 `\|denom‖∞`） | `UnitaryEvalPoint::initial` |
| sqff 微调 | 若 `p(x₀)` 对 `main` 非 sqff → `x₀ += 1`（有界） | `bump_sqff`，上限 `UNITARY_MAX_TRY = 48` |
| 失败换点 | `x₀ ← ⌊x₀ · 73794 / 27011⌋ + 1` | `UnitaryEvalPoint::advance` |
| stream 终止 | `x₀` 位长 `> 256` bit → `EvalBaseStream::next()` 返回 false | 与外层 `ntry > 48` 取先触发者 |

**line 25 轨迹前几项**（`‖p‖∞` 小，初值 `38`）：

```text
38 → 104 → 285 → 779 → 2129 → 5817 → …
```

每步约乘以 `73794/27011 ≈ 2.732`（giac-rs 在整数商后再 `+1`）。

**与上游 `unitaryfactor` 的细微差：** `gausspol.cc` L6761 为 `x0 = iquo(x0·73794, 27011)`（**无** `+1`）；giac-rs `advance` **含** `+1`，测试 `eval_base_stream_upstream_only` 按后者验收。`gcdheu`（L4703）同样无 `+1`。

### 3.2 数学原理：GCDHEU 赋值–重构

算法族来自 **Char–Geddes–Gonnet** 的启发式 GCD（**GCDHEU**）：

1. 选大整数赋值点 `n`
2. 计算 `a(n)`、`b(n)`（或多元情形下递归代入）
3. 算整数 GCD `g(n) = igcd(a(n), b(n))`
4. 用 **对称模 `n` 数位展开** 把 `g(n)` 重构为多项式 `g*(x)`（Maple 内核称 `genpoly`；giac 称 **`pzadic`**）
5. **试除验证**：`g*` 须同时整除 `a`、`b`；否则说明 `n` 太小或遭遇 spurious factor，增大 `n` 重试

这是典型的 **Las Vegas** 启发式：循环可能很久，但一旦试除通过，结果正确。

**`pzadic` 与 GCDHEU 的关系：** `pzadic(f, B)` 正是第 4 步的 n-adic 抬升——把赋值后常数系数 `c_e` 按 base `B` 展开为 `C_e(y)` 的多项式 digit，见 §4。`unitaryfactor` 不单独发明新的 advance 理论，而是把 GCDHEU 的「换 `n` 再试」借来做多元分解的外层赋值搜索。

### 3.3 初值 `2·‖p‖∞ + 2` 的含义

设真 GCD（或真因子）系数绝对值上界约为 `α`。要保证单点 digit 重构不混淆相邻系数，需要赋值点 `n` 足够大，使 `n` 大于重构多项式系数的两倍量级（Gonnet 讲义中的下界论证；最坏情形与 `deg`、`α` 相关）。

实践中不先解出精确下界，而取与输入规模相关的**保守初值**：

```text
n₀ = 2 · max_coeff_magnitude(p) + 2
```

与 Maple `gcdheu` / giac `listmax` 一致。多元时 `unitaryfactor` 对当前 sqff 块 `unitaryp` 取 `‖·‖∞`。

### 3.4 步进因子 `73794/27011` 的来源

| 出处 | 失败后的 `n` 更新 |
|------|-------------------|
| Char–Geddes–Gonnet 伪代码 / Gonnet 讲义 | `n := 2·n + 1` |
| **Maple 内核 `gcdheu`** | `n := ⌊n · 73794 / 27011⌋` |
| giac `gausspol.cc`（`gcdheu`、`unitaryfactor`） | 同 Maple（`iquo`） |
| giac-rs `advance` | `⌊n · 73794 / 27011⌋ + 1` |

`73794/27011` **不是** Char–Geddes–Gonnet 原文公式，而是 **Maple 对 GCDHEU 的工程常数**（约在 `2` 与 `3` 之间的有理倍率，使 `n` 较快增大又避免 `2n+1` 式过快膨胀）。有据可查的文献描述见 **Liao & Fateman**, *Evaluation of the Heuristic Polynomial GCD*, ISSAC 1995（伪代码 `β ← β * 73794 / 27011`）。

giac 从 Maple 实现链路继承该常数；`unitaryfactor` 与 `gcdheu` 共用，故分解外层轨迹与 GCD 启发式同源。

### 3.5 有界性与代价

- **外层 peel 轮数：** `ntry ≤ UNITARY_MAX_TRY`（48）→ 最多 48 次 `advance` 级换点（生产路径常在较早基成功即停，如 line 25 约 `104`）。
- **位长熔断：** `EvalBaseStream::next()` 要求 `x₀.bits() ≤ 256`；典型小范数输入先触达 48 轮上限，而非 bit 上限。
- **P2a 局部窗**（§5）在**固定外层 `B`** 上另扫 `b ∈ [B-(need-1), …]`，不计入 upstream `advance` 轨迹；大 `B`（如 `779`）时代价主要在赋值后大整数一元分解，见测试 `p2a_line25_upstream_trajectory` 性能分析。

sqff 微调发生在**每一轮**外层尝试的内层；`EvalBaseStream` 只负责 `initial` 与 `advance`。

---

## 4. 主路径：`pzadic`（单点 digit 抬升）

### 4.1 数学

设赋值因子 `f(x) = Σ_{e=0}^n c_e x^e`，其中每个 `c_e ∈ ℚ` 为常数（已代入 `eval_var = B`）。

真因子形如：

```text
f̂(x, y) = Σ_{e=0}^n C_e(y) · x^e
```

其中 `y` 即 `eval_var`，`C_e(y) ∈ ℚ[y]`，且 `C_e(B) = c_e`。

**`pzadic` 的做法：** 对每个常数 `c_e`，把有理数 `c_e` 按 base `B` 做**对称数位展开**（upstream `smod`），用 `y^j` 编码第 `j` 位：

```text
c_e = d₀ + d₁·B + d₂·B² + …   →   C_e(y) = d₀ + d₁·y + d₂·y² + …
```

（实现：`PzadicLift::pzadic`，`giac-rs/crates/giac-poly/src/factor/unitary.rs`）

### 4.2 适用条件

- 每个 `C_e(y)` 的次数较低，且单点 base-`B` 展开能忠实表示（如 `C_e(y) = y+3`、纯常数等）
- 与上游 `pzadic` 使 `dim+1`、在反序变量下展开 digit 的语义一致

### 4.3 失败情形

当 `C_e(y)` 次数较高或 digit 展开在**当前** `B` 上不能代表真多项式时，单点 `pzadic` 给出错误 `f̂`，`f̂ ∤ p`。

**典型：** line 25 第二因子含 `y³`；在 `B = 38` 上 slot-1 的 pzadic 抬升不整除（slot-0 往往仍可 pzadic 成功）。

---

## 5. 增强路径：P2a（局部窗 + monic + Lagrange）

当 `pzadic` 失败时，`try_lift_and_peel` 调用 `lift_factor_multi_eval`（P2a）。

### 5.1 核心恒等式

将真因子写成：

```text
f̂(x, y) = Σ_{e=0}^n C_e(y) · x^e ,   C_e(y) ∈ ℚ[y]
```

在赋值点 `y = b_i` 上，赋值因子（未归一前）满足：

```text
f(b_i, x) = Σ_e C_e(b_i) · x^e
```

若在**多个** distinct 的 `b_i` 上取得**同一语义槽位**的因子，则 `(b_i, coeff_e(f(b_i,x)))` 是 `C_e(y)` 的样本 → 可插值恢复 `C_e(y)`。

### 5.2 三步流程

```mermaid
flowchart TD
    A["pzadic(f, B) 失败"] --> B["局部窗采样 b ∈ [B-(need-1), …]"]
    B --> C["y↦b，sqff 微调，factor_univariate(ev)"]
    C --> D["sort 槽位 + monic 归一（x 首项系数=1）"]
    D --> E["记录 (b, coeff_x^e) 样本"]
    E --> F{样本数 ≥ need?}
    F -->|否| B
    F -->|是| G["对每个 e：Lagrange → C_e(y)"]
    G --> H["组装 f̂ = Σ C_e(y)x^e + x^n"]
    H --> I{f̂ | p ?}
    I -->|是| J["peel"]
    I -->|否| K["放弃该槽位/该 B"]
```

#### （1）局部窗（非全局扫描）

- 外层 `B` 来自 upstream `initial` / `advance`
- P2a **仅在 pzadic 失败时**启动，在 `B` 附近采样：

```text
b ∈ [ B - (need-1), B - (need-2), … , B + sqff微调 ]
```

其中 `need = min(deg_y(p)+1, MULTI_EVAL_MAX_SAMPLES)`（默认上限 8）。

**line 25 例：** `B=38`，`need=6` → 窗含 `33…38`，覆盖有效邻域点 `34`，但**不**扫全局 `2..N`。

#### （2）monic 归一

在 `ℚ[x]` 上，因子相差非零常数倍。不同 `b_i` 上 `factor_univariate` 的 leading coefficient 可能不同，直接对系数插值会错位。

**做法：** 对每个样本因子 `g`，除以 `x` 的首项系数，使最高次项系数恒为 `1`：

```text
g_monic(x) = g(x) / lc_x(g)
```

只记录 `e < deg_x` 的系数样本；最高次项固定为 `x^n`（系数 `1`），不插值。

#### （3）Lagrange 插值

对每个次数 `e`，样本 `(b_0,v_0),…,(b_{k-1},v_{k-1})`：

```text
C_e(y) = Σ_i v_i · L_i(y) ,   L_i(y) = Π_{j≠i} (y - b_j)/(b_i - b_j)
```

在 `ℚ[y]` 上精确构造（有理系数算术）。组装：

```text
f̂(x,y) = Σ_{e=0}^{n-1} C_e(y)·x^e + x^n
```

### 5.3 槽位对齐

多采样点须对应**同一个**多元因子。giac-rs 对 `fz` 按 `(deg_x, leading_coeff)` 排序，用固定 `factor_slot` 索引；各 `b_i` 上 `deg_x` 不一致的样本丢弃。

### 5.4 验证闭环

插值只产生**候选**；必须 `UnivariateIn(f̂).divides(p)` 才接受。失败则换槽位或换外层 `B`。

**注意：** 生产路径是**顺序 peel**（剥一个因子后更新 `unitaryp` 再继续），不是要求单点 batch 一次剥光所有槽位。

---

## 6. 与上游 giac 的对应

| 环节 | 上游 `gausspol.cc` | giac-rs |
|------|-------------------|---------|
| 外层赋值 | `x0=2np+2`，sqff bump，`x0←iquo(x0·73794,27011)`（giac-rs：`+1`） | `UnitaryEvalPoint` + `EvalBaseStream` |
| 抬升 | `pzadic(*f_it, x0)` | `PzadicLift::pzadic` |
| 失败补救 | 主要靠递归 / 多轮剥离 | **显式 P2a**（局部窗 + monic + Lagrange） |
| 非 monic 首项 | `unitarize` → `unitaryfactor` → `ununitarize` | `unitarize` / `ununitarize` 第二路径 |
| `main` 次数归零 | `trunc1` + 递归 | `factor_constant_tail` |

P2a 是 giac-rs 在保持 upstream 外层轨迹前提下的**增强 lift**，不是替代 `unitaryfactor` 的独立分解算法。

---

## 7. 完整管线位置（FAC-G1）

```text
factor_sqff_over_coeff_ring_ctx
  → try_sparse_factor / try_sparse_factor_bi
  → try_hensel_lift_bivariate
  → try_unitary_factor → unitary_factor_rev
       ├─ pzadic peel（主）
       └─ P2a lift_factor_multi_eval（fallback）
```

---

## 8. 回归锚点（line 25）

合成多项式（Hensel / sparse 均失败，unitary 门禁）：

```text
f1 = 3*x - y^2 + y - 5
f2 = x*y + 3*x - y^2 - 1 + y^3
p  = f1 * f2
```

| 路径 | 结果 |
|------|------|
| `try_sparse_factor` | None |
| `try_hensel_lift_bivariate` | None |
| `try_unitary_factor` | 2 因子（P2a + 顺序 peel） |
| 测试 | `unitary_factor_line25_l22_y3`、`testfactor_line25_unitaryfactor_gate` |

---

## 9. 实现索引

| 函数 / 类型 | 文件 | 角色 |
|-------------|------|------|
| `unitary_factor_rev` | `factor/unitary.rs` | 主循环：赋值 → 分解 → peel |
| `PzadicLift::pzadic` | 同上 | 单点 digit 抬升 |
| `lift_factor_multi_eval` | 同上 | P2a 局部窗 + 插值 |
| `try_lift_and_peel` | 同上 | pzadic 优先，P2a fallback |
| `lagrange_interp_coeff` | 同上 | `ℚ[y]` 插值核 |
| `monic_wrt_main` | 同上 | 跨采样点系数对齐 |
| `EvalBaseStream` | 同上 | upstream-only 外层轨迹 |

---

## 10. 参考

### 代码

- 上游：`giac/giac-2.0.0/src/gausspol.cc` — `pzadic` L3893、`unitaryfactor` L6701、`gcdheu` L4662、`unitarize` L6783
- giac-rs：`giac-rs/crates/giac-poly/src/factor/unitary.rs` — `UnitaryEvalPoint`、`EvalBaseStream`、`PzadicLift`

### 文献（赋值轨迹 / GCDHEU）

| 文献 | 说明 |
|------|------|
| **Char, Geddes, Gonnet** — *GCDHEU: Heuristic Polynomial GCD Algorithm based on Integer GCD Computation* | EUROSAM'84 (LNCS 174); 正式版 **J. Symbolic Computation 7(1):31–48, 1989** — 赋值、整数 GCD、n-adic 重构、试除验证 |
| **Geddes, Czapor, Labahn** — *Algorithms for Computer Algebra* (Kluwer, 1992) | 教材第 6 章 GCD；GCDHEU 系统化叙述 |
| **Gonnet** — [Heuristic Algorithms](https://people.inf.ethz.ch/gonnet/CAII/HeuristicAlgorithms/node1.html) (1999) | 在线讲义；伪代码用 `n := 2·n+1` 换点 |
| **Liao & Fateman** — *Evaluation of the Heuristic Polynomial GCD*, ISSAC 1995 | Maple `gcdheu` 实现剖析；**`73794/27011`** 步进常数 |
| **Kaltofen** — 对多元 GCDHEU 正确性的补充证明（arXiv cs/0206032 等） | 多元情形赋值假设的修补 |

### 项目文档

- 缺口与优先级：[GIAC-poly-unitaryfactor-gaps](issues/GIAC-poly-unitaryfactor-gaps.md)
- API 分层：[giac-poly-api-stability.md](giac-poly-api-stability.md) §FAC-G1
- 嵌套环除法：[GIAC-poly-nested-ring-types](issues/GIAC-poly-nested-ring-types.md)
