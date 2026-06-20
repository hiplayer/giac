# GIAC-AlgExt — `AlgExtData` / `Expr::AlgExt` 落地与算法适配清单

**状态:** open（Phase A 骨架已合入 giac-rs）  
**类型:** AFK（实现向）  
**相关:** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md)（giac-poly 分阶段待办）、[rust-migration-plan.md](../rust-migration-plan.md) §3.2、`cas-long-term-vision.md` §5 Phase B、`GIAC-205`（`rootof` 浅层对接）  
**上游参考:** `giac/giac-1.5.0/src/alg_ext.cc`、`alg_ext.h`、`gen.h` `ref_algext`  
**Rust 落点:** `giac-core::algebra::{alg_ext, ext_tower, alg_ext_c}`、`Expr::AlgExt` / `Expr::AlgExtC`、`giac-poly::Poly<AlgExtC>`

---

## 问题陈述

giac-rs 原先仅用 `Func(RootOf, …)` 表示代数数，**没有**扩域算术类型。与 upstream `_EXT` 不对齐，导致：

- `solve` / `roots` 无法对 `rootof` 结果做 `+ − × ÷`
- `normal` / `factor` / `sturm` / `integrate` 遇到代数系数时无法继续
- `testcas`、`testgeo`、`flanex` 大量 `rootof` 用例只能符号打印，不能参与运算

本 issue 跟踪：**已落地的最小 `AlgExt` 内核** + **各 crate 算法适配优先级**。

---

## 数学原理

本节记录 `AlgExt` / `rootof` 的数学模型、表示约定、设计动机，以及与本工程算法的对应关系（normative 背景，供实现与 review 参照）。

### 1. 适用范围：代数数，不是全体实数

从有理数 **ℚ** 出发，添加一个（或有限个）满足 **有理系数多项式** 的数 α，得到代数扩域 **ℚ(α) ⊆ ℂ**（实代数数落在 **ℝ ∩ ℚ(α)**）。

```text
P(α) = 0，P ∈ ℚ[t]，P 不可约（或取最小多项式）
ℚ(α) = { f(α) | f ∈ ℚ[t] }，维数 = deg(P)
```

**能精确表示的：** √2、∛2、黄金比 `(1+√5)/2`、某五次方程的特定根等 **代数数**。

**不能表示的：** π、e、ln 2、sin(1) 等 **超越数**（不满足任何非零有理系数多项式）。它们需 `Func(Sin/Exp/Ln, …)`、级数或 `Numeric` 浮点层，**不属于 `AlgExt`**。

与 [cas-long-term-vision.md](../cas-long-term-vision.md) **L3** 一致：精确 CAS 允许结果超出纯 ℚ 有理式，但仍在 **代数闭包** 或 **未求值符号** 内；超越结果走 `Unevaluated` / 数值层。

### 2. 数据结构 ↔ 数学对象

`AlgExtData` 编码扩域 **ℚ(α) 中的一个元素**：

| 字段 | 数学含义 |
|------|----------|
| `min_poly` | α 的 **最小多项式** `P`（giac `poly1`，**高次系数在前**） |
| `coords` | 元素在基 `{1, α, α², …, α^{d−1}}` 下的坐标：`c₀ + c₁α + …`（同 poly1 顺序） |
| `root_index` | 可选：`P` 有多个根时，选定哪一支（对标 `select_root` / `proot`） |

**恒等式：**

```text
元素值 = coords 所代表的多项式在 α 处取值
约束   = min_poly(α) = 0
运算   = 对 coords 做多项式 +、×，再 mod min_poly（`ext_add` / `ext_mul`）
```

**与 `rootof(num, minpoly)` 的对应：**

```text
rootof([c_{d−1}, …, c₀], poly1[P])  ↔  min_poly = P，coords = [c_{d−1}, …, c₀]
```

示例（giac-rs 当前 MVP，系数 ∈ ℚ）：

| 数 | `min_poly`（高→低） | `coords` | 说明 |
|----|---------------------|----------|------|
| **√2** | `[1, 0, -2]`（t²−2） | `[1, 0]` | α=√2，元素即 α |
| **−√2** | 同上 | `[-1, 0]` | 同域另一根 |
| **∛2** | `[1, 0, 0, -2]`（t³−2） | `[1, 0, 0]` | β=∛2 |

### 3. 同域与跨域：`common_EXT`

**同域**（相同 `min_poly`）：`ext_add` / `ext_mul` 直接在 `coords` 上 mod `min_poly` 即可（giac-rs Phase A 已实现）。

**跨域**（不同 `min_poly`）：元素分别落在 ℚ(√2) 与 ℚ(∛2)，**不能直接相加**。须先 **合并扩域**（upstream `common_EXT` / `common_minimal_POLY`）：

```text
√2  ∈ ℚ(√2)          min_poly: t²−2
∛2  ∈ ℚ(∛2)          min_poly: u³−2
√2+∛2 ∈ ℚ(√2,∛2)     公共 min_poly: M(γ)=0，deg(M) ≤ 6
```

**√2 + ∛2 的精确表示** 不是两个 `rootof` 并列，而是 **一个** `AlgExt`：

1. `common_EXT` 构造包含 √2、∛2 的 ℚ(γ)，γ 的最小多项式 **M** 在 ℚ 上次数为 **6**（经典结果）：

   ```text
   M(t) = t⁶ − 6t⁴ − 4t³ + 12t² − 24t − 8
   poly1: [1, 0, -6, -4, 12, -24, -8]
   ```

2. 将 √2、∛2 **嵌入** 新基 `{1, γ, …, γ⁵}`，再 `ext_add`，得到 6 个有理坐标 `coords`（非直观的 `[1,1]`）。

giac-rs **已实现** B-05（T4a 塔 compositum 默认、`subfield_common_pair`、cache）；跨域 `add` 经 `align_elements` → `common_over_q`。

### 4. 为何用「ℚ 多项式 + ℚ 系数」表示

| 好处 | 说明 |
|------|------|
| **精确** | 无浮点误差；`√2·√2` 恒为 `2` |
| **封闭性** | 代数数对 `+ − × ÷`（分母非零）封闭；结果仍可用 `(min_poly, coords)` 表示 |
| **可判定相等** | 同域内 `coords` mod `min_poly` 相同 ⟺ 数学相等 |
| **与 `Poly<Rat>` 统一** | 系数环从 ℚ 推广到 ℚ(α)；gcd / factor / resultant 可扩到无理系数 |
| **对齐 upstream** | Giac `_EXT`、`testcas`、`testgeo` 均建立在此模型上 |

**与仅保留 `Func(RootOf, …)` 符号的区别：**

| | `RootOf` 符号 | `Expr::AlgExt` + 运算 |
|--|---------------|------------------------|
| 打印 / 解析 | ✅ | ✅（可转回 `rootof`） |
| `+ − ×` | ❌ | ✅（同域；跨域需 `common_EXT`） |
| `/`、`normal` 约分 | ❌ | 需 `inv_EXT`（A-02） |
| 多项式 gcd / factor | ❌ | 需 **`Poly<AlgExtC>`**（§8.4） |
| 浮点近似 | ❌ | 需 `alg_evalf` → `Numeric`（D-03），**不**扩 `Expr` 变体 |

### 5. 三层架构（精确 vs 数值）

```text
Expr::AlgExtC    — 精确复代数数（目标标量，§8.3）
Poly<AlgExtC>    — 系数在代数闭包子集上的多项式算法（§8.4）
Numeric / evalf  — AlgExtC::evalf → F64，见 cas-long-term-vision CasResult::Numeric
```

`evalf` **不应** 把结果存回 `AlgExt`；精确与数值分层，避免污染符号层。

### 6. 对本工程哪些算法有用

**直接依赖（无 AlgExt 则无法精确完成或需登记 DIV）：**

| 算法 | 用途 | 典型式子 |
|------|------|----------|
| `solve` / `roots` | 多项式无理根 | `solve(t²−2=0,t)` |
| `factor` | 不可约二次/三次因子 | `factor(x²−2)` |
| `normal` / `ratnormal` | 含根式分式化简 | `(√2+1)/(√2−1)` |
| `partfrac` / `integrate` | 对数项参数为代数数 | `∫1/(x²−2)dx` |
| `sturm` / `realroot` | 代数多项式实根隔离 | `realroot(x²−2)` |
| `resultant` / gcd（扩域） | 含根式多项式 | 几何消元 |

**间接依赖：** 精确几何（`testgeo`）、`charpoly` / `egv`（特征值为代数数）、扩域上 Groebner。

**不依赖 AlgExt：** `sin`/`exp`/`ln`、一般 transcendental 极限/积分（走 calculus / MRV / Risch 等其它管线）。

### 7. MVP 与过渡态限制的数学含义

| 限制 | 数学含义 | 解除方向（见 §8 目标架构） |
|------|----------|---------------------------|
| 系数仅 **ℚ**（Int/Rat） | `coords`、`min_poly` 系数不能含符号 **A,B,x** | 阶段 3：`PolyCoeff` + 参数字段；与 `assume` 结合见 **§8.9** |
| **`i` 为符号而非生成元** | `Complex(re, im)` 中 `im·i` 的 `i` 不在任何 `min_poly` 内 | 阶段 2b：`AlgExtC` 或塔层 `Adj(i, t²+1)` |
| **双形态** `AlgExt` / `Complex(0, AlgExt)` | 同一复代数数两种 `Expr` 形状，算法需分支 | `canonicalize → AlgExtC` |
| **跨域 re/im** | 实部、虚部可能落在不同 ℚ(α) | `ExtensionTower::common` |
| **`common_EXT` MVP** | 仅 primitive element 试探，未缓存、未对齐 upstream 全算法 | 阶段 2：`ExtensionField` + 缓存 |
| **无 `evalf`** | 不能从精确根式得到小数近似 | 阶段 4：`AlgExtC::evalf` → `Numeric` |

### 8. 目标架构：塔式 `_EXT` + `Poly<AlgExtC>`（normative）

本节取代原「阶段 1–5 平铺表」，作为 **Phase B 完成后的 normative 数据模型**。当前 giac-rs 处于 **过渡态**（§8.1）；最终形态对齐 upstream `ref_algext` 的 **塔式扩域** + 系数环 **`Poly<AlgExtC>`**。

#### 8.1 三层与过渡态（当前 → 目标）

```text
                    ┌─────────────────────────────────────────┐
  标量（根的值）     │  AlgExtC  z = re + im·i，re,im ∈ K      │
                    │  K = 塔顶域 ExtensionTower              │
                    └─────────────────┬───────────────────────┘
                                      │ 系数环
                    ┌─────────────────▼───────────────────────┐
  多项式算法层       │  Poly<AlgExtC>   gcd / factor / roots   │
                    └─────────────────┬───────────────────────┘
                                      │ 求值 / 显示
                    ┌─────────────────▼───────────────────────┐
  Expr 层（用户式）  │  AlgExtC | AlgExt(im=0) | rootof(...) │
                    └─────────────────────────────────────────┘

过渡态（已部分落地）：
  Expr::AlgExt          — 实代数数 ℚ(α)
  Expr::Complex(re,im)  — re,im 可为 AlgExt；i 仍为 Context 符号
  fold_complex_algext_* — 在 eval 中 patch，非一等类型

目标态：
  Expr::AlgExtC(Arc<AlgExtCData>)  — 复代数数唯一标量形态
  canonicalize: Complex/AlgExt → AlgExtC；im=0 可降回 AlgExt 显示
```

#### 8.2 塔式 `ExtensionTower`（对标 upstream `_EXT`）

upstream 每次 `rootof` / 开方实质是 **在基域上再 adjoin 一个代数元**。Rust 用显式塔记录每次扩张及嵌入：

```rust
/// 扩域塔：K₀=ℚ → K₁=ℚ(α₁) → … → Kₙ（塔顶）
pub enum ExtensionTower {
    /// K = ℚ
    Base,
    /// K = parent(γ)，γ 满足 min_poly(γ)=0
    Adj {
        parent: Arc<ExtensionTower>,
        /// 新生成元最小多项式，系数 ∈ parent 的可 rationalize 部分（Phase C 前 ∈ ℚ）
        min_poly: Vec<ExprArc>,
        /// parent 的标准基向量嵌入塔顶基：dim(parent) 个长度 dim(top) 的坐标列
        /// embed_parent[j] = parent 的第 j 个基向量在塔顶下的 coords
        embed_parent: Vec<Vec<ExprArc>>,
        degree: usize, // = deg(min_poly)
    },
}

/// 塔顶域的共享句柄（Arc 比较 + 可选 canonical hash）
pub struct ExtensionField {
    pub tower: Arc<ExtensionTower>,
    /// 预计算：塔顶乘法表（coords × coords → coords mod min_poly），可选
    mul_table: Option<Arc<MulTable>>,
}
```

**运算语义：**

| 操作 | 行为 |
|------|------|
| `tower.adj(P)` | 若 `P` 在 parent 上不可约 → 新塔层；否则可能降次或报错 |
| `tower.common(a, b)` | 两子塔的最小公共超塔；对标 `common_EXT` / `common_minimal_POLY` |
| `embed(x: Kᵢ, target: Kⱼ)` | 沿塔嵌入矩阵把 coords 提升到 j≥i 的基 |
| 元素 `+ − × ÷` | 先 `common` 到同一 `ExtensionField`，再 mod 塔顶 `min_poly` |

**Lazy `common`（Phase S0/S1，[GIAC-lazy-common-tower-plan.md](GIAC-lazy-common-tower-plan.md)）：**

| 层 | 何时 | `common`? |
|----|------|-----------|
| L0 `Expr` / display | 解析、未 eval、`to_rootof_expr` | 否 |
| L1 `AlgExtData` | `from_rootof`、同域 `+ − ×` | 否（最小域） |
| L2 对齐后 | 跨域 `add`/`mul`/`eq_mod`、经 [`align_elements`](../../giac-rs/crates/giac-core/src/algebra/ext_tower.rs) | 运算时可能 `common` |

决策树（同域 / 子域嵌入 / `common_cache`）见塔计划 **§12.5**；`fold_algext_sum` Split 语义见 **§12.6 E**。S0 已落地：`embedding_for` + 统一 `align_elements` 入口。

**示例塔（`solve(t⁴−2=0)` 自然生长）：**

```text
K₀ = ℚ
K₁ = ℚ(α)     min_poly: t²−2          α = √2
K₂ = K₁(β)    min_poly: u²−α（或 u⁴−2 一次 adj） β = 2^(1/4)
K₃ = K₂(i)    min_poly: v²+1           i 为形式生成元

±2^(1/4)           ∈ AlgExtC { field: K₂, re: ±β, im: 0 }
±i·2^(1/4)         ∈ AlgExtC { field: K₃, re: 0, im: β }   // β 嵌入 K₃
```

对比过渡态：`algext_sqrt_branches` 对 `u<0` 返回 `Complex(0, AlgExt)`，**未**把 `i` 写入塔；目标态开方失败时 **`tower.adj(u²+1)` 或 `tower.adj(i)`** 再返回 `AlgExtC`。

#### 8.3 `AlgExtC`：复代数标量

```rust
/// z ∈ K[i] / (i²+1)，K = ExtensionField 的塔顶
pub struct AlgExtCData {
    pub field: Arc<ExtensionField>,
    /// re, im ∈ K，在塔顶基下各一组 coords（长度 = field.degree）
    pub re: Vec<ExprArc>,
    pub im: Vec<ExprArc>,
    pub root_index: Option<u32>,
}

impl AlgExtCData {
    // (a+bi) ± (c+di)
    fn add(&self, other: &Self) -> Result<Self, EvalError>;
    // (a+bi)(c+di) = (ac−bd) + (ad+bc)i，ac,bd,... 为 K 内 ext_mul
    fn mul(&self, other: &Self) -> Result<Self, EvalError>;
    // z⁻¹ = (a−bi) / N(z)，N(z)=a²+b² ∈ K
    fn inv(&self) -> Result<Self, EvalError>;
    fn eq_mod(&self, other: &Self) -> Result<bool, EvalError>;
    fn evalf(&self, prec: u32) -> Result<(f64, f64), EvalError>; // → Numeric，不进 Expr
}
```

**与 `Expr` 的关系：**

| 形态 | 目标 |
|------|------|
| `Expr::AlgExt(a)` | `AlgExtC { re: a.coords, im: 0, field: a.field }` |
| `Expr::Complex(re, im)` | `AlgExtC { re, im, field: common(re.field, im.field) }` |
| 显示 | `im=0` → `rootof(...)`；否则 `(re)+(im)*i` 或 `rootof(...)*i` |
| `eval` | 优先产出 `AlgExtC`；兼容期可降 `AlgExt` / `Complex` |

#### 8.4 `Poly<AlgExtC>`：多项式算法闭环

系数环从 `Poly<Rat>` 提升到 `Poly<AlgExtC>` 后，B-02 / B-03 / B-04 走 **同一条管线**：

```rust
// giac-poly（或 giac-core::algebra）扩展
pub struct Poly<C> { /* 现有稀疏/ dense 结构 */ }

impl Poly<AlgExtC> {
    fn gcd(&self, other: &Self) -> Self;
    fn factor(&self) -> Vec<(Self, u32)>;
    fn roots(&self, var: &Var) -> Vec<AlgExtCData>;
    // Sturm / realroot：系数在 K 内，序列在 Poly<AlgExt> 或 Poly<AlgExtC> 上
}
```

**`solve` 目标流程（取代 biquadratic 特判）：**

```text
1. expr → Poly<Rat>（或 Poly<Var>）
2. factor / square-free 分解
3. 对每个不可约因子 f：
     若 deg(f) 可根式 / rootof → Poly<AlgExtC>::roots(f)
     实根隔离 → realroot 用 Poly<AlgExt>（im=0 子集）
4. 根列表 Vec<AlgExtC> → Expr::List（canonicalize 显示）
```

| 式子 | Poly<Rat> | Poly<AlgExtC> roots |
|------|-----------|---------------------|
| `t²−2=0` | 一次因子 | ±√2 |
| `t⁴−2=0` | 双二次 | ±2^(1/4), ±i·2^(1/4)（4 个 AlgExtC） |
| `x²−2` factor | 不可约二次 | rootof 因子 |

#### 8.5 精确 vs 数值（不变）

```text
Expr::AlgExtC    — 精确复代数数（本 issue 目标标量）
Poly<AlgExtC>    — 精确多项式算法
Numeric / evalf  — AlgExtC::evalf → (re, im) f64，不进 Expr
```

#### 8.6 中期 vs 长期（Horizon）

阶段编号（§8.7）描述**交付顺序**；本节用 **中期 / 长期** 描述**目标层次**，避免把「类型重构」与「多项式管线」混为一谈。

| | **中期** | **长期** |
|---|----------|----------|
| **核心问题** | 复代数数**怎么存、怎么算** | 含代数根的多项式**怎么 factor / gcd / roots** |
| **交付重心** | `ExtensionTower` + `ExtensionField` + **`AlgExtC`** | **`Poly<AlgExtC>`** + `evalf` + 参数扩域 |
| **对应阶段** | **2**（塔 + common）、**2b**（AlgExtC + canonicalize） | **5**（Poly 算法）、**4**（evalf）、**3**（PolyCoeff） |
| **完成后** | 标量 `+−×÷` 闭环；`i` 进塔；跨域 `common` 可缓存 | 通用 `solve` / `factor` / sturm；无 biquadratic 等特判 |
| **仍可保留** | `Poly<AlgExtC>::roots` 仅二次/双二次；`solve` 部分特判 | — |
| **验收** | `AlgExtC` 四则与 `canonicalize`；`(i·α)²=−α²` | `testcas` / `testgeo` 子集；`Poly<AlgExtC>::factor` |

```text
现在（过渡态）          中期                    长期
AlgExt + Complex patch  ExtensionTower          Poly<AlgExtC>
biquadratic 特判    →   AlgExtC 标量闭环    →   gcd/factor/roots 通用
i 为符号                i 进塔 / K[i]           evalf；参数 A,B
```

**判断标准：**

- **中期完成** = 类型设计正确、标量运算正确；算法可以「够用但不通用」。
- **长期完成** = 多项式算法在 `Poly<AlgExtC>` 上闭环；过渡态 patch 可删除。

#### 8.7 实施阶段（交付顺序）

| 阶段 | Horizon | 交付物 | 能力 | 与过渡态关系 |
|------|---------|--------|------|--------------|
| **A** ✅ | 现在 | `AlgExtData`, `ext_*`, `inv`, P0/P1 MVP | 实代数数同域运算；`Complex+AlgExt` patch | 当前主分支 |
| **2** | 中期 | `ExtensionField` + `ExtensionTower::Adj` + `common` 缓存 | 跨域合并、塔顶统一；`AlgExt.min_poly` → `field` | 重构 `AlgExtData` 字段 |
| **2b** | 中期 | `AlgExtCData`, `Expr::AlgExtC`, `canonicalize` | 复代数数 `+−×÷`；`i` 进塔或 K[i]/(i²+1) | 替代 `fold_complex_algext_*` |
| **3** | 长期 | `PolyCoeff` / 参数 min_poly | 符号参数 A,B | Phase C |
| **4** | 长期 | `AlgExtC::evalf` | 浮点近似 | D-03 |
| **5** | 长期 | `Poly<AlgExtC>` gcd/factor/roots | 通用 solve/factor/sturm | B-02/B-03/B-04 终态 |

**giac-poly 分任务清单：** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md)（Phase 0–4、P2-3 partfrac、里程碑 M1–M4）。

**推荐 PR 顺序（自阶段 2 起）：**

```text
ExtensionTower + ExtensionField::common（重构 B-05）
  → AlgExtData.field 迁移（兼容 min_poly 访问器）
  → AlgExtC + canonicalize + inv
  → giac-poly Poly<AlgExtC> 骨架 + roots 二次/双二次
  → 删除 biquadratic_rootof_roots / algext_sqrt_branches 特判
  → Poly<AlgExtC> factor/gcd
  → evalf
```

#### 8.8 与 upstream C++ 的对应

| upstream | giac-rs 目标 |
|----------|--------------|
| `ref_algext` / `_EXT` | `AlgExtCData` + `ExtensionField` |
| `common_EXT` | `ExtensionField::common` |
| `ext_add` / `ext_mul` / `inv_EXT` | K 内 coords 运算；`AlgExtC` 复乘法 |
| `horner_rootof` / `proot` | `AlgExtC::evalf` |
| `gausspol` + `_EXT` 系数 | `Poly<AlgExtC>` |
| `complex_mode` + `_EXT` | 统一 `AlgExtC`，不再分支 `ComplexVal` 与 `AlgExt` |
| `sym2poly.cc::check_assume` | `Context::assumptions` + `check_assume` 在 e2r / factor / solve 前 |

#### 8.9 参数 A,B 与 `assume`（阶段 3 / Phase C）

阶段 3（`PolyCoeff`、参数 `min_poly`）与 [cas-long-term-vision.md](../cas-long-term-vision.md) §5.2–5.3 的 **Phase C 参数化** 共用同一套 `Context` 约定。参数 **A,B** 与 **`assume`** 解决不同问题，经 `Context` 一并传给各算法。

**分工：**

| 机制 | 回答的问题 | 典型内容 |
|------|------------|----------|
| **`SymbolRole::Parameter`**（A,B,C,D） | 这个符号**是什么角色**？ | 不是积分/求导的活跃变量 |
| **`assume`** | 对这个符号**额外知道什么**？ | `A>0`、`A` 为整数、`A≠0` |
| **`PolyCoeff` / 参数字段**（阶段 3） | 系数环**能含什么**？ | `min_poly`、`coords` 系数 ∈ ℚ(A,B,…) |

```text
Context
  ├─ symbol_roles:  A → Parameter,  x → Variable
  ├─ assumptions:   Positive(A), Integer(B), …
  └─ active_var:    x   // diff / integrate / limit 默认变量
```

**A,B 本身不是假设**：未写 `assume` 时，`A` 仍可在式中出现，算法须 **保守**（见下）。

**各层配合：**

| 层 | `Parameter` 作用 | `assume` 作用 |
|----|------------------|---------------|
| **微积分** | `diff(...,x)` 时 A,B 视为常数 | 选分支，如 `integrate(abs(A*x),x)` 在 `assume(A>0)` 下化简 |
| **化简** | 系数环可含 A,B | 无 assume：`sqrt(A^2)` 不化为 `A`；`assume(B>0)` 可化简 |
| **solve** | 未知数 x，A,B 留系数 | `solve(A*x+B=0,x)` 需 `assume(A≠0)` 或 `Failed(NeedsAssume)` |
| **AlgExt / Poly\<AlgExtC\>** | A,B 进 `PolyCoeff` | 选 `rootof` 分支、实根、`factor` 可判定性；**不写入 coords** |

**参数字段 + assume 示例：**

```text
min_poly: t² - A              // 根 ±√A；系数环 ℚ(A)
assume(A>0)                   // 选正根 / 实根分支
assume(A, integer)            // 可走整数 factor 路径
```

**数据流（概念）：**

```text
Expr:  factor(A*t^2 - B)  或  rootof(..., poly1[1,0,-A])
           │
           ▼
Context::symbol_roles  →  A,B 为 Parameter
Context::assumptions   →  A>0 ?  A integer ?  A≠0 ?
           │
           ▼
expr_to_poly / PolyCoeff  →  Poly over ℚ(A,B)
           │
           ▼
factor / solve / AlgExtC
  ├─ 分支明确     →  Exact(...)
  ├─ 缺 assume    →  Failed(NeedsAssume)     // 非静默错答
  └─ 需分情形     →  piecewise（Phase D）
```

**Normative 设计原则：**

1. **角色与假设分离** — `A` 是参数 ≠ `A>0`；两者独立配置。
2. **默认保守** — 无 `assume` 时不猜符号正负、非零（[cas-long-term-vision.md](../cas-long-term-vision.md) §5.3）。
3. **assume 不进 AlgExt 坐标** — 假设在 `Context`，不在 `coords` / `min_poly` 里编码「A>0」。
4. **失败显式** — `CasResult::Failed(NeedsAssume(...))` 优于静默约分/错分支。
5. **`purge` 清假设** — `purge(A)` 清除 `vars` 绑定与相关 assumption（对标 upstream；见 **GIAC-204b**、[parser-token-map.md](../parser-token-map.md) §6）。

**giac-rs 现状与缺口：**

| 项 | 状态 |
|----|------|
| `Context::assumptions`（`Integer` / `Real` / `Positive`） | ✅ 字段存在 |
| `symbol_roles` | ❌ 未实现 |
| 关系假设 `A≠0`、`A>B` | ❌ 未实现 |
| 语句级 `assume` / `purge` 栈 | ❌ **GIAC-204b** |
| 算法侧 `check_assume`（e2r / factor / solve） | ❌ 未接 |

**与 AlgExt Horizon 的关系：**

| Horizon | A,B | assume |
|---------|-----|--------|
| **现在** | 不得出现在 `min_poly`/`coords` | 几乎未接算法 |
| **中期**（塔 + AlgExtC） | 仍 ℚ 系数；A,B 仅在普通 `Expr` 有理式 | 与 AlgExt 无直接耦合 |
| **长期**（阶段 3 + `Poly<AlgExtC>`） | 参数进系数环 | **必须**结合 assume 做分支、实根、约分 |

阶段 3 的 `PolyCoeff` 扩展「能写什么式子」；Phase C 的 assume 引擎扩展「在什么条件下可化简/求解」——二者在 `Poly<AlgExtC>::factor/roots` 处汇合。

---

## 已落地（Phase A + P1 过渡态）

| 项 | 状态 | 说明 |
|----|------|------|
| `AlgExtData { min_poly, coords, root_index }` | ✅ | `giac-core/src/algebra/alg_ext.rs` |
| `Expr::AlgExt(Arc<AlgExtData>)` | ✅ | `giac-core/src/expr.rs` |
| `from_rootof` / `to_rootof_expr` | ✅ | 有理 `poly1` 系数；显示仍走 `rootof(...)` |
| `ext_add` / `ext_sub` / `ext_mul` / `inv`（同域） | ✅ | P0 A-01/A-02 |
| `eval(rootof(...))` 提升 | ✅ | 可构造则 → `AlgExt` |
| `fold_algext_*` / `ratnormal` / `assert_equiv` | ✅ | P0 A-03/A-04 |
| `common_EXT` MVP | ✅ | primitive element 试探（B-05 子集） |
| `algext_square_roots` / `algext_sqrt_branches` | ✅ | 实根 + `Complex(0, AlgExt)` 复根 |
| `fold_complex_algext_*` | ✅ | **过渡态**：`Complex` 与 `AlgExt` 混合加减乘 |
| `biquadratic_rootof_roots` | ✅ | `solve(t^4-2=0)` 四根（2 实 + 2 复） |
| `quadratic_rootof` / `realroot` rootof 端点 | ✅ | B-03/B-04 子集 |
| `factor(x^2-2)` rootof 钩子 | ✅ | B-02 子集 |

**过渡态边界（§8.1，待阶段 2/2b/5 消除）：**

- `i` 为 `Context` 符号，非塔生成元
- 复代数数为 `Complex(0, AlgExt)` 与 `AlgExt` 双形态
- `solve` / `roots` 依赖 `biquadratic` / `algext_sqrt_branches` 特判，非 `Poly<AlgExtC>`
- `common_EXT` 未缓存、未完全对齐 upstream
- 系数环仍 **ℚ**；`evalf` 未接

---

## 算法适配清单

按 **阻塞 conformance / Phase B** 优先级排序。每项标注：必须识别的 `Expr` 形态、目标行为、上游对标、验收。

### P0 — 核心求值与化简（阻塞几乎所有代数路径）

| ID | 模块 | 算法 / 入口 | 须适配的行为 | 上游 | 验收 |
|----|------|-------------|--------------|------|------|
| A-01 | `giac-core::eval` | `eval_add` / `eval_mul` | `AlgExt ± AlgExt`、`AlgExt * AlgExt`、与 `Int`/`Rat` 混合 | `ext_add` / `ext_mul` | ✅ 同域 `fold_algext_*`；跨域 / `inv` 未做 |
| A-02 | `giac-core::eval` | `eval_frac` | 分母为 `AlgExt` 时 `inv_EXT` 或 `NotImplemented` | `inv_EXT` | 明确错误，不 panic |
| A-03 | `giac-simplify::normal` | `ratnormal` | `AlgExt` 视为原子常数；`Frac` 中含 `AlgExt` 的约分 | `ext_reduce` | `normal(rootof²-2)` 相关式 |
| A-04 | `giac-simplify::equiv` | `assert_equiv` | `AlgExt` 同域相等比较（`coords` + `min_poly`） | — | conformance `assert_equiv` |
| A-05 | `giac-core::display` | `format_expr` | 已走 `to_rootof_expr`；需稳定排序 | `symb_rootof` | golden 字符串不变 |

### P1 — 多项式与求根（阻塞 GIAC-205 深化、GIAC-232）

| ID | 模块 | 算法 | 须适配的行为 | 上游 | 验收 |
|----|------|------|--------------|------|------|
| B-01 | `giac-core::algebra::poly` | `expr_to_poly` | **拒绝**代数系数；**提升**走 `poly_alg_from_expr`（P1 stub）；契约 [expr-poly-conversion.md](../expr-poly-conversion.md) | `e2r` + `_EXT` | ✅ `TypeError` + `expr_contains_alg_coeff` |
| B-02 | `giac-poly` | `gcd` / `factor` / `roots` | 系数环 **`Poly<AlgExtC>`**；过渡态二次 rootof 钩子 | `gausspol` `algext_convert` | `factor(x^2-2)`；终态任意次数 |
| B-03 | `giac-solve` | `solve` / `roots` | 根为 **`AlgExtC`**；`Poly<AlgExtC>::roots` | `solve.cc` + `rootof` | `solve(t^4-2=0,t)` 四根 ✅ 过渡态 |
| B-04 | `giac-solve` | `sturm` / `realroot` | 实根：`Poly<AlgExt>`（`AlgExtC.im=0`） | `alg_ext.cc` `sturm` | `realroot(x^2-2)` |
| B-05 | `giac-core` | **`ExtensionTower::common`** | 塔式 `common_EXT` + 缓存 | `alg_ext.cc` L51–52 | ✅ `(√2)+(∛2)`；[lazy-common-tower](GIAC-lazy-common-tower-plan.md) T4a/T4b |
| B-06 | `giac-core` | **`AlgExtC`** + `canonicalize` | 复代数数一等类型；`i` 进塔 | — | `i·rootof(...)` 可 `+−×÷` |

### P2 — 微积分与极限（Phase B）

| ID | 模块 | 算法 | 须适配的行为 | 验收 |
|----|------|------|--------------|------|
| C-01 | `giac-calculus::diff` | `diff` | `AlgExt` 对 `x` 为常数（0） | `diff(rootof(...),x)=0` |
| C-02 | `giac-calculus::integrate` | 有理式 / `partfrac` | 分母含不可约二次 → `AlgExt` 对数项 | `∫1/(x^2-2)dx` |
| C-03 | `giac-calculus::limit` | 有限点代入 | `AlgExt` 作为常数代入 | `limit(rootof(...),x,0)` |
| C-04 | `giac-calculus::risch` | 塔扩展 | 代数塔层识别 `AlgExt` | GIAC-229+ |

### P3 — 线代、几何、数值（后期）

| ID | 模块 | 算法 | 说明 |
|----|------|------|------|
| D-01 | `giac-linalg` | `egv` / `jordan` / `charpoly` | 特征值可为 `AlgExt` |
| D-02 | `giac-geo` | 距离、交点 | `testgeo` 依赖 |
| D-03 | `giac-core` | `evalf` / `alg_evalf` | `approx_rootof`、`proot` 浮点逼近 |
| D-04 | `giac-solve` | `froot` / `froots` | GIAC-232；因子根列表 |

### P4 — 解析与程序（低优先）

| ID | 模块 | 说明 |
|----|------|------|
| E-01 | `giac-parse` | 解析结果直接构造 `AlgExt`（若 `convert_rootof`） |
| E-02 | `giac-prog` | `keep_algext` 伪变量行为 |
| E-03 | `giac-simplify::factor` | `factor` 输出含 `AlgExt` 因子 |

---

## 推荐实施顺序

**已完成（过渡态）：**

```text
P0 (A-01…A-05) → P1 子集 (B-01…B-05 MVP, Complex+AlgExt, t^4-2)
```

**目标架构（§8.6 中期/长期，§8.7 交付顺序）：**

```text
ExtensionTower + ExtensionField::common（B-05 重构）
  → AlgExtC + canonicalize（B-06）
  → Poly<AlgExtC> roots 二次/双二次（B-03 终态）
  → 移除 biquadratic / algext_sqrt_branches 特判
  → Poly<AlgExtC> factor/gcd（B-02 终态）
  → sturm / realroot on Poly<AlgExt>（B-04）
  → C-01…C-02 → AlgExtC::evalf（D-03）
```

---

## 数据模型

### 当前（过渡态，Phase A + P1）

```rust
pub struct AlgExtData {
    pub min_poly: Vec<ExprArc>,
    pub coords: Vec<ExprArc>,
    pub root_index: Option<u32>,
}
// 复根：Expr::Complex(Expr::int(0), AlgExt(...))，i 为符号
```

### 目标（§8 normative）

```rust
pub enum ExtensionTower { Base, Adj { parent, min_poly, embed_parent, degree } }
pub struct ExtensionField { pub tower: Arc<ExtensionTower>, /* mul_table */ }

pub struct AlgExtData {
    pub field: Arc<ExtensionField>,
    pub coords: Vec<ExprArc>,
    pub root_index: Option<u32>,
}

pub struct AlgExtCData {
    pub field: Arc<ExtensionField>,
    pub re: Vec<ExprArc>,
    pub im: Vec<ExprArc>,
    pub root_index: Option<u32>,
}
// Expr::AlgExtC(Arc<AlgExtCData>)
// giac-poly: Poly<AlgExtC>
```

**与 `Func(RootOf)` 的关系：**

- **输入：** `rootof(num, minpoly)` 求值 → 构造/扩展 `ExtensionTower` → `AlgExtC` 或 `AlgExt`（`im=0`）
- **输出：** `format_expr(AlgExtC)` → `rootof(...)` 或 `(re)+(im)*i`（兼容 golden）
- **退化：** 符号参数、元数不足 → 保留 `Func(RootOf, …)`

---

## 测试计划

| 层级 | 命令 / 用例 |
|------|-------------|
| 单元 | `cargo test -p giac-core alg_ext`；`AlgExtC` 加减乘除逆 |
| 塔 | `ExtensionTower::adj`；`common(ℚ(√2), ℚ(∛2))` 次数 6 |
| solve | `cargo test -p giac-solve rootof`；终态 `Poly<AlgExtC>::roots` |
| conformance | `batch2_giac205_solve_quadratic_rootof`；`solve(t^4-2=0,t)` 四根 |
| 属性 | `α²=2`；`(i·α)²=−2`；`AlgExtC` 范数除法 |
| 参数 + assume | `normal((A*x)^2/A)` 分 assume 有无；`solve(A*x+B=0,x)` + `assume(A≠0)`（阶段 3 / Phase C，§8.9） |

---

## 开放问题（HITL）

1. **系数环：** Phase C / 阶段 3 前 `min_poly` 不含 A,B（§8.9）；参数化与 assume 同期交付。
2. **`keep_algext`：** 默认保留 `rootof` 还是 `evalf` 展开？
3. **`CasResult::Algebraic`：** WASM 边界暴露 `AlgExtC` 还是 `Exact(Expr)`？
4. **塔 vs 扁平 `common`：** 每次 `rootof` 是否必须 `Adj` 新层，还是允许 flatten 到单个 `min_poly`（upstream 两种路径并存）？
5. **`i` 层时机：** `AlgExtC` 引入时是否默认 adjoin `t²+1`，还是按需 `adj`（影响 `solve` 塔深度）？
6. **assume 与 AlgExt 分支：** `assume(A>0)` 是否足够选 `rootof` 正支，还是需要显式 `root_index` / `proot`（§8.9）？

---

## 参考

- C++：`ext_add`、`ext_mul`、`inv_EXT`、`common_EXT`、`ext_reduce`、`horner_rootof`；`sym2poly.cc::check_assume`
- 文档：`rust-migration-plan.md` L133–L134、L180–L185；[cas-long-term-vision.md](../cas-long-term-vision.md) §5.2–5.3（参数 A,B 与 assume）；[parser-token-map.md](../parser-token-map.md) §6
- 已关闭浅层：`GIAC-205`（二次 `rootof` 符号）；本 issue 跟踪 **真扩域类型** 与各算法适配
