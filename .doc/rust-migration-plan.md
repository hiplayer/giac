# Giac CAS → Rust 实现方案

基于当前 **53 个 CTest / ~220 种 API / 19 个功能域**（见 `functional-coverage.md`）界定迁移范围。原则：**测试即规格**，先达到现有 golden 回归等价，再扩展未覆盖功能。

> **长远愿景：** 数学边界、初等函数能力分层、参数 A/B/C/D 语义见 [cas-long-term-vision.md](cas-long-term-vision.md)。本文 Phase 0–5 对应彼处 **Phase A（基线）**。

---

## 1. 范围界定

### 1.1 MVP（必须实现）

| 阶段 | 功能域 | 代表测试 | API 约数 |
|------|--------|----------|----------|
| P0 | 基础设施 | — | Expr/解析/上下文 |
| P1 | 基础算术与复数 | `test_cas_basic`, `testcas` | ~15 |
| P2 | 化简与三角 | `test_trig`, `testnormalize` | ~23 |
| P3 | 多项式 | `test_poly*`, `testfactor`, `test_modular` | ~40 |
| P4 | 线性代数 | `test_linalg*`, `test_linalg_decomp` | ~35 |
| P5 | 方程与 Sturm | `test_solve*`, `test_sturm*` | ~8 |
| P6 | 微积分 | `test_integrate*`, `testintegrate`, `testother` | ~25 |
| P7 | 级数/极限/分式 | `test_series`, `testlimit`, `testpartfrac*` | ~8 |
| P8 | ODE | `test_desolve*` | 1 |
| P9 | 数论/置换 | `test_ifactor*`, `test_permu*` | ~16 |
| P10 | 数值/变换/向量 | `test_numerical`, `test_laplace`, `test_vector_calc` | ~20 |
| P11 | 特殊函数/概率 | `test_special`, `test_probstat` | ~15 |
| P12 | 程序与替换 | `test_subst`, `test_prog` | ~6 |
| P13 | 导出 | `test_export` | 2 |
| P14 | 综合回归 | `flanex`, `testcas` | 交叉验证 |

### 1.2 明确不做（当前无测试或已排除）

- GUI：`icas.cc`, `Graph*.cc`, FLTK 相关
- 兼容层：`maple.cc`, `ti89.cc`, `rpn.cc`, `pari.cc`
- 外部依赖：`cocoa.cc`（CoCoA Groebner 后端）
- 稀疏/3D：`sparse.cc`, `plot3d.cc`
- `gbasis` 完整求基（仅 `greduce` 有测试且可用）

---

## 2. Rust 工程结构

```
giac-rs/
├── Cargo.toml                 # workspace
├── crates/
│   ├── giac-core/             # Expr, Ident, Domain, EvalError, Context
│   ├── giac-parse/            # 词法/语法（兼容 giac 脚本子集）
│   ├── giac-simplify/         # normal, ratnormal, texpand, trig 化简
│   ├── giac-poly/             # 多项式环、gcd/factor/resultant/mod
│   ├── giac-linalg/           # 矩阵、rref/svd/lu/qr、特征值
│   ├── giac-solve/            # solve, linsolve, fsolve, sturm, realroot
│   ├── giac-calculus/         # diff, integrate, limit, series, partfrac
│   ├── giac-ode/              # desolve
│   ├── giac-transform/        # laplace, ilaplace, fourier_*
│   ├── giac-num/              # ifactors, primes, perm, mod 运算
│   ├── giac-special/          # gamma, Beta, orthogonal poly, 分布
│   ├── giac-vector/           # curl, div, grad, potential
│   ├── giac-prog/             # subst, sum, product, for/while
│   ├── giac-emit/             # latex, mathml
│   ├── giac-groebner/         # greduce（Phase 3+，gbasis 后置）
│   └── giac-ffi/                # 可选：C ABI 供 C++ giac 对照
├── apps/
│   └── giac-cli/              # 兼容 `giac script` 行为
└── tests/
    └── conformance/           # 消费 bin/ + check/ golden 的集成测试
```

**WASM 输出目标**

项目最终需编译为 WebAssembly，在浏览器端运行完整的 CAS 功能。所有 crate 选型必须满足
`wasm32-unknown-unknown` 目标开箱编译，禁止引入含 C FFI 的传递依赖。具体约束：

| 约束 | 说明 |
|------|------|
| 数值 linalg | 使用 `nalgebra`（纯 Rust，no_std 可选），**禁止** `faer`（SIMD intrinsics 不兼容 WASM） |
| 大整数 | `num-bigint` 可编译到 WASM；`rug`（GMP FFI）**不可用**，WASM 场景降级为 `num-bigint` |
| 任意精度浮点 | WASM 下暂用 `f64`，后续可选 `soft-float` 纯 Rust MPFR 替代 |
| CLI / FFI crate | `giac-cli`、`giac-ffi` 不参与 WASM 构建，通过 `cfg` feature gate 隔离 |
| 输出格式 | WASM 层暴露 `eval_to_string(input: &str) -> String`，格式化走 `giac-core::format_expr` |

```
┌──────────────────────────────────────────────┐
│  Browser (WASM)                              │
│  giac-wasm ──┐                               │
│              ├── giac-core ── giac-linalg     │
│              ├── giac-parse                   │
│              └── giac-simplify ── giac-poly   │
│  (giac-cli / giac-ffi 不参与 WASM 编译)       │
└──────────────────────────────────────────────┘
```

| 能力 | Crate | 说明 |
|------|-------|------|
| 大整数 | `num-bigint` / `rug` | 对标 libtommath；模运算密集处可用 `rug` |
| 有理数 | `num-rational` | 精确系数；与 giac `fraction` 对齐 |
| 浮点 | `f64` + 规范化；可选 `rug::Float`（仅 native） | `evalf`/`fsolve` 需可配置精度；WASM 下无 `rug` |
| 线性代数 | `nalgebra` | `f64` 矩阵；WASM 输出必需；符号矩阵自研 |
| 随机 | `rand` + `rand_distr` | `randNorm` 等 |
| 序列化/测试 | `insta`, `serde_json` | golden 快照与 diff |
| 解析 | `logos` + `lalrpop`/`chumsky` | 替代 flex/bison |

---

## 3. 核心类型设计（替代 `gen`）

Giac 的 `gen` 是 tagged union + 引用计数（`dispatch.h` 中 22 种 `type` + `subtype`），所有 API 统一传递。Rust 侧**不做一个 Rust 版 `gen`**，而是分层强类型 + 显式转换。

### 3.1 `gen` → Rust 类型对照表

| C++ `gen` (`type`) | giac 含义 | 建议 Rust 类型 | 所在 crate | MVP |
|--------------------|-----------|----------------|------------|-----|
| `_INT_` | 小整数 | `Expr::Int(i32)` 或 `BigInt` | `giac-core` | ✓ |
| `_ZINT` | GMP 大整数 | `Expr::Int(BigInt)` | `giac-core` | ✓ |
| `_DOUBLE_` | 机器浮点 | `Numeric::F64(f64)` | `giac-core` | ✓ |
| `_FLOAT_` | 即时浮点 | `Numeric::F64(f64)` | `giac-core` | ✓ |
| `_REAL` | MPFR 任意精度 | `Numeric::Float(PrecFloat)` | `giac-core` | ✓ |
| `_FRAC` | `Tfraction<gen>`，分子分母可为符号 | `Expr::Frac(num, den)` | `giac-core` | ✓ |
| `_CPLX` | `gen[2]` 复数 | `Expr::Complex(re, im)` | `giac-core` | ✓ |
| `_IDNT` | 标识符 | `Expr::Symbol(Ident)` | `giac-core` | ✓ |
| `_SYMB` | `symbolic{sommet, feuille}` | `Add`/`Mul`/`Pow`/`Func` | `giac-core` | ✓ |
| `_VECT` + `_SEQ__VECT` | 序列 `[a,b,c]` | `Expr::Seq(Vec<Expr>)` | `giac-core` | ✓ |
| `_VECT` + `_LIST__VECT` | 列表 | `Expr::List(Vec<Expr>)` | `giac-core` | ✓ |
| `_VECT` + `_SET__VECT` | 集合 | `Expr::Set(Vec<Expr>)` | `giac-core` | 后期 |
| `_VECT` + `_MATRIX__VECT` | 矩阵 | `Expr::Matrix(Vec<Vec<Expr>>)` | `giac-core` | ✓ |
| `_VECT` + `_POLY1__VECT` | 一元多项式系数向量 | `Poly1<D>`（内部） | `giac-poly` | ✓ |
| `_VECT` + `_INTERVAL__VECT` | 区间 | `Expr::Interval(lo, hi)` | `giac-core` | 后期 |
| `_VECT` + `_ASSUME__VECT` | 假设 | `Context::assumptions` | `giac-core` | ✓ |
| `_VECT` + `_PNT__VECT` 等 | 几何对象 | `GeomObj` | `giac-geo` | 后期 |
| `_VECT` + `_PRG__VECT` | 程序语句 | `Stmt` / `Program` | `giac-prog` | ✓ |
| `_POLY` | `tensor<gen>` 稀疏多项式 | `Poly<Coeff>` | `giac-poly` | ✓ |
| `_SPOL1` | 稀疏一元多项式 | `SparsePoly1<Coeff>` | `giac-poly` | ✓ |
| `_EXT` / `_ROOT` | 代数扩域 / `rootof` | `Expr::AlgExt(AlgExtData)` | `giac-core` | ✓ |
| `_MOD` | `n mod m` | `Expr::Mod(val, m)` + `ModInt` | `giac-core`/`giac-poly` | ✓ |
| `_STRNG` | 字符串值 | `Expr::Str(String)` 或 `Emit` 层 | `giac-core` | 后期 |
| `_FUNC` | 一等函数指针 | `Callable::Builtin(FuncKind)` / `UserFn` | `giac-prog` | 后期 |
| `_MAP` | `map<gen,gen>` | `Expr::Map(BTreeMap<Expr,Expr>)` | `giac-core` | 后期 |
| `_USER` | 用户自定义类型 | `UserValue` trait object | — | 不做 |
| `_EQW` | 排版公式树 | `EqwNode` | `giac-emit` | 后期 |
| `_GROB` / `_POINTER_` | 图形 / 裸指针 | — | — | 不做 |

**MVP 缺口说明**（现有 53 项测试会触发，§3 初版 `Expr` 未列出）：

- **`_EXT` / `_ROOT`**：`testcas`、`testgeo`、`test_integrate_more` 大量 `rootof(...)`；需扩域算术，不能仅用 `Func(RootOf, …)` 打印。
- **`_MOD`**：`test_modular`、`flanex` 的 `x % 13`、`chinrem`；系数环为 `ℤ/mℤ`，非普通 `Mod` 函数调用。
- **`_REAL`**：`test_numerical`、`testcas` 的 `evalf`/`float2rational`；`f64` 不够。
- **`_FRAC`**：有理式 `(x+1)/(x^2-1)` 在 `normal` 前常保持分式形态；`Ratio<BigInt>` 仅覆盖数字。
- **`Stmt`**：`test_prog` 的 `for`/`while`；表达式 AST 不能承载语句。
- **`Seq` vs `List`**：giac 序列与列表语义不同；解析输出必须区分。

### 3.2 建议 Rust 类型（完整草图）

MVP 约定：**节点统一用 `Arc<Expr>`**（见 §3.5）。下列草图子节点写 `Arc<Expr>`；复合变体用 `Arc` 包裹整棵子树。

```rust
// ── giac-core ──────────────────────────────────────────

/// 顶层符号表达式（用户可见 CAS 值）
pub enum Expr {
    Int(BigInt),
    Rat(Ratio<BigInt>),                      // 仅数字有理数；符号有理式用 Frac
    Frac(Arc<Expr>, Arc<Expr>),              // 对标 _FRAC
    Mod(Arc<Expr>, Arc<Expr>),               // 对标 _MOD（表面语法 a % m）
    Complex(Arc<Expr>, Arc<Expr>),
    Symbol(Ident),
    Add(Vec<Arc<Expr>>),
    Mul(Vec<Arc<Expr>>),
    Pow(Arc<Expr>, Arc<Expr>),
    Func(FuncKind, Vec<Arc<Expr>>),          // 对标 _SYMB
    Seq(Vec<Arc<Expr>>),                     // [a,b,c]  对标 _SEQ__VECT
    List(Vec<Arc<Expr>>),                    // 列表     对标 _LIST__VECT
    Matrix(Vec<Vec<Arc<Expr>>>),             // 对标 _MATRIX__VECT
    Relation(RelOp, Arc<Expr>, Arc<Expr>),
    AlgExt(Arc<AlgExtData>),                 // 对标 _EXT / rootof
    Interval(Arc<Expr>, Arc<Expr>),          // 对标 _INTERVAL__VECT
    Str(String),                             // 对标 _STRNG
    Undefined,
}

/// 代数扩域元素（对标 ref_algext / rootof）
pub struct AlgExtData {
    pub min_poly: Vec<Arc<Expr>>,     // 最小多项式系数 [c0,c1,...]
    pub coords: Vec<Arc<Expr>>,       // 在扩域基下的坐标
    pub root_index: Option<u32>,      // rootof 选取第 k 个根
}

/// 任意精度数值（evalf / fsolve 数值路径）
pub enum Numeric {
    Rat(Ratio<BigInt>),
    F64(f64),
    Float(PrecFloat),                  // rug / 任意精度，对标 _REAL
}

pub enum FuncKind { Sin, Cos, Integrate, Factor, Gcd, /* ... 220+ */ }

pub struct Context {
    pub vars: HashMap<Ident, Arc<Expr>>,
    pub assumptions: Vec<Assumption>,  // assume[integer] 等
    pub complex_mode: bool,
    pub epsilon: f64,
    pub series_order: u32,
    pub float_digits: u32,            // evalf 精度
}

pub enum EvalError {
    TooFewArgs, NotDifferentiable, GroebnerFailed,
    DivisionByZero, TypeError, /* ... */
}

// ── giac-poly ──────────────────────────────────────────

pub struct Poly<C> { /* 稀疏/稠密，内部表示 */ }
pub struct SparsePoly1<C> { /* 对标 _SPOL1 */ }
pub struct Poly1<C>(Vec<C>);         // 对标 _POLY1__VECT

/// 模整数（系数环）
pub struct ModInt { pub val: BigInt, pub modulus: BigInt }

/// 模多项式环 ℤ/mℤ[x]
pub struct PolyMod {
    pub poly: Poly<ModInt>,
    pub modulus: BigInt,
}

// Expr ↔ Poly 桥接（对标 e2r / r2e）
impl Expr {
    pub fn to_poly(&self, vars: &[Ident]) -> Result<Poly<Ratio<BigInt>>, ConvError>;
}
pub fn from_poly(p: &Poly<Ratio<BigInt>>, vars: &[Ident]) -> Arc<Expr>;

// 算法签名：不可变语义，返回新 Arc（子树未变则 Arc::clone 复用）
pub fn normal(expr: &Expr, ctx: &Context) -> Result<Arc<Expr>, EvalError>;

// ── giac-prog ──────────────────────────────────────────

pub enum Stmt {
    Assign(Ident, Arc<Expr>),
    For { var: Ident, from: Arc<Expr>, to: Arc<Expr>, body: Vec<Stmt> },
    While { cond: Arc<Expr>, body: Vec<Stmt> },
    Print(Vec<Arc<Expr>>),
    ExprStmt(Arc<Expr>),
}
pub struct Program(Vec<Stmt>);

pub enum Callable {
    Builtin(FuncKind),
    UserFn { params: Vec<Ident>, body: Program },
}

// ── giac-linalg（数值快路径，不进 Expr）────────────────

// svd / fsolve / newton 使用 f64 矩阵，与符号 Matrix<Expr> 分离
```

### 3.3 转换边界（不进 `Expr` 的类型）

| 转换 | 方向 | 说明 |
|------|------|------|
| `e2r` / `r2e` | `Expr ↔ Poly<Rat>` | gcd/factor/resultant 入口；失败返回 `ConvError` |
| `mod` 语法 | `Expr → PolyMod` | `factor(p) mod 13` 在模多项式环运算 |
| `evalf` | `Expr → Numeric` | 符号层下沉到任意精度浮点 |
| `rootof` 运算 | `AlgExtData` 内部 | `ext_add/mul/div` 对标 `alg_ext.cc` |
| 语句执行 | `Program → Expr` | `giac-cli` 逐句 `eval`，变量写入 `Context` |
| 排版 | `Expr → EqwNode → String` | `latex`/`mathml` 独立 AST，不复用计算 `Expr` |

### 3.4 设计要点

1. **符号层与数值层分离**：`Expr` 用于 CAS；`Numeric` / `Matrix<f64>` 用于 `evalf`/`fsolve`/`svd`。
2. **多项式快路径独立**：`Poly`/`PolyMod`/`SparsePoly1` 不进用户-facing `Expr` 枚举；通过 `to_poly()` 进入（对标 `_POLY`/`_SPOL1`）。
3. **代数数是一等公民**：`AlgExt` 变体 + `giac-core::alg_ext` 模块，不可推迟到「后期」（`testcas`/`geo` 依赖）。
4. **模算术是一等公民**：`Mod`/`ModInt`/`PolyMod`，`flanex` 与 `test_modular` 依赖。
5. **语句与表达式分离**：`Stmt`/`Program` 在 `giac-prog`；`Expr` 不承载 `for`/`while`。
6. **错误类型**：`Result<Expr, EvalError>`，避免 giac 把错误塞进 `gen` 字符串。
7. **不可变 + 结构共享**：变换返回新 `Arc<Expr>`，替代 `gen` 浅拷贝 + 原地修改（详见 §3.5）。

### 3.5 MVP 内存策略：仅 `Arc<Expr>`

**已定**：MVP 阶段不引入 arena（`bumpalo` 等），节点统一为 `Arc<Expr>`。

| 层面 | MVP 做法 |
|------|----------|
| **算法语义** | 不可变：`fn algo(expr: &Expr, ctx: &Context) -> Result<Arc<Expr>, EvalError>` |
| **结构共享** | 子树未变 → `Arc::clone`；子树变了 → 只重建路径上的父节点 |
| **会话状态** | `Context.vars: HashMap<Ident, Arc<Expr>>` |
| **多项式/数值内核** | `Expr → Poly` / `Matrix<f64>` 后可在专用类型上可变运算，不绑 AST 策略 |
| **arena** | **不做**；待 profiling 证明分配是瓶颈后再评估（如 `expand`/`series` 热点） |

算法编写约定：

```rust
// ✓ MVP：借用输入，返回新 Arc
fn diff(expr: &Expr, var: &Ident, ctx: &Context) -> Result<Arc<Expr>, EvalError>;

// ✗ 避免：原地修改 AST
fn diff_in_place(expr: &mut Expr, ...) -> ...;

// ✓ 子树复用示例
if Arc::ptr_eq(&old_child, &new_child) {
    return Ok(Arc::clone(expr_arc));  // 整棵子树不变
}
```

---

## 4. 解析与脚本兼容

### 4.1 目标

- 能直接执行现有 `bin/test_*` 与 `check/test*` 脚本（giac 语法子集）。
- 一行一条语句，`;` 分隔；支持 `:=`、`^`、`**`、`mod`、矩阵 `[[...]]`。

### 4.2 实现策略

```
源代码 → Tokenizer(logos) → Parser(lalrpop) → Vec<Stmt> → Evaluator
```

- **Stmt** 枚举：`ExprStmt(Expr)`, `Assign(Ident, Expr)`, `FuncDef(...)`, `For/While`。
- **内置函数表**：`HashMap<&'static str, BuiltinFn>`，与 C++ `at_*` 注册表对应（MVP ~210 项，见 [`.doc/builtin-api-map.md`](builtin-api-map.md)）。
- **Token / 优先级 / 多模式 / Context**：见 [`.doc/parser-token-map.md`](parser-token-map.md)。

### 4.3 与 C++ 对照

初期可保留 `giac` 二进制，conformance 测试双跑：

```
giac script.in  > cpp.out
giac-cli script.in > rust.out
diff cpp.out rust.out
```

---

## 5. 分阶段实施路线

> **各阶段通用门禁**（与 golden 验收并列）：`cargo test --workspace` + **`cargo ci-clippy`** 全绿。见 §8.4.1、[supplement §7](rust-migration-supplement.md#7-工程门禁)。

### Phase 0 — 骨架（2–3 周）

- [ ] workspace + `giac-core` + `giac-parse`
- [ ] `Expr` 运算：`+ - * / ^`，`simplify` 基础
- [ ] `giac-cli` 读 stdin/文件，打印结果
- [ ] conformance harness：遍历 `bin/test_cas_basic`
- [ ] `clippy.toml` + `[workspace.lints]` + `cargo ci-clippy` 可跑通

### Phase 1 — 基础域（3–4 周）

对标：`test_cas_basic`, `testcas` 前半（复数、gcd、normal）

| 模块 | 实现 |
|------|------|
| `giac-core` | `abs`, `arg`, `conj`, `re`, `im`, `sign`, `gcd`（整数） |
| `giac-simplify` | `normal`, `ratnormal`, `expand` |
| `giac-num` | 整数 `gcd`, `lcm` |

**验收**：`giac_check_cas` 前 50 行 golden 子集通过。

### Phase 2 — 多项式环（4–6 周）

对标：`test_poly*`, `test_factor`, `giac_check_factor`, `test_modular`

| 模块 | 实现 |
|------|------|
| `giac-poly` | `Poly` 表示、单项序、`gcd`/`lcm`/`quo`/`rem` |
| | `factor`（有限域 + 整数）、`egcd`, `resultant`, `roots` |
| | `chinrem`, 模多项式 `gcd`/`rref` |
| `giac-groebner` | **仅** `greduce`（多项式长除法链） |

**验收**：`test_poly_ext`, `test_modular`, `giac_check_factor` 全通过。

### Phase 3 — 线性代数（3–4 周）

对标：`test_linalg*`, `test_linalg_decomp`, `test_gauss_ext`

| 模块 | 实现 |
|------|------|
| `giac-linalg` | 符号矩阵、`rref`, `det`, `linsolve`, `ker`, `image` |
| | 数值：`lu`, `qr`, `svd`（`nalgebra`）——**必须使用 `nalgebra`，项目要求 WASM 输出** |
| | `jordan`, `egv`, `charpoly`, `gramschmidt` |

> **WASM 输出 & 为什么选 `nalgebra` 而非 `faer`**
>
> 项目最终目标之一是将 CAS 编译为 WebAssembly，在浏览器端运行。`nalgebra` 是纯 Rust
> 实现、无 C FFI 依赖，`wasm-pack` / `wasm32-unknown-unknown` 目标开箱即用；而 `faer`
> 内部使用 `std::simd` 与平台相关 intrinsics，当前在 WASM 目标下编译存在兼容性问题且
> 性能退化严重。因此 Phase 3 数值线性代数一律使用 `nalgebra`：
>
> - `lu` → `nalgebra::linalg::LU`
> - `qr` → `nalgebra::linalg::QR`
> - `svd` → `nalgebra::linalg::SVD`
> - 向量/矩阵基本运算 → `nalgebra::DMatrix<f64>` / `DVector<f64>`
>
> 如未来 `faer` 的 WASM 支持成熟，可作为可选后端切换，但 MVP 以 `nalgebra` 为准。

### Phase 4 — 求解与微积分（6–8 周，难度最高）

对标：`test_solve*`, `test_integrate*`, `testintegrate`, `testother`, `testlimit`

| 模块 | 实现 |
|------|------|
| `giac-solve` | 多项式 `solve`, `linsolve`, `fsolve`/`newton`, `sturm`/`realroot` |
| `giac-calculus` | `diff`/`derive`, `limit`, `series`/`taylor`, `partfrac` |
| | `integrate`：先规则表 + 部分分式，再逐步 port `risch` |
| `giac-ode` | `desolve` 线性常系数 ODE |

**验收**：`giac_check_integrate`, `giac_check_limit`, `giac_check_other`。

### Phase 5 — 扩展域（4–5 周）

对标：其余 bin 脚本 + `flanex` 子集

- `giac-transform`：`laplace`/`ilaplace`/`fourier_*`
- `giac-vector`：`curl`, `divergence`, `laplacian`, `potential`
- `giac-num`：`ifactors`, `isprime`, `permu*`
- `giac-special`：`Beta`, `gamma`, `zeta`, `rand*`
- `giac-prog`：`subst`, `sum`, `product`, `for`/`while`
- `giac-emit`：`latex`, `mathml`（模板输出）

**验收**：全量 53 CTest conformance（Rust 侧等价测试）。

---

## 6. 测试与等价策略

### 6.1 三层测试

```
单元测试 (crate内)     →  每个 API 的最小用例
黄金回归 (conformance) →  check/*.out diff
覆盖追踪 (optional)    →  cargo-tarpaulin 对标 gcov
```

### 6.2 Golden 规范化

移植 `cmake/giac_check_normalize.sh` 逻辑到 Rust：

- `cas_floats`：浮点格式统一到 10 位有效数字
- `geo`：绘图 ID 替换（若实现 plot）

### 6.3 已知差异处理

| 差异 | 策略 |
|------|------|
| 随机函数 `rand*` | 测试只验证类型/范围，或 seed 固定 |
| `gbasis` 失败 | MVP 不实现；`greduce` 用自研约化 |
| 表达式排序不同但等价 | 可选 `canonicalize` 后再 diff |
| `factor` segfault 用例 | Rust 应返回 `Err` 而非 panic |

### 6.4 算法偏离原则（不必照抄 giac）

**已定**：Rust 实现**不绑定** giac 源码的算法选择与实现细节；绑定的是**数学语义 + 测试规格**。giac C++ 仅作参考实现与回归对照，不是权威。

#### 何时可以改算法

| 情形 | 做法 | 示例 |
|------|------|------|
| giac 有 bug（崩溃、错答、静默 UB） | 用正确算法；**不**复现 bug | `factor` segfault → 返回 `Err` |
| 有更优算法（复杂度、稳定性） | 可替换；验证见下文 | F4 替代 F5、更稳的 `rref` 主元策略 |
| 输出形式不同但数学等价 | `canonicalize` 后 diff，或等价判定 | `x+1` vs `1+x`，`sqrt(2)/2` vs `1/sqrt(2)` |
| 输出形式不同且无法规范化 | 更新 golden + 记录偏离 | 主动选择更简单闭式解 |
| 无测试覆盖的行为 | 自订语义，补测试 | 新 `integrate` 规则 |

#### 验证层次（由强到弱，组合使用）

```
① 数学等价判定     diff(new, giac) 失败时，用 CAS 自检：subs、expand、normal 差为零
② 属性测试         proptest：交换律、gcd*lcm=n*m、integrate 再 diff 还原等
③ Golden 回归      check/*.out、bin 脚本 — MVP 门禁，但不是唯一真理
④ 单元测试         每个 API 最小用例 + 边界（除零、空矩阵、不可积）
⑤ 对照 giac 双跑   conformance harness 可选模式：同输入比 giac 与 giac-rs
⑥ 第三方参照       Sage/SymPy/文献算法 — 改算法时优先用独立来源佐证
```

**MVP 门禁**：

1. **测试**：`cargo test --workspace` 全绿
2. **Clippy**：`cargo ci-clippy` 全绿（见 §8.4.1）
3. **Golden**：53 项 golden 全绿 **或** 已记录的已知偏离项全部有等价证明 / 更新后的 golden

#### 偏离处理流程

1. **发现**：Rust 与 giac golden 不一致，或 giac 本身崩溃/明显错答。
2. **归类**：bug 修复 / 等价不同形 / 真语义分歧。
3. **验证**：
   - bug 修复、真分歧 → 独立证明正确（手算、`expand(sub(a,b))=0`、文献）。
   - 等价不同形 → 实现 `canonicalize` 或 `assert_equiv(a, b, ctx)`。
4. **落盘**：写入 [`.doc/known-divergences.md`](known-divergences.md)（输入、giac 输出、Rust 输出、理由、验证方式）。
5. **更新测试**：修 golden **或** conformance 对偏离项走等价判定而非字面 diff。

#### `assert_equiv` 示意（conformance 增强）

完整规格见 [`.doc/conformance-testing.md` §3](conformance-testing.md)。

```rust
/// 判定两个 Expr 是否数学等价（比字面 golden 更权威）
fn assert_equiv(a: &Expr, b: &Expr, ctx: &Context) -> bool {
    let d = normal(&sub(a, b)?, ctx)?;
    is_zero(&d)
}
```

积分、求解等多解 API：允许 Rust 返回不同但等价的闭式，只要 `assert_equiv` 通过；必要时 golden 改用「等价类」或更新为更简答案。

#### 不必照抄的具体项

| 领域 | giac 做法 | Rust 可自选 |
|------|-----------|-------------|
| 化简 | `normal` 启发式链 | 自己的 canonical form + 规则优先级 |
| 积分 | 规则表 + 部分 Risch | 先完备规则表，Risch 子集按需 |
| Groebner | CoCoA F5 + 原生 F4 混用 | 自研 `greduce`（`gbasis` 已后置） |
| 因式分解 | `ezgcd` + 启发式 | Zassenhaus / Wang 等标准流程 |
| 浮点 | 混用 `double` / MPFR | 统一 `PrecFloat`，边界更可预测 |

**结论**：算法可以自由选型；**验证靠数学等价 + 测试套件 +  documented 偏离**，不是靠与 `solve.cc` 行级一致。

---

## 7. 模块映射（C++ → Rust）

| C++ 源文件 | Rust crate | 已测 API |
|------------|------------|----------|
| `gen.cc`, `symbolic.cc` | `giac-core` | 基础运算符 |
| `input_lexer.cc`, `input_parser.yy` | `giac-parse` | 全部脚本 |
| `usual.cc` | `giac-simplify` | trig, normal |
| `sym2poly.cc`, `ezgcd.cc`, `modpoly.cc` | `giac-poly` | gcd, factor, mod |
| `lin.cc`, `vecteur.cc` | `giac-linalg` | 矩阵全套 |
| `solve.cc`, `csturm.cc` | `giac-solve` | solve, sturm |
| `intg.cc`, `derive.cc`, `risch.cc`, `series.cc` | `giac-calculus` | 微积分 |
| `desolve.cc` | `giac-ode` | desolve |
| `intgab.cc` | `giac-transform` | laplace, fourier |
| `ifactor.cc`, `permu.cc` | `giac-num` | 数论、置换 |
| `moyal.cc`, `misc.cc` | `giac-special` | 特殊函数 |
| `subst.cc`, `prog.cc` | `giac-prog` | subst, 控制流 |
| `tex.cc`, `mathml.cc` | `giac-emit` | latex, mathml |
| `solve.cc`/`cocoa.cc` | `giac-groebner` | greduce only |
| `libtommath` | `num-bigint`/`rug` | 大整数 |

---

## 8. 内存安全与 unsafe 纪律

Rust 重写的主要收益之一是**默认内存安全**。本计划要求：业务逻辑用 safe Rust 实现，unsafe 仅作为可选边界层，且受 CI 与评审约束。

### 8.1 总体原则

| 原则 | 说明 |
|------|------|
| **Safe by default** | 符号运算、化简、求解、微积分等核心算法全部用 safe Rust |
| **unsafe 隔离** | 仅 `giac-ffi`（及未来明确的 shim crate）允许 `unsafe` |
| **不移植 C++ 指针语义** | 禁止把 `gen` 的引用计数/别名模型直译为裸指针或 `unsafe` 共享可变状态 |
| **纯 Rust 优先** | 大整数/有理数用 `num-bigint`、`num-rational`（或 `rug`），避免直接绑 GMP |
| **可证明边界** | 若必须 `unsafe`，需文档化不变量 + 单元测试 + Miri 覆盖 |

### 8.2 Crate 级 `#![deny(unsafe_code)]`

以下 crate **默认禁止** `unsafe`（在 `lib.rs` 顶部声明）：

- `giac-core`, `giac-parse`, `giac-simplify`, `giac-poly`
- `giac-linalg`, `giac-solve`, `giac-calculus`, `giac-ode`
- `giac-transform`, `giac-num`, `giac-special`, `giac-vector`
- `giac-prog`, `giac-emit`, `giac-groebner`
- `apps/giac-cli`, `tests/conformance`

**唯一例外**：`giac-ffi`（可选，用于 C ABI 对照或渐进替换）。该 crate 内：

```rust
// giac-ffi/src/lib.rs
#![warn(unsafe_op_in_unsafe_fn)]
// 不 deny：此处集中承载 FFI，但每个 unsafe 块需 SAFETY 注释
```

若未来为性能引入 `unsafe` 优化（如零拷贝缓冲区），应新建独立 crate（如 `giac-poly-unsafe`），**不得**污染上述 core crate。

### 8.3 依赖与 FFI 策略

| 场景 | 策略 |
|------|------|
| 大整数/模运算 | `num-bigint` / `rug`，不直接 `libc` 调 GMP |
| 数值线性代数 | `nalgebra` / `faer`（safe API） |
| 与 C++ giac 对照 | 仅 `giac-ffi` 暴露 `extern "C"`；conformance 测试优先纯 Rust 路径 |
| `rug` 内部 unsafe | 视为依赖边界，不在业务 crate 再包一层裸 FFI |

### 8.4 静态分析与 CI

每个 PR 对 workspace 执行（**均为必过，除标注可选**）：

```bash
cd giac-rs
cargo test --workspace          # P0 必过
cargo ci-clippy                 # P0 必过；见 .cargo/config.toml
cargo miri test -p giac-core    # P1：有 unsafe 或 FFI 变更时
```

**Clippy 安装**（任选其一，版本须匹配 `rustc`）：

- `rustup component add clippy`
- Ubuntu/Debian 系统 Rust：`sudo apt install rust-clippy`（勿 `cargo install clippy`）

- **Clippy**：`clippy.toml` + `[workspace.lints]`；详见 §8.4.1 与 [supplement §6.4 / §7](rust-migration-supplement.md#64-clippy-检查移植规范)
- **Miri**（可选）：至少覆盖 `giac-core`、`giac-poly`、`giac-simplify` 的属性测试与 golden 用例
- **`cargo deny`**（可选）：审计含 `unsafe` 的传递依赖

### 8.4.1 Clippy 检查规范（PR 必过）

| 项 | 要求 |
|----|------|
| **命令** | `cd giac-rs && cargo ci-clippy`（等价 `cargo clippy --workspace --all-targets -- -D warnings`） |
| **生效范围** | workspace 全部 member（`giac-core`、`giac-parse`、`giac-cli`、`tests/conformance` 等） |
| **配置文件** | `giac-rs/clippy.toml`（阈值、`disallowed-*`）；`Cargo.toml` `[workspace.lints]`（级别） |
| **新 crate** | `Cargo.toml` 须含 `[lints] workspace = true` |
| **unwrap/expect** | 生产路径见 [supplement §6](rust-migration-supplement.md#6-错误处理与-unwrap-规范)；crate 根 `cfg_attr(not(test), warn(...))` |
| **豁免** | `#[allow(clippy::…)]` 须 PR 说明；不得静默下调 workspace **deny** 规则 |

**deny 级（不可绕过）**：`unsafe_code`、`undocumented_unsafe_blocks`、`missing_safety_doc`、`mem_forget`，及 `clippy.toml` 中 `disallowed-types` / `disallowed-macros` / `disallowed-methods`。

**暂 allow（待代码收紧后可升为 warn）**：`indexing_slicing`、`integer_division`、`cast_*` 等——CAS/矩阵/解析中的合法模式；见 supplement §6.4 Lint 分级表。

### 8.5 代码评审检查项

合并前确认：

1. 新增 `unsafe` 是否仅在 `giac-ffi`（或已批准的 shim crate）？
2. 每个 `unsafe` 块是否有 `// SAFETY:` 说明前置条件与不变量？
3. 是否可用 safe 抽象替代（`Vec`、`Cow`、`Arc`）？MVP 阶段不引入 arena。
4. 是否避免 `gen` 式别名可变共享（优先 `Expr` 不可变 + 显式 `subst`/`map`）？
5. 错误路径是否返回 `Result`/`EvalError`，而非 `unwrap`/`panic`（对照 [supplement §6](rust-migration-supplement.md#6-错误处理与-unwrap-规范)）？
6. **`cargo ci-clippy` 是否全绿？** 新 crate 是否继承 `[lints] workspace = true`（对照 §8.4.1）？

### 8.6 与 C++ `gen` 的对照

| C++ `gen` 习惯 | Rust 纪律 |
|----------------|-----------|
| 引用计数 + 原地修改 | `Arc<Expr>` 或 `Expr` 拥有所有权；变换返回新 `Expr` |
| `ptr` 别名与 `ref_count` | 禁止；用类型系统表达共享（`Arc`）或借用（`&Expr`） |
| 全局 `context` 可变状态 | `Context` 显式传参；`thread_local` 仅用于兼容层 |
| 未定义行为容忍（如 factor segfault） | 返回 `Err`；测试断言错误而非崩溃 |

**结论**：按当前已覆盖 API 范围，**绝大多数代码应为 safe Rust**；`unsafe` 体量应 < 1% 行数，且集中在可选 FFI 边界。

---

## 9. 风险与决策

| 决策点 | 建议 | 理由 |
|--------|------|------|
| 全符号 vs 混合数值 | 混合 | 测试已覆盖 `evalf`/`fsolve`/`svd` 数值路径 |
| Groebner | 后置，先 greduce | `gbasis` 无稳定测试 |
| 积分算法 | 规则表 → Risch 子集 | `testintegrate` 67 条是主要难度 |
| 解析兼容性 | 先脚本子集，非完整 Xcas | 测试均为 `.in` 脚本，非 GUI session |
| 并行 | 后期 | `threaded.cc` 无专项测试 |
| 许可证 | GPL-3.0 继承 | giac 源码 GPL；Rust 重写需合规 |

---

## 10. 里程碑与完成定义

| 里程碑 | 标志 | 预估 |
|--------|------|------|
| M0 | `giac-cli` 跑通 `test_cas_basic` | 1 月 |
| M1 | 多项式 + 因式分解 golden 全过 | +1.5 月 |
| M2 | 线代 + 求解 | +1 月 |
| M3 | 微积分 check 套件 | +2 月 |
| M4 | 53 项 conformance 全绿 | +1.5 月 |
| **合计** | 可替换 headless giac CAS | **~7 月**（2 人） |

**完成定义**：

1. `tests/conformance` 对 `bin/` + `check/` 全部输入与 golden 等价（含规范化），不依赖 C++ `giac` 二进制
2. `cargo test --workspace` 与 **`cargo ci-clippy`** 持续全绿（工程门禁，见 [supplement §7](rust-migration-supplement.md#7-工程门禁)）

---

## 11. 建议的第一步

```bash
# 1. 初始化 workspace
cargo new giac-rs --workspace
cd giac-rs && mkdir -p crates/giac-core/src crates/giac-parse/src tests/conformance

# 2. 实现最小闭环
#    giac-parse: "1+2;" → Expr
#    giac-core:  eval → 3
#    giac-cli:   读 bin/test_cas_basic，打印每行结果

# 3. 接入 conformance
#    tests/conformance/src/harness.rs 读取 ../bin/* 与 ../giac/giac-1.5.0/check/*
#    双跑或单跑对比 golden
```

优先 port 顺序：**core → poly → simplify → linalg → solve → calculus**，与测试 ROI 和依赖关系一致。
