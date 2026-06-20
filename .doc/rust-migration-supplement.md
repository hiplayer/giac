# Giac → Rust 移植补充分析

本文档基于 [`rust-migration-plan.md`](rust-migration-plan.md)、[`functional-coverage.md`](functional-coverage.md)、[`test-inventory.md`](test-inventory.md) 及 giac-1.5.0 源码，补充迁移计划尚未展开的内容：**计划缺口**、**依赖库清单**、**文法解析架构**、**非图像 CAS 功能全景**。

---

## 1. 迁移计划深化文档索引

以下子文档已落盘，对应原 `rust-migration-plan.md` 中的待补充项。

| 主题 | 文档 |
|------|------|
| Token / 优先级 / 多模式 / `cas_setup` / `assume` | [`parser-token-map.md`](parser-token-map.md) |
| 测试 API → `at_*` → Rust crate（210 项） | [`builtin-api-map.md`](builtin-api-map.md) |
| 算法源文件分工（normalize/factor/partfrac 等） | [`module-division.md`](module-division.md) |
| Golden / `assert_equiv` / 随机 seed | [`conformance-testing.md`](conformance-testing.md) |
| 已知偏离登记 | [`known-divergences.md`](known-divergences.md) |

### 1.1 解析（摘要）

- **80+ token**、18 级优先级、隐式乘 `T_IMPMULT` → 见 `parser-token-map.md` §2–3。
- **MVP 模式：** 仅 `xcas_mode(0)`；TI（`xcas_mode==3`）、Maple、RPN、Python 兼容 **不做**。
- **Harness 注入：** check 首行 `cas_setup(0,0,0,1,0,[1e-10,1e-17],12,[1,50,0,25],0,0,0),xcas_mode(0)` → 12 字段 `Context` 映射表见 `parser-token-map.md` §5–6。

### 1.2 类型与模块（摘要）

- giac 全库 **~1852** 个 `at_*`；MVP 仅需 **~210** 个测试覆盖 API → `builtin-api-map.md`。
- **`normalize.cc` / `factor.cc` / `partfrac.cc` / `integrate.cc` 是 profiler 可执行文件**，非库；算法在 `usual.cc`、`ezgcd.cc`、`sym2poly.cc`、`intg.cc` → `module-division.md` §1–4。

### 1.3 测试（摘要）

- **53 CTest** = 41 bin + 10 check + 1 tommath + 1 segfault → `conformance-testing.md` §1。
- **`assert_equiv`：** `normal(sub(a,b))` 为零；积分/求解用 `diff` 还原 → §3。
- **随机数：** tinymt32 + Box-Muller；测试验类型/有限性，可选 `srand` 固定 seed → §4。
- **DIV-001：** factor segfault → `known-divergences.md`。

### 1.4 工程（摘要）

| 项 | 说明 |
|----|------|
| flex/bison | 源：`input_lexer.ll`、`input_parser.yy`；CMake **用预生成 .cc**，改 grammar 需本地 regen → `parser-token-map.md` §1 |
| libtommath | 静态链入；→ `num-bigint` / `rug` |
| gettext | `translate.cc`；integrate/geo 测试设 `LANG=en_US.UTF-8`；MVP 可硬编码英文 |
| FLTK / gl2ps | 可选 GUI；headless 不移植 |
| 并行 | `threaded.cc` 无测试；Rust `rayon` 后置 |
| GPL-3.0 | Rust 重写需 NOTICE + 源码提供说明（待 legal 评审） |

---

## 2. 移植依赖库

分两层：**giac C++ 现有依赖**（理解移植源）与 **Rust 侧建议依赖**（目标实现）。

### 2.1 giac C++ 构建依赖（当前 upstream）

| 依赖 | 类型 | 主要功能 | 移植策略 |
|------|------|----------|----------|
| **libtommath** 0.39 | 静态库 | 大整数 `mp_int`；GMP 不可用时的整数后端 | → `num-bigint` 或 `rug::Integer` |
| **flex** 2.5.x | 构建时 | 生成 `input_lexer.cc` 词法分析器 | → `logos` |
| **bison** | 构建时 | 生成 `input_parser.cc` LALR 语法 | → `lalrpop` 或 `chumsky` |
| **C++ 标准库** | 系统 | `vector`、`map`、异常、iostream | → Rust 标准库 + `thiserror` |
| **librt** | 系统 (Linux) | 时钟等 | 按需 |
| **FLTK** | 可选 | GUI：`icas.cc`、`Graph*.cc` 等 | **不做**（headless CAS） |
| **gl2ps** | 源码编入 | PostScript/PDF 图形导出 | **不做** |
| **CoCoA** | 外部（已排除） | F5 Groebner；`cocoa.cc` 在 CAS 构建中 stub | 自研 `greduce`；`gbasis` 后置 |
| **tinymt32** | 内嵌源码 | 梅森旋转随机数；`randNorm` 等 | → `rand` + `rand_distr` |
| **TmpFGLM / TmpLESystemSolver** | 内嵌 C++ | 原生 Groebner/FGLM 辅助 | 参考算法，Rust 自研 |
| **gettext** | 可选 | 错误/帮助信息 i18n | MVP 英文/中文硬编码或 `rust-i18n` 后置 |

> 注：giac 历史上亦支持 GMP/MPFR（`config.h`）；本仓库 CMake 构建以 **libtommath** 为主，浮点路径混用 `double` 与 `_REAL`（MPFR 类）。

### 2.2 Rust 移植目标依赖（建议 workspace）

#### 2.2.1 核心数学与数据结构

| Crate | 功能 | 对标的 giac 模块 |
|-------|------|------------------|
| `num-bigint` | 任意精度整数 | `_ZINT`、libtommath、`ifactor` |
| `num-rational` | 精确有理数 | 数字 `_FRAC`、`Ratio` |
| `num-traits` | 数值 trait 抽象 | 泛型系数环 |
| `rug`（可选） | GMP/MPFR 绑定；模运算、任意精度浮点 | `_REAL`、`evalf`、`fsolve` 高精度路径 |
| `num-complex`（可选） | 复数 `f64` 快路径 | 数值 `evalf` |

#### 2.2.2 线性代数与数值

| Crate | 功能 | 对标的 giac API |
|-------|------|-----------------|
| `nalgebra` 或 **faer** | `f64` 矩阵 LU/QR/SVD | `lu`, `qr`, `svd`, `fsolve`, `newton` |
| 自研符号矩阵 | `Expr::Matrix` 上的 `rref`/`det`/`linsolve` | `lin.cc`, `vecteur.cc`, `gauss.cc` |

#### 2.2.3 解析与 CLI

| Crate | 功能 |
|-------|------|
| `logos` | 词法分析（替代 flex） |
| `lalrpop` 或 `chumsky` | 语法分析（替代 bison） |
| `clap` | `giac-cli` 命令行 |
| `thiserror` / `anyhow` | 错误类型 |

#### 2.2.4 测试与质量

| Crate | 功能 |
|-------|------|
| `insta` | golden 快照 |
| `serde` / `serde_json` | 测试 fixture、偏离记录 |
| `proptest` | 属性测试（gcd、交换律、等价性） |
| **`cargo clippy`**（**必过**） | 静态分析；`cargo ci-clippy`，见 §6.4、§7 |
| `cargo-tarpaulin`（CI） | 覆盖率 |
| `cargo-deny`（可选） | 依赖审计 |

#### 2.2.5 随机与特殊函数

| Crate | 功能 |
|-------|------|
| `rand` + `rand_distr` | `randNorm`, `randchisquare` 等 |
| 自研或 `special`（按需） | `gamma`, `Beta`, `zeta`, 正交多项式 |

#### 2.2.6 后置 / 可选

| Crate | 功能 |
|-------|------|
| `rayon` | 并行（`threaded.cc` 无测试，后置） |
| `bumpalo` | arena（profiling 后评估） |
| `giac-ffi` + `libc` | C ABI 对照层（唯一允许大量 `unsafe` 的 crate） |

### 2.3 依赖选型决策摘要

| 决策点 | 建议 | 理由 |
|--------|------|------|
| 大整数 | MVP `num-bigint`；模/GCD 热点可换 `rug` | 纯 Rust、Miri 友好 |
| 浮点 | 符号层 `PrecFloat`（`rug::Float`）；快路径 `f64` | 对齐 `evalf`/`float2rational` 测试 |
| 线性代数 | 符号自研 + 数值 `faer` | 测试已覆盖双路径 |
| 解析 | `logos` + `lalrpop` | 可维护、类型安全 |
| GMP 直连 | 避免（除 `rug` 边界） | §8 unsafe 纪律 |

---

## 3. giac 文法解析架构

### 3.1 总体流水线

```
输入字符串
  → python2xcas() 等预处理（gen.cc::try_parse）
  → set_lexer_string() 绑定 flex scanner + context
  → giac_yyparse()  bison 语法分析
  → parsed_gen()    存入 thread-local 解析结果
  → eval / 脚本执行
```

**源文件：**

| 文件 | 工具 | 职责 |
|------|------|------|
| `input_lexer.cc` | flex 2.5 生成 | 词法：数字、标识符、运算符、容器定界符 |
| `input_parser.yy` → `input_parser.cc` | bison 生成 | 语法：表达式、语句、程序 |
| `lexer.h` | flex 头 | scanner 类型、`yylex` 声明 |
| `gen.cc` | — | `try_parse` / `gen(const string&)` 入口 |
| `global.cc` | — | `context`：模式标志、精度、假设 |
| `prog.cc` | — | 程序执行、`for`/`while` 运行时 |

语法栈深度：`YYINITDEPTH=4000`，`YYMAXDEPTH=20000`（桌面）；嵌入式平台更小。

### 3.2 词法层（flex）

词法器为 **reentrant**（`yyscan_t`），通过 `YY_EXTRA_TYPE` 携带 `const context *`，使同一 token 在不同模式下产生不同 token 类型。

**主要 token 类别：**

| 类别 | 代表 token | 说明 |
|------|------------|------|
| 字面量 | `T_NUMBER`, `T_LITERAL`, `T_STRING` | 整数/浮点、π/i/∞、字符串 |
| 标识符 | `T_SYMBOL`, `T_EXPRESSION` | 变量与函数名 |
| 算术 | `T_PLUS`, `T_MOINS`, `T_FOIS`, `T_DIV`, `T_MOD`, `T_POW`, `T_SQ` | `^` 与 `**` |
| 关系 | `T_TEST_EQUAL`, `T_EQUAL` | `==`, `=`（方程） |
| 集合 | `T_UNION`, `T_INTERSECT`, `T_MINUS`, `T_INTERVAL` | 区间与集合运算 |
| 容器 | `T_VECT_BEGIN/END`, `T_VECT_DISPATCH`, `T_INDEX_BEGIN` | `[ ]`, `{ }`, 下标 |
| 矩阵 | `T_MATRICE_BEGIN/END` | `[[ ]]` |
| 特殊 | `T_ROOTOF_BEGIN/END`, `T_SPOLY1_*`, `T_POLY1_*` | 代数数、稀疏/稠密多项式 |
| 语句 | `T_SEMI`, `T_AFFECT`, `TI_STO` | `;`, `:=`, TI 存储 |
| 控制流 | `T_FOR`, `T_WHILE`, `T_IF`, `T_RETURN`, … | 完整 Xcas 程序语法 |
| 模式 | `TI_*`, `T_RPN_*`, `T_MAPLELIB` | TI-89 / RPN / Maple 方言 |

**容器 dispatch（`T_VECT_DISPATCH`）** 由词法器设置 subtype，语法层构造不同 `gen` 向量类型：

- `_SEQ__VECT` — 序列 `[a,b,c]`
- `_LIST__VECT` — 列表
- `_SET__VECT` — 集合
- `_MATRIX__VECT` — 矩阵
- `_POLY1__VECT` — 一元多项式系数向量
- `_ASSUME__VECT` — 假设
- `_PNT__VECT`, `_CURVE__VECT`, … — 几何/绘图对象（非 MVP）

**模式相关行为示例：**

- `xcas_mode==3`：`;` → `TI_SEMI`；`: ` → `TI_DEUXPOINTS`
- `xcas_mode==0` 且无 `i:=…`：`i` 解析为 `sqrt(-1)`（解析后二次扫描）
- `rpn_mode`：反引号等 RPN 定界符
- `python_compat`：圆括号列表 → `python_list`

### 3.3 语法层（bison）

- **起始符号：** `input` → `correct_input`（一条或多条 `exp;`）
- **语义值类型：** `YYSTYPE = gen`（直接在语法动作中构造 AST）
- **表达式 `exp`：** 运算符、函数调用、下标、程序、控制流均在此层
- **语句序列：** `correct_input` 递归合并为 `_SEQ__VECT`

**运算符优先级（从低到高，节选）：**

```
T_AFFECT (:>)  >  TI_SEMI  >  T_VIRGULE  >  T_AND_OP  >  T_EQUAL  >
T_TEST_EQUAL  >  T_UNION/INTERSECT  >  T_PLUS/MOINS  >
T_FOIS/T_IMPMULT  >  T_DIV  >  T_MOD  >  T_POW  >  T_FACTORIAL  >
T_UNARY_OP  >  T_COMPOSE
```

**关键语法构造 → AST：**

| 语法 | AST 形式 |
|------|----------|
| `a+b` | `symbolic(at_plus, [a,b])` |
| `a^b` | `symb_pow(a,b)` |
| `a mod b` | `normalmod(a,b)` 或 `symbolic(mod,...)` |
| `f(x)` | `symbolic(at_of, [f, args])` |
| `a:=b` | `symb_sto(b, a)` |
| `x' ` | `symbolic(at_derive, x)` |
| `rootof(...)` | `algebraic_EXTension(...)` → `_EXT` |
| `for i from a to b do ...` | `symbolic(at_for, [init, cond, step, bloc])` |
| `[[1,2],[3,4]]` | `gen(matrix, _MATRIX__VECT)` |

**解析入口（`gen.cc`）：**

1. 去前导空白、`python2xcas` 转换
2. `set_lexer_string(s, scanner, contextptr)`
3. `giac_yyparse(scanner)` → 填充 `parsed_gen`
4. 若 `xcas_mode==0`，检查 `i` 是否被用户绑定；否则重解析并将 `i` 视为 `sqrt(-1)`

### 3.4 Rust 移植要点（`giac-parse`）

MVP 脚本子集应对齐 **bin/** + **check/** 测试，而非完整 Xcas：

**必须实现：**

- 行尾 `;`、`:=`、表达式/赋值语句
- 算术、`^`/`**`、`mod`、关系运算符
- 函数调用 `f(x,y)`、嵌套
- 向量/矩阵/序列字面量
- 下标 `a[i]`、`a[i,j]`
- `for`/`while`/`sum`/`product`（`test_prog`）
- 特殊：`rootof`、`assume`（测试中出现）

**可推迟：**

- TI / Maple / RPN 全套 token
- `switch`/`try_catch`/`piecewise` 完整方言
- 单位换算 `T_UNIT`、Logo 语法
- 帮助 `T_HELP`、spreadsheet `$` 语法

---

## 4. 非图像 CAS 功能全景

以下按 **源码模块 → 功能域** 列举，**排除** GUI/绘图/3D/兼容层；标注测试覆盖与 Rust crate 映射。

图例：**✓** = 有专项或 check 回归；**△** = 仅 flanex 等综合用例；**✗** = 无稳定测试 / 计划不做

### 4.1 核心运行时

| 模块 | 功能 | 测试 | Rust crate |
|------|------|------|------------|
| `gen.cc` | `gen` 类型、运算、解析入口、打印 | ✓ | `giac-core` |
| `symbolic.cc` | 符号对象 `symbolic{sommet,feuille}` | ✓ | `giac-core` |
| `index.cc` | 单项序、指数向量 | ✓ | `giac-poly` |
| `identificateur.cc` | 标识符表 | ✓ | `giac-core` |
| `first.cc` | 初始化、常量 | ✓ | `giac-core` |
| `global.cc` | 全局 context、`cas_setup` | △ | `giac-core` |
| `unary.cc` | 一元函数 dispatch | ✓ | `giac-core` |
| `threaded.cc` | 线程参数 | ✗ | 后置 |

### 4.2 解析与程序

| 模块 | 功能 | 测试 | Rust crate |
|------|------|------|------------|
| `input_lexer.cc` | 词法 | ✓（全脚本） | `giac-parse` |
| `input_parser.yy` | 语法 | ✓ | `giac-parse` |
| `prog.cc` | 程序执行、控制流 | ✓ `test_prog` | `giac-prog` |
| `subst.cc` | 替换 `subst` | ✓ | `giac-prog` |
| `find_global_var.cc` | 全局变量查找 | △ | `giac-prog` |

### 4.3 化简与基础运算

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `usual.cc` | 三角/双曲化简、expand、normal | `normal`, `ratnormal`, `texpand`, `tlin`, `halftan`, `trig2exp`, … | ✓ |
| `normalize.cc` | 非递归规范化 | `non_recursive_normal` | ✓ check |
| `misc.cc` | 杂项化简、概率、特殊函数、向量算子 | `simplify`, `Beta`, `gamma`, `evalf`, `fsolve`, … | ✓ |
| `softmath.cc` | 软数学辅助 | △ | `giac-simplify` / `giac-special` |

### 4.4 多项式与代数

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `sym2poly.cc` | 表达式 ↔ 多项式 | `e2r`, `r2e`, `valuation` | ✓ |
| `ezgcd.cc` | GCD/LCM/因式分解 | `gcd`, `lcm`, `factor`, `egcd`, `resultant` | ✓ |
| `modpoly.cc` | 模多项式 | `chinrem`, `smod`, 模 `gcd`/`rref` | ✓ |
| `modfactor.cc` | 模因式分解 | △ | ✓ modular |
| `gausspol.cc` | 多项式高斯消元 | △ | ✓ |
| `gauss.cc` | 多项式矩阵 | `gauss`, `content`, `quo`, `rem` | ✓ |
| `factor.cc` | 因式分解入口 | `factor`, `factors` | ✓ |
| `partfrac.cc` | 部分分式 | `partfrac`, `propfrac` | ✓ |
| `alg_ext.cc` | 代数扩域 `rootof` | `rootof`, 扩域算术 | ✓ |
| `TmpFGLM.C` | FGLM Groebner | `greduce`（部分） | △ |
| `TmpLESystemSolver.C` | 线性方程组求解辅助 | △ | `giac-groebner` |
| `cocoa.cc` | CoCoA Groebner | `gbasis`（stub） | ✗ |
| `sparse.cc` | 稀疏矩阵 | ✗ | 不做 |

### 4.5 线性代数

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `lin.cc` | 符号矩阵 | `rref`, `det`, `rank`, `ker`, `image`, `jordan`, `egv`, `charpoly`, … | ✓ |
| `vecteur.cc` | 向量运算 | `cross`, `dot`, `tran`, `hadamard`, … | ✓ |

### 4.6 方程求解与根隔离

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `solve.cc` | 代数/数值求解 | `solve`, `linsolve`, `froots`, `realroot`, `fsolve` | ✓ |
| `csturm.cc` | Sturm 序列 | `sturm`, `sturmab` | ✓ |

### 4.7 微积分

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `derive.cc` | 微分 | `diff`, `derive` | ✓ |
| `intg.cc` | 积分主模块 | `integrate`, `int` | ✓ |
| `integrate.cc` | 积分辅助 | △ | ✓ |
| `risch.cc` | Risch 算法 | `risch` | ✓ |
| `series.cc` | 级数 | `series`, `taylor`, `limit` | ✓ |
| `intgab.cc` | 拉普拉斯/傅里叶 | `laplace`, `ilaplace`, `fourier_an`, `fourier_cn` | ✓ |

### 4.8 微分方程

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `desolve.cc` | ODE | `desolve` | ✓ |

### 4.9 数论与组合

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `ifactor.cc` | 整数分解、素数 | `ifactors`, `isprime`, `nextprime`, `divisors`, `euler`, `nCr`, `nPr`, `bernoulli` | ✓ |
| `permu.cc` | 置换 | `permu2cycles`, `cycles2permu`, `permuorder`, … | ✓ |

### 4.10 特殊函数与四元数

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `moyal.cc` | Moyal 分布等 | `moyal`, `UTPN` | ✓ |
| `quater.cc` | 四元数 | `quaternion` | ✓ |
| `misc.cc` | 正交多项式、随机 | `legendre`, `hermite`, `laguerre`, `randNorm`, … | ✓ |

### 4.11 几何（符号，非绘图 GUI）

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `plot.cc` | 符号几何构造、部分绘图 | `triangle`, `mkisom` 相关 | ✓ geo/flanex |
| `isom.cc` | 几何同构 | `mkisom` | ✓ |
| `plot3d.cc` | 3D 几何（CAS 构建 stub） | `plan`, `sphere` | ✗ stub |

> `testgeo` 的 133 条 golden 含几何构造与代数化简，**非** FLTK 绘图；Rust 需实现符号几何子集或等价代数路径。

### 4.12 变换与向量分析

| 功能域 | 代表 API | 测试 | Rust crate |
|--------|----------|------|------------|
| 拉普拉斯/傅里叶 | `laplace`, `ilaplace`, `fourier_*` | ✓ | `giac-transform` |
| 向量微积分 | `grad`, `curl`, `div`, `laplacian`, `potential`, `vpotential`, `hessian` | ✓ | `giac-vector` |

### 4.13 导出

| 模块 | 功能 | 代表 API | 测试 |
|------|------|----------|------|
| `tex.cc` | LaTeX | `latex` | ✓ |
| `mathml.cc` | MathML | `mathml` | ✓ |
| `hevea2mml.cc` | Hevea→MathML | ✗ | 不做 |

### 4.14 辅助 / 不移植

| 模块 | 功能 | 移植 |
|------|------|------|
| `help.cc`, `translate.cc` | 帮助与 i18n | 后置 |
| `maple.cc`, `ti89.cc`, `rpn.cc`, `pari.cc` | 兼容方言 | **不做** |
| `signalprocessing.cc`, `graphtheory.cc`, `optimization.cc`, `lpsolve.cc` | 信号/图论/优化/LP | 无测试，不做 |
| `luabridge.cc` | Lua 嵌入 | 不做 |
| `renee.cc` | 内部实验 | 不做 |
| `tinymt32.cc` | 随机数引擎 | → `rand` |

### 4.15 按测试规格归纳的功能域（MVP 边界）

与 `functional-coverage.md` 一致，**220+ API** 归入 **19 个功能域**（不含 GUI）：

1. 基础设施 — 解析、Context、`eval`、`sto`、`quote`
2. 基础算术与复数 — `abs`, `arg`, `conj`, `re`, `im`, `sign`, `gcd`
3. 化简与三角 — `normal`, `ratnormal`, `texpand`, `tlin`, …
4. 多项式 — `factor`, `gcd`, `quo`, `rem`, `resultant`, `roots`, `e2r`, …
5. 模运算 — `chinrem`, `smod`, 模多项式运算
6. 线性代数 — `rref`, `det`, `linsolve`, `svd`, `lu`, `qr`, …
7. 方程与 Sturm — `solve`, `sturm`, `realroot`
8. 微积分 — `diff`, `integrate`, `limit`, `series`, `partfrac`, `risch`
9. ODE — `desolve`
10. 级数/极限/分式 — 与微积分重叠，独立 golden
11. 数论 — `ifactors`, `isprime`, `nCr`, …
12. 排列 — `permu*`, `cycles*`
13. 数值 — `evalf`, `fsolve`, `float2rational`
14. 变换 — `laplace`, `fourier_*`
15. 向量分析 — `curl`, `div`, `grad`, `potential`, …
16. 特殊函数/概率 — `gamma`, `Beta`, `zeta`, `rand*`
17. 程序与替换 — `subst`, `sum`, `product`, `for`, `while`
18. 导出 — `latex`, `mathml`
19. Groebner — **仅** `greduce`（`gbasis` 不做）
20. 几何符号 — `triangle`, 构造线（`testgeo`）
21. 同构 — `mkisom`
22. 四元数 — `quaternion`

### 4.16 明确排除（非图像但也不移植）

- GUI 会话：`icas.cc`, `Xcas1.cc`, `Editeur.cc`, `kdisplay.cc`, …
- 2D/3D 绘图引擎：`Graph.cc`, `Graph3d.cc`, `gl2ps.c`, `opengl.cc`
- 电子表格 UI：`Tableur.cc`
- HTML 会话：`cas2html.cc`
- Maple/TI/RPN/Pari 兼容层
- CoCoA `gbasis` 完整实现
- 稀疏矩阵 `sparse.cc`（无测试）

---

## 5. 参考路径

| 路径 | 说明 |
|------|------|
| `giac/giac-1.5.0/src/input_parser.yy` | 语法定义 |
| `giac/giac-1.5.0/src/input_lexer.cc` | 词法（flex 生成） |
| `giac/giac-1.5.0/src/gen.cc` | 解析入口 ~11681 行起 |
| `giac/CMakeLists.txt` | 构建源列表与测试注册 |
| `giac/cas_excluded_stubs.cc` | cocoa/plot3d stub |
| `.doc/functional-coverage.md` | API ↔ 测试映射 |
| `.doc/rust-migration-plan.md` | 主迁移方案 |
| 本文 §6 | 错误处理与 unwrap 规范 |
| 本文 §7 | **工程门禁（含 Clippy 必过）** |
| `giac-rs/clippy.toml` | Clippy 阈值与安全 deny-list |

---

## 6. 错误处理与 unwrap 规范

Rust 移植与 C++ `gen` 的重要差异之一：用户输入与算法失败应返回 `Result`，而不是 panic 或把错误编码进表达式字符串。

### 6.1 分层规则

| 层级 | 策略 |
|------|------|
| **公开 API**（`eval`, `parse_program`, `exec_stmt`, CLI） | 必须 `Result<_, EvalError>` / `ParseError`；调用方用 `?` |
| **库内部求值/化简/解析** | 优先 `?` 传播 `EvalError`；禁止对用户可控输入 `unwrap` |
| **可证明不变量** | 极少数 `expect("…")` 可接受，须注释为何不会失败 |
| **`#[cfg(test)]` / `tests/`** | `unwrap()` / `expect()` 允许；测试边界可 panic |
| **Conformance harness** | 内部 `Result` + 带上下文的 `map_err`；测试函数返回 `Result<(), String>` |

### 6.2 数值转换

`BigInt` → 有界整数（指数、矩阵维数等）统一经 `giac-core::num_util`：

```rust
use crate::num_util::{bigint_to_nonneg_u32, bigint_to_u32_abs};

let e_u = bigint_to_nonneg_u32(e)?;   // 非负 u32
let e_u = bigint_to_u32_abs(e)?;      // 绝对值 u32（负指数）
```

禁止在生产路径写 `n.to_string().parse().unwrap()`。

### 6.3 示例

```rust
// 好：公开入口
pub fn eval(expr: &Expr, ctx: &Context) -> Result<ExprArc, EvalError> { … }

// 好：内部传播
let p = expr_to_poly(a.as_ref())?;
let e_u = bigint_to_nonneg_u32(e)?;

// 好：测试
let r = eval(e.as_ref(), &ctx).unwrap();

// 差：生产路径
let n = as_int(x).unwrap();
```

### 6.4 Clippy 检查（移植规范）

Clippy 自 **Phase 0 起为必过项**，与 `cargo test --workspace` 并列；未通过不得合并 PR。完整门禁见 [§7](#7-工程门禁)。

**配置文件**（`giac-rs/`）：

| 文件 | 作用 |
|------|------|
| `clippy.toml` | MSRV、复杂度阈值、数学单字母绑定、`disallowed-*` 硬禁止列表 |
| `Cargo.toml` `[workspace.lints]` | 全 workspace lint 级别（`unsafe_code = deny`、浮点/索引告警等） |
| `.cargo/config.toml` | `cargo ci-clippy` 别名（`-D warnings`） |

**生产代码** `unwrap` / `expect`：`giac-core` / `giac-parse` / `giac-simplify` crate 根：

```rust
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]
```

`cargo clippy --lib` 对生产代码告警；`cargo test` / 集成测试编译带 `cfg(test)` 时不告警。

**`clippy.toml` 要点**（CAS 特性）：

- `msrv = "1.75.0"` 与工具链一致
- 放宽 `cognitive-complexity` / `too-many-lines`（`eval` 等算法体量大）
- `allowed-idents-below-min-chars` 允许 `x`, `y`, `i` 等数学符号
- `disallowed-types`：`CString` / `CStr` / `NonNull`（仅 `giac-ffi` 可用）
- `disallowed-macros`：`todo!` / `unimplemented!` → 应返回 `EvalError`
- `disallowed-methods`：`abort` / `unreachable_unchecked`
- **不在** `disallowed-methods` 中封 `unwrap`（测试边界允许，见 §6.1）
- Clippy **1.75** 不支持 `check-private-items`（勿写）；`single_char_binding_names` 亦非 1.75 lint 名

**Lint 分级**（`Cargo.toml` `[workspace.lints]` + `cargo ci-clippy` = `-D warnings`）：

| 级别 | 示例 |
|------|------|
| **deny** | `unsafe_code`, `undocumented_unsafe_blocks`, `mem_forget` |
| **allow（暂）** | `indexing_slicing`, `integer_division`, `cast_*` — CAS/矩阵/解析大量合法索引与整除，待逐步收紧 |
| **crate 根 warn** | `unwrap_used`, `expect_used`（`cfg_attr(not(test))`） |

**CI**（需先安装 Clippy，**不要** `cargo install clippy`——crates.io 上的包已废弃）：

| 环境 | 安装 |
|------|------|
| **rustup** | `rustup component add clippy` |
| **Ubuntu/Debian**（系统 `rustc` 1.75，如本机 `/usr/bin/rustc`） | `sudo apt install rust-clippy` |
| 版本须与 `rustc --version` 一致；勿装 `rust-1.8x-clippy` 若 compiler 仍是 1.75 |

```bash
cd giac-rs
cargo ci-clippy    # 或: cargo clippy --workspace --all-targets -- -D warnings
```

**PR 约束**：

- 新增/修改 Rust 代码后本地跑通 `cargo ci-clippy`（退出码 0）
- 新 crate 须在 `Cargo.toml` 加 `[lints] workspace = true`
- 禁止为消告警随意 `#[allow(clippy::…)]`；确需豁免须在 PR 说明理由
- 勿改 `clippy.toml` / `[workspace.lints]` 放宽 **deny** 级规则（`unsafe_code` 等）除非 ADR 记录

### 6.5 Conformance harness 模式

```rust
fn run_script_lines(lines: &[&str]) -> Result<Vec<String>, String> {
    let stmts = parse_program(&input, &ctx).map_err(|e| format!("parse: {e}"))?;
    for (idx, stmt) in stmts.iter().enumerate() {
        exec_stmt(stmt, &mut ctx).map_err(|e| format!("eval stmt {idx}: {e}"))?;
    }
    Ok(out)
}

#[test]
fn cas_tst_first_25_batch() -> Result<(), String> {
    let got = run_script_lines(&input_slice)?;
    assert!(passed >= 23);
    Ok(())
}
```

失败时测试输出含 `eval stmt 12: TypeError(...)`，便于定位 golden 回归。

---

## 7. 工程门禁

移植代码合并与阶段验收时，**除 golden 外**须满足：

| 优先级 | 检查 | 命令 | 通过标准 |
|--------|------|------|----------|
| **P0** | 测试 | `cargo test-timeout`（`giac-rs/`；见 [conformance-testing §5.1](conformance-testing.md#51-单测超时推荐-cargo-test-timeout)） | 全绿 |
| **P0** | **Clippy** | `cargo ci-clippy` | 退出码 0（`-D warnings`） |
| P1 | Miri（有 `unsafe` 时） | `cargo miri test -p giac-core` | 无 UB |
| P2 | 覆盖率（核心 crate） | `cargo tarpaulin -p giac-core` | ≥70%（目标，非阻塞 MVP） |

**Clippy 安装**（版本须与 `rustc --version` 一致；**勿** `cargo install clippy`）：

| 环境 | 命令 |
|------|------|
| rustup | `rustup component add clippy` |
| Ubuntu/Debian 系统 Rust 1.75 | `sudo apt install rust-clippy` |

**最小本地自检**（提交 PR 前）：

```bash
cd giac-rs
cargo test-timeout          # 推荐；并行 + 单测超时（需 cargo-nextest）
cargo ci-clippy
```

**测试通过后再做（算法 crate 必做）：** 复审 diff — 临时形状匹配是否净减少、新增/变更 `fn` 是否已标 tier 并更新 `.doc/*-api-stability.md`（见 [algorithm-expr-api.md §7.2](algorithm-expr-api.md#72-测试通过后提交--合入前复审)）。

全量无超时包装（不推荐日常/CI）：`cargo test --workspace`。说明见 [conformance-testing.md §5.1](conformance-testing.md#51-单测超时推荐-cargo-test-timeout) 与 `giac-rs/README.md`。

配置细节：§6.4；unsafe 纪律：[`rust-migration-plan.md` §8](rust-migration-plan.md#8-内存安全与-unsafe-纪律)。
