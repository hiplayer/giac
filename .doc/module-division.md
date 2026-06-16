# 算法模块分工（C++ 源文件）

澄清迁移计划中易混淆的文件：**同名 `.cc` 可能是独立 profiler 可执行文件，而非库模块**。

---

## 1. _profiler 可执行文件（非库，勿移植）_

以下文件仅含 `main()`，用于本地性能剖析，**不在** `giac` CMake 源列表中：

| 文件 | 用途 | 实际算法所在 |
|------|------|--------------|
| `normalize.cc` | 对 stdin 表达式跑 `normal()` | `usual.cc`, `sym2poly.cc` |
| `partfrac.cc` | 对 stdin 跑 `partfrac()` | `sym2poly.cc` |
| `factor.cc` | 对 stdin 跑 `factor()` | `ezgcd.cc`, `modfactor.cc` |
| `integrate.cc` | 对 stdin 跑 `integrate()` | `intg.cc`, `risch.cc` |

---

## 2. 化简与规范化

| 文件 | 职责 | 代表 API | Rust crate |
|------|------|----------|------------|
| **`usual.cc`** | 主化简：`normal`, `ratnormal`, `expand`, `factor` 入口、三角化简、`assume` | `normal`, `texpand`, `tlin`, `halftan`, `chinrem` | `giac-simplify` |
| **`sym2poly.cc`** | 表达式↔多项式、`partfrac` 算法、`non_recursive_normal` 注册 | `e2r`, `r2e`, `partfrac`, `non_recursive_normal` | `giac-poly`, `giac-simplify` |
| **`misc.cc`** | 杂项：`simplify`, `reorder`, `canonical_form`, `evalf`, 向量算子 | `simplify`, `reorder`, `evalf`, `curl` | `giac-simplify`, `giac-special` |
| **`symbolic.cc`** | 符号对象运算、quoted 函数 | 内部 | `giac-core` |

**要点：**

- `non_recursive_normal` 在 `sym2poly.cc` 注册，实现委托给 `normal()`（`usual.cc` 链）。
- `check/testnormalize`（188 条）测的是 **non_recursive_normal**，不是独立算法。

---

## 3. 多项式与因式分解

| 文件 | 职责 | 代表 API | Rust crate |
|------|------|----------|------------|
| **`ezgcd.cc`** | GCD/LCM、一元/多元 gcd、**因式分解核心** | `gcd`, `lcm`, `factor`, `egcd`, `resultant` | `giac-poly` |
| **`modfactor.cc`** | 模因式分解 | 模 `factor` | `giac-poly` |
| **`modpoly.cc`** | 模多项式环运算 | `chinrem`, `smod`, 模 gcd/rref | `giac-poly` |
| **`gausspol.cc`** | 多项式表示、单项序底层 | 内部 | `giac-poly` |
| **`gauss.cc`** | 多项式矩阵、`gauss` 消元 | `gauss`, `content`, `quo`, `rem` | `giac-poly`, `giac-linalg` |
| **`ifactor.cc`** | 整数分解、`abcuv` 等 | `ifactors`, `isprime`, `abcuv` | `giac-num`, `giac-poly` |
| **`alg_ext.cc`** | 代数扩域 `rootof` | `rootof`, 扩域算术 | `giac-core` |
| **`TmpFGLM.C`** | Groebner FGLM 辅助 | `greduce` 原生路径 | `giac-groebner` |
| **`TmpLESystemSolver.C`** | 线性方程组（Groebner 子步骤） | 内部 | `giac-groebner` |
| **`cocoa.cc`** | CoCoA F5（**CAS 构建 stub**） | `gbasis` 不可用 | 不做 |

**要点：**

- **`factor.cc` 不是库**；`factor` API 在 `ezgcd.cc` + `usual.cc` 入口。
- `gbasis` 依赖 CoCoA，当前 stub；MVP 只做 `greduce`（`bin/test_groebner`）。

---

## 4. 微积分

| 文件 | 职责 | 代表 API | Rust crate |
|------|------|----------|------------|
| **`derive.cc`** | 微分 | `diff`, `derive` | `giac-calculus` |
| **`intg.cc`** | **积分主入口**、规则表、定积分 | `integrate`, `int`, `ibp`, `ibpu` | `giac-calculus` |
| **`risch.cc`** | Risch 扩展积分 | `risch` | `giac-calculus` |
| **`series.cc`** | 级数与极限 | `series`, `taylor`, `limit` | `giac-calculus` |
| **`intgab.cc`** | 拉普拉斯、傅里叶 | `laplace`, `ilaplace`, `fourier_*` | `giac-transform` |
| **`integrate.cc`** | profiler 可执行文件 | — | 不移植 |

**要点：**

- 所有积分测试（`testintegrate` 67 条）入口是 **`intg.cc::_integrate`**，按需调用 `risch.cc`。
- `partfrac` 在 **`sym2poly.cc`**，被积分/heaviside 等调用。

---

## 5. 方程、线代、其他

| 文件 | 职责 | Rust crate |
|------|------|------------|
| `solve.cc` | `solve`, `linsolve`, `fsolve` | `giac-solve` |
| `csturm.cc` | `sturm`, `sturmab`, `realroot` | `giac-solve` |
| `desolve.cc` | `desolve` | `giac-ode` |
| `lin.cc`, `vecteur.cc` | 符号/数值矩阵 | `giac-linalg` |
| `subst.cc`, `prog.cc` | 替换、程序、`cas_setup` | `giac-prog`, `giac-core` |
| `permu.cc`, `ifactor.cc` | 置换、数论 | `giac-num` |
| `moyal.cc`, `quater.cc` | 特殊函数、四元数 | `giac-special` |
| `tex.cc`, `mathml.cc` | 导出 | `giac-emit` |
| `isom.cc`, `plot.cc`（符号几何） | 几何 | `giac-geo`（或 core 子集） |

---

## 6. `at_*` 注册规模

| 范围 | 数量 | 说明 |
|------|------|------|
| giac 全库 `define_unary_function_ptr5(at_*)` | **~1852** | 含 GUI/兼容 |
| 测试覆盖 API | **~210** | 见 [`builtin-api-map.md`](builtin-api-map.md) |
| MVP 需实现 | **~210** | 以测试为准，非全库 |

注册机制：`define_unary_function_ptr5(at_foo, …, &__foo, …)` 将字符串 `"foo"` 绑定到 `at_foo` 函数指针；解析器 `T_UNARY_OP` 直接产生 `gen(at_foo, arity)`。

Rust 侧建议：

```rust
pub enum FuncKind { Integrate, Factor, /* … */ }

pub static BUILTINS: phf::Map<&'static str, FuncKind> = /* 210 项 MVP 表 */;
```

完整对照表：[`builtin-api-map.md`](builtin-api-map.md)。
