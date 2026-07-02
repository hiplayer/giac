# GIAC — 任意精度浮点（MPFR 等价 / `Digits` / `evalf` 高精度）

**状态:** open（plan，待立项实现）
**类型:** 功能补齐 / upstream 等价
**根因:** giac C++ 上游有完整 MPFR/MPFI 浮点栈（`real_object` + `real_interval` + `Digits` + 全套超越函数），giac-rs 当前仅有 `f64`，无法跑 `evalf(sin(1/10), 50)` 等 golden
**Rust 落点:** 新建 `giac-core::float`（或 `giac-core::real`）模块 + `giac-core::context::decimal_digits` 字段 + `Expr::Real` variant
**相关:** [GIAC-poly-f5-fglm-complete-square-decision](GIAC-poly-f5-fglm-complete-square-decision.md)（**不依赖**——本 issue 与之并行；LLL/完备平方判定用 `Ratio<BigInt>` 不用任意精度浮点）、[GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md)、[GIAC-expr-api-test-audit](GIAC-expr-api-test-audit.md)
**快照:** 2026-06-29

---

## 0. 背景与目标

### 0.1 upstream 有什么

giac C++ 的 `gen` 类型有一个 `real_object` variant（`giac/giac-2.0.0/src/gen.h:318-390`）：

```cpp
class real_object {
#ifdef HAVE_LIBMPFR
    mpfr_t inf;          // 优先 MPFR
#else
    mpf_t inf;           // 回退 GMP mpf
#endif
    ...
    virtual gen sqrt() const;
    virtual gen exp() const;
    virtual gen log() const;
    virtual gen sin() const;  // ... 全套 sinh/cosh/tanh/asin/acos/atan/...
};
```

外加 `real_interval`（`gen.h:399`，基于 MPFI，区间算术）。

用户层：`Digits` 命令（`prog.cc:8178`）+ `decimal_digits(context)` 全局精度（`global.cc:2524`）。

### 0.2 giac-rs 现状

| 子件 | 状态 |
|---|---|
| 浮点类型 | **仅 `f64`**（`crates/giac-core/src/float_format.rs` 只格式化 f64） |
| 任意精度浮点 | **无** |
| 区间算术 | **无** |
| `Digits` 命令 | **无** |
| 超越函数到 N 位 | **无**（仅 `f64` 标准库 `sin`/`exp`/...） |
| `evalf` 到 N 位 | **无** |
| `f64_to_rational` | 有（`num_util.rs:80`，反向桥；仅 f64 精度） |

`rg -n 'mpfr\|mpf\|bigfloat\|RealObject\|Digits\|decimal_digits' crates/` 全部 zero matches。

### 0.3 目标

把 giac upstream 的「任意精度浮点 + 全套超越函数 + `Digits` + `evalf` 高精度」对齐到 giac-rs，使 `evalf(sin(1/10), 50)` 等 golden 能跑通。

### 0.4 不做的事

- **不**为 [GIAC-poly-f5-fglm-complete-square-decision](GIAC-poly-f5-fglm-complete-square-decision.md) 的 LLL / 完备平方判定引入——该场景用 `Ratio<BigInt>`，本 issue 与之独立
- **不**实现 MPFI 区间算术（多数 golden 不需要，留作 follow-up）
- **不**形式化证明浮点算法正确性（ Lean 后续）
- **不**移植 upstream 的 `mpfr_t` C ABI——纯 Rust 重写

---

## 1. 数学 / 工程骨架

### 1.1 类型设计

```rust
// crates/giac-core/src/real.rs
pub struct RealObject {
    value: FBig,         // dashu 任意精度浮点
    precision: u32,      // bits of precision（对应 upstream decimal_digits 换算）
}

pub enum Expr {
    ...
    /// 任意精度浮点（对应 upstream gen::type=_REAL_）
    Real(Arc<RealObject>),
    ...
}
```

### 1.2 精度模型

| upstream | giac-rs 对应 | 备注 |
|---|---|---|
| `decimal_digits(ctx) = N`（十进制位数） | `ctx.decimal_digits = N` | 用户接口 |
| `mpfr_t` 内部精度 = `bits2digits⁻¹(N)` | `FBig::with_precision(bits)` | `dashu` 用 bits |
| `evalf(x, N)` → `RealObject` precision N | 同 | 主入口 |

### 1.3 `evalf` 分流

```rust
fn evalf(e: &Expr, ctx: &Context) -> Result<Expr, EvalError> {
    let n = ctx.decimal_digits;
    if n <= 15 {
        // 走现有 f64 路径（保持性能 + 与 upstream Digits≤13 默认一致）
        evalf_f64(e).map(Expr::Double)
    } else {
        // 走任意精度路径
        evalf_real(e, n).map(|r| Expr::Real(Arc::new(r)))
    }
}
```

### 1.4 超越函数对齐

| upstream `real_object::` | `dashu::FBig::` | 备注 |
|---|---|---|
| `sqrt` | `SquareRoot::sqrt` | dashu_base trait |
| `exp` | `exp` | 内置 |
| `log` | `ln` | 内置 |
| `sin/cos/tan` | `sin/cos/tan` | dashu 0.4 起（PR #60） |
| `sinh/cosh/tanh` | （需自写或加依赖） | `exp(x) ± exp(-x)` 简单组合 |
| `asin/acos/atan/atan2` | `asin/acos/atan/atan2` | dashu 0.4 起 |
| `asinh/acosh/atanh` | （需自写） | log 形式简单 |
| `abs/inv/is_zero/is_inf/is_nan` | `abs/rec/is_zero/...` | FBig 不支持 NaN（panic 代替） |

**关键差异**：dashu `FBig` **无 NaN**——溢出/除零会 panic 而非返回 NaN。giac-rs 工程规范要求 `Result<_, EvalError>`，所以需要包一层：把 panic 转 `EvalError::FloatDomain`。

---

## 2. 缺口盘点

### 2.1 已有（无需新写）

| 子件 | 位置 |
|---|---|
| `f64` 浮点格式化 | `giac-core::float_format` |
| `f64_to_rational` 反向桥 | `giac-core::num_util:80` |
| `Context.epsilon: f64` | `giac-core::context:34` |
| `Expr` 类型框架 | `giac-core::expr` |
| `evalf` 入口 | 待确认（应该在 `giac-calculus` 或 `giac-simplify`） |

### 2.2 缺件（需新写）

| 子件 | 工作量估 | 依赖 |
|---|---|---|
| **A1** 引入 `dashu` workspace 依赖 + WASM 编译验证 | 小（~5 行 Cargo + 1 验证脚本） | `dashu = "0.4"` |
| **A2** `RealObject` 结构 + `Expr::Real` variant | 小（~50 行） | A1 |
| **A3** `Context::decimal_digits` 字段 + `Digits` 命令解析 | 小（~30 行） | A2 |
| **A4** `evalf` 分流（f64 vs Real）+ 基础算术 (`+ - * /`) | 中（~200 行） | A2 + A3 |
| **A5** 超越函数 `sqrt/exp/ln/sin/cos/tan/asin/acos/atan/atan2` 接 dashu | 中（~300 行） | A4 |
| **A6** 双曲 / 反双曲（`dashu` 缺，自写） | 小中（~150 行，exp(x)±exp(-x) / log 形式） | A5 |
| **A7** `real2int` / `real2double` / `real2rational` 桥（对应 upstream `gen.h:397-398`） | 小（~80 行） | A4 |
| **A8** panic → `EvalError::FloatDomain` 包装 | 小（~30 行 + 测试） | A4-A6 |
| **A9** `print` / `dbgprint` 格式化（对齐 upstream `printmpf_t`） | 中（~150 行） | A4 |
| **A10** `evalf` golden 测试接入（conformance） | 中（~300 行测试 + 跑 N 位 golden） | A5 + A9 |

### 2.3 依赖决策

| 选项 | 评价 |
|---|---|
| **`dashu = "0.4"`**（推荐） | 纯 Rust、no_std、WASM-safe、MIT/Apache 双许可、MSRV 1.68（giac-rs MSRV 1.75 兼容）、有 sin/cos/tan/exp/ln/sqrt/cbrt/nth_root/asin/acos/atan/atan2/sin_cos |
| `rug` (GMP FFI) | **禁**（WASM 场景禁止，giac-rust-engineering.mdc 明确） |
| `faer` | **禁**（WASM 不兼容） |
| `rust_decimal` | 仅 128-bit decimal，精度上限不够 |
| `bigdecimal` | 同上 |
| `malachite` | 任意精度整数强，浮点支持弱 |
| 自写 | 工作量 2000+ 行，重复造轮子 |

**结论**：选 `dashu`。它就是「Rust native MPFR 替代品」的明确目标（README 原话："intended to be a Rust native alternative to GNU GMP + MPFR"）。

---

## 3. 实现优先级（三阶段）

### Phase A — 基础设施 + `Digits` 命令（先做，~1 周）

| ID | 内容 | 工作量 |
|---|---|---|
| **A1** | `dashu` workspace 依赖 + `cargo build --target wasm32-unknown-unknown` 验证 | 0.5 天 |
| **A2** | `RealObject` 结构 + `Expr::Real` variant + 基础 `Display` | 1 天 |
| **A3** | `Context::decimal_digits` 字段 + `Digits` 命令解析（读写） | 0.5 天 |
| **A4** | `evalf` 分流框架 + `+ - * /` 算术 + `f64 ↔ Real` 互转 | 2 天 |
| **A7** | `real2int/double/rational` 桥 | 1 天 |
| **A8** | panic → `EvalError` 包装 + 测试 | 1 天 |

**DoD**：
- `Digits := 50; evalf(1/3)` 返回 `0.33333...`（50 位）
- `evalf(2.5 + 1.25)` 返回 `3.75`（精度保持）
- WASM 编译绿（`cargo build -p giac-wasm --target wasm32-unknown-unknown`）
- `cargo test-timeout` 全绿（不破坏现有 f64 测试）

### Phase B — 超越函数对齐（核心 golden 解锁，~2 周）

| ID | 内容 | 工作量 |
|---|---|---|
| **A5** | `sqrt/exp/ln/sin/cos/tan/asin/acos/atan/atan2` 接 dashu | 5 天 |
| **A6** | 双曲 / 反双曲自写 | 2 天 |
| **A9** | `print` 格式化对齐 upstream `printmpf_t`（N 位科学计数法） | 3 天 |
| **A10** | `evalf` golden 测试接入（跑 `giac_check_cas` 里的 N 位用例） | 2 天 |

**DoD**：
- `evalf(sin(1/10), 50)` 与 upstream giac 输出 diff ≤ 1 ulp（或字面一致，看 golden 类型）
- `evalf(exp(1), 100)` 同
- `evalf(pi, 50)` 同（`pi` 用 dashu 内置 `FBig::PI` 或自写 AGM）
- conformance golden 里的高精度 `evalf` 用例 100% 通过（或登记 known-divergence）

### Phase C — 完备化与 follow-up（独立排期）

| ID | 内容 | 价值 |
|---|---|---|
| **C1** MPFI 区间算术（`real_interval`） | rigor computation / 数值保证；多数 golden 不需要 |
| **C2** `evalf` 复数路径（`Complex<Real>`） | `exp(i*pi)` 等；可推迟 |
| **C3** 性能优化（dashu 内部精度切换 + 缓存 `pi`/`e` 常数） | 大规模 `Digits=1000` 用例性能 |
| **C4** 浮点 ↔ 有理的高精度桥（`real_to_rational` 给 LLL 用） | 跨域桥；目前 LLL 走 `Ratio<BigInt>`，不需要 |

---

## 4. 落地建议

### 4.1 推荐执行顺序

```
Phase A（基础设施，~1 周）
  ├─ A1 dashu 引入 + WASM 验证
  ├─ A2 RealObject + Expr::Real
  ├─ A3 decimal_digits + Digits 命令
  ├─ A4 evalf 分流 + 算术
  ├─ A7 real2int/double/rational 桥
  └─ A8 panic→EvalError 包装

Phase B（超越函数 + golden，~2 周）
  ├─ A5 sqrt/exp/ln/sin/cos/tan/...
  ├─ A6 双曲 / 反双曲自写
  ├─ A9 print 对齐 upstream
  └─ A10 evalf golden 接入

Phase C（完备化，独立）
  ├─ C1 MPFI 区间算术
  ├─ C2 复数 evalf
  ├─ C3 性能优化
  └─ C4 跨域桥
```

### 4.2 与现有 issue 的关系

- **独立于** [GIAC-poly-f5-fglm-complete-square-decision](GIAC-poly-f5-fglm-complete-square-decision.md)：完备平方判定走 `Ratio<BigInt>`，不用任意精度浮点。两者可完全并行。
- **依赖** [GIAC-expr-api-tech-debt](GIAC-expr-api-tech-debt.md) 的 `Expr` variant 扩展机制（待确认）
- **解锁** [conformance 测试](../conformance-testing.md) 中所有 `evalf(..., N>15)` golden 用例

### 4.3 风险登记

| 风险 | 缓解 |
|---|---|
| `dashu` panic 替代 NaN 与 giac-rs `Result` 规范冲突 | A8 显式包装 + `EvalError::FloatDomain` |
| `dashu` 超越函数精度与 upstream MPFR 不完全一致 | A10 用 upstream golden 做对照，diff > 1 ulp 登记 known-divergence |
| `dashu` 性能比 MPFR 慢 2-5× | Phase C C3 优化；接受 Phase B 性能差（golden 不卡时间） |
| WASM 编译失败 | A1 先验证；`dashu` README 明确 no_std 支持，风险低 |
| `Expr::Real` 与 `Expr::Double` 共存导致下游 `match` 漏 arm | 用 exhaustive match + clippy `match_same_arms` 强制 |
| `decimal_digits` 默认值与 upstream 不一致（upstream 默认 13） | A3 显式设默认 13 + 测试覆盖 |

### 4.4 何时考虑性能优化

只有当：
1. Phase B 完成后某 `Digits=1000` golden 用例 > 5s，且
2. profile 显示 >50% 时间在 `dashu` 内部

才进入 Phase C C3。否则 `dashu` 默认性能对 golden 套件足够。

---

## 5. 参考

### 5.1 upstream 实现

| 子件 | 位置 |
|---|---|
| `real_object` 类定义 | `giac/giac-2.0.0/src/gen.h:318-390` |
| `real_interval` 区间算术 | `giac/giac-2.0.0/src/gen.h:399-` |
| `Digits` 命令注册 | `giac/giac-2.0.0/src/prog.cc:8178` |
| `decimal_digits(ctx)` 全局精度 | `giac/giac-2.0.0/src/global.cc:2524` |
| `bits2digits` 换算 | `giac/giac-2.0.0/src/prog.cc:7950` |
| `set_decimal_digits` | `giac/giac-2.0.0/src/prog.cc:7958` |
| `real2int` / `real2double` 桥 | `giac/giac-2.0.0/src/gen.h:397-398` |
| `evalf` 主入口 | `giac/giac-2.0.0/src/usual.cc`（待确认行号） |

### 5.2 Rust 依赖

| crate | 版本 | 评价 |
|---|---|---|
| **`dashu`** | 0.4.4（2026-05-31） | 首选；纯 Rust、no_std、WASM-safe、MIT/Apache、MSRV 1.68、有全套超越函数（PR #60 加 sin/cos/tan/asin/acos/atan/atan2/sin_cos） |
| `dashu-float` API 参考 | https://docs.rs/dashu/latest/dashu/float/struct.FBig.html | FBig 是主类型，`with_precision / exp / ln / sin / cos / ...` |
| `dashu` GitHub | https://github.com/cmpute/dashu | README 明确「Rust native alternative to GMP + MPFR」 |

### 5.3 性能 benchmark（待 Phase B 完成后实测）

未实测。Phase A1 落地后跑：
- `evalf(pi, 100)` 用时 vs upstream `Digits=100; evalf(pi)`
- `evalf(sin(1/10), 50)` 同
预期 `dashu` 比 MPFR 慢 2-5×（通用经验），对 golden 套件足够。

### 5.4 giac-rs 内部参考

| 子件 | 位置 |
|---|---|
| 当前 `f64` 浮点格式化 | `giac-rs/crates/giac-core/src/float_format.rs` |
| `f64_to_rational` 反向桥 | `giac-rs/crates/giac-core/src/num_util.rs:80` |
| `Context.epsilon: f64` | `giac-rs/crates/giac-core/src/context.rs:34` |
| `Expr` 类型框架 | `giac-rs/crates/giac-core/src/expr.rs`（待确认） |
| workspace `Cargo.toml` | `giac-rs/Cargo.toml`（依赖加在这里） |
| 工程规范（WASM / 依赖） | `.cursor/rules/giac-rust-engineering.mdc` |

### 5.5 工程规范约束

引用 `giac-rust-engineering.mdc`：

> - 业务 crate **默认 `#![deny(unsafe_code)]`**；仅 `giac-ffi` 允许 FFI
> - **WASM**：`wasm32-unknown-unknown` 须可编译；**禁止** C FFI 传递依赖
> - WASM 场景：**禁止** `rug`（GMP FFI）；用 `num-bigint`
> - 数值 linalg：**`nalgebra`**；**禁止** `faer`（WASM 不兼容）

`dashu` 满足全部约束（纯 Rust、no_std、无 FFI、WASM-safe）。本 issue 引入它符合规范。
