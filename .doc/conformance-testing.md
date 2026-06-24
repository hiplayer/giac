# Conformance 测试规格

Golden 回归、等价判定、随机数与归一化规则。配合 [`rust-migration-plan.md` §6](rust-migration-plan.md) 使用。

**算法 crate 内单元/通路单测分层（A/B/C）与审计表：** [test-writing-spec.md](test-writing-spec.md)、[issues/GIAC-expr-api-test-audit.md](issues/GIAC-expr-api-test-audit.md)。

---

## 1. 测试计数（统一口径）

| 类别 | 数量 | CTest 前缀 |
|------|------|------------|
| libtommath 冒烟 | 1 | `tommath_smoke_test` |
| bin 覆盖脚本 | **41** | `giac_bin_*` |
| check golden diff | **10** | `giac_check_*` |
| 已知 segfault | 1 | `giac_check_factor_segfault_known` |
| **合计** | **53** | |

**输入文件：** 41（bin）+ 10（check）= **51** 个；表达式约 **1210** 条；API **~210** 种（见 `functional-coverage.md`）。

> 旧文档「42 bin」为过时计数；以本表为准。

---

## 2. Golden 归一化（移植 `cmake/giac_check_normalize.sh`）

### 2.1 `cas_floats`（`giac_check_cas`）

对所有裸浮点字面量应用 **10 位有效数字** `%.10g`：

```perl
s/(?<![\w.])(-?\d+\.\d+)(?![\w.])/sprintf("%.10g", $1)/ge
```

Rust 实现：

```rust
fn normalize_cas_floats(s: &str) -> String {
    // 正则: (?<![\w.])(-?\d+\.\d+)(?![\w.])
    // 替换为 format!("{:.10g}", val)
}
```

### 2.2 `geo`（`giac_check_geo`）

对含 `pnt(pnt` 的行替换不稳定 ID：

| 模式 | 替换 |
|------|------|
| `[536870xxx]` | `[ID]` |
| `,536870xxx,` | `,ID,` |
| `,[0-9]{1,10}]` | `,L]` |
| `,[0-9]{1,10},` | `,L,` |

### 2.3 其他 check 测试

| CTest | 归一化 | 备注 |
|-------|--------|------|
| `giac_check_integrate` | `LANG=en_US.UTF-8` | 区域设置 |
| `giac_check_cas` | `cas_floats` + unset LANG | 删 `session.tex` |
| `giac_check_geo` | `geo` + UTF-8 | |
| 其余 | 无 | 字面 diff |

---

## 3. `assert_equiv` 规格

比字面 golden 更权威的**数学等价**判定；用于积分/求解多解、化简排序差异。

### 3.1 核心算法

```rust
/// 判定 a 与 b 在 ctx 下是否数学等价
pub fn assert_equiv(a: &Expr, b: &Expr, ctx: &Context) -> Result<bool, EvalError> {
    let d = simplify::normal(&sub(a, b)?, ctx)?;
    Ok(is_zero(&d, ctx))
}
```

### 3.2 `sub` / `is_zero`

| 函数 | 语义 | 边界 |
|------|------|------|
| `sub(a,b)` | 构造 `a + (-b)` 或代数减 | 类型不匹配 → `TypeError` |
| `is_zero(e, ctx)` | `normal(e)==0` 或 `e` 在假设下恒为 0 | 见 §3.3 |
| `normal` | 与 giac `normal` 同规格目标，算法可不同 | 必须终止 |

### 3.3 等价判定边界

| 情形 | 处理 |
|------|------|
| 浮点结果 | 先 `evalf` 到 `ctx.float_digits`，比较差 < `ctx.epsilon` |
| 复数 | 比较实部/虚部分别 `is_zero` |
| 矩阵/向量 | 逐元素 `assert_equiv` |
| 方程解（集合） | 转为集合差 `A \\ B` 与 `B \\ A` 均为空（后置） |
| 积分常数 `+C` | 剥离常数项再比，或 `diff` 还原被积函数 |
| 超时 | `max_recursion` /  wall-clock 上限 → 回退字面 diff 并记 DIV-* |

### 3.4 多解 API 策略

| API | 策略 |
|-----|------|
| `integrate` | `diff(rust) == integrand` 且 `diff(giac) == integrand`；或 `assert_equiv(rust, giac)` |
| `solve` | 解集代入原方程；或有限解逐一 `assert_equiv` |
| `factor` | 展开乘积等于原式 |
| `normal` / 化简 | 默认 `assert_equiv` |

### 3.5 命令 I/O 契约（L1 normative）

L1 的**数学属性**见 §3.4；**输入/输出形态**写在各 crate 的 `*-api-stability.md`（[algorithm-expr-api.md §4](algorithm-expr-api.md#4-各-crate-规范入口索引) 索引）。

| 命令族 | 契约章节 | Conformance 门禁 | L1 属性（不可改脚本弱化） |
|--------|----------|------------------|---------------------------|
| `factor` | [giac-simplify-api-stability.md §5.1](giac-simplify-api-stability.md#51-factorexpr--io-契约normative) | `giac_check_factor` | `expand(factor(p)) ≡ expand(p)`；因子 ∈ ℚ[…] |
| `integrate` / `int` | `giac-calculus-api-stability.md` | `giac_check_integrate` | 微分还原 / `assert_equiv` |
| `limit` | `limit-engine-expr-api.md` | `giac_check_limit` | 极限值属性 |
| `solve` / `roots` | `giac-solve-api-stability.md` | `test_solve` 等 | 解代入原方程 |
| `normal` / `expand` | `giac-simplify-api-stability.md` | `giac_check_cas` 等 | `assert_equiv` |

**动到某命令的实现时，PR 描述须引用上表对应契约章节**（或新建的 `*-expr-api.md`）。

**禁止为让 L1 绿而修改：**

- `tests/conformance/scripts/sympy_verify.py` 中该命令的 `verify_property` / 属性判定
- 把 L1 改成弱断言（子串 `contains`、非空输出、删 failing 用例）

改 L1 定义本身：先改契约文档 → 人工确认与 giac 语义一致 → 再改脚本/测试。

### 3.6 L1 失败处理（须人工确认）

Conformance L1（SymPy 属性 / `giac_check_*` 逐行测）失败时，**不得**为求 CI 绿而直接改功能实现或弱化门禁。流程：

```text
1. 查因 — 对照 §3.5 契约章节 + 上游 golden（check/*.out）+ 必要时 giac C++（module-division.md）
2. 能修且符合契约 — 改实现；PR 引用契约章节；跑对应 giac_check_* 全量
3. 短期不能修 — 双轨（见 [test-writing-spec.md §8](test-writing-spec.md#8-l1-conformance-失败须人工确认)），须登记
4. 改 L1 规格 — 仅人工确认后：先改契约文档，再改 sympy_verify / 测试
```

**双轨（L1 专用）：**

| 轨道 | 做法 | CI |
|------|------|-----|
| **目标 L1** | 原 `factor_sympy_line_XX` 等加 `#[ignore = "L1-*: 原因"]` | 默认跳过 |
| **smoke-until** | 保留或新增弱测（非空、结构、crate 内 A 层）；注释 `smoke-until L1-*: delete when …` | 默认须绿 |
| **登记** | [GIAC-expr-api-test-contains-cleanup.md §8 L1 表](issues/GIAC-expr-api-test-contains-cleanup.md#81-l1-conformance-阻塞l1-) | 必做 |

示例（conformance 单条 L1）：

```rust
// smoke-until L1-FACTOR-15: delete when `factor_sympy_line_15` green
#[test]
fn factor_smoke_line_15_non_empty() {
    let out = giac_conformance::run_line("factor((x^2-3*x+1)*(x^2+x+1))").unwrap();
    assert!(!out.is_empty());
}

#[test]
#[ignore = "L1-FACTOR-15: stack overflow in factor_poly_form / FieldSession Δ<0"]
fn factor_sympy_line_15() -> Result<(), String> {
    giac_conformance::assert_factor_line_sympy(15)
}
```

解除：`#[ignore]` 去掉 → L1 绿 → **删除** smoke-until → 删 §8 L1 行。

本地验 L1：`cargo nextest run --release -p giac-conformance --test giac_check_factor factor_sympy_line_15`；验 ignore 目标：`cargo test -p giac-conformance --test giac_check_factor -- --ignored factor_sympy_line_15`。

### 3.7 Conformance 集成

```rust
enum CheckOutcome {
    LiteralMatch,
    EquivMatch,      // assert_equiv 通过
    KnownDivergence, // known-divergences.md 登记
    Fail,
}

fn check_line(input: &str, expected: &str, ctx: &Context) -> CheckOutcome {
    let got = eval_and_print(input, ctx)?;
    if normalize(&got) == normalize(expected) {
        return LiteralMatch;
    }
    if assert_equiv(&parse(&got)?, &parse(expected)?, ctx)? {
        return EquivMatch;
    }
    if is_known_divergence(input, &got, expected) {
        return KnownDivergence;
    }
    Fail
}
```

---

## 4. 随机数 seed 策略

### 4.1 giac 实现链

| 组件 | 位置 | 说明 |
|------|------|------|
| PRNG | `global.cc::giac_rand` | 桌面：**tinymt32**（`tinymt32.cc`） |
| seed | `prog.cc::srand` / `rand_seed` | `srand(t)` → `rand_seed(t, ctx)` |
| 正态 | `moyal.cc::randNorm` | Box-Muller，消耗 2× `giac_rand` |
| χ² | `moyal.cc::randchisquare` | 2k 个 uniform |

`rand_max2` = `0x7FFFFFFF`（`prog.h`）。

### 4.2 测试策略

| 测试 | 策略 |
|------|------|
| `bin/test_special`, `test_probstat` | 脚本 **不** 断言具体随机序列；只验 **类型**、有限性、分布矩（可选） |
| 需可重复 diff | harness 在脚本前注入 `srand(固定值);` |
| Rust 默认 | `StdRng::seed_from_u64(0xGIAC_TEST_SEED)`，**不**复现 tinymt32 位级序列 |
| 位级兼容（可选） | 移植 `tinymt32` 或 crate 包装；仅当 golden 含具体随机输出时需要 |

**推荐 MVP：**

```rust
const TEST_RNG_SEED: u64 = 42424242;

pub fn test_context_with_fixed_rng() -> Context {
    let mut ctx = Context::default_check_setup(); // cas_setup 等价
    ctx.rng = TestRng::seed(TEST_RNG_SEED);
    ctx
}
```

对 `randNorm` / `randchisquare`：

1. **不** 与 giac 浮点输出逐字 diff；
2. 断言 `is_finite`、实数类型；
3. 大样本检验均值≈μ、方差≈σ²（属性测试，容差 `proba_epsilon`）。

### 4.3 若未来需要 giac 序列对齐

1. 在 conformance 开头统一 `srand(12345);`
2. Rust 侧实现 `tinymt32_generate_uint32` 兼容层；
3. 记录于 `known-divergences.md`（浮点 Box-Muller 顺序须一致）。

---

## 5. Harness 初始化清单

每个 check/bin 脚本运行前：

```text
1. Context ← cas_setup 测试默认值（见 parser-token-map.md §5）
2. xcas_mode ← 0
3. （可选）srand 固定 seed
4. LANG ← 测试要求（integrate/geo 用 en_US.UTF-8）
5. 逐行 parse → eval → print
6. normalize → diff golden
7. 失败时尝试 assert_equiv → 仍失败查 known-divergences.md

### 5.1 单测超时（推荐 `cargo test-timeout`）

单元测试与 conformance 均可能因逻辑死循环挂起。giac-rs 全 workspace **800+** 测，**默认用 `cargo test-timeout`**，不要习惯性跑裸 `cargo test --workspace`。

**两层超时（职责不同）：**

| 层 | 机制 | 作用 |
|----|------|------|
| **整测** | `giac-rs/.config/nextest.toml` → `slow-timeout`（默认 **10s**） | 挂死的 Rust `#[test]` 进程被 nextest 终止，套件继续 |
| **子进程** | `GIAC_CHECK_TIMEOUT_SECS`（默认 **10**）→ `subprocess_timeout()` | SymPy `python3` 校验：`command_with_timeout` 轮询 + `kill`，避免 `wait` 无限阻塞 |

Rust `eval` / `run_line` / `run_lines` **不再**包线程级 per-line 超时；挂死由 nextest 整测上限兜底。多行脚本合在一个 `#[test]` 里时，任一行挂死会耗尽该测的 nextest 预算。

**为什么更快：**

| | `cargo test-timeout` | `cargo test --workspace` |
|--|----------------------|---------------------------|
| 执行器 | cargo-nextest，多核并行 + 单测超时杀进程 | libtest，无单测超时 |
| 挂死 | 超时后跳过，套件继续 | 整个进程卡死 |
| 无 nextest 时 | `./scripts/test-with-timeout.sh` 逐测 GNU `timeout`（极慢） | — |

**安装与运行：**

```bash
cd giac-rs
cargo install cargo-nextest --locked --version 0.9.85   # rustc 1.75；0.9.86+ 需 1.91
cargo test-timeout                      # 推荐：workspace 全量
# 或
./scripts/test-with-timeout.sh
```

超时策略见 `giac-rs/.config/nextest.toml`（`cargo test-timeout` 经 `.cargo/config.toml` 使用 **--release**；默认 slow-timeout **10s**）。SymPy 子进程上限：`GIAC_CHECK_TIMEOUT_SECS`（默认 10），见 `giac_conformance::subprocess_timeout()`。

**调试子集**可用 debug `cargo nextest run`（较慢，重算路径可能触发 slow-timeout）：

```bash
cargo test -p giac-calculus ck_int_61
cargo test -p giac-conformance --test giac_check_integrate giac_check_integrate_ck_int_60
```

超时后脚本会打印 `RUST_BACKTRACE=1 cargo test -p …` 便于定位；典型根因：`poly_divrem` 索引错误导致 `inv_EXT` 无限循环。

---

## 6. 外部开源语料

Maxima / SymPy / Rubi 等第三方测试的目录布局、JSON schema、`extract_maxima_rtest.py` 规范见 **[external-test-resources.md](external-test-resources.md)**。

giac-rs 路径：`giac-rs/tests/conformance/fixtures/{maxima,sympy,rubi}/`。

---

## 7. 代数 `ext_tower` 与可审计的正确性验证

Golden 覆盖端到端表达式；**域扩张 / compositum** 另有一套分层证据，不替代 check golden，与之互补。

| 文档 | 内容 |
|------|------|
| [giac-tower-common-math.md](giac-tower-common-math.md) | compositum 数学参考、Lean 分层引理、**§4 可审计验证四层做法** |
| [GIAC-lazy-common-tower-plan.md](issues/GIAC-lazy-common-tower-plan.md) | T4a/T4b 实施与坐标基 §12.9 |

**日常命令：**

```bash
cargo test -p giac-core                                    # 默认 T4a 塔 common
cargo test -p giac-core --no-default-features            # bisect：Phase 0 flatten
cargo test -p giac-core tower_common_matches_flatten -- --ignored  # flatten 慢对照
```

**原则：** 数学正确性以 oracle + 不变量为主；Lean 证引理、Rust 测实例（见 giac-tower-common-math §4.2 层 D）。

**并发：** R4/R5 后无进程级 `FieldRegistry`；`ext_tower` 单元测可并行（`cargo nextest` 默认多核）。
