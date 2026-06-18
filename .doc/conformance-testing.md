# Conformance 测试规格

Golden 回归、等价判定、随机数与归一化规则。配合 [`rust-migration-plan.md` §6](rust-migration-plan.md) 使用。

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

### 3.5 Conformance 集成

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
```

---

## 6. 外部开源语料

Maxima / SymPy / Rubi 等第三方测试的目录布局、JSON schema、`extract_maxima_rtest.py` 规范见 **[external-test-resources.md](external-test-resources.md)**。

giac-rs 路径：`giac-rs/tests/conformance/fixtures/{maxima,sympy,rubi}/`。
