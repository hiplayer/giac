# GIAC limit @ +∞：对齐 upstream `unidirectional_limit` 与 Maxima 回归

## 背景

Maxima `tests/rtest_limit*.mac` 中抽取了 14 条 `limit` 用例（`giac-calculus::limit::tests::maxima_rtest`）。此前 10 条标记为 `#[ignore]`，根因是 giac-rs 在 `+infinity` 上依赖多项式快路径 / 倒数 L'Hôpital，未走 upstream `series.cc` 的 **MRV 主路径**。

本次按 upstream 完整解法重构 `+∞` 极限管线，并记录仍待实现的缺口。

## Upstream 对照（`giac-1.5.0/src/series.cc`）

| Upstream | giac-rs 对应 |
|----------|----------------|
| `_pow2exp` / `pow2expln` 预处理 | `risch::pow2expln` + `preprocess::merge_exp_quotients` + `fold_exp_zero_linear` |
| `unidirectional_limit` → `mrv_lead_term` | `limit_unidirectional_plus_infinity` → `mrv_lead_term_plus_infinity` → `limit_from_mrv_lead_term`（coeff 递归） |
| `mrv_compare`（指数快慢比较） | **部分**：`choose_mrv_w` + `linear_coeff_in_var` + `is_negative_const_expr`（常数 `ln` 用 `f64` 估值） |
| 有限点 `x + 1/x` 换元 | `limit_at_plus_infinity_fallback`：`x=1/u` + `series_at_zero` / Laurent |
| `series__SPOL1` ordre 循环 | `mrv_series_lead_loop`（GIAC-216e，已有） |

## `+∞` 调用顺序（`limit_at_plus_infinity`）

> **计划变更**：步骤 2–4 为写死快路径，拟按本文 [「去写死」方案](#去写死sqrt代数式-极限对齐方案) 移除，统一为 MRV + 倒数级数。

1. **`limit_unidirectional_plus_infinity`**（主路径，对齐 upstream）
2. ~~有理式首项（`limit_rational_leading_at_infinity`）~~ → 并入 fallback Laurent
3. ~~代数共轭（`sqrt` 差分 / `sqrt(S)-var`）~~ → **拟删除**
4. ~~有理式 / `sqrt` 商（CK-INT-58 形状）~~ → **拟删除**，改 `surd2pow` + 级数
5. **倒数换元 + 稀疏级数**（`limit_at_plus_infinity_fallback`）

不再把「仅含 `exp` 的嵌套式」作为 MRV 门槛；`mrv_limit_eligible` 仅做节点/深度界。

## 主要代码改动

| 模块 | 改动摘要 |
|------|----------|
| `preprocess.rs` | `merge_exp_quotients`、`fold_exp_zero_linear` |
| `mrv.rs` | `Pow(const,var)` 进 MRV；`choose_mrv_w` 支持一般负线性指数；`linear_coeff_in_var` |
| `mrv_lead_term.rs` | `mrv_limit_eligible`；`limit_unidirectional_plus_infinity`；`omega_tends_to_zero` 扩展 |
| `asymptotic.rs` | MRV 优先；共轭 `sqrt(S)-x`；`limit_at_zero_rational_lead` |
| `sparse_series.rs` | `atan(1/u)`、`sqrt` 二项级数 |
| `mod.rs` | `limit_via_reciprocal` 改走 `limit_at_zero_fallback`（避免 atan L'Hôpital 挂起） |

## Maxima 回归状态（2026-06-17，续）

| 用例 | 期望 | 状态 |
|------|------|------|
| `sin(x)/x` @ 0 | 1 | ✅ |
| `(1-cos)/x²` @ 0 | 1/2 | ✅ |
| `(1+1/n)^n` @ +∞ | `exp(1)` | ✅ |
| `a/n` @ +∞ | 0 | ✅ |
| `7^n/8^n` @ +∞ | 0 | ✅ MRV |
| `4^n/2^(2n)` @ +∞ | 1 | ✅ |
| `x*(sqrt(1+x²)-x)` @ +∞ | 1/2 | ✅ 共轭 / fallback |
| `x/(x^ln(x))` @ +∞ | 0 | ✅ 形状识别 |
| `(1+1/x)*(sqrt(x+1)+1)` @ +∞ | +∞ | ✅ `limit_poly_over_sqrt` |
| `x*atan(x)/(x+1)` @ +∞ | `pi/2` | ✅ 倒数级数 + `series_atan_of_inv` |
| `(3^x+5^x)^(1/x)` | 5 | ✅ `limit_exp_sum_nth_root` |
| gruntz `exp*(exp(...)-exp(...))` | -1 | ❌ MRV 级数 |
| CK-INT-60 比值 | 1 | ❌ MRV 比值收敛 |
| gruntz 嵌套 exp 差 | 1 | ❌ MRV 级数 |

**通过：11/14**；**仍 ignore：3 条 gruntz**。

### Phase 0 fallback-only（`limit_at_plus_infinity_fallback`）

| 用例 | 期望 | 状态 |
|------|------|------|
| `(x+1)/(x-1)` | 1 | ✅ |
| `x*(sqrt(1+x²)-x)` | 1/2 | ✅ |
| CK-INT-58 | +∞ | ✅ |
| `x*atan(x)/(x+1)` | `pi/2` | ✅ |
| `(1+1/x)*(sqrt(x+1)+1)` | +∞ | ✅ |

**fallback 矩阵：5/5 全绿。**

## 后续（对齐 upstream 完整实现）

1. **`mrv_compare`**：`ln(a)/ln(b)` 的 MRV 主项比较（`3^x` vs `5^x`、`x^ln(x)` 等）
2. **`series_div`**：支持含 `pi` 的常数主项（`atan` @ +∞）
3. **CK-60 比值**：扩展 `remove_lnexp` / peel（无 `exp(x)` 差分形状）
4. **gruntz 嵌套 exp**：216e 级数路径覆盖更多 `exp` 差分

## 测试

```bash
# 全部 maxima 子集（含 ignore）
cargo test -p giac-calculus --lib maxima_rtest -- --include-ignored

# 仅 CI 用例（无 ignore）
cargo test -p giac-calculus --lib maxima_rtest
```

---

## 去写死：`sqrt`/代数式 `+∞` 极限对齐方案

> 目标：弱化或移除 `rationalize_sqrt_in_expr` 等形状快路径，把 `sqrt`/代数式的 `+∞` 极限统一交给 **MRV** 或 **倒数换元 + 稀疏级数**，与 upstream 一致。

### 1. 当前与 upstream 的差异

#### 1.1 upstream `limit` @ `+∞` 实际管线（`series.cc`）

```
limit(e, x, +inf)
  ├─ 直接代入 / partfrac @ lim_point（有理式常可在此结束）
  ├─ surd2pow（含 sqrt / surd → 有理指数幂，再递归 limit）
  ├─ limit_symbolic_preprocess + pow2expln
  ├─ 双向极限且解析：series__SPOL1(e, x, +inf, ordre)（ordre 递增）
  └─ unidirectional_limit：
        +∞ 时 **不做** x=1/u，直接 mrv_lead_term → 由 exponent/coeff 判 0/有限/±∞
```

要点：

- **没有** `rationalize_sqrt_in_expr`、`limit_conjugate_sqrt_at_infinity` 这类独立分支。
- `sqrt` 通过 **surd2pow** 纳入统一幂级数 / MRV 体系，而非共轭有理化 + 多项式次数比较。
- 倒数换元 `lim_point ± inv(x)` 用于**有限点**单侧极限；`+∞` 的 unidirectional 路径**不依赖**共轭形状识别。
- 级数引擎 `series__SPOL1` / `in_series__SPOL1` 在 `lim_point = ±∞` 下对 `lvx` 做渐近展开，覆盖 `sqrt`、`sin`、`ln` 等，不仅限于 `exp` MRV 集合。

#### 1.2 giac-rs 当前 `limit_at_plus_infinity` 管线

```
limit_at_plus_infinity
  1. limit_unidirectional_plus_infinity     ← 对齐 upstream（MRV）
  2. limit_rational_leading_at_infinity     ← 写死：多项式次数比较
  3. limit_conjugate_sqrt_at_infinity      ← 写死：sqrt 差分 / sqrt(S)-x 共轭
  4. limit_rational_over_sqrt_quotient_*    ← 写死：CK-INT-58 分式+sqrt 商
  5. limit_at_plus_infinity_fallback      ← 较通用：x=1/u + Laurent / sparse_series
```

| 能力 | upstream | giac-rs 现状 | 差距 |
|------|----------|--------------|------|
| MRV @ `+∞` | `unidirectional_limit` 主路径 | 已有 `limit_unidirectional_plus_infinity` | gruntz / `mrv_compare` 未全 |
| `sqrt` @ `+∞` | `surd2pow` + `series__SPOL1` 或递归 `limit` | **共轭快路径** + 部分 `series_sqrt` | 快路径与 upstream 架构不一致 |
| 有理式 @ `+∞` | `partfrac` / 级数主项 | `limit_rational_leading`（快路径） | 功能重叠，可并入 fallback |
| `atan` @ `+∞` | 级数 / `mrv_lead_term` | 倒数级数未闭环（`series_div`+`pi`） | 无共轭捷径，但级数未够 |
| `surd2pow` 预处理 | `limit` 入口统一 | **未实现** | 根式处理分散在共轭里 |

#### 1.3 「写死」函数清单（拟移除或降级）

| 函数 | 行数级 | 识别形状 | upstream 等价物 |
|------|--------|----------|-----------------|
| `rationalize_sqrt_in_expr` | ~30 | `x*(sqrt±…)`、`sqrt(a)±sqrt(b)` | 无（应走级数） |
| `rationalize_sqrt_minus_var` | ~25 | `sqrt(S)-x` | 无 |
| `rationalize_sqrt_difference` | ~30 | `sqrt(a)-sqrt(b)` | 无 |
| `limit_conjugate_sqrt_at_infinity` | ~40 | 上式 + 多项式首项 | `series__SPOL1` @ `+∞` |
| `limit_rational_over_sqrt_quotient_at_infinity` | ~35 | `(poly)/sqrt(frac)` | `surd2pow` + `limit` / 级数 |
| `limit_rational_leading_at_infinity` | ~25 | 有理式次数比 | `partfrac` @ `+∞` 或倒数 Laurent |

共轭快路径能算对 `x*(sqrt(1+x²)-x)→1/2`，是因为该形状**恰好**落在手写规则里；换形（如 `sqrt(1+x²)-x/2`、`2*sqrt(x²+1)-2x`）即失效，而 upstream 级数路径更稳。

#### 1.4 通用路径已具备但未单独兜住的能力

倒数换元后（`x=1/u`，求 `u→0`）：

| 原式 @ `+∞` | 换元后 @ `0` | `sparse_series` 需求 |
|-------------|--------------|------------------------|
| `x*(sqrt(1+x²)-x)` | `(sqrt(1+u²)-1)/u²` | `series_sqrt` + `series_div`（**已有雏形**） |
| `(x+1)/(x-1)` | `(1+u)/(1-u)` | 有理 Laurent（**已有**） |
| `(x+1)/sqrt((x+1)/(x-1))` CK-58 | 有理式组合 | `sqrt` 级数 + 除法（需加强） |
| `x*atan(x)/(x+1)` | `atan(1/u)/(1+u)` | `series_atan_of_inv` + `series_div`（**部分**） |

结论：**去掉共轭快路径的前提**是补强 `limit_at_zero_fallback` 的级数链，而非再加新形状表。

---

### 2. 目标架构（对齐 upstream）

```
limit @ +∞
  │
  ├─ [A] preprocess
  │     pow2expln, merge_exp_quotients, fold_exp_zero_linear
  │     surd2pow（新增，对齐 limit 入口）
  │
  ├─ [B] limit_unidirectional_plus_infinity   ← 主路径（含 exp / 超越）
  │     失败且非「纯代数」→ 继续
  │
  └─ [C] limit_at_plus_infinity_fallback      ← 代数 / 超越的统一兜底
        x = 1/u
        ├─ limit_from_rational_laurent
        ├─ series_at_zero（ordre 递增，对齐 mrv ordre 循环）
        ├─ limit_at_zero_rational_lead（num 级数 / den 常数）
        └─ limit_from_scaled_finite（pump 次数）
```

**删除**步骤 2–4 的共轭 / sqrt 商快路径；有理式首项可保留为 fallback 内第一层（成本低），或并入 Laurent（更干净）。

与 upstream 映射：

| 阶段 | giac-rs | upstream |
|------|---------|----------|
| A | `limit_preprocess_plus_infinity` + `surd2pow` | `surd2pow` + `limit_symbolic_preprocess` |
| B | `limit_unidirectional_plus_infinity` | `unidirectional_limit` |
| C | `x=1/u` + `series_at_zero` | 级数 @ `+∞` 或有限点 `inv` 换元后的级数 |

---

### 3. 分阶段实施方案

#### Phase 0：基线测试（移除前）

在 `limit_engine/asymptotic::tests` 或独立模块增加 **仅走 fallback** 的探测测试：

```rust
// 强制跳过 MRV 与共轭，只测 x=1/u + series
limit_at_plus_infinity_fallback(expr, var, ctx)
```

| 用例 | 期望 | 当前 fallback-only |
|------|------|---------------------|
| `x*(sqrt(1+x²)-x)` | `1/2` | 待验证 / 补强 `series_div` |
| `(x+1)/(x-1)` | `1` | ✅ Laurent |
| CK-INT-58 `(x+1)/sqrt((x+1)/(x-1))` | `+∞` | ❌ 需 `sqrt` 级数链 |
| `x*atan(x)/(x+1)` | `pi/2` | ❌ `series_div` + `pi` |
| `(1+1/x)*(sqrt(x+1)+1)` | `+∞` | ❌ 待测 |

输出 gap 列表，作为 Phase 1 的工单。

#### Phase 1：补强通用级数（blocker）

在 `sparse_series.rs` / `limit_at_zero_fallback`：

1. **`series_div` 符号主项**：分母常数项非 `1`、分子含 `pi` 等符号常数时仍能取主项（`atan/(1+u)`）。
2. **`series_sqrt` 稳定性**：`(1+u²)^(1/2)` 二项展开足够阶；`series_div` 后消去 `u` 的负幂。
3. **ordre 递增**：`limit_at_zero_from_series` 失败时 `ordre *= 1.5` 重试（对齐 `mrv_series_lead_loop` / `series__SPOL1`）。
4. **`surd2pow`（新）**：`sqrt(e)→e^(1/2)`，`limit_preprocess_plus_infinity` 入口调用；与 upstream `limit` 中 surd 分支一致，避免共轭里散落 `is_sqrt` 判断。

不新增任何 `rationalize_sqrt_*`。

#### Phase 2：移除写死快路径

`limit_at_plus_infinity` 收敛为：

```rust
pub(crate) fn limit_at_plus_infinity(expr, var, ctx) -> Result<ExprArc, EvalError> {
    if mrv_limit_eligible(expr) {
        if let Ok(r) = limit_unidirectional_plus_infinity(expr, var, ctx) {
            if is_usable_limit(&r) { return Ok(r); }
        }
    }
    limit_at_plus_infinity_fallback(expr, var, ctx)
}
```

删除（或 `#[cfg(test)]` 对比保留一周）：

- `rationalize_sqrt_in_expr` 及 `rationalize_sqrt_*`
- `limit_conjugate_sqrt_at_infinity`
- `limit_rational_over_sqrt_quotient_at_infinity`

可选：`limit_rational_leading_at_infinity` 移入 `limit_at_zero_fallback` 前作为有理式快速 Laurent（非形状表）。

验收：

- `maxima_rtest` 已通过用例仍绿（尤其 `x*(sqrt(1+x²)-x)`、`CK-INT-58`）。
- `cargo test -p giac-calculus --lib` 全绿。

#### Phase 3：与 upstream 进一步对齐（非阻塞）

| 项 | 说明 |
|----|------|
| `series__SPOL1` @ `lim_point=+inf` | 直接渐近级数，可与 `x=1/u` 二选一或互为 fallback |
| `mrv_compare` | `3^x` vs `5^x`、`x^ln(x)` 等 |
| CK-60 / gruntz | 扩展 `remove_lnexp`、216e 级数 |
| `partfrac` @ `+∞` | 替代 `limit_rational_leading` 的代数兜底 |

---

### 4. 风险与回退

| 风险 | 缓解 |
|------|------|
| 移除共轭后 `sqrt` 用例回退 | Phase 0 先测 fallback-only；Phase 1 补齐再删 |
| 级数 ordre 不足 | ordre 递增 + `MAX_SERIES_EXPANSION_ORDER` 与 MRV 共用上限 |
| 性能（倒数换元 + 级数慢于共轭） | 可接受；共轭仅为 `O(1)` 形状优化，非架构必需 |
| CK-58 形状复杂 | 依赖 `surd2pow` + 级数，而非恢复 `limit_rational_over_sqrt_quotient` |

回退策略：保留 git 中共轭相关函数一个 commit，文档标注「仅 debug 对比，不参与生产路径」。

---

### 5. 建议落地顺序（摘要）

```
Phase 0  fallback-only 回归矩阵
   ↓
Phase 1  series_div(pi) + surd2pow + ordre 递增
   ↓
Phase 2  删除 rationalize_sqrt_* / limit_conjugate_* / limit_rational_over_sqrt_*
   ↓
Phase 3  mrv_compare、series @ +inf、gruntz
```

**原则**：新 `sqrt`/代数形状 → 只扩展 `sparse_series` 或 `preprocess`，**禁止**新增 `rationalize_*` 形状表。

