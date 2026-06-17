# GIAC-216e — `series_lead_at_zero` 未收敛到统一 `SparseSeries` 路径

**状态:** closed  
**类型:** AFK（实现向）/ 部分 HITL（算法选型若需调整 MRV 语义）  
**Blocked by:** GIAC-216c+（`remove_lnexp` 已接入 `SparseSeries::add` / `normalize_map`）  
**相关:** GIAC-215 / GIAC-216、`mrv_lead_term.rs`、`sparse_series.rs`、`mrv_series_lead.rs`  
**验收用例:** CK-INT-61；`giac-calculus` `engine_ck_int_61` / `mrv_lead_ck_int_61`

---

## 问题陈述

`mrv_lead_term_plus_infinity` 在 MRV 换元后调用 `series_lead_at_zero`，期望 **仅通过** `SparseSeries` + upstream 风格 `remove_lnexp`（`padd` 系数合并）得到主项 `(exponent, coeff)`。

**现状：** 对嵌套 `exp(inner)`（`inner` 为含 `exp` 的分式）表达式，该路径 **不能** 稳定给出正确主项，最终 **回退** 到 `mrv_series_lead::mrv_lead_term_at_zero`（lead-only 特例）。

```text
series_lead_at_zero
  ├─ loop: series_at_zero + remove_lnexp(lead.coeff)
  │     └─ lead_coeff_ready?  (无 ln(w)、无 w 依赖)
  └─ 失败 → mrv_lead_fallback → mrv_lead_term_at_zero  ← CK-INT-61 仍走此路
```

这表示 **216c+ 的统一目标尚未完成**：`remove_lnexp` 已在 `padd` 就位，但 **MRV 主项提取** 与 upstream `mrv_lead_term` 仍有结构性差距。

---

## 失败机理（以 CK-INT-61 为例）

MRV 换元后（`w = exp(-x)`，`x → -ln(w)`）形状近似：

```text
(exp(inner') - w^-1) * (-ln(w))^-1
```

其中 `inner'` 仍为 `(-ln(w))*w / (w + exp(...))` 类嵌套分式。

| 步骤 | 现象 |
|------|------|
| 逐项 `series_at_zero` 再 `Add` | `exp(inner')` 与 `-w^-1` 同阶 `w^-1` 系数 `1` 与 `-1` 相加为 0，**丢失下一阶** |
| `remove_lnexp` 代数恒等式 `exp(f)-w^-1 = w^-1(exp(f+ln(w))-1)` | 对复杂 `inner'`，`f+ln(w)` 无法化成小量 ε，表达式几乎不变 |
| `series_at_zero` 完整展开 | `NotImplemented("series")`，或主项系数仍含 `(-ln(w))^-1` |
| `lead_coeff_ready` | 拒绝含 `ln(w)` 的主项 → 触发回退 |
| lead-only 回退 | `try_exp_minus_w_inv` / `second_term_inner_plus_ln` 等 **形状特例** 给出 `exp=-1, coeff=-exp(2)` |

---

## 已知实现缺口

1. **`try_order` 循环无效：** `series_lead_at_zero` 将 `try_order` 增至 `2*MAX_SERIES_ORDER`，但 `series_at_zero` 内 `order.min(MAX_SERIES_ORDER)` 恒 cap 为 10，**增高阶重试无实际效果**。
2. **`(-ln(w))^-1` 语义：** 即 `x^-1` 换元结果，非普通 `w^k` Laurent 项；当前以符号常数系数保留，未纳入完整 `padd` 消去链。
3. **符号 `exp` 保留不完整：** 虽对含 `ln(w)` 的 arg 保留 `exp(arg.to_expr(w))`，嵌套分式 arg 仍会在子展开中过早数值化。
4. **upstream 前置改写未移植：** `mrv_lead_term` 在 `series__SPOL1` 前有 `ln(exp(g)^k)` 重写、`upscale`、`ln(w)→±g` 代回、系数 `undef` 时 **升阶循环** 与 `spdiv` 求逆等（`series.cc` ~2853–2975）。

---

## 背景：Gruntz / MRV 标准理论

**Gruntz/MRV 理论**是符号计算中求 **x → ±∞ 极限** 的标准算法框架，由 Bruno Gruntz 在 1990 年代系统化（博士论文与后续论文，常称 **Gruntz algorithm**）。GIAC、Maple、Mathematica 的 `limit` 在 `+∞` 情形基本都沿此思路实现。

### 核心想法：比较谁“增长更快”

当 `x → +∞` 时，不同函数趋于无穷的速度不同，例如：

```text
exp(x)  ≫  x^100  ≫  ln(x)  ≫  1
```

算法先回答：**表达式里哪个子式是“最快变化”的（most rapidly varying, MRV）？** 然后用它做换元，把 `+∞` 问题化成 **在 0 附近的级数/主项问题**。

### MRV（Most Rapidly Varying）

对表达式 `f(x)`，在 `x → +∞` 时：

1. 找出所有显著依赖 `x` 且趋于无穷/0 的子式（如 `exp(-x)`、`exp(x)`、`x` 等）
2. 按渐近快慢排序
3. 选最慢趋于 0 的作为辅助变量 **w**（GIAC 里常为 `w = exp(-x)` 这类）

giac-rs 对应模块：

| 步骤 | 模块 |
|------|------|
| 收集 MRV 集 | `mrv.rs` → `mrv_at_plus_infinity()` |
| 选 w | `choose_mrv_w()` |
| 换元重写 | `mrv_lead_term.rs` → `rewrite_in_mrv_w()` |
| w=0 级数 | `sparse_series.rs` → `series_at_zero()` |
| 取主项 | `series_lead_at_zero()` |

### Gruntz 算法主流程（简化）

```text
limit(f(x), x, +∞)
  → 找 MRV 集
  → 选 w → 0
  → 换元: x, exp(±x), ... 用 w, ln(w) 表示
  → 在 w=0 做级数展开
  → 取主项 lead term
  → 主项系数仍依赖 x? → 递归再跑 MRV
  → 否则由主项符号/阶数得极限
```

要点：

- 不是一次 L'Hôpital，而是 **反复** MRV 换元 + 主项提取
- 级数系数可暂时含 `ln(w)`，最后 `rewrite_ln_w` 代回 `x`
- 两主导项同阶相消时，需 **提高展开阶** 或 `remove_lnexp` 类代数化简（CK-INT-61 卡点）

### 与普通 Taylor 级数的区别

| | 普通 Taylor | Gruntz/MRV 级数 |
|--|-------------|-----------------|
| 展开点 | 有限点（如 0） | `w → 0`，`w` 为 MRV 辅助变量 |
| 变量 | 原变量 `x` | 换元后的 `w` |
| 系数 | 数字或 `x` 的多项式 | 可含 `ln(w)`、`exp(...)` 等 |
| 目的 | 逼近函数值 | 比较阶数、求极限 |

`SparseSeries` 在 `_mrv_w` 处做 `series_at_zero`，即 Gruntz 框架中 **“MRV 换元后的局部展开”** 这一步。

### 为何称“标准理论”

Gruntz 证明：对 **有限组合初等函数**（`+,-,*,/,^,exp,ln,sin,cos,...` 等组成的塔），在适当条件下，MRV 比较 + 换元 + 主项提取 **能判定极限**（存在、无穷、或不存在）。

实务参照：

- **GIAC** `series.cc` 的 `mrv()`、`mrv_lead_term()` 即该理论的工程实现
- **Maple** `limit` 文档明确引用 Gruntz
- **Mathematica** 内部有类似 MRV/级数机制（实现细节不公开）

### 与 giac-rs `limit_engine` 的对应关系

| Gruntz 步骤 | giac-rs 现状 |
|-------------|--------------|
| MRV 集合 | `mrv.rs` ✅ |
| 选 w、换元 | `rewrite_in_mrv_w` ✅ |
| w=0 级数 | `SparseSeries` + `remove_lnexp` ⚠️ 部分 |
| 取主项、回代 | `series_lead_at_zero` ⚠️ 复杂式回退 lead-only |
| 升阶/消去循环 | upstream 有，Rust 未完整 ❌ |

**本 issue（GIAC-216e）要补的是 Gruntz 理论中“级数 + 主项 + 消去”那一段**，而非重新实现 MRV 集合本身。路线 **A**（对齐 upstream `mrv_lead_term`）是 Gruntz 理论的 **工程子集**；路线 **C** 指按论文完整实现含 `ln(w)` 系数域的 MRV-series。

**文献:** Gruntz, *On Computing Limits in a Symbolic Manipulation System* (1995)；van der Hoeven, transseries（超渐近/对数-指数级数扩展）。

---

## 候选算法与路线（按推荐优先级）

### A. 对齐 upstream `mrv_lead_term` 完整循环（推荐，风险低）

**来源:** `giac/giac-1.5.0/src/series.cc` `mrv_lead_term` + `series__SPOL1` + `padd`/`remove_lnexp`。

**要点:**

- 主项系数为 `undef` 或 lead 消去时，`ordre = ordre*1.5+1` **迭代升阶** 重算级数（非仅增大 `try_order` 参数名）。
- `series_flags & 0x1`（limit 模式）下 `padd` 对复杂系数调用 `remove_lnexp`（Rust 已部分实现）。
- 级数后 `subst(ln(w), ±g)` 将辅助变量 ln 换回真实 MRV 元。
- 系数仍 `undef` 时对级数做 `spdiv(1, p)` 再 `pnormal`（分式主项情形）。
- `exp` 为 MRV 元时：`ln(exp(g)^k*...) → k*g + ln(...)` 预重写。

**验收:** CK-INT-61 不经 `mrv_series_lead` 回退；删除或收缩 `mrv_lead_fallback`。

---

### B. 强化 `remove_lnexp` + 主项专用级数（中等工作量）

在 `remove_lnexp.rs` 下沉 lead-only 中已验证的代数：

- `exp(f) - w^-1` 恒等式（**已有**）。
- `f + ln(w)` 的 **次主项** 提取（`second_term_inner_plus_ln` 逻辑泛化）：对 `inner` 为分式时，从分母扰动 `den - w` 读下一阶。
- `divide_lead_coeffs` 型 `ln(w)` 幂次在分子分母间对消（`mrv_series_lead` 已有）。

配合 `SparseSeries`：**Add 不逐项展开**，先对整体 `remove_lnexp`，再单次 `series_at_zero`。

---

### C. Gruntz / MRV 理论标准路径（长期对照）

见上文 **「背景：Gruntz / MRV 标准理论」**。在已有 MRV 骨架上，补全系数域 `R(w, ln(w), exp(...))` 上的 MRV-series 与递归主项提取。

**优点:** 理论完备；**缺点:** 实现量大，超出当前 bounded `MAX_SERIES_ORDER=10` 引擎时需重新定义级数环。

---

### D. 对数-指数级数 / Transseries（过重，作备选）

对 `1/ln(w)`、`ln(w)^k` 等 **非整数幂 / 对数尺度** 项，扩展稀疏级数为：

```text
Σ c_{k,j} · w^k · (ln w)^j
```

或 Ecalle 超渐近级数。

**适用:** 极限比较阶精细分析；**超出** 当前 bounded `MAX_SERIES_ORDER=10` 引擎范围。

---

### E. 保持 lead-only 模块、缩小特例面（保守）

保留 `mrv_series_lead` 为 **主项提取器**，但：

- 将 `try_exp_minus_w_inv` 等改写为 `remove_lnexp` 的 **规则表**；
- `series_lead_at_zero` 显式文档化为「先尝试 SparseSeries，失败则 lead-only」；
- 不为 CK-INT-61 单独写形状路由。

**优点:** 快速稳定；**缺点:** 与 upstream 单路径 `mrv_lead_term` 仍有架构差。

---

## 建议实施切片（tracer bullets）

| 切片 | 内容 | 验收 |
|------|------|------|
| **216e-1** | 修复 `try_order`：升阶时真正提高 `MAX_SERIES_ORDER` 或内部 `ordre` 乘子循环 | 单元：人造 lead 消去后次项可现 |
| **216e-2** | 移植 upstream `mrv_lead_term` 升阶循环 + `ln(w)→g` 代回 | CK-INT-61 无回退 |
| **216e-3** | `ln(exp(g)^k)` 预重写 + `upscale` 路径 | 嵌套 `exp(-x)` 族回归 |
| **216e-4** | 删除 `mrv_series_lead` 中重复形状检测；收敛测试 103+ | `cargo test -p giac-calculus` |
| **216e-5**（可选） | `ln(w)` 混合系数域或 transseries 调研 spike | 设计笔记 / ADR |

---

## 验收标准（关闭本 issue）

- [x] `mrv_series_lead_loop` 对 CK-INT-61 **不调用** `mrv_lead_fallback` 即返回正确 `MrvLeadTerm`（`coeff=-exp(2)`）。
- [x] `mrv_series_lead.rs` 体积显著缩小，仅保留 `normalize_expr_quotients`。
- [x] `giac-calculus --lib` 全绿；无新增 unbounded `expand`。
- [x] `phase4-issues.md` 更新：CK-INT-61 不再标注为 lead-only 回退。

---

## 参考实现位置

| 组件 | Rust | Upstream giac |
|------|------|----------------|
| 稀疏级数 | `sparse_series.rs` | `sparse_poly1`, `series__SPOL1` |
| padd + remove_lnexp | `normalize_map`, `remove_lnexp.rs` | `padd`, `remove_lnexp` (`series.cc` ~263–328, ~2261) |
| MRV 主项 | `mrv_lead_term.rs`, `mrv_series_lead.rs` | `mrv_lead_term` (~2853) |
| MRV 集 | `mrv.rs` | `mrv()` |
| 有界参数 | `bounds.rs` `MAX_SERIES_ORDER=10` | `max_series_expansion_order` |

---

## 备注

提交 `5154682`（GIAC-216c++）已接入 `remove_lnexp`，本 issue 跟踪 **后续统一收敛** 工作，与 216c+「消除 CK-INT-61 硬编码形状路由」目标一致。
