# GIAC — F5 FGLM / F-module √-reduction 调试 postmortem

**状态:** open（解决方案待落地）
**类型:** 调试复盘 / 重构提案
**相关:** [GIAC-poly-f5-fglm-over-coefficient-field](GIAC-poly-f5-fglm-over-coefficient-field.md)（P3 主体）、[GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md)（约定二元根因）、[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)
**快照:** 2026-06-29
**触发提交:** `giac-rs 6693e6b`（P3e 接线 + `LIFT_BUDGET` 回归防护）

---

## 1. 本次调试遇到的错误清单

按出现顺序,7 个独立错误（P3c→P3e）:

| # | 症状 | 根因类别 | 修复 |
|---|------|---------|------|
| 1 | p-adic 提升结果错 | 表示 | `u` 只 `mod p` 一次、`mod p^k` 算术里零填充。改为每步重算 `u mod p^k = numer·denom⁻¹ mod p^k` |
| 2 | p-adic 提升错 / 除零 | 表示 | 模多项式 `m_int` 只 `mod p`、首系数 `denom_lcm` 可能非 `p` 单元。改为每步 `mod pk2` + 过滤整除 `denom_lcm` 的素数 |
| 3 | 编译错 `Monomial::pow` 缺失 | 机械 | 加 `mon_x_pow` 辅助（重复乘法构造 `x^i`） |
| 4 | 编译错 `cannot move out of *c` | 机械 | BigInt 算术改 `&*c` / `.clone()` |
| 5 | **`padic_sq=true` 但 `verified=false` 悖论** | **约定错配** | `field.element_mul` 读 HighFirst,F-module 管线写 LowFirst,对同一个 `Vec` 解读相反。验证改用独立 low-first `Poly<Ratio>` 算术 |
| 6 | 10 个 release 测 100s 超时 | 搜索无上限 | `sqrt_base_case` 对非平方穷举 `60 素数 × 8 组合`。加 `LIFT_BUDGET=8` 硬上限 |
| 7 | norm 预过滤误杀 A₄ 真平方 | **约定错配（同 #5）** | `N_{K/ℚ}(u)` 早退数学正确,但递归调用处 low/high-first 错配 → krylov 看到「转置」元素 → norm 算错。移除该过滤 |

**关键观察:#5 与 #7 是同一个根因踩了两次。** 第一次花最多时间（远端 `verified=false` 症状),第二次又差点回归。

---

## 2. 根因分类

| 类别 | 占调试时间 | 性质 | 命中错误 |
|------|-----------|------|---------|
| **A. 约定不可见** | ~45% | 设计缺陷（主因） | #5, #7 |
| **B. 隐式双义参数** | ~15% | 接口设计 | 见 §3.2 |
| **C. 搜索函数无代价上限** | ~20% | 工程契约 | #6 |
| **D. 两套算术世界无适配层** | ~15% | 边界设计 | #5 的另一面 |
| **E. 机械 Rust 摩擦** | ~5% | 正常日常 | #3, #4 |

### 2.1 主因 A — 约定没编码进类型

`CoordsQ = Vec<Ratio<BigInt>>`（`field_arith.rs:21`),域的顺序 `POLY1_Q = Poly1Order::HighFirst`(`:24`),F-module 管线全程 low-first。**ordering 只在注释里**(如 `poly_roots.rs:1140` `// low-first, monic`、`:1545`）。

后果:
- `field.element_mul`(HighFirst)与 F-module `Poly`(LowFirst)对同一 `Vec` 解读相反,编译器抓不到;
- 只能靠远端症状(`verified=false`、norm 误杀)暴露,定位链长;
- 同一坑踩两次(#5、#7)。

这与 [GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md) §1 记录的「高次在前 ↔ 升幂」二元是**同一根因**——那份重构计划尚未落地,F-module 便在未统一的约定上踩坑。

### 2.2 B — 隐式双义参数

`sqrt_fmodule(field, u_coords, m_gen_low: Option<&CoordsQ>)` 的 `Option` 同时编码:
- 「我是不是递归顶层」(控制流);
- 「有没有生成元极小多项式」(数据)。

两个正交含义挤在一个 `Option` 里,使「哪条路径会执行」难以推理;#7 的 norm 过滤放错位置(放进了递归也走到的基例)即由此而来。

`sqrt_base_case(field, u_coords, m_gen_low)`:`field` 已隐含生成元极小多项式,又单独传 `m_gen_low`,**两者可不一致却无校验**——冗余且危险。

### 2.3 C — 搜索伪装成纯函数

`sqrt_base_case` 本质是搜索(穷举素数×符号组合),却伪装成纯函数,签名上看不出「可能很慢/不终止」。#6 的 100s 超时是必然:非平方在模 p 经常「是平方」,搜索不收敛。这与宽窄无关,是「搜索 vs 计算」没区分。

### 2.4 D — 无边界适配层

`field.element_mul`(HighFirst)与 `Poly<Ratio>`(LowFirst)在 F-module 内部直接混用,边界没有显式转换函数。#5 的修复(验证改独立 low-first 算术)是**绕过**而非**根治**——约定错配仍潜伏在 `sqrt_fmodule` 的 `f_basis`/`uk`/`dp_kcoords` 组装处(只是恰好自洽,见 diag 通过)。

---

## 3. 解决方案（按性价比排序）

### 3.1 ★ 根治 A+D:给 CoordsQ 加 order newtype

最小版:

```rust
pub struct LowFirstQ(pub Vec<Ratio<BigInt>>);
pub struct HighFirstQ(pub Vec<Ratio<BigInt>>);
// 边界转换:rev()
impl LowFirstQ { pub fn from_high(h: &HighFirstQ) -> Self { … } }
impl HighFirstQ { pub fn from_low(l: &LowFirstQ) -> Self { … } }
```

- `field.element_mul` 只吃 `HighFirstQ`;F-module 管线只用 `LowFirstQ`;
- 编译期抓所有错配,#5/#7 类复发一次消除;
- 边界一对 `to_low`/`to_high`,把 D 的隐患局部化(F-module 内部纯 LowFirst,外部只从边界进出)。

**与 [GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md) §2.1「显式 poly1 顺序」一致**——可作为该重构的第一个落地切片(只做 CoordsQ 的 order newtype,不动 dense poly1 全量)。

### 3.2 修 B:`Option` → 显式 mode

```rust
enum FmoduleMode<'a> { Top, Recurse { m_gen: &'a CoordsQ } }
fn sqrt_fmodule(field, u, mode: FmoduleMode) -> Option<CoordsQ>
```

- 「顶层(禁基例)」与「递归(允许基例)」成为显式状态;
- norm 过滤这类「只在某模式合法」的逻辑有明确归属;
- `sqrt_base_case` 不再收 `m_gen_low`,改从 `field` 取生成元极小多项式,消除冗余不一致。

### 3.3 修 C:搜索函数显式 budget

```rust
fn sqrt_base_case(field, u, m_gen, budget: u32) -> Option<CoordsQ>
```

- 把「这是搜索」写进签名,调用方控代价;
- 或统一接 `Fuel`(与 `EvalError` 管线一致),让超限走 `Err` 而非静默穷举;
- 当前 `LIFT_BUDGET=8` 是硬编码兜底,应升级为参数。

### 3.4 不改:E 类

#3(`Monomial::pow`)、#4(BigInt move)是正常 Rust 摩擦,不该当设计问题改。

---

## 4. 落地建议

- **P3 已交付**(`6693e6b`,A₄→dim 12 绿、回归封顶),**本 postmortem 不阻塞 P4**。
- §3.1(order newtype)是性价比最高项,建议作为 `GIAC-dense-poly1-refactor` 的**首切片**单独排期:影响面 `field_arith` + `ext_tower` + `poly_roots::sqrt_fmodule`,编译器驱动,可小步验证。
- §3.2、§3.3 可与 P4 顺手做(都在 `sqrt_fmodule`/`sqrt_base_case`,改动局部)。
- 凡新写域元素算术,**禁止再裸用 `CoordsQ` 跨约定边界**——要么走 newtype,要么边界显式 `rev()` + 注释,否则归入本 postmortem §A 复发类。
