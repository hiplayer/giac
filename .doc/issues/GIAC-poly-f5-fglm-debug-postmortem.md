# GIAC — F5 FGLM / F-module √-reduction 调试 postmortem

**状态:** open（解决方案待落地）
**类型:** 调试复盘 / 重构提案
**相关:** [GIAC-poly-f5-fglm-over-coefficient-field](../issues_resolved/GIAC-poly-f5-fglm-over-coefficient-field.md)（P3 主体）、[GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md)（约定二元根因）、[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)
**快照:** 2026-06-29
**触发提交:** `giac-rs 6693e6b`（P3e 接线 + `LIFT_BUDGET` 回归防护）

---

## 0. TL;DR — 预防规则与 revisit trigger

### 0.1 预防规则（review / lint 可检查）

1. **任何 `fn` 取 `CoordsQ` 必须在签名或紧邻注释声明 order**（`// order: HighFirst` / `LowFirst`）；过渡期由 `scripts/lint-coordsq-order.sh` 扫描强制，newtype 落地后由类型本身强制。
2. **搜索伪装成纯函数的 `fn` 签名必须带 `Fuel` / `budget`**（§3.3）；裸穷举不进生产路径。
3. **`Option<…>` 参数不得同时编码控制流与数据**（§3.2）；递归 fn 的 base-case / 递归分支用显式 `enum` mode 分发。
4. **跨算术世界（sparse `Poly` ↔ dense `CoordsQ`）必须经 dense-poly1-refactor §4.2 转换 API**，禁止裸坐标互灌（D4 防线）。

### 0.2 revisit trigger（本 postmortem 何时关闭）

| 条件 | 动作 |
|---|---|
| `GIAC-dense-poly1-refactor` §3.1 CoordsQ 双边 newtype 落地 + #5/#7 反向误用编译失败验证 | 关闭 §A、§D |
| `sqrt_base_case` 接 `Fuel` 参数（§3.3 局部） | 关闭 §C |
| `sqrt_fmodule` 改 `FmoduleMode`（§3.2） | 关闭 §B |
| 全链 fuel 审计独立 issue 闭合（§3.3 follow-up F2） | 关闭 §C 残余 |

---

## 1. 本次调试遇到的错误清单

按出现顺序,7 个独立错误（P3c→P3e）:

| # | 症状 | 根因类别 | 定位耗时¹ | 修复 |
|---|------|---------|----------|------|
| 1 | p-adic 提升结果错 | 表示 | 30min–4h | `u` 只 `mod p` 一次、`mod p^k` 算术里零填充。改为每步重算 `u mod p^k = numer·denom⁻¹ mod p^k` |
| 2 | p-adic 提升错 / 除零 | 表示 | 30min–4h | 模多项式 `m_int` 只 `mod p`、首系数 `denom_lcm` 可能非 `p` 单元。改为每步 `mod pk2` + 过滤整除 `denom_lcm` 的素数 |
| 3 | 编译错 `Monomial::pow` 缺失 | 机械 | <30min | 加 `mon_x_pow` 辅助（重复乘法构造 `x^i`） |
| 4 | 编译错 `cannot move out of *c` | 机械 | <30min | BigInt 算术改 `&*c` / `.clone()` |
| 5 | **`padic_sq=true` 但 `verified=false` 悖论** | **约定错配** | >4h | `field.element_mul` 读 HighFirst,F-module 管线写 LowFirst,对同一个 `Vec` 解读相反。验证改用独立 low-first `Poly<Ratio>` 算术 |
| 6 | 10 个 release 测 100s 超时 | 搜索无上限 | 30min–4h | `sqrt_base_case` 对非平方穷举 `60 素数 × 8 组合`。加 `LIFT_BUDGET=8` 硬上限 |
| 7 | norm 预过滤误杀 A₄ 真平方 | **约定错配（同 #5）** | 30min–4h | `N_{K/ℚ}(u)` 早退数学正确,但递归调用处 low/high-first 错配 → krylov 看到「转置」元素 → norm 算错。移除该过滤 |

¹ **定位耗时**：从症状出现到定位根因的时间；修复写码 <30s 全列等时，不计入。档位 <30min / 30min–4h / >4h。

**关键观察:#5 与 #7 是同一个根因踩了两次。** 第一次花最多时间（远端 `verified=false` 症状),第二次又差点回归。

---

## 2. 根因分类

| 类别 | 占调试时间 | 性质 | 命中错误 |
|------|-----------|------|---------|
| **A. 约定不可见（intra-world）** | ~60% | 设计缺陷（主因） | #5, #7 |
| **B. 隐式双义参数** | ~15% | 接口设计 | 见 §3.2 |
| **C. 搜索函数无代价上限** | ~20% | 工程契约 | #6 |
| **D. 跨算术世界无桥接（inter-world）** | 0%（本轮未触发） | 边界设计 | —（dense-poly1-refactor D4 防线） |
| **E. 机械 Rust 摩擦** | ~5% | 正常日常 | #3, #4 |

> **A 与 D 正交**：A 是同一 `Vec` 内部两种解读（HighFirst ↔ LowFirst），本轮主因；D 是 sparse `Poly` ↔ dense `CoordsQ` 两种载体之间无转换函数，本轮未触发但为 dense-poly1-refactor D4 处理的独立维度。原版 §2 把 #5 的"另一面"计入 D 致两类重叠，现拆分：#5 全归 A，D 重定义为 inter-world。

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

### 2.4 D — 跨算术世界无桥接（inter-world，本轮未触发）

sparse `Poly`（升幂、`giac-poly::univariate`）与 dense `CoordsQ`（高次在前、`field_arith`）两套载体之间无显式转换函数；F-module 内部 `Poly<Ratio>` 与 `CoordsQ` 直接互灌处即此类隐患点。本轮 #5 是同一 `CoordsQ` 内部的 HighFirst/LowFirst 解读反转（intra-world，归 A），未触及 sparse↔dense 跨载体边界，故 D 本轮 0% 调试时间。但 `sqrt_fmodule` 的 `f_basis`/`uk`/`dp_kcoords` 组装处一旦经 sparse `Poly` 中转，即落入 D；dense-poly1-refactor §4.2 转换 API（D4）是 D 的根治防线，本 postmortem 不重复处理。

### 2.5 约定边界 inventory（§3.1 落地的影响面实证）

按 `field_arith.rs` + `poly_roots.rs:1458+` 核对的 HighFirst/LowFirst 生产者/消费者：

| 函数 | 位置 | 消费 order | 产出 order | newtype 必触点 |
|---|---|---|---|---|
| `field.element_mul` | ext_tower | HighFirst | HighFirst | 是 |
| `poly_add/sub/mul/neg/scale/div_rem/ext_gcd` | field_arith:45–93 | HighFirst | HighFirst | 是 |
| `poly_reduce` / `poly_inv_mod` | field_arith:71/76 | HighFirst | HighFirst | 是 |
| `poly_*_with_coeffs_in_field` / `poly_*_blocks` | field_arith:516–568 | HighFirst（`ParentBlockRing`） | HighFirst | 是 |
| `krylov_minpoly_coords` | poly_roots:1907 | （读 field） | **LowFirst, monic** | 是（关键 producer） |
| `poly_from_low` | poly_roots:1479 | **LowFirst** | `Poly`（sparse） | 是（函数名已编码约定） |
| `poly_to_low_bigint` | poly_roots:1539 | `PolyMod` | **LowFirst** | 是 |
| `sqrt_base_case(field, u_coords, m_gen_low)` | poly_roots:1458 | u/m_gen 均 **LowFirst** | **LowFirst** | 是 |
| `sqrt_fmodule(field, u, m_gen_low: Option)` | poly_roots:1887 | u **LowFirst**；`m_gen_low` LowFirst | **LowFirst** | 是 |
| `crt_combine_poly` / `ff_sqrt` | poly_roots | PolyMod 内部 | PolyMod | 否（纯 𝔽_p） |

**冒烟实例（#5 精确触发点）：** `sqrt_fmodule` 在 `poly_roots:1926` 把 `field.element_mul(&pw, u_coords)` 的返回值（HighFirst）push 进 `f_basis` 后当 LowFirst 用——`f_basis` 类型若改为 `LowFirstQ`，此行即编译失败。这是 §3.1 双边 newtype 的最小验证锚点。

---

## 3. 解决方案（按性价比排序）

### 3.1 ★ 根治 A:给 CoordsQ 加双边 order newtype

**方案（双边 newtype + 过渡 guardrail）：**

**第 0 步（过渡 guardrail，✅ 已落地 P1，见 [remediation](GIAC-poly-f5-fglm-postmortem-remediation.md)）：** `field.element_mul` / `krylov_minpoly_coords` / `poly_from_low` 标注 `// order:` 注释；`scripts/lint-coordsq-order.sh` 扫描 `fn.*CoordsQ` 要求相邻 `// order: HighFirst|LowFirst`，未标注的进 `scripts/coordsq-order-allowlist` ratchet（仅缩）。**注：** 原草案提的 `debug_assert_eq!(order, …)` 不可行——`CoordsQ = Vec<Ratio>` 无运行期 order 值，order 是隐式约定，运行期 assert 无从比较；这正是 P2 newtype 要解决的（把 order 编进类型，编译期抓）。故 P1 guardrail 是 lint + 注释，runtime assert 并入 P2。

**第 1 步（双边 newtype）：**

```rust
pub struct LowFirstQ(pub Vec<Ratio<BigInt>>);
pub struct HighFirstQ(pub Vec<Ratio<BigInt>>);
// 边界转换:rev()
impl LowFirstQ { pub fn from_high(h: &HighFirstQ) -> Self { … } }
impl HighFirstQ { pub fn from_low(l: &LowFirstQ) -> Self { … } }
```

- `field.element_mul` 只吃 `HighFirstQ`；F-module 管线只用 `LowFirstQ`；
- 编译期抓所有 HighFirst↔LowFirst 错配，#5/#7 类复发一次消除；
- 边界一对 `to_low`/`to_high`，把跨约定隐患局部化（F-module 内部纯 LowFirst，外部只从边界进出）。

**为什么双边而非单边 `LowFirstQ`：** 单边能抓本轮 100% 已知错配（方向固定：F-module 写 low、field 读 high），但 HighFirst 侧若未来新写代码出现 intra-world 错配则抓不到。双边 newtype 把 order 编码进**类型**，与 `field_arith` 现有 `POLY1_Q: Poly1Order = HighFirst` 常量参数（运行期传）**互补而非重复**——常量管 `dense::*` 内部一致性，newtype 管跨函数边界一致性，两者分层不冲突。churn 据表 §2.5 inventory 约 30–60 处 `CoordsQ` 裸用点，作为 `GIAC-dense-poly1-refactor` 首切片单独排期。

**影响面：** 见 §2.5 inventory；`krylov_minpoly_coords` 是 LowFirst monic minpoly 的唯一生产者，改其返回类型为 `LowFirstQ` 一个点即锁死下游一半约定。

**与 [GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md) §2.1「显式 poly1 顺序」一致**——可作为该重构的第一个落地切片（只做 CoordsQ 的 order newtype，不动 dense poly1 全量）。

**DoD（复用现有测，不重造 minimal reproducer）：**
- #5 场景 = `diag_fmodule_sqrt_recovery`（已 green，`6693e6b`）；
- #7 场景 = `quartic_a4_galois_dim_le_12`（已 unignore green）；
- newtype 落地后：两者仍 green **且** 反向误用（把 HighFirst `CoordsQ` 喂给 `sqrt_fmodule` 的 `f_basis`，即 §2.5 冒烟实例 `poly_roots:1926`）编译失败——用 `#[cfg(test)]` `compile_fail` 测或手工 `git revert` 锚点行验证。

### 3.2 修 B:`Option` → 显式 mode

```rust
enum FmoduleMode<'a> { Top, Recurse { m_gen: &'a CoordsQ } }
fn sqrt_fmodule(field, u, mode: FmoduleMode) -> Option<CoordsQ>
```

- 「顶层(禁基例)」与「递归(允许基例)」成为显式状态;
- norm 过滤这类「只在某模式合法」的逻辑有明确归属;
- `sqrt_base_case` 不再收 `m_gen_low`,改从 `field` 取生成元极小多项式,消除冗余不一致。

### 3.3 修 C:搜索函数显式 Fuel

**拍板：用 `Fuel`（非 `budget: u32`）。** crate 已有 `Fuel` 抽象（`EvalError` 管线、S7 `Fuel::new(8)`），无理由造第二个。

```rust
fn sqrt_base_case(field, u, m_gen, fuel: &mut Fuel) -> Option<CoordsQ>
```

- 把「这是搜索」写进签名，调用方控代价；
- 超限走 `None`/`Err` 而非静默穷举；
- 当前 `poly_roots:1476` `const LIFT_BUDGET: u32 = 8` 升级为调用方传入 `Fuel::new(N)`，`const` 删除。

**局部 vs 全链：** 本 postmortem 只做局部（`sqrt_base_case` 加 `Fuel`，封住已爆发的 #6）。全链 fuel 审计——60 素数外循环（`small_primes(60)`，magic number 应参数化）、Newton lift 步数（`target_bits=320` 决定，但 `p^k > 2^320` 的 k 无硬上限）、`factor_mod_irreducibles` 内部——也是无界搜索，但改 `sqrt_fmodule` 递归结构超出了"postmortem 修复"范畴，**记为 follow-up F2 开独立 issue**。

### 3.4 不改:E 类

#3(`Monomial::pow`)、#4(BigInt move)是正常 Rust 摩擦,不该当设计问题改。

---

## 4. 落地建议

- **P3 已交付**（`6693e6b`，A₄→dim 12 绿、回归封顶），**本 postmortem 不阻塞 P4**。
- §3.1（双边 order newtype + 过渡 guardrail）是性价比最高项，建议作为 `GIAC-dense-poly1-refactor` 的**首切片**单独排期：影响面见 §2.5 inventory，编译器驱动，DoD 见 §3.1 末。
- §3.2、§3.3 可与 P4 顺手做（都在 `sqrt_fmodule`/`sqrt_base_case`，改动局部）。
- 凡新写域元素算术，**禁止再裸用 `CoordsQ` 跨约定边界**——要么走 newtype，要么边界显式 `rev()` + 注释，否则归入本 postmortem §A 复发类。
- **关闭条件**：见 §0.2 revisit trigger 表。

### 4.1 Follow-ups

| ID | 内容 | 触发条件 |
|----|------|----------|
| F1 | **norm 预过滤复活**：#7 移除的 `N_{K/ℚ}(u)` 早退数学正确，仅因 §A 约定错配误杀。§3.1 双边 newtype 落地后约定错配消失，norm 预过滤可安全复活，是纯 perf 收益（早退非平方 u，省 `LIFT_BUDGET` 次 lift）。 | §3.1 合并后评估 |
| F2 | **全链 fuel 审计**：`sqrt_fmodule` 递归链上的 60 素数外循环、Newton lift 步数、`factor_mod_irreducibles` 内部均无界搜索，§3.3 局部只封 `sqrt_base_case`。 | 开独立 issue |
