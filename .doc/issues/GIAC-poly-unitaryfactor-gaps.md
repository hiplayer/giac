# GIAC-poly-unitaryfactor — `unitaryfactor` / `pzadic` 缺口与优先级

**状态:** open（P0 已落地，P1 待做）  
**类型:** 算法 / FAC-G1 尾部  
**上游基线:** `giac/giac-2.0.0` `gausspol.cc` `unitaryfactor` / `pzadic` / `unitarize` / `do_factor_hensel` L7044–7077  
**实现:** `giac-rs/crates/giac-poly/src/factor/unitary.rs`  
**相关:** [GIAC-simplify-poly-upstream-gaps](GIAC-simplify-poly-upstream-gaps.md) §2.2、[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)、[giac-poly-api-stability.md](../giac-poly-api-stability.md) §FAC-G1  
**快照日期:** 2026-06-19

---

## 1. 问题陈述

giac-rs 已接入 `factor/unitary.rs`（`UnitaryEvalPoint`、`PzadicLift`、`PzadicDraft`/`LiftedFactor`、`unitary_factor_rev`），并在 `poly_uni` 的 sqff 链中排在 **sparse → Hensel 之后**。

**现状：** conformance / testfactor 主路径已绿（L20–L24 等由 Tower / Hensel 覆盖）；**P0 架构统一已完成**（单一入口 `unitary_factor_rev`、统一 `lift_candidates`）；**P1 算法尾链**（忠实 `pzadic`、`unitarize`、`trunc1`、`reverse`）仍为 FAC-G1 尾部最大能力缺口。

**目标：** 对齐上游 **有界递归启发式**（非 closed-form），使 `unitaryfactor` 成为多元分解的真实兜底，而非形状特判集合。

---

## 2. 已交付（MVP）

| 项 | 说明 |
|----|------|
| 类型分层 | `UnitaryEvalPoint`、`PzadicLift`、`PzadicDraft`、`LiftedFactor`、`UnivariateIn::divides` |
| 管线接线 | `poly_uni`：`|others|≥1` 统一 `try_unitary_factor` → `unitary_factor_rev`（无二元旁路） |
| 赋值搜索 | `x0 = 2‖p‖∞+2`；sqff 微调；`x0 ← x0·73794/27011`；base 位长 ≤256 |
| 仿射特判 | `deg_x=1` 时 `±eval_var` 两候选 + `pzadic_digits` fallback |
| 测试 | `testfactor_unitary_bilinear_gate` ✅；147 passed / 6 ignored |

---

## 3. 开放问题（按根因）

### U1 — `PzadicLift` 非上游 `pzadic`（**架构**）

| 上游 | giac-rs |
|------|---------|
| `pzadic` 使 `dim+1`，在反序变量下做 base-`x₀` 数位展开 | flat `Poly` 上手搓 digit；无扩维语义 |
| 与 `reverse()` + 递归变量序咬合 | 固定 `(main, eval_var)` 二元视图 |

**后果：** 无法从赋值像恢复 `y²`、`y³`、`xy` 等系数结构；line 25 类用例必然失败。

**禁止方向：** 再加 `quadratic template`、`y³ template` 等 per-case 补丁。

---

### U2 — 二元 / 多元路径分裂（**架构**）— ✅ P0 已闭合

| 路径 | 行为 |
|------|------|
| ~~`try_unitary_factor_bivariate`~~ | 已删除；batch peel 内联为 `try_peel_all_at_eval` |
| `unitary_factor_rev` | 唯一入口；`lift_candidates` 统一 peel |

~~`poly_uni` 在 `|others|=1` 时**不走**完整递归~~ → 现与 `≥2` 相同，经 `try_unitary_factor`。

---

### U3 — 缺 `unitarize` / `ununitarize`（**上游链**）

上游 `do_factor_hensel` 在 `unitaryfactor` 失败后对 **非 monic 首项系数** 做：

```text
unitarize → reverse → unitaryfactor → reverse → ununitarize
```

典型：`3x - y² + y - 5`（L22 因子）。giac-rs **未实现**。

---

### U4 — 缺 `trunc1` 常数项分支（**上游链**）

上游 `unitaryfactor` 在 `main` 次数归零后，对常数项 `trunc1` 再递归分解。giac-rs 递归壳未实现此尾部。

---

### U5 — 缺 `reverse()` 统一入口（**管线**）

上游对 sqff 块先 `pcur.reverse()` 再 `unitaryfactor`，因子再 `reverse` 回来。giac-rs 在 `factor_multivariate` 边界未做等价反序，仅 `vars_rev` 在 `try_unitary_factor` 内部。

---

### U6 — 单点赋值信息不足（**数学**）

在 `eval_var = B` 时，真因子若含 `y²`、`y³`，赋值像为 **ℚ[x] 上常数系数**；**单点无法唯一恢复** `ℚ[x,y]` 中的高次项。

上游靠 **递归降维 + 多轮剥离** 逐步消化；若仍不足，需 **多点 coeff 插值**（增强，非主路径）。

---

### U7 — `sparse_bi` sum-coeff 仍缺（**FAC-G1 并行缺口**）

与 unitary 独立，但同属 `do_factor_hensel` 尾部能力：`sparse_bi` 重建仅单项式系数 MVP；sum-coeff / 完整 dilation 未补。部分多元式在到达 unitary 前即失败或极慢。

---

## 4. 回归门禁

| 测试 | 状态 | 意图 |
|------|------|------|
| `testfactor_unitary_bilinear_gate` | ✅ | 仿射提升最小覆盖 |
| `testfactor_line25_unitaryfactor_gate` | `#[ignore]` | Hensel 失败、需真实 pzadic/递归 |
| `unitary_factor_line25_l22_y3` | `#[ignore]` | 同上（单元测试） |

**Line 25 多项式（合成）：**

```text
f1 = 3*x - y^2 + y - 5
f2 = x*y + 3*x - y^2 - 1 + y^3
p  = f1 * f2
```

| 路径 | 结果 |
|------|------|
| `try_sparse_factor` (x,y)/(y,x) | None |
| `try_hensel_lift_bivariate` | None |
| `try_unitary_factor` / `unitary_factor_rev` | None（base 2..39 无 peel） |
| `factor_multivariate` | 1 因子（误判不可约） |

---

## 5. 优先级路线图

与 [GIAC-simplify-poly-upstream-gaps](GIAC-simplify-poly-upstream-gaps.md) §5 对齐；**unitary 属 P2**（不阻塞 conformance golden）。

### P0 — 类型与单一入口（架构）— ✅ 2026-06-19

| ID | 工作 | 验收 |
|----|------|------|
| **U-P0a** | `PzadicDraft` / `LiftedFactor` IR（`nested.rs`）；`lifted_factor_peel_vs_div_rem` 误用测 | ✅ |
| **U-P0b** | 删除 `try_unitary_factor_bivariate`；`|others|=1` 统一 `try_unitary_factor` | ✅ bilinear gate |
| **U-P0c** | `unitary_factor_rev` 统一 `lift_candidates`（无 `lift()` 单候选） | ✅ |

---

### P1 — 上游尾链（算法）— ✅ 2026-06-19（line 25 门禁仍 P2）

| ID | 工作 | 验收 |
|----|------|------|
| **U-P1a** | 忠实 `pzadic`（centered `smod` digit lift） | ✅ bilinear；`pzadic_lift_linear_via_digits` |
| **U-P1b** | `unitarize` / `ununitarize` | ✅ `unitarize_extracts_leading_coeff`；二次路径 |
| **U-P1c** | `trunc1` 常数项尾部 + `reverse_var_order` 辅助 | ✅ `factor_constant_tail`；`vars_rev` 接线（无符号翻转） |
| **U-P1d** | 移除 `±eval_var` 仿射特判 | ✅ 仅 `PzadicLift::pzadic` |

**未达里程碑：** `testfactor_line25_unitaryfactor_gate` 仍 `#[ignore]` — `y³` 交叉项需 **P2a** 多点插值（单点 `pzadic` + 尾链仍超时/失败）。

---

### P2 — 增强与并行缺口

| ID | 工作 | 验收 |
|----|------|------|
| **U-P2a** | 多点赋值 coeff 插值（仅当 P1 单点仍失败时） | line 25 或文档登记已知偏离 |
| **U-P2b** | `sparse_bi` sum-coeff 重建（FAC-G1，与 unitary 并行） | 减少进入 unitary 前的超时/误判 |
| **U-P2c** | `giac-poly-api-stability.md` / `unitary.rs` tier 复审 | API 表与实现一致 |

**估时：** 1–2 周（可与 P1 后半并行）

---

### 不做 / 低优先级

| 项 | 原因 |
|----|------|
| 更多仿射 / 二次 / 三次 template | 违背「一种环一种算法」；维护爆炸 |
| 展开六次积多元分解（line 23 结构层） | 属 simplify `factor(Mul)` + FAC-G2，非 unitary |
| 将 unitary 提至 P0 | conformance 已绿；Hensel/Tower 先覆盖 |

---

## 6. 验证命令

```bash
cd giac-rs
# 全量（含 6 ignored）
cargo test --release -p giac-poly --lib

# unitary 相关
cargo test --release -p giac-poly unitary --lib
cargo test --release -p giac-poly testfactor_unitary --lib

# line 25 门禁（实现 P1 后取消 ignore）
cargo test --release -p giac-poly line25 --lib -- --include-ignored

cargo test --release -p giac-conformance --test giac_check_factor
```

---

## 7. 参考

- 上游：`giac/giac-2.0.0/src/gausspol.cc` L3893 `pzadic`、L6701 `unitaryfactor`、L6783 `unitarize`
- 实现：`giac-rs/crates/giac-poly/src/factor/unitary.rs`
- 接线：`giac-rs/crates/giac-poly/src/factor/poly_uni.rs` `factor_sqff_over_coeff_ring_ctx`
- 嵌套环类型：[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)
