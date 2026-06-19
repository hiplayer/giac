# GIAC-poly-unitaryfactor — `unitaryfactor` / `pzadic` 缺口与优先级

**状态:** open（P0/P1/P2a/P2b 已落地；P2c API 复审待做）  
**类型:** 算法 / FAC-G1 尾部  
**上游基线:** `giac/giac-2.0.0` `gausspol.cc` `unitaryfactor` / `pzadic` / `unitarize` / `do_factor_hensel` L7044–7077  
**实现:** `giac-rs/crates/giac-poly/src/factor/unitary.rs`  
**相关:** [GIAC-simplify-poly-upstream-gaps](GIAC-simplify-poly-upstream-gaps.md) §2.2、[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)、[giac-poly-api-stability.md](../giac-poly-api-stability.md) §FAC-G1  
**快照日期:** 2026-06-19

---

## 1. 问题陈述

giac-rs 已接入 `factor/unitary.rs`（`UnitaryEvalPoint`、`PzadicLift`、`PzadicDraft`/`LiftedFactor`、`unitary_factor_rev`），并在 `poly_uni` 的 sqff 链中排在 **sparse → Hensel 之后**。

**现状：** conformance / testfactor 主路径已绿；**P0 架构统一**、**P1 上游尾链**、**P2a 多点 coeff 插值**、**P2b sparse_bi sum-coeff** 均已落地。line 25 门禁（`L22+y³`）已取消 `#[ignore]`。

**剩余缺口：** `reverse()` 在 `factor_multivariate` 边界未接线（U5）；`ununitarize` 往返非精确（y² 符号）；3+ aux 的 `try_sparse_factor_bi` 仍仅 2-aux。

---

## 2. 已交付

| 项 | 说明 |
|----|------|
| 类型分层 | `UnitaryEvalPoint`、`PzadicLift`、`PzadicDraft`、`LiftedFactor`、`UnivariateIn::divides` |
| 管线接线 | `poly_uni`：`|others|≥1` 统一 `try_unitary_factor` → `unitary_factor_rev`（无二元旁路） |
| 赋值搜索 | `x0 = 2‖p‖∞+2`；sqff 微调；`x0 ← x0·73794/27011`；**无** `2..N` 小整数扫描；base 位长 ≤256 |
| P1 尾链 | 忠实 `pzadic`；`unitarize`/`ununitarize`；`trunc1` 常数项尾部 |
| P2a | `lift_factor_multi_eval`（局部窗 `[base0-(need-1),…]` + monic 插值）；`try_lift_and_peel` fallback |
| P2b | `reconstruct_factor_dual_embed` sum-coeff；`sparse_factor_tri_var_sum_coeff` ✅ |
| 测试 | line25 gate ✅；P2a 专项 3 + e2e 1（见 §6）；154 passed / 4 ignored |

---

## 3. 开放问题（按根因）

### U1 — `PzadicLift` 与上游 `pzadic` 语义 — **部分闭合（P1+P2a）**

| 上游 | giac-rs |
|------|---------|
| `pzadic` 使 `dim+1`，base-`x₀` 数位展开 | `PzadicLift::pzadic` + P2a 多点 Lagrange fallback |
| 与 `reverse()` + 递归变量序咬合 | `vars_rev` 内部反序；边界 `reverse_var_order` 未接线（U5） |

单点 `pzadic` 对 `y³` 交叉项仍不足；P2a 在单点失败时用多点插值恢复（line 25 验收）。

---

### U2 — 二元 / 多元路径分裂 — ✅ P0 已闭合

单一入口 `unitary_factor_rev`；`lift_candidates` 统一 peel。

---

### U3 — `unitarize` / `ununitarize` — ✅ P1 已闭合

`try_unitary_factor` 第二路径：`unitarize → unitary_factor_rev → ununitarize`。`unitarize_extracts_leading_coeff` 覆盖 L22 首项。

**已知：** `ununitarize` 往返非精确（y² 符号）；仅提取 `an` 已测。

---

### U4 — `trunc1` 常数项分支 — ✅ P1 已闭合

`factor_constant_tail` / `factor_constant_tail_into` 在 `main` 次数归零后递归。

---

### U5 — 缺 `reverse()` 统一入口（**管线**）

上游对 sqff 块先 `pcur.reverse()` 再 `unitaryfactor`。giac-rs 仅用 `vars_rev`，3+ 元符号翻转风险未消。

---

### U6 — 单点赋值信息不足 — ✅ P2a 已闭合（增强路径）

`lift_factor_multi_eval`：在 `pzadic` 不整除时，对赋值因子各 `main` 系数做 Lagrange 插值（`MULTI_EVAL_MAX_SAMPLES=8`）。

---

### U7 — `sparse_bi` sum-coeff — ✅ P2b 已闭合

`reconstruct_factor_dual_embed` + monomial loop；`sparse_factor_tri_var_sum_coeff` ✅。`try_dilation_sparse_bi` 随机 dilation 已接线。

---

## 4. 回归门禁

| 测试 | 状态 | 意图 |
|------|------|------|
| `testfactor_unitary_bilinear_gate` | ✅ | 仿射提升最小覆盖 |
| `testfactor_line25_unitaryfactor_gate` | ✅ | Hensel 失败、P2a multi-eval |
| `unitary_factor_line25_l22_y3` | ✅ | 单元测试同上 |

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
| `try_unitary_factor` | ✅（b≈34 + multi-eval） |
| `factor_multivariate` | 2 因子 |

---

## 5. 优先级路线图

### P0 — 类型与单一入口 — ✅ 2026-06-19

### P1 — 上游尾链 — ✅ 2026-06-19

### P2 — 增强与并行缺口

| ID | 工作 | 验收 |
|----|------|------|
| **U-P2a** | 多点赋值 coeff 插值 | ✅ line 25 gate |
| **U-P2b** | `sparse_bi` sum-coeff 重建 | ✅ `sparse_factor_tri_var_sum_coeff` |
| **U-P2c** | `giac-poly-api-stability.md` / `unitary.rs` tier 复审 | 待做 |

---

### P2a 回归（无小整数扫描）

| 测试 | 验收 |
|------|------|
| `eval_base_stream_upstream_only` | 首基 = `initial`，次基 = `advance(initial)` |
| `p2a_sample_window_anchors_below_base0` | 局部采样窗低于 `initial` |
| `p2a_line25_upstream_trajectory` | 前 4 个 upstream 基 ≥1 次 partial peel |
| `unitary_factor_line25_l22_y3` + `testfactor_line25_unitaryfactor_gate` | 端到端 2 因子 |

---

## 6. 验证命令

```bash
cd giac-rs
cargo test --release -p giac-poly --lib
cargo test --release -p giac-poly line25 --lib
cargo test --release -p giac-poly sparse_factor_tri_var_sum_coeff --lib
cargo test --release -p giac-conformance --test giac_check_factor
```

---

## 7. 参考

- 上游：`giac/giac-2.0.0/src/gausspol.cc` L3893 `pzadic`、L6701 `unitaryfactor`、L6783 `unitarize`
- 实现：`giac-rs/crates/giac-poly/src/factor/unitary.rs`
- 接线：`giac-rs/crates/giac-poly/src/factor/poly_uni.rs` `factor_sqff_over_coeff_ring_ctx`
- 嵌套环类型：[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)
