# GIAC-poly P0  backlog — 架构 + 算法剩余项

**状态:** open  
**类型:** 索引 / AFK  
**相关:** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md)、[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)、[GIAC-poly-unitaryfactor-gaps](GIAC-poly-unitaryfactor-gaps.md)、[GIAC-simplify-poly-upstream-gaps](GIAC-simplify-poly-upstream-gaps.md) §2  
**快照:** 2026-06-19

---

## 1. 表示层 / 架构（nested-ring）

| 项 | 目标 | 状态 | 验收 |
|----|------|------|------|
| `factor/*` 热路径禁裸 `Poly::div_rem` | 嵌套环用 `UnivariateIn::divides` / `FlatUni` | **✅** | sparse/hensel/unitary 已迁；zassenhaus/fpx/cyclotomic 仅 `FlatUni` 或模域 |
| `BivariateEmbed` / `EmbedFactorDraft` | embed 阶段 IR，无 `Poly` 往返 | **✅** | `nested.rs` + `sparse.rs` reconstruct |
| `SqffRingCtx` / `PolyFactorTower` | sqff 上下文统一 | **✅** | `ctx.rs` + `tower.rs` + `poly_uni.rs` |
| `FlatUni` / `MultivariatePoly` | 平坦一元 vs 多元展示边界 | **✅** | `nested.rs`；`univariate.rs` / `util.rs` 分流 |

细节与 Phase 清单：[GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)（Phase 1–3 已勾选）。

---

## 2. P0 算法

| ID | 工作 | 状态 | 验收 / 备注 |
|----|------|------|-------------|
| **U5** | `factor_multivariate` 边界 `reverse()` 再 `unitaryfactor` | **Partial** | `|vars|≥3` 时对 `p`+因子做 `reverse_var_order`；二元仅用 `vars_rev`（reverse `p` 会破坏 line25 赋值点） |
| **sparse_bi 3+ aux** | `|aux|≥3` 时非直接 `None` | **Partial** | 已实现 **aux 两两配对** 调用 `try_sparse_factor_bi_two_aux`；完整三 aux 同 embed 仍缺 |
| **partfrac 高次** | 不可约次数 >3、实二次 disc>0 | **✅ disc>0** | deg≤3 不可约 ✅；disc>0 K 路由 ✅（DIV-084）；deg>3 单项式仍 `NotImplemented` |

### partfrac 缺口明细

| 情形 | 现状 |
|------|------|
| sqff 因子 `deg > 3` | `NotImplemented` |
| sqff 因子 `deg = 3` 不可约 | ✅（ℚ 单项式，如 `x/(x³+2)`） |
| 实二次 `disc > 0` 且 factor 不分裂 | ✅ K 路由（P4-2/3）；DIV-084 |
| 线性 / 重根 / `disc≤0` 二次 | ✅ |
| 仿射幂次线性方程组 | ✅（次数 ≤3） |

---

## 3. 优先级（剩余）

```text
P1  sparse_bi — 三 aux 同 embed sum-coeff（非仅配对）
P2  nested-ring — clippy 门禁禁 factor/* 嵌套环 div_rem（可选）
P2  F4′ — 四次 Galois σ(κ) 优化（见 algext-backlog §2.1）
```

---

## 4. 验证

```bash
cd giac-rs
cargo test -p giac-poly --lib unitary sparse_factor testfactor_line
cargo test -p giac-poly --lib partfrac
```
