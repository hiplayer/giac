# GIAC-poly — 扩域系数 gcd/factor 任务优先级（上游对齐）

**状态:** open  
**类型:** 实施计划 / AFK  
**父项:** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) §6 **P3-1…P3-3**、§7 **P4-***  
**相关:** [GIAC-algext-adoption](GIAC-algext-adoption.md) §8.4、[GIAC-poly-p3-6-quartic-roots-gaps](GIAC-poly-p3-6-quartic-roots-gaps.md)（P3-6 ✅）、[expr-poly-conversion.md](../expr-poly-conversion.md)、[giac-poly-api-stability.md](../giac-poly-api-stability.md)  
**上游基线:** `giac/giac-2.0.0` — `gausspol.cc`（`gcd`→`gcd_ext`、`ext_factor`/`ext_factor_nodegck`）、`threaded.cc`（`mod_gcd_ext`）、`sym2poly.cc`（partfrac + `_EXT`）  
**Rust 落点:** `giac-poly::Poly<AlgExtCPolyCoeff>`、`giac-core::algebra::{field_session, poly, poly_alg_coeff}`  
**快照:** 2026-06-23（**T0-1…T0-3 ✅**、**T1-1…T1-3 ✅**）

---

## 1. 问题陈述

`Poly<AlgExtC>` 表示层（P1）与 deg≤4 求根（P3-6）已落地，但 **K 上 gcd/factor 主路径仍缺**，导致：

- `factor(x²−2)`、`gcd(x²−2, x−√2)` 无法走多项式管线
- `partfrac(1/(x²−2), x)` 卡在 disc>0 二次分裂
- `solve` / `integrate` 仍依赖 `giac-simplify` / `giac-solve` 的 rootof 形状特判

**目标：** 与 upstream `_EXT` 系数多项式语义对齐——**不在 `Poly<ℚ>` 或 per-case 形状表上打补丁**。

**已就绪（不重复做）：**

| 组件 | 状态 |
|------|------|
| `PolyCoeff for AlgExtCPolyCoeff`（含 `inv`） | ✅ |
| `Poly<C>` 稀疏环、`try_add`/`try_mul` | ✅ |
| `FieldSession` + `poly_algext_roots` deg≤4 | ✅ P3-6 |
| `poly_alg_from_expr` / `algext_poly_to_expr` | ✅ P1-4 |

---

## 2. 上游算法锚点

| 能力 | 上游入口 | 要点 |
|------|----------|------|
| gcd over K | `gausspol.cc` `gcd` → `pt==_EXT` → **`gcd_ext`** | 模 gcd / 子结果式；系数、余式 ∈ 同一 `_EXT` |
| factor over K | **`ext_factor`** → **`ext_factor_nodegck`** | 先 `common_EXT`/`ext_reduce`；sqff → d=1 线性 / d=2 二次分裂 / 高次子路径 |
| 二次分裂 | `ext_factor_nodegck` + `addtov` | `algebraic_EXTension` + √Δ；复根走 `AlgExtC` |
| partfrac | `sym2poly.cc` + `ext_factor` | 分母在 K 上 factor 后线性方程组 |
| **明确不做** | — | `_EXT` 系数 **无** Hensel / sparse_bi / unitaryfactor（与 [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) §6 一致） |

---

## 3. 任务优先级表

| 优先级 | ID | 任务 | 上游对标 | 依赖 | 验收 |
|:--:|:--:|---|---|---|---|
| **P0** | **T0-1** | **系数域对齐契约**：`Poly<AlgExtC>` 运算前经 `FieldSession::align`；同塔顶 K | `common_EXT` + `ext_reduce`（`gausspol.cc` ~6086–6115） | P1 ✅ | ✅ `align_algext_polys` / `ensure_common_field_for_polys` |
| **P0** | **T0-2** | **一元 `div_rem_wrt` / `quo_exact_wrt` 泛型化**：`Poly<C>` + `FlatUni<C>`，`C=AlgExtCPolyCoeff` | `_EXT` 系数 `quo`/`rem` | T0-1 | ✅ `rem(x²−2, x−√2)=0` |
| **P0** | **T0-3** | **首项归一化 `monic_wrt`**：leading coeff ∈ K 用 `inv` | `ext_factor_nodegck` ~6144–6148 | T0-2 | ✅ `monic_wrt_algext` lc=1 |
| **P1** | **T1-1** | **`egcd` / `gcd` over K**（一元优先） | `gcd_ext`（`threaded.cc` `mod_gcd_ext`） | T0-2 | ✅ `gcd(x²−2, x−√2)=x−√2`；`gcd(x²−2, x+√2)=x+√2`；`gcd(x²−2, x+1)=1` |
| **P1** | **T1-2** | **`content` / `primitive_part` / `square_free_part` over K** | `ext_factor` 前 sqff + `lcmdeno` + pp | T1-1 | ✅ sqff 链系数均在 K |
| **P1** | **T1-3** | **二次分裂 `split_quadratic_factor`**（disc≤0 / disc>0 / 复根） | `ext_factor_nodegck` d=2 + `addtov` | T0-1, P3-6 deg2 ✅ | ✅ `factor(x²−2)→(x−√2)(x+√2)` |
| **P2** | **T2-1** | **`factor_univariate_over_k` 主路径**：sqff → 一次 → 二次分裂 → K 内有理根 → 不可约 `[g]` | `ext_factor` + 次数校验（~6345–6351） | T1-1…T1-3 | `factor(x⁴−4)`；`(x²+1)(x²−2)` 降次 |
| **P2** | **T2-2** | **`factor_into` / `eval_factor` 接线** | `usual.cc` factor + `algext_convert` | T2-1 | `eval(factor(x²−2))` 展开 = 原式 |
| **P2** | **T2-3** | **solve 降次共用**（backlog P3-7 / P4-6） | `solve.cc` + `ext_factor` 递归 | T2-1, P3-6 ✅ | `(x²+1)(x³−x+1)` 五根；删 solve 形状特判 |
| **P2** | **T2-4** | **partfrac 消费 K 上 factor**（backlog P2-3/P2-4） | `sym2poly.cc` partfrac | T1-3, T2-1 | `partfrac(1/(x²−2),x)` 两项一次 |
| **P3** | **T3-1** | **多元 gcd（低维）**：子结果式 PRS 系数环 `PolyCoeff` 泛型 | `gcdheu` / `gcd_ext` 多元 | T1-1 稳定 | `gcd(x²−2·y, x−√2·y)` 类 |
| **P3** | **T3-2** | **高次不可约 witness**：sqff 后不可约 → `[g]`（或 K 上 Zassenhaus 子集） | `ext_factor` 高次（无 EXT-Hensel） | T2-1 | 与 upstream「不 silent 错分」一致 |
| **P3** | **T3-3** | **跨 crate 清理** | — | T2-2/T2-3 | 删 `try_factor_quadratic_rootof`、`biquadratic_rootof` 等；见 [GIAC-rs-four-crates-dedup-architecture](GIAC-rs-four-crates-dedup-architecture.md) C3-2 |

**与 backlog ID 映射：**

| 本 issue | [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) |
|----------|-----------------------------------------------------------|
| T0-* | P0-D 收尾、P1 表示层契约 |
| T1-* | **P3-1**、**P3-2** |
| T1-3, T2-* | **P2-5**、**P3-3**、**P3-7** |
| T2-4 | **P2-3**、**P2-4** |
| T2-2/T2-3/T3-3 | **P4-1…P4-7** |

---

## 4. 推荐实施波次

```text
波次 1（1–2 周）：T0-1 → T0-2 → T0-3 → T1-1
波次 2（1–2 周）：T1-2 → T1-3 → T2-1
波次 3（1 周）：  T2-2 → T2-3 → T2-4
波次 4（按需）：  T3-* + known-divergences 登记
```

依赖图：

```text
FieldSession / PolyCoeff ✅
    └─ T0-1 align
         └─ T0-2 div_rem_wrt
              └─ T0-3 monic
                   └─ T1-1 gcd
                        ├─ T1-2 sqff/pp
                        │    └─ T1-3 split_quadratic
                        │         └─ T2-1 factor_over_k
                        │              ├─ T2-2 eval factor
                        │              ├─ T2-3 solve 降次
                        │              └─ T2-4 partfrac
                        └─ T3-1 多元 gcd
```

---

## 5. 明确不做

1. **`Poly<AlgExtC>` 上 Hensel / sparse_bi / unitaryfactor** — upstream `_EXT` 亦不走（见 backlog §6）。
2. **通用五次根式** — deg≥5 不可约 → `rootof(α, P)` 一支（Abel–Ruffini）。
3. **参变元 A,B 塔** — Phase C / P3-5；不阻塞本计划 M2/M3。
4. **用 `Poly::div_rem`（多元 leading）验 K 上一元整除** — 须 `quo_exact_wrt` / `FlatUni<C>`（见 [GIAC-poly-nested-ring-types](GIAC-poly-nested-ring-types.md)）。

---

## 6. 纪律（避免走弯路）

1. **禁止** 在 `Poly<ℚ>` 或 `giac-simplify` 增 `try_factor_*_rootof` 形状表。
2. **禁止** 在 `partfrac.rs` 手写 rootof 分裂后仍用 `Poly<ℚ>` 解方程。
3. **必须** 含代数系数走 `poly_alg_from_expr`（路径 B）；`expr_to_poly` 拒绝 AlgExt（[expr-poly-conversion.md](../expr-poly-conversion.md)）。
4. 新增 API 标 tier；合入前按 [algorithm-expr-api.md](../algorithm-expr-api.md) §7.2 复审。

---

## 7. 验证门禁

```bash
cd giac-rs

# 波次 1
cargo test -p giac-core alg_ext poly_alg_coeff
cargo test -p giac-poly --lib

# 波次 2 起（示例单测名，落地时补）
cargo test -p giac-poly gcd_algext
cargo test -p giac-poly factor_x_squared_minus_2

# 波次 3 跨 crate
cargo test -p giac-conformance --test giac_check_cas
cargo test -p giac-conformance --test giac_check_integrate  # partfrac / C-02
```

| 波次 | 硬门禁 |
|------|--------|
| 1 | `rem(x²−2, x−√2)=0`；`gcd(x²−2, x−√2)=x−√2` |
| 2 | `factor(x²−2)` 乘积 = 原式；`factor(x⁴−4)` 分裂 |
| 3 | `eval(factor(x²−2))` conformance；`partfrac(1/(x²−2),x)` |

---

## 8. 任务 ID 速查

| 优先级 | IDs |
|--------|-----|
| P0 基础除法 | T0-1, T0-2, T0-3 |
| P1 核心环算法 | T1-1, T1-2, T1-3 |
| P2 管线闭环 | T2-1, T2-2, T2-3, T2-4 |
| P3 扩展与清理 | T3-1, T3-2, T3-3 |
