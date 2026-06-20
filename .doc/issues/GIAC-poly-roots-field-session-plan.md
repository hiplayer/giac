# GIAC — `poly_algext_roots` FieldSession 改造方案

**状态:** draft（P3-6 续作）  
**依据:** [algorithm-expr-api.md](../algorithm-expr-api.md) §6（通用原则）、§6.3（`FieldSession`）；P3-6 review（2026-06）  
**父项:** [GIAC-poly-p3-6-quartic-roots-gaps.md](GIAC-poly-p3-6-quartic-roots-gaps.md)  
**落点:** `giac-core::algebra::{field_session,poly_roots}`  
**快照:** 2026-06-20

---

## 1. 目标

| 项 | 说明 |
|----|------|
| **G5 类型化** | ambient **K** + working **L** 由 `FieldSession` 承载，子算法不再裸传 `&Arc<ExtensionField>` |
| **消 align 补丁** | Cardano / resolvent / split 内禁止散落 `align_coeff` + `ring_int(K,…)` |
| **P3-6 验收** | `t⁴+t+1=0` 四根、`eq_mod` 验证，单测 **≤10s** |
| **回归不退化** | 现有 7 个 `poly_roots` 单测保持绿 |

**非目标（本方案不做）：** P4-6 solve 接线、P3-7 可约降次、`PolyCoeff::coeff_one` 全 crate 重构。

---

## 2. 现状基线（已实现，勿回退）

| 模块 | 状态 |
|------|------|
| T3+ `adjoin_irreducible_parent_coeffs` | ✅ |
| `generator_coords`（ℚ 扩张 `v[dim-2]=1`） | ✅ |
| `poly_algext_roots` 入口 + deg 1–4 骨架 | ✅ |
| 二次 / 双二次 + K₁ 混合系数 normalize | ✅（补丁式） |
| 纯三次 `t³−2` | ⚠️ 仅 1 实根 |
| 一般四次 `t⁴+t+1` | ❌ 超时，测试 ignore |

**技术债：** `infer_field` / `normalize_coeffs` / `lift_to_field` / 手写 `align_coeff` 分散在 `poly_roots.rs`；三次 Cardano、四次 split 仍绑入口 **K**。

---

## 3. 架构：两层类型

### 3.1 `PolyInK`（入口规范化，一次性）

```text
poly_algext_roots(p, var)
  → K = infer_field(p)           // 最大维非 base 域
  → p_k = normalize_coeffs(p, K) // 所有项系数 ∈ K
  → p_m = monic_univariate(p_k)
  → FieldSession::new(K)
  → roots_dispatch(session, p_m, var)
```

**规则：** 进入 `roots_dispatch` 后，多项式系数 **语义上全在 `session.ambient()`**；`Poly` 树本身仍不携带 K（类型层用调用约定 + 私有模块边界约束）。

可选：新增 `struct PolyInK { poly, ambient: Arc<ExtensionField> }` 私有包装，禁止从外部构造未 normalize 的多项式。

### 3.2 `FieldSession`（算法全程）

```rust
// giac-core/src/algebra/field_session.rs（新文件）

pub struct FieldSession {
    ambient: Arc<ExtensionField>,  // K，不变
    working: Arc<ExtensionField>,   // L，单调扩大
}

impl FieldSession {
    pub fn new(ambient: Arc<ExtensionField>) -> Self;
    pub fn ambient(&self) -> &Arc<ExtensionField>;
    pub fn working(&self) -> &Arc<ExtensionField>;

    // 常数（总在 L 上）
    pub fn zero(&self) -> AlgExtCPolyCoeff;
    pub fn one(&self) -> AlgExtCPolyCoeff;
    pub fn int(&self, n: i64) -> Result<AlgExtCPolyCoeff, EvalError>;
    pub fn half(&self) -> Result<AlgExtCPolyCoeff, EvalError>;

    // 嵌入 / 对齐
    pub fn lift(&self, c: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError>;
    pub fn align(&self, a: &AlgExtCPolyCoeff, b: &AlgExtCPolyCoeff)
        -> Result<(AlgExtCPolyCoeff, AlgExtCPolyCoeff), EvalError>;

    // 扩域 primitive（内部 bump working）
    pub fn adjoin_sqrt(&mut self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError>;
    pub fn adjoin_cbrt(&mut self, u: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError>;

    // 算术（自动 lift + align 到 working）
    pub fn add(&self, a: &AlgExtCPolyCoeff, b: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError>;
    pub fn mul(&self, a: &AlgExtCPolyCoeff, b: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError>;
    pub fn div(&self, a: &AlgExtCPolyCoeff, b: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError>;
    pub fn neg(&self, a: &AlgExtCPolyCoeff) -> Result<AlgExtCPolyCoeff, EvalError>;

    fn bump_to(&mut self, field: &Arc<ExtensionField>);
}
```

**`bump_to` 策略：**

1. 若 `ExtensionField::is_subfield_of(working, new)` → `working = new`
2. 否则若 `is_subfield_of(new, working)` → 不变
3. 否则 **优先 T2 链 embed**，避免立即 `common_over_q`（四次性能关键）

**`adjoin_sqrt` / `adjoin_cbrt`：** 包装现有 `algext_square_roots` / `algext_cube_root`（T3+ 路径），取生成元后 `bump_to(β.field)`，返回 **已 lift 到 L 的** `AlgExtCPolyCoeff`。

---

## 4. PR 竖切（修订）

在 [P3-6 PR 顺序](GIAC-poly-p3-6-quartic-roots-gaps.md#4-算法竖切建议-pr-顺序) 中 **插入 PR-B′**，后续 PR 均基于 session。

```text
PR-A   T3+ + generator_coords          ✅ 已完成
PR-B′  FieldSession + 二次/verify 迁移  ← 下一步
PR-C′  三次 Cardano on session + ω 分支
PR-D′  resolvent 公式 golden + 构造 on session
PR-E′  四次 split on session + t⁴+t+1 去 ignore
PR-F   P4-6 接线（不变）
```

### PR-B′：FieldSession 骨架 + 二次迁移

| 任务 | 文件 |
|------|------|
| 新增 `field_session.rs` | `int/one/zero/half/lift/align/bump_to` |
| `poly_roots` 二次路径改用 `&FieldSession` | 删 `ring_int(field,…)`、`lift_to_field`、`align_coeff` 自由函数 |
| 测试 helper `verify_root(session, p, root)` | 与入口共用 `normalize+monic` |
| 单测 | `field_session` 常数在 L 上；K₁ 二次回归 |

**DoD：** 现有 4 个二次/双二次相关测试绿；无 PR 内新增手写 `align_coeff`。

### PR-C′：三次 on session

| 任务 | 说明 |
|------|------|
| `one_cubic_root(&mut session, …)` | 每步 sqrt/cbrt 后 `session.adjoin_*` |
| 纯三次 `t³+a` | `cube_roots_via_unity`：在 L 上构造 ω，返回 3 根 |
| 一般三次 | deflate → `quadratic_roots(&mut session, …)` |
| 复数 `sqrt_disc` | `session` 上 `mul_i`（改走 `AlgExtCData::from_coords_q`） |

**DoD：** `roots_cubic_t3_minus_2` 恢复 **3 根** 且全 `verify_root`；Cardano 单步单测 ≤10s。

### PR-D′：resolvent golden

| 任务 | 说明 |
|------|------|
| 手算 / upstream 对照 | `t⁴+t+1` depressed 系数 → resolvent 系数写入单测 fixture |
| `build_resolvent_cubic` | 仅 session 算术，禁止裸 `field` |
| 可选 O2 | parent-block 稠密向量（**仅**当 sparse 仍超时） |

### PR-E′：四次 split + 性能

| 任务 | 说明 |
|------|------|
| `quartic_roots(&mut session, …)` | resolvent → `cubic_roots(session)` → 对每个 z `split_depressed_quartic` |
| `split_*` 内 `m = session.adjoin_sqrt(…)` | split 二次在 **同一 session**（L 已含 z,m） |
| dedup | `eq_mod` 前 `session.align` 到同一 L |
| 性能 | resolvent 只解 **1 个** z（或按 upstream 选判别式非零分支）；缓存 adjoin registry |
| 去 `#[ignore]` | `roots_quartic_t4_plus_t_plus_1` ≤10s |

**DoD：** 四根 verify；`cargo test -p giac-core poly_roots` 全绿无 ignore。

---

## 5. `poly_roots.rs` 函数签名迁移表

| 现签名 | 迁移后 |
|--------|--------|
| `quadratic_roots(p, var, &field)` | `quadratic_roots(session, p, var)` |
| `cubic_roots(p, var, &field)` | `cubic_roots(session, p, var)` |
| `one_cubic_root(p, var, &field)` | `one_cubic_root(session, p, var)` |
| `quartic_roots(p, var, &field)` | `quartic_roots(session, p, var)` |
| `split_depressed_quartic(..., &field)` | `split_depressed_quartic(session, ...)` |
| `ring_int/half/lift/align_coeff` | **删除**，迁入 `FieldSession` |
| `coeff_at(..., &field)` | `coeff_at(..., session.ambient())` 或 `session.zero()` |

**保留私有：** `infer_field`、`normalize_coeffs`、`monic_univariate`（仅入口）。

---

## 6. 四次超时专项（PR-E′ 内）

| 根因 | 对策 |
|------|------|
| Cardano 链多次 `common_over_q` flatten | session `bump_to` 优先 T2；同链 adjoin 复用 registry |
| 对 3 个 resolvent 根全 split | 先实现 upstream 等价「选一个 z + 共轭分支」；3 根全试作 follow-up |
| `align` 在不等域上触发 flatten | split 前 `session.lift` 全部 depressed 系数到当前 L |

**门禁：** 每个 primitive 单测（`adjoin_sqrt`、`one_cubic_root` on 固定 resolvent）独立 ≤1s。

---

## 7. 测试策略

| 层 | 内容 |
|----|------|
| **L0 塔** | `generator_coords`、T3+ adjoin、subfield embed（已有） |
| **L1 session** | `K→L bump`；`int(2)` 在 adjoin 后仍 `eq_mod` 正确 |
| **L2 primitive** | `adjoin_sqrt(4α)` on K₁；`adjoin_cbrt(2)` on ℚ |
| **L3 公式** | 二次/三次/四次 end-to-end + `verify_root` |
| **L4 golden** | resolvent 系数向量；`t⁴+t+1` 四根 |

**超时规则：** 单测 `timeout 10`；超则改算法/域策略，不加 ignore（除已知 giac 上游 bug 登记 `known-divergences.md`）。

**共享 verify：**

```text
verify_root(session, p_raw, root):
  p = monic(normalize(p_raw, session.ambient()))
  eval in root.field via session.lift
```

---

## 8. 风险与缓解

| 风险 | 缓解 |
|------|------|
| `FieldSession` API 膨胀 | 仅暴露 roots 所需；算术 `add/mul/div` 可 Phase 2 |
| `&mut session` 传递链长 | 子函数返回 `(session, value)` 或 `SessionGuard` 仅 PR-B′ 用 `&mut` |
| 与 `PolyCoeff::coeff_one` 冲突 | 本模块禁止 `Poly::ring_one` 构造 K 上 poly；用 `session.one()` 造常数项 |
| resolvent 公式错 | PR-D′ 先 golden 再写 split |

---

## 9. 完成定义（P3-6 + G5）

- [ ] `field_session.rs` + 单测
- [ ] `poly_roots` 无自由 `align_coeff` / `ring_int(ambient, …)`（grep 为 0）
- [ ] `roots_cubic_t3_minus_2`：3 根
- [ ] `roots_quartic_t4_plus_t_plus_1`：4 根，无 ignore，≤10s
- [ ] [GIAC-poly-p3-6-quartic-roots-gaps.md](GIAC-poly-p3-6-quartic-roots-gaps.md) G0–G5 行更新
- [ ] [algorithm-expr-api.md](../algorithm-expr-api.md) §6.3 草图改为 **Stable API** 表（实现后）

---

## 10. 交叉引用

| 文档 | 关系 |
|------|------|
| [algorithm-expr-api.md §6.3](../algorithm-expr-api.md) | `FieldSession` 设计规则来源 |
| [GIAC-poly-p3-6-quartic-roots-gaps.md](GIAC-poly-p3-6-quartic-roots-gaps.md) | 缺口 / 验收 |
| [expr-poly-conversion.md](../expr-poly-conversion.md) | path B 混合系数 |
| [giac-tower-common-math.md](../giac-tower-common-math.md) | T2/T4a common 策略 |
