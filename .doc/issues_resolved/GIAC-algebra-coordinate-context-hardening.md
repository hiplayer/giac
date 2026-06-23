# GIAC algebra — 坐标上下文与 embedding hardening

**状态:** **resolved**（2026-06-23；P0–P2 已验收）  
**类型:** correctness / API hardening  
**相关:** [GIAC-poly-quartic-roots-F1-F5.md](../issues/GIAC-poly-quartic-roots-F1-F5.md)、[giac-tower-common-math.md](../giac-tower-common-math.md)、[giac-core-algebra-api-stability.md](../giac-core-algebra-api-stability.md)  
**落点:** `giac-core::algebra::{ext_tower,field_session,poly_roots}`

**长期规则（仍有效）：** §2 开发规则适用于新增 embedding / 坐标入口；新 PR 改 `ext_tower` / `poly_roots` 时对照 §2 复审。

---

## 1. 问题陈述

`ExtensionField + CoordsQ` 是一个过宽 API：`CoordsQ` 不携带所属 field、basis/layout、是否已按目标 field 规范化等上下文。调用者一旦用 `dimension()` 或下标手写坐标，就容易把 primitive 幂基假设带进 tower block layout。

本次 M3/F3 根因：

```text
rational_subfield_embedding(Q, K_tower)
```

旧实现把 ℚ 的 `1` 放到目标坐标最后一维；这只适合 primitive over ℚ。parent-coeff tower 的常数应在 block0，应由 `target.embed_rational(1)` 生成。错误嵌入导致 `eval_vanishes` 把最高项系数 `1` lift 到 dim24 后落在错误坐标，四次根误判为不零化。

---

## 2. 开发规则（长期有效）

1. **目标 field 负责构造坐标**  
   常数、零、一、生成元只能由 `target.{embed_rational,zero_coords,one_coords,generator_coords}` 生成。禁止用 `dim - 1`、block offset 等手写“常数在哪”。

2. **embedding 必须保持环结构**  
   每个 embedding 至少测试：
   - `emb(0) == target.zero`
   - `emb(1) == target.embed_rational(1)`
   - generator/minpoly 关系
   - 小样本 `emb(a+b)=emb(a)+emb(b)`、`emb(ab)=emb(a)emb(b)`

3. **L2/L3 不传裸 `CoordsQ` 作为业务参数**  
   裸坐标只留在 L1；跨 field 进入算法必须包在 `AlgExtCPolyCoeff` 或经 `FieldSession::{lift,align}`。

4. **按 `LayerArithMode` 分派，不按维数猜**  
   `dimension()==24` 不能说明 primitive/tower/flatten。危险函数要么统一调用 field constructor，要么显式 match `layer_arith_mode()`。

5. **fallback 不能返回错根**  
   所有 adjoin/deflate fallback 必须以 `roots_all_vanish` / `eq_mod(0)` 作出口门禁。

---

## 3. 优先级（已全部完成）

| 优先级 | 项 | 状态 | 验收 |
|--------|----|------|------|
| **P0** | 修 `rational_subfield_embedding` 使用 `target.embed_rational(1)` | ✅ | `rational_embedding_into_tower_uses_constant_block` 绿；F3 e2e 无 ignore |
| **P0** | `eval_vanishes` / `verify_root` 用 `eq_mod(0)` 判零 | ✅ | dim24 上不再依赖坐标全零 |
| **P0** | Quartic fallback 出口加 `roots_all_vanish` guard | ✅ | fallback 不返回错根 |
| **P1** | embedding ring-homomorphism 测试矩阵 | ✅ | `embedding_ring_hom_*` 三条 + constant_block 两条 |
| **P1** | API stability 登记 | ✅ | `giac-core-algebra-api-stability.md` |
| **P2** | `CoordsInFieldForTest` 轻量包装（L1 测试试点） | ✅ | `ext_tower` 测试内 |
| **P2** | 扫描手写坐标下标 | ✅ | `embed_rationals_into` / `alg_ext_c` 已收敛 |

---

## 4. 已落地修复（归档）

- `ext_tower::{rational_subfield_embedding,embed_rationals_into}` 改为以 target `embed_rational(1)` 填 embedding matrix。
- 新增 `rational_embedding_into_tower_uses_constant_block`、`common_rational_to_tower_embedding_uses_constant_block`。
- 新增 embedding ring-homomorphism 测试矩阵：rational→tower、parent→child、composite chain。
- 新增 `CoordsInFieldForTest`（测试专用，非生产 API）。
- `poly_roots::{eval_vanishes,verify_root}` 改为 `eq_mod(0)`；`verify_root` 与 `PolyInK::prepare` 同路径。
- `quartic_roots_by_adjoin_deflate`：仅 `roots_all_vanish` 时返回。
- `euler_four_roots_vanish`、`roots_quartic_t4_plus_t_plus_1` 无 `#[ignore]` 且全绿。

---

## 5. 未纳入本 issue 的后续（可选）

| 方向 | 说明 |
|------|------|
| 生产 `CoordsInField` | P2 仅测试试点；全面 L1 包装另开 issue |
| 新 embedding 入口 | 按 §2 补 ring-hom 测试 + stability 登记 |
