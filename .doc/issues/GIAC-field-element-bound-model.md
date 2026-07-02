# GIAC field-element bound model — 代数数绑定型表示演进

**状态:** open（P2a/P2b/P3a/P3b-1/P3c ✅ 落地 2026-07；P3b-3 同域换皮按决策跳过；P1 使用点硬化已落地，见 §1）
**类型:** architecture / correctness hardening
**相关:** [issues_resolved/GIAC-algebra-coordinate-context-hardening.md](../issues_resolved/GIAC-algebra-coordinate-context-hardening.md) §5、[giac-core-algebra-api-stability.md](../giac-core-algebra-api-stability.md)、[GIAC-algext-adoption.md](GIAC-algext-adoption.md) §7（flatten bisect fallback）
**落点:** `giac-core::algebra::{ext_tower,field_arith,common_minimal,field_session,poly_roots,alg_ext,alg_ext_c}`

---

## 1. 背景与问题陈述

`ExtensionField + CoordsQ` 是**分离表示**：域对象和坐标是两个独立值，没有语言层面绑定。

- `ExtensionField` 不知道某个 `CoordsQ` 的存在；`CoordsQ = Vec<Ratio>` 也不知道自己属于哪个域。
- 只在 `element_mul(field, &coords)` 这种**使用点**配对，靠约定（`coords.len() == field.dimension()` 且语义对）对齐。
- 坐标来源很多（域方法、字面量、`pad_to_len`、`factor_to_parent_blocks` 投影、cache 反序列化、跨域传递）；只有"从域方法取坐标"天然一致，外部塞的靠约定。
- 约定破了就错配：长度错（`block.len() != parent.dim`，被 `pad_to_len` 静默兜底成垃圾）、语义错（K₁ 坐标配 K₂ 域，dim 同但无意义）。

**已落地的使用点硬化**（前序 issue + 本次 P1）：
- `rational_subfield_embedding` 走 `target.embed_rational(1)`（归档 issue P0）。
- `is_embedded_rational` 修正蕴含方向，识别 parent-block 塔的嵌套有理常数位（本次）。
- `verify_embedded_generator` ParentBlocks 分支复合 `try_subfield_embedding(parent, operand) ∘ embedding_for(operand, pair)`，支持两 parent-block 塔 compositum（本次）。
- `build_adjoin_parent_coeffs` 形态校验（monic / block 长 = parent_dim / degree≥1，本次）。
- `CoordsInFieldForTest` 测试试点（归档 issue P2）。
- **P2a `ParentBlock` newtype（存储层先行，2026-07 落地）**：`min_poly_parent_blocks: Option<Vec<ParentBlock>>`，`ParentBlock::new(parent, coords)` 构造时挡 `len != parent.dim`；`build_adjoin_parent_coeffs` 为唯一持久存储写入门（裸 `Vec<CoordsQ>` → `Vec<ParentBlock>` 转换在此）；`LayerMinPolyForAdjoin::ParentBlocks` 升 `&[ParentBlock]` 零分配借存储；瞬态 compositum blocks（`MbSpec`、`factor_to_parent_blocks`、`common_minimal_poly_over_parent`）保持 `&[CoordsQ]` 由消费方把守。`semantic_key_bytes`/`layer_minpoly_parent_coeffs`/`flatten_min_poly_over_q_cold`/`verify_embedded_generator`/`common_adjoin_sibling_over` 解包 `as_coords()`。单测 `parent_block_new_rejects_wrong_length` + `build_adjoin_parent_coeffs_rejects_wrong_length_block`。

这些把错配挡在**使用点 / 登记点**，但坐标仍裸、配对仍靠约定——错配的**结构性源头**（独立对象可乱配）未消除。

## 2. 目标模型：FieldElement 绑定整体

对齐 upstream `gen`/`_EXT`：每个代数数从诞生起就是 `(field, coords)` 绑定整体，所有运算通过整体入口，内部 align 后产出新绑定整体。**无裸坐标流动 → 错配无入口**。

```rust
// 概念目标（非实现）
pub struct FieldElement {
    field: Arc<ExtensionField>,
    coords: HighFirstQ,   // 构造时校验 len == field.dimension()
}

impl FieldElement {
    pub fn new(field: Arc<ExtensionField>, coords: HighFirstQ) -> Result<Self, EvalError> {
        if coords.as_slice().len() != field.dimension() { return Err(...); }
        Ok(Self { field, coords })
    }
    pub fn mul(&self, other: &Self) -> Result<Self, EvalError> {
        // 内部：同域直接算；异域先 align 到公共域，结果包成新 FieldElement
        ...
    }
}
```

- **长度错配**：构造时挡。
- **语义错配**：`a.field != b.field` 拒绝直接算，强制走 align——比 `ParentBlock`（只挡长度）更彻底。
- **对齐 upstream**：`_EXT` 自带 `(coords, minpoly)`，`operator+` 内部 `common_EXT` 合并域、embed、运算、包新 `_EXT`。Rust `align_pair`/`common_cache` 已做同样的事，只是对齐责任在**调用方**；FieldElement 把责任移到**元素自身**。

## 3. 分层优先级

| 优先级 | 项 | 内容 | 代价 | 阻塞 conformance? |
|--------|----|------|------|-------------------|
| ~~**P2a**~~ ✅ | `ParentBlock` newtype（存储层先行） | `min_poly_parent_blocks: Option<Vec<ParentBlock>>`，构造时绑 parent + 校验 `len == parent.dim`；字段私有，唯一构造路径。**不动运行时元素表示**，只把"配对"从 `build_adjoin_parent_coeffs` 一个函数提到类型构造器，覆盖全路径（`factor_to_parent_blocks`、`common_adjoin_sibling_over`、test fixtures 等）。 | 局部类型重构，无 API 大改 | 否 — 2026-07 落地 |
| ~~**P2b**~~ ✅ | `MonicLayerMinPoly` newtype | 在 P2a 上再加一层：把层 minpoly 全部不变量（长度 + monic + degree≥1）集中到 `MonicLayerMinPoly::new(parent, blocks)`；`build_adjoin_parent_coeffs` 退化成一行。塔字段 `Option<Vec<ParentBlock>>` → `Option<MonicLayerMinPoly>`。 | 局部 | 否 — 2026-07 落地 |
| **P3** | `FieldElement` 运行时全栈 | `FieldElement { field, coords }` 作为运算单元；`ExtensionField::element_*` 方法族重构成 `FieldElement` 方法，align 内化；`AlgExtCPolyCoeff` 等系数类型对齐。消除所有裸坐标流动。 | **大改**（algebra 层入口换型）+ Arc 引用计数开销（upstream 同样接受）+ 与现有分离表示桥接 | 否 |
| ~~**P3a**~~ ✅ | `FieldElement` 类型 + 核心运算（align 内化） | 定义 `FieldElement { field: Arc<ExtensionField>, coords: HighFirstQ }`，构造挡长度错配；`add/sub/neg/mul/inv/div/eq_mod` 异域经 `align_pair` 内化 align（同域 ptr-eq/语义 eq 快路径）；`element_*` **保留不删**（增量非破坏）；ring-hom 单测（同域 √2·√2=2 / √2·inv=1 / (√2+1)(√2−1)=1；跨域 √2·√3→compositum dim 4，(√6)²=6）。 | 局部新增 | 否 — 2026-07 落地 |
| **P3b** | `element_*` 调用点迁移 | 把 ~285 处 `element_*(field, &a, &b)` 调用点逐文件迁移到 `a.mul(&b)`（异域自动 align）；`element_*` API 稳定性登记更新。 | 中改（机械迁移，分文件 PR） | 否 |
| ~~**P3c**~~ ✅ | `AlgExtData`/`AlgExtC` 双轨对齐到 `FieldElement` | `AlgExtData::add/sub/mul/eq_mod`+`add_aligned`/`mul_aligned` 经 `as_field_element()` 路由到 `FieldElement`（删私有 align_pair）；`AlgExtC` 跨域 align（`from_complex_parts`/`align_with`/`align_with_in_cache`）经 `FieldElement::align_into*`；`FieldElement` 新增 raw-align API；`poly_roots::embed_real_coords_for_parent` 跨域 align 迁移。re/im 同域算术保留 `field.element_*`（决策 C）。 | 中改 | 否 — 2026-07 落地 |

**建议次序：** P2a → P2b → P3。P2a/P2b 是存储层局部硬化，可独立 PR；P3 是全栈架构演进，依赖 P2a/P2b 的类型基建。

## 4. 代价与难点

1. **Arc 引用计数开销**：每个 `FieldElement` 带 `Arc<ExtensionField>`，海量元素（poly 系数、中间结果）有构造/销毁的引用计数增减。upstream `gen` 也是引用计数指针，同样开销，可接受。
2. **无循环引用**：`block.field = parent_field`，parent 链终止于 Base，是 DAG（向下指不回指 self）。✓（存储层 `Vec<FieldElement>` 安全）
3. **API 大改**（仅 P3）：`element_mul(field, &a, &b)` → `a.mul(&b)`，整个 algebra 层入口换型。主要工作量。
4. **域无关字面量保留**：`minpoly_x2_plus_c(-3)` 等 **adjoin 之前** 的 ℚ 系数字面量无域可绑（域是 adjoin 创建的），保持 `CoordsQ`。分层：`min_poly_q`（ℚ 字面量）裸；`min_poly_parent_blocks`（parent 域元素，adjoin 时 parent 已存在）可绑。
5. **align 性能**：运算时查"两边域是否同"（`Arc::ptr_eq` / `==`，便宜）；异域走 `common_cache`（已有），不每次重算。

## 5. 与现状的关系（分层对照）

| 层 | 现状 | FieldElement 模型 |
|---|---|---|
| ℚ 字面量（adjoin 前） | `CoordsQ`（裸，无域） | `CoordsQ`（不变，无域可绑） |
| 层 minpoly parent 系数 | `Vec<CoordsQ>` + 运行时校验 | `Vec<ParentBlock>` / `Vec<FieldElement>`（绑 parent，构造校验） |
| 运行时元素 | `ExtensionField` 方法 + `HighFirstQ`，调用方 align | `FieldElement` 方法，元素自身 align |
| 错配防线 | 使用点约定 + `build_adjoin_parent_coeffs` 校验 | 结构性消除（无裸坐标） |

## 6. 验收（按优先级）

- ~~**P2a**~~ ✅（2026-07）：`min_poly_parent_blocks` 类型升 `Vec<ParentBlock>`；持久存储唯一写入门 `build_adjoin_parent_coeffs` 经 `ParentBlock::new(parent, coords)`；构造时挡长度错配；`LayerMinPolyForAdjoin::ParentBlocks` 升 `&[ParentBlock]` 零分配借存储；瞬态 compositum blocks 保持 `&[CoordsQ]` 由消费方把守；全 workspace 1205 测试全绿；新增 `parent_block_new_rejects_wrong_length` + `build_adjoin_parent_coeffs_rejects_wrong_length_block` 单测；clippy 干净。
- ~~**P2b**~~ ✅（2026-07）：塔字段 `Option<Vec<ParentBlock>>` → `Option<MonicLayerMinPoly>`；`MonicLayerMinPoly::new(parent, blocks)` 集中全部不变量（degree≥1 + monic + per-block 长度经 `ParentBlock::new`）；`build_adjoin_parent_coeffs` 退一行（`let layer = MonicLayerMinPoly::new(parent, layer_blocks)?;`）；读取点（`semantic_key_bytes`/`layer_min_poly_exprs`/`layer_minpoly_parent_coeffs`/`flatten_min_poly_over_q_cold`/`layer_minpoly_coords_for_adjoin`）经 `.as_blocks()`；`LayerMinPolyForAdjoin::ParentBlocks` 仍 `&[ParentBlock]`；全 workspace 1207 测试全绿；新增 `monic_layer_min_poly_rejects_non_monic` + `monic_layer_min_poly_rejects_degree_zero` 单测；clippy 干净。
- **P2b**：`MonicLayerMinPoly::new` 挡 monic + 长度 + degree；`build_adjoin_parent_coeffs` 退一行；现有测试全绿。
- **P3**：`FieldElement` 作为运算单元；`a.mul(b)` 异域自动 align；无裸 `CoordsQ` 业务参数（L2/L3）；`element_*` API 稳定性登记更新；全 workspace + conformance 全绿。
  - ~~**P3a**~~ ✅（2026-07）：`FieldElement { field, coords }` 类型 + 核心运算（`add/sub/neg/mul/inv/div/eq_mod/is_zero/is_one`），异域 `align_pair` 内化（同域 ptr-eq/语义 eq 快路径），构造挡长度错配；`element_*` 保留不删（增量非破坏）；`field_element_new_rejects_wrong_length` + `field_element_same_field_ring_hom` + `field_element_cross_field_mul_aligns` 单测；nextest --workspace 1210 passed；clippy 干净。
  - **P3b**（in progress）：~285 处 `element_*` 调用点逐文件迁移到 `FieldElement` 方法；`element_*` API 稳定性登记更新。
    - ~~**P3b-1**~~ ✅（2026-07）：`FieldElement` session-aware 运算族（`add/sub/mul/div/eq_mod_with_session(&FieldSession)`）——经 `FieldSession::align_pair` 走 session `common_cache`（dedup），补齐 P3a 静态 `align_pair`（ephemeral cache 无 dedup）的缺口，为迁移 session-using 消费者铺路；`field_element_session_aware_mul_matches_static_and_dedups` 单测（结果等价 + 重复 align 命中 cache）；nextest --workspace 1211 passed；clippy 干净。
    - **P3b-2**（pending → 部分并入 P3c）：跨域 align 消费者点已随 P3c 迁移（`poly_roots::embed_real_coords_for_parent` 经 `FieldElement::align_into_with_session`）。
    - **P3b-3**（pending）：剩余 ~280 同域 `element_*` 调用点——经评估为同域无 align 收益换皮（每 op 增 Arc+coords clone），按决策跳过；如需绑定防错配再逐文件 idiomatic 提升。
  - ~~**P3c**~~ ✅（2026-07）：`AlgExtData`/`AlgExtC` 双轨对齐到 `FieldElement`。
    - `AlgExtData::add/sub/mul/eq_mod` + `add_aligned`/`mul_aligned`（alg_ext.rs）经 `as_field_element()` 路由到 `FieldElement`（静态 + session-aware），删私有 `align_pair`/`align_pair_with_session`。
    - `AlgExtC::from_complex_parts`/`align_with`/`align_with_in_cache`（alg_ext_c.rs）的跨域 align 经 `FieldElement::align_into`/`align_into_with_cache`；re/im 同域算术保留 `field.element_*`（按决策 C 跳过同域）。
    - `FieldElement` 新增 raw-align API：`align_into` / `align_into_with_session` / `align_into_with_cache`（暴露 aligned 坐标供 re/im split 等场景）。
    - `poly_roots::embed_real_coords_for_parent` 跨域 align 经 `FieldElement::align_into_with_session`。
    - nextest --workspace 1211 passed（AlgExtData/AlgExtC 算术被 poly_alg_*/conformance 重度覆盖，无回归）；clippy 干净。

## 7. 候选路径：AlgExt\* 薄包装融合（已采用，PR1/PR2 落地）

P3c 落地后 `AlgExtData`/`AlgExtC` 的算术已 delegate 到 `FieldElement`（桥接）。进一步"融合"= 让 `AlgExt*` 退成 `FieldElement` 的薄包装：存储从 `Vec<ExprArc>` 换成 `HighFirstQ`（`Vec<Ratio>`），Expr 互转 API 保留为视图。评估如下（2026-07，原评估未采用；后续 profile/架构复审后按 PR 切分推进，PR1/PR2 已落地——见 §7.6/§7.7）。

### 7.1 目标形态

- `AlgExtData = { field: Arc<ExtensionField>, coords: HighFirstQ, root_index: Option<u32> }`（`coords` 从 `Vec<ExprArc>` → `HighFirstQ`）≈ `FieldElement` + `root_index`。
- `AlgExtCData = { field, re: HighFirstQ, im: HighFirstQ, root_index }`（`re`/`im` 从 `Vec<ExprArc>` → `HighFirstQ`）= 两个 `FieldElement` + 共享 field。
- `AlgExtCPolyCoeff(pub AlgExtCData)` 薄 newtype，自动跟进。
- 算术已 delegate（P3c）；Expr 互转 API（`into_expr`/`from_rootof`/`from_field_coords`）保留，内部改 `coords_to_expr` 按需重建。

### 7.2 改动点统计（grep 实测）

- **`AlgExtData.coords` 直接访问 ~7 处**（`alg_ext.rs` 内 4：`is_zero`/`is_one`/`canonical_poly1_expr`/`coords_q`；`alg_ext_c.rs` 3：`rationalize_poly1(&a.coords)`）+ 构造 `from_coords_q`/`from_field_coords`/`into_expr`。
- **`AlgExtCData.re`/`.im` 直接访问 ~50 处跨 6 文件**（`poly_roots.rs` ~13、`alg_ext_c.rs` ~20 内部、`field_session.rs` ~13、`common_minimal.rs` 2、`poly_alg_factor.rs` 2、`eval.rs` 1、`alg_ext.rs` 2）。
- 模式集中为 3 种：(a) `rationalize_poly1(&x.re)` → 改 `x.re.0.clone()`（**省一次 rationalize，收益**）；(b) `x.im.iter().all(|c| c.is_zero())` → `x.im.as_slice().iter().all(...)`；(c) `x.re.clone()` 塞回 `Expr::Add`/`from_field_coords` → 改 `coords_to_expr(&x.re)` 重建。
- **`Expr::AlgExt`/`Expr::AlgExtC` 消费面**（eval/normal/显示读 `.coords`/`.re`/`.im` 渲染）：尚未 grep，需单独摸。

### 7.3 收益

1. 消除存储冗余：存 `Ratio`，省每次运算的 `rationalize_poly1`（`coords_q`/`re_q`/`im_q` 从"rationalize 存的 ExprArc"变成"clone 已是 Ratio 的 HighFirstQ"）。
2. `as_field_element()` 路径更短：从"rationalize + 构造 FieldElement"变成"clone HighFirstQ + 包 FieldElement"。
3. 融合度：`AlgExtData ≈ FieldElement + root_index`，双轨收窄为视图层。

### 7.4 代价

1. **`into_expr` 重建**：从 O(1) 包 `Arc`（复用存的 `ExprArc`）→ 每次 `coords_to_expr`（`Ratio→Expr::int/rat`，每系数一次 `Expr` 构造 + `Arc`）。eval/normal 把代数数塞回 Expr 树时调用，可能热——需 profile。
2. **`into_expr` fallback 语义变**：`AlgExtC` 的 `into_expr` 现有 `Err(_) => Expr::Add(self.re.clone())` fallback 保留原始符号 `ExprArc`；改存 `HighFirstQ` 后 fallback 要 `coords_to_expr` 重建，且重建出的是有理 `Expr`（丢失"原始符号形式"——但实测 `rationalize_poly1` 要求有理系数，原始本就只是有理包装，无真正符号损失）。
3. **~60 处机械改动跨 6+ 文件**，量大、易错。
4. **回归面大**：`AlgExtC` 算术被 `poly_alg_*`/conformance 重度覆盖，50 处改动任一处错都炸。

### 7.5 风险与决策

- `root_index` **非阻碍**：实测全 crate 无 `Some(_)` 赋值、无读取消费、无 consistency 校验——纯透传预留字段（全 `None`），融合加个透传字段即可。
- 真正阻碍已消除（coords 实际只是有理包装、root_index 无逻辑），**剩代价 = into_expr 重建 + 60 处机械改动 + 回归面**。
- **未采用理由**：收益（省 rationalize + 路径缩短）是性能边际 + 架构整洁度；代价是 60 处改动 + into_expr 重建潜在热路径回退 + 大回归面。P3c 桥接已让两层协同（Expr 层保 Expr 互转 + 计算层保绑定 align），融合的边际收益不抵代价。
- **采用条件**：若 profile 显示 `rationalize_poly1`（`coords_q`/`re_q`/`im_q`）成瓶颈，或 Expr round-trip 性能不敏感（conformance 实测无回退），再按 PR 切分推进：PR1 `AlgExtData`（~7 处，小）→ PR2 `AlgExtC`（~50 处，大，分文件）。

### 7.6 PR1 落地（AlgExtData.coords → HighFirstQ）

**已落地**（2026-07）。`AlgExtData.coords` 从 `Vec<ExprArc>` 改存 `HighFirstQ`（`Vec<Ratio>`，存时 `pad_to_len` 到 dim）。

- `alg_ext.rs`：结构体字段 `coords: HighFirstQ`；`from_coords_q` 存 `HighFirstQ::new(pad_to_len(...))`（省 `coords_to_expr`）；`coords_q` 返 `self.coords.0.clone()`（省 `rationalize_poly1`，O(1)）；`as_field_element` 用 `self.coords.clone()`（省 `coords_q` + `HighFirstQ::new` 包装）；`is_zero`/`is_one` 读 `Ratio`；`to_rootof_expr` 内联 trim leading zeros + `ratio_to_expr_arc` 逐系数重建（infallible，保持 `-> ExprArc` 签名）；移除 `coords_to_expr` import。
- `alg_ext_c.rs` 3 处外部 `rationalize_poly1(&a.coords)` → `a.coords.0.clone()`/`a.coords.clone()`（存时已 pad，长度与原 `rationalize_poly1` 保持一致，无行为变化）。
- `giac-calculus/expr_util.rs` `depends_on_var` 的 `Expr::AlgExt` arm 删 `a.coords.iter().any(...)` 子句（coords 现为 `Ratio`，绝不含变量；原 `rationalize` 后 `Int/Rat` 该检查亦恒 false，无行为变化）。
- 验证：nextest release --workspace 1211 passed + clippy 干净。
- 未简化：`mul_rational`/`inv`/`common_ext`/`algext_square_roots`/`algext_cube_root`/`fold_algext_sum` 内 `HighFirstQ::new(self.coords_q()?)` 模式保留（`coords_q` 现已 O(1)，等价开销，改动面 vs 收益不抵）。
- 待 PR2：`AlgExtCData.re`/`.im` `Vec<ExprArc>` → `HighFirstQ`（~50 处跨 6 文件）。

### 7.7 PR2 落地（AlgExtCData.re/.im → HighFirstQ）

**已落地**（2026-07）。`AlgExtCData.re`/`.im` 从 `Vec<ExprArc>` 改存 `HighFirstQ`（`Vec<Ratio>`，存时 `pad_to_len` 到 dim）。

- `alg_ext_c.rs`：结构体 `re: HighFirstQ`/`im: HighFirstQ`；`from_coords_q` 存 `HighFirstQ::new(pad_to_len(...))`（省 `coords_to_expr`）；`re_q`/`im_q` 返 `self.re.clone()`/`self.im.clone()`（省 `rationalize_poly1`，O(1)）；`coords_as_rational` 签名改 `&[Ratio<BigInt>]`（直接消费存的高理坐标，省逐 `ExprArc` rationalize）；`re_to_legacy_expr`/`re_to_legacy_arc`/`im_i_term` 用 `self.re.as_slice()` + `AlgExtData::from_field_coords_q` + `coords_to_expr(...).unwrap_or_default()`（fallback `Expr::Add` infallible）。
- 外部消费者（6 文件）：
  - `field_session.rs`：`rationalize_poly1(&inner.re)` → `inner.re.0.clone()`（5 处：`try_sqrt_in_field`/`try_sqrt_in_field_shallow`/`adjoin_sqrt`/`adjoin_cbrt`/`embed_coeff_into` + `is_negative_rational` let-else 去 rationalize）；`from_field_coords(.., inner.re.clone())` → `from_field_coords_q(.., inner.re.0.clone())`；`mul_formal_i` 字面构造 `re: coords_to_expr(&zero_coords())?` → `re: inner.field.zero_coords()`（HighFirstQ）。
  - `common_minimal.rs`：`rationalize_poly1(&inner.re)` → `inner.re.0.clone()`（`as_rational_coeff`）。
  - `poly_roots.rs`：`rationalize_poly1(&inner.re)` → `inner.re.0.clone()`（`embed_real_coords_for_parent`/`build_min_poly_parent_blocks`/`build_min_poly_parent_blocks_v2`）；`AlgExtCData` 字面构造 `im: coords_to_expr(&zero_coords())?` → `im: inner.field.zero_coords()`；diag test 内 `rationalize_poly1(&dq_inner.re).unwrap()` → `dq_inner.re.0.clone()`（4 处）。
  - `poly_alg_factor.rs`：`expr_to_ratio(alg_c.re.first()?.as_ref()).ok()?` → `alg_c.re.first().cloned()?`（re 现 `HighFirstQ`，首系数即 `Ratio`，省 `expr_to_ratio`）。
  - `giac-solve/realroot.rs`：`expr_is_real` 的 `AlgExtC` arm `z.im.iter().all(expr_is_real)` → `z.im.iter().all(|c| c.is_zero())`（im 现为 `Ratio`，原 `expr_is_zero` 只 match `Int(0)` 漏 `Rat(0)`，改后等价且更直接）。
  - `giac-calculus/expr_util.rs`：`depends_on_var` 的 `Expr::AlgExtC` arm `z.re/z.im.iter().any(depends_on_var)` → `Expr::AlgExtC(_) => false`（re/im 现为 `Ratio`，绝不含变量；原 `rationalize` 后 `Int/Rat` 该检查亦恒 false，无行为变化）。
- **底层 bug 修复（PR2 暴露 → 根因修复，非 test workaround）**：
  - **症状链**：`realroot(x⁴−1)` PR2 后返 4 项含 `rootof([-1],poly1[1,0,poly1[-1/2,0]])`（嵌套 `poly1`，conformance parser 不认）+ 有理 ±1，而非期望的 2 项 `[[-1,1],[1,1]]`。
  - **根因**（`ext_tower.rs` `try_square_root_in_field_impl`）：`try_sqrt_basis_squares` 只探 ±basis 元素（`e[i]=1`，cand `{[1],[-1]}`，查 `coords_square_eq_mod`），**不探有理倍数**。ℚ 上 `sqrt(4)=2·basis[0]` 非 basis 元素 → 返 `None`。biquadratic 对 `x⁴−1`（disc=4）走 `adjoin_sqrt_new(4)` → `adjoin x²−4`（**reducible，另一 bug**）生成 `Q(√4,u)` 域，根 ±1 表示为 `AlgExt`（`to_rootof_expr` → 嵌套 `poly1`）而非 `Int`。
  - **为何 ℚ 需特判**：`try_sqrt_subfield_descent` 遍历 `proper_subfields_chain(field)` 递归 `try_square_root_in_field_impl` 到子域。ℚ 是 Base，`proper_subfields_chain(ℚ)` 为空（`parent_field()==None`），无子域可 descent——ℚ 必须靠自身 base case。
  - **修复**（正确根因，非 hack）：`try_square_root_in_field_impl` 开头加 ℚ Base fast path：`rational_sqrt(u[0])`（`num_util::integer_sqrt` 检查 numer/denom 完全平方，负数/非完全平方返 `None`）→ 返 `[sqrt]`。这是 `try_sqrt_subfield_descent` 递归的 **base case**：任意非 Base 域的有理平方根（如 `Q(√2)` 中 `sqrt(4)`）经 descent 到 ℚ 走此特判，`emb.apply` 回原域——**递归机制自动覆盖全塔，无需每层重复**。ceil 注释已写明（只 ℚ；deeper towers fall through，实指 descent 递归到 ℚ）。
  - **非 workaround 的判据**：不是在 `realroot` 层 hack 简化 `rootof`→`Int`（那会掩盖根因），而是在 `ext_tower` sqrt 求解的数学正确层补全 ℚ Base 缺失的有理平方判定。修复后 biquadratic 在 ℚ 直接返 `Int ±1`，realroot 2 项，`rootof` 路径不再触发。
  - **P3c 隐藏机制**：P3c `sturmab_count_rational_poly(x⁴−1)` 返 `real_count=0`（x⁴−1 实根 ±1，应 = 2——另一 P3c bug），`real_algebraic_roots` 在 `real_count==0` 早返 `vec![]`，从未走到 biquadratic，隐藏了 ℚ sqrt 漏 + adjoin reducible 两个 bug。PR2 改 `AlgExtC` 算术后 sturmab 走通（`real_count=2`，根因未深挖——疑 P3c `rationalize_poly1` 路径健壮性，PR2 `re_q` infallible 修正之），暴露上述 bug 链。
  - **遗留 gap（登记 follow-up，非 PR2 scope）**：
    1. `adjoin` reducible minpoly 未拦截：`adjoin x²−4`（x²−4=(x−2)(x+2)）应失败却生成 `Q(√4)=Q(2)=ℚ` 的伪 dim-2 域。ℚ sqrt 特判绕开了此路径，但 reducible adjoin 本身是 bug，需在 adjoin 入口加 irreducibility 检查。
    2. `to_rootof_expr` 对嵌套塔输出嵌套 `poly1`（layer minpoly 系数是 parent 域 `Poly1` 元素）——不符合 giac `rootof` minpoly over-ℚ 语法。roundtrip 内部认（`try_as_algext_data`/`from_rootof`），但 `giac-parse` 外部 parse 不认。P3c 已有，ℚ sqrt 修复后 realroot 不再触发，但其他路径可能触发。应改用 `top_min_poly_exprs`（flatten over ℚ）——曾试但 `assert_equiv` 不识别 `rootof([c],flatten_minpoly)` 常数（`c` mod minpoly 求值未简化），回退；需配套 `assert_equiv`/eval 简化 `rootof` 常数。
    3. `sturmab_count_rational_poly(x⁴−1)` P3c 返 0 的根因（PR2 算术修正掩盖）——若 P3c `rationalize` 路径确有健壮性 bug，应独立定位而非依赖 PR2 算术掩盖。
- 移除 unused import：`rationalize_poly1`（`alg_ext_c.rs`/`common_minimal.rs`/`field_session.rs`/`poly_roots.rs` 顶层 + `poly_roots.rs` test 内 2 处）、`expr_to_ratio`（`poly_alg_factor.rs`）。
- 验证：nextest release --workspace 1211 passed + clippy 干净。
- fusion 完成：`AlgExtData.coords`/`AlgExtCData.re`/`.im` 全部存 `HighFirstQ`，算术 delegate `FieldElement`（P3c），Expr 互转 API 退成视图层（`coords_to_expr`/`ratio_to_expr_arc` 按需重建）。`root_index` 透传预留（全 `None`）。

### 7.8 遗留 GAP 修复（2026-07）

§7.7 三个遗留 gap 全部修复，nextest release --workspace 1216 passed + conformance 207 passed + clippy 干净。

1. **GAP-1 `adjoin` reducible minpoly 漏拦 → 已拦截**（`ext_tower.rs` `build_adjoin_irreducible`）。
   - 加 `ensure_minpoly_irreducible(parent, min_poly_q)`：`parent==ℚ` 时 `minpoly_coords_to_poly` → `factor_into`（ℚ 上完全），过滤 degree-0 content factor（`factor_into` 对 `x²−1/2` 返 `[1/2, 2x²−1]`，常数 content 非 reducible 证据），非常数 factor 数 >1 → `Err`。`degree<2` 放行（trivially irreducible）。
   - `parent` 非 Base **跳过**：`factor_into_algext` 当前误判 `u²+3 over Q(√2)` 为 reducible（factorer 不健全），若启用会破坏所有合法嵌塔 adjoin（实测 9+ test 失败）。ceil：非 Base reducible adjoin 仍为潜在风险，待 algext factorer 健全后启用；realroot 的 `x²−4 over ℚ` bug 由 ℚ 检查全覆盖。
   - 测试 `adjoin_irreducible_rejects_reducible_over_q`：`x²−4` 拒、`x²−2` 放。
2. **GAP-2 `to_rootof_expr` 嵌套 `poly1` → 用 ℚ-flatten minpoly**（`alg_ext.rs` `to_rootof_expr`）。
   - minpoly 从 `self.min_poly()`（layer，parent 域 `Poly1` 系数 → 嵌套）改 `self.field.top_min_poly_exprs()`（ℚ-flatten，纯 `Int/Rat` 系数，符合 giac `rootof` over-ℚ 语法）。
   - 单层域（parent==ℚ）flatten==layer，无变化；嵌塔输出 flatten（degree=dim，无嵌套 `poly1`）。
   - num 保持 layer coords（display 语义；嵌塔 roundtrip 无 test 覆盖，单层 roundtrip `algext_to_rootof_roundtrip_display` 不受影响）。
   - 测试 `to_rootof_expr_nested_tower_uses_flatten_minpoly`：嵌塔 `Q(√2,u²=3)` 的 rootof minpoly 系数全 `Int/Rat`（无嵌套 `poly1`）。
3. **GAP-3 `sturmab(x⁴−1)` P3c real_count=0 根因 → 确认 PR2 根因消除（非掩盖）**。
   - 根因：P3c `AlgExtCPolyCoeff` 算术走 `rationalize_poly1`（`Vec<ExprArc>`↔`CoordsQ` 往返），该路径在 Sturm 序列 `div_rem_wrt_algext`/`sign_of_algext_value` 求值中出错（具体 round-trip 步骤未逐一 bisect，因 PR2 已移除该路径）。
   - PR2 fusion（`re`/`im` 存 `HighFirstQ`）**移除算术的 rationalize 路径**（`re_q`/`im_q` O(1) 直返 `Ratio`），不是在 sturmab 层加特判绕过——根因路径不复存在，属根因消除而非掩盖。
   - 无独立代码修复（rationalize 算术路径已随 PR2 删除）；加 regression test `sturmab_count_x4_minus_one_two_real_roots`（`sturmab_count_rational_poly(x⁴−1,−10,10)==2`）锁住，防止 rationalize 算术路径回退。

## 8. 不在本 issue 范围

- 新 embedding 入口的 ring-hom 测试 + stability 登记（归档 issue §2 长期规则，随 PR 走）。
- 嵌套环类型清理（见 [GIAC-poly-nested-ring-types.md](GIAC-poly-nested-ring-types.md)）。
- flatten bisect fallback 的去除（见 [GIAC-algext-adoption.md](GIAC-algext-adoption.md) §7）。
