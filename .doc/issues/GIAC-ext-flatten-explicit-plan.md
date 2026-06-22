# GIAC — 去掉 `ExtensionDesc.min_poly_over_q: OnceLock`，显式 flatten API

**状态:** draft（R0–R5 ✅ 之后的收束项）  
**类型:** 实施计划 / AFK  
**父项:** [GIAC-ext-registry-removal-plan.md](GIAC-ext-registry-removal-plan.md)（Phase R6 建议名）  
**相关:** [GIAC-algext-adoption.md](GIAC-algext-adoption.md) §8.2；[giac-tower-common-math.md](../giac-tower-common-math.md)；[GIAC-poly-algext-backlog.md](GIAC-poly-algext-backlog.md) §5.1（deg≥5 solve 策略）  
**上游参考:** `giac/giac-2.0.0/src/alg_ext.cc` — `algebraic_EXTension`、`common_EXT`、`common_minimal_POLY`  
**Rust 落点:** `ext_tower.rs`、`field_session.rs`（可选 cache）、`poly_roots.rs`  
**快照:** 2026-06-22

---

## 1. 问题陈述

R0 用 `OnceLock<CoordsQ>` 缓存 **ℚ-flatten 本原元多项式**，避免 nested adjoin 在 register 时 `compose` / `char_poly`。R3 已将 **`AlgExtData::min_poly()`** 切到 **layer**（`layer_min_poly_exprs`），与 upstream `*(_EXTptr+1)` 对齐。

**仍存在的架构债：**

1. `ExtensionDesc` 上仍挂 **惰性 flatten 字段**（`OnceLock`），语义像「每个域句柄都可能有第二套 minpoly」。
2. **热路径误触发 flatten**：`compute_common_tower`（T4a）对 sibling 调 `min_poly_over_q()`；嵌套塔会无谓 `compose`。
3. `element_*_primitive` 经 `ensure_min_poly_over_q()`，对 adjoin over ℚ 本可读 `tower.min_poly_q`。
4. `PartialEq` / 测试曾围绕「register 不算 flatten」设计；`OnceLock` 仍让 flatten **隐式**、难审计。

**目标：** flatten 成为 **显式、可选、可定位** 的算法步骤；`ExtensionDesc` 只描述塔 + layer（对标 upstream adjoin），不缓存派生 flatten。

**非目标：**

- 删除 flatten 算法（`--no-default-features` bisect、`compute_common_flatten` 保留）
- 删除 `ExtensionField::rational()` 的 `static OnceLock`（MSRV 1.75 单例初始化；与 flatten 无关）
- 改 `Expr::AlgExt` AST

---

## 2. 目标架构（normative）

### 2.1 现状（R3 后）

```text
ExtensionDesc {
  tower: Arc<ExtensionTower>     // layer minpoly 真源
  min_poly_over_q: OnceLock<…>   // 惰性 flatten 缓存（派生）
  parent_field: Option<Arc<…>>
}

adjoin / 塔运算 ──► 不依赖 OnceLock
隐式 accessor ──► min_poly_over_q() / top_min_poly_exprs() 可能 compose
```

### 2.2 目标（R6）

```text
ExtensionDesc {
  tower: Arc<ExtensionTower>     // layer 真源（唯一 stored minpoly）
  parent_field: Option<Arc<…>>
  // 无 min_poly_over_q 字段
}

layer 运算 / display / 同域判定 ──► tower + layer_min_poly_exprs
flatten 需要时 ──► flatten_min_poly_over_q(&field, session?) -> CoordsQ   // 冷算或 session memo
可选 session ──► FieldSession::flatten_cache: HashMap<sem_key, CoordsQ>

adjoin over ℚ / flatten compositum 句柄：
  tower.Adj.min_poly_q 即为 ℚ 上 primitive poly（构造时已知，非 lazy）
```

### 2.3 与 upstream 对照

| 概念 | upstream `_EXT` | giac-rs 目标（R6） |
|------|-----------------|-------------------|
| 元素 minpoly | `*(_EXTptr+1)` layer `v` | `layer_min_poly_exprs()` / `tower` |
| adjoin 构造 | `algebraic_EXTension(a,v)` O(1) | `build_adjoin_*` O(1)，无 flatten |
| 同域 | `*(_EXTptr+1)` 相等 | `tower` / `semantic_key` / `Arc` 相等 |
| 跨域合并 | `common_EXT` → `common_minimal_POLY` | `common_over_q`（tower 或 flatten 路径） |
| ℚ-flatten poly | **不**挂在每个 `_EXT` 上 | **不**挂在 `ExtensionDesc` 上；显式函数或 session cache |
| 域标识 | `*(_EXTptr+2)` | `semantic_key` + `FieldSession` cache（Rust 无第三 gen 槽） |

**Rust 有意保留、upstream 无直接对应：**

- `ExtensionTower` 显式塔链 + `parent_field`（比 C++ 指针链更可审计）
- `FieldSession` per `Context` 的 `common_cache` / `adjoin_cache`（进程内 dedup，无全局 `Mutex`）
- `tower-common`（T4a）相对 flatten 的 compositum 快路径

---

## 3. 优势与代价

### 3.1 优势

| 优势 | 说明 |
|------|------|
| **语义清晰** | 真源 = layer；flatten = 派生算法，与 adoption §8.2 / upstream 心理模型一致 |
| **触发面可审计** | `rg flatten_min_poly_over_q` 即全量调用点；无隐式 `OnceLock::set` |
| **热路径更便宜** | T4a `compute_common_tower` 用 layer adjoin，嵌套塔不再误 compose |
| **不可变 `Arc<ExtensionDesc>`** | 去掉 per-field 内部可变性；`Eq` 无副作用顾虑更少 |
| **测试更简单** | R0 测改为「flatten 函数输入输出」；删除 `min_poly_is_cached` |
| **与 R3 一致** | `AlgExtData::min_poly` 已 layer；句柄上不再存第二套 minpoly |

### 3.2 代价 / 风险

| 代价 | 缓解 |
|------|------|
| 重复 flatten 计算 | 可选 `FieldSession.flatten_cache`（`semantic_key` 键） |
| `compose_min_poly_over_q` 递归曾依赖 parent `ensure_min_poly_over_q` | 改为递归调 `flatten_min_poly_over_q(parent)` |
| API 破坏：`min_poly_over_q()` 公开 accessor | 标 **deprecated** 一 PR → 删或改为 `flatten_min_poly_over_q` |
| flatten bisect 测 | 仍显式调 flatten；`#[ignore]` 慢测保留 |

---

## 4. FlattenCache 规范（PR-R6d，normative）

R6 去掉 `ExtensionDesc` 上的 `OnceLock` 后，**可选**在 `FieldSession` 内 memo flatten 结果。规范与现有 `common_cache` / `adjoin_cache` **同构**，避免第三套并发/生命周期模型。

### 4.1 存储形态

```rust
type FlattenCache = HashMap<Vec<u8>, CoordsQ>;

pub struct FieldSession {
    // … ambient, working, common_cache, adjoin_cache, base_by_min_poly …
    flatten_cache: Rc<RefCell<FlattenCache>>,
}
```

| 字段 | 类型 | 含义 |
|------|------|------|
| **键** | `Vec<u8>` | `field.semantic_key()`（tower + layer minpoly 字节；**禁止** `field.id()`） |
| **值** | `CoordsQ` | ℚ 上 flatten 后的本原元多项式（`poly1` 高次在前，`Ratio<BigInt>` 系数） |

**生命周期：**

| 事件 | 行为 |
|------|------|
| `FieldSession::new` | 空 `HashMap` |
| `Context::clone` / `fork_ambient` | `Rc::clone(&flatten_cache)` — 与 `common_cache` 共享 |
| 独立 `Context::new()` | 各自冷 cache |
| `#[cfg(test)] clear_caches_for_test` | `flatten_cache.clear()` |

**不存储：** 进程级全局表、`ExtensionDesc` 内字段、layer minpoly（真源在 `tower`）。

### 4.2 键规范

与 [GIAC-ext-registry-removal-plan.md](GIAC-ext-registry-removal-plan.md) §7 `semantic_key` 一致：

```text
cache_key = field.semantic_key()
          = ExtensionTower::semantic_key_bytes()   // parent 链 + layer kind + min_poly_key
```

| 可用作键？ | 原因 |
|-----------|------|
| `semantic_key` | **是**（normative）；同塔 duplicate `Arc`、session adjoin dedup 后必命中 |
| `Arc::ptr_eq` | 否 — 同塔不同 `id` 会重复冷算 |
| flatten 结果的 `min_poly_key_bytes` | 否 — 不同 tower 构造可数学等价；缓存语义应对 **域描述符** 去重 |
| `field.id()` | 否（R1 已废弃作数学/cache 键） |

**正确性：** 固定 `ExtensionDesc`（不可变 `tower`）下，flatten 结果由塔结构唯一确定；session 内无失效/驱逐。

### 4.3 API 与查找顺序

**冷算（无 memo）：**

```rust
// ext_tower.rs — 迁自 compute_lazy_min_poly_over_q
pub(crate) fn flatten_min_poly_over_q_cold(
    field: &ExtensionField,
) -> Result<CoordsQ, EvalError>;
```

**统一入口（带可选 session）：**

```rust
pub(crate) fn flatten_min_poly_over_q(
    field: &ExtensionField,
    session: Option<&FieldSession>,
) -> Result<CoordsQ, EvalError>;
```

**查找顺序（必须按此实现）：**

```text
1. 快路径（不写 flatten_cache）
   - Base → 空 / 零维约定
   - is_simple_over_q() 或 parent.is_base()
     → clone tower.Adj.min_poly_q（layer = flatten，O(1)）

2. session = Some(s)
   - key = field.semantic_key()
   - if let Some(hit) = s.flatten_cache.borrow().get(&key) → return hit.clone()
   - poly = flatten_min_poly_over_q_cold(field)?
   - s.flatten_cache.borrow_mut().insert(key, poly.clone())
   - return poly

3. session = None（静态 API / 无 Context 单测）
   - flatten_min_poly_over_q_cold(field)   // 与 R5b ephemeral common 同级
```

**`FieldSession` 薄封装（可选，便于 roots 热路径）：**

```rust
impl FieldSession {
    pub(crate) fn flatten_min_poly_over_q(
        &self,
        field: &ExtensionField,
    ) -> Result<CoordsQ, EvalError> {
        flatten_min_poly_over_q(field, Some(self))
    }
}
```

### 4.4 与 layer / primitive 的分工

| 需求 | 来源 | 是否用 FlattenCache |
|------|------|-------------------|
| `AlgExtData::min_poly()` / `to_rootof` | `layer_min_poly_exprs()` | **否** |
| `element_*_primitive`（adjoin over ℚ） | `tower.min_poly_q` | **否** |
| `element_*_tower`（嵌套塔） | parent 系数环 + layer blocks | **否** |
| `compute_common_tower`（T4a） | sibling **layer** minpoly adjoin | **否** |
| `compute_common_flatten` | `flatten_min_poly_over_q(·, session?)` | **是**（嵌套域） |
| `poly_roots::approx_real_sign` | flatten 数值嵌入 | **是**（建议 `ctx.session()`） |
| `top_min_poly_exprs` | flatten 元数据 | **是**（显式调用方） |

### 4.5 数据流

```text
                    ┌──────────────────────────────────┐
                    │ FieldSession                     │
                    │  flatten_cache: RefCell<HashMap>   │
                    │    sem_key(√2)      → [1,0,-2]   │
                    │    sem_key(K₁(u²−3)) → degree-4  │
                    └───────────────┬──────────────────┘
                                    │
poly_roots / flatten common ────────┤ flatten_min_poly_over_q(f, Some(session))
                                    │
              miss ─────────────────┼──► flatten_min_poly_over_q_cold(f)
                                    │         compose / char_poly_matrix
                                    └──► insert(key, CoordsQ); return
```

### 4.6 与 upstream 的差异（刻意）

| | upstream | giac-rs FlattenCache |
|--|----------|----------------------|
| flatten 存放位置 | 不算在 `_EXT` 上；`common_minimal_POLY` 现场算 | session 内 `HashMap` memo |
| 作用域 | 单次 `common_EXT` 调用 |  per `Context` / `fork_ambient` 共享 |
| 数学 | 每次重算 | 同 `semantic_key` 结果相同；等价 |

### 4.7 实现注意

1. **值必须 owned clone**：`borrow` 结束后仍可返回 `CoordsQ`。
2. **`fork_ambient` 必须** `Rc::clone(&self.flatten_cache)`，与 `common_cache` 一致。
3. **R6a/R6b 可不引入** FlattenCache；R6d 再加字段，避免与删 `OnceLock` 混在同一 PR。
4. **体积**：每域一条 `CoordsQ`，度 ≤ 各层次数之积；单次 eval 触达域数量有限，可接受。
5. **递归 compose**：`flatten_min_poly_over_q_cold` 内对 parent 调 `flatten_min_poly_over_q(parent, session)`（同一 session 可命中 parent 项）。

### 4.8 测试

| 测项 | 断言 |
|------|------|
| `flatten_cache_hit_same_semantic_key` | 同塔 duplicate `Arc`，第二次不增 compose 计数（若保留计数器） |
| `context_clone_shares_flatten_cache` | `ctx.clone()` 后 flatten 命中 |
| `context_new_independent_flatten_cache` | 两 `Context::new()` 不共享 |
| `simple_over_q_skips_cache` | ℚ(√2) 走快路径，`flatten_cache_len` 不变 |
| `flatten_compose_fail_no_cache_insert` | `compose` 返回 `NotImplemented` 时 cache 无脏条目 |
| `nested_simple_over_q_fast_path` | 嵌套塔 `!is_simple_over_q()` 不走 §4.3 快路径 |

### 4.9 `compose` 失败、体积与 k 搜索上限

| 情形 | 规范 |
|------|------|
| `compose_min_poly_over_q` 在 `k=1..12` 内无本原元 | 返回 `EvalError::NotImplemented`；**禁止**向 `flatten_cache` insert |
| 嵌套塔冷算 | 度可达各层次数之积；单次 eval 域数量有限，**无驱逐**（ponytail：后期可加 LRU） |
| 高维 parent 递归 compose | 指数级搜索成本属 flatten 算法本身；R6 后每次冷算付全价，**R6d** 用 `semantic_key` 去重 |
| `is_simple_over_q()` 快路径 | **仅** adjoin over ℚ 或 flatten compositum 句柄；**嵌套** parent≠Base 必须 `flatten_min_poly_over_q_cold` |

### 4.10 T4a 资格与 flatten fallback

`tower_common_eligible` 要求两边 `is_simple_over_q()` 且非子域关系（见 `ext_tower.rs`）。

| 场景 | common 路径 | R6 后 |
|------|-------------|-------|
| ℚ(√2) × ℚ(√3) | T4a layer adjoin | R6a 用 layer，不 flatten |
| 两不同五次不可约 adjoin（并列因子） | **不满足** T4a → `compute_common_flatten` | 仍走显式 `flatten_min_poly_over_q`；慢但可审计 |
| T3+ sibling（parent 系数层 minpoly） | T4a 用 **layer blocks**，非 flatten | 见 §5.1 `layer_minpoly_coords_for_adjoin`（**P0**） |

---

## 5. 分阶段实施（建议 PR 顺序）

### 5.0 修改优先级总表

| 优先级 | 含义 | 落点 PR |
|--------|------|---------|
| **P0** | 阻塞热路径 / T4a 正确性；R6a 必须 | PR-R6a |
| **P1** | 删 `OnceLock` 前 API 与调用点收束 | PR-R6b、PR-R6c |
| **P2** | 性能 memo；可与 R6c 并行 | PR-R6d（可选） |
| **P3** | 与 R6 解耦；记入 plan，**不**挡 R6 合并 | 后续 P3-6 / P3-7 / 数值 |

| ID | 任务 | 优先级 | PR |
|----|------|--------|-----|
| M1 | `layer_minpoly_coords_for_adjoin`（rational + T3+ blocks） | **P0** | R6a |
| M2 | `compute_common_tower` 改 layer adjoin（含 T3+） | **P0** | R6a |
| M3 | `element_*_primitive` 读 `tower.min_poly_q` | **P0** | R6a |
| M4 | 统一 `flatten_min_poly_over_q(f, session?)` 签名（§4.3） | **P1** | R6b |
| M5 | `compose` / `compute_common_flatten` / `top_min_poly_exprs` 迁移 | **P1** | R6b |
| M6 | `approx_real_sign` → `flatten_min_poly_over_q(..., session)` | **P1** | R6b |
| M7 | `min_poly_over_q()` deprecated → 删 | **P1** | R6b→R6c |
| M8 | `duplicate_field_arc_for_test` 去掉 flatten 字段复制 | **P1** | R6c |
| M9 | R0 测试改写 + api-stability | **P1** | R6c |
| M10 | `FieldSession.flatten_cache` + §4.8 测 | **P2** | R6d |
| M11 | T3+ `compute_common_tower` 回归测 | **P1** | R6a 或 R6b |
| M12 | registry 父计划链到本 issue | **P1** | 文档 |
| M13 | deg≥5 数值 / `all_real_roots_minpoly` 扩展 | **P3** | 非 R6；见 §10 |
| M14 | `alg_ext_c::align_with` session 对齐 | **P3** | 独立 issue |

### 5.1 `layer_minpoly_coords_for_adjoin`（normative，P0）

R6a **禁止**在 T4a 对 sibling 调 `min_poly_over_q()` / flatten。统一 helper：

```rust
// ext_tower.rs — pipeline private
enum LayerMinPolyForAdjoin<'a> {
    Rational(&'a CoordsQ),           // tower.min_poly_q（T1/T2/flatten compositum）
    ParentBlocks(&'a [CoordsQ]),    // T3+ min_poly_parent_blocks
}

fn layer_minpoly_coords_for_adjoin(
    field: &ExtensionField,
) -> LayerMinPolyForAdjoin<'_>;
```

| `tower` 形态 | 返回 | 用于 |
|--------------|------|------|
| `AdjOverQ` / simple over ℚ | `Rational(tower.min_poly_q)` | `compute_common_tower` sibling adjoin |
| `AdjParentCoeffs`（T3+） | `ParentBlocks(...)` | 在 parent 上 `adjoin_irreducible_parent_coeffs` |
| Base | 错误 / 不可 adjoin | — |

**验收：** `adjoin(K₁, u²−α)` common 路径 register 阶段无 `compose`；与 `min_poly_is_layer_not_flatten_over_nested_adjoin` 一致。

### PR-R6a — 修热路径 + primitive 不碰 OnceLock（**P0**，小 diff，先合）

| 任务 | 优先级 | 文件 |
|------|--------|------|
| 实现 `layer_minpoly_coords_for_adjoin`（§5.1） | P0 | `ext_tower.rs` |
| `compute_common_tower`：sibling → **layer**（含 T3+ `ParentBlocks`） | P0 | `ext_tower.rs` |
| `element_*_primitive`：`ensure_min_poly_over_q` → `tower.min_poly_q`（仅 `!uses_tower_arithmetic()`） | P0 | `ext_tower.rs` |
| 新增 T3+ common 回归测（M11） | P1 | `ext_tower.rs` tests |
| 回归：`common_sqrt2_sqrt3_*`、`t4a::*` 全绿 | P0 | tests |

**验收：** T4a common(√2,√3) 不调用 `compose_min_poly_over_q`（`#[cfg(test)]` 计数器或 mock 可选）。

### PR-R6b — 显式 flatten API，窄化 accessor（**P1**）

| 任务 | 优先级 | 文件 |
|------|--------|------|
| `flatten_min_poly_over_q_cold` + `flatten_min_poly_over_q(f, session?)`（**签名以 §4.3 为准**） | P1 | `ext_tower.rs` |
| `compute_common_flatten`、`top_min_poly_exprs` 改调 flatten（`session` 有则传） | P1 | `ext_tower.rs` |
| `poly_roots::approx_real_sign` → `flatten_min_poly_over_q(&field, ctx.session())` | P1 | `poly_roots.rs` |
| `compose_min_poly_over_q` 内 parent：`flatten_min_poly_over_q(parent, session)?` | P1 | `ext_tower.rs` |
| `min_poly_over_q()` 标 `#[deprecated]`，thin wrapper 调 `flatten_min_poly_over_q(f, None)` | P1 | `ext_tower.rs` |
| `tower_common_matches_flatten_*` 两边显式 flatten | P1 | `ext_tower.rs` tests |

**验收：** `cargo test -p giac-core`；`--no-default-features` flatten smoke 仍绿。

### PR-R6c — 删除 `OnceLock` 字段（**P1**）

| 任务 | 优先级 | 文件 |
|------|--------|------|
| 删 `min_poly_over_q: OnceLock`、`min_poly_lock`、`ensure_min_poly_over_q`、`min_poly_is_cached` | P1 | `ext_tower.rs` |
| `build_adjoin_*` / `build_base_extension_uncached`：仅写 `tower.min_poly_q` | P1 | `ext_tower.rs` |
| `duplicate_field_arc_for_test`：只复制 `tower` + `parent_field`（**不**复制 flatten） | P1 | `ext_tower.rs` |
| R0 测试改写：`flatten_min_poly_over_q` 度 4；register 后无 flatten 字段 | P1 | `ext_tower.rs` tests |
| 删 `min_poly_over_q()` 公开 accessor；更新 [giac-core-algebra-api-stability.md](../giac-core-algebra-api-stability.md) | P1 | `.doc` |
| [GIAC-ext-registry-removal-plan.md](GIAC-ext-registry-removal-plan.md) 增 Phase **R6** 行链到本文 | P1 | `.doc` |

**保留：** `ExtensionField::rational()` 的 `static OnceLock<Arc<…>>`（ℚ 单例）。

### PR-R6d（可选）— Session flatten cache（**P2**）

见 **§4 FlattenCache 规范**。任务摘要：

| 任务 | 优先级 | 文件 |
|------|--------|------|
| `flatten_cache: Rc<RefCell<FlattenCache>>`；`fork_ambient` / `clear_caches_for_test` 接线 | P2 | `field_session.rs` |
| R6b 已落地的 `flatten_min_poly_over_q(f, session)` 接入 cache 查找（§4.3 顺序） | P2 | `ext_tower.rs` |
| `poly_algext_roots_for_ctx` 等改 `ctx.session().flatten_min_poly_over_q(…)` | P2 | `poly_roots.rs` |

**验收：** §4.8 测项 + 同一 session 内同 `semantic_key` 不重复 compose。

---

## 6. 调用点迁移清单

| 当前 | R6 后 | 优先级 |
|------|--------|--------|
| `element_*_primitive` → `ensure_min_poly_over_q` | `tower.min_poly_q`（adjoin over ℚ / flatten 域） | P0 |
| `compute_common_tower` → `sibling.min_poly_over_q()` | `layer_minpoly_coords_for_adjoin(sibling)` | P0 |
| `compute_common_flatten` → `a.min_poly_over_q()` | `flatten_min_poly_over_q(a, session?)` | P1 |
| `compose_min_poly_over_q` → parent `ensure_min_poly_over_q` | `flatten_min_poly_over_q(parent, session?)` | P1 |
| `top_min_poly_exprs` | `coords_to_expr(&flatten_min_poly_over_q(self, None)?)` | P1 |
| `poly_roots::approx_real_sign` | `flatten_min_poly_over_q(&field, ctx.session())` | P1 |
| `AlgExtData::min_poly` | 不变（`layer_min_poly_exprs`） | — |
| `tower_common_matches_flatten_*` | 两边显式 `flatten_min_poly_over_q(..., None)` | P1 |
| `duplicate_field_arc_for_test` | 不复制 flatten；仅 `tower` + `parent_field` | P1 |
| `min_poly_over_q()` 公开 API | deprecated → 删除；索引见 api-stability | P1 |
| `alg_ext_c::align_with`（静态 ephemeral common） | **不在 R6**；独立 session 对齐（M14） | P3 |

---

## 7. 测试策略

1. **PR-R6a 后：** `cargo test -p giac-core`；T4a 无回归；**新增** T3+ sibling common 不测 flatten（M11）。
2. **PR-R6b 后：** `compose` 递归经 `flatten_min_poly_over_q`；`approx_real_sign` 仍仅 deg≤3 数值（见 §10.3）。
3. **PR-R6c 后：** 删 `min_poly_is_cached`；`flatten_min_poly_over_q_nested_degree_four`；`duplicate_field_arc` 无 flatten 字段。
4. **PR-R6d 后：** §4.8 全表；`compose` 失败不污染 cache。
5. **慢测：** `tower_common_matches_flatten_*` 仍 `#[ignore]`；flatten bisect `--no-default-features`。
6. **合并门禁：** `cargo test-timeout` + `cargo ci-clippy`。

---

## 8. DoD

- [x] `ExtensionDesc` 无 `min_poly_over_q` / `OnceLock`（flatten 字段）
- [x] `flatten_min_poly_over_q` 为唯一 flatten 入口（+ 可选 session 版）
- [x] T4a `compute_common_tower` 不触发 flatten
- [x] `AlgExtData::min_poly` 仍为 layer（R3 不退化）
- [x] `ExtensionField::rational()` 单例仍可用（允许保留一处 `static OnceLock`）
- [x] FlattenCache（R6d）符合 §4 键/生命周期/快路径规范
- [x] `layer_minpoly_coords_for_adjoin` 覆盖 T3+；T4a 不测 flatten
- [x] api-stability 文档更新；registry 计划 Phase R6 链到本 issue
- [ ] `cargo test-timeout` + `cargo ci-clippy` 全绿（合并门禁）

---

## 9. 架构 review 摘要（vs upstream）

**对齐：** adjoin 只存 layer；flatten 仅在 common/roots 等「贵路径」显式出现 — 与 `algebraic_EXTension` + `common_EXT` 分工一致。

**Rust 增强（保留）：** 显式塔类型、`FieldSession` 本地化 cache、语义键 — 弥补 upstream 无全局 field 表但 `_EXTptr+2` 等隐式约定难审计的问题。

**刻意不同：** flatten 结果可选 session 缓存（upstream 每次 `common_minimal_POLY` 重算）；数学等价，性能可优于冷路径重复 compose。

**不建议：** 为删光文件内 `OnceLock` 而引入 `once_cell` 或升 MSRV 仅换 `LazyLock` — ℚ 单例一处 `OnceLock` 成本可接受。

---

## 10. 与 deg≥5 / solve 管线（P3-7）的关系

**normative 引用：** [GIAC-poly-algext-backlog.md](GIAC-poly-algext-backlog.md) §5.1 — 通用五次根式**不做**；不可约 deg≥5 → **`rootof(α, Pᵢ)`** + factor 降次；数值互补 `fsolve`，不替代精确 `rootof`。

### 10.1 R6 对五次主路径的影响

| 方面 | 影响 | 优先级 |
|------|------|--------|
| 不可约 P₅ → `rootof` + `adjoin_irreducible(ℚ, P₅)` | **正面**：用 **layer** minpoly，**不**需 flatten | — |
| `poly_algext_roots` deg=5 闭式 | **无影响**：仍 `NotImplemented`；属 P3-7/P4-6，非 R6 | P3 |
| factor 降次后递归 roots | **轻微正面**：R6a 减误 flatten，session 更干净 | P3 |
| 两不同五次因子进公共域 | T4a **不适用** → 仍 `compute_common_flatten`；R6 只使 flatten **显式** | P3 |
| 可解五次 / 分裂域特判 | 应优先 **tower common + layer adjoin**；R6 是前提，**不**实现 Galois | P3 |

**结论：** R6 是 **ext_tower 基础设施**；五次精确求根依赖 **P3-7 factor + `rootof` +（可选）T3+ adjoin**，**不**应把 flatten 后的 ℚ-本原元多项式当作 deg≥5 的工作模数。

### 10.2 与四次基建（P3-6）的衔接

一般四次在 K 上 resolvent 依赖 **T3+ layer adjoin**（backlog §5.1「与塔计划阻塞」）。R6a **P0** 的 `layer_minpoly_coords_for_adjoin` 直接服务 P3-6，避免 common 时误 flatten sibling。

### 10.3 数值分支上限（**P3**，非 R6 范围）

`poly_roots::approx_real_sign` 在 R6b 后仍经 `flatten_min_poly_over_q`，但：

- `all_real_roots_minpoly` 仅解 **deg 1–3**；deg≥4 返回空；
- `dim > 6` 时 `approx_real_sign` 直接 `None`。

故 R6 **不**开通五次欧拉式数值启发；若后期需要实根隔离，应走 **sturm / fsolve** 或扩展数值内核（M13），与删 `OnceLock` 正交。

### 10.4 FlattenCache 与五次

纯 `rootof` 路径通常 **不**重复 flatten；R6d cache 主要惠及 `compute_common_flatten` 与 roots 数值嵌入。dim=5 域缓存一条 degree-5 `CoordsQ` 合理；嵌套 flatten 度可能 >5，与五次求根策略无关。
