# giac-core / algebra API 稳定性分层

**规范来源:** [algorithm-expr-api.md](algorithm-expr-api.md)  
**Expr ↔ Poly 边界:** [expr-poly-conversion.md](expr-poly-conversion.md)（normative）  
**扩域求根:** [issues/GIAC-poly-roots-field-session-plan.md](issues/GIAC-poly-roots-field-session-plan.md)  
**Registry 移除 / session 接线:** [issues/GIAC-ext-registry-removal-plan.md](issues/GIAC-ext-registry-removal-plan.md)（R4–R5）  
**代码:** `giac-rs/crates/giac-core/src/algebra/*`

---

## 1. 源码注释格式

| 标记 | 用于 | 可见性 |
|------|------|--------|
| `/// **Stable** — …` | 有 I/O 契约、可跨 crate 调用 | `pub` |
| `/// **Stable (bounded)** — …` | 稳定但次数/域/算法阶段受限 | `pub` |
| `/// **Partial** — …` | 窄路径；注释写扩展或退役 | `pub` |
| `// **Pipeline private** — …` | 模块内辅助 | `fn` 私有 |
| `// **Temporary** — …` | shim / 待 FieldSession 吸收 | `fn` 私有 |

**未标注的私有函数:** 运行 `giac-rs/scripts/annotate_api_tiers.py` 补 `// **Pipeline private**`；Per-file 全表见本文末尾。

**提交前复审:** [algorithm-expr-api.md §7.2](algorithm-expr-api.md#72-测试通过后提交--合入前复审)

---

## 2. Crate 公开 API（`giac-core` re-export）

### 2.1 Expr ↔ Poly（path A / B）

| 函数 | 层级 | 模块 | 边界 |
|------|------|------|------|
| `expr_to_poly` | **Stable** | `poly` | 系数 ∈ ℚ；含 `AlgExt`/`rootof` → `TypeError` |
| `poly_alg_from_expr` | **Stable** | `poly` | 系数 ∈ K / `AlgExtC` |
| `poly_to_expr` | **Stable** | `poly` | ℚ 系数保持有理 |
| `algext_poly_to_expr` | **Stable** | `poly` | `PolyAlgExt` → `Expr` |
| `univariate_poly_to_poly1_expr` | **Stable** | `poly` | 一元 `poly1[…]` 高次在前 |
| `expr_contains_alg_coeff` | **Stable** | `poly_conv` | 调用方选 path A vs B |

### 2.2 扩域运算

| 函数 | 层级 | 模块 | 边界 |
|------|------|------|------|
| `AlgExtData` 算术 | **Stable** | `alg_ext` | 同域 ±×、`inv`、`eq_mod` |
| `fold_algext_sum` / `fold_algext_product` | **Stable** | `alg_ext` | 规范合并 |
| `fold_algext_sum_for_ctx` / `fold_algext_product_for_ctx` | **Stable (bounded)** | `alg_ext` | 经 `Context::session()` 的 `common_cache`（R5 eval） |
| `common_ext` | **Stable** | `alg_ext` | 两元公共扩域 |
| `algext_square_roots` / `algext_cube_root` | **Stable (bounded)** | `alg_ext` | 扩域内开方 |
| `canonicalize_to_algext_c` | **Stable** | `alg_ext_c` | Expr → `AlgExtCData` |
| `ExtensionField` / `ExtensionTower` / `ExtensionDesc` | **Stable** | `ext_tower` | R3：`ExtensionField` = `ExtensionDesc`；`layer_min_poly_exprs` / `top_min_poly_exprs` |
| `AlgExtData::min_poly` | **Stable** | `alg_ext` | R3：塔顶 **layer** minpoly（非 ℚ-flatten） |
| `poly_algext_roots` | **Stable (bounded)** | `poly_roots` | deg 1–4；独立 session（测试 / 无 Context 调用方） |
| `poly_algext_roots_for_ctx` | **Stable (bounded)** | `poly_roots` | 同上；`fork_ambient` 共享 `ctx.session()` extension cache（R5 solve） |
| `FieldSession` | **Stable (bounded)** | `field_session` | ambient K + working L；extension cache 可 `fork_ambient` 共享 |
| `Context::session` / `session_mut` | **Stable (bounded)** | `context` | 每 `Context` 一份 `Rc<FieldSession>`；`Clone` 浅拷贝共享 cache（`!Send`） |

### 2.3 Extension cache 与线程（R4–R5）

| 存储 | 位置 | 生命周期 | 键 |
|------|------|----------|-----|
| `common_cache` | `FieldSession`（`RefCell<…>`，可 `fork_ambient` 共享） | per `Context`；`Clone` 共享 | 域对 `semantic_key` 字典序 |
| `adjoin_cache` / `base_by_min_poly` | 同上 | 同上 | adjoin 语义键 / ℚ minpoly 字节 |
| 进程级 `FieldRegistry` | — | **已删除**（R4） | — |

静态 `ExtensionField::{adjoin_irreducible, common_over_q, align_elements}`：**无 session 时**每调用 ephemeral cache（R5b，无跨调用 dedup）。热路径用 `FieldSession` / `Context::session()` 做 adjoin 与 common dedup。

---

## 3. I/O 契约（§3 模板）

### `expr_to_poly`

| 字段 | 说明 |
|------|------|
| **输入** | `&Expr`；多项式形状；系数 ∈ ℚ |
| **输出** | `Ok(Poly)` 或 `TypeError`（含代数系数） |
| **禁止** | 含 `rootof` / `AlgExt` 静默落入 |

### `poly_alg_from_expr`

| 字段 | 说明 |
|------|------|
| **输入** | `&Expr`；变元 + `AlgExtC` 系数 |
| **输出** | `Ok(PolyAlgExt)` |
| **上下文** | 入口 normalize 到 ambient K（`poly_roots` 内） |

### `poly_algext_roots`

| 字段 | 说明 |
|------|------|
| **输入** | `&PolyAlgExt`, `&Var` |
| **输出** | `Vec<AlgExtCPolyCoeff>` |
| **上下文** | 内部 `FieldSession::new(infer_field(p))`；与 `Context` 无关 |
| **边界** | deg≤4；纯三次 `t³+a` 三根（ω 分支）；resolvent 公式见 `resolvent_golden_*`；一般四次 e2e 见 `#[ignore]` |

### `poly_algext_roots_for_ctx`

| 字段 | 说明 |
|------|------|
| **输入** | `&PolyAlgExt`, `&Var`, `&Context` |
| **输出** | 同 `poly_algext_roots` |
| **上下文** | `ctx.session().fork_ambient(K)`：独立 ambient/working，**共享** extension cache（`common` / adjoin dedup） |
| **调用方** | `giac-solve::eval_solve` → `try_algext_or_rootof_roots` |
| **边界** | 同 `poly_algext_roots` |

### `Context::session` / `session_mut`

| 字段 | 说明 |
|------|------|
| **输入** | `&Context` |
| **输出** | `Rc<FieldSession>` |
| **语义** | 该 `Context` 拥有的 session；`Context::new()` 各自独立 cache |
| **`Clone`** | 浅拷贝 `Rc` → 兄弟 handle 共享 cache |
| **线程** | `Context: !Send`（`Rc<FieldSession>`） |

### `fold_algext_sum_for_ctx` / `fold_algext_product_for_ctx`

| 字段 | 说明 |
|------|------|
| **输入** | `&[ExprArc]`, `&Context` |
| **输出** | 合并后的 `ExprArc` |
| **上下文** | 跨域 `align` 走 `ctx.session().common_cache`（`eval` Add/Mul 路径） |
| **对比** | 无 Context 的 `fold_algext_sum` 用 ephemeral `align_elements`（无 session cache） |

### `ExtensionField::common_over_q` / `align_elements`

| 字段 | 说明 |
|------|------|
| **输入** | 两个 `Arc<ExtensionField>`（+ coords for align） |
| **输出** | `CommonFieldPair` / `AlignedElements` |
| **缓存** | 无跨调用 cache；`FieldSession::common_over_q` / `common_over_q_with_cache` |
| **显式 session** | `common_over_q_with_cache` / `align_elements_with_cache`；`FieldSession::common_over_q` |

### `expr_contains_alg_coeff`

| 字段 | 说明 |
|------|------|
| **输入** | `&Expr` |
| **输出** | `bool` |
| **用途** | 调用方在 `expr_to_poly` 前分支 |

---

## 4. Pipeline private（代表性）

| 模块 | 函数 | 说明 |
|------|------|------|
| `poly_roots` | `infer_field`, `normalize_coeffs`, `pure_cubic_roots`, `build_resolvent_cubic` | 算术在 `FieldSession`；入口 normalize 到 K |
| `field_session` | `fork_ambient`, `adjoin_irreducible*`, `adjoin_sqrt` | extension cache + working L（R4–R5） |
| `field_arith` | `poly_*_with_coeffs_in_field` | 坐标多项式算术 |
| `ext_tower` | `align_elements_in_cache`, `common_over_q_in_cache`, `build_adjoin_*` | R4：无 `FieldRegistry`；dedup 在 session |

## 5. 跨 crate 契约

| 调用方 | Stable API | 注意 |
|--------|------------|------|
| `giac-simplify` | `expr_to_poly`, `poly_to_expr` | 超越叶子失败时原样返回 |
| `giac-solve` | `expr_to_poly`, `poly_algext_roots_for_ctx`, `quadratic_rootof_roots` | solve 经 `Context` session cache；一般四次 → `poly_algext_roots_for_ctx` 或 NotImplemented |
| `giac-calculus` | `expr_to_poly`（partfrac） | 仅 ℚ 有理分式 |
| `giac-core::eval` | `fold_algext_*_for_ctx` | Add/Mul 中 AlgExt 合并用 `ctx.session()` |
| `giac-poly` | — | 本模块提供 Expr 边界，算法在 `giac-poly` |

---

## 6. 维护

```bash
cd giac-rs
python3 scripts/annotate_api_tiers.py              # 补 tier 注释（幂等）
python3 scripts/annotate_api_tiers.py --inventory   # 刷新本文 Per-file 表
```

新增 `pub` 函数时：源码 tier 注释 → 更新 §2–§3 → 跑 inventory。

---

## Per-file function inventory (generated)

Regenerate: `python3 scripts/annotate_api_tiers.py --inventory`

### `alg_ext.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `min_poly` | **Stable** | `Poly::min_poly` |
| `zero` | **Stable** | Poly zero |
| `one` | **Stable** | Poly one |
| `from_field_coords` | **Stable** | `from_field_coords` |
| `from_coords_q` | **Pipeline private** | `from_coords_q` |
| `into_expr` | **Stable** | `into_expr` |
| `from_rootof` | **Stable** | `from_rootof` |
| `from_rootof_over` | **Stable** | `from_rootof_over` |
| `add` | **Stable** | Poly addition |
| `sub` | **Stable** | Poly subtraction |
| `mul` | **Stable** | Poly multiplication |
| `inv` | **Stable** | `inv` |
| `eq_mod` | **Stable** | `eq_mod` |
| `neg` | **Stable** | Poly negation |
| `mul_rational` | **Stable** | `mul_rational` |
| `is_zero` | **Stable** | Poly is zero |
| `is_one` | **Stable** | Poly is one |
| `to_rootof_expr` | **Stable** | `to_rootof_expr` |
| `coords_q` | **Pipeline private** | `coords_q` |
| `fields_same` | **Pipeline private** | `fields_same` |
| `align_pair` | **Pipeline private** | `align_pair` |
| `fold_algext_sum` | **Stable** | canonical sum of AlgExt terms |
| `fold_algext_sum_for_ctx` | **Stable (bounded)** | fold sum with Context session cache (R5 eval) |
| `fold_algext_sum_mode` | **Stable** | fold_algext_sum with mode |
| `complex_algext_parts` | **Pipeline private** | `complex_algext_parts` |
| `complex_algext_to_expr` | **Pipeline private** | `complex_algext_to_expr` |
| `fold_complex_algext_sum` | **Stable** | sum with complex AlgExtC parts |
| `fold_complex_algext_product` | **Stable** | product with complex AlgExtC parts |
| `complex_algext_mul_parts` | **Pipeline private** | `complex_algext_mul_parts` |
| `fold_algext_product` | **Stable** | canonical product of AlgExt terms |
| `fold_algext_product_for_ctx` | **Stable (bounded)** | fold product with Context session cache (R5 eval) |
| `try_rootof_to_algext` | **Stable** | Func(RootOf) → AlgExt Expr |
| `contains_algext` | **Stable** | subtree contains AlgExt or rootof |
| `try_as_algext_data` | **Stable** | view Expr as AlgExtData if present |
| `algext_square_roots` | **Stable (bounded)** | square roots in extension field |
| `algext_cube_root` | **Stable (bounded)** | cube root in extension field |
| `cube_minpoly_via_matrix` | **Pipeline private** | `cube_minpoly_via_matrix` |
| `algext_sqrt_branches` | **Stable (bounded)** | sqrt branches as Expr list |
| `common_ext` | **Stable** | common extension for two AlgExt values |
| `sqrt_minpoly_via_matrix` | **Pipeline private** | `sqrt_minpoly_via_matrix` |
| `algext_mul_squares_to_two` | **Pipeline private** | `algext_mul_squares_to_two` |
| `algext_add_neg_cancels` | **Pipeline private** | `algext_add_neg_cancels` |
| `algext_sqrt_of_sqrt2` | **Pipeline private** | `algext_sqrt_of_sqrt2` |
| `algext_cube_root_of_two` | **Pipeline private** | `algext_cube_root_of_two` |
| `algext_sqrt_of_neg_sqrt2_is_complex` | **Pipeline private** | `algext_sqrt_of_neg_sqrt2_is_complex` |
| `complex_algext_sum_cancels` | **Pipeline private** | `complex_algext_sum_cancels` |
| `algext_to_rootof_roundtrip_display` | **Pipeline private** | `algext_to_rootof_roundtrip_display` |
| `algext_inv_divides_to_one` | **Pipeline private** | `algext_inv_divides_to_one` |
| `subfield_embed_rational_into_sqrt2` | **Pipeline private** | `subfield_embed_rational_into_sqrt2` |
| `common_ext_sqrt2_cbrt2` | **Pipeline private** | `common_ext_sqrt2_cbrt2` |
| `algext_add_reverse_order_after_common_cache` | **Pipeline private** | `algext_add_reverse_order_after_common_cache` |
| `fold_algext_sum_merges_equal_fields_without_ptr_eq` | **Pipeline private** | `fold_algext_sum_merges_equal_fields_without_ptr_eq` |
| `from_coords_q_preserves_k2_tensor_coords_roundtrip` | **Pipeline private** | `from_coords_q_preserves_k2_tensor_coords_roundtrip` |
| `fold_algext_sum_rat_on_k2_merges_via_embed_rational` | **Pipeline private** | `fold_algext_sum_rat_on_k2_merges_via_embed_rational` |
| `embed_rational_on_k2_places_constant_in_block_u0` | **Pipeline private** | `embed_rational_on_k2_places_constant_in_block_u0` |
| `fold_algext_sum_canonical_order_independent` | **Pipeline private** | `fold_algext_sum_canonical_order_independent` |
| `fold_algext_sum_split_two_fields_merges_via_lazy_common` | **Pipeline private** | `fold_algext_sum_split_two_fields_merges_via_lazy_common` |
| `fold_algext_sum_rat_on_k1_still_merges_into_algext` | **Pipeline private** | `fold_algext_sum_rat_on_k1_still_merges_into_algext` |
| `algext_frac_via_eval` | **Pipeline private** | `algext_frac_via_eval` |

### `alg_ext_c.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `zero` | **Stable** | Poly zero |
| `one` | **Stable** | Poly one |
| `from_alg_ext` | **Stable** | `Poly::from_alg_ext` |
| `from_complex_parts` | **Stable** | `from_complex_parts` |
| `from_coords_q` | **Pipeline private** | `from_coords_q` |
| `re_q` | **Pipeline private** | `re_q` |
| `im_q` | **Pipeline private** | `im_q` |
| `ensure_same_field` | **Pipeline private** | `ensure_same_field` |
| `add` | **Stable** | Poly addition |
| `sub` | **Stable** | Poly subtraction |
| `mul` | **Stable** | Poly multiplication |
| `neg` | **Stable** | Poly negation |
| `inv` | **Stable** | `inv` |
| `eq_mod` | **Stable** | `eq_mod` |
| `is_zero` | **Stable** | Poly is zero |
| `is_one` | **Stable** | Poly is one |
| `into_expr` | **Stable** | `into_expr` |
| `to_expr` | **Stable** | `to_expr` |
| `re_to_legacy_expr` | **Pipeline private** | `re_to_legacy_expr` |
| `im_to_legacy_expr` | **Pipeline private** | `im_to_legacy_expr` |
| `align_pair` | **Stable** | `align_pair` |
| `align_with` | **Pipeline private** | `align_with` |
| `canonicalize_to_algext_c` | **Stable** | Expr → canonical AlgExtCData |
| `expr_to_field_element` | **Pipeline private** | `expr_to_field_element` |
| `algext_c_from_algext_is_real` | **Pipeline private** | `algext_c_from_algext_is_real` |
| `i_times_sqrt2_squared_is_minus_two` | **Pipeline private** | `i_times_sqrt2_squared_is_minus_two` |
| `canonicalize_algext_roundtrip` | **Pipeline private** | `canonicalize_algext_roundtrip` |
| `canonicalize_complex_with_algext_im` | **Pipeline private** | `canonicalize_complex_with_algext_im` |

### `ext_tower.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `next_field_id` | **Pipeline private** | `next_field_id` |
| `dimension` | **Stable** | `Poly::dimension` |
| `is_simple_over_q` | **Stable** | `Poly::is_simple_over_q` |
| `eq` | **Stable** | `Poly::eq` |
| `is_base` | **Stable** | `Poly::is_base` |
| `rational` | **Stable** | `Poly::rational` |
| `id` | **Stable** | `id` |
| `tower` | **Stable** | `tower` |
| `dimension` | **Stable** | `dimension` |
| `parent_field` | **Stable** | `parent_field` |
| `adjoin_irreducible` | **Stable** | `adjoin_irreducible` |
| `adjoin_irreducible_parent_coeffs` | **Stable** | `adjoin_irreducible_parent_coeffs` |
| `adjoin_irreducible_over_q` | **Stable** | `adjoin_irreducible_over_q` |
| `is_subfield_of` | **Stable** | `is_subfield_of` |
| `try_subfield_embedding` | **Partial** | optional algorithm path `try_subfield_embedding` |
| `top_min_poly_exprs` | **Stable** | `top_min_poly_exprs` |
| `zero_coords` | **Stable** | `zero_coords` |
| `one_coords` | **Stable** | `one_coords` |
| `generator_coords` | **Stable** | `generator_coords` |
| `uses_tower_arithmetic` | **Pipeline private** | `uses_tower_arithmetic` |
| `embed_rational` | **Stable** | `embed_rational` |
| `layer_ext_degree` | **Pipeline private** | `layer_ext_degree` |
| `layer_minpoly_parent_coeffs` | **Pipeline private** | `layer_minpoly_parent_coeffs` |
| `unflatten_layer_blocks` | **Pipeline private** | `unflatten_layer_blocks` |
| `flatten_layer_element` | **Pipeline private** | `flatten_layer_element` |
| `element_add` | **Stable** | `element_add` |
| `element_add_tower` | **Pipeline private** | `element_add_tower` |
| `element_add_primitive` | **Pipeline private** | `element_add_primitive` |
| `element_sub` | **Stable** | `element_sub` |
| `element_sub_tower` | **Pipeline private** | `element_sub_tower` |
| `element_sub_primitive` | **Pipeline private** | `element_sub_primitive` |
| `element_neg` | **Stable** | `element_neg` |
| `element_neg_tower` | **Pipeline private** | `element_neg_tower` |
| `element_neg_primitive` | **Pipeline private** | `element_neg_primitive` |
| `element_mul` | **Stable** | `element_mul` |
| `element_mul_tower` | **Pipeline private** | `element_mul_tower` |
| `element_mul_primitive` | **Pipeline private** | `element_mul_primitive` |
| `element_inv` | **Stable** | `element_inv` |
| `element_inv_tower` | **Pipeline private** | `element_inv_tower` |
| `element_inv_primitive` | **Pipeline private** | `element_inv_primitive` |
| `element_eq_mod` | **Stable** | `element_eq_mod` |
| `element_is_zero` | **Stable** | `element_is_zero` |
| `element_is_one` | **Stable** | `element_is_one` |
| `ensure_same_field_len` | **Pipeline private** | `ensure_same_field_len` |
| `common_over_q` | **Stable** | ephemeral per-call cache；`FieldSession` 显式 dedup |
| `common_over_q_with_cache` | **Stable (bounded)** | explicit `CommonCache` |
| `embedding_for` | **Stable** | `embedding_for` |
| `align_elements` | **Stable** | ephemeral per-call cache（R5b） |
| `align_elements_with_cache` | **Stable (bounded)** | explicit session cache |
| `apply` | **Stable** | `Poly::apply` |
| `identity` | **Stable** | `Poly::identity` |
| `new` | **Pipeline private** | `new` |
| `min_poly_key_bytes` | **Pipeline private** | adjoin / base dedup key (R4) |
| `adjoin_cache_key_rational` | **Pipeline private** | session adjoin dedup key |
| `build_adjoin_irreducible` | **Pipeline private** | construct adjoin without dedup |
| `build_adjoin_parent_coeffs` | **Pipeline private** | parent-coeff adjoin construct |
| `build_base_extension_uncached` | **Pipeline private** | ℚ(θ) field without dedup |
| `get_or_create_base_by_min_poly` | **Pipeline private** | flatten compositum field dedup |
| `common_over_q_in_cache` | **Pipeline private** | `common_over_q_in_cache` |
| `align_elements_in_cache` | **Pipeline private** | `align_elements_in_cache` |
| `compute_common_flatten_for_test` | **Stable** | `compute_common_flatten_for_test` |
| `common_cache_len_for_test` | **Stable** | `common_cache_len_for_test` |
| `duplicate_field_arc_for_test` | **Stable** | `duplicate_field_arc_for_test` |
| `layer_min_poly_exprs` | **Stable** | R3：塔顶 layer minpoly（`AlgExtData::min_poly`） |
| `top_min_poly_exprs` | **Stable** | ℚ-flatten primitive poly（roots 元数据） |
| `layer_minpoly_rational_constants` | **Pipeline private** | `layer_minpoly_rational_constants` |
| `flatten_min_poly_over_q` | **Pipeline private** | explicit ℚ-flatten (+ optional session memo) |
| `flatten_min_poly_over_q_cold` | **Pipeline private** | cold flatten without memo |
| `compose_min_poly_over_q` | **Pipeline private** | `compose_min_poly_over_q` |
| `rational_subfield_embedding` | **Pipeline private** | `rational_subfield_embedding` |
| `flatten_layer_blocks` | **Pipeline private** | `flatten_layer_blocks` |
| `direct_adjoin_parent_embedding` | **Pipeline private** | `direct_adjoin_parent_embedding` |
| `compose_field_embeddings` | **Pipeline private** | `compose_field_embeddings` |
| `compute_common_dispatch` | **Pipeline private** | `compute_common_dispatch` |
| `subfield_common_pair` | **Pipeline private** | `subfield_common_pair` |
| `tower_common_eligible` | **Pipeline private** | `tower_common_eligible` |
| `pick_tower_adjoin_parent` | **Pipeline private** | `pick_tower_adjoin_parent` |
| `embedding_for_common_operand` | **Pipeline private** | `embedding_for_common_operand` |
| `compute_common_tower` | **Pipeline private** | `compute_common_tower` |
| `simple_over_q_embedding` | **Pipeline private** | `simple_over_q_embedding` |
| `embed_simple_over_q_coords` | **Pipeline private** | `embed_simple_over_q_coords` |
| `tower_adjoin_parent_for_test` | **Stable** | `tower_adjoin_parent_for_test` |
| `compute_common_flatten` | **Pipeline private** | `compute_common_flatten` |
| `embed_rationals_into` | **Pipeline private** | `embed_rationals_into` |
| `common_primitive_sum` | **Pipeline private** | `common_primitive_sum` |
| `embedding_matrix_from_theta` | **Pipeline private** | `embedding_matrix_from_theta` |
| `embedding_matrix_from_theta_block` | **Pipeline private** | `embedding_matrix_from_theta_block` |
| `embed_in_gamma_vector` | **Pipeline private** | `embed_in_gamma_vector` |
| `embed_coords` | **Stable** | `embed_coords` |
| `rational_field_dimension_one` | **Pipeline private** | `rational_field_dimension_one` |
| `t3a_adjoin_k1_u2_minus_sqrt2_has_dimension_four` | **Pipeline private** | `t3a_adjoin_k1_u2_minus_sqrt2_has_dimension_four` |
| `adjoin_cbrt2_generator_cubes_to_two` | **Pipeline private** | `adjoin_cbrt2_generator_cubes_to_two` |
| `adjoin_sqrt2_dimension_two` | **Pipeline private** | `adjoin_sqrt2_dimension_two` |
| `common_sqrt2_cbrt2_has_degree_six` | **Pipeline private** | `common_sqrt2_cbrt2_has_degree_six` |
| `common_cache_identity_is_fast` | **Pipeline private** | `common_cache_identity_is_fast` |
| `common_cache_hits_same_pair` | **Pipeline private** | `common_cache_hits_same_pair` |
| `field_registry_dedup_same_minpoly` | **Pipeline private** | renamed `adjoin_cache_dedup_same_minpoly` (R4) |
| `embedding_for_matches_source_not_operand_order` | **Pipeline private** | `embedding_for_matches_source_not_operand_order` |
| `align_elements_reverse_order_after_cache_warm` | **Pipeline private** | `align_elements_reverse_order_after_cache_warm` |
| `t1a_adjoin_base_sqrt2_matches_legacy` | **Pipeline private** | `t1a_adjoin_base_sqrt2_matches_legacy` |
| `t1b_adjoin_k1_u2_minus_3_has_dimension_four` | **Pipeline private** | `t1b_adjoin_k1_u2_minus_3_has_dimension_four` |
| `t1_try_subfield_embedding_k1_into_k2` | **Pipeline private** | `t1_try_subfield_embedding_k1_into_k2` |
| `t4b_common_subfield_is_superfield_not_compositum` | **Pipeline private** | `t4b_common_subfield_is_superfield_not_compositum` |
| `t2_align_subfield_does_not_grow_common_cache` | **Pipeline private** | `t2_align_subfield_does_not_grow_common_cache` |
| `t2_algext_add_same_generator_in_superfield_no_common_cache` | **Pipeline private** | `t2_algext_add_same_generator_in_superfield_no_common_cache` |
| `k2_embedded_sqrt2_squared_is_two` | **Pipeline private** | `k2_embedded_sqrt2_squared_is_two` |
| `t3_k1_sqrt2_squared_is_two` | **Pipeline private** | `t3_k1_sqrt2_squared_is_two` |
| `t3_k2_sqrt2_beta_times_beta_is_three_sqrt2` | **Pipeline private** | `t3_k2_sqrt2_beta_times_beta_is_three_sqrt2` |
| `common_q_sqrt2_without_tower_common` | **Pipeline private** | `common_q_sqrt2_without_tower_common` |
| `common_sqrt2_sqrt3_has_flatten_parent_none` | **Pipeline private** | `common_sqrt2_sqrt3_has_flatten_parent_none` |
| `common_sqrt2_sqrt3_has_tower_parent` | **Pipeline private** | `common_sqrt2_sqrt3_has_tower_parent` |
| `align_sqrt2_plus_sqrt3_reverse_order` | **Pipeline private** | `align_sqrt2_plus_sqrt3_reverse_order` |
| `common_cache_hits_after_tower_common` | **Pipeline private** | `common_cache_hits_after_tower_common` |
| `align_sqrt2_cbrt2_dim_six` | **Pipeline private** | `align_sqrt2_cbrt2_dim_six` |
| `tower_common_invariants_sqrt2_sqrt3` | **Pipeline private** | `tower_common_invariants_sqrt2_sqrt3` |
| `tower_common_matches_flatten_minpoly_on_sqrt2_sqrt3` | **Pipeline private** | `tower_common_matches_flatten_minpoly_on_sqrt2_sqrt3` |

### `field_arith.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `dense_q` | **Pipeline private** | `dense_q` |
| `trim_leading_zero` | **Stable** | `trim_leading_zero` |
| `poly_degree` | **Stable** | `poly_degree` |
| `poly_add` | **Stable** | `poly_add` |
| `poly_sub` | **Stable** | `poly_sub` |
| `poly_mul` | **Stable** | `poly_mul` |
| `poly_neg` | **Stable** | `poly_neg` |
| `poly_scale` | **Stable** | `poly_scale` |
| `poly_reduce` | **Stable** | `poly_reduce` |
| `poly_inv_mod` | **Stable** | `poly_inv_mod` |
| `poly_ext_gcd` | **Stable** | `poly_ext_gcd` |
| `poly_divrem` | **Stable** | `poly_divrem` |
| `pad_to_len` | **Stable** | `pad_to_len` |
| `generator_coords` | **Stable** | `generator_coords` |
| `expr_to_ratio` | **Stable** | `expr_to_ratio` |
| `rationalize_poly1` | **Stable** | `rationalize_poly1` |
| `coords_to_expr` | **Stable** | `coords_to_expr` |
| `ratio_to_expr_arc` | **Stable** | `ratio_to_expr_arc` |
| `min_poly_exprs_to_q` | **Stable** | `min_poly_exprs_to_q` |
| `poly1_coeffs` | **Stable** | `poly1_coeffs` |
| `canonical_poly1_expr` | **Stable** | `canonical_poly1_expr` |
| `minpoly_at_square` | **Stable** | `minpoly_at_square` |
| `minpoly_at_cube` | **Stable** | `minpoly_at_cube` |
| `embed_in_cube_extension` | **Stable** | `embed_in_cube_extension` |
| `embed_in_square_extension` | **Stable** | `embed_in_square_extension` |
| `identity_matrix` | **Stable** | `identity_matrix` |
| `mat_mul` | **Stable** | `mat_mul` |
| `kron_left` | **Stable** | `kron_left` |
| `mult_matrix_of_element` | **Stable** | `mult_matrix_of_element` |
| `parent_basis_vector` | **Stable** | `parent_basis_vector` |
| `tower_flat_index` | **Pipeline private** | `tower_flat_index` |
| `tower_blocks_to_flat` | **Pipeline private** | `tower_blocks_to_flat` |
| `mul_tower_basis_by_layer_generator` | **Pipeline private** | `mul_tower_basis_by_layer_generator` |
| `mult_matrix_of_adjoin_generator` | **Stable** | `mult_matrix_of_adjoin_generator` |
| `char_poly_matrix` | **Stable** | `char_poly_matrix` |
| `newton_char_poly` | **Pipeline private** | `newton_char_poly` |
| `apply_linear_map` | **Stable** | `apply_linear_map` |
| `coords_all_zero` | **Stable** | `coords_all_zero` |
| `coords_is_one` | **Stable** | `coords_is_one` |
| `new` | **Stable** | `new` |
| `zero` | **Stable** | Poly zero |
| `one` | **Stable** | Poly one |
| `is_zero` | **Stable** | Poly is zero |
| `add` | **Stable** | Poly addition |
| `sub` | **Stable** | Poly subtraction |
| `neg` | **Stable** | Poly negation |
| `mul` | **Stable** | Poly multiplication |
| `inv` | **Pipeline private** | `inv` |
| `poly_add_with_coeffs_in_field` | **Stable** | `poly_add_with_coeffs_in_field` |
| `poly_mul_with_coeffs_in_field` | **Stable** | `poly_mul_with_coeffs_in_field` |
| `poly_reduce_with_coeffs_in_field` | **Stable** | `poly_reduce_with_coeffs_in_field` |
| `poly_sub_with_coeffs_in_field` | **Stable** | `poly_sub_with_coeffs_in_field` |
| `poly_neg_with_coeffs_in_field` | **Stable** | `poly_neg_with_coeffs_in_field` |
| `poly_inv_mod_with_coeffs_in_field` | **Stable** | `poly_inv_mod_with_coeffs_in_field` |

### `poly.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `var` | **Stable** | Poly univariate generator |
| `poly_from_expr_shape` | **Pipeline private** | `poly_from_expr_shape` |
| `rational_leaf` | **Pipeline private** | `rational_leaf` |
| `rational_poly_ring` | **Pipeline private** | `rational_poly_ring` |
| `algext_leaf` | **Pipeline private** | `algext_leaf` |
| `algext_poly_ring` | **Pipeline private** | `algext_poly_ring` |
| `expr_to_poly` | **Stable** | path A in [expr-poly-conversion.md](../../../../.doc/expr-poly-conversion.md). |
| `poly_alg_from_expr` | **Stable** | path B in [expr-poly-conversion.md](../../../../.doc/expr-poly-conversion.md). |
| `assemble_poly_expr` | **Pipeline private** | `assemble_poly_expr` |
| `poly_to_expr` | **Stable** | `Poly` over ℚ → `Expr` (coefficients remain rational). |
| `algext_poly_to_expr` | **Stable** | `Poly<AlgExtC>` → `Expr` (coefficients as `AlgExt` / `rootof` / `AlgExtC`). |
| `univariate_poly_to_poly1_expr` | **Stable** | univariate `Poly` w.r.t. `var` → `poly1[coeffs…]` Expr (giac high-degree-first order). |
| `poly_mod_to_expr` | **Stable** | PolyMod → Expr |
| `u64_to_expr_int` | **Stable** | `u64_to_expr_int` |
| `ratio_to_expr` | **Stable** | `Ratio<BigInt>` → `Expr` (`Int` / `Rat`). |
| `monomial_to_expr` | **Pipeline private** | `monomial_to_expr` |
| `vars_from_expr` | **Stable** | sorted variables in Expr |
| `collect_vars` | **Pipeline private** | `collect_vars` |
| `gcd_linear_bridge` | **Pipeline private** | `gcd_linear_bridge` |
| `expr_to_poly_matches_expanded_power` | **Pipeline private** | `expr_to_poly_matches_expanded_power` |
| `expr_to_poly_high_degree_within_limit` | **Pipeline private** | `expr_to_poly_high_degree_within_limit` |
| `expr_to_poly_rejects_algext` | **Pipeline private** | `expr_to_poly_rejects_algext` |
| `expr_to_poly_rejects_nested_rootof` | **Pipeline private** | `expr_to_poly_rejects_nested_rootof` |
| `poly_alg_from_expr_rejects_rational` | **Pipeline private** | `poly_alg_from_expr_rejects_rational` |
| `poly_alg_from_expr_constant_rootof` | **Pipeline private** | `poly_alg_from_expr_constant_rootof` |
| `poly_alg_from_expr_x_squared_minus_two_over_k` | **Pipeline private** | `poly_alg_from_expr_x_squared_minus_two_over_k` |
| `flat_uni_algext_degree` | **Pipeline private** | `flat_uni_algext_degree` |
| `poly_alg_from_expr_univariate_roundtrip` | **Pipeline private** | `poly_alg_from_expr_univariate_roundtrip` |
| `univariate_poly_to_poly1_expr_quadratic` | **Pipeline private** | `univariate_poly_to_poly1_expr_quadratic` |

### `poly_alg_coeff.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `new` | **Stable** | `new` |
| `as_inner` | **Stable** | `Poly::as_inner` |
| `into_inner` | **Stable** | `Poly::into_inner` |
| `fmt` | **Stable** | `Poly::fmt` |
| `from` | **Stable** | `Poly::from` |
| `coeff_zero` | **Stable** | `Poly::coeff_zero` |
| `coeff_one` | **Stable** | `Poly::coeff_one` |
| `coeff_is_zero` | **Pipeline private** | `coeff_is_zero` |
| `coeff_is_one` | **Pipeline private** | `coeff_is_one` |
| `coeff_add` | **Pipeline private** | `coeff_add` |
| `coeff_sub` | **Pipeline private** | `coeff_sub` |
| `coeff_neg` | **Pipeline private** | `coeff_neg` |
| `coeff_mul` | **Pipeline private** | `coeff_mul` |
| `coeff_div` | **Pipeline private** | `coeff_div` |
| `algext_c_coeff_mul_sqrt2` | **Pipeline private** | `algext_c_coeff_mul_sqrt2` |

### `poly_conv.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `expr_contains_alg_coeff` | **Stable** | predicate for callers choosing `expr_to_poly` vs `poly_alg_from_expr`. |
| `rational_poly_reject` | **Stable** | `rational_poly_reject` |
| `first_alg_coeff_message` | **Pipeline private** | `first_alg_coeff_message` |
| `detects_rootof_in_sum` | **Pipeline private** | `detects_rootof_in_sum` |
| `detects_algext_in_product` | **Pipeline private** | `detects_algext_in_product` |
| `rational_poly_has_no_alg_coeff` | **Pipeline private** | `rational_poly_has_no_alg_coeff` |

### `field_session.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `FieldSession` | **Stable (bounded)** | ambient K + working L；interior `RefCell` on caches / working |
| `new` | **Stable** | independent extension caches |
| `fork_ambient` | **Stable (bounded)** | fresh K/L; share `common` / adjoin caches with parent session (R5) |
| `adjoin_irreducible` | **Stable (bounded)** | session adjoin dedup (R4) |
| `adjoin_irreducible_parent_coeffs` | **Stable (bounded)** | T3+ parent-coeff adjoin dedup |
| `common_over_q` | **Stable (bounded)** | `common_over_q` with session cache |
| `align_elements` | **Stable (bounded)** | align coords with session cache |
| `common_cache_len` | **Stable (bounded)** | diagnostics / tests |
| `ambient` | **Stable** | ambient K (normalized poly coefficients) |
| `working` | **Stable** | current working L (`Arc` clone) |
| `zero` / `one` / `int` / `half` | **Stable** | constants in L |
| `lift` | **Stable** | embed coeff into working L |
| `align` | **Stable** | align pair; bump L to common field |
| `adjoin_sqrt` | **Stable (bounded)** | adjoin √u (real); return generator in L |
| `adjoin_cbrt` | **Stable (bounded)** | adjoin ∛u (real); return generator in L |
| `adjoin_primitive_cube_root_of_unity` | **Stable (bounded)** | adjoin ω with ω²+ω+1=0; pure cubic roots |
| `add` / `mul` / `div` / `neg` | **Stable** | arithmetic with auto-align |
| `bump_to` | **Pipeline private** | grow working L toward superfield or align target |
| `with_common_cache` | **Pipeline private** | reentrant cache closure (static API) |
| `get_or_create_base_by_min_poly` | **Pipeline private** | flatten path ℚ(θ) dedup |
| `set_working` | **Pipeline private** | Euler branch checkpoint |

### `poly_roots.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `poly_algext_roots` | **Stable (bounded)** | exact AlgExtC roots deg 1–4; own `FieldSession::new` |
| `poly_algext_roots_for_ctx` | **Stable (bounded)** | roots with `ctx.session()` shared extension cache |
| `infer_field` | **Pipeline private** | infer ambient K from PolyAlgExt coefficients |
| `coeff_at` | **Stable** | univariate coefficient at exponent (uses session.zero) |
| `normalize_coeffs` | **Pipeline private** | lift all coeffs to ambient K via session |
| `monomial_to_poly` | **Pipeline private** | `monomial_to_poly` |
| `monic_univariate` | **Pipeline private** | `monic_univariate` |
| `linear_root` | **Pipeline private** | `linear_root` |
| `quadratic_roots` | **Pipeline private** | quadratic via session + sqrt_disc |
| `sqrt_disc` | **Pipeline private** | √Δ via session.adjoin_sqrt; imaginary via mul_i |
| `mul_i` | **Pipeline private** | formal i·z on working field |
| `cubic_roots` | **Pipeline private** | pure cubic or Cardano+deflate |
| `pure_cubic_roots` | **Pipeline private** | t³+a₀: β=∛(−a₀), roots β·ω^k |
| `dedup_roots` | **Pipeline private** | merge roots modulo eq_mod |
| `one_cubic_root` | **Pipeline private** | one Cardano root (session adjoin) |
| `quartic_roots` | **Pipeline private** | resolvent cubic + split |
| `biquadratic_roots` | **Pipeline private** | `biquadratic_roots` |
| `depress_quartic` | **Pipeline private** | `depress_quartic` |
| `build_resolvent_cubic` | **Pipeline private** | R(z)=z³−pz²−4rz+(4pr−q²) on session |
| `split_depressed_quartic` | **Pipeline private** | split via session.adjoin_sqrt |
| `deflate_monic` | **Pipeline private** | `deflate_monic` |
| `algext_c_sqrt` | **Pipeline private** | legacy sqrt for biquadratic u-roots |
| `algext_c_cube_root` | **Pipeline private** | legacy cbrt (biquadratic path) |
| `verify_root` | **Pipeline private** | monic+normalize eval via session |
| `cubic_one_root_vanishes` | **B** | one_cubic_root Cardano; verify_root |
| `quadratic_sqrt_four_times_sqrt2_over_k1` | **B** | session.adjoin_sqrt on K₁ |
| `quadratic_x2_minus_sqrt2_roots_vanish` | **B** | quadratic_roots over K₁ |
| `roots_quadratic_x2_minus_2` | **B** | poly_algext_roots deg=2 over ℚ |
| `roots_quadratic_x2_minus_sqrt2_over_k` | **B** | poly_algext_roots over ℚ(√2) |
| `roots_cubic_t3_minus_2` | **B** | pure cubic: 3 roots β·ω^k; verify_root each |
| `resolvent_golden_t4_plus_t_plus_1` | **B** | PR-D′: depressed (0,1,1) → R(z)=z³−4z−1 |
| `depressed_t4_plus_t_plus_1_coeffs` | **B** | PR-D′: t⁴+t+1 depression coeffs |
| `roots_quartic_t4_plus_t_plus_1` | **B** | quartic e2e (ignored: PR-E′ perf) |
| `roots_biquadratic_t4_minus_2` | **B** | biquadratic four roots |
