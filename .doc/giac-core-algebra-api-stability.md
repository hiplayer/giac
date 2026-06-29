# giac-core::algebra API 稳定性分层

**规范来源:** [algorithm-expr-api.md](algorithm-expr-api.md)  
**Expr ↔ Poly:** [expr-poly-conversion.md](expr-poly-conversion.md)  
**扩域求根:** [issues/GIAC-poly-roots-field-session-plan.md](issues/GIAC-poly-roots-field-session-plan.md)  
**代码:** `giac-rs/crates/giac-core/src/algebra/*`

---

## 1. 源码注释格式

| 标记 | 用于 | 可见性 |
|------|------|--------|
| `/// **Stable** — …` | 转换 / 谓词 / 环运算 | `pub` |
| `/// **Stable (bounded)** — …` | 次数或维数受限（roots、adjoin） | `pub` |
| `/// **Partial** — …` | 窄子集或启发式 | `pub` |
| `// **Pipeline private** — …` | 管线步骤、形状检测 | `fn` 私有 |

**未标注私有 fn:** `python3 scripts/annotate_api_tiers.py`；Per-file 表见本文末尾。

---

## 2. 模块公开 API（代表性）

### `poly.rs` / `poly_conv.rs` — Expr ↔ Poly

| 函数 | 层级 | 说明 |
|------|------|------|
| `expr_to_poly` | **Stable** | 路径 A：ℚ 多项式；拒绝 `AlgExt` / `rootof` |
| `poly_to_expr` | **Stable** | `Poly` → `Expr` |
| `poly_alg_from_expr` | **Stable** | 路径 B：`Poly<AlgExtCPolyCoeff>` |
| `algext_poly_to_expr` | **Stable** | `PolyAlgExt` → `Expr` |
| `expr_contains_alg_coeff` | **Stable** | 选择路径 A vs B |
| `poly_algext_from_poly` | **Stable** | ℚ `Poly` 嵌入 `PolyAlgExt` |

### `poly_roots.rs` — 扩域求根

| 函数 | 层级 | 说明 |
|------|------|------|
| `poly_algext_roots` | **Stable (bounded)** | deg 1–4 精确根；deg ≥ 5 → `NotImplemented` |
| `poly_algext_roots_for_ctx` | **Stable (bounded)** | 同上 + `Context` session cache（R5） |

**入口管线（1B）：** `infer_field → normalize_coeffs → monic_univariate → FieldSession`；仅经私有 `PolyInK::prepare` 进入 `roots_dispatch`。

### `poly_alg_sturm.rs` — Sturm over K (P3-4 / B-04)

| 函数 | 层级 | 说明 |
|------|------|------|
| `sturm_sequence_wrt_algext` | **Stable** | sqff(p), P₁=p′, 链式 `div_rem_wrt_algext` |
| `sturm_sign_variations_at_algext` | **Stable** | V(a) at rational a ∈ ℚ ⊂ K |
| `sturmab_count_wrt_algext` | **Stable** | root count in (a, b] |
| `sturmab_count_rational_poly` | **Stable** | ℚ[var] lift → ambient K Sturm count |

**边界：** 求值点限于 ℚ 嵌入；`realroot` 全 Sturm 区间隔离仍待接。

### `poly_alg_resultant.rs` — Resultant over K (P3-4 / B-04)

| 函数 | 层级 | 说明 |
|------|------|------|
| `resultant_wrt_algext` | **Stable** | 一元 Sylvester  resultant；标量 ∈ K |

### `field_session.rs` — 显式 K / L

| 函数 | 层级 | 说明 |
|------|------|------|
| `FieldSession::new` | **Stable** | ambient **K** + working **L** |
| `FieldSession::fork_ambient` | **Stable (bounded)** | 新 K，共享 extension cache |
| `.int` / `.half` / `.zero` / `.one` | **Stable** | working **L** 上常数 |
| `.lift` / `.align` / `.add` / `.mul` / `.div` | **Stable** | 对齐后算术 |
| `.adjoin_sqrt` / `.adjoin_cbrt` / `.adjoin_irreducible` | **Stable (bounded)** | 单调扩域 |

### `alg_ext.rs` — Expr 层 AlgExt

| 函数 | 层级 | 说明 |
|------|------|------|
| `fold_algext_sum` / `fold_algext_product` | **Stable** | 同域合并 |
| `fold_algext_*_for_ctx` | **Stable** | 带 `Context` session |
| `try_as_algext_data` | **Stable** | 分解 `Expr::AlgExt` |
| `contains_algext` | **Stable** | 子树含 AlgExt |
| `common_ext` | **Stable (bounded)** | 两元公共扩域 |

---

## 3. I/O 契约

### `expr_to_poly` / `poly_alg_from_expr`

见 [expr-poly-conversion.md](expr-poly-conversion.md) §2–§3。含代数系数 **必须** 走路径 B 或保持 `Expr` 层运算。

### `poly_algext_roots`

| 字段 | 说明 |
|------|------|
| **输入** | `p: &PolyAlgExt`, `var: &Var`；`p` 不必预 normalize（入口 `PolyInK::prepare`） |
| **输出** | `Vec<AlgExtCPolyCoeff>` 根；去重 |
| **上下文** | `infer_field(p)` → ambient **K**；`FieldSession` 上 working **L** 单调扩大 |
| **边界** | deg ≤ 4；分裂域 `dim(L) ≤ 24`（`POLY_ROOTS_DIM_HARD`）；resolvent 阶段 `dim(L) ≤ 6` |
| **禁止** | 裸 `align_coeff` 链；未 normalize 多项式直入 `roots_dispatch` |

### `poly_algext_roots_for_ctx`

| 字段 | 说明 |
|------|------|
| **输入** | `p`, `var`, `ctx: &Context` |
| **缓存** | `ctx.session().fork_ambient(K)` 复用 adjoin / common cache |
| **输出** | 同 `poly_algext_roots` |

### `FieldSession`

| 不变量 | 规则 |
|--------|------|
| **ambient K** | 创建后不变 |
| **working L** | `bump_to` / adjoin 仅增大 L |
| **align** | 二元运算前对齐到公共 **L** |

---

## 4. Pipeline private（代表性）

| 模块 | 类型 / 函数 | 说明 |
|------|-------------|------|
| `poly_roots` | `PolyInK` | normalize + monic 包装；唯一进 dispatch 路径 |
| `poly_roots` | `normalize_coeffs`, `monic_univariate` | 系数 embed **K**、首一化 |
| `poly_roots` | `verify_root` | 测试用；与 `PolyInK::prepare` 同路径 |
| `ext_tower` | `ExtensionField`, adjoin 构造 | 塔与嵌入 |
| `field_arith` | `coords_to_expr`, `poly_reduce` | ℚ 坐标与多项式 gcd 辅助 |

---

## 5. 维护

```bash
cd giac-rs
python3 scripts/annotate_api_tiers.py
python3 scripts/annotate_api_tiers.py --inventory
```

PR 改动 `giac-core/src/algebra` 时更新本文 §2–§3 与 Per-file 表（[algorithm-expr-api.md §7.2](algorithm-expr-api.md#72-测试通过后提交--合入前复审)）。

---

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
| `align_pair_with_session` | **Pipeline private** | align with optional Context session (R5) |
| `add_aligned` | **Pipeline private** | add after align (optional session) |
| `mul_aligned` | **Pipeline private** | mul after align (optional session) |
| `fold_algext_sum` | **Stable** | canonical sum of AlgExt terms |
| `fold_algext_sum_for_ctx` | **Stable (bounded)** | fold sum with Context session |
| `fold_algext_sum_mode` | **Stable** | fold_algext_sum with mode |
| `fold_algext_sum_mode_impl` | **Pipeline private** | `fold_algext_sum_mode_impl` |
| `complex_algext_parts` | **Pipeline private** | `complex_algext_parts` |
| `complex_algext_to_expr` | **Pipeline private** | `complex_algext_to_expr` |
| `fold_complex_algext_sum` | **Stable** | sum with complex AlgExtC parts |
| `fold_complex_algext_product` | **Stable** | product with complex AlgExtC parts |
| `complex_algext_mul_parts` | **Pipeline private** | `complex_algext_mul_parts` |
| `fold_algext_product` | **Stable** | canonical product of AlgExt terms |
| `fold_algext_product_for_ctx` | **Stable (bounded)** | fold product with Context session |
| `fold_algext_product_impl` | **Pipeline private** | `fold_algext_product_impl` |
| `rootof_from_minpoly` | **Stable** | `rootof([n₀, n₁, …], minpoly)` as an `Expr`. |
| `quadratic_rootof_branches` | **Stable (bounded)** | two `rootof` branches `±α` for quadratic irrational roots. |
| `try_rootof_to_algext` | **Stable** | Func(RootOf) → AlgExt Expr |
| `contains_algext` | **Stable** | subtree contains AlgExt or rootof |
| `try_as_algext_data` | **Stable** | view Expr as AlgExtData if present |
| `algext_square_roots` | **Stable (bounded)** | square roots in extension field (blind adjoin) |
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
| `min_poly_is_layer_not_flatten_over_nested_adjoin` | **Pipeline private** | `min_poly_is_layer_not_flatten_over_nested_adjoin` |
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
| `align_pair_in_cache` | **Stable** | `align_pair_in_cache` |
| `align_with` | **Pipeline private** | `align_with` |
| `align_with_in_cache` | **Pipeline private** | `align_with_in_cache` |
| `canonicalize_to_algext_c` | **Stable** | Expr → canonical AlgExtCData |
| `expr_to_field_element` | **Pipeline private** | `expr_to_field_element` |
| `algext_c_from_algext_is_real` | **Pipeline private** | `algext_c_from_algext_is_real` |
| `i_times_sqrt2_squared_is_minus_two` | **Pipeline private** | `i_times_sqrt2_squared_is_minus_two` |
| `canonicalize_algext_roundtrip` | **Pipeline private** | `canonicalize_algext_roundtrip` |
| `canonicalize_complex_with_algext_im` | **Pipeline private** | `canonicalize_complex_with_algext_im` |

### `common_minimal.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `new` | **Pipeline private** | `new` |
| `min_g` | **Pipeline private** | `min_g` |
| `k` | **Pipeline private** | `k` |
| `mat_theta` | **Pipeline private** | `mat_theta` |
| `w_a` | **Pipeline private** | `w_a` |
| `w_b` | **Pipeline private** | `w_b` |
| `na` | **Pipeline private** | `na` |
| `nb` | **Pipeline private** | `nb` |
| `dim` | **Pipeline private** | `dim` |
| `side_spec_for` | **Pipeline private** | `side_spec_for` |
| `tensor_index_exp` | **Pipeline private** | `tensor_index_exp` |
| `tensor_index` | **Pipeline private** | `tensor_index` |
| `add_tensor_term_rational` | **Pipeline private** | bivariate reduction mod rational `ma`, `mb` (U1a). |
| `mul_by_gen_a_rational` | **Pipeline private** | `mul_by_gen_a_rational` |
| `mul_by_gen_b_rational` | **Pipeline private** | `mul_by_gen_b_rational` |
| `mul_by_theta_rational` | **Pipeline private** | `mul_by_theta_rational` |
| `build_compositum_matrix_rational` | **Pipeline private** | `build_compositum_matrix_rational` |
| `mat_b_from_spec` | **Pipeline private** | `mat_b_from_spec` |
| `build_mat_theta` | **Pipeline private** | `build_mat_theta` |
| `tensor_unit_1` | **Pipeline private** | `tensor_unit_1` |
| `tensor_gen_a` | **Pipeline private** | `tensor_gen_a` |
| `tensor_gen_b` | **Pipeline private** | `tensor_gen_b` |
| `build_compositum_matrix` | **Pipeline private** | `build_compositum_matrix` |
| `compositum_rref_solvable` | **Pipeline private** | `compositum_rref_solvable` |
| `minpoly_from_rref` | **Pipeline private** | `minpoly_from_rref` |
| `embed_poly_from_rref_col` | **Pipeline private** | `embed_poly_from_rref_col` |
| `common_minimal_poly_core` | **Pipeline private** | `common_minimal_poly_core` |
| `common_minimal_poly_over_q` | **Pipeline private** | `common_minimal_poly_over_q` |
| `common_minimal_poly_over_parent` | **Pipeline private** | `common_minimal_poly_over_parent` |
| `element_pow_coords` | **Pipeline private** | `element_pow_coords` |
| `invert_rational_matrix` | **Pipeline private** | invert square matrix over ℚ (dim ≤ compositum bound). |
| `power_basis_invertible` | **Pipeline private** | `power_basis_invertible` |
| `nested_primitive_from_k` | **Pipeline private** | `nested_primitive_from_k` |
| `embedding_passes_generator_gate` | **Pipeline private** | `embedding_passes_generator_gate` |
| `power_basis_to_ops` | **Pipeline private** | `power_basis_to_ops` |
| `try_nested_tensor_ops_embed` | **Pipeline private** | `try_nested_tensor_ops_embed` |
| `nested_embedding_from_compositum` | **Pipeline private** | `nested_embedding_from_compositum` |
| `flat_primitive_candidates` | **Pipeline private** | `flat_primitive_candidates` |
| `embed_side` | **Pipeline private** | `embed_side` |
| `embed_simple_power` | **Pipeline private** | `embed_simple_power` |
| `simple_layer_minpoly` | **Pipeline private** | `simple_layer_minpoly` |
| `primitive_minpoly_for_nested` | **Pipeline private** | `primitive_minpoly_for_nested` |
| `charpoly_of_element_op` | **Pipeline private** | `charpoly_of_element_op` |
| `minpoly_and_prim_for_common` | **Pipeline private** | `minpoly_and_prim_for_common` |
| `common_operand_degree` | **Pipeline private** | `common_operand_degree` |
| `order_operands_for_common` | **Pipeline private** | `order_operands_for_common` |
| `is_embedded_rational` | **Pipeline private** | `is_embedded_rational` |
| `coeff_as_embedded_rational` | **Pipeline private** | `coeff_as_embedded_rational` |
| `factor_to_minpoly_q` | **Pipeline private** | `factor_to_minpoly_q` |
| `factor_to_parent_blocks` | **Pipeline private** | `factor_to_parent_blocks` |
| `compute_compositum_upstream` | **Pipeline private** | `compute_compositum_upstream` |
| `compute_compositum_upstream_inner` | **Pipeline private** | `compute_compositum_upstream_inner` |
| `compute_common_minimal_pair` | **Pipeline private** | `compute_common_minimal_pair` |
| `compute_common_minimal_pair_inner` | **Pipeline private** | `compute_common_minimal_pair_inner` |
| `compute_compositum_upstream_for_test` | **Pipeline private** | `compute_compositum_upstream_for_test` |
| `compute_common_minimal_pair_for_test` | **Pipeline private** | `compute_common_minimal_pair_for_test` |
| `common_minimal_poly_over_parent_for_test` | **Pipeline private** | `common_minimal_poly_over_parent_for_test` |

### `compositum_session.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `with_compositum_session` | **Pipeline private** | `with_compositum_session` |

### `ext_tower.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `next_field_id` | **Pipeline private** | `next_field_id` |
| `layer_minpoly_block_to_expr` | **Pipeline private** | one layer-minpoly coefficient (parent-field element) as `Expr` |
| `dimension` | **Stable** | `Poly::dimension` |
| `is_simple_over_q` | **Stable** | `Poly::is_simple_over_q` |
| `min_poly_key` | **Pipeline private** | `min_poly_key` (R1 semantic_key / adjoin cache) |
| `parent_blocks_key` | **Pipeline private** | `parent_blocks_key` |
| `semantic_key_bytes` | **Pipeline private** | `semantic_key_bytes` |
| `eq` | **Stable** | `Poly::eq` (R0: tower only; lazy flatten not compared) |
| `is_base` | **Stable** | `Poly::is_base` |
| `rational` | **Stable** | `Poly::rational` |
| `id` | **Stable** | `id` |
| `tower` | **Stable** | `tower` |
| `dimension` | **Stable** | `dimension` |
| `semantic_key` | **Stable (bounded)** | `semantic_key` |
| `primitive_modulus` | **Pipeline private** | `primitive_modulus` |
| `parent_field` | **Stable** | `parent_field` |
| `adjoin_irreducible` | **Stable** | `adjoin_irreducible` |
| `adjoin_irreducible_parent_coeffs` | **Stable** | `adjoin_irreducible_parent_coeffs` |
| `adjoin_irreducible_over_q` | **Stable** | `adjoin_irreducible_over_q` |
| `is_subfield_of` | **Stable** | `is_subfield_of` |
| `try_subfield_embedding` | **Partial** | optional algorithm path `try_subfield_embedding` |
| `layer_min_poly_exprs` | **Stable** | `layer_min_poly_exprs` |
| `top_min_poly_exprs` | **Stable** | `top_min_poly_exprs` |
| `zero_coords` | **Stable** | `zero_coords` |
| `one_coords` | **Stable** | `one_coords` |
| `generator_coords` | **Stable** | `generator_coords` |
| `uses_tower_arithmetic` | **Pipeline private** | `uses_tower_arithmetic` |
| `layer_arith_mode` | **Pipeline private** | `layer_arith_mode` |
| `parent_coeff_ring_capable` | **Pipeline private** | `parent_coeff_ring_capable` |
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
| `try_square_root_in_field` | **Stable (bounded)** | try square root in extension field |
| `try_square_root_in_field_shallow` | **Pipeline private** | shallow sqrt probe |
| `element_is_zero` | **Stable** | `element_is_zero` |
| `element_is_one` | **Stable** | `element_is_one` |
| `ensure_same_field_len` | **Pipeline private** | `ensure_same_field_len` |
| `compositum` | **Stable** | `compositum` |
| `compositum_with_cache` | **Stable (bounded)** | compositum with session cache |
| `embedding_for` | **Stable** | `embedding_for` |
| `align_pair` | **Stable** | `align_pair` |
| `align_pair_with_cache` | **Stable (bounded)** | align with session cache |
| `with_ephemeral_common_cache` | **Pipeline private** | ephemeral cache for static API (R5b: no cross-call dedup). |
| `compositum_in_cache` | **Pipeline private** | `compositum_in_cache` |
| `align_pair_in_cache` | **Pipeline private** | `align_pair_in_cache` |
| `apply` | **Stable** | `Poly::apply` |
| `field` | **Pipeline private** | `field` |
| `arc` | **Pipeline private** | `arc` |
| `arc_ref` | **Pipeline private** | `arc_ref` |
| `into_operand` | **Pipeline private** | `into_operand` |
| `from_arc` | **Pipeline private** | `from_arc` |
| `deref` | **Pipeline private** | `deref` |
| `identity` | **Stable** | `Poly::identity` |
| `new` | **Pipeline private** | `new` |
| `strategy` | **Pipeline private** | `strategy` |
| `pair` | **Pipeline private** | `pair` |
| `identity` | **Pipeline private** | `identity` |
| `deref` | **Pipeline private** | `deref` |
| `verify_common_pair` | **Pipeline private** | `verify_common_pair` |
| `element_pow_coords` | **Pipeline private** | `element_pow_coords` |
| `eval_rational_poly1_at_field` | **Pipeline private** | Horner eval of monic poly1 `[c_deg,…,c_0]` at field element. |
| `verify_embedded_one` | **Pipeline private** | `verify_embedded_one` |
| `verify_embedded_generator` | **Pipeline private** | layer generator + minpoly vanishing in common field. |
| `common_cache_key` | **Pipeline private** | `common_cache_key` |
| `min_poly_key_bytes` | **Pipeline private** | byte key for ℚ minpoly dedup (R4 session cache). |
| `adjoin_cache_key_rational` | **Pipeline private** | adjoin dedup key (R4). |
| `adjoin_cache_key_parent_blocks` | **Pipeline private** | parent-coeff adjoin dedup key (R4). |
| `proper_subfields_chain` | **Pipeline private** | collect proper subfields along parent chain (immediate → … → ℚ). |
| `try_preimage_under_embedding` | **Pipeline private** | solve M·v = u for embedding matrix (tgt × src); None if u ∉ Im(M). |
| `gauss_elim_rref` | **Pipeline private** | ℚ Gaussian elimination; returns pivot column per row, inconsistent flag. |
| `try_sqrt_subfield_descent` | **Pipeline private** | F5 S3: u ∈ F ⊂ L → sqrt in F, embed back. |
| `try_sqrt_quadratic_top_layer` | **Pipeline private** | F5 S4: top deg-2 layer u = a + b·g, (a+bg)² = u (standard 2×pd layout). |
| `quadratic_layer_generators` | **Pipeline private** | collect deg-2 layer generators embedded in `field`. |
| `try_sqrt_layer_generator_algebra` | **Pipeline private** | F5 S5: products / ratios of quadratic layer generators. |
| `try_sqrt_pairwise_fallback` | **Pipeline private** | S6: dim≤6 pairwise ±1 ponytail fallback. |
| `try_sqrt_layer_generators` | **Pipeline private** | `try_sqrt_layer_generators` |
| `try_sqrt_basis_squares` | **Pipeline private** | operational basis ±eᵢ. |
| `try_square_root_in_field_impl` | **Pipeline private** | F5 core: structural probes S0–S6 (no dim³ enum). |
| `try_square_root_in_field_shallow_impl` | **Pipeline private** | shallow: S1–S3 + S5 products (Euler hot path). |
| `try_sqrt_small_combo_with_quadratic_gens` | **Pipeline private** | F4: e_i ± k·g combinations for last quadratic adjoin layers (dim≤12). |
| `embed_coords_in` | **Pipeline private** | `embed_coords_in` |
| `coords_square_eq_mod` | **Pipeline private** | `coords_square_eq_mod` |
| `algext_from_coords` | **Pipeline private** | `AlgExtData` from operational coords in `field`. |
| `layer_minpoly_rational_constants` | **Pipeline private** | `layer_minpoly_rational_constants` |
| `build_adjoin_irreducible` | **Pipeline private** | construct adjoin field without session dedup (R4). |
| `build_adjoin_parent_coeffs` | **Pipeline private** | construct parent-coeff adjoin without session dedup (R4). |
| `get_or_create_base_by_min_poly` | **Pipeline private** | construct base extension without session dedup (R4/R5b). |
| `build_base_extension_uncached` | **Pipeline private** | construct base extension without dedup (R4). |
| `duplicate_field_arc_for_test` | **Stable** | `duplicate_field_arc_for_test` |
| `layer_minpoly_coords_for_adjoin` | **Pipeline private** | `layer_minpoly_coords_for_adjoin` |
| `flatten_min_poly_over_q_cold` | **Pipeline private** | `flatten_min_poly_over_q_cold` |
| `flatten_min_poly_over_q` | **Pipeline private** | `flatten_min_poly_over_q` |
| `compose_min_poly_via_charpoly` | **Pipeline private** | compositum minpoly via multiplication-matrix charpoly (no k-search). |
| `compose_min_poly_over_q` | **Pipeline private** | `compose_min_poly_over_q` |
| `rational_subfield_embedding` | **Pipeline private** | `rational_subfield_embedding` |
| `flatten_layer_blocks` | **Pipeline private** | `flatten_layer_blocks` |
| `fields_same_parent` | **Pipeline private** | `fields_same_parent` |
| `direct_adjoin_parent_embedding` | **Pipeline private** | `direct_adjoin_parent_embedding` |
| `compose_field_embeddings` | **Pipeline private** | `compose_field_embeddings` |
| `mult_matrix_layer_gen` | **Pipeline private** | `mult_matrix_layer_gen` |
| `compute_common_dispatch` | **Pipeline private** | `compute_common_dispatch` |
| `subfield_common_pair` | **Pipeline private** | `subfield_common_pair` |
| `simple_over_q_embedding` | **Pipeline private** | `simple_over_q_embedding` |
| `embed_simple_over_q_coords` | **Pipeline private** | `embed_simple_over_q_coords` |
| `common_adjoin_sibling_over` | **Pipeline private** | `common_adjoin_sibling_over` |
| `try_common_adjoin_one_simple` | **Pipeline private** | nested + simple-over-ℚ compositum when flatten k-search fails |
| `common_primitive_sum` | **Pipeline private** | `common_primitive_sum` |
| `embed_coords` | **Stable** | `embed_coords` |
| `new` | **Pipeline private** | `new` |
| `zero` | **Pipeline private** | `zero` |
| `one` | **Pipeline private** | `one` |
| `rational` | **Pipeline private** | `rational` |
| `generator` | **Pipeline private** | `generator` |
| `add` | **Pipeline private** | `add` |
| `mul` | **Pipeline private** | `mul` |
| `embedded_by` | **Pipeline private** | `embedded_by` |
| `eq_mod` | **Pipeline private** | `eq_mod` |
| `assert_embedding_ring_hom` | **Pipeline private** | `assert_embedding_ring_hom` |
| `r6_nested_adjoin_layer_two_flatten_explicit_four` | **Pipeline private** | `r6_nested_adjoin_layer_two_flatten_explicit_four` |
| `compose_minpoly_charpoly_matches_nested_sqrt2_sqrt3` | **Pipeline private** | `compose_minpoly_charpoly_matches_nested_sqrt2_sqrt3` |
| `common_adjoin_nested_sqrt2_sqrt3_with_sqrt5` | **Pipeline private** | `common_adjoin_nested_sqrt2_sqrt3_with_sqrt5` |
| `t4a_sqrt2_sqrt3_passes` | **Pipeline private** | `t4a_sqrt2_sqrt3_passes` |
| `subfield_k1_in_k2_passes` | **Pipeline private** | `subfield_k1_in_k2_passes` |
| `adjoin_nested_sqrt2_sqrt3_with_sqrt5_passes` | **Pipeline private** | `adjoin_nested_sqrt2_sqrt3_with_sqrt5_passes` |
| `sqrt2_sqrt3_passes_verify` | **Pipeline private** | `sqrt2_sqrt3_passes_verify` |
| `compositum_reentry_returns_budget_err` | **Pipeline private** | `compositum_reentry_returns_budget_err` |
| `sqrt2_cbrt2_passes_verify` | **Pipeline private** | `sqrt2_cbrt2_passes_verify` |
| `dispatch_prefers_minimal_over_t4a` | **Pipeline private** | `dispatch_prefers_minimal_over_t4a` |
| `over_parent_sqrt2_u2_minus_alpha_degree_four` | **Pipeline private** | `over_parent_sqrt2_u2_minus_alpha_degree_four` |
| `simple_layer_minpoly_from_field` | **Pipeline private** | `simple_layer_minpoly_from_field` |
| `compositum_plan_tensor_dims_consistent` | **Pipeline private** | `compositum_plan_tensor_dims_consistent` |
| `verify_k5_flat_k2_core_embeddings_pass_v2` | **Pipeline private** | `verify_k5_flat_k2_core_embeddings_pass_v2` |
| `flatten_k2_sqrt2_sqrt3_is_fast` | **Pipeline private** | `flatten_k2_sqrt2_sqrt3_is_fast` |
| `nested_sqrt5_upstream_passes_verify` | **Pipeline private** | `nested_sqrt5_upstream_passes_verify` |
| `dispatch_nested_sqrt5_prefers_upstream` | **Pipeline private** | `dispatch_nested_sqrt5_prefers_upstream` |
| `factor_x4_minus_4_over_sqrt2_selects_quadratic` | **Pipeline private** | `factor_x4_minus_4_over_sqrt2_selects_quadratic` |
| `r6_parent_coeff_adjoin_flatten_explicit_four` | **Pipeline private** | `r6_parent_coeff_adjoin_flatten_explicit_four` |
| `rational_field_dimension_one` | **Pipeline private** | `rational_field_dimension_one` |
| `t3a_adjoin_k1_u2_minus_sqrt2_has_dimension_four` | **Pipeline private** | `t3a_adjoin_k1_u2_minus_sqrt2_has_dimension_four` |
| `rational_embedding_into_tower_uses_constant_block` | **Pipeline private** | `rational_embedding_into_tower_uses_constant_block` |
| `embedding_ring_hom_rational_to_tower` | **Pipeline private** | `embedding_ring_hom_rational_to_tower` |
| `common_rational_to_tower_embedding_uses_constant_block` | **Pipeline private** | `common_rational_to_tower_embedding_uses_constant_block` |
| `embedding_ring_hom_parent_to_child` | **Pipeline private** | `embedding_ring_hom_parent_to_child` |
| `embedding_ring_hom_composite_chain` | **Pipeline private** | `embedding_ring_hom_composite_chain` |
| `adjoin_cbrt2_generator_cubes_to_two` | **Pipeline private** | `adjoin_cbrt2_generator_cubes_to_two` |
| `adjoin_sqrt2_dimension_two` | **Pipeline private** | `adjoin_sqrt2_dimension_two` |
| `try_square_root_sqrt2_in_q_sqrt2_sqrt3` | **Pipeline private** | `try_square_root_sqrt2_in_q_sqrt2_sqrt3` |
| `try_square_root_sqrt6_in_q_sqrt2_sqrt3` | **Pipeline private** | `try_square_root_sqrt6_in_q_sqrt2_sqrt3` |
| `try_square_root_sqrt8_in_q_sqrt2` | **Pipeline private** | `try_square_root_sqrt8_in_q_sqrt2` |
| `try_square_root_sqrt2_in_q_sqrt2` | **Pipeline private** | `try_square_root_sqrt2_in_q_sqrt2` |
| `try_square_root_in_parent_coeff_tower` | **Pipeline private** | `try_square_root_in_parent_coeff_tower` |
| `try_square_root_after_align_in_common` | **Pipeline private** | `try_square_root_after_align_in_common` |
| `common_sqrt2_cbrt2_has_degree_six` | **Pipeline private** | `common_sqrt2_cbrt2_has_degree_six` |
| `common_cache_identity_is_fast` | **Pipeline private** | `common_cache_identity_is_fast` |
| `common_cache_hits_same_pair` | **Pipeline private** | `common_cache_hits_same_pair` |
| `adjoin_cache_dedup_same_minpoly` | **Pipeline private** | `adjoin_cache_dedup_same_minpoly` |
| `r1_duplicate_field_arc_subfield_embedding_parent_semantic_match` | **Pipeline private** | `r1_duplicate_field_arc_subfield_embedding_parent_semantic_match` |
| `r1_duplicate_field_arc_common_cache_semantic_key` | **Pipeline private** | `r1_duplicate_field_arc_common_cache_semantic_key` |
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
| `align_sqrt2_plus_sqrt3_reverse_order` | **Pipeline private** | `align_sqrt2_plus_sqrt3_reverse_order` |
| `common_cache_hits_after_dispatch` | **Pipeline private** | `common_cache_hits_after_dispatch` |
| `align_sqrt2_cbrt2_dim_six` | **Pipeline private** | `align_sqrt2_cbrt2_dim_six` |
| `dispatch_invariants_sqrt2_sqrt3` | **Pipeline private** | `dispatch_invariants_sqrt2_sqrt3` |

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
| `rational_rref` | **Pipeline private** | `rational_rref` |

### `field_session.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `new` | **Stable** | `FieldSession::new` |
| `fork_ambient` | **Stable (bounded)** | fork ambient with shared extension cache |
| `with_common_cache` | **Pipeline private** | run closure on common cache (reentrant with adjoin caches). |
| `clear_caches_for_test` | **Pipeline private** | `clear_caches_for_test` |
| `flatten_min_poly_over_q` | **Pipeline private** | explicit ℚ-flatten with session memo (R6). |
| `flatten_cache_len` | **Pipeline private** | `flatten_cache_len` |
| `adjoin_irreducible` | **Stable (bounded)** | adjoin irreducible |
| `adjoin_irreducible_parent_coeffs` | **Stable (bounded)** | adjoin with parent-coeff layer minpoly |
| `adjoin_parent_coeff_layer` | **Stable (bounded)** | parent-coeff adjoin on working field |
| `adjoin_cache_len` | **Pipeline private** | `adjoin_cache_len` |
| `get_or_create_base_by_min_poly` | **Pipeline private** | `get_or_create_base_by_min_poly` |
| `ambient` | **Stable** | ambient K (normalized poly coefficients) |
| `working` | **Stable** | current working L |
| `common_cache_len` | **Stable (bounded)** | common cache length |
| `compositum` | **Stable (bounded)** | compositum with session cache |
| `align_pair` | **Stable (bounded)** | align coords on session cache |
| `zero` | **Stable** | zero in L |
| `one` | **Stable** | one in L |
| `int` | **Stable** | integer constant in L |
| `half` | **Stable** | 1/2 in L |
| `lift` | **Stable** | lift coeff to working field |
| `align` | **Stable** | align pair on working field |
| `mul_formal_i` | **Stable (bounded)** | multiply by formal i |
| `sqrt_principal` | **Stable (bounded)** | principal sqrt on working field |
| `sqrt_in_field` | **Stable (bounded)** | sqrt in field or adjoin |
| `try_sqrt_in_field` | **Stable (bounded)** | try existing square root |
| `try_sqrt_in_field_shallow` | **Pipeline private** | Euler second-sqrt fast path before blind adjoin |
| `adjoin_sqrt` | **Stable (bounded)** | adjoin sqrt primitive |
| `adjoin_sqrt_new` | **Pipeline private** | blind adjoin sqrt |
| `adjoin_cbrt` | **Stable (bounded)** | adjoin cbrt primitive |
| `adjoin_primitive_cube_root_of_unity` | **Stable (bounded)** | adjoin ω for pure cubic roots |
| `checkpoint` | **Stable (bounded)** | working-field checkpoint |
| `restore` | **Pipeline private** | restore working field |
| `set_working` | **Pipeline private** | `set_working` |
| `add` | **Stable** | add with auto-align |
| `mul` | **Stable** | multiply with auto-align |
| `div` | **Stable** | divide with auto-align |
| `neg` | **Stable** | negate (no align) |
| `bump_to` | **Pipeline private** | `bump_to` |
| `rat` | **Pipeline private** | rational constant in field |
| `coeff_in_field` | **Pipeline private** | `coeff_in_field` |
| `is_negative_rational` | **Pipeline private** | negative constant in ℚ ⊂ K (for Δ<0 guard) |
| `embed_coeff_into` | **Pipeline private** | `embed_coeff_into` |
| `coeff_from_coords` | **Pipeline private** | embed coords as coeff in `field`. |
| `r2_session_common_cache_hit_on_second_common` | **Pipeline private** | `r2_session_common_cache_hit_on_second_common` |
| `restore_checkpoint_discards_later_adjoin` | **Pipeline private** | `restore_checkpoint_discards_later_adjoin` |
| `int_and_one_on_rational` | **Pipeline private** | `int_and_one_on_rational` |
| `lift_and_align_bump_working_on_k1` | **Pipeline private** | `lift_and_align_bump_working_on_k1` |
| `adjoin_sqrt_matches_sqrt_layer` | **Pipeline private** | `adjoin_sqrt_matches_sqrt_layer` |
| `adjoin_parent_coeff_quadratic_with_linear_term` | **Pipeline private** | `adjoin_parent_coeff_quadratic_with_linear_term` |
| `adjoin_quadratic_over_parent_coeff_cubic_parent` | **Pipeline private** | `adjoin_quadratic_over_parent_coeff_cubic_parent` |
| `r6_flatten_cache_hit_same_semantic_key` | **Pipeline private** | `r6_flatten_cache_hit_same_semantic_key` |
| `r6_simple_over_q_skips_flatten_cache` | **Pipeline private** | `r6_simple_over_q_skips_flatten_cache` |
| `set_working_restores_adjoin` | **Pipeline private** | `set_working_restores_adjoin` |
| `embed_coeff_into_rejects_non_subfield_no_grow` | **Pipeline private** | `embed_coeff_into_rejects_non_subfield_no_grow` |

### `galois_automorphism.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `solve_linear_system` | **Pipeline private** | solve `A·c = b` over ℚ; `A` is `rows × cols`, returns `c` or `None`. |
| `embed_parent_coords` | **Pipeline private** | embed parent coords into child (lower block). |
| `try_conjugate_by_power_basis` | **Pipeline private** | if `x = Σ c_k src^k` in `field`, return `Σ c_k dst^k`. |
| `try_apply_sigma_on_parent` | **Pipeline private** | σ(x) via conjugation in a cubic (or smaller) subfield, embed back. |
| `conjugate_map_sends` | **Pipeline private** | verify σ(src)=dst and σ(u)=v. |
| `coords_square_eq_mod` | **Pipeline private** | `coords_square_eq_mod` |
| `try_sigma_kappa` | **Pipeline private** | σ(κ) with κ²=u, σ(u)=v, σ\|_P from subfield conjugation. |
| `try_galois_sqrt_second` | **Pipeline private** | F4′: √v = σ(√u) when σ(α)=β on resolvent splitting field. |

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
| `poly_algext_from_poly` | **Stable** | lift rational polynomial for `poly_algext_roots`. |
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
| `poly_algext_from_poly_embeds_rational` | **Pipeline private** | `poly_algext_from_poly_embeds_rational` |
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
| `groebner_lex_over_qsqrt2` | **Pipeline private** | `groebner_lex_over_qsqrt2` |
| `groebner_grevlex_over_qsqrt2` | **Pipeline private** | `groebner_grevlex_over_qsqrt2` |
| `fglm_over_qsqrt2` | **Pipeline private** | `fglm_over_qsqrt2` |

### `poly_alg_factor.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `common_minpoly_var` | **Pipeline private** | `common_minpoly_var` |
| `factor_into_algext` | **Stable (bounded)** | irreducible factors over K for univariate `p`; `None` if not split. |
| `univariate_main_var` | **Pipeline private** | single main variable when `p` is univariate in it. |
| `expand_factor_pairs` | **Pipeline private** | `expand_factor_pairs` |
| `factor_into_via_algext` | **Stable** | factor `Poly` over ℚ via K[var] pipeline when univariate split exists. |
| `factor_univariate_pairs_over_k` | **Stable (bounded)** | sqff factor pairs over K for univariate `Poly` over ℚ (T2-3). |
| `factor_univariate_over_k` | **Stable** | factor `flat` in K[var] as `(factor, multiplicity)` pairs. |
| `factor_univariate_flat_over_k` | **Stable** | flat list of factors (with repetition for multiplicity). |
| `aligned_flat` | **Pipeline private** | normalize and re-wrap as [`FlatUni`]. |
| `try_as_rational_poly` | **Pipeline private** | embed `FlatUni` as ℚ `Poly` when all coeffs lie in ℚ. |
| `factor_square_free_over_k` | **Pipeline private** | square-free factor in K[var]. |
| `try_factor_by_roots` | **Pipeline private** | full linear split when all roots lie in K. |
| `try_factor_via_x_squared` | **Pipeline private** | p(x)=h(x²); factor h quadratically and lift to x²−r. |
| `is_even_only` | **Pipeline private** | `is_even_only` |
| `halve_exponents` | **Pipeline private** | `halve_exponents` |
| `linear_root_coeff` | **Pipeline private** | `linear_root_coeff` |
| `product_divides` | **Pipeline private** | `product_divides` |
| `minpoly_q_to_poly_over_field` | **Pipeline private** | embed monic `poly1` / ℚ into `PolyAlgExt` over `base`. |
| `minpoly_coords_to_poly` | **Pipeline private** | `poly1` minpoly over ℚ → `Poly` in `var`. |
| `factor_minpoly_over_field` | **Pipeline private** | `factor_minpoly_over_field` |
| `select_factor_for_common` | **Pipeline private** | `select_factor_for_common` |
| `push_factor` | **Pipeline private** | `push_factor` |
| `x_var` | **Pipeline private** | `x_var` |
| `flat` | **Pipeline private** | `flat` |
| `x_squared_minus_2` | **Pipeline private** | `x_squared_minus_2` |
| `x_fourth_minus_4` | **Pipeline private** | `x_fourth_minus_4` |
| `factor_into_algext_x_squared_minus_2` | **Pipeline private** | `factor_into_algext_x_squared_minus_2` |
| `factor_into_via_algext_from_rational_poly` | **Pipeline private** | `factor_into_via_algext_from_rational_poly` |
| `factor_x_squared_minus_2_over_k1` | **Pipeline private** | `factor_x_squared_minus_2_over_k1` |
| `factor_irreducible_deg5_witness_over_k` | **Pipeline private** | `factor_irreducible_deg5_witness_over_k` |
| `factor_quartic_reducible_over_q_yields_four_linear` | **Pipeline private** | `factor_quartic_reducible_over_q_yields_four_linear` |
| `factor_minpoly_x4_minus_4_over_sqrt2` | **Pipeline private** | `factor_minpoly_x4_minus_4_over_sqrt2` |
| `factor_x_fourth_minus_4_over_q` | **Pipeline private** | `factor_x_fourth_minus_4_over_q` |

### `poly_alg_ops.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `infer_ambient_field` | **Stable** | infer ambient **K** from polynomial coefficients (max-dimension field). |
| `normalize_algext_poly` | **Stable** | lift all coefficients into session working field **L**. |
| `align_algext_polys` | **Stable** | align two polynomials to a common working field **L** (monotone bump). |
| `ensure_common_field_for_polys` | **Stable** | bump session **L** to contain all coefficient fields in `polys`. |
| `div_rem_wrt_algext` | **Stable** | `(q, r)` with `a = q*b + r` in K[var] after coefficient alignment. |
| `quo_exact_wrt_algext` | **Stable** | exact quotient in K[var]; `Err` if remainder nonzero. |
| `monic_wrt_algext` | **Stable** | monic normalize w.r.t. `var` in K[var]. |
| `gcd_wrt_algext` | **Stable** | gcd in K[var] after coefficient alignment (monic). |
| `gcd_algext` | **Stable** | multivariate gcd in K[x₁,…,xₙ] via subresultant PRS (T3-1). |
| `egcd_wrt_algext` | **Stable** | extended gcd in K[var]; `(g, s, t)` with monic `g`. |
| `content_wrt_algext` | **Stable** | scalar content in K. |
| `primitive_part_wrt_algext` | **Stable** | primitive part in K[var]. |
| `square_free_part_wrt_algext` | **Stable** | square-free part in K[var]. |
| `split_quadratic_factor` | **Stable** | split monic quadratic into linear factors `(var − root)` over K. |
| `vars_in_algext` | **Pipeline private** | variable union for multivariate gcd |
| `flat_aligned` | **Pipeline private** | aligned [`FlatUni`] for one polynomial. |
| `aligned_pair` | **Pipeline private** | aligned pair in K[var]. |
| `monomial_to_poly` | **Pipeline private** | `monomial_to_poly` |
| `eval_univariate_algext_at` | **Pipeline private** | `eval_univariate_algext_at` |
| `x_var` | **Pipeline private** | `x_var` |
| `sqrt2_coeff` | **Pipeline private** | `sqrt2_coeff` |
| `x_squared_minus_2` | **Pipeline private** | `x_squared_minus_2` |
| `x_minus_sqrt2` | **Pipeline private** | `x_minus_sqrt2` |
| `rem_x_squared_minus_2_mod_x_minus_sqrt2` | **Pipeline private** | `rem_x_squared_minus_2_mod_x_minus_sqrt2` |
| `monic_wrt_leading_one` | **Pipeline private** | `monic_wrt_leading_one` |
| `align_scales_do_not_change_divisibility` | **Pipeline private** | `align_scales_do_not_change_divisibility` |
| `x_plus_sqrt2` | **Pipeline private** | `x_plus_sqrt2` |
| `y_var` | **Pipeline private** | `y_var` |
| `x_squared_minus_2_y_squared` | **Pipeline private** | `x_squared_minus_2_y_squared` |
| `x_minus_sqrt2_y` | **Pipeline private** | `x_minus_sqrt2_y` |
| `gcd_bivariate_x_squared_minus_2y_squared_and_x_minus_sqrt2_y` | **Pipeline private** | `gcd_bivariate_x_squared_minus_2y_squared_and_x_minus_sqrt2_y` |
| `gcd_x_squared_minus_2_and_x_minus_sqrt2` | **Pipeline private** | `gcd_x_squared_minus_2_and_x_minus_sqrt2` |
| `gcd_x_squared_minus_2_and_x_plus_sqrt2` | **Pipeline private** | `gcd_x_squared_minus_2_and_x_plus_sqrt2` |
| `gcd_x_squared_minus_2_and_x_plus_one_is_one` | **Pipeline private** | `gcd_x_squared_minus_2_and_x_plus_one_is_one` |
| `split_quadratic_x_squared_minus_2` | **Pipeline private** | `split_quadratic_x_squared_minus_2` |
| `square_free_part_x_squared_minus_2_is_self` | **Pipeline private** | `square_free_part_x_squared_minus_2_is_self` |
| `flat_uni_quo_exact_matches_div_rem` | **Pipeline private** | `flat_uni_quo_exact_matches_div_rem` |

### `poly_alg_partfrac.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `quadratic_factor_needs_k_in_q` | **Pipeline private** | irreducible quadratic over ℚ with disc > 0 (splits in K). |
| `partfrac_needs_k_split` | **Stable** | ℚ partfrac cannot split an irreducible quadratic (disc>0) → use K[var]. |
| `partfrac_sqff_factor_needs_k` | **Stable** | sqff quadratic factor irreducible in ℚ but splits in K (disc > 0). |
| `partfrac_rational_terms_over_k` | **Stable (bounded)** | `(poly_part, terms)` with denominators factored in K[var]. |
| `linear_root_coeff` | **Pipeline private** | `-c0/lc` for monic-linear `lc*x + c0`. |
| `eval_wrt` | **Pipeline private** | Horner evaluation in K. |
| `x` | **Pipeline private** | `x` |
| `partfrac_one_over_x_squared_minus_two` | **Pipeline private** | `partfrac_one_over_x_squared_minus_two` |
| `partfrac_x_over_x_squared_minus_two` | **Pipeline private** | `partfrac_x_over_x_squared_minus_two` |
| `partfrac_needs_k_split_false_for_x_squared_minus_one` | **Pipeline private** | `partfrac_needs_k_split_false_for_x_squared_minus_one` |

### `poly_alg_resultant.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `resultant_wrt_algext` | **Stable** | Sylvester resultant of univariate `a`, `b` ∈ K[var]; scalar in K. |
| `leading_coeff_wrt` | **Pipeline private** | leading coefficient in K[var]. |
| `univariate_coeffs_asc` | **Pipeline private** | ascending coeffs c₀ + c₁ var + … in K. |
| `sylvester_det2_wrt` | **Pipeline private** | Res(ax+b, cx+d) = ad − bc. |
| `sylvester_det_field` | **Pipeline private** | Sylvester matrix determinant in K. |
| `det_field` | **Pipeline private** | Gaussian elimination determinant in K. |
| `coeff_pow` | **Pipeline private** | c^exp in K. |
| `x_var` | **Pipeline private** | `x_var` |
| `session_for` | **Pipeline private** | `session_for` |
| `sqrt2_coeff` | **Pipeline private** | `sqrt2_coeff` |
| `x_squared_minus_2` | **Pipeline private** | `x_squared_minus_2` |
| `x_minus_sqrt2` | **Pipeline private** | `x_minus_sqrt2` |
| `from_q` | **Pipeline private** | `from_q` |
| `resultant_two_linear_polys` | **Pipeline private** | `resultant_two_linear_polys` |
| `resultant_shared_factor_is_zero` | **Pipeline private** | `resultant_shared_factor_is_zero` |
| `resultant_x_squared_minus_2_and_x_minus_sqrt2` | **Pipeline private** | `resultant_x_squared_minus_2_and_x_minus_sqrt2` |
| `resultant_x_squared_minus_2_and_x_plus_one` | **Pipeline private** | `resultant_x_squared_minus_2_and_x_plus_one` |
| `resultant_quadratic` | **Pipeline private** | `resultant_quadratic` |

### `poly_alg_sturm.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `sturm_sequence_wrt_algext` | **Stable** | classical Sturm chain in K[var]: P₀ = sqff(p), P₁ = P₀′, Pᵢ₊₁ = −rem(Pᵢ₋₁, Pᵢ). |
| `sturm_sign_variations_at_algext` | **Stable** | sign-variation count V(a) for a Sturm sequence at rational `a` ∈ ℚ ⊂ K. |
| `sturmab_count_wrt_algext` | **Stable** | root count in `(a, b]` for `p` ∈ K[var] (odd-multiplicity sqff convention on ℚ embed). |
| `sturmab_count_rational_poly` | **Stable** | `sturmab_count_wrt_algext` for ℚ[var] via `poly_algext_from_poly` lift. |
| `eval_algext_at_rational` | **Pipeline private** | Horner evaluation at x ∈ ℚ ⊂ L. |
| `rational_coeff` | **Pipeline private** | embed ℚ constant into current L. |
| `negate_algext_poly` | **Pipeline private** | −p in K[var]. |
| `sign_of_algext_value` | **Pipeline private** | sign for Sturm (rational or positive real constant in ℚ ⊂ K). |
| `sign_variations_i8` | **Pipeline private** | sign changes in Sturm sequence values. |
| `x_var` | **Pipeline private** | `x_var` |
| `session_for` | **Pipeline private** | `session_for` |
| `sturm_sequence_x_squared_minus_two` | **Pipeline private** | `sturm_sequence_x_squared_minus_two` |
| `sturmab_count_x_squared_minus_two` | **Pipeline private** | `sturmab_count_x_squared_minus_two` |

### `poly_conv.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `expr_contains_alg_coeff` | **Stable** | predicate for callers choosing `expr_to_poly` vs `poly_alg_from_expr`. |
| `rational_poly_reject` | **Stable** | `rational_poly_reject` |
| `first_alg_coeff_message` | **Pipeline private** | `first_alg_coeff_message` |
| `detects_rootof_in_sum` | **Pipeline private** | `detects_rootof_in_sum` |
| `detects_algext_in_product` | **Pipeline private** | `detects_algext_in_product` |
| `rational_poly_has_no_alg_coeff` | **Pipeline private** | `rational_poly_has_no_alg_coeff` |

### `poly_roots.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `prepare` | **Pipeline private** | infer K → normalize → monic; sole entry to roots_dispatch |
| `prepare_with_session` | **Pipeline private** | prepare with fresh session from coefficient field of p |
| `poly_algext_roots` | **Stable (bounded)** | exact AlgExtC roots deg 1–4; quartic resolvent gap |
| `poly_algext_roots_for_ctx` | **Stable (bounded)** | roots with Context session |
| `poly_algext_roots_in_session` | **Pipeline private** | shared roots dispatch |
| `roots_dispatch` | **Pipeline private** | degree dispatch on prepared monic input |
| `infer_field` | **Pipeline private** | infer ambient K from PolyAlgExt coefficients |
| `coeff_at` | **Stable** | univariate coefficient at exponent |
| `normalize_coeffs` | **Pipeline private** | lift all coeffs to ambient K |
| `monomial_to_poly` | **Pipeline private** | `monomial_to_poly` |
| `monic_univariate` | **Pipeline private** | `monic_univariate` |
| `linear_root` | **Pipeline private** | `linear_root` |
| `quadratic_roots_formula` | **Pipeline private** | quadratic roots via √Δ (no x²+bx+c adjoin layer) |
| `sqrt_disc` | **Pipeline private** | sqrt(Δ) via session; imaginary branch when Δ<0 in ℚ ⊂ K |
| `mul_i` | **Pipeline private** | formal i times real z on session working field |
| `cubic_depressed_parts` | **Pipeline private** | √γ = −q / (√α·√β) for depressed x⁴+px²+qx+r (q≠0) |
| `finish_cubic_roots` | **Pipeline private** | `finish_cubic_roots` |
| `f2_pure_cubic_roots_in_session` | **Pipeline private** | F2 A′: monic t³+a₀ → β·ω^k (ω adjoin once). |
| `f2_one_cubic_root_for_deflate` | **Pipeline private** | F2 deflate path: one root via irreducible t³+p_dep·t+q (no Cardano stack). |
| `f2_split_cubic_via_deflate_in_session` | **Pipeline private** | F2: z₀ + deflate quadratic; √Δ only via sqrt_in_field (F1). |
| `split_monic_cubic_roots_in_session` | **Pipeline private** | unified cubic split (P3-6 §3 / F2). |
| `resolvent_cubic_roots_in_session` | **Pipeline private** | resolvent three roots in L (split_monic_cubic_roots_in_session) |
| `cubic_roots` | **Pipeline private** | `cubic_roots` (delegates to F2 resolvent API) |
| `pure_cubic_roots` | **Pipeline private** | monic t³+a₀=0: β=∛(−a₀), roots β·ω^k |
| `root_coords_lex_key` | **Pipeline private** | lex sort roots by coords in working field (F3 canonical α,β,γ) |
| `canonical_sort_roots` | **Pipeline private** | `canonical_sort_roots` |
| `dedup_roots` | **Pipeline private** | merge roots modulo minpoly |
| `one_cubic_root` | **Pipeline private** | `one_cubic_root` |
| `casus_adjoin_cubic_root` | **Pipeline private** | casus: adjoin one root of monic t³+pt+q (Δ<0) directly |
| `embed_real_coords_for_parent` | **Pipeline private** | `embed_real_coords_for_parent` |
| `quartic_roots` | **Pipeline private** | `quartic_roots` |
| `quartic_roots_by_adjoin_deflate` | **Pipeline private** | fallback: adjoin one quartic root, then solve deflated cubic. |
| `literal_monic_univariate` | **Pipeline private** | `literal_monic_univariate` |
| `adjoin_one_root_of_monic` | **Pipeline private** | `adjoin_one_root_of_monic` |
| `biquadratic_roots` | **Pipeline private** | `biquadratic_roots` |
| `depress_quartic` | **Pipeline private** | `depress_quartic` |
| `build_resolvent_cubic` | **Pipeline private** | Ferrari resolvent R(z)=z³−pz²−4rz+(4pr−q²) on session |
| `sqrt_in_field_euler_second` | **Pipeline private** | F4′: pick β among conjugates with √(β/α) ∈ L(√α); else Galois σ; else blind adjoin. |
| `try_galois_sqrt_second_coeff` | **Pipeline private** | F4′: embed coeffs in parent P and call `galois_automorphism`. |
| `try_euler_one_alpha` | **Pipeline private** | F4′: one Euler branch with α = sorted[alpha_idx]. |
| `euler_depressed_quartic_roots` | **Pipeline private** | Euler resolvent: four roots from three resolvent zeros α,β,γ (F3 single path) |
| `euler_derived_sqrt_gamma` | **Pipeline private** | √γ = −q / (√α·√β) for depressed x⁴+px²+qx+r (q≠0) |
| `euler_four_roots_from_triple` | **Pipeline private** | `euler_four_roots_from_triple` |
| `roots_all_vanish` | **Pipeline private** | `roots_all_vanish` |
| `root_vanishes` | **Pipeline private** | verify root vanishes (fresh session, like verify_root) |
| `eval_vanishes` | **Pipeline private** | quick vanishing check on working session |
| `deflate_monic` | **Pipeline private** | `deflate_monic` |
| `real_from_pure_imag` | **Pipeline private** | real coords of z = i·r (formal i, r ∈ K real) |
| `sqrt_branches_pure_imag` | **Pipeline private** | ±√(i·r) for formal i and real r (cyclotomic biquadratic) |
| `sqrt_branches` | **Pipeline private** | ±√z via session (replaces blind `algext_square_roots`) |
| `coeff_inv` | **Pipeline private** | `coeff_inv` |
| `eq_mod` | **Pipeline private** | `eq_mod` |
| `coeff_inv` | **Stable** | `Poly::coeff_inv` |
| `eq_mod` | **Stable** | `Poly::eq_mod` |
| `verify_root` | **Pipeline private** | verify root vanishes mod minpoly (same prepare path as roots) |
| `krylov_minpoly_coords` | **Pipeline private** | `krylov_minpoly_coords` |
| `poly_mul_mod` | **Pipeline private** | `poly_mul_mod` |
| `poly_from_low` | **Pipeline private** | `poly_from_low` |
| `poly_to_low` | **Pipeline private** | `poly_to_low` |
| `poly_inv_mod` | **Pipeline private** | `poly_inv_mod` |
| `ff_mul` | **Pipeline private** | `ff_mul` |
| `ff_pow` | **Pipeline private** | `ff_pow` |
| `ff_degree` | **Pipeline private** | `ff_degree` |
| `ff_sqrt` | **Pipeline private** | `ff_sqrt` |
| `rational_reconstruct` | **Pipeline private** | `rational_reconstruct` |
| `ff_extgcd` | **Pipeline private** | `ff_extgcd` |
| `ff_inv` | **Pipeline private** | `ff_inv` |
| `crt_combine_poly` | **Pipeline private** | `crt_combine_poly` |
| `sqrt_base_case` | **Pipeline private** | `sqrt_base_case` |
| `padic_sqrt_lift` | **Pipeline private** | `padic_sqrt_lift` |
| `pn_reduce` | **Pipeline private** | `pn_reduce` |
| `pn_mulmod` | **Pipeline private** | `pn_mulmod` |
| `pn_inv` | **Pipeline private** | `pn_inv` |
| `lcm` | **Pipeline private** | `lcm` |
| `mon_x_pow` | **Pipeline private** | `mon_x_pow` |
| `ff_reduce` | **Pipeline private** | `ff_reduce` |
| `poly_to_low_bigint` | **Pipeline private** | `poly_to_low_bigint` |
| `crt_combine` | **Pipeline private** | `crt_combine` |
| `mod_inv` | **Pipeline private** | `mod_inv` |
| `small_primes` | **Pipeline private** | `small_primes` |
| `is_prime_i64` | **Pipeline private** | `is_prime_i64` |
| `mod_pow_i64` | **Pipeline private** | `mod_pow_i64` |
| `large_primes` | **Pipeline private** | `large_primes` |
| `isqrt` | **Pipeline private** | `isqrt` |
| `sqrt_fmodule` | **Pipeline private** | `sqrt_fmodule` |
| `q_session` | **Pipeline private** | `q_session` |
| `rat_coeff` | **Pipeline private** | `rat_coeff` |
| `cubic_one_root_vanishes` | **Pipeline private** | `cubic_one_root_vanishes` |
| `quadratic_sqrt_four_times_sqrt2_over_k1` | **Pipeline private** | `quadratic_sqrt_four_times_sqrt2_over_k1` |
| `quadratic_x2_minus_sqrt2_roots_vanish` | **Pipeline private** | `quadratic_x2_minus_sqrt2_roots_vanish` |
| `roots_quadratic_x2_minus_2` | **Pipeline private** | `roots_quadratic_x2_minus_2` |
| `roots_quadratic_x2_minus_sqrt2_over_k` | **Pipeline private** | `roots_quadratic_x2_minus_sqrt2_over_k` |
| `roots_cubic_t3_minus_2` | **Pipeline private** | `roots_cubic_t3_minus_2` |
| `resolvent_cubic_all_roots_vanish` | **Pipeline private** | `resolvent_cubic_all_roots_vanish` |
| `resolvent_one_cubic_root_z3_minus_4z_minus_1` | **Pipeline private** | `resolvent_one_cubic_root_z3_minus_4z_minus_1` |
| `resolvent_golden_t4_plus_t_plus_1` | **Pipeline private** | `resolvent_golden_t4_plus_t_plus_1` |
| `depressed_t4_plus_t_plus_1_coeffs` | **Pipeline private** | `depressed_t4_plus_t_plus_1_coeffs` |
| `adjoin_sqrt_squares_one_resolvent_root` | **Pipeline private** | `adjoin_sqrt_squares_one_resolvent_root` |
| `adjoin_sqrt_after_resolvent_deflate_only` | **Pipeline private** | `adjoin_sqrt_after_resolvent_deflate_only` |
| `adjoin_sqrt_after_quadratic_roots_formula` | **Pipeline private** | `adjoin_sqrt_after_quadratic_roots_formula` |
| `adjoin_sqrt_after_sqrt_disc_on_resolvent` | **Pipeline private** | `adjoin_sqrt_after_sqrt_disc_on_resolvent` |
| `resolvent_dim_bound_t4_plus_t_plus_1` | **Pipeline private** | `resolvent_dim_bound_t4_plus_t_plus_1` |
| `eval_vanishes_direct_root_dim4` | **Pipeline private** | `eval_vanishes_direct_root_dim4` |
| `euler_gamma_relation_holds` | **Pipeline private** | `euler_gamma_relation_holds` |
| `euler_four_roots_vanish` | **Pipeline private** | `euler_four_roots_vanish` |
| `field_session_dimension_bound_quartic` | **Pipeline private** | `field_session_dimension_bound_quartic` |
| `field_session_dimension_bound_quartic_tight` | **Pipeline private** | `field_session_dimension_bound_quartic_tight` |
| `quartic_adjoin_deflate_t4_plus_t_plus_1` | **Pipeline private** | `quartic_adjoin_deflate_t4_plus_t_plus_1` |
| `quartic_a4_galois_dim_le_12` | **Pipeline private** | `quartic_a4_galois_dim_le_12` |
| `roots_quartic_t4_plus_t_plus_1` | **Pipeline private** | `roots_quartic_t4_plus_t_plus_1` |
| `roots_x3_minus_x_plus_1_vanish` | **Pipeline private** | `roots_x3_minus_x_plus_1_vanish` |
| `field_session_dimension_bound_cubic_x3_minus_x_plus_1` | **Pipeline private** | `field_session_dimension_bound_cubic_x3_minus_x_plus_1` |
| `f2_resolvent_split_no_cardano_stack` | **Pipeline private** | `f2_resolvent_split_no_cardano_stack` |
| `roots_biquadratic_x4_plus_1` | **Pipeline private** | `roots_biquadratic_x4_plus_1` |
| `roots_biquadratic_t4_minus_2` | **Pipeline private** | `roots_biquadratic_t4_minus_2` |
| `diag_a4_sqrt_probe_gap` | **Pipeline private** | `diag_a4_sqrt_probe_gap` |
| `diag_adjoin_collapse_recovery` | **Pipeline private** | `diag_adjoin_collapse_recovery` |
| `diag_fmodule_step1` | **Pipeline private** | `diag_fmodule_step1` |
| `diag_fmodule_sqrt_recovery` | **Pipeline private** | `diag_fmodule_sqrt_recovery` |
| `proto_subst` | **Pipeline private** | `proto_subst` |
| `proto_rational_roots` | **Pipeline private** | `proto_rational_roots` |
| `generic_vandermonde` | **Pipeline private** | `generic_vandermonde` |
| `substitute_linear` | **Pipeline private** | `substitute_linear` |
| `proto_try_sqrt_flat_over_q` | **Pipeline private** | `proto_try_sqrt_flat_over_q` |
| `verify` | **Pipeline private** | `verify` |
| `dfs` | **Pipeline private** | `dfs` |
| `diag_norm_of_disc_c_is_square` | **Pipeline private** | `diag_norm_of_disc_c_is_square` |
| `diag_bruteforce_sqrt_of_resolvent_core` | **Pipeline private** | `diag_bruteforce_sqrt_of_resolvent_core` |
| `integer_sqrt` | **Pipeline private** | `integer_sqrt` |
