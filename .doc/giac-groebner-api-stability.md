# giac-groebner API 稳定性分层

**规范来源:** [algorithm-expr-api.md](algorithm-expr-api.md)  
**表示层:** 全程 `giac_poly::Poly` / `PolyMod`；变量序显式传入  
**代码:** `giac-rs/crates/giac-groebner/src/lib.rs`

---

## 1. 源码注释格式

| 标记 | 用于 | 可见性 |
|------|------|--------|
| `/// **Stable (bounded)** — …` | 约化入口 | `pub` |

---

## 2. Crate 公开 API

| 函数 | 层级 | 说明 |
|------|------|------|
| `greduce` | **Stable (bounded)** | 多元多项式对 Gröbner 基约化（lex，`vars` 显式） |
| `greduce_grevlex` | **Stable (bounded)** | 多元多项式对 Gröbner 基约化（grevlex，`vars` 显式） |
| `greduce_mod` | **Stable (bounded)** | 模 `p` 约化 |
| `groebner_basis_lex` | **Stable (bounded)** | lex Gröbner 基（Buchberger，泛型 `C: FieldCoeff`） |
| `groebner_basis_grevlex` | **Stable (bounded)** | grevlex Gröbner 基（Buchberger，0-dim 快速，泛型 `C: FieldCoeff`） |
| `fglm` | **Stable (bounded)** | FGLM grevlex→lex 转换（0-dim，泛型 `C: FieldCoeff`） |

**未实现（规划）:** F4/F5 选择策略；当前 GB 构造为 Buchberger（lex/grevlex）。泛型入口均 `C: FieldCoeff`，ℚ 经包装兼容。

---

## 3. I/O 契约 — `greduce`

| 字段 | 说明 |
|------|------|
| **输入** | `poly: &Poly`, `basis: &[Poly]`, `vars: &[Var]` |
| **输出** | 约化 remainder `Poly` |
| **上下文** | `vars` 定义 lex 序；非 Gröbner 基时结果非规范 normal form |
| **禁止** | 假设 `Poly` 树隐含变量序 |

---

## 4. 维护

```bash
cd giac-rs
python3 scripts/annotate_api_tiers.py
python3 scripts/annotate_api_tiers.py --inventory
```

---

## Per-file function inventory (generated)

Regenerate: `python3 scripts/annotate_api_tiers.py --inventory`

### `lib.rs`

| Function | Tier | Description |
|----------|------|-------------|
| `pmul` | **Pipeline private** | `pmul` |
| `psub` | **Pipeline private** | `psub` |
| `padd` | **Pipeline private** | `padd` |
| `fadd` | **Pipeline private** | `fadd` |
| `fsub` | **Pipeline private** | `fsub` |
| `fmul` | **Pipeline private** | `fmul` |
| `fdiv` | **Pipeline private** | `fdiv` |
| `lt` | **Pipeline private** | leading-term dispatch |
| `greduce` | **Stable (bounded)** | multivariate reduce mod ideal (lex) |
| `greduce_grevlex` | **Stable (bounded)** | multivariate reduce mod ideal (grevlex) |
| `greduce_order` | **Pipeline private** | `greduce_order` |
| `greduce_mod` | **Stable (bounded)** | reduce PolyMod mod basis |
| `make_monic_order` | **Pipeline private** | `make_monic_order` |
| `spoly_order` | **Pipeline private** | `spoly_order` |
| `autoreduce_order` | **Pipeline private** | `autoreduce_order` |
| `buchberger_order` | **Pipeline private** | `buchberger_order` |
| `groebner_basis_lex` | **Stable (bounded)** | lex Gröbner basis over a field via Buchberger |
| `groebner_basis_grevlex` | **Stable (bounded)** | grevlex Gröbner basis over a field via Buchberger |
| `monom_pow` | **Pipeline private** | FGLM helper |
| `monomials_of_degree` | **Pipeline private** | FGLM helper |
| `rec` | **Pipeline private** | `rec` |
| `fglm_solve` | **Pipeline private** | FGLM linear solve |
| `fglm_matvec` | **Pipeline private** | FGLM helper |
| `fglm` | **Stable (bounded)** | FGLM grevlex→lex over a field |
| `half` | **Pipeline private** | `half` |
| `greduce_xy_minus_1` | **Pipeline private** | `greduce_xy_minus_1` |
| `greduce_circle` | **Pipeline private** | `greduce_circle` |
| `greduce_mod_linear` | **Pipeline private** | `greduce_mod_linear` |
| `groebner_lex_x2_1_y_minus_x` | **Pipeline private** | `groebner_lex_x2_1_y_minus_x` |
| `groebner_lex_unit_ideal` | **Pipeline private** | `groebner_lex_unit_ideal` |
| `groebner_grevlex_unit_inconsistent_under_100ms` | **Pipeline private** | `groebner_grevlex_unit_inconsistent_under_100ms` |
| `groebner_grevlex_x2_2_y2_3_d4` | **Pipeline private** | `groebner_grevlex_x2_2_y2_3_d4` |
| `gb_set_eq` | **Pipeline private** | `gb_set_eq` |
| `fglm_x2_1_y_minus_x_matches_lex` | **Pipeline private** | `fglm_x2_1_y_minus_x_matches_lex` |
| `fglm_x2_plus_y2_5_xy_2_d4_generic` | **Pipeline private** | `fglm_x2_plus_y2_5_xy_2_d4_generic` |
| `fglm_unit_ideal_returns_none` | **Pipeline private** | `fglm_unit_ideal_returns_none` |
