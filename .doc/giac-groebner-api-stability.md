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
| `greduce_mod` | **Stable (bounded)** | 模 `p` 约化 |

**未实现（规划）:** 完整 Gröbner 基构造（F4/F5/Buchberger）；本 crate 当前仅 **约化** 步骤。

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
| `greduce` | **Stable (bounded)** | multivariate reduce mod ideal (lex) |
| `greduce_mod` | **Stable (bounded)** | reduce PolyMod mod basis |
| `half` | **Pipeline private** | `half` |
| `greduce_xy_minus_1` | **Pipeline private** | `greduce_xy_minus_1` |
| `greduce_circle` | **Pipeline private** | `greduce_circle` |
| `greduce_mod_linear` | **Pipeline private** | `greduce_mod_linear` |
