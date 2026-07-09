# HNF 单位 arch 列：`exp=0` → `vecslice(C,1,zc)`

**状态:** implemented（`hnf_spec.rs`，与 Pari `hnffinal` 同路）  
**上游:** Pari/GP `basemath/hnf_snf.c` + `basemath/buch2.c`  
**源码树:** `/home/kanli.hu/upstream/pari/src/`（`gp` 为 `/home/kanli.hu/upstream/pari/gp` → 运行时二进制；算法以 `src/` 为准）  
**Rust:** `giac-rs/crates/giac-core/src/algebra/hnf_spec.rs`、`class_group.rs`

---

## 数学对象

| 符号 | 含义 |
|------|------|
| **IdealFactorizationRelation** `(a, γ)` | `∏ 𝔭_j^{a_j} = (γ)`，`a` 不全为 0；贡献类格 + arch |
| **UnitArchRelation** `(0, ε)` | `ε ∈ O_K^×`，`a ≡ 0`；理想因子分解平凡，**只贡献** `fixarch(ε)` |
| **`zc`** | `|C| − |W| − |B|`；`vecslice(C,1,zc)` 为单位格 arch 列 |
| **虚假 `zc`** | `zc>0` 但单位槽 arch 全零 → `pari_need=0` 不可信 |

代数单位满足 `(ε) = O_K`，故在因子基上 **`v_𝔭(ε)=0` ∀𝔭** → 指数向量 **恒为 0**。这不是「无关系」，而是 **arch-only 单位关系**。

---

## Pari 上游行为（对照）

### 1. `add_rel_i` — 仍收录 `exp=0`（`buch2.c:2284`，`nz==n+1` 于 `:2289`）

```c
if (nz == n+1) { k = 0; goto ADD_REL; }  /* nz = 首个非零指标；n+1 ⇒ 全零 */
```

`nz == KC+1` 时跳过 mod-p 秩增长，**仍写入** `REL_t`（`rel->m = ε`，`rel->emb` 后续填 arch）。

### 2. `hnfadd_i` — `C` 矩阵布局（`hnf_snf.c:612`）

```
co = |C|, lH = |W|, lB = |B|, col = co - lB
C = [ unit slots: cols 1..col-lH | W arch: cols col-lH+1..co ]
Cnew = concat(extraC, vecslice(C, col-lH+1, co))
hnffinal(...) → 重排
*ptC = concat(vecslice(C,1,col-lH), Cnew)
```

新 arch 列先与 **W 尾**拼接，再经 `hnffinal`；**旧单位槽** `vecslice(C,1,col-lH)` 原样保留。

### 3. `hnffinal` — 零列 → 单位槽（`hnf_snf.c:112-116`）

```c
H = ZM_hnflll(matgen, &U, 0);
H += lg(H)-1 - lnz; H[0] = evaltyp(t_MAT) | _evallg(lnz+1);
zc = col - lnz; /* # of 0 columns, correspond to units */
```

- `ZM_hnflll(..., remove=0)` 保留前导零列；指针前移只留尾部 `lnz×lnz`
- `zc = col - lnz`；`C*U` 后前 `zc` 列为单位 arch
- 零整数列经 `ZM_rowrankprofile` 归为 dependent（非 pivot），对应 arch 进 `zc`（非 W）

### 4. `matbotid` / `matbotidembs`（`buch2.c:3887`，`flag` 代数模式）

代数 `bnfinit` 用 **嵌入单位阵块**（非裸 `fixarch`）标记关系行；浮点模式用 `get_embs` → `rel_embed` → `fixarch`。  
giac-rs GRH 路径对标 **浮点 + fixarch**（`relation_fixarch_column`），与 `flag=0` 一侧一致。

---

## Rust 实现

**根因（旧）：** `hnffinal` 漏做 Pari 的 `H += lg(H)-1-lnz` 切片 → 对角/`C` 重排错位；零列 arch 落在 W 尾。曾用 `hnf_merge_unit_arch_columns` 显式填 `zc`（DIV-102）。

**现行为（与 Pari 同路）：**

1. `hnfadd_i` / `hnfspec`：**不拆** `exp=0`；全列进 `rowrankprofile` + `hnffinal`
2. `hnffinal`：`zm_hnflll` 后若 `ncol > lnz`，只取**末尾** `lnz` 列；`zc = col - lnz`
3. 证书仍用 `vecslice(C,1,zc)` → `compute_R`；**禁止**用 W 列冒充 `Ar`
4. **GRH grow**（`GrhRelCache`）：`matbotid` / `matbotidembs` + `c_lift_embs`（对标 Pari `bnfinit(,1)`）

```text
hnfadd_i / hnfspec
  └─ rowrankprofile + hnffinal（含 H 前导零列切片）→ C[0..zc)
```

**不再使用** `hnf_merge_unit_arch_columns` / `class_group` 事后 patch。

---

## 关系生成侧（`class_group.rs`）

| 函数 | 作用 |
|------|------|
| `add_unit_arch_relation` | `(0, ε)` 绕过 mod-p；独立 `seen_unit_gamma` |
| `enumerate_algebraic_unit_candidates` | `ℤ[α]` 有界枚举 `|N|=1`，按 arch 范数排序 |
| `enumerate_relations_unit_oriented` | FB 分解单位 + arch-only 单位 |
| `spurious_unit_zc` / `effective_pari_need` | 虚假满秩侦测 → 强制 `RELAT_BURST` |

---

## 测试锚

```bash
cargo test -p giac-core --release --lib hnfadd_zero_exp
cargo test -p giac-core --release --lib x3_11_unit_arch
cargo test -p giac-core --release --lib class_number_general_cert_grh_x3_11_h2_probe
```

- `hnfadd_zero_exp_arch_fills_unit_slot`：类格满 + `hnfadd(0,ε)` → `C[0]` 非零  
- `x3_11_unit_arch_relation_fills_hnf_unit_slot`：ℚ(∛11) 诊断批 + `ε=[-2,4,1]`  
- `class_number_general_cert_grh_x3_11_h2_probe`：端到端 GRH 证书

---

## 后续（可选）

- 类型化 `BuchmannRelation::{Ideal, UnitArch}` 替代 `(Vec<BigInt>, HighFirstQ)` 隐式约定  
- 精度提升时只重算 `embs` 顶行、不重跑 HNF（Pari `flag=1` + `myprecdbl` 循环；`hnf_prec_snaps` 仍走全量 `fixarch` stub）
