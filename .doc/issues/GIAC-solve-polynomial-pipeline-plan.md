# GIAC-solve — 全多项式方程管线计划

**状态:** open  
**类型:** 计划 / AFK 索引  
**父项:** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) §5.1（P3-7、P4-6、P4-7）  
**相关:** [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md)、[GIAC-rs-crate-dedup-plan](GIAC-rs-crate-dedup-plan.md) D3-1  
**Rust 落点:** `giac-solve::{solve,froot,rootof}`、`giac-poly::factor`  
**快照:** 2026-06-23（S0/S2/S3 + P3-6 C1–F3 绿后复审）

**P3-6 缺口索引:** [GIAC-poly-p3-6-quartic-roots-gaps](GIAC-poly-p3-6-quartic-roots-gaps.md)（含 S0 暴露的三次 **C1–C4**）

---

## 0. 目标

实现 normative **全多项式方程** solve 管线（[algext-backlog §5.1](GIAC-poly-algext-backlog.md)）：

```text
solve(P)   P ∈ Poly<ℚ>, 变元 t
  ├─ sqff / factor 降次
  ├─ 有理根 → ℚ
  └─ 不可约因子 Pᵢ
        ├─ deg(Pᵢ) ≤ 4 → poly_algext_roots
        └─ deg(Pᵢ) ≥ 5 → rootof(α, Pᵢ) 一支（无根式闭式）
```

**非目标：** 通用五次及以上根式公式（Abel–Ruffini）；`Poly<AlgExtC>` 上 Hensel/Zassenhaus 全移植。

---

## 1. 基线（已完成，勿重复立项）

| 能力 | 状态 | 证据 |
|------|------|------|
| `poly_algext_roots` deg 1–4 | ✅ | `poly_roots::tests` **25/25** |
| 四次塔 F1–F3 | ✅ | `sqrt_in_field`、split resolvent、Euler 单路径 |
| 三次 C1–C4 | ✅ | `x³−x+1` verify + solve 代入 |
| `solve(t⁴+t+1=0,t)` e2e | ✅ | `solve_quartic_t4_plus_t_plus_1` |
| **F4 维数** | ✅ | `t⁴+t+1` 实测 **dim=24** ≤ `D_HARD` |

**文档债务：** [F1–F5](GIAC-poly-quartic-roots-F1-F5.md) §0.1 进度表过期；[D3-1](GIAC-rs-crate-dedup-plan.md) 部分 AC 未勾选。

---

## 2. 缺口症状（驱动优先级）

> **2026-06-23（复审）：** P3-6 核心与 S0/S2/S3 已绿；F4 硬顶 24 ✅；F4′ σ(κ) open。

| 输入 | 现状 | 根因 / issue |
|------|------|--------------|
| `solve(x⁵-x+1=0,x)` | ✅ 1 支 `rootof` | S0 |
| `solve((x²+1)(x³-x+1)=0,x)` | ✅ 5 根 + 代入 | C1 已修复 |
| `solve(x³−x+1=0,x)` | ✅ | C1/C4 |
| `solve(t⁴+t+1=0,t)` e2e | ✅ | S2 |
| `froot` deg≤4 | ✅ | S3 |
| **F4 维数** | ✅ dim=24 | `D_HARD=24`；无未证之 12 目标 |
| `realroot` 一般代数 | ❌ | S6 |

---

## 3. 优先级总表

| 优先级 | ID | 内容 | 估时 | 解除条件 |
|--------|-----|------|------|----------|
| **P0** | **S0** | 共享内核 `solve_factor_roots` + deg≥5 `rootof` 一支 | 2–3d | — |
| **P0** | **S1** | `eval_solve` 接 factor/sqff 递归（P4-6 骨架） | 2–3d | S0 |
| **P1** | **S2** | 验收测 + 删 stale ignore / 更新 D3-1 文档 | 0.5d | S1 |
| **P1** | **S3** | `froot` / `froots` 共用 S0 内核（P4-7） | 1d | S0 |
| **P2** | **S4** | 删 `rootof.rs` 二次/双二次 fallback（P4-1） | 1d | S1 绿 |
| **P2** | **S5** | factor 覆盖加固（P3-7 子集：solve 触达路径） | 2–4d | S1 |
| **P3** | **S6** | `realroot` + Sturm 隔离（P3-4 子集） | 3–5d | S0 |
| **P3** | **S7** | `proot` / `fsolve` 与精确根互补 | 2–3d | 无 |
| **P4** | **S8** | partfrac disc>0（P2-3）、factor over K（P3-3） | 长期 | 独立竖切 |

```text
P0  S0 共享内核 + deg≥5 rootof
     ↓
P0  S1 eval_solve factor 递归
     ↓
P1  S2 验收 / 文档    P1  S3 froot 接线
     ↓
P2  S4 删 rootof 旁路   P2  S5 factor 加固
     ↓
P3  S6 realroot        P3  S7 数值互补
     ↓
P4  S8 代数系数生态（partfrac / factor over K）
```

---

## 4. 分阶段规格

### S0 — 共享因子求根内核（P0）

**落点:** `giac-solve/src/solve.rs`（或新 `solve_poly.rs`）

**API（管线私有）:**

```text
solve_univariate_over_q(poly, var, ctx) -> Vec<ExprArc>
  1. if deg=0: TypeError / []
  2. rational roots (现有 giac_poly 或 factor_into 线性因子)
  3. sqff_factorization 或 factor_into
     - 失败且 deg≥5 → 单因子 [poly]
  4. for each sqff factor F with multiplicity m:
       append solve_irreducible_factor(F, var, ctx) × m
  5. dedup（eq_mod 或字面）; return

solve_irreducible_factor(F, var, ctx):
  d = deg(F)
  if d ≤ 4: poly_algext_roots_for_ctx → Expr
  if d ≥ 5: irreducible_rootof_branch(F, var)  // 一支 rootof([1,0,...], minpoly)
```

**`irreducible_rootof_branch`:** `univariate_poly_to_poly1_expr` + 现有 `rootof_expr([1,0], minpoly)`；**禁止**走 `biquadratic_rootof_roots`。

**验收:**

- [ ] `solve(x⁵-x+1=0,x)` → 1 支 `rootof`，非 `expected quartic`
- [ ] `solve((x²+1)(x³-x+1)=0,x)` → 5 根（2+3），`assert_equiv` 代入原式
- [ ] `solve(x⁴-1=0,x)` → 4 根（factor 降次，非整式 quartic 硬算）
- [ ] 单元测 `solve_irreducible_deg5_one_branch`

**不做:** Galois 可解五次特判；`poly_algext_roots` deg>4 闭式。

---

### S1 — `eval_solve` 统一入口（P0）

**改动:** `poly_roots_as_exprs` 改为调用 S0；删除「整式 deg 直撞 fallback」路径。

```text
eval_solve:
  equation_to_poly → solve_univariate_over_q
```

保留 `giac_poly::roots` 仅作 **deg≤1 有理根快路径**（可选）；deg≥2 统一进 S0。

**验收:** S0 全套 + 既有 `solve_biquadratic` / `solve_quadratic` 仍绿。

---

### S2 — 验收与文档清扫（P1）

| 任务 | 说明 |
|------|------|
| 删 `solve_quartic_t4_plus_t_plus_1` 的 `#[ignore]` | roots 层已绿 |
| 增 `solve_reducible_deg5_product` | `(x²+1)(x³-x+1)` |
| 增 `solve_irreducible_deg5` | `x⁵-x+1` |
| 更新 [F1–F5](GIAC-poly-quartic-roots-F1-F5.md) §0.1 → M3 达成 | |
| 勾选 [D3-1](GIAC-rs-crate-dedup-plan.md) AC | |

---

### S3 — `froot` 共用内核（P1, P4-7）

**改动:** `froot.rs::solve_factor_roots` 改调 `solve_irreducible_factor`（带重数），deg≤4 不再 `NotImplemented`。

**验收:** `froot((x²+1)(x-1))` 给出 3 个有理/代数根对。

---

### S4 — 删 `rootof.rs` 形状表（P2, P4-1）

**范围:**

| 删/收窄 | 保留 |
|---------|------|
| `try_algext_or_rootof_roots` 末尾 `quadratic_rootof` / `biquadratic` fallback | `irreducible_rootof_branch`（deg≥5） |
| `biquadratic_rootof_roots` 生产路径 | `rootof_expr` 构造器（S0 用） |

**验收:** `giac-solve` 全测绿；`rootof.rs` 仅 deg≥5 + 测试对照。

---

### S5 — factor 触达路径加固（P2, P3-7 子集）

仅加固 **solve 竖切会碰到的** `factor_into` 失败情形；不全移植 Zassenhaus over `AlgExt`。

| 加固 | 例 |
|------|-----|
| 可约五次显式分解 | `(x-1)(x⁴+x³+…)` |
| `x^n ± 1` / cyclotomic 已有路径确认 | `x⁶-1` |
| factor 失败 → 当作不可约单因子（S0 兜底） | 不 panic |

**验收:** solve 不因 factor 失败而误报 `expected quartic`。

---

### S6 — `realroot` 代数实根（P3, P3-4 子集）

**路径:** Sturm（已有 `sturm.rs`）→ 区间隔离 → 区间端点 `rootof` 或有理端点。

**验收:** `realroot(x²-2)` 两实根；`realroot(x³-x-1)` 一支实区间（可选 `#[ignore]` 精度）。

**阻塞:** 无 S0；可与 S4 并行。

---

### S7 — 数值互补（P3）

| API | 最小交付 |
|-----|----------|
| `fsolve` | 文档化边界；代数系数式子返回明确 `Err` 或数值近似 |
| `proot` | 委托 `fsolve` 或区间 Newton；登记 `known-divergences` |

**不替代** S0 精确 `rootof`。

---

### S8 — 代数系数生态（P4，独立）

|  backlog ID | 内容 |
|-------------|------|
| P2-3 | `partfrac(1/(x²-2))` disc>0 分裂 |
| P3-1/2/3 | gcd、sqff、factor over `Poly<AlgExtC>` |
| P3-5 | 嵌套环 `Poly<AlgExtC>[main]` |

与 solve 管线 **解耦**；S0–S4 完成后按需竖切。

---

## 5. 里程碑

- [ ] **M1 — P0 闭环:** S0+S1；`x⁵-x+1`、`(x²+1)(x³-x+1)` CLI 绿
- [ ] **M2 — P1 清扫:** S2+S3；`froot` deg≤4；无 stale ignore
- [ ] **M3 — P2 去旁路:** S4+S5；`rootof.rs` 无二次/双二次生产路径
- [ ] **M4 — P3 实根/数值:** S6 可选 + S7 登记
- [ ] **M5 — P4 生态:** S8 按 partfrac/integrate 需求立项

---

## 6. AFK 抓取顺序

| 顺序 | PR 范围 | 估时 | 备注 |
|------|---------|------|------|
| 1 | **S0** 内核 + 3 个红灯测 | 2–3d | 先 `solve_irreducible_deg5` 红灯 |
| 2 | **S1** `eval_solve` 接线 | 1–2d | 同 PR 或紧随 |
| 3 | **S2** unignore + 文档 | 0.5d | |
| 4 | **S3** froot | 1d | |
| 5 | **S4** 删 fallback | 1d | 单独 PR 便于 review |
| 6 | **S5** factor 加固 | 按需 | 遇红再加 |

**禁止:**

- 在 S0 前为 `x⁵-x+1` 加 `rootof.rs` 新形状表
- 在 S1 前 unignore 却不接 factor 递归（可约高次仍挂）
- 实现通用五次根式闭式

---

## 7. 测试矩阵

| 测试名 | 阶段 | 输入 | 期望 |
|--------|------|------|------|
| `solve_irreducible_deg5` | S0 | `x⁵-x+1=0` | 1× `rootof` |
| `solve_reducible_deg5_product` | S0 | `(x²+1)(x³-x+1)=0` | 5 根，`assert_equiv` |
| `solve_quartic_t4_plus_t_plus_1` | S2 | `t⁴+t+1=0` | 4 根，**无 ignore** |
| `solve_x4_minus_1` | S1 | `x⁴-1=0` | 4 根（factor 路径） |
| `froot_mixed_degrees` | S3 | `(x-1)(x²-2)` | 3 根对 |
| `solve_quadratic_double_root` | 回归 | 既有 | 仍绿 |

---

## 8. 与 upstream giac 的差异

**原则（[rust-migration-plan](../rust-migration-plan.md) §6.4）：** 绑定数学语义 + 测试规格，**不绑定** `solve.cc` 行级算法与 `gen` 旁路表。

### 8.1 管线骨架 — 一致

upstream `solve_cleaned`（`solve.cc` ~2437）与 `froot`/`addfactors`（`misc.cc`）的核心模式相同：

```text
factor / sqff → 对每个因子递归 solve → 合并根
```

本计划 S0/S1 是在 **补齐 Rust 侧尚未接线的这一骨架**；不是发明新数学。

### 8.2 关键差异

| 维度 | upstream giac | 本计划（giac-rs） |
|------|---------------|-------------------|
| **代码形态** | 单文件 `solve.cc` + `gen` 全局递归 | `giac-poly` factor + `giac-core` roots + `giac-solve` 接线 |
| **deg≤4 不可约** | 多路径：`in_solve` 三次 ALT/resultant、四次 factor+`rootof`、二次 `algebraic_EXTension` | **统一** `poly_algext_roots` + `FieldSession`（Cardano/Ferrari/Euler） |
| **deg≥5 不可约** | `solve(vecteur)` 返回 **一支** `algebraic_EXTension`（~L1242） | **一致**：一支 `rootof(α, minpoly)` |
| **数值兜底** | `has_num_coeff` → `proot`/`realroot`/`evalf` 早且多 | S6/S7 后置；精确根优先 |
| **超越/含参** | `abs`、`tan`、分式幂、复合式、`fsolve`、RUR 等大量特判 | **本计划不覆盖**；`eval_solve` 仅多项式 + 极简 `sin(x)=0` |
| **可解五次特判** | 分圆 `is_cyclotomic`、复合 `translate_gcddeg` 等 | 标为 P4 可选，首版不做 |
| **输出形态** | `gen`/`rootof`/`_EXT` 混用；`complex_mode` 影响分支 | 统一 `AlgExtC` / `Expr::AlgExt`；`Context` 管 assume |

### 8.3 我们比 upstream **更强** 之处

- **deg≤4**：显式四根（`poly_algext_roots`），非 dense `solve(vecteur)` 的一支扩展。
- **塔语义**：`FieldSession` + `sqrt_in_field` 可审计维数；upstream `_EXT` common 隐式且难测。
- **类型分层**：`Poly<ℚ>` / `Poly<AlgExtC>` / `FieldSession` 分离，便于 partfrac、integrate 复用。

### 8.4 我们比 upstream **更弱**（已知、可接受）

- 无 RUR / 区间算术实根隔离（upstream `solve` RUR 分支 ~8240+）。
- 无 `modsolve`、不等式、`assume` 分支过滤。
- 数值 `proot`/`fsolve` 未与精确根闭环。

偏离登记：字面 golden 不同时走 `assert_equiv` + [known-divergences](../known-divergences.md)。

---

## 9. 共享工具层（为其他 CAS 模块预留）

S0 **不应**只活在 `giac-solve`；应落成 **`giac-core::algebra` 稳定（bounded）API**，供多模块共用。

### 9.1 推荐落点与 tier

| API | 建议 crate | Tier | 消费者 |
|-----|------------|------|--------|
| `univariate_factor_roots_over_q(poly, var, ctx)` | `giac-core::algebra::poly_roots` | **Stable (bounded)** | `solve`, `froot`, `roots` |
| `irreducible_factor_roots(f, var, session)` | 同上 | **Stable (bounded)** | 上项内部；partfrac、integrate |
| `irreducible_rootof_branch(f, var)` | `giac-core::algebra` | **Stable (bounded)** | deg≥5；`eval_rootof` 文档化 |
| `poly_algext_roots_for_ctx` | 已有 | **Stable (bounded)** | 同上 |
| `FieldSession` | 已有 | **Stable** | roots、将来 factor/gcd over K |
| `factor_into` / `square_free_factorization` | `giac-poly` | **Stable (bounded)** | solve 只 **调用**，不复制 |

`giac-solve::eval_solve` 仅做：`equation_to_poly` → `univariate_factor_roots_over_q` → `List<Expr>`。

### 9.2 未来模块接线（S0 完成后即解锁）

```text
univariate_factor_roots_over_q
  ├─ giac-solve     solve / froot / roots
  ├─ giac-calculus  Rothstein–Trager（特征根）、partfrac 分母零点
  ├─ giac-simplify  factor 后有理根探测（子集）
  ├─ giac-linalg    egv / charpoly deg≤4
  └─ giac-core      eval_partfrac（P2-3 分裂二次因子）
```

**纪律（[algorithm-expr-api](../algorithm-expr-api.md) §6）：**

1. 禁止在 `partfrac.rs` / `integrate.rs` 再写 `rootof` 形状表 — 调 `irreducible_factor_roots`。
2. 扩域计算带 `FieldSession`（或 `poly_algext_roots_for_ctx`），禁止裸 `Poly<AlgExtC>` + 散落 `align`。
3. deg≥5 统一一支 `rootof`；要小数走 `fsolve`/`proot`（S7），不另开根式旁路。

### 9.3 S0 实现约束（预留扩展）

- 返回 `Vec<ExprArc>`（用户层）+ 可选 `Vec<AlgExtCPolyCoeff>`（算法层）薄封装，避免 integrate 再 `Expr` 往返。
- `session` 从 `Context` fork，便于多根列表共享 `common_cache`（upstream 无此显式契约）。
- factor 失败 → 单因子 `[poly]` 兜底（与 upstream 数值 fallback 前行为一致）。

---

## 10. 索引

| 主题 | 文档 |
|------|------|
| 本计划 | [GIAC-solve-polynomial-pipeline-plan](GIAC-solve-polynomial-pipeline-plan.md) |
| AlgExt backlog | [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) §5.1 |
| 四次塔（已完成） | [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) |
| dedup D3-1 | [GIAC-rs-crate-dedup-plan](GIAC-rs-crate-dedup-plan.md) |
| 表示层 API | [algorithm-expr-api](../algorithm-expr-api.md) §6.3 |
| upstream 参考 | `giac-2.0.0/src/solve.cc`、`misc.cc` `addfactors` |
