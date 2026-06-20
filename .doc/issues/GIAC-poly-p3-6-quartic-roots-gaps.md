# GIAC — P3-6 通用四次 `Poly<AlgExtC>::roots` 缺口清单

**状态:** open  
**类型:** 索引 / AFK 竖切  
**父项:** [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) §6 **P3-6**  
**相关:** [GIAC-lazy-common-tower-plan](GIAC-lazy-common-tower-plan.md) §Phase T3+、§11；[GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md) §D5；[giac-tower-common-math.md](../giac-tower-common-math.md)  
**上游参考:** `giac/giac-1.5.0/src/gausspol.cc`（`_EXT` 系数 `roots` / resolvent）、`alg_ext.cc`  
**Rust 落点:** `giac-core::algebra::poly`（`PolyAlgExt`）、`giac-poly`（泛型 `roots`）；**非** `giac-solve/rootof.rs` 新特判  
**快照:** 2026-06-20

---

## 1. 目标（P3-6 是什么）

在基域 **K = ℚ**（及 solve 会话内已登记的扩域）上，对 **`Poly<AlgExtCPolyCoeff>`** 一元四次多项式求**全部代数根**：

| 算法骨架 | 说明 |
|----------|------|
| **Resolvent cubic** | 由 \(P(t)\in K[t]\)（\(\deg P=4\)）构造三次 resolvent \(R(u)\in K[u]\) |
| **P2-6 内层** | 在 K 上解 \(R(u)=0\)（至多 3 个 \(u\in\overline K\)） |
| **T3+ adjoin** | 对每个 resolvent 根 \(u\)，在 K 上 **register** \(u^2-\alpha\) 层（\(\alpha\in K\) 非常数 embed） |
| **K 上二次 split** | 在 adjoin 后的域上解两个二次因子 → 共 4 个 \(t\) 根 |
| **输出** | `Vec<AlgExtCPolyCoeff>`（或经 `algext_poly_to_expr` 的 `List<Expr>`） |

**验收（normative）：** `solve(t⁴+t+1=0, t)` 四根，两两 `eq_mod` 对齐；**不得** `NotImplemented("general quartic rootof")`。

**明确不做：** Ferrari 在 `giac-solve` 旁路；通用五次根式；Hensel / unitaryfactor over `Poly<AlgExtC>`。

---

## 2. 现状快照

| 组件 | 状态 | 说明 |
|------|------|------|
| `Poly<AlgExtC>` 表示 + `PolyCoeff` | ✅ P1-1/2/4 | `poly_alg_coeff.rs`、`poly.rs` 桥接 |
| `AlgExtC` 同域 ±×、`inv`、`common` | ✅ P0-A/B + T4a | 塔默认 `tower-common` |
| T3 parent 系数 `element_*` | ✅ | T1b \(u^2-3\) 等；**不含** parent 系数层 **登记** |
| **T3+** `adjoin(K, u²−α)` | ❌ **硬阻塞** | `adjoin_irreducible` 仍拒 `layer_minpoly_rational_constants`（`ext_tower.rs` ~164） |
| `giac_poly::roots` | ❌ deg≤2（ℚ） | `resultant.rs`；无 `Poly<C>` 泛化 |
| **`Poly<AlgExtC>::roots`** | ❌ **未实现** | 无二次/三次/四次 AlgExt 路径 |
| P2-1 二次 roots over K | ❌ | 仍靠 `giac-solve/rootof.rs` `quadratic_rootof_roots` |
| P2-2 双二次 | ⚠️ 过渡态 | `biquadratic_rootof_roots`；**非** `Poly<AlgExtC>::roots` |
| P2-6 三次 roots over K | ❌ | resolvent 内层缺失 |
| P3-1 gcd/quo/rem over K | ❌ | resolvent / split 可能需要；至少需可测 remainder |
| P3-3 factor over K | ❌ | 可约四次应先降次（P3-7）；一般式仍要 P3-6 |
| 一般四次 `rootof` | ❌ | `rootof.rs:26` → `NotImplemented("general quartic rootof")` |
| dense poly1（D1–D4） | ✅ | T3 `ParentBlockRing` 已收敛；**不挡** P3-6 主路径（sparse + `PolyCoeff`） |
| D5 `FlatUni<AlgExtC>` ↔ dense | ☐ 可选 | 见 §5.3 |

---

## 3. 缺口清单（按阻塞级别）

### 3.1 硬依赖 — 无则 P3-6 无法验收

| ID | 缺口 | 为何需要 | 落点 / 参考 |
|----|------|----------|-------------|
| **G0** | **T3+** `adjoin(K, u²−α)` | 一般四次 resolvent 步 minpoly 系数 \(\in K\setminus\mathbb{Q}\)（如 \(u^2-\sqrt2\)）；现无法 **register** 新层 | [tower-plan §Phase T3+](GIAC-lazy-common-tower-plan.md)；`ext_tower.rs` |
| **G1** | **`Poly<AlgExtC>::roots` 三次（P2-6）** | Resolvent cubic 是 P3-6 **内层**；缺则只能 stub 四次 | 新模块 `giac-poly` 或 `giac-core::algebra::poly_roots`；对标 `gausspol.cc` `_EXT` |
| **G2** | **`Poly<AlgExtC>::roots` 二次（P2-1）** | Resolvent 后 K-adjoin 上的二次 split | 同上；替代 `quadratic_rootof_roots` |
| **G3** | **泛型 `roots(p: &Poly<C>, var)` 或 `PolyAlgExt::roots`** | 入口 API；`C=AlgExtCPolyCoeff` 分发 deg 1/2/3/4 | `giac-poly` + `giac-core` 接线 |
| **G4** | **四次 resolvent 构造 + split 编排** | P3-6 本体：\(P\to R\to\) roots\(_3\)\(\to\) adjoin \(\to\) roots\(_2\) \(\times 2\) | 新代码；数学见 §4 |
| **G5** | **系数域 K 与 `AlgExtCData.field` 一致** | resolvent 系数、adjoin 层须在同一会话塔顶；多根 `eq_mod` 不依赖 embed 顺序 | S0 align + T2 子域嵌入；**算法设计：** [`algorithm-expr-api.md`](../algorithm-expr-api.md) §6.3 `FieldSession`（ambient K + working L），禁止裸 `Poly` 树隐式承载分裂域 |

### 3.2 软依赖 — 可竖切内先做 stub，验收前须补齐

| ID | 缺口 | 说明 | 可延期到 |
|----|------|------|----------|
| **G6** | P3-1 gcd / rem over `Poly<AlgExtC>` | 子结果式 gcd 或 pseudo-remainder；split 前约化 | P3-6 PR 内最小 subset |
| **G7** | P3-2 sqff / content / primitive | 可约四次、重根；避免 resolvent 在非 sqff 上失效 | 与 P3-7 共用 |
| **G8** | P3-3 / P3-7 factor 降次 | \((t^2+1)(t^2+2)\) 等**不应**撞四次未实现 | P4-6 solve 管线；P3-6 可先只处理 **不可约** 四次 + 文档登记 |
| **G9** | P2-2 双二次迁入 `Poly<AlgExtC>::roots` | `t⁴−2` 可走特判或 resolvent；与 P4-1 删 `biquadratic_rootof_roots` 对齐 | P4-1 |
| **G10** | P4-6 solve 统一管线 | `solve` 调 `Poly<AlgExtC>::roots` 而非 `rootof.rs` 形状表 | P3-6 算法绿后再接 |

### 3.3 可选优化 — 不挡 P3-6 立项

| ID | 优化 | 收益 | 参考 |
|----|------|------|------|
| **O1** | D5 `FlatUni<AlgExtCPolyCoeff>` ↔ dense | resolvent 若走 parent 系数稠密向量，复用 `giac-poly::dense` | [dense-poly1-refactor §D5](GIAC-dense-poly1-refactor.md) |
| **O2** | resolvent 在 T3 已登记塔顶用 `ParentBlockRing` | 少 flatten；与 T1b 路径一致 | `field_arith.rs` + dense |
| **O3** | 可根式四次特判（biquadratic、\(t^4+a\)） | 快路径；**不能**替代一般式验收 | 现有 `biquadratic_rootof_roots` 迁移 |

---

## 4. 算法竖切（建议 PR 顺序）

```text
PR-A  T3+ adjoin(K, u²−α) + T3+a 验收
        ↓
PR-B  Poly<AlgExtC>::roots 二次 (P2-1) + 单测 roots(x²−2)
        ↓
PR-C  Poly<AlgExtC>::roots 三次 (P2-6) + solve(t³−2) 预备
        ↓
PR-D  resolvent cubic 构造（Poly<AlgExtC> 系数算术）
        ↓
PR-E  四次 split + PolyAlgExt::roots 入口 + solve(t⁴+t+1) 验收
        ↓
PR-F  P4-6 接线；删 rootof 一般四次 NotImplemented；P3-7 可约降次
```

**PR-D/E 核心步骤（实现 checklist）：**

1. [ ] 输入：monic 或 primitive 四次 `P` over K（`PolyAlgExt` + 变元 `t`）
2. [ ] 构造 resolvent `R(u) ∈ K[u]`（系数为 `AlgExtCPolyCoeff`；公式与 upstream 对齐或文献登记）
3. [ ] `roots_3 = R.roots()` → 0–3 个 `u`（P2-6）
4. [ ] 对每个 `u`：T3+ `adjoin(K, u²−α)`（\(\alpha\) 由 resolvent 分支确定）
5. [ ] 在 adjoin 域上解两个二次 → 合并为 4 个根（dedup / align）
6. [ ] 输出经 `eq_mod` 与手算 / golden 对齐

---

## 5. 与相邻工作的边界

### 5.1 与 T3+

- **T3 alone 不够：** 可在塔顶对 parent 系数做 `element_mul`，但 **不能** register \(u^2-\sqrt2\) 这类层。
- **P3-6 验收硬绑 T3+：** 在 T3+ 前仅允许双二次 + 可约降次 **临时** 路径，须登记 [known-divergences](../known-divergences.md)，**不得** 标 P3-6 完成。

### 5.2 与 dense poly1（D1–D4 ✅）

- P3-6 **主路径**：sparse `Poly<AlgExtCPolyCoeff>` + `PolyCoeff`（与 backlog §6 一致）。
- T3 `ParentBlockRing` 已可用于 K 上稠密运算；若 resolvent 实现选用 dense，走 **O1/O2**，非必需。

### 5.3 与 giac-solve

| 现状 | P3-6 后 |
|------|---------|
| `quadratic_rootof_roots` / `biquadratic_rootof_roots` | 迁入 `Poly<AlgExtC>::roots`（P4-1） |
| `general quartic rootof` → NotImplemented | 删除；改调 poly 层 |
| `solve` 形状表 | P4-6 factor + 递归 `roots` |

---

## 6. 测试与 DoD

### 6.1 单元 / 集成

| 用例 | 期望 |
|------|------|
| `roots(x²−2, x)` over K=ℚ | 两分支 `AlgExtC`，`eq_mod` ±√2 |
| `roots(t³−2, t)` | 三实根/一支实 + 复共轭（按域设定） |
| **resolvent 中间量** | 对已知四次手算 resolvent 系数 |
| **`solve(t⁴+t+1=0, t)`** | 四根；无 NotImplemented |
| `solve(t⁴−2=0, t)` | 四根（可与 P2-2 双二次对照） |
| T3+a | `adjoin(ℚ(√2), u²−√2)` 维数 4 + `eq_mod` |
| 可约 `(t²+1)(t²+2)` | factor 降次后递归（P3-7；可标 follow-up） |

### 6.2 DoD（P3-6 关单）

- [ ] **G0–G5** 全部 ✅
- [ ] `cargo test` 含 `PolyAlgExt` / solve 四次用例绿
- [ ] 无新增 `giac-solve` 四次特判
- [ ] 更新 [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) P3-6 行；偏离登记 `known-divergences`

---

## 7. 风险

| 风险 | 缓解 |
|------|------|
| T3+ registry 与 `min_poly_over_q` 张量基错误 | 先 T3+a 小域手算 + `eq_mod` |
| resolvent 分支 / adjoin 层爆炸 | 缓存 adjoin；同 splitting field 用 T2 嵌入 |
| `AlgExtCPolyCoeff::coeff_zero/one` 硬绑 rational 域 | roots 入口显式传入 K 的 field handle；或从 poly 项推断 field |
| 与 upstream `_EXT` 输出形状不一致 | 测 `eq_mod` 而非字符串；登记已知偏离 |

---

## 8. 交叉引用

| 文档 | 章节 |
|------|------|
| [GIAC-poly-algext-backlog](GIAC-poly-algext-backlog.md) | §5.1 四次策略、§6 P3-6、§7 P4-6 |
| [GIAC-lazy-common-tower-plan](GIAC-lazy-common-tower-plan.md) | §Phase T3+、§11 solve |
| [GIAC-dense-poly1-refactor](GIAC-dense-poly1-refactor.md) | D5 可选 |
| [GIAC-algext-adoption](GIAC-algext-adoption.md) | §8.4 目标架构 |

**AFK 建议：** 先 **PR-A（T3+）** 与 **PR-B/C（P2-1/6）** 可并行两人；**PR-E** 依赖前三项全部绿后再开。  
**续作（FieldSession）：** 见 [GIAC-poly-roots-field-session-plan.md](GIAC-poly-roots-field-session-plan.md) — 在 PR-B 与 PR-E 之间插入 **PR-B′…E′**，消除裸 `field` + 手工 align，解决四次超时与 G5。
