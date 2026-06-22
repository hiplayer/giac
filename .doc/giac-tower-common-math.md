# GIAC 塔式 `common` — 数学背景与可审计正确性

**状态:** living doc（T4a/T4b 实现参考）  
**代码:** `giac-core::algebra::ext_tower`（`compute_common_tower` / `compute_common_flatten` / `subfield_common_pair`）  
**计划:** [GIAC-lazy-common-tower-plan.md](issues/GIAC-lazy-common-tower-plan.md)

---

## 1. 我们在算什么

| 问题 | 数学对象 | 实现 |
|------|----------|------|
| 单步扩张 | \(K = F(\alpha) \cong F[x]/(m(x))\) | `adjoin_irreducible` |
| 并列域公共域 | compositum \(K_1 \vee K_2\)（含两者最小扩张） | T4a `compute_common_tower` 或 flatten |
| 子域包含 | \(K_1 \subseteq K_2\) 时 common = 提升嵌入 | T4b `subfield_common_pair` |
| 跨域运算 | 固定 \(\iota_i : K_i \hookrightarrow K\) 后在同一 ambient 做 `element_*` | `align_elements` + `FieldEmbedding` |

**giac-rs 约定：** 测试规格与 golden 等价优先；实现可偏离 C++ 行级路径，但 compositum / 嵌入须数学正确（见 `.cursor/rules` 与 `known-divergences.md`）。

---

## 2. 算法与可读参考

### 2.1 塔式 compositum（T4a，推荐）

**构造：** 设 \(K_1 = \mathbb{Q}(\alpha)\)，\(K_2 = \mathbb{Q}(\beta)\) 为次数 \(n,m\) 的简单扩张，\(m_\beta \in \mathbb{Q}[x]\) 在 \(K_1\) 上仍不可约（或至少 `adjoin` 给出次数 \(nm\) 扩张）时：

\[
K_1 \vee K_2 \;\cong\; K_1(\beta) \;\cong\; \mathbb{Q}(\alpha,\beta)
\]

**实现：** 选较小次数域为 parent，`adjoin_irreducible(parent, minpoly(sibling))`；  
`embed_parent` = 子域链嵌入；`embed_sibling` = \(\mathbb{Q}(\beta) \to K_1(\beta)\) 由 \(\beta \mapsto\) 新生成元 决定。

**经典教材（域论基础）：**

- S. Lang, *Algebra* — §V.1 代数扩张、有限扩张次数公式 \([K_1(\beta):K_1]=\deg m_\beta\)。
- D. S. Dummit & R. M. Foote, *Abstract Algebra* — §13 域扩张、复合域。
- I. Stewart, *Galois Theory* — compositum 与 Galois 群对应（概念层）。

### 2.2 Flatten 本原元（Phase 0 fallback）

**构造（Kronecker / 本原元试探）：** 在 \(K_1 \otimes_\mathbb{Q} K_2\) 的矩阵模型中，对整数 \(k\) 令 \(\theta = \alpha + k\beta\)，计算 \(\theta\) 在 \(\mathbb{Q}\) 上的最小多项式；当 \(\deg m_\theta = nm\) 时 \( \mathbb{Q}(\theta) \cong K_1 \vee K_2\)。

**实现：** `common_primitive_sum` + `char_poly_matrix`（Newton 恒等式），\(k=1..12\)。

**算法文献（可计算域论 / CAS）：**

- H. Cohen, *A Course in Computational Algebraic Number Theory* — §4.5 相对扩张、矩阵表示；与 giac 本原元 common 同类思想。
- M. Pohst & H. Zassenhaus, *Algorithmic Algebraic Number Theory* — 数域表示与嵌入。
- J. von zur Gathen & J. Gerhard, *Modern Computer Algebra* — §16 格罗埃纳基与理想运算（间接用于 minpoly）。

**局限（为何 T4a 替换默认）：** 特征多项式 + \(k\) 搜索贵；登记为扁平 `Adj{parent:Base}`，丢失塔结构；坐标为 \(\theta\) 幂基，与 T3 张量基不一致（plan §12.9）。

### 2.3 子域嵌入（T2/T4b）

**事实：** 若 \(K = K_{\text{parent}}(\gamma)\) 为一步 adjoin，则 \(K_{\text{parent}} \hookrightarrow K\) 在 operational 基下为块对角（parent 坐标进 block \(u^0\)）。

**参考：** 任何标准域论教材的「塔定理」\( [K:\mathbb{Q}] = [K:F][F:\mathbb{Q}] \)；矩阵形式见 Cohen §2.2（数域的矩阵表示）。

### 2.4 与 upstream giac 的关系

giac C++ `common_EXT` / 代数扩展在 `modpoly` / `algext` 模块（见 `.doc/module-division.md`）。Rust 侧 **绑定数学语义 + 测试**，不绑定其 profiler 或具体 \(k\) 上界；flatten 路径保留作 bisect fallback。

---

## 3. 以后要形式化证明：推荐思路

形式化目标不是一次证明「整个 giac-core」，而是 **分层证书（layered certificates）**，与 Rust 代码一一对应。

### 3.1 工具链选型

| 层级 | 工具 | 已有库 |
|------|------|--------|
| 代数与域论 | **Lean 4 + Mathlib4** | `Field`, `IntermediateField`, `Algebra`, `FiniteDimensional`, `Polynomial` |
| 可选验证辅助 | Isabelle/HOL (*Algebra*) | 较成熟但 CAS 对接少 |
| Rust 代码 | **Kani / Prusti** | 仅适合内存安全、小循环；**不适合**整个 minpoly 算法 |
| 规格测试 | `assert_equiv` + Sage/SymPy oracle | 项目现有 conformance 路线 |

**建议主路径：** Mathlib 证明数学引理 → Rust 测试 + 文档引用引理编号 → 关键小函数用 `#[cfg(test)]` 性质测试（proptest）。

### 3.2 证明分层（自底向上）

```text
L0  数据表示
    CoordsQ 长度 = dim(K)；poly1 顺序 ↔ Polynomial (reverse …)

L1  域公理（相对 F）
    element_add/mul 封闭；embed_rational 是环同态 F → K

L2  单步 adjoin
    adjoin_irreducible 给出 K ≅ F[x]/(m)；generator_coords = x̄

L3  塔运算
    element_*_tower 等价于 mod 层 minpoly 的多项式算术（张量基）

L4  嵌入
    FieldEmbedding.matrix 是 F-线性；compose = 矩阵乘
    try_subfield_embedding 沿 parent 链 = 块对角复合

L5  compositum（核心）
    引理 A：m_β 在 K_1 不可约 ⇒ [K_1(β):Q] = [K_1:Q]·deg m_β
    引理 B：K_1(β) 同时包含 K_1 与 Q(β) 的拷贝
    定理：compute_common_tower 的 ambient 是 compositum（在同构意义下唯一）

L6  flatten fallback
    引理 C：本原元定理 + 一般元素 θ=α+kβ 对「几乎所有 k」生成全扩张
    定理：common_primitive_sum 当 deg=m 时给出 K_1∨K_2 的 primitive 表示

L7  align / common 接口
    align_elements 三分支覆盖互斥情况；embedding_for 与 cache id 排序无关
```

**Rust 侧 DoD：** 每层至少一条 **golden + 小域手算** 回归；L5/L6 用同一对 \((\sqrt2,\sqrt3)\)、\((\sqrt2,\sqrt[3]2)\) 交叉验证 `min_poly_over_q` 次数与 `element_eq_mod`。

### 3.3 Mathlib 可落地的第一批引理（示例）

在 Lean 4 中可优先形式化（名称示意）：

1. `FiniteDimensional.finrank_mul` — 塔次数乘积（L5 引理 A 的特例）。
2. `IntermediateField.adjoin_simple_to_polynomial` — 简单扩张与商环同构（L2）。
3. `Algebra.IsAlgebraic.compose` — 代数元复合仍代数（保证 adjoin 良定）。
4. **Compositum 唯一性** — 用「最小同时包含 \(K_1,K_2\) 的域」定义，证明同构唯一（L5 定理框架）。

Rust `compute_common_tower` 的正确性陈述（伪规范）：

```text
∀ (K1 K2 : SimpleExt Q) (parallel : ¬ subfield K1 K2),
  let K := adjoin_tower pick_parent K1 K2 in
  ∃ ι1 ι2, ringHom K1 K ∧ ringHom K2 K ∧ finrank Q K = finrank Q K1 * finrank K1 K2_layer
```

### 3.4 工程实践（短期不必等证明）

1. **双路径等价测** — `tower-common` 开：`min_poly_over_q` 次数、生成元平方、align 后和（快测见 `t4a::tower_common_invariants_*`；完整 flatten 对照见 `#[ignore]` 的 `tower_common_matches_flatten_*`）。
2. **已知偏离登记** — 坐标基不同但 `eq_mod` 同；写入 `known-divergences.md`。
3. **不可约性** — 当前 `adjoin` 不验证 irreducible；形式化前应在 Rust 加可选 `debug_assert` 或文档假设。
4. **CI 分层** — 默认 T4a（`default = ["tower-common"]`）；bisect 用 `--no-default-features` 跑 flatten；证明脚本与 CI 解耦。

---

## 4. 可审计的正确性验证（giac-rs 口径）

本节所称 **可审计的正确性验证** 指：用文档、测试与交叉对照积累可复查的证据链；**不是**狭义 machine-checked 全库形式化证明（后者见 §3 层 D / Lean 引理，长期可选）。

**不把整个 CAS 写进 Lean 时，工程上如何积累这类证据：**

### 4.1 证据光谱（由重到轻）

```text
Lean/Coq 全证明  →  Refinement/提取  →  性质测试  →  Golden oracle  →  手算小例
  （极少全做）        （关键核心）       (proptest)    (giac check)      (单元测试)
```

**giac-rs 当前主战场：** Golden + conformance + 小域代数回归（见 [conformance-testing.md](conformance-testing.md)）。  
这与 Sage、PARI、upstream giac 的常见做法一致：**测试即规格**，不是证明助手即规格。

**Rust 内存安全工具（Kani / Prusti / MIRI）** 只管 panic/越界，**不**覆盖 `element_mul` 的域论语义；代数正确性不走这条线。

### 4.2 四层工程做法（推荐落地顺序）

#### 层 A — 规格冻结 + Oracle（已在用）

| 做法 | 说明 |
|------|------|
| Golden / `assert_equiv` | C++ giac check 为 oracle；数学等价即可，字符串可 normalize |
| 小域手算 | √2、√3、K₁(β) 张量坐标表（plan §12.9） |
| 文档陈述 | 每个算法 1 段数学构造 + 1 个最小例子（本文 §1–2） |
| **dense ↔ sparse 审计（D4）** | `giac-poly::dense::convert`：`reverse_coeffs`、`sparse_ascending_to_dense_high_first`、`dense_high_first_to_sparse`；单测 `dense::tests::{sparse_dense_high_first_roundtrip,dense_div_rem_matches_*}` |

**DoD：** 改 `align_elements` / `common_*` 必带 S0 逆序或子域/并列用例。

#### 层 B — 双实现 / 交叉验证（T4a 核心）

| 做法 | giac-rs 实例 |
|------|----------------|
| 同算法两实现 | `compute_common_tower`（默认） vs `compute_common_flatten`（`--no-default-features`） |
| 不变量对齐 | 次数乘积、α²=2、embed 后 `element_add` 非零 |
| 外部 oracle | 可选 Sage/SymPy snippet 对 dim≤4 算 minpoly（未强制进 CI） |

**注意：** 两路径 **coords 字面可不同**（flatten θ 基 vs 塔张量基）；比较用 `min_poly_over_q` 次数、域内 `element_eq_mod`，勿比坐标向量逐分量（除非同一 `field` 句柄）。

**慢测策略：** flatten 特征多项式对照放 `#[ignore]`，文档注明 `cargo test -p giac-core -- --ignored`；日常 CI 跑默认塔路径快测。

#### 层 C — 性质测试 + 有限模型（下一步，性价比高）

| 做法 | 适用 |
|------|------|
| **proptest** | 随机 `(K1,K2)` simple-over-ℚ、逆序 `align_elements`、`embedding_for` 与 `add` 交换 |
| **GF(p) 镜像** | 在有限域复刻 adjoin/common，次数≤3 时可 brute force 验证同态 |
| **cheap 不变量** | `dim`、`trace`、eval at rational point；避免比全表达式 |

**DoD：** 每个新 Layer（L4–L7）至少 1 条性质或 GF(p) 测，登记在 `ext_tower` tests 或独立 `tests/algebra_props.rs`。

#### 层 D — 证明助手（长期，只证引理）

| 谁证什么 | 工具 |
|----------|------|
| 塔次数、adjoin 同构、compositum 存在性 | Lean 4 + Mathlib4 |
| Rust 实现 | 注释链接定理名；**不**要求逐行提取 |
| CI | `cargo test` 与 `lake build` 分 job |

**行业参照：** CompCert/seL4 全证内核；HACL* 证 crypto 原语；**开放 CAS 几乎不 full proof**。Mathlib 作语义权威，不作 Rust 逐行证书。

### 4.3 模块级验证优先级（ext_tower）

| 优先级 | 模块 / API | 证据类型 |
|--------|------------|----------|
| P0 | `align_elements`, `embedding_for` | S0 逆序 + cache 预热 |
| P0 | `fold_algext_sum` + `embed_rational` | K₂ 有理项 merge 回归 |
| P1 | `compute_common_tower` | 塔结构 + dim；√2+∛2 快测（默认） |
| P1 | `subfield_common_pair` (T4b) | `common(K₁,K₂)=K₂`、不 flatten |
| P2 | `compute_common_flatten` | `#[ignore]` 与 tower 对照；bisect fallback |
| P3 | `poly_reduce` monic 假设 | 文档 + `debug_assert`（形式化前）；T3 见 [GIAC-dense-poly1-refactor §4.4](issues/GIAC-dense-poly1-refactor.md) |
| P3 | poly1 坐标 layout | `giac-poly::dense::convert` roundtrip 测（D4） |

### 4.4 CI 与 feature 约定

```bash
# 默认：T4a 塔 compositum + 快测
cargo test -p giac-core

# bisect / Phase 0 flatten 回滚
cargo test -p giac-core --no-default-features

# 可选：flatten 与 tower minpoly 慢对照
cargo test -p giac-core tower_common_matches_flatten -- --ignored
```

Lean 证明仓库（若另建 `giac-proofs/`）**不**阻塞 Rust CI。

**并发：** R4/R5 后无进程级 `field_registry`；`ext_tower` 单元测可并行（`cargo nextest` 默认多核）。

### 4.5 与项目其它文档的关系

- **测试规格 authoritative：** [conformance-testing.md](conformance-testing.md) §7  
- **已知数学/表示偏离：** [known-divergences.md](known-divergences.md)  
- **迁移原则（测试即规格、不复现 C++ bug）：** [rust-migration-plan.md](rust-migration-plan.md) §6.4  

---

## 5. 索引

| 主题 | 文档 |
|------|------|
| 实施阶段 | [GIAC-lazy-common-tower-plan.md](issues/GIAC-lazy-common-tower-plan.md) |
| 坐标基对照 | 同上 §12.9 |
| dense poly1 转换 | [GIAC-dense-poly1-refactor.md](issues/GIAC-dense-poly1-refactor.md) §4.2；`giac-poly::dense::convert` |
| adoption B-05 | [GIAC-algext-adoption.md](issues/GIAC-algext-adoption.md) |
| 四次 roots Lean 4 | [GIAC-poly-quartic-lean4-verification.md](issues/GIAC-poly-quartic-lean4-verification.md) |
| 测试规格 | [conformance-testing.md](conformance-testing.md) §7 |
| 已知偏离 | [known-divergences.md](known-divergences.md) |
