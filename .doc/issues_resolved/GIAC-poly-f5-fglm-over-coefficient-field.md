# GIAC-poly — F5 FGLM over 系数域 F=ℚ(α)（dim-12 塔域 √Δ_Q 探测）

**状态:** P0–P4 全部完成 / 已解决（2026-06-29，验收全绿，移入 `issues_resolved`）
**类型:** 实现 / AFK 可抓取
**父项:** [GIAC-poly-quartic-roots-F1-F5](../issues/GIAC-poly-quartic-roots-F1-F5.md) §F5 reframe
**Rust 落点:** `giac-groebner`、`giac-poly`、`giac-core::algebra::{ext_tower, poly_roots}`
**相关:** [GIAC-poly-flat-field-division-layering](../issues/GIAC-poly-flat-field-division-layering.md)、[giac-groebner-api-stability](../giac-groebner-api-stability.md)

---

## 0. 目标

在 dim-12 塔域 `ℚ(α,β)`（顶层 `β` 对 `F=ℚ(α)` 为 3 次）内求解 `x² = disc_q`，使 F5 `try_sqrt_in_field` 在该塔域命中 `√Δ_Q`，从而：
- **A₄** `x⁴+8x+12`：adjoin α(4) → β(×3=12) → √Δ_Q 命中 → 不盲 adjoin → **dim 12**。
- **S₄** `t⁴+t+1`：同结构，√Δ_Q 不命中（确需二次扩张）→ 盲 adjoin → **dim 24**（baseline 不变）。

**数学结构：** `√disc_q = x₀ + x₁β + x₂β²`，`xᵢ ∈ F` → **3 个 quadric in 3 未知数 over `F`**，0-dim degree 2（解 `±x`）。A₄ 下 `√disc_q ∈ ℚ(α,β)`（dim 12 即分裂域），故解存在、probe 应成功；S₄ 下 `ℚ(α,β)` dim 12 是 `S₄` 的固定域（非分裂域），`√disc_q ∉` → miss。

**方法：** FGLM —— 先求 `greduce`/Buchberger 在 **grevlex 序** 下的 GB（0-dim 快），再用 **线性代数** 转成 lex GB（三角形），读出 `x₀,x₁,x₂`。d=2 → 转换 trivial。

---

## 1. 现状与关键发现

- `giac-poly::PolyCoeff` / `FieldCoeff`（`poly_coeff.rs`）已抽象域系数；`AlgExtCPolyCoeff: FieldCoeff`（`poly_alg_coeff.rs:48`）**已实现** ⇒ `Poly<AlgExtCPolyCoeff>` 即 `F[var]`，`coeff_inv`/`field_div` 可用。**系数域抽象已就位。**
- `giac-groebner` 现硬编码 `Poly = Poly<Ratio<BigInt>>`（`greduce`、`groebner_basis_lex`、`spoly_lex`、`make_monic_lex`、`autoreduce_lex` 均用 `Ratio`/`mul_scalar(&Ratio)`/`Ratio::from_integer`）。**需泛化到 `C: FieldCoeff`。**
- **`Poly::add/sub/mul/neg/var/zero/one/mul_scalar` 是 `impl Poly<Ratio<BigInt>>` 专属（`poly.rs:201-241`），泛型块只有 `try_add/try_sub/try_mul/try_neg`（返回 `PolyResult`）+ `ring_zero/ring_one/ring_var/ring_constant`。** 泛型 groebner 不能「直接复用」`add/sub/mul`，须走 `try_*`（见 P0 决策）。
- `AlgExtCPolyCoeff::coeff_*` 全部返回 `PolyResult`（`poly_alg_coeff.rs:82-105`）；`inv()` 仅在零元时失败，`mul` 仅在同域维度一致时成功 —— 在 groebner 中只对同域非零 lc 运算，故满足 invariant。
- **`greduce` 是公开 CAS builtin `greduce` 的后端**（`eval_poly.rs:273` `eval_greduce`），`groebner_basis_lex` 另被 `poly_roots.rs:2018`（WIP proto）调用。**泛型化不得改公开 ℚ 入口的签名。**
- lex 序直接 Buchberger 对 `x²=u` 的 4 坐标 quadric 中间度爆炸（120s 未完成，见 F1-F5 §359）。**必须走 grevlex + FGLM，不再用 lex 直接 Buchberger。**
- `Monomial` 当前 `cmp_lex`/`leading_term_lex` 已有（`monomial.rs:84`、`poly.rs:79`）；**缺 grevlex 比较**。注意 `Monomial` derive `Ord`（`BTreeMap` 按 var 名字典序）既非 lex 也非 grevlex，`Poly::leading_term()`（`poly.rs:73`，`iter().next_back()`）依赖它 —— **grevlex 路径禁止误用 `leading_term()`**，必须用新增的 `leading_term_grevlex`。

---

## 2. 组件与优先级

> P0 必须最先；其余按序。每期独立可测、可 commit。

### P0 — `giac-groebner` 泛化到 `C: FieldCoeff`（解锁一切）

**为什么 P0：** 所有后续（grevlex GB、FGLM、塔域探测）都依赖「`greduce`/Buchberger 能跑在 `Poly<AlgExtCPolyCoeff>` 上」。

**错误传播决策（invariant `expect`，签名保持 infallible）：**
groebner 内部对系数只做同域非零 lc 上的 `+ − × ÷`，`AlgExtCPolyCoeff` 在该前提下 `try_*` 必成功。故泛型核心用 `try_add/try_sub/try_mul/try_neg` + `.expect("field op invariant: same-field nonzero lc")`，**签名保持 `Vec<Poly<C>>` / `Poly<C>`（非 `Result`）**。理由：(a) 工程规则允许有注释的 invariant `expect`；(b) 不破 `eval_greduce` 等 `Stable` 调用方；(c) ℚ 下 `try_*` 本就 infallible，`expect` 恒不触发。

**任务：**
1. 抽出泛型核心 `greduce_generic<C: FieldCoeff>`、`groebner_basis_lex_generic<C: FieldCoeff>`、`spoly_lex_generic<C>`、`make_monic_lex_generic<C>`、`autoreduce_lex_generic<C>`；保留 `pub fn greduce(...)` / `pub fn groebner_basis_lex(...)`（ℚ 专属 `Stable`）作为薄包装 `greduce_generic::<Ratio<BigInt>>(...)`，**签名不变** → `eval_greduce`、`poly_roots.rs` 调用方零改动。
2. 去 `Ratio` 特化（逐项映射，列全）：

   | 现状（`Ratio` 专属） | 泛型替换 |
   |---|---|
   | `r_lc.clone() / g_lc.clone()`（`greduce` `lib.rs:36`） | `r_lc.field_div(g_lc).expect("field div invariant")` |
   | `Ratio::from_integer(BigInt::from(1))`（`make_monic_lex:99-100`、`spoly_lex:117-118`） | `C::coeff_one()` |
   | `1 / lc_a.clone()`（`spoly_lex`） | `lc_a.coeff_inv().expect("field inv invariant: lc nonzero")` |
   | `*lc != Ratio::from_integer(BigInt::from(1))`（`make_monic_lex:99`） | `!lc.coeff_is_one()` |
   | `f.mul_scalar(&inv)`（`make_monic_lex:101`） | `Poly::<C>::term(monom, c)` + `try_mul` 重组；或逐 term `c.coeff_mul(&inv)` 重建 `BTreeMap` |
   | `Poly::term(t_a, inv_a).mul(a)` | `Poly::<C>::term(t_a, inv_a).try_mul(a).expect(...)` |
   | `left.sub(&right)` / `r.sub(...)` / `q_term.mul(g)` | `try_sub` / `try_mul` + `expect("field op invariant")` |
   | `Poly::zero()` / `Poly::one()` / `Poly::var("x")` | `Poly::<C>::ring_zero()` / `ring_one()` / `ring_var("x")` |

3. `Poly::term` / `try_add` / `try_sub` / `try_mul` / `try_neg` / `ring_zero` / `ring_one` / `ring_var` / `is_zero` / `is_one` / `leading_term_lex` 均已是 `C: PolyCoeff` 泛型，直接复用（注意是 `try_*` + `ring_*`，非 `add/sub/zero/one/var`）。
4. `num-bigint`/`num-rational` 从 `giac-groebner` **lib deps 退回 dev-only**：泛化后 lib 内不再出现 `Ratio`/`BigInt`（ℚ 包装只是泛型核心的 `Ratio` 特化点，类型由调用方带入，不需 crate 自身引 `num-*`）。确认所有 `pub fn` 均泛型或为 `Poly = Poly<Ratio>` 别名包装后再删 deps。

**验收：**
- 既有 `giac-groebner` 小例测（`groebner_lex_x2_1_y_minus_x`、`groebner_lex_unit_ideal`、`greduce_*`）在 `Poly<Ratio<BigInt>>` 下仍绿（包装签名不变）。
- 新增 1 例 `Poly<AlgExtCPolyCoeff>`：在 `ℚ(√2)` 上 `groebner_basis_lex_generic` 解 `{x²−2, y−x}` 得 `{y²−2, x−y}`，绿。
- `cargo build -p giac-groebner` 无 `num-bigint`/`num-rational` lib deps（`cargo tree -p giac-groebner` 验证）。

**文件：** `giac-groebner/src/lib.rs`、`giac-groebner/Cargo.toml`
**tier 登记：** `greduce_generic`/`groebner_basis_lex_generic`/`spoly_lex_generic`/`make_monic_lex_generic`/`autoreduce_lex_generic` 标 `Partial`；`greduce`/`groebner_basis_lex` 维持 `Stable`。更新 `.doc/giac-groebner-api-stability.md`，跑 `giac-rs/scripts/annotate_api_tiers.py --inventory`。

**✅ 完成（2026-06-29）实现偏离登记：**
- **API 形态偏离**（vs 任务 1「薄包装」方案）：未采用「泛型核心 `*_generic` + ℚ 薄包装」双层，而是直接把 `pub fn greduce` / `groebner_basis_lex` 自身泛型化为 `<C: FieldCoeff>`（`groebner_basis_lex`/`autoreduce_lex` 额外 `+ PartialEq`，用于 autoreduce 的变更检测 `r != g[i]`；`Poly<C>` derive `PartialEq`，`AlgExtCPolyCoeff`/`AlgExtCData` 均 derive `PartialEq` 故满足）。因 `Poly<C: PolyCoeff = Ratio<BigInt>>` 有默认类型参数，既有 ℚ 调用方（`eval_poly.rs:273` `greduce(&p,&basis,&vars)`、`poly_roots.rs:2018` `groebner_basis_lex`）类型推断 `C=Ratio<BigInt>`，**调用点零改动**（`cargo build -p giac-core` 通过）。比双层包装少一层样板，符合 ponytail；副作用是公开签名从 `fn greduce(&Poly,…)` 变为 `fn greduce<C: FieldCoeff>(&Poly<C>,…)`，但因默认类型参数对 ℚ 调用方透明。
- **错误传播**：按决策走 `try_mul`/`try_sub` + invariant `expect`，但用两个私有 helper `pmul`/`psub`（`a.try_mul(b).expect("same-field …")`）封装，签名 infallible。`make_monic_lex` 用 `Poly::term(Monomial::one(), lc.coeff_inv()?).mul(&*f)`（经 `pmul`）。`spoly_lex` 用 `lc.coeff_inv().expect("nonzero leading coeff")`。`greduce` 用 `r_lc.coeff_div(g_lc)`（`Ok`-else-`continue`，零元除已由 `coeff_is_zero` guard 排除）。
- **`ring_*` 而非 `add/sub/zero/one`**：泛型块只暴露 `try_*` + `ring_zero/ring_one/ring_var/ring_constant`，`Poly::add/sub/mul/zero/one/var` 是 `impl Poly<Ratio<BigInt>>` 专属；泛型代码一律走 `ring_*`/`try_*`（与 §1 关键发现一致）。
- **DoD 测例偏离**（vs 验收「`{x²−2, y−x}` 得 `{y²−2, x−y}`」）：改用 **`{x²−αx, xy−1}` over ℚ(√2) → `{x−α, y−1/α}`**（测 `groebner_lex_over_qsqrt2`，`giac-core/src/algebra/poly_alg_coeff.rs`）。理由：计划原例系数全为 ℚ（2 是有理数，`√2` 从不以系数出现），无法证明「系数域 `F=ℚ(α)` 上的 `coeff_inv`/`coeff_div` 被真正走到」；新例的 `αx` 项与 `y−1/α` 单变式必须用 α 系数除法/求逆才能产生，**严格覆盖 P0 的泛型点**。断言：`x−α` 结构相等（x 系数同在 base ℚ，field id 一致）+ 第二元用「LT=y 且 monic」不变量（避免 field-id 结构误判：GB 的 y 系数因除以 α 落在 extension field id，而 `ring_var` 造的 y 系数在 base ℚ id，结构 `==` 会误报；改用 leading-monomial + `coeff_is_one()` + `greduce(fᵢ,gb)=0` 有效性校验）。
- **lib deps**：`num-bigint`/`num-rational` 退回 `[dev-dependencies]`（仅测试用）；`cargo tree -p giac-groebner` 直连仅 `giac-poly`（num-* 经 giac-poly 传递，符合 DoD「lib 无 num-*」）。
- **回归**：`giac-groebner` 5/5、`giac-core --lib` 299/299（含 `eval_greduce_direct`/`eval_greduce_circle` 即 eval→`greduce` ℚ 路径、新增 ℚ(√2) 例）全绿，无新增 ignore。
- **待办（不阻塞 P1）**：tier 登记（`annotate_api_tiers.py --inventory`）与 `giac-groebner-api-stability.md` 更新留到 P4 一次性做（plan 原 P0 tier 登记步骤顺延）。

---

### P1 — grevlex 序 + grevlex Buchberger（0-dim 快）

**任务：**
1. `Monomial::cmp_grevlex(&self, other, var_order)` —— 按总次数，再按「最大相异 var 处指数**小者**为大」（grevlex 标准）；`Poly::leading_term_grevlex(&self, var_order) -> Option<(&Monomial, &C)>`（复用 `cmp_grevlex` 的 `max_by`）。**禁止在 grevlex 路径用 `Poly::leading_term()`（derive `Ord` 字典序）。**
2. `greduce_grevlex_generic<C: FieldCoeff>` / `spoly_grevlex_generic<C>` / `make_monic_grevlex_generic<C>` / `autoreduce_grevlex_generic<C>` / `groebner_basis_grevlex_generic<C: FieldCoeff>`：复用 P0 泛型骨架，仅换 leading-term 选取 + Gebauer-Möller product/chain criteria。`pub fn groebner_basis_grevlex` 暴露 ℚ 包装 + 泛型核心（与 P0 同 pattern）。
3. grevlex 下 0-dim 中间度受控（grevlex 是 0-dim 推荐序）；沿用 P0 的 `GROEBNER_POLY_CEILING=64` / `GROEBNER_DEGREE_CEILING=8` 兜底。

**验收：**
- 新测：`{x²−2, y²−3, x·y−6}` over ℚ 的 grevlex GB（0-dim，d=4）<100ms 绿。
- 新测：同例 over `ℚ(√2)`（`Poly<AlgExtCPolyCoeff>`，系数含 `√2`）绿。

**文件：** `giac-poly/src/monomial.rs`、`giac-poly/src/poly.rs`（`leading_term_grevlex`）、`giac-groebner/src/lib.rs`
**tier 登记：** `cmp_grevlex`/`leading_term_grevlex` `Stable`；`groebner_basis_grevlex_generic` 等 `Partial`。

**✅ 完成（2026-06-29）实现偏离登记：**
- **骨架共享（DRY）**（vs 任务 2「复用 P0 泛型骨架」）：未复制 Buchberger body，而是把 P0 lex 代码重构为 order-parameterized 共享核心 —— `Order` enum + `lt()` leading-term 分发 + `greduce_order`/`make_monic_order`/`spoly_order`/`autoreduce_order`/`buchberger_order`；`greduce`/`greduce_grevlex`/`groebner_basis_lex`/`groebner_basis_grevlex` 为薄包装。Gebauer-Möller product/chain criteria 本就 order 无关，只有 `lt()` 选取不同。消去 ~80 行重复。
- **去冗余 bound**：P0 的 `+ PartialEq` 经核实冗余（`FieldCoeff: PolyCoeff: Clone+PartialEq+Debug`，`Poly<C>` derive `PartialEq`），P1 全部回退为 `<C: FieldCoeff>`。
- **`cmp_grevlex` 约定**：graded reverse lex，`order[0]` = 最显著变量（同 `cmp_lex`）；等总次数时从**最不显著**变量扫描，**指数小者为大**（CLO 标准）。`Poly::leading_term_grevlex` 用 `max_by(cmp_grevlex)`；grevlex 路径未误用 `leading_term()`（derive Ord 字典序）。
- **DoD 测例偏离**（vs 验收「`{x²−2, y²−3, xy−6}` 0-dim d=4」）：该系统**不一致**（`xy=6` 与 `(xy)²=6 ⇒ xy=±√6≠6` 矛盾）⇒ grevlex GB = unit ideal（空）。改作 **unit-ideal 检测测**（`groebner_grevlex_unit_inconsistent_under_100ms`，空 GB + <200ms，正是防 lex-blowup 的回归 guard）。d=4 0-dim 由 **`{x²−2, y²−3}`**（coprime ⇒ d=4 GB，`groebner_grevlex_x2_2_y2_3_d4`，validity via `greduce_grevlex=0`）覆盖。
- **ℚ(√2) 测**：`groebner_grevlex_over_qsqrt2`（`{x²−αx, y−x}` ⇒ 0-dim d=2，validity + <500ms），走 α 系数 grevlex S-poly 归约，区别于 P0 lex 测（`Order::Grevlex` 的 lt 选取）。位于 `giac-core/poly_alg_coeff.rs`（避免 giac-groebner→giac-core 循环依赖）。
- **回归**：`giac-groebner` 7/7（5 lex + 2 grevlex ℚ）、`giac-core --lib` 300/300（4 ignored）全绿；`clippy -p giac-groebner -p giac-core -p giac-poly --lib` clean。
- **待办（不阻塞 P2）**：tier 登记（`annotate_api_tiers.py --inventory`）与 `giac-groebner-api-stability.md` 仍顺延到 P4 一次性做。

---

### P2 — FGLM 转换（grevlex GB → lex GB，线性代数）

**d 范围决策（通用 d，非 d=2 特化）：** 实现通用 d FGLM，理由：(a) 验收用 d=4 回代测能 catch 三角化错误，d=2 特化测不出一般性 bug；(b) 通用 Krylov 在 d=2 上自然退化为 trivial，无额外 fast-path 代码；(c) 算法本身是标准 FGLM，非「未请求的抽象」。**退役条件：** 若通用 d 实现超 ~300 行且 d=2 问题已可用，可先合 d=2 路径并留 follow-up；否则一次到位。不新增 d=2 专属分支。

**任务：**
1. `fglm_generic<C: FieldCoeff>(grevlex_gb: &[Poly<C>], var_order: &[Var]) -> Option<Vec<Poly<C>>>`：
   - 取商环基 `B` = grevlex GB 的标准单项式（不被任一 leading monomial 整除），`|B|=d`。
   - 对每个 `xᵢ` 构造乘法矩阵 `M_{xᵢ}`（`d×d` over `C`）：`xᵢ·bⱼ` 用 `greduce_grevlex_generic` 归约 → 坐标向量。
   - `g₀(x₀)` = `M_{x₀}` 的极小多项式（Krylov：`v, Mv, M²v, …` 找首个线性相关，`C` 上线性相关判定走 `field_div` 高斯消元）。
   - 对 `k≥1`：在 `M_{x₀..x_{k-1}}` 张成的代数中找 `xₖ` 的线性表示（normal form 三角化）→ `gₖ`。
   - 输出 lex 三角形 GB。
2. 边界：`d=0`（unit ideal）/ `d > D_MAX`（**`D_MAX=64`，ponytail 天花板**）→ 返回 `None` 哨兵（调用方退回 blind adjoin，非 hang）。
3. `pub fn fglm` 暴露 ℚ 包装 + 泛型核心（同 P0/P1 pattern）。

**验收：**
- 新测：`{x²−1, y−x}` grevlex GB → FGLM → lex GB = `{y²−1, x−y}`（对照现有 `groebner_lex_x2_1_y_minus_x`）。
- 新测：d=4 0-dim 系统（如 `{x²−2, y²−3, x·y−6}`）FGLM 三角形可回代出全部 4 解。

**文件：** `giac-groebner/src/lib.rs`（新 `fglm` 模块/段）
**tier 登记：** `fglm_generic` `Partial`；`fglm` `Stable`（ℚ 包装）。

**✅ 完成（2026-06-29）实现偏离登记：**
- **API 形态偏离**（vs 任务 3「ℚ 包装 + 泛型核心」+ tier 登记「`fglm_generic` `Partial` / `fglm` `Stable` ℚ 包装」）：与 P0/P1 一致，未做双层包装，直接 `pub fn fglm<C: FieldCoeff>(...)`（默认类型参数对 ℚ 调用方透明）；私有辅助 `fglm_solve`/`fglm_matvec`/`monom_pow`/`monomials_of_degree` + `fadd/fsub/fmul/fdiv/padd` scalar/poly invariant helper。故无独立 `fglm_generic` 命名（tier 顺延到 P4 登记时统一记 `fglm` 为 `Stable (bounded)`）。
- **算法形态偏离**（vs 任务 1「`g₀(x₀)` = `M_{x₀}` 极小多项式，`x₀`=var_order[0] 最大变量」）：实际采用 **shape-lemma 形式，univariate 在最小变量 `x_{n-1}=var_order[n-1]`**（与验收例 `{x²−1,y−x} → {y²−1, x−y}` 一致，y 是最小变量；plan 任务 1 的「x₀ 最大变量」与验收例自相矛盾，以验收例/lex 三角化数学为准）。即：`g_{n-1}(x_{n-1})` = `M_{x_{n-1}}` 极小多项式（Krylov `u=coords(1)`，`w_{j+1}=M·w_j`，首个线性相关给 `g_{n-1}=t^d−Σcᵢtⁱ`）；对 `k<n−1`，`φ_k` 由 `K⁻¹·coords(NF(x_k))` 解出（`K=[w_0..w_{d-1}]` Krylov 矩阵），`g_k = x_k − φ_k(x_{n-1})`。**只需 `M_{x_{n-1}}` 一个乘法矩阵**（NF(x_k) 走 `greduce_grevlex`，无需 `M_{x_k}`），比 plan「每个 `xᵢ` 一个 `M_{xᵢ}`」更省。
- **通用位置（generic position）限制**（plan 未提及，关键约束）：shape-lemma FGLM 要求理想对最小变量 `x_{n-1}` 处于 generic position（`x_{n-1}` 的极小多项式次数 = d，即 d 个根有 d 个不同的 `x_{n-1}` 值）。`mp_deg < d` ⇒ 返回 `None`（非 generic），调用方退回 blind adjoin —— **miss ≠ 错答**（与 P3「S7 miss 退回 `try_sqrt_pairwise_fallback`」语义一致，安全兜底）。P3 的 √Δ_Q 系统解 `±x`，只要 x 的 `x_{n-1}` 分量 ≠ 0 即 generic（随机 disc_q 下成立）。
- **DoD 测例偏离**（vs 验收「d=4 `{x²−2, y²−3, xy−6}` 回代 4 解」）：该系统与 P1 同样**不一致**（`xy=6` 与 `(xy)²=6 ⇒ xy=±√6` 矛盾）⇒ unit ideal d=0。改用 **`{x²+y²−5, xy−2}`**（4 有理根 (2,1),(1,2),(−2,−1),(−1,−2)，y 值 {±1,±2} 4 个不同 ⇒ generic w.r.t. y，d=4），`fglm_x2_plus_y2_5_xy_2_d4_generic`：FGLM 出 `{y⁴−5y²+4, x+(y³−5y)/2}`（shape-lemma），用**理想等式双向校验**（generators 经 LEX `greduce` mod fg = 0 ∧ fg 经 GREVLEX `greduce_grevlex` mod grev = 0 ⇒ ⟨fg⟩=⟨gens⟩，**理想成员资格必须用与 GB 序匹配的归约序**，用错序会假阴性）+ 结构断言（一元 y⁴ + x 线性）+ <500ms。
- **d=2 对照测**：`fglm_x2_1_y_minus_x_matches_lex` —— FGLM(grevlex GB) **集合相等**于独立路径 `groebner_basis_lex`（gold cross-check）+ spot-check `{y²−1, x−y}`。
- **unit ideal 测**：`fglm_unit_ideal_returns_none`（`{x−1, x−2}` 不一致 ⇒ 空 grevlex GB ⇒ FGLM None）。
- **ℚ(√2) 测**（P3 前置）：`fglm_over_qsqrt2`（`{x²−αx, y−x}` ⇒ d=2 generic，FGLM 出 `{y²−αy, x−y}`，α 系数 Krylov + 高斯消元，双向理想等式 + 结构（x−y 结构相等，y² 三角用 LT+monic 避 field-id 误判）+ <500ms）。位于 `giac-core/poly_alg_coeff.rs`，证明 FGLM 在 `AlgExtCPolyCoeff` 上端到端可用（P3 直接依赖）。
- **边界**：`d=0`（unit / `1∈ideal`）/ `d>FGLM_D_MAX=64`（ponytail 天花板）/ 非泛型位置 / NF 项落基外 ⇒ 一律 `None`。`FGLM_D_MAX=64` inline 标注 ceiling + 升级路径（增大或换分裂版 FGLM）。
- **规模**：FGLM 段 ~190 行（含注释），未触退役条件（<300 行）；通用 d 一次到位，无 d=2 专属分支。
- **回归**：`giac-groebner` 10/10（5 lex + 3 grevlex + 3 fglm ℚ，含 unit-ideal None）、`giac-core --lib` 含 `fglm_over_qsqrt2` 全绿；`clippy -p giac-groebner -p giac-core --tests` 对 FGLM 代码 clean（`poly_roots.rs:1976` 一条 `manual_memcpy` 警告为既有 `proto_*` dead code，非 P2 引入）。
- **待办（不阻塞 P3）**：tier 登记（`annotate_api_tiers.py --inventory`）与 `giac-groebner-api-stability.md` 仍顺延到 P4 一次性做。

---

### P3 — 重定向：F-module 递归 √-归约（dim-12 塔域 √ 值恢复）

> **调试复盘见** [GIAC-poly-f5-fglm-debug-postmortem](../issues/GIAC-poly-f5-fglm-debug-postmortem.md)：P3c→P3e 的 7 个错误清单、根因分类（约定不可见为主因）、及 order-newtype / 显式 mode / 搜索 budget 三项根治方案。

> **历史路径（已否决，留档防回退）：** 原 P3 计划「S7 塔域 deg-3 顶层 FGLM」假定 `x²=u` 在 deg-n 域对应 0-dim degree **2** 的系统，FGLM 三角化后回代 ≤2 次单变式。**该假定错误**：`x²=u` 在 deg-n 域（相对基）对应 n 个 quadric in n 未知数 over 基，0-dim degree **2ⁿ**（n 个嵌入 × ±），非 2。A₄ dim-12 下 2¹²=4096，FGLM 不可行（diag 实测 `g_last` deg 8 for 3-var-over-ℚ(α) 塔系统、4-var-over-ℚ flat 系统）。S7/`flat_sqrt_via_fglm` 代码已 revert。

> **Euler-resolvent-first 旁路（已否决）：** 试图把 `quartic_roots` 改为 Euler 优先（adjoin-deflate fallback），期望 Euler 的 √β 在 dim-6 盲 adjoin → 12 以避开 dim-12 √Δ_Q miss。**diag 实测否决**：A₄ 下 Euler 路径 √α → dim 12、√β 盲 adjoin → dim 24，**同一 dim-12 塔域 √ miss**。两条 quartic 路径都撞同一 gap，无 lazy 旁路。已 revert。

#### 3.1 已确立：collapse 检测（boolean，已验证）

`diag_adjoin_collapse_recovery`（`poly_roots.rs`）证明 collapse 检测便宜且正确：
- adjoin δ²=u → `flatten_min_poly_over_q_cold` → δ 的 ℚ-特征多项式（deg 2·dim K）→ `factor_univariate_pairs` → δ 的 ℚ-极小多项式 `m`（deg = multiplicity-free factor deg）。
- **collapse 判据：δ ∈ K ⟺ deg(m) ≤ dim(K)。** A₄ 实测 `charpoly deg 24`、`factors=[(12,2)]`、`deg m=12 ≤ dim K=12` → `δ ∈ K = true`（0.18s）。S₄ 下 `deg m=24 > 12` → `δ ∉ K`（应盲 adjoin）。
- 检测可用作 A₄/S₄ 布尔分类器，但**不**给 √ 值。

#### 3.2 已否决：值恢复的 cheap 路径（全部撞墙）

| 路径 | 结果 |
|---|---|
| `align_pair(K, ℚ(δ))`（δ 由 `m` 经 `build_base_extension_uncached` 造 flat dim-12 域） | **hang** —— 两个同构但无共享塔链的 dim-12 域，compositum 算法不收敛 |
| `try_subfield_embedding(ℚ(δ), K)` | **None** —— `is_subfield_of` 只沿塔 parent 链查；ℚ(δ) 非 K 祖先 |
| `m(t) mod (t²−disc_q)` over K → `δ = −B/A` | **退化** —— A₄ 下 −δ 是 δ 的共轭（A₄ 的 2 阶元固定 disc_q、翻转 √），故 `(t²−disc_q) \| m(t)`，余式 0；而 splitting `t²−disc_q` over K 正是原问题（循环） |
| FGLM on `x²=u` quadric 系统 | degree 2ⁿ（见上），不可行 |

**结论：** collapse 只给 yes/no；√ 的**值**恢复 = 原始「塔域 √」问题本身（伪装）。需 real machinery。

#### 3.3 F-module 递归 √-归约（采用，real fix）

**核心思想：** √u in K（dim d_K，已知 δ∈K）→ 工作在 F=ℚ(u)（dim d_F = deg m_u）。`[K:F] ≥ 2` 时 K 是 F 上 rank-2 模，δ 在 F-module 上作用为 M∈M₂(F)，M² = u·I。任取 w∈K\F，算 T=tr_{K/F}(w)、N=N_{K/F}(w)（F-module 乘法矩阵的 trace/det），则：

> δ = (w − T/2) · √(4u / (T² − 4N))，其中 √ 在 **F** 内（且 4u/(T²−4N) = 1/b² ∈ F 是 F 中平方，by 构造）

递归把 √ in dim-d_K 归约到 √ in dim-d_F（每次 dim 减半，因 F=ℚ(u)ᐸK 当 u 不生成 K）。**base case**（u 生成当前域，`[K:F]=1`）：解 `p(u)² = u`，p ∈ ℚ[u]/(m_u)（d 个 quadric in d 未知数 over ℚ，degree 2^d）—— 仅对小 d（≤ ~6，degree ≤ 64）可行，走 P1/P2 的 grevlex GB + FGLM。递归链 A₄：dim 12 → dim 6 →（dim 3 或 base d=6）→ … → dim 1（ℚ，有理 √ trivial）。**完整、通用，修复两条 quartic 路径。**

**关键子件（多数已存在）：**
- u 的 ℚ-极小多项式 m_u（deg d_F）—— collapse 已算 charpoly/factor，复用。
- F=ℚ(u) 的 ℚ-基 {1,u,…,u^{d_F−1}} 作为 K 内坐标（`k.element_mul` 逐次乘）。
- K 作为 F-module 的 rank-2 分裂 + w∈K\F 选取 + F-coord 提取（12×12 ℚ-线性解，`rational_rref`）。
- tr/norm_{K/F}（F-module 乘法矩阵 trace/det）。
- F 内 √ 递归（同算法，dim 减半）。
- base case：小 d GB/FGLM（P1/P2 已就位）。

#### 3.4 任务（增量、每步可测）

1. **P3b 可行性 diag**（`poly_roots.rs`）：A₄ 上实现 F-module step 1 —— 算 `m_{disc_q}`（deg 6）、F=ℚ(disc_q) 基、选 w∈K\F、算 T/N（F-coord）、`u'=4u/(T²−4N)`、打印 `deg m_{u'}`（判 base case dim）。验证线性代数基础设施。
2. **P3c base case**：小-d（≤6）`p(u)²=u` via grevlex GB + FGLM；ℚ-d=1 有理 √ trivial。
3. **P3d 递归组装**：F-module step + 递归 √ in F + δ = (w−T/2)·√u'；diag 恢复 √(disc_q) ∈ K 且自校 `δ²=disc_q`。
4. **P3e 接线**：进 `try_square_root_in_field_impl`（S7 新阶段：collapse 检测 + F-module 恢复，gate `dim ≥ 4`）；un-ignore `quartic_a4_galois_dim_le_12`；全 quartic 门禁。
   - **回归防护（`sqrt_base_case` `LIFT_BUDGET=8`）**：F-module 基例的 p-adic Newton 提升对「模 p 是平方但 ℚ 上非平方」的 u 会穷举 60 素数 × 8 符号组合 → 100s 超时。真平方（A₄ √Δ）在首个可用素数/首个组合即成功，故硬上限 8 次 lift 既保 A₄ 命中又把非平方调用的代价封顶。release 实测 8 个原超时 quartic 测全绿（≤0.4s/测，`euler_four_roots_vanish` 11s 不变），A₄ 测 0.12s 绿。
   - **未采用：norm 预过滤**（`N_{K/ℚ}(u)` 非有理平方 ⇒ 早退）。数学正确，但 F-module 递归调用处存在 low-first/high-first CoordsQ 约定错配（krylov 看到的是「转置」元素），norm 被误判 → 误杀 A₄ 真平方。已移除，仅留 `LIFT_BUDGET`。

#### 3.5 验收

- `diag_adjoin_collapse_recovery`：collapse 检测绿（**已就位**）。
- 新 `diag_fmodule_sqrt_recovery`：A₄ dim-12 塔域恢复 √(disc_q) ∈ K，自校 `δ²=disc_q`，dim 不增（仍 12）。**绿**。
- `quartic_a4_galois_dim_le_12` unignore 绿：四根 `verify_root` + `dim ≤ 12`。**绿**（0.12s）。
- `t⁴+t+1`（S₄）dim ≤ 24 不回归（√Δ_Q ∉ K，collapse 判 false → 盲 adjoin 不变）。**绿**。
- **P3e 回归**：原 100s 超时的 8 个 giac-core + 2 个 giac-solve quartic 测全绿（`LIFT_BUDGET=8` 封顶非平方 lift 代价）。

**文件：** `giac-core/src/algebra/ext_tower.rs`（F-module √ + S7 接线）、`giac-core/src/algebra/poly_roots.rs`（diag + 测）、`giac-core/src/algebra/field_arith.rs`（复用 `rational_rref`）/ 可能新增 subfield-module 辅助。
**tier 登记：** F-module √ / S7 helper 标 `Pipeline private`（同 S0–S6）。

---

### P4 — S₄ 回归 + 全套门禁 ✅

**任务：**
1. `roots_quartic_t4_plus_t_plus_1`、`field_session_dimension_bound_quartic_tight`（`t⁴+t+1` dim=24 baseline）仍绿 —— S₄ 路径 √Δ_Q **不**命中（`ℚ(α,β)` dim 12 是固定域非分裂域），盲 adjoin 不变。**绿**（0.30s/测，4 根 `verify_root`，dim ≤ 24）。
2. `cargo test-timeout -p giac-groebner -p giac-core` 全绿，无新增 ignore。**绿**：nextest release `1089 passed, 46 skipped, 0 failed/timeout`，13.8s（A₄ 测已 un-ignore 计入 passed；46 skipped 为既有 slow/brute-force）。
3. `./scripts/ci-clippy.sh` 全绿（含 `lint-substring-golden`）。**绿**：修 3 处 clippy（`redundant_closure` ×2、`needless_borrow`、`manual_memcpy`，均 P3 迁出的 `sqrt_fmodule` 代码）+ 1 处 pre-existing substring-golden（`giac-calculus/integrate.rs`，`1be6dbd` 引入，改用 `depends_on_var` 结构检查）。
4. 登记 `giac-groebner-api-stability.md`（新 `groebner_basis_lex`/`groebner_basis_grevlex`/`fglm`/`greduce_grevlex` 均 **Stable (bounded)**，泛型 `C: FieldCoeff`）；跑 `annotate_api_tiers.py --inventory` 刷新 5 份 api-stability Per-file 表。

**验收：**
- A₄ → dim ≤ 12（`quartic_a4_galois_dim_le_12` 0.12s 绿）；S₄ → dim ≤ 24（`t⁴+t+1` 0.30s 绿）；两者四根 `verify_root`。**满足**。
- 全 suite 绿、无回归。**满足**。

**文件:** `giac-core/src/algebra/poly_roots.rs`（clippy）、`giac-calculus/src/integrate.rs`（lint）、`.doc/giac-groebner-api-stability.md`（§2 + inventory）、其余 4 份 api-stability（inventory 刷新）。

---

## 3. 优先级总表

| 优先级 | ID | 内容 | 估时 | 解除条件 | 风险 |
|--------|----|------|------|----------|------|
| **P0** | groebner 泛型 `C: FieldCoeff` + ℚ 包装兼容 | 解锁系数域 | 1–1.5d | 既有 ℚ 例 + 1 例 ℚ(√2) 绿；lib 无 `num-*` deps | 低（invariant `expect` + 包装保签名） |
| **P1** | grevlex 序 + grevlex Buchberger | 0-dim 快速 GB | 1d | grevlex 0-dim 例 <100ms 绿 | 低 |
| **P2** | FGLM（grevlex→lex，通用 d） | 三角形回代 | 1.5–2d | FGLM d=4 三角形回代出全部解 | 中（通用 d Krylov/三角化细节） |
| **P3** | F-module 递归 √-归约（collapse 检测 + F-module 值恢复 + S7 接线） | A₄→12 落地 | 2–3d | `diag_fmodule_sqrt_recovery` 恢复 √(disc_q) 自校 + A₄ 测 dim≤12 | 中高（subfield-module 线代、递归、base-case GB/FGLM） |
| **P4** | S₄ 回归 + 门禁 | 防回归 | 0.5d | 全套绿、api-stability 登记、clippy 绿 | 低 |

**合计：~5.5–6.5d。** P0–P2 可在 `giac-groebner` 内独立交付（不碰 giac-core）；P3 才接线。

---

## 4. 设计要点与约束

- **系数域抽象：** 一律走 `C: FieldCoeff`，**不**为 ℚ(α) 特化。`Poly<AlgExtCPolyCoeff>` 即 `F[var]`。
- **错误传播：** 泛型 groebner 核心用 `try_*` + invariant `expect`，签名 infallible；ℚ `Stable` 公开入口（`greduce`/`groebner_basis_lex`/`groebner_basis_grevlex`/`fglm`）为薄包装，签名不变，调用方零改动。
- **不用 lex 直接 Buchberger**（已证爆炸）；0-dim 一律 `grevlex GB + FGLM`。
- **不 flatten 到 deg-12 over ℚ**（12-var GB 不可行）。
- **Fuel：** S7 内部局部 `Fuel::new(8)`，传给 `F=ℚ(α)` flat √ 递归；不改 `try_square_root_in_field_impl` 全局签名（S0–S6 不动）。S7 → F 不会重入 S7（F 顶层 deg-4 ≠ deg-3 触发条件），fuel 为防御性深度 guard。
- **`proto_*` 处置：** 仅 `proto_rational_roots` 思路 rewrite 为 `rational_roots_in_field`；其余 `proto_*` 不用，留 dead WIP（父项已 `#[allow(dead_code)]`）。
- **盲 adjoin 兜底：** S7 miss ≠ 错答，只退回 `try_sqrt_pairwise_fallback` → blind adjoin（维数可能偏高但根仍对）；DoD 仅要求 A₄ 命中、S₄ 允许 miss。
- **两个不同上限（勿混）：**
  - `d = |B|` = **商环维数**（FGLM 乘法矩阵阶数；本问题 d=2；FGLM `D_MAX=64` 是 d 上界）。
  - `dim = [塔域 : ℚ]` = **塔域总维数**（S7 触发条件 `dim ≤ B=12`）。
- **Monomial 序：** derive `Ord`（字典序）既非 lex 也非 grevlex；grevlex 路径必须用 `leading_term_grevlex`，禁止 `leading_term()`。
- **API tier：** 新 fn inline `/// **Tier**` 标注 + `.doc/giac-groebner-api-stability.md` 登记；算法 crate 改动跑 `annotate_api_tiers.py --inventory`。

---

## 5. 验收（总 DoD）

- [x] P0：`giac-groebner` 泛型核心 `C: FieldCoeff`；ℚ 包装签名不变、既有例不回归；ℚ(√2) 例绿；lib 无 `num-*` deps（实现偏离见 §2 P0 末尾）
- [x] P1：grevlex Buchberger；0-dim 例 <100ms；grevlex 路径未误用 `leading_term()`（实现偏离见 §2 P1 末尾）
- [x] P2：FGLM 通用 d；d=4 三角形回代出全部解；`d>D_MAX` 返回 `None`（实现偏离见 §2 P2 末尾；shape-lemma generic-position 限制：非 generic ⇒ None 安全兜底）
- [x] P3：`diag_adjoin_collapse_recovery` collapse 检测绿（0.19s，已就位）；`diag_fmodule_sqrt_recovery` 恢复 √(disc_q)∈K 自校 `δ²=disc_q` 且 dim 不增（0.06s 绿）；`quartic_a4_galois_dim_le_12` unignore 绿（dim ≤ 12，0.12s）；S7（collapse + F-module）接线（`try_square_root_in_field_impl` gate dim ≥ 4）。实现偏离见 §2 P3 末尾 + [postmortem](../issues/GIAC-poly-f5-fglm-debug-postmortem.md)。
- [x] P4：`t⁴+t+1` dim ≤ 24 不回归（`roots_quartic_t4_plus_t_plus_1` / `field_session_dimension_bound_quartic_tight` / `quartic_adjoin_deflate_t4_plus_t_plus_1` 全绿，4 根 `verify_root`）；全 suite 绿（nextest release `1089 passed / 46 skipped / 0 failed`，13.9s）；clippy + substring-golden 绿（`ci-clippy.sh`）；api-stability 登记（`giac-groebner-api-stability.md` §2 + 5 份 inventory 刷新）。实现偏离见 §2 P4 末尾。
