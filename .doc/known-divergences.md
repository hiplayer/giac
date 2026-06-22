# 已知偏离记录（Rust vs giac golden）

Rust 实现与 giac C++ golden **字面不一致**但可能数学等价，或 **故意修正 giac bug** 的条目。格式见 [`rust-migration-plan.md` §6.4](rust-migration-plan.md)。

**状态：** `open` = 待处理；`accepted` = 已验证并更新 golden/测试；`wontfix-giac` = giac 侧问题，Rust 不跟随。

---

## 模板

```markdown
### DIV-NNN: 简短标题

- **状态:** open | accepted | wontfix-giac
- **输入:** `表达式;`
- **giac 输出:** `...`
- **Rust 输出:** `...`
- **归类:** bug 修复 | 等价不同形 | 真语义分歧
- **理由:** …
- **验证:** expand(sub(a,b))=0 / 手算 / 文献
- **测试处理:** 更新 golden | assert_equiv | 跳过
```

---

## 已记录条目

### DIV-001: factor 大指数有理式 segfault

- **状态:** wontfix-giac
- **输入:** `factor((x^202+x^101+1)/(x^2+x+1));`（前置 `cas_setup(0,0,0,1,0,[1e-10,1e-17],12,[1,50,0,25],0,0,0),xcas_mode(0);`）
- **giac 输出:** 进程 **SIGSEGV**（`giac_check_factor_segfault_known`, WILL_FAIL）
- **Rust 输出:** `Err(EvalError::…)` 或有限时间内返回正确因子分解
- **归类:** bug 修复
- **理由:** giac 未定义行为；Rust 不得 panic/segfault
- **验证:** 数学上因子分解存在且唯一（代数验证）
- **测试处理:** Rust 专项测试断言 `Result::Err` 或正确 `Ok`；不纳入 giac golden diff

---

## 预期常见偏离（尚未触发，预留）

| ID | 领域 | 说明 | 策略 |
|----|------|------|------|
| DIV-010 | 化简 | `x+1` vs `1+x` 排序 | `canonicalize` 后 diff |
| DIV-020 | 积分 | 不同但等价的原函数 | `assert_equiv` + diff 还原 |
| DIV-030 | 求解 | 三角方程多解分支 | 集合等价或更新 golden |
| DIV-076–083 | 代数扩展 / `poly_roots` | `i`、√ 分支、二次/双二次/四次根形、塔顺序、列表序、显示 | 见下文 §Phase B；`assert_equiv` / 集合等价 |
| DIV-040 | 浮点 | `evalf` 末位差异 | `cas_floats` 10 位有效数字 |
| DIV-050 | 随机 | `randNorm` 序列不同 | 固定 seed 或只验类型/范围 |
| DIV-060 | Groebner | `gbasis` 不可用 | MVP 不实现；`greduce` 自研路径 |

---

## MVP 实施限制（giac-rs）

与 giac 默认上界对齐的显式边界；超出返回 `TypeError` 而非未定义行为。

| 常量 | 值 | 范围 | giac 参考 | Rust 位置 |
|------|-----|------|-----------|-----------|
| `MAX_POLY_EXPONENT` | `1024` | `expr_to_poly` 中 `Pow` 指数 ≤ 1024 | `series.cc` 级数阶、`GBASISF4_MAX_TOTALDEG` | `giac-core::limits` |

说明：`expr_to_poly` 用 `u64` 解析并校验 `<= MAX_POLY_EXPONENT`；`giac-poly` 单项式指数同为 `u64`。`BigInt::pow` 仍仅接受 `u32`，超大指数经 `giac-poly::exp::bigint_pow` 桥接。

---

## Phase 3 线性代数（DIV-070–073）

### DIV-070: 3×3 SVD 有理数近似与重构

- **状态:** accepted
- **输入:** `svd([[1,2,1],[3,4,1],[1,5,6]]);`
- **giac 输出:** 浮点 `U, Σ, Vᵀ` 三元组（giac 对符号矩阵给出数值警告）
- **Rust 输出:** 经 `f64_to_rational` / `cas_floats` 规范化后的 `U, Σ, Vᵀ`；大分母有理数以 10 位有效数字浮点打印
- **归类:** 等价不同形
- **理由:** 固定分母 `1_000_000` 有理数近似无法稳定重构 `A`；改为 Stern-Brocot 重构 + SymPy 验证 `U·Σ·Vᵀ ≈ A`（容差 ≤ 1e-5）
- **验证:** `sympy_verify.py` SVD 重构；`giac-linalg` 3×3 单元测试
- **测试处理:** `phase3_triple::test_linalg_decomp_triple` 不再跳过 3×3 SVD；`is_known_numeric_gap` 已移除

---

### DIV-071: `gramschmidt` 内积 / 开方 / 输出形

- **状态:** accepted
- **输入:** `gramschmidt([1,1+x],(p,q)->integrate(p*q,x,-1,1));`
- **giac 输出:** `[1/sqrt(2),(1+x-1)/(sqrt(6)/3)]`（等价于 `[(sqrt(2))^-1, x*((sqrt(2/3))^-1)]`）
- **Rust 输出:** `[(sqrt(2))^-1,x*((sqrt(2/3))^-1)]`
- **归类:** 等价不同形
- **理由:** 定积分边界、lambda 保持符号、多项式积展开等路径修复后，正交归一性成立但打印形式与 giac 不同
- **验证:** SymPy 属性：`⟨v_i,v_j⟩=0`（i≠j）、`⟨v_i,v_i⟩=1`
- **测试处理:** `phase3_triple` / `sympy_all` 已移除 `gramschmidt` skip

---

### DIV-072: `egv` / `jordan` 3×3 与一般情形

- **状态:** accepted（MVP 子集）
- **输入:** `egv([[4,1,-2],[1,2,-1],[2,1,0]]);`、`jordan([[1,1],[0,1]]);`、`jordan([[1,1,0],[0,1,0],[0,0,2]]);`
- **giac 输出:** `"Not diagonalizable at eigenvalue 2"`；Jordan 对 `(J,P)` 与输入矩阵相关（giac 约定）
- **Rust 输出:** 不可对角化时同信息字符串；已是 Jordan 形时 `(J=A, P=I)`，满足 `P⁻¹JP=A`；整数 3×3 可对角化时返回特征值序列
- **归类:** 等价不同形（Jordan 规范形选取）+ MVP 范围限制
- **理由:** 整数特征多项式有理根 + 几何重数判定覆盖 conformance 脚本；一般符号 3×3 / 完整 Jordan 基未实现
- **验证:** SymPy：`det(A-λI)=0`（`egv`）；`P⁻¹JP ≈ A`（`jordan`）
- **测试处理:** `test_linalg_ext` / `phase3_triple` 已移除 `egv`/`jordan` skip；一般情形返回 `NotImplemented` 不 panic

---

### DIV-073: `gauss` 合同变换 vs 字面等价

- **状态:** accepted
- **输入:** `gauss(x^2+y^2-1,[x,y]);`、`gauss(x*y,[x,y]);`、`gauss(x^2-y^2,[x,y]);`（`bin/test_gauss_ext`）
- **giac 输出:** 合同对角形（与输入二次型合同，未必逐字相等）
- **Rust 输出:** 合同对角形；与 giac 在部分用例上字面不同
- **归类:** 等价不同形
- **理由:** `gauss` 实现合同变换 `PᵀAP`，非 `normal(sub(a,b))=0` 字面相等；giac 与 Rust 对角形选取可不同但二次型等价
- **验证:** SymPy 弱检查（随机代入验良定义）；`assert_equiv` 用于其他代数输出
- **测试处理:** `phase3_triple::test_gauss_ext_triple` 验 SymPy 属性，不要求 giac-rs 与 giac 字面一致

---

### DIV-075: `texpand(cos(3*x))` 幂次括号

- **状态:** accepted
- **输入:** `texpand(cos(3*x));`
- **giac 输出:** `4*cos(x)^3-3*cos(x)`
- **Rust 输出:** `4*((cos(x))^3)-3*cos(x)`（display 括号差异）
- **归类:** 等价不同形
- **理由:** 数学等价；`expand` 保留 `Pow` 时 display 加括号
- **验证:** SymPy `expand_trig` 属性验证
- **测试处理:** `trig_format_diff`；SymPy 验证不要求字面一致

---

## Phase B 代数扩展 / `poly_roots`（DIV-076–083）

与 giac C++ `alg_ext.cc` / `solve.cc` / `rootof` 在 **符号约定** 上的系统性差异。数学上多解集等价，但 **min_poly、coords、打印、根列表顺序** 常与 golden 字面不同。背景见 [GIAC-algext-adoption.md](issues/GIAC-algext-adoption.md) §7–8、[giac-tower-common-math.md](giac-tower-common-math.md) §4。

### DIV-076: 虚数单位 `i` 的表示（形式 `AlgExtC.im` vs `cst_i` vs 塔 adjoin）

- **状态:** open（过渡态；目标态 adjoin \(v^2+1\) 见 adoption §8.2）
- **输入:** `sqrt(-2)`；`rootof([1,0,1])`；`solve(t^2+1=0,t)`；扩域内 `sqrt(-u)`（\(u>0\) 在 \(K\) 中）
- **giac 输出:** 上下文符号 `cst_i`；`rootof([1,0,1])` **化为** `i`；负判别式 `cst_i*sqrt(-delta)`；`common_EXT` 对 `sqrt(X)` 与 `sqrt(-X)` 合并时 adjoin `i`（`alg_ext.cc`）
- **Rust 输出:** `FieldSession::mul_formal_i` — \(i\cdot z\) 写入 `AlgExtC.im`，**不把** \(v^2+1\) adjoin 进塔；`algext_sqrt_branches` 对 \(u<0\) 返回 `Complex(0, AlgExt)`；`sqrt_principal` / `sqrt_disc` 对 ℚ 负元 `adjoin_sqrt(|u|)` 再 `mul_formal_i`
- **归类:** 等价不同形（过渡态）+ 目标态与 giac 仍可能字面不同
- **理由:** Rust 显式塔 + `AlgExtC` 一等类型；giac 混用符号 `i`、`_EXT` 与 `common_EXT` 特判；过渡态未实现 adoption 目标 \(K_3=K_2(i)\)
- **验证:** `(i·√2)² = -2` 的 `eq_mod`；`canonicalize_to_algext_c` 往返；与 giac 比 **集合/等价**，不比 `i` 是否进 min_poly
- **测试处理:** `poly_roots` / `field_session` 单元测 `eq_mod`；conformance 用 `assert_equiv`；不追 giac `rootof([1,0,1])→i` 字面

---

### DIV-077: 平方根主分支（`sqrt` / `normalize_sqrt` vs `adjoin_sqrt` / `pick_sqrt_euler`）

- **状态:** open
- **输入:** `solve(x^2-2=0,x)`；`sqrt(Δ)`（二次/三次 resolvent）；Euler 四次中 `pick_sqrt_euler(α,β,γ)`
- **giac 输出:** `fastsign` / `complex_mode` / `normalize_sqrt(sqrt(·))` 隐式选支；二次 `delta_prime=sqrt(delta_prime)` 或 `i*sqrt(-delta_prime)`（`solve.cc` `in_solve`）
- **Rust 输出:** `adjoin_sqrt` → `algext_square_roots` 取 **一个** 生成元（`pop()`）；`pick_sqrt_euler` 对 \(K\) 中可能为负的元用 `approx_real_sign`（deg≤6 数值嵌入）在实 adjoin 与 `i·√|u|` 间切换；无 giac 式 `normalize_sqrt`
- **归类:** 等价不同形
- **理由:** giac 分支绑定符号 `sqrt` 化简链；Rust 分支绑定 **新 adjoin 元** 及其符号翻转；同一 \(\sqrt{u}\) 的 coords 可差整体负号或不同塔层
- **验证:** `beta^2 eq_mod u`；代入原多项式 `root_vanishes`；必要时 `evalf` 同象限
- **测试处理:** `assert_equiv` / 零化子验证；golden 字面 diff 记本条；PR-F1 `sqrt_in_field` 落地后复审是否收窄

---

### DIV-078: 二次根外层形状（`rootof([±1,0],P)` vs Vieta vs 公式根）

- **状态:** open
- **输入:** `solve(t^2-2=0,t)`；`roots(x^2+b*x+c,x)` over ℚ
- **giac 输出:** `solve(vecteur)`：`algebraic_EXTension([1,0,-b/2], u²-Δ)` 与 `([-1,0,-b/2], …)`，或 adjoin 整式 `v` 时 `([1,0,0],v)` / `([-1,-b,0],v)`（`solve.cc` 1214–1238）；`rootof([1,0],P)` / `rootof([-1,0],P)`
- **Rust 输出:** **三条路径并存：** (1) `poly_roots::quadratic_roots` — adjoin \(x^2+bx+c\)，\((\alpha,-b-\alpha)\)；(2) `quadratic_roots_formula` — \((-b\pm\sqrt\Delta)/2\)；(3) `giac-solve::quadratic_rootof_roots` — `rootof([1,0],P)` 与 `[-1,0]`
- **归类:** 等价不同形
- **理由:** 同一方程可选 **最小 \(\Delta\)-扩张**、**原式扩张** 或 **公式根**；min_poly 与 coords 无统一 canonical 形
- **验证:** 两根均 `root_vanishes`；Vieta \(r_1+r_2=-b\)，\(r_1 r_2=c\)
- **测试处理:** P4 统一 `solve` 管线后收敛为一种对外形；此前各路径单测 `eq_mod`，golden 不要求字面一致

---

### DIV-079: 双二次 \(t^4+bt^2+c\) 四根（实 `algext_square_roots` vs `Complex(0,·)`）

- **状态:** open
- **输入:** `solve(t^4-2=0,t)`；`biquadratic_rootof_roots` / `poly_roots::biquadratic_roots`
- **giac 输出:** 经 `solve` / `rootof` + `sqrt` 链；负 \(u=t^2\) 时常为 `i*sqrt(...)` 形状；adoption 文档期望四根 \(\pm 2^{1/4}\)、\(\pm i\cdot 2^{1/4}\) 为 `AlgExtC`（§8.2）
- **Rust 输出:** `poly_roots::biquadratic_roots` → `algext_c_sqrt`（仅 **im=0**）→ 四个 **实部** 根（如 \(t^4-2\) 得 \(\pm 2^{1/4}\)）；`giac-solve::biquadratic_rootof_roots` → `algext_sqrt_branches` → \(u<0\) 时 `Complex(0,±β)`。**crate 内两路径不一致**
- **归类:** 等价不同形 + Rust 内部未统一
- **理由:** `algext_c_sqrt` 不做负扩域元的虚分支；`algext_sqrt_branches` 用 `Complex` 补丁且不 adjoin \(i\) 进塔
- **验证:** 四根均零化 \(t^4+bt^2+c\)；集合 \(\{r\}\) 与 giac 四解 `assert_equiv`（若 giac 给出四根）
- **测试处理:** `biquadratic_t_fourth_minus_two` / `roots_biquadratic_t4_minus_2` 验零化与个数；合并路径时更新本条

---

### DIV-080: 一般四次方程（Euler + resolvent vs giac 单根 `algebraic_EXTension`）

- **状态:** open（`poly_roots` 四次 e2e 仍 `#[ignore]`，塔 `adjoin_sqrt` 嵌入问题）
- **输入:** `solve(t^4+t+1=0,t)`；一般 depressed quartic + Ferrari resolvent
- **giac 输出:** `solve(vecteur v)` 在 **deg>2**（除三次特判）时通常只返回 **一个** `algebraic_EXTension([1,0], v)`（`solve.cc` 1241–1243）；无与 Rust 对等的「resolvent 三次 + Euler 四根 + ε 符号表」管线
- **Rust 输出:** `quartic_roots`：depress → resolvent \(R(z)=z^3-pz^2-4rz+(4pr-q^2)\) → `one_cubic_root` + `quadratic_roots_formula` → `euler_depressed_quartic_roots`（6 置换 × flip 搜索）；约定 \(\sqrt\gamma=-q/(\sqrt\alpha\sqrt\beta)\)；四符号 \((\varepsilon_1,\varepsilon_2,\varepsilon_3)\in\{(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)\}\)，\(\varepsilon_1\varepsilon_2\varepsilon_3=+1\)
- **归类:** 真语义分歧（**解的个数与算法**）+ Rust 自研规格（四根零化）
- **理由:** giac 对一般四次 **不承诺** 四根 closed form 列表；Rust 以 **四根均在 \(K(\cdots)\) 且零化** 为测试规格，不以 giac 字面为准
- **验证:** 四根 `root_vanishes`；`euler_gamma_relation_holds`（\(\sqrt\alpha\sqrt\beta\sqrt\gamma+q=0\)）；与 giac 仅比 **是否存在代数解**，不比四根 coords
- **测试处理:** `euler_four_roots_vanish`、`roots_quartic_t4_plus_t_plus_1`（ignore 至 F1/F2 塔修复）；不纳入 giac golden 四根字面 diff

---

### DIV-081: 扩域塔构造顺序与 `common_EXT`（adjoin 顺序 / `sqrt_in_field`）

- **状态:** open
- **输入:** resolvent 后连续 `adjoin_sqrt`；`align_elements` 跨会话；`common(√2,√3)` 类合并
- **giac 输出:** `algebraic_EXTension` + `common_EXT`（含 `sqrt(X)`/`sqrt(-X)`、`sqrt` 同对象分数幂等特判，`alg_ext.cc` 910–940）；`convert_rootof` / `keep_algext` 控制是否保持符号 `rootof`
- **Rust 输出:** 显式 `ExtensionField` 塔，每次 `adjoin_sqrt` 叠层；**尚无** `sqrt_in_field`（PR-F1）；resolvent 后再 `adjoin_sqrt` 可破坏 \(\varepsilon^2=u\)（`adjoin_sqrt_after_resolvent_deflate_only` 回归）；T4a `compute_common_tower` 与 flatten **coords 可不同**（[giac-tower-common-math.md](giac-tower-common-math.md) §4.2）
- **归类:** 等价不同形（同域 `eq_mod`）+ 实现缺口（错误 adjoin）
- **理由:** 构造顺序不同 → 同一数的 coords 不同；重复 adjoin 未识别已有平方根时 Rust 行为错误（非 giac 约定问题）
- **验证:** 同域 `element_eq_mod`；`min_poly_over_q` 次数；塔计划 S0 逆序 / `tower_common_invariants_*`
- **测试处理:** `adjoin_sqrt_squares_one_resolvent_root` 等回归；F1/F2 完成后复审本条「实现缺口」是否关闭

---

### DIV-082: 代数根列表顺序与去重

- **状态:** accepted
- **输入:** `solve(t^4-2=0,t)`；`roots` / `poly_algext_roots` 任意次数
- **giac 输出:** `solve` 合并顺序依赖因子分解、`mergevecteur`、`isolate_mode`、分支构造顺序
- **Rust 输出:** `dedup_roots` 按 **插入顺序** 保留首个；Euler 搜索 **先命中** 的置换/flip 即返回；`quadratic_roots` 固定 `[r1, -b-r1]`
- **归类:** 等价不同形
- **理由:** 无 canonical 排序约定；多解集为集合而非有序列表
- **验证:** 集合等价：每个 giac 根与某 Rust 根 `eq_mod`；个数一致
- **测试处理:** conformance 用集合比较或 `assert_equiv` 排序无关；不要求列表下标对齐

---

### DIV-083: 显示层（`rootof` / `AlgExt` / `AlgExtC` / `Complex`）

- **状态:** open
- **输入:** 纯实代数数；复代数数；`to_rootof_expr` 往返
- **giac 输出:** `rootof([c…], poly1[…])` 未求值或经 `normal`；`rootof([1,0,1])→i`；复数 `re+im*i`
- **Rust 输出:** `im=0` → `AlgExt` / `rootof` 显示；`im≠0` → `AlgExtC` 或过渡态 `Complex(0, AlgExt)`；**不**自动把 `rootof([1,0,1])` 折叠为 `i`
- **归类:** 等价不同形
- **理由:** 表示层策略不同；Rust 优先塔顶 coords + `AlgExtC` 规范形（adoption §8.3）
- **验证:** `to_rootof_expr` / `try_as_algext_data` 往返 `eq_mod`；打印 diff 记 `display` 测而非数学失败
- **测试处理:** `algext_to_rootof_roundtrip_display`；golden 字符串 diff 引用本条

---

## 维护规则

1. conformance 失败时先归类，再决定修 Rust / 更新 golden / 记本条。
2. 每条必须有 **验证方式**；不能只写「Rust 更对」。
3. `accepted` 后同步更新 Rust 测试或 golden 子集。
