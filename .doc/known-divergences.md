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

## 维护规则

1. conformance 失败时先归类，再决定修 Rust / 更新 golden / 记本条。
2. 每条必须有 **验证方式**；不能只写「Rust 更对」。
3. `accepted` 后同步更新 Rust 测试或 golden 子集。
