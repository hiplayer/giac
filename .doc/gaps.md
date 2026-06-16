# 测试覆盖缺口

基于当前 **53** 个 CTest 目标与 giac-1.5.0 源码模块对比，以下功能域**无专项测试或覆盖薄弱**。

## 未纳入构建 / 无测试

| 模块 | 说明 |
|------|------|
| `sparse.cc` | 稀疏矩阵，已从部分构建排除 |
| `cocoa.cc` | 依赖 CoCoA 外部库 |
| `plot3d.cc` | 3D 绘图，无 GUI 难触发 |
| `maple.cc` / `ti89.cc` | Maple/TI 兼容模式 |
| `rpn.cc` | RPN 逆波兰模式 |
| `pari.cc` | Pari/GP 接口 |

## 有源码但测试极少

| 模块 | 当前覆盖 | 缺口 |
|------|---------|------|
| `help.cc` | 间接 | 无 help 查询测试 |
| `icas.cc` / GUI 相关 | 低 | 无 FLTK/GUI 测试 |
| `moyal.cc` | ~5% | 仅 `test_special` 中 moyal/Beta 等 |
| `misc.cc` | 分散 | 特殊函数不全 |
| `prog.cc` | `test_subst` | 无 for/while 控制流 |
| `global.cc` | 间接 | 无 cas_setup 变体测试 |

## 已有 API 但未单独成测

以下 API 仅出现在 `flanex` 等综合用例中，无独立回归：

- `gbasis` / `greduce`（Groebner，且 flanex 中可能报错）
- `gramschmidt`、`potential`、`vpotential`（向量分析）
- `lu`、`qr`、`svd`（矩阵分解，flanex 中有）
- `partfrac`（除 `check/testpartfrac` 单条外）
- `evalf`、`float2rational`（数值计算）

## 建议补充方向

| 方向 | 状态 | 新增测试 |
|------|------|----------|
| check 套件扩展 | 部分完成 | `giac_check_partfrac_ext`（4 条 partfrac 黄金 diff） |
| 模运算 | 已补充 | `bin/test_modular` |
| 数值求解 | 已补充 | `bin/test_numerical` |
| 概率统计 | 已补充 | `bin/test_probstat` |
| 矩阵分解 | 已补充 | `bin/test_linalg_decomp`（lu/qr/svd/gramschmidt） |
| 向量分析 | 已补充 | `bin/test_vector_calc` |
| Groebner | 部分完成 | `bin/test_groebner`（greduce；gbasis 仍不可用） |
| 程序控制流 | 已补充 | `bin/test_prog`（for/while/sum/product） |

仍待补充：`chk_fhan*` 系列、`gbasis`、Maple/RPN 兼容模式、GUI 测试。

