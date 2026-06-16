# 功能覆盖分析

按 CAS 功能域列出测试覆盖的 API、关联源码模块与测试来源。

## 微积分
**源码模块：** `intg.cc`, `derive.cc`, `risch.cc`, `series.cc`, `intgab.cc`
**覆盖 API（23/23）：** `integrate`, `int`, `diff`, `derive`, `limit`, `series`, `taylor`, `partfrac`, `risch`, `ibp`, `ibpu`, `ibpdv`, `proot`, `rootof`, `laplacian`, `hessian`, `divergence`, `curl`, `potential`, `vpotential`, `preval`, `fxnd`, `simplify`
**测试来源：**
- `bin/input`
- `bin/test_diff`
- `bin/test_factor`
- `bin/test_integrate`
- `bin/test_integrate_ext`
- `bin/test_integrate_more`
- `bin/test_limit`
- `bin/test_series`
- `check/flanex`
- `check/testcas`
- `check/testintegrate`
- `check/testlimit`
- `check/testother`
- `check/testpartfrac`

## 微分方程
**源码模块：** `desolve.cc`
**覆盖 API（1/1）：** `desolve`
**测试来源：**
- `bin/test_desolve`
- `bin/test_desolve_ext`

## 方程求解
**源码模块：** `solve.cc`
**覆盖 API（4/4）：** `solve`, `linsolve`, `froots`, `realroot`
**测试来源：**
- `bin/3`
- `bin/test_linalg`
- `bin/test_solve`
- `bin/test_solve_ext`
- `bin/test_sturm_ext`
- `check/flanex`
- `check/testcas`

## 多项式与因式分解
**源码模块：** `sym2poly.cc`, `modpoly.cc`, `ezgcd.cc`, `gausspol.cc`, `ifactor.cc`, `gauss.cc`
**覆盖 API（33/33）：** `factor`, `ifactors`, `factors`, `divis`, `content`, `quo`, `rem`, `irem`, `iquo`, `idivis`, `gcd`, `lcm`, `egcd`, `abcuv`, `simp2`, `horner`, `resultant`, `roots`, `peval`, `e2r`, `r2e`, `valuation`, `divpc`, `ptayl`, `fcoeffs`, `propfrac`, `chinrem`, `cyclotomic`, `canonical_form`, `gbasis`, `greduce`, `partfrac`, `float2rational`
**测试来源：**
- `bin/test_cas_basic`
- `bin/test_factor`
- `bin/test_ifactor`
- `bin/test_ifactor_ext`
- `bin/test_poly`
- `bin/test_poly_ext`
- `check/flanex`
- `check/testcas`
- `check/testfactor`
- `check/testnormalize`
- `check/testpartfrac`

## 线性代数
**源码模块：** `lin.cc`, `vecteur.cc`, `gauss.cc`
**覆盖 API（32/32）：** `rref`, `det`, `tran`, `trn`, `ker`, `image`, `jordan`, `egv`, `egvl`, `pcar`, `charpoly`, `linsolve`, `svd`, `lu`, `qr`, `rank`, `cross`, `dot`, `idn`, `ranm`, `hadamard`, `hilbert`, `vandermonde`, `adjoint_matrix`, `changebase`, `pmin`, `signature`, `gramschmidt`, `syst2mat`, `axq`, `qxa`, `gauss`
**测试来源：**
- `bin/test_gauss_ext`
- `bin/test_linalg`
- `bin/test_linalg_ext`
- `bin/test_poly`
- `check/flanex`
- `check/testcas`

## 三角与双曲
**源码模块：** `usual.cc`
**覆盖 API（17/17）：** `texpand`, `tlin`, `halftan`, `lin`, `tsimplify`, `tcollect`, `lncollect`, `trig2exp`, `trigcos`, `trigsin`, `tan2sincos`, `acos2asin`, `asin2acos`, `atan2asin`, `asin2atan`, `acos2atan`, `atan2acos`
**测试来源：**
- `bin/test_trig`
- `check/flanex`
- `check/testcas`

## 拉普拉斯与傅里叶
**源码模块：** `intgab.cc`, `misc.cc`
**覆盖 API（4/4）：** `laplace`, `ilaplace`, `fourier_an`, `fourier_cn`
**测试来源：**
- `bin/test_laplace`
- `check/flanex`

## 数论
**源码模块：** `ifactor.cc`
**覆盖 API（8/8）：** `isprime`, `nextprime`, `prevprime`, `euler`, `divisors`, `nCr`, `nPr`, `bernoulli`
**测试来源：**
- `bin/test_ifactor`
- `bin/test_ifactor_ext`
- `bin/test_permu`
- `bin/test_permu_ext`
- `check/flanex`
- `check/testcas`

## 排列与置换
**源码模块：** `permu.cc`
**覆盖 API（8/8）：** `permu2cycles`, `cycles2permu`, `permuorder`, `permu2mat`, `is_permu`, `est_permu`, `est_cycle`, `p1op2`
**测试来源：**
- `bin/test_permu`
- `bin/test_permu_ext`
- `check/flanex`

## Sturm 根隔离
**源码模块：** `csturm.cc`
**覆盖 API（2/2）：** `sturm`, `sturmab`
**测试来源：**
- `bin/test_sturm`
- `bin/test_sturm_ext`
- `check/flanex`

## 几何与绘图
**源码模块：** `plot.cc`
**覆盖 API（3/9）：** `plotfunc`, `triangle`, `xyztrange`
**测试来源：**
- `check/flanex`
- `check/testgeo`

## 几何同构
**源码模块：** `isom.cc`
**覆盖 API（1/1）：** `mkisom`
**测试来源：**
- `bin/test_isom`
- `check/flanex`

## 特殊函数
**源码模块：** `moyal.cc`, `misc.cc`
**覆盖 API（14/14）：** `gamma`, `zeta`, `psi`, `Beta`, `polygamma`, `UTPN`, `legendre`, `hermite`, `laguerre`, `tchebyshev1`, `tchebyshev2`, `moyal`, `randNorm`, `randchisquare`
**测试来源：**
- `bin/test_special`
- `check/flanex`

## 四元数
**源码模块：** `quater.cc`
**覆盖 API（1/1）：** `quaternion`
**测试来源：**
- `bin/test_special`

## 化简与规范化
**源码模块：** `usual.cc`, `misc.cc`
**覆盖 API（6/6）：** `normal`, `ratnormal`, `non_recursive_normal`, `epsilon2zero`, `truncate`, `reorder`
**测试来源：**
- `bin/test_factor`
- `check/flanex`
- `check/testcas`
- `check/testgeo`
- `check/testintegrate`
- `check/testnormalize`

## 复数与基础
**源码模块：** `usual.cc`, `gen.cc`
**覆盖 API（6/6）：** `abs`, `arg`, `conj`, `re`, `im`, `sign`
**测试来源：**
- `bin/test_cas_basic`
- `check/testcas`

## 公式导出
**源码模块：** `tex.cc`, `mathml.cc`
**覆盖 API（2/2）：** `latex`, `mathml`
**测试来源：**
- `bin/test_export`

## 替换与程序
**源码模块：** `subst.cc`, `prog.cc`
**覆盖 API（3/3）：** `subst`, `sum`, `product`
**测试来源：**
- `bin/test_subst`
- `check/flanex`
- `check/testcas`

## 模运算
**源码模块：** `modpoly.cc`
**覆盖 API（1/1）：** `chinrem`
**测试来源：**
- `check/flanex`

## 其他（仅少量用例）

以下 API 出现在测试中，但未归入上述主功能域：

- `acsc` — check/flanex
- `append` — check/flanex
- `asec` — check/flanex
- `assume` — check/flanex, check/testgeo, check/testintegrate
- `cas_setup` — check/testfactor, check/testgeo, check/testintegrate
- `coeff` — check/testcas
- `color` — check/testgeo
- `cot` — check/flanex
- `csc` — check/flanex
- `degree` — check/flanex
- `equal` — check/testcas
- `eval` — check/testcas
- `evalf` — check/flanex, check/testcas
- `exlr` — check/flanex
- `exp2pow` — check/flanex
- `f` — bin/test_subst
- `fdistrib` — check/testcas
- `ichinrem` — check/testcas
- `inv` — check/testcas
- `is_element` — check/testgeo
- `is_prime` — check/testcas
- `jacobi_symbol` — check/testcas
- `lagrange` — check/flanex
- `ln` — bin/1, bin/2
- `lname` — check/flanex, check/testcas
- `lvar` — check/flanex, check/testcas
- `nodisp` — check/testgeo
- `pcoeff` — check/testcas
- `purge` — check/flanex, check/testgeo
- `quote` — check/testcas
- `sec` — check/flanex
- `size` — check/flanex
- `smod` — check/testcas
- `sort` — check/testcas
- `sto` — check/testcas
- `suppress` — check/flanex
- `tan` — check/flanex
- `trace` — check/flanex

