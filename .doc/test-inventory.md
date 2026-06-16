# 测试用例清单

## bin/ 覆盖驱动脚本

| 文件 | CTest | 表达式数 | API 数 | 说明 |
|------|-------|---------|--------|------|
| `1` | `giac_bin_1` | 5 | 1 | 对数化简冒烟 |
| `2` | `giac_bin_2` | 5 | 1 | 对数化简冒烟 |
| `3` | `giac_bin_3` | 1 | 1 | 超越方程冒烟 |
| `input` | `giac_bin_input` | 5 | 2 | 积分化简冒烟 |
| `test_cas_basic` | `giac_bin_test_cas_basic` | 3 | 3 | abs/gcd/conj 等基础 |
| `test_desolve` | `giac_bin_test_desolve` | 3 | 1 | ODE desolve |
| `test_desolve_ext` | `giac_bin_test_desolve_ext` | 3 | 1 | 重根/非齐次/强迫 ODE |
| `test_diff` | `giac_bin_test_diff` | 2 | 1 | 单变量微分 diff |
| `test_export` | `giac_bin_test_export` | 2 | 2 | latex/mathml |
| `test_factor` | `giac_bin_test_factor` | 3 | 3 | 符号因式分解、展开 |
| `test_gauss_ext` | `giac_bin_test_gauss_ext` | 3 | 1 | gauss 多项式消元 |
| `test_ifactor` | `giac_bin_test_ifactor` | 3 | 3 | ifactors/isprime/nextprime |
| `test_ifactor_ext` | `giac_bin_test_ifactor_ext` | 7 | 7 | divisors/euler/idivis/iquo/irem |
| `test_integrate` | `giac_bin_test_integrate` | 3 | 2 | 基础不定积分 |
| `test_integrate_ext` | `giac_bin_test_integrate_ext` | 4 | 3 | 分部积分、定积分、多元 derive、proot |
| `test_integrate_more` | `giac_bin_test_integrate_more` | 5 | 3 | exp*cos、ln、risch、rootof |
| `test_isom` | `giac_bin_test_isom` | 3 | 1 | mkisom 旋转/反射 |
| `test_laplace` | `giac_bin_test_laplace` | 3 | 3 | laplace/ilaplace/fourier_an |
| `test_limit` | `giac_bin_test_limit` | 2 | 1 | 极限 |
| `test_linalg` | `giac_bin_test_linalg` | 5 | 4 | 矩阵幂、rref、linsolve、det、charpoly |
| `test_linalg_ext` | `giac_bin_test_linalg_ext` | 6 | 6 | jordan/egv/ker/image/pcar/tran |
| `test_permu` | `giac_bin_test_permu` | 3 | 3 | 置换循环、nCr |
| `test_permu_ext` | `giac_bin_test_permu_ext` | 5 | 4 | is_permu/permuorder/permu2mat/nPr |
| `test_poly` | `giac_bin_test_poly` | 5 | 5 | gcd/quo/rem/content/gauss |
| `test_poly_ext` | `giac_bin_test_poly_ext` | 7 | 7 | egcd/abcuv/simp2/lcm/horner/resultant/roots |
| `test_series` | `giac_bin_test_series` | 3 | 2 | taylor/series |
| `test_solve` | `giac_bin_test_solve` | 2 | 1 | 方程求解 |
| `test_solve_ext` | `giac_bin_test_solve_ext` | 3 | 1 | 三角方程、方程组、无理根 |
| `test_special` | `giac_bin_test_special` | 7 | 7 | moyal/Beta/polygamma/UTPN/随机分布/quaternion |
| `test_modular` | `giac_bin_test_modular` | 8 | 8 | modp/smod/chinrem/模多项式 gcd/rref |
| `test_numerical` | `giac_bin_test_numerical` | 7 | 4 | fsolve/newton/evalf/float2rational |
| `test_probstat` | `giac_bin_test_probstat` | 10 | 10 | Beta/UTPN/随机分布/gamma/zeta/psi |
| `test_linalg_decomp` | `giac_bin_test_linalg_decomp` | 6 | 5 | lu/qr/svd/gramschmidt/trace |
| `test_vector_calc` | `giac_bin_test_vector_calc` | 6 | 6 | potential/vpotential/div/curl/laplacian/hessian |
| `test_partfrac_ext` | `giac_bin_test_partfrac_ext` | 4 | 1 | 扩展 partfrac |
| `test_prog` | `giac_bin_test_prog` | 5 | 3 | for/while/sum/product |
| `test_groebner` | `giac_bin_test_groebner` | 2 | 1 | greduce（Groebner 约化） |
| `test_sturm` | `giac_bin_test_sturm` | 3 | 2 | sturm/sturmab |
| `test_sturm_ext` | `giac_bin_test_sturm_ext` | 3 | 3 | 重根 sturm、sturmab、realroot |
| `test_subst` | `giac_bin_test_subst` | 5 | 4 | subst、函数定义、sum、product |
| `test_trig` | `giac_bin_test_trig` | 4 | 3 | texpand/halftan/lin |

## check/ 回归测试（golden diff）

| 输入文件 | CTest | Golden | 归一化 | 表达式数 | 说明 |
|---------|-------|--------|--------|---------|------|
| `flanex` | `giac_check_flanex` | `flanex.out` | — | 284 | 284 条跨模块综合回归 |
| `testpartfrac_ext` | `giac_check_partfrac_ext` | `partfrac_ext.out` | — | 4 | 扩展 partfrac 黄金 diff |
| `testcas` | `giac_check_cas` | `cas.out` | cas_floats | 271 | 272 条 CAS 综合：复数、gcd、因式分解、积分、矩阵 |
| `testfactor` | `giac_check_factor` | `factor.out` | — | 31 | 31 条因式分解专项 |
| `testgeo` | `giac_check_geo` | `geo.out` | geo | 133 | 133 条几何：三角形、中线、垂直平分线、角平分线 |
| `testintegrate` | `giac_check_integrate` | `integrate.out` | en_US.UTF8 | 67 | 67 条高难度积分+级数+极限（黄金 diff） |
| `testlimit` | `giac_check_limit` | `limit.out` | — | 52 | 52 条极限 |
| `testnormalize` | `giac_check_normalize` | `normalize.out` | — | 188 | 188 条 non_recursive_normal 规范化 |
| `testother` | `giac_check_other` | `other.out` | — | 4 | 4 条高难度积分（含参/无穷限） |
| `testpartfrac` | `giac_check_partfrac` | `partfrac.out` | — | 1 | 部分分式分解 |

## libtommath

| 文件 | CTest | 说明 |
|------|-------|------|
| `libtommath/tests/smoke_test.c` | `tommath_smoke_test` | mp_init/mp_add/mp_cmp 冒烟 |

## 已知问题

| CTest | 说明 |
|-------|------|
| `giac_check_factor_segfault_known` | `factor((x^202+x^101+1)/(x^2+x+1))` 已知 segfault，标记 WILL_FAIL |
