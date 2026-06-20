# 测试覆盖 API → at_* → Rust crate 对照

## MVP API 是什么

**MVP API** = **M**inimum **V**iable **P**roduct **API**（最小可行产品 API）。

giac-rs Phase 4 要实现的 **builtin 函数子集**，以 **`bin/` + `check/` 黄金回归** 能 eval 为准，**不是** upstream giac-2.0.0 全库 1:1 移植。

| 范围 | 数量 | 说明 |
|------|------|------|
| giac 全库 `at_*` 注册 | **~1852** | 含 GUI、TI/Maple 兼容、极少用函数 |
| **MVP API**（本表） | **~210** | `functional-coverage.md` / 测试规格归纳；Rust 侧见 `FuncKind` / `BUILTINS` |
| 本表行数 | **~250** | 含别名、多 crate 映射；与「~210 种 API」为同一批能力的不同计数口径 |

**判定：** 某函数在 MVP 内 ⟺ 出现在本表且被 conformance / check 覆盖；缺口以「黄金行能否 eval」为准，而非 upstream 是否有对应 C++ 实现。

**不做 MVP：** 全库其余 ~1600+ 个 `at_*`（如 `gbasis` 依赖 CoCoA stub、TI 方言等）— 见 [module-division.md](module-division.md) §6、[rust-migration-plan.md](rust-migration-plan.md) §1.1。

---

来源：`functional-coverage.md` 中 **250** 个 API；giac 全库约 **1852** 个 `at_*` 注册，MVP 仅实现本表。

| API | C++ `at_*` | 主注册文件 | 功能域 | Rust crate |
|-----|------------|------------|--------|------------|
| `Beta` | `at_Beta` | `moyal.cc` | 特殊函数 | `giac-special` |
| `UTPN` | `at_UTPN` | `moyal.cc` | 特殊函数 | `giac-special` |
| `abcuv` | `at_abcuv` | `ifactor.cc` | 多项式与因式分解 | `giac-poly` |
| `abs` | `at_abs` | `usual.cc` | 复数与基础 | `giac-core` |
| `acos2asin` | `at_acos2asin` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `acos2atan` | `at_acos2atan` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `acsc` | `at_acsc` | `misc.cc` | 其他 | `giac-core` |
| `adjoint_matrix` | `at_adjoint_matrix` | `misc.cc` | 线性代数 | `giac-linalg` |
| `append` | `at_append` | `prog.cc` | 其他 | `giac-core` |
| `arg` | `at_arg` | `usual.cc` | 复数与基础 | `giac-core` |
| `asec` | `at_asec` | `misc.cc` | 其他 | `giac-core` |
| `asin2acos` | `at_asin2acos` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `asin2atan` | `at_asin2atan` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `assume` | `at_assume` | `usual.cc` | 其他 | `giac-core` |
| `atan2acos` | `at_atan2acos` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `atan2asin` | `at_atan2asin` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `axq` | `at_axq` | `?` | 线性代数 | `giac-linalg` |
| `bernoulli` | `at_bernoulli` | `intg.cc` | 数论 | `giac-num` |
| `canonical_form` | `at_canonical_form` | `misc.cc` | 多项式与因式分解 | `giac-poly` |
| `cas_setup` | `at_cas_setup` | `prog.cc` | 其他 | `giac-core` |
| `changebase` | `at_changebase` | `misc.cc` | 线性代数 | `giac-linalg` |
| `charpoly` | `at_charpoly` | `misc.cc` | 线性代数 | `giac-linalg` |
| `chinrem` | `at_chinrem` | `usual.cc` | 模运算 | `giac-poly` |
| `coeff` | `at_coeff` | `misc.cc` | 其他 | `giac-core` |
| `color` | `at_color` | `plot.cc` | 其他 | `giac-core` |
| `conj` | `at_conj` | `usual.cc` | 复数与基础 | `giac-core` |
| `content` | `at_content` | `misc.cc` | 多项式与因式分解 | `giac-poly` |
| `cot` | `at_cot` | `misc.cc` | 其他 | `giac-core` |
| `cross` | `at_cross` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `csc` | `at_csc` | `misc.cc` | 其他 | `giac-core` |
| `curl` | `at_curl` | `permu.cc` | 微积分 | `giac-calculus` |
| `cycles2permu` | `at_cycles2permu` | `permu.cc` | 排列与置换 | `giac-num` |
| `cyclotomic` | `at_cyclotomic` | `usual.cc` | 多项式与因式分解 | `giac-poly` |
| `degree` | `at_degree` | `misc.cc` | 其他 | `giac-core` |
| `derive` | `at_derive` | `derive.cc` | 微积分 | `giac-calculus` |
| `desolve` | `at_desolve` | `desolve.cc` | 微分方程 | `giac-ode` |
| `det` | `at_det` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `diff` | `at_diff` | `derive.cc` | 微积分 | `giac-calculus` |
| `divergence` | `at_divergence` | `permu.cc` | 微积分 | `giac-calculus` |
| `divis` | `at_divis` | `ifactor.cc` | 多项式与因式分解 | `giac-poly` |
| `divisors` | `at_divisors` | `misc.cc` | 数论 | `giac-num` |
| `divpc` | `at_divpc` | `misc.cc` | 多项式与因式分解 | `giac-poly` |
| `dot` | `at_dot` | `misc.cc` | 线性代数 | `giac-linalg` |
| `e2r` | `at_e2r` | `sym2poly.cc` | 多项式与因式分解 | `giac-poly` |
| `egcd` | `at_egcd` | `usual.cc` | 多项式与因式分解 | `giac-poly` |
| `egv` | `at_egv` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `egvl` | `at_egvl` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `epsilon2zero` | `at_epsilon2zero` | `misc.cc` | 化简与规范化 | `giac-simplify` |
| `equal` | `at_equal` | `?` | 其他 | `giac-core` |
| `est_cycle` | `at_est_cycle` | `?` | 排列与置换 | `giac-num` |
| `est_permu` | `at_est_permu` | `?` | 排列与置换 | `giac-num` |
| `euler` | `at_euler` | `ifactor.cc` | 数论 | `giac-num` |
| `eval` | `at_eval` | `usual.cc` | 其他 | `giac-core` |
| `evalf` | `at_evalf` | `usual.cc` | 其他 | `giac-core` |
| `exlr` | `at_exlr` | `?` | 其他 | `giac-core` |
| `exp2pow` | `at_exp2pow` | `subst.cc` | 其他 | `giac-core` |
| `expand` | `at_expand` | `lin.cc` | 其他 | `giac-core` |
| `f` | `at_f` | `?` | 其他 | `giac-core` |
| `factor` | `at_factor` | `sym2poly.cc` | 多项式与因式分解 | `giac-poly` |
| `factors` | `at_factors` | `ifactor.cc` | 多项式与因式分解 | `giac-poly` |
| `fcoeffs` | `at_fcoeffs` | `?` | 多项式与因式分解 | `giac-poly` |
| `fdistrib` | `at_fdistrib` | `?` | 其他 | `giac-core` |
| `float2rational` | `at_float2rational` | `misc.cc` | 多项式与因式分解 | `giac-poly` |
| `for` | `at_for` | `prog.cc` | 其他 | `giac-core` |
| `fourier_an` | `at_fourier_an` | `intg.cc` | 拉普拉斯与傅里叶 | `giac-transform` |
| `fourier_cn` | `at_fourier_cn` | `intg.cc` | 拉普拉斯与傅里叶 | `giac-transform` |
| `froots` | `at_froots` | `?` | 方程求解 | `giac-solve` |
| `fsolve` | `at_fsolve` | `solve.cc` | 其他 | `giac-core` |
| `fxnd` | `at_fxnd` | `ifactor.cc` | 微积分 | `giac-calculus` |
| `gamma` | `at_gamma` | `?` | 特殊函数 | `giac-special` |
| `gauss` | `at_gauss` | `gauss.cc` | 线性代数 | `giac-linalg` |
| `gbasis` | `at_gbasis` | `solve.cc` | 多项式与因式分解 | `giac-poly` |
| `gcd` | `at_gcd` | `usual.cc` | 多项式与因式分解 | `giac-poly` |
| `gramschmidt` | `at_gramschmidt` | `misc.cc` | 线性代数 | `giac-linalg` |
| `greduce` | `at_greduce` | `solve.cc` | 多项式与因式分解 | `giac-poly` |
| `hadamard` | `at_hadamard` | `permu.cc` | 线性代数 | `giac-linalg` |
| `halftan` | `at_halftan` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `hermite` | `at_hermite` | `permu.cc` | 特殊函数 | `giac-special` |
| `hessian` | `at_hessian` | `permu.cc` | 微积分 | `giac-calculus` |
| `hilbert` | `at_hilbert` | `permu.cc` | 线性代数 | `giac-linalg` |
| `horner` | `at_horner` | `modpoly.cc` | 多项式与因式分解 | `giac-poly` |
| `ibp` | `at_ibp` | `?` | 微积分 | `giac-calculus` |
| `ibpdv` | `at_ibpdv` | `intg.cc` | 微积分 | `giac-calculus` |
| `ibpu` | `at_ibpu` | `misc.cc` | 微积分 | `giac-calculus` |
| `ichinrem` | `at_ichinrem` | `usual.cc` | 其他 | `giac-core` |
| `idivis` | `at_idivis` | `ifactor.cc` | 多项式与因式分解 | `giac-poly` |
| `idn` | `at_idn` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `ifactors` | `at_ifactors` | `ifactor.cc` | 多项式与因式分解 | `giac-poly` |
| `ilaplace` | `at_ilaplace` | `desolve.cc` | 拉普拉斯与傅里叶 | `giac-transform` |
| `im` | `at_im` | `usual.cc` | 复数与基础 | `giac-core` |
| `image` | `at_image` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `int` | `at_int` | `ti89.cc` | 微积分 | `giac-calculus` |
| `integrate` | `at_integrate` | `intg.cc` | 微积分 | `giac-calculus` |
| `inv` | `at_inv` | `usual.cc` | 其他 | `giac-core` |
| `iquo` | `at_iquo` | `usual.cc` | 多项式与因式分解 | `giac-poly` |
| `irem` | `at_irem` | `usual.cc` | 多项式与因式分解 | `giac-poly` |
| `is_element` | `at_is_element` | `plot.cc` | 其他 | `giac-core` |
| `is_permu` | `at_is_permu` | `permu.cc` | 排列与置换 | `giac-num` |
| `is_prime` | `at_is_prime` | `usual.cc` | 其他 | `giac-core` |
| `isprime` | `at_isprime` | `ti89.cc` | 数论 | `giac-num` |
| `jacobi_symbol` | `at_jacobi_symbol` | `usual.cc` | 其他 | `giac-core` |
| `jordan` | `at_jordan` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `ker` | `at_ker` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `lagrange` | `at_lagrange` | `misc.cc` | 其他 | `giac-core` |
| `laguerre` | `at_laguerre` | `permu.cc` | 特殊函数 | `giac-special` |
| `laplace` | `at_laplace` | `desolve.cc` | 拉普拉斯与傅里叶 | `giac-transform` |
| `laplacian` | `at_laplacian` | `permu.cc` | 微积分 | `giac-calculus` |
| `latex` | `at_latex` | `tex.cc` | 公式导出 | `giac-emit` |
| `lcm` | `at_lcm` | `usual.cc` | 多项式与因式分解 | `giac-poly` |
| `legendre` | `at_legendre` | `permu.cc` | 特殊函数 | `giac-special` |
| `limit` | `at_limit` | `series.cc` | 微积分 | `giac-calculus` |
| `lin` | `at_lin` | `lin.cc` | 三角与双曲 | `giac-simplify` |
| `linsolve` | `at_linsolve` | `solve.cc` | 线性代数 | `giac-linalg` |
| `ln` | `at_ln` | `usual.cc` | 其他 | `giac-core` |
| `lname` | `at_lname` | `prog.cc` | 其他 | `giac-core` |
| `lncollect` | `at_lncollect` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `lu` | `at_lu` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `lvar` | `at_lvar` | `sym2poly.cc` | 其他 | `giac-core` |
| `mathml` | `at_mathml` | `mathml.cc` | 公式导出 | `giac-emit` |
| `mkisom` | `at_mkisom` | `isom.cc` | 几何同构 | `giac-geo` |
| `modp` | `at_modp` | `maple.cc` | 其他 | `giac-core` |
| `moyal` | `at_moyal` | `moyal.cc` | 特殊函数 | `giac-special` |
| `nCr` | `at_nCr` | `ti89.cc` | 数论 | `giac-num` |
| `nPr` | `at_nPr` | `ti89.cc` | 数论 | `giac-num` |
| `newton` | `at_newton` | `solve.cc` | 其他 | `giac-core` |
| `nextprime` | `at_nextprime` | `usual.cc` | 数论 | `giac-num` |
| `nodisp` | `at_nodisp` | `prog.cc` | 其他 | `giac-core` |
| `non_recursive_normal` | `at_non_recursive_normal` | `sym2poly.cc` | 化简与规范化 | `giac-simplify` |
| `normal` | `at_normal` | `sym2poly.cc` | 化简与规范化 | `giac-simplify` |
| `p1op2` | `at_p1op2` | `permu.cc` | 排列与置换 | `giac-num` |
| `partfrac` | `at_partfrac` | `sym2poly.cc` | 多项式与因式分解 | `giac-poly` |
| `pcar` | `at_pcar` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `pcoeff` | `at_pcoeff` | `vecteur.cc` | 其他 | `giac-core` |
| `permu2cycles` | `at_permu2cycles` | `permu.cc` | 排列与置换 | `giac-num` |
| `permu2mat` | `at_permu2mat` | `permu.cc` | 排列与置换 | `giac-num` |
| `permuorder` | `at_permuorder` | `permu.cc` | 排列与置换 | `giac-num` |
| `peval` | `at_peval` | `vecteur.cc` | 多项式与因式分解 | `giac-poly` |
| `plotfunc` | `at_plotfunc` | `plot.cc` | 几何与绘图 | `giac-geo` |
| `pmin` | `at_pmin` | `misc.cc` | 线性代数 | `giac-linalg` |
| `polygamma` | `at_polygamma` | `moyal.cc` | 特殊函数 | `giac-special` |
| `potential` | `at_potential` | `misc.cc` | 微积分 | `giac-calculus` |
| `preval` | `at_preval` | `misc.cc` | 微积分 | `giac-calculus` |
| `prevprime` | `at_prevprime` | `usual.cc` | 数论 | `giac-num` |
| `product` | `at_product` | `ti89.cc` | 替换与程序 | `giac-prog` |
| `proot` | `at_proot` | `vecteur.cc` | 微积分 | `giac-calculus` |
| `propfrac` | `at_propfrac` | `ifactor.cc` | 多项式与因式分解 | `giac-poly` |
| `psi` | `at_psi` | `?` | 特殊函数 | `giac-special` |
| `ptayl` | `at_ptayl` | `misc.cc` | 多项式与因式分解 | `giac-poly` |
| `purge` | `at_purge` | `rpn.cc` | 其他 | `giac-core` |
| `qr` | `at_qr` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `quaternion` | `at_quaternion` | `quater.cc` | 四元数 | `giac-special` |
| `quo` | `at_quo` | `usual.cc` | 多项式与因式分解 | `giac-poly` |
| `quote` | `at_quote` | `usual.cc` | 其他 | `giac-core` |
| `qxa` | `at_qxa` | `?` | 线性代数 | `giac-linalg` |
| `r2e` | `at_r2e` | `sym2poly.cc` | 多项式与因式分解 | `giac-poly` |
| `randNorm` | `at_randNorm` | `moyal.cc` | 特殊函数 | `giac-special` |
| `randchisquare` | `at_randchisquare` | `moyal.cc` | 特殊函数 | `giac-special` |
| `rank` | `at_rank` | `misc.cc` | 线性代数 | `giac-linalg` |
| `ranm` | `at_ranm` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `ratnormal` | `at_ratnormal` | `maple.cc` | 化简与规范化 | `giac-simplify` |
| `re` | `at_re` | `usual.cc` | 复数与基础 | `giac-core` |
| `realroot` | `at_realroot` | `csturm.cc` | 方程求解 | `giac-solve` |
| `rem` | `at_rem` | `usual.cc` | 多项式与因式分解 | `giac-poly` |
| `reorder` | `at_reorder` | `misc.cc` | 化简与规范化 | `giac-simplify` |
| `resultant` | `at_resultant` | `sym2poly.cc` | 多项式与因式分解 | `giac-poly` |
| `risch` | `at_risch` | `risch.cc` | 微积分 | `giac-calculus` |
| `rootof` | `at_rootof` | `alg_ext.cc` | 微积分 | `giac-calculus` |
| `roots` | `at_roots` | `misc.cc` | 多项式与因式分解 | `giac-poly` |
| `rref` | `at_rref` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `sec` | `at_sec` | `misc.cc` | 其他 | `giac-core` |
| `series` | `at_series` | `series.cc` | 微积分 | `giac-calculus` |
| `sign` | `at_sign` | `usual.cc` | 复数与基础 | `giac-core` |
| `signature` | `at_signature` | `permu.cc` | 线性代数 | `giac-linalg` |
| `simp2` | `at_simp2` | `ifactor.cc` | 多项式与因式分解 | `giac-poly` |
| `simplify` | `at_simplify` | `subst.cc` | 微积分 | `giac-calculus` |
| `size` | `at_size` | `vecteur.cc` | 其他 | `giac-core` |
| `smod` | `at_smod` | `usual.cc` | 其他 | `giac-core` |
| `solve` | `at_solve` | `solve.cc` | 方程求解 | `giac-solve` |
| `sort` | `at_sort` | `prog.cc` | 其他 | `giac-core` |
| `sto` | `at_sto` | `usual.cc` | 其他 | `giac-core` |
| `sturm` | `at_sturm` | `alg_ext.cc` | Sturm 根隔离 | `giac-solve` |
| `sturmab` | `at_sturmab` | `alg_ext.cc` | Sturm 根隔离 | `giac-solve` |
| `subst` | `at_subst` | `usual.cc` | 替换与程序 | `giac-prog` |
| `sum` | `at_sum` | `usual.cc` | 替换与程序 | `giac-prog` |
| `suppress` | `at_suppress` | `misc.cc` | 其他 | `giac-core` |
| `svd` | `at_svd` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `syst2mat` | `at_syst2mat` | `permu.cc` | 线性代数 | `giac-linalg` |
| `tan` | `at_tan` | `usual.cc` | 其他 | `giac-core` |
| `tan2sincos` | `at_tan2sincos` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `taylor` | `at_taylor` | `ti89.cc` | 微积分 | `giac-calculus` |
| `tchebyshev1` | `at_tchebyshev1` | `permu.cc` | 特殊函数 | `giac-special` |
| `tchebyshev2` | `at_tchebyshev2` | `permu.cc` | 特殊函数 | `giac-special` |
| `tcollect` | `at_tcollect` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `texpand` | `at_texpand` | `lin.cc` | 三角与双曲 | `giac-simplify` |
| `tlin` | `at_tlin` | `lin.cc` | 三角与双曲 | `giac-simplify` |
| `trace` | `at_trace` | `vecteur.cc` | 其他 | `giac-core` |
| `tran` | `at_tran` | `vecteur.cc` | 线性代数 | `giac-linalg` |
| `triangle` | `at_triangle` | `plot.cc` | 几何与绘图 | `giac-geo` |
| `trig2exp` | `at_trig2exp` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `trigcos` | `at_trigcos` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `trigsin` | `at_trigsin` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `trn` | `at_trn` | `permu.cc` | 线性代数 | `giac-linalg` |
| `truncate` | `at_truncate` | `misc.cc` | 化简与规范化 | `giac-simplify` |
| `tsimplify` | `at_tsimplify` | `subst.cc` | 三角与双曲 | `giac-simplify` |
| `valuation` | `at_valuation` | `misc.cc` | 多项式与因式分解 | `giac-poly` |
| `vandermonde` | `at_vandermonde` | `permu.cc` | 线性代数 | `giac-linalg` |
| `vpotential` | `at_vpotential` | `misc.cc` | 微积分 | `giac-calculus` |
| `while` | `at_while` | `?` | 其他 | `giac-core` |
| `xyztrange` | `at_xyztrange` | `plot.cc` | 几何与绘图 | `giac-geo` |
| `zeta` | `at_zeta` | `?` | 特殊函数 | `giac-special` |
