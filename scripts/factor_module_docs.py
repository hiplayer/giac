#!/usr/bin/env python3
"""Replace generic module headers in giac-poly/src/factor/*.rs with role-specific docs."""

from pathlib import Path

FACTOR = Path(__file__).resolve().parents[1] / "crates/giac-poly/src/factor"

HEADERS = {
    "mod.rs": """\
//! Polynomial factorization over ℚ (and ℤ/pℤ for `factor_poly_mod`).
//!
//! **Upstream:** `gausspol.cc` / `ezgcd.cc` — `do_factor_hensel`, `factor`.
//! **缺口:** FAC-G1 `try_sparse_factor`、FAC-G2 参系数塔、FAC-G3 混合次数二元 Hensel。
//!
//! | Tier | 入口 |
//! |------|------|
//! | **Stable (bounded)** | `factor_into`, `factor_poly`, `factor_into_by_rational_roots` |
//! | **Stable** | `factor_poly_mod`, `factor_mod_irreducibles` |
//!
//! 子模块：`multivariate` → `univariate` / `hensel` / `zassenhaus` / `patterns`。
//! 全函数 tier 见 `.doc/giac-poly-api-stability.md` § `factor/mod.rs`。
""",
    "multivariate.rs": """\
//! Multivariate factorization over ℚ: pattern table → univariate → main-variable tower.
//!
//! **Pipeline:** `factor_multivariate` → `factor_multivariate_rec` → `factor_wrt_main_var`.
//! **Stable (bounded):** `factor_into_poly`, `factor_multivariate`.
//! **Pipeline private:** `factor_multivariate_rec`, `factor_wrt_main_var`.
""",
    "univariate.rs": """\
//! Univariate factorization over ℚ: rational roots, quadratics, Zassenhaus, low-degree patterns.
//!
//! **Stable (bounded):** `factor_univariate_flat`, `factor_univariate_pairs`, `factor_power_pairs`.
//! **Partial:** `try_factor_biquadratic`, `try_factor_two_cubics`.
//! **Pipeline private:** `find_rational_root`, `factor_square_free`, `factor_quadratic`, …
""",
    "hensel.rs": """\
//! Bivariate factorization in ℚ[y][x] via Hensel lifting @ y=0 + interpolation fallback.
//!
//! **Upstream:** `gausspol.cc` `try_hensel_lift_factor`, `hensel_lift`.
//! **Partial:** `try_hensel_lift_bivariate` (FAC-G3 混合次数仍可能 None).
//! **Pipeline private:** `hensel_lift_at_zero`, `try_hensel_lift_interp`, rat-vector helpers.
""",
    "zassenhaus.rs": """\
//! Univariate Zassenhaus + modular Hensel lifting (integer poly → factors over ℚ).
//!
//! **Partial:** `try_zassenhaus_factor`.
//! **Pipeline private:** modular egcd, Hensel lift, factor combination recovery.
""",
    "fpx.rs": """\
//! Irreducible factorization of polynomials over finite fields (Cantor–Zassenhaus).
//!
//! **Stable:** `factor_fpx`, `degree`.
//! **Pipeline private:** Yun square-free, distinct-degree, CZ block split, Berlekamp-style linear.
""",
    "modular.rs": """\
//! Display-level factorization mod p: `Poly` → `PolyMod` → `factor_fpx` → lift back.
//!
//! **Stable:** `factor_poly_mod`.
//! **Pipeline private:** `modpoly_to_poly`.
""",
    "poly_uni.rs": """\
//! Polynomials as univariate in main var with coefficients in nested Poly ring (bivariate steps).
//!
//! **Stable / Partial:** `coeff_wrt_poly`, `content_wrt`, `factor_sqff_over_coeff_ring`, …
//! **Pipeline private:** Kronecker / eval lift fallbacks (`try_kronecker_bivariate`, …).
""",
    "power.rs": """\
//! Perfect-power detection: `(base)^exp` and `(linear)^n` shapes before general factor.
//!
//! **Stable:** `as_perfect_power`, `try_linear_power`.
//! **Pipeline private:** `try_nth_root`, `try_binomial_square`.
""",
    "patterns.rs": """\
//! Fast-path pattern factorization (x^n±1, x^n±y^n, cyclotomic hooks) before Hensel.
//!
//! **Partial:** `try_factor_patterns`, `factor_xn_minus_one_display`.
//! **Pipeline private:** `factor_xn_minus_yn`, `is_binomial_diff_power`, …
""",
    "cyclotomic.rs": """\
//! Cyclotomic polynomials Φ_n and specialized x^n ± 1 factorization.
//!
//! **Stable:** `cyclotomic_poly`.
//! **Partial:** `factor_xn_minus_one`, `try_factor_xn_minus_one`, `try_factor_xn_plus_one`, …
""",
    "sqrt.rs": """\
//! Quadratic factorization with sqrt display (Expr string helpers for giac-simplify).
//!
//! **Partial:** `quadratic_sqrt_factor_exprs`.
""",
    "util.rs": """\
//! Shared utilities for factor pipeline: vars, content, primitive part, nth roots, …
//!
//! **Stable:** `vars_in`, `ratio_perfect_sqrt`, `coeff_wrt`, `main_var`, …
""",
    "tracer.rs": """\
//! Regression tests mapped from upstream `check/testfactor` (Issue 2.3/2.4).
//!
//! L20/L22 仍 `#[ignore]`（FAC-G2/G3）。测试辅助 `x()`, `y()`, `assert_factors` 为 Pipeline private。
""",
}


def strip_old_header(text: str) -> str:
    lines = text.splitlines(keepends=True)
    i = 0
    if lines and lines[0].startswith("#!"):
        i = 1
    while i < len(lines) and (lines[i].startswith("//!") or lines[i].strip() == ""):
        i += 1
    return "".join(lines[i:])


def main() -> None:
    for name, header in HEADERS.items():
        path = FACTOR / name
        if not path.exists():
            continue
        body = strip_old_header(path.read_text())
        path.write_text(header + "\n" + body)
        print("updated", name)


if __name__ == "__main__":
    main()
