#!/usr/bin/env python3
"""Add API tier comments to algorithm crate sources.

Idempotent: skips functions that already have a tier marker in the preceding lines.
Run from giac-rs/:
  python3 scripts/annotate_api_tiers.py              # annotate missing tiers
  python3 scripts/annotate_api_tiers.py --inventory    # refresh *-api-stability.md Per-file tables
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TIER_MARKERS = (
    "**Stable**",
    "**Stable (bounded)**",
    "**Stable (crate-internal)**",
    "**Partial**",
    "**Pipeline**",
    "**Pipeline private**",
    "**Temporary**",
)

# (src_root relative to giac-rs/crates, stability doc basename, skip filenames)
CRATE_TARGETS: list[tuple[str, str, frozenset[str]]] = [
    ("giac-simplify/src", "giac-simplify-api-stability.md", frozenset()),
    ("giac-poly/src", "giac-poly-api-stability.md", frozenset({"tests_phase2.rs"})),
    ("giac-calculus/src", "giac-calculus-api-stability.md", frozenset()),
    ("giac-core/src/algebra", "giac-core-algebra-api-stability.md", frozenset({"test_fixtures.rs"})),
    ("giac-solve/src", "giac-solve-api-stability.md", frozenset()),
    ("giac-ode/src", "giac-ode-api-stability.md", frozenset()),
    ("giac-groebner/src", "giac-groebner-api-stability.md", frozenset()),
]

# name -> (tier, short description)
FN_TIERS: dict[str, tuple[str, str]] = {
    # --- giac-simplify ---
    "expand": ("Stable", "expand with Full policy"),
    "expand_with_policy": ("Stable", "expand with explicit ExpandPolicy"),
    "expand_polynomial": ("Stable", "expand without distributing through exp/ln"),
    "expr_contains_exp_ln": ("Pipeline private", "predicate: subtree contains exp or ln"),
    "expand_mul_pair": ("Pipeline private", "distribute one mul factor over add"),
    "expand_pow": ("Pipeline private", "expand integer powers and mod-poly powers"),
    "repeated_mul": ("Pipeline private", "repeated multiply for small integer power"),
    "expand_binomial": ("Pipeline private", "binomial power via poly or Expr::pow"),
    "normal": ("Stable", "expand then polynomial collect"),
    "modulus_from_expr": ("Pipeline private", "coerce Expr modulus to i64"),
    "mod_err": ("Pipeline private", "removed — poly returns EvalError directly"),
    "ratnormal": ("Stable", "single fraction in lowest terms"),
    "ratnormal_algext": ("Temporary", "AlgExt shim via eval; retire A-03"),
    "rational_parts": ("Pipeline private", "Expr to (num Poly, den Poly)"),
    "rational_add": ("Pipeline private", "add rationals with common denominator"),
    "reduce_fraction": ("Pipeline private", "gcd-reduce num/den Poly pair"),
    "factor": ("Stable (bounded)", "structural then giac-poly factor"),
    "factor_expr": ("Pipeline private", "recursive factor on Mul/Pow/Frac"),
    "flatten_mul": ("Pipeline private", "flatten Mul to factor vec"),
    "factor_poly_form": ("Pipeline private", "normal→poly→factor_into chain"),
    "try_factor_quadratic_rootof": ("Temporary", "quadratic→rootof linear factors"),
    "try_factor_quadratic_sqrt": ("Temporary", "quadratic sqrt factors when with_sqrt"),
    "rational_num_den": ("Pipeline private", "Expr leaf to (num,den) Poly"),
    "ifactor": ("Stable", "integer prime factorization display"),
    "is_probable_prime_u64": ("Pipeline private", "Miller-Rabin for u64"),
    "mod_mul": ("Pipeline private", "u64 modular multiply"),
    "mod_pow": ("Pipeline private", "u64 modular exponentiation"),
    "pollard_rho": ("Pipeline private", "Pollard rho on BigInt"),
    "pollard_rho_u64": ("Pipeline private", "Pollard rho on u64"),
    "sub": ("Stable", "build a-b Expr"),
    "is_zero": ("Stable", "normal(e)==0"),
    "assert_equiv": ("Stable", "normal(a-b)==0 with narrow radical drift"),
    "canonical_radical": ("Temporary", "drift: 1/sqrt(n)↔sqrt(n)/n for assert_equiv"),
    "inv_sqrt_to_mul": ("Temporary", "drift helper for canonical_radical"),
    "texpand": ("Partial", "trig/exp/ln arg expand then algebraic expand"),
    "texpand_rec": ("Pipeline private", "recursive texpand on Expr tree"),
    "expand_sin_arg": ("Pipeline private", "sin angle-sum and n*x rules"),
    "expand_cos_arg": ("Pipeline private", "cos angle-sum and n*x rules"),
    "expand_exp_arg": ("Pipeline private", "exp of sum → product of exp"),
    "expand_ln_arg": ("Pipeline private", "ln of product → sum of ln"),
    "expand_sin_nx": ("Pipeline private", "sin(nx) for small integer n"),
    "expand_cos_nx": ("Pipeline private", "cos(nx) for small integer n"),
    "halftan": ("Partial", "half-angle tan on narrow pattern"),
    "halftan_half_angle_rational": ("Pipeline private", "Weierstrass tan(v/2) form"),
    "tan_half": ("Pipeline private", "tan(v/2) Expr builder"),
    "lin": ("Partial", "linearize exp products and (exp+1)^2"),
    "lin_rec": ("Pipeline private", "recursive lin on Expr tree"),
    "contains_exp": ("Pipeline private", "subtree contains exp"),
    "expand_integer_pow": ("Pipeline private", "expand exp-base integer power"),
    "sin_expr": ("Pipeline private", "build sin Expr"),
    "cos_expr": ("Pipeline private", "build cos Expr"),
    "integer_multiple": ("Pipeline private", "detect n*x integer multiple"),
    "detect_halftan_tan": ("Pipeline private", "detect sin(2x)/(1+cos(2x))"),
    "as_frac": ("Pipeline private", "view Expr as Frac pair"),
    "unwrap_unit_mul_owned": ("Pipeline private", "peel unit coefficient from mul"),
    "sin_double_angle": ("Pipeline private", "detect sin(2k*x) in halftan"),
    "cos_double_angle": ("Pipeline private", "detect cos(2k*x) in halftan"),
    "lin_exp_plus_one_pow": ("Pipeline private", "expand (exp+1)^2 only"),
    "install_simplify": ("Stable", "register DefaultAlgebraPlugin"),
    "xcas_default": ("Stable", "Context with simplify plugin"),
    # --- giac-poly core ---
    "zero": ("Stable", "Poly zero"),
    "one": ("Stable", "Poly one"),
    "constant": ("Stable", "Poly scalar constant"),
    "var": ("Stable", "Poly univariate generator"),
    "is_zero": ("Stable", "Poly is zero"),
    "is_one": ("Stable", "Poly is one"),
    "leading_term": ("Stable", "leading term by total degree"),
    "leading_term_lex": ("Stable", "leading term with variable order"),
    "term": ("Stable", "monomial × coefficient"),
    "degree": ("Stable", "total degree"),
    "add": ("Stable", "Poly addition"),
    "sub": ("Stable", "Poly subtraction"),
    "neg": ("Stable", "Poly negation"),
    "mul": ("Stable", "Poly multiplication"),
    "mul_scalar": ("Stable", "scale Poly by rational"),
    "pow": ("Stable", "Poly integer power"),
    "content": ("Stable", "integer content of Poly"),
    "primitive_part": ("Stable", "divide out content"),
    "monic": ("Stable", "divide by leading coeff"),
    "div_rem": ("Stable", "multivariate division with remainder"),
    "div_exact": ("Stable", "exact division if remainder zero"),
    "gcd": ("Stable", "Poly gcd via subresultant"),
    "lcm": ("Stable", "Poly lcm"),
    "horner": ("Stable", "Horner eval at rational point"),
    "integer_content_gcd": ("Pipeline private", "gcd of rational coeffs as Ratio"),
    "quo": ("Stable", "exact quotient Poly/ Poly"),
    "rem": ("Stable", "remainder Poly/ Poly"),
    "egcd": ("Stable", "extended gcd (s,t,g)"),
    "simp2": ("Stable", "reduce fraction pair by gcd"),
    "abcuv": ("Stable", "Bezout coeffs for au+bv=c"),
    "gauss": ("Stable", "Gauss elimination on Poly rows"),
    "smod": ("Stable", "symmetric mod for i64"),
    "irem": ("Stable", "integer remainder"),
    "modp": ("Stable", "Poly → PolyMod mod p"),
    "resultant": ("Stable", "univariate resultant"),
    "coeff_at": ("Stable", "univariate coefficient at exponent"),
    "univariate_degree": ("Stable", "degree w.r.t. var"),
    "roots": ("Stable (bounded)", "low-degree exact roots as Poly factors"),
    "chinrem": ("Stable", "Chinese remainder two residues"),
    "chinrem_lists": ("Stable", "CRT fold over lists"),
    "partfrac_terms": ("Stable (bounded)", "partial fraction terms"),
    "partfrac_rational_terms": ("Stable (bounded)", "partfrac with poly part"),
    "bigint_pow": ("Pipeline private", "BigInt pow with overflow check"),
    # univariate
    "univariate_derivative": ("Stable", "derivative w.r.t. var"),
    "square_free_factorization": ("Stable", "Yun square-free factors"),
    "square_free_part": ("Stable", "product of square-free factors"),
    "substitute_univariate": ("Stable", "substitute var → Poly"),
    "odd_multiplicity_part": ("Stable", "odd multiplicity factor"),
    "gcd_univariate": ("Stable", "univariate gcd"),
    "sturm_sequence": ("Stable", "Sturm chain"),
    "eval_univariate_at": ("Stable", "Horner eval"),
    "sign_variations": ("Stable", "sign change count in sequence"),
    "sturm_sign_variations_at": ("Stable", "Sturm sign count at point"),
    "sturmab_count": ("Stable", "root count in (a,b)"),
    "subresultant_gcd": ("Stable", "multivariate gcd subresultant"),
    # factor entry
    "factor_into": ("Stable (bounded)", "irreducible factors over Q or None"),
    "factor_poly": ("Stable (bounded)", "legacy factor display product"),
    "factor_into_by_rational_roots": ("Stable (bounded)", "univariate via rational roots"),
    "factor_poly_mod": ("Stable", "factor mod p display"),
    "factor_mod_irreducibles": ("Stable", "irreducibles over F_p"),
    "as_perfect_power": ("Stable", "detect base^exp decomposition"),
    "try_linear_power": ("Partial", "detect (linear)^n"),
    "quadratic_sqrt_factor_exprs": ("Partial", "quadratic sqrt factor strings"),
    "factor_power_pairs": ("Stable (bounded)", "factors with multiplicities"),
    "factor_univariate_flat": ("Partial", "flat irreducible list"),
    "factor_univariate_pairs": ("Partial", "pairs with multiplicity"),
    "try_zassenhaus_factor": ("Partial", "Zassenhaus+Hensel lift"),
    "try_hensel_lift_bivariate": ("Partial", "bivariate Hensel at y=0 + interp fallback; FAC-G3 gap"),
    "lift_factor_from_aux_evals": ("Pipeline private", "Hensel lift factor from auxiliary eval tracks"),
    "try_lift_factors_in_aux_var": ("Pipeline private", "lift bivariate factors via aux variable"),
    "factor_multivariate_rec": ("Pipeline private", "multivariate factor recursion (patterns→uni→main var)"),
    "find_rational_root": ("Pipeline private", "rational root via rational root theorem"),
    "modpoly_to_poly": ("Pipeline private", "PolyMod → Poly over ℤ/pℤ for display"),
    "try_factor_patterns": ("Partial", "cyclotomic/binomial pattern table"),
    "factor_xn_minus_one_display": ("Partial", "x^n-1 display factorization"),
    "factor_into_poly": ("Stable (bounded)", "factor_multivariate ok→Some"),
    "factor_multivariate": ("Stable (bounded)", "multivariate factorization"),
    "factor_fpx": ("Stable", "factor irreducible PolyMod"),
    "vars_in": ("Stable", "sorted variables in Poly"),
    "ratio_perfect_sqrt": ("Stable", "detect perfect square Ratio"),
    "coeff_wrt": ("Stable", "coefficient Poly w.r.t. var^exp"),
    "coeff_wrt_poly": ("Stable", "Poly coeff as Poly in nested var"),
    "content_wrt": ("Stable", "content w.r.t. main var"),
    "primitive_part_wrt": ("Stable", "primitive part w.r.t. var"),
    "term_with_var": ("Stable", "coeff * var^exp as Poly"),
    "derivative_wrt": ("Stable", "derivative w.r.t. var into coeff ring"),
    "square_free_wrt": ("Stable", "square-free factors w.r.t. var"),
    "substitute_poly": ("Stable", "substitute var → Poly"),
    "factor_sqff_over_coeff_ring": ("Partial", "sqff factor over coeff ring"),
    # tresultant
    "num_minus_t_derivative": ("Stable", "RT numerator derivative"),
    "tresultant_eliminate_x": ("Stable", "eliminate x via t-resultant"),
    "eval_param_poly": ("Stable", "substitute parameter in Poly"),
    "rational_roots_in_t": ("Stable", "rational roots in parameter t"),
    "biquadratic_res_conjugate_pairs": ("Partial", "RT biquadratic resolvent"),
    "biquartic_conjugate_pairs": ("Partial", "RT biquartic resolvent"),
    # cyclotomic
    "divisors_u64": ("Pipeline private", "divisors of n"),
    "cyclotomic_poly": ("Stable", "n-th cyclotomic polynomial"),
    "factor_xn_minus_one": ("Partial", "factor x^n-1 via cyclotomic"),
    "factor_x2n_plus_xn_plus_1": ("Partial", "factor x^2n+x^n+1 pattern"),
    "try_factor_x2n_plus_xn_plus_1": ("Partial", "detect and factor x^2n+x^n+1"),
    "try_factor_xn_minus_one": ("Partial", "detect x^n-1"),
    "is_xn_minus_one_poly": ("Pipeline private", "shape test x^n-1"),
    "try_factor_xn_plus_one": ("Partial", "detect x^n+1"),
    # factor/unitary.rs — FAC-G1 unitaryfactor / pzadic / P2a
    "try_unitary_factor": ("Partial", "FAC-G1 last-resort; sparse/Hensel fallback; bounded GCDHEU eval stream"),
    "unitary_factor_rev": ("Partial", "core unitaryfactor loop on vars_rev; pzadic peel + P2a fallback"),
    "reverse_var_order": ("Pipeline private", "upstream tensor reverse on variable indices"),
    "trunc1_drop_var": ("Pipeline private", "drop eval_var tail (upstream trunc1)"),
    "untrunc1_insert_var": ("Pipeline private", "reinsert eval_var with zero exp (upstream untrunc1)"),
    "unitarize": ("Pipeline private", "scale to unitary leading coeff w.r.t. eval_var"),
    "ununitarize": ("Pipeline private", "undo unitarize scaling factor"),
    "initial": ("Pipeline private", "GCDHEU eval base 2·‖p‖∞+2"),
    "base": ("Pipeline private", "unitary eval base accessor"),
    "set_base": ("Pipeline private", "set unitary eval base"),
    "bump_sqff": ("Pipeline private", "sqff micro-bump base += 1"),
    "advance": ("Pipeline private", "upstream eval step ⌊base·73794/27011⌋+1"),
    "as_ratio": ("Pipeline private", "eval base as Ratio<BigInt>"),
    "current": ("Pipeline private", "EvalBaseStream current point"),
    "current_mut": ("Pipeline private", "EvalBaseStream mutable current point"),
    "next": ("Pipeline private", "advance outer eval base or stop when bits > 256"),
    "upstream_bases": ("Pipeline private", "first N bases on EvalBaseStream (tests)"),
    "eval_coeff_groups": ("Pipeline private", "group terms by eval_var exponent"),
    "pow_poly": ("Pipeline private", "integer exponentiation in Poly ring"),
    "draft_from": ("Pipeline private", "build PzadicDraft from eval factor"),
    "lift_candidates": ("Pipeline private", "pzadic lift candidates from draft"),
    "pzadic": ("Pipeline private", "faithful base-B digit lift (dim+1 via eval_var)"),
    "factor_sort_key": ("Pipeline private", "sort key (deg, lc) for eval factors"),
    "sort_eval_factors": ("Pipeline private", "stable sort eval factor slots"),
    "lagrange_interp_coeff": ("Pipeline private", "P2a Lagrange coeff in eval_var"),
    "monic_wrt_main": ("Pipeline private", "normalize factor monic w.r.t. main"),
    "lift_factor_multi_eval": ("Pipeline private", "P2a local-window multi-point coeff lift"),
    "try_lift_and_peel": ("Pipeline private", "pzadic peel then P2a fallback + divides check"),
    "sym_mod_digit": ("Pipeline private", "symmetric mod digit for pzadic expansion"),
    "factor_constant_tail_into": ("Pipeline private", "recurse constant tail into factor list"),
    "factor_constant_tail": ("Pipeline private", "factor tail when main degree → 0"),
    "try_peel_all_at_eval": ("Pipeline private", "batch peel all slots at one eval base"),
    "factor_at_eval": ("Pipeline private", "factor eval image w.r.t. main"),
    "verified_product": ("Pipeline private", "check factor product equals orig"),
    "is_sqff_wrt_main": ("Pipeline private", "sqff test w.r.t. main var"),
    "linfnorm": ("Pipeline private", "L∞ norm of Poly coefficients"),
    # --- giac-core / algebra ---
    "expr_to_poly": ("Stable", "path A: Expr → Poly over Q"),
    "poly_alg_from_expr": ("Stable", "path B: Expr → PolyAlgExt over K"),
    "poly_to_expr": ("Stable", "Poly over Q → Expr"),
    "algext_poly_to_expr": ("Stable", "PolyAlgExt → Expr"),
    "univariate_poly_to_poly1_expr": ("Stable", "univariate Poly → poly1 Expr (high-degree-first)"),
    "poly_mod_to_expr": ("Stable", "PolyMod → Expr"),
    "ratio_to_expr": ("Stable", "Ratio<BigInt> → Expr"),
    "vars_from_expr": ("Stable", "sorted variables in Expr"),
    "expr_contains_alg_coeff": ("Stable", "predicate: choose expr_to_poly vs poly_alg_from_expr"),
    "poly_algext_roots": ("Stable (bounded)", "exact AlgExtC roots deg 1–4; quartic resolvent gap"),
    "contains_algext": ("Stable", "subtree contains AlgExt or rootof"),
    "try_as_algext_data": ("Stable", "view Expr as AlgExtData if present"),
    "try_rootof_to_algext": ("Stable", "Func(RootOf) → AlgExt Expr"),
    "fold_algext_sum": ("Stable", "canonical sum of AlgExt terms"),
    "fold_algext_sum_mode": ("Stable", "fold_algext_sum with mode"),
    "fold_algext_product": ("Stable", "canonical product of AlgExt terms"),
    "fold_complex_algext_sum": ("Stable", "sum with complex AlgExtC parts"),
    "fold_complex_algext_product": ("Stable", "product with complex AlgExtC parts"),
    "algext_square_roots": ("Stable (bounded)", "square roots in extension field"),
    "algext_cube_root": ("Stable (bounded)", "cube root in extension field"),
    "algext_sqrt_branches": ("Stable (bounded)", "sqrt branches as Expr list"),
    "common_ext": ("Stable", "common extension for two AlgExt values"),
    "canonicalize_to_algext_c": ("Stable", "Expr → canonical AlgExtCData"),
    "infer_field": ("Pipeline private", "infer ambient K from PolyAlgExt coefficients"),
    "normalize_coeffs": ("Pipeline private", "lift all coeffs to ambient K"),
    "align_coeff": ("Pipeline private", "align two AlgExtCPolyCoeff to common field; retire FieldSession"),
    "lift_to_field": ("Pipeline private", "embed coeff into target ExtensionField"),
    # --- giac-solve ---
    "eval_solve": ("Stable (bounded)", "solve via poly roots, rootof, or linsolve"),
    "eval_froot": ("Stable (bounded)", "rational roots of univariate poly"),
    "eval_realroot": ("Stable (bounded)", "real roots via Sturm isolation"),
    "eval_fsolve": ("Partial", "Newton numeric solve stub"),
    "eval_sturm": ("Stable", "Sturm sequence for univariate poly"),
    "eval_sturmab": ("Stable", "root count in (a,b) via Sturm"),
    "quadratic_rootof_roots": ("Stable (bounded)", "two rootof branches for quadratic"),
    "biquadratic_rootof_roots": ("Partial", "biquadratic rootof; general quartic NotImplemented"),
    "install_solve": ("Stable", "register DefaultSolvePlugin"),
    # --- giac-ode ---
    "eval_desolve": ("Stable (bounded)", "linear constant-coefficient ODE subset"),
    "install_ode": ("Stable", "register DefaultOdePlugin"),
    # --- giac-groebner ---
    "greduce": ("Stable (bounded)", "multivariate reduce mod ideal (lex)"),
    "greduce_mod": ("Stable (bounded)", "reduce PolyMod mod basis"),
    # --- giac-calculus (stable API not yet inline-marked) ---
    "diff": ("Stable", "symbolic derivative"),
    "eval_diff": ("Stable", "builtin diff evaluator"),
    "integrate": ("Stable", "symbolic integration"),
    "eval_integrate": ("Stable", "builtin integrate evaluator"),
    "eval_limit": ("Stable", "limit evaluator"),
    "eval_series": ("Stable", "series expansion evaluator"),
    "eval_risch": ("Partial", "Risch integration subset"),
    "hermite_reduce": ("Stable", "Hermite reduction on rational tower"),
    "pow2expln": ("Stable", "pow → exp/ln tower rewrite"),
    "risch_tower": ("Stable", "build Risch integration tower"),
    "rlvarx": ("Stable", "Risch log extension variable"),
    "rothstein_trager_integrate": ("Partial", "Rothstein–Trager narrow path"),
    "depends_on_var": ("Stable", "Expr depends on variable"),
    "is_const_wrt": ("Stable", "Expr constant w.r.t. variable"),
    "install_calculus": ("Stable", "register DefaultCalculusPlugin"),
    "canonical_mrv_coeff": ("Stable", "MRV coefficient canonical form"),
    "decompose_mrv_coeff": ("Stable", "decompose MRV coeff into parts"),
    "canonical_exp_diff": ("Stable", "canonical exp-times-(exp-1) form"),
    "match_exp_times_exp_minus_one": ("Stable", "match exp*(exp-1) pattern"),
    "is_exp_minus_one_factor": ("Stable", "predicate: factor is exp(ε)-1"),
}

FN_RE = re.compile(
    r"^(\s*)((?:pub\s*(?:\(\s*crate\s*\))?\s+)?(?:async\s+)?fn\s+)(\w+)"
)


def has_tier(lines: list[str], idx: int) -> bool:
    for j in range(idx - 1, max(idx - 12, -1), -1):
        s = lines[j].strip()
        if not s:
            continue
        if any(m in s for m in TIER_MARKERS):
            return True
        if s.startswith("///") or s.startswith("// **退役"):
            continue
        if s.startswith("#"):
            continue
        break
    return False


def find_test_regions(text: str) -> list[tuple[int, int]]:
    """Return line ranges (start,end) inside `mod tests` blocks."""
    regions: list[tuple[int, int]] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if re.match(r"\s*#\[cfg\(test\)\]", lines[i]) and i + 1 < len(lines):
            if "mod tests" in lines[i + 1]:
                start = i
                brace = 0
                started = False
                for j in range(i + 1, len(lines)):
                    for ch in lines[j]:
                        if ch == "{":
                            brace += 1
                            started = True
                        elif ch == "}":
                            brace -= 1
                    if started and brace == 0:
                        regions.append((start, j + 1))
                        i = j + 1
                        break
                else:
                    i += 1
                continue
        i += 1
    return regions


def in_region(line_no: int, regions: list[tuple[int, int]]) -> bool:
    return any(s <= line_no < e for s, e in regions)


def tier_for(name: str, line: str, is_impl_method: bool) -> tuple[str, str]:
    if name in FN_TIERS:
        return FN_TIERS[name]
    if is_impl_method and name not in ("new", "default"):
        return ("Stable", f"`Poly::{name}`")
    if "pub fn" in line or "pub(crate) fn" in line:
        if name.startswith("try_"):
            return ("Partial", f"optional algorithm path `{name}`")
        return ("Stable", f"`{name}`")
    if name.startswith("try_"):
        return ("Pipeline private", f"optional fallback `{name}`")
    return ("Pipeline private", f"`{name}`")


def comment_line(indent: str, tier: str, desc: str, is_pub: bool) -> str:
    body = f"**{tier}** — {desc}"
    if is_pub and "pub(crate)" not in indent:
        return f"{indent}/// {body}\n"
    return f"{indent}// {body}\n"


def process_file(path: Path) -> int:
    text = path.read_text()
    lines = text.splitlines(keepends=True)
    test_regions = find_test_regions(text)
    added = 0
    out: list[str] = []
    for i, line in enumerate(lines):
        m = FN_RE.match(line.rstrip("\n"))
        if m and not in_region(i, test_regions):
            if not has_tier(lines, i):
                indent, _, name = m.group(1), m.group(2), m.group(3)
                is_impl = i > 0 and "impl " in "".join(
                    l.strip() for l in lines[max(0, i - 15) : i] if l.strip().startswith("impl ")
                )
                tier, desc = tier_for(name, line, is_impl)
                is_pub = re.search(r"\bpub\s+fn\b", line) and "pub(crate)" not in line
                out.append(comment_line(indent, tier, desc, is_pub))
                added += 1
        out.append(line)
    if added:
        path.write_text("".join(out))
    return added


def add_module_inventory_header(path: Path, doc_name: str) -> None:
    text = path.read_text()
    if "**API inventory:**" in text:
        return
    header = (
        f"//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;\n"
        f"//! full module index in `.doc/{doc_name}`.\n//!\n"
    )
    if text.startswith("#!"):
        # insert after first line block of #!
        end = 0
        for i, line in enumerate(text.splitlines(keepends=True)):
            if not line.startswith("//!") and not line.startswith("#!["):
                end = sum(len(l) + 1 for l in text.splitlines(keepends=True)[:i] if True)
                break
        else:
            end = len(text)
        # simpler: prepend if no module doc
        pass
    lines = text.splitlines(keepends=True)
    insert_at = 0
    if lines and lines[0].startswith("#!["):
        insert_at = 1
    elif lines and lines[0].startswith("//!"):
        while insert_at < len(lines) and (
            lines[insert_at].startswith("//!") or lines[insert_at].strip() == ""
        ):
            insert_at += 1
        header = "//!\n" + header
    else:
        pass  # will prepend
    if insert_at == 0 and not text.startswith("//!"):
        path.write_text(header + text)
    elif "**API inventory:**" not in text:
        new_lines = lines[:insert_at] + [header] + lines[insert_at:]
        path.write_text("".join(new_lines))


def write_inventory_docs() -> None:
    doc_root = ROOT.parent / ".doc"
    tier_re = re.compile(
        r"\*\*(Stable \(crate-internal\)|Stable \(bounded\)|Stable|Partial|Temporary|Pipeline private|Pipeline)\*\* — (.+)"
    )
    for src_rel, doc_name, skip in CRATE_TARGETS:
        src_root = ROOT / "crates" / src_rel
        if not src_root.is_dir():
            print(f"skip missing {src_rel}", file=sys.stderr)
            continue
        by_file: dict[str, list[tuple[str, str, str]]] = {}
        for path in sorted(src_root.rglob("*.rs")):
            if path.name in skip:
                continue
            rel = str(path.relative_to(src_root))
            items: list[tuple[str, str, str]] = []
            lines = path.read_text().splitlines()
            for i, line in enumerate(lines):
                m = FN_RE.match(line)
                if not m:
                    continue
                name = m.group(3)
                tier, desc = "Pipeline private", f"`{name}`"
                for j in range(i - 1, max(i - 10, -1), -1):
                    s = lines[j].strip()
                    if not s:
                        continue
                    tm = tier_re.search(s)
                    if tm:
                        tier, desc = tm.group(1), tm.group(2)
                        break
                    if s.startswith("///") or s.startswith("// **退役"):
                        continue
                    break
                items.append((name, tier, desc))
            if items:
                by_file[rel] = items
        out = [
            "## Per-file function inventory (generated)",
            "",
            "Regenerate: `python3 scripts/annotate_api_tiers.py --inventory`",
            "",
        ]
        for f, items in by_file.items():
            out.append(f"### `{f}`")
            out.append("")
            out.append("| Function | Tier | Description |")
            out.append("|----------|------|-------------|")
            for name, tier, desc in items:
                out.append(f"| `{name}` | **{tier}** | {desc.replace('|', chr(92)+'|')} |")
            out.append("")
        inv = "\n".join(out)
        doc = doc_root / doc_name
        if not doc.exists():
            print(f"skip inventory (no doc): {doc_name}", file=sys.stderr)
            continue
        text = doc.read_text()
        marker = "## Per-file function inventory"
        if marker in text:
            text = text[: text.index(marker)] + inv.rstrip() + "\n"
        else:
            text = text.rstrip() + "\n\n---\n\n" + inv + "\n"
        doc.write_text(text)
        n = sum(len(v) for v in by_file.values())
        print(f"inventory {doc_name}: {n} functions")


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "--inventory":
        write_inventory_docs()
        return 0
    total = 0
    for src_rel, doc_name, skip in CRATE_TARGETS:
        src_root = ROOT / "crates" / src_rel
        if not src_root.is_dir():
            continue
        for path in sorted(src_root.rglob("*.rs")):
            if path.name in skip:
                continue
            add_module_inventory_header(path, doc_name)
            total += process_file(path)
    print(f"annotated {total} functions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
