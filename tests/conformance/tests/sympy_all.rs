//! SymPy third-party verification for Phase 1–3 (all conformance scripts).

use giac_conformance::{
    assert_sympy_script, assert_sympy_script_lines, assert_testcas_sympy_range, run_line,
    sympy_verify_line, verify_sympy, ALL_SYMPTY_SCRIPTS, PHASE1_SCRIPTS, PHASE2_SCRIPTS,
    PHASE3_SCRIPTS, upstream_root, run_script,
};

#[test]
fn sympy_available() {
    let script = giac_conformance::sympy_script();
    assert!(script.exists(), "sympy_verify.py required at {}", script.display());
}

macro_rules! sympy_script_test {
    ($fn_name:ident, $script:literal, $strict:expr) => {
        #[test]
        fn $fn_name() -> Result<(), String> {
            if $strict {
                assert_sympy_script($script, |_| false)
            } else {
                assert_sympy_script($script, is_known_sympy_gap)
            }
        }
    };
}

sympy_script_test!(phase1_test_cas_basic_sympy, "test_cas_basic", true);

sympy_script_test!(phase2_test_poly_sympy, "test_poly", false);
sympy_script_test!(phase2_test_poly_ext_sympy, "test_poly_ext", false);
sympy_script_test!(phase2_test_modular_sympy, "test_modular", false);
sympy_script_test!(phase2_test_factor_sympy, "test_factor", false);
sympy_script_test!(phase2_test_groebner_sympy, "test_groebner", false);

macro_rules! phase3_script_test {
    ($fn_name:ident, $script:literal) => {
        #[test]
        fn $fn_name() -> Result<(), String> {
            assert_sympy_script_lines($script, is_phase3_skip, is_known_sympy_gap)
        }
    };
}

phase3_script_test!(phase3_test_linalg_sympy, "test_linalg");
phase3_script_test!(phase3_test_linalg_ext_sympy, "test_linalg_ext");
phase3_script_test!(phase3_test_linalg_decomp_sympy, "test_linalg_decomp");
phase3_script_test!(phase3_test_gauss_ext_sympy, "test_gauss_ext");

#[test]
fn test_cas_basic_each_line_sympy() -> Result<(), String> {
    let cases = [
        ("abs(1+2*i)", "sqrt(5)"),
        ("gcd(45,75)", "15"),
        ("conj((1+2*i)^2)", "-3-4*i"),
    ];
    for (line, want) in cases {
        let got = run_line(line)?;
        assert_eq!(got, want, "{line}");
        verify_sympy(line, &got)?;
    }
    Ok(())
}

macro_rules! testcas_sympy_chunk {
    ($fn_name:ident, $start:expr, $end:expr, $min:expr) => {
        #[test]
        fn $fn_name() -> Result<(), String> {
            assert_testcas_sympy_range($start, $end, is_known_testcas_gap, $min)
        }
    };
}

testcas_sympy_chunk!(testcas_sympy_00_09, 0, 10, 8);
testcas_sympy_chunk!(testcas_sympy_10_19, 10, 20, 8);
testcas_sympy_chunk!(testcas_sympy_20_29, 20, 30, 8);
testcas_sympy_chunk!(testcas_sympy_30_39, 30, 40, 8);
testcas_sympy_chunk!(testcas_sympy_40_49, 40, 50, 8);

macro_rules! bin_script_smoke {
    ($fn_name:ident, $script:literal) => {
        #[test]
        fn $fn_name() -> Result<(), String> {
            let path = upstream_root().join("bin").join($script);
            let outputs = run_script(&path)?;
            let lines = giac_conformance::script_lines(&path)?;
            assert_eq!(lines.len(), outputs.len());
            for (line, out) in lines.iter().zip(outputs.iter()) {
                let r = sympy_verify_line(line, out)?;
                if !r.ok && !is_known_sympy_gap(line) && !is_known_testcas_gap(line) {
                    return Err(format!("SymPy: `{line}` -> `{out}`"));
                }
            }
            Ok(())
        }
    };
}

bin_script_smoke!(smoke_test_cas_basic_sympy, "test_cas_basic");
bin_script_smoke!(smoke_test_poly_sympy, "test_poly");
bin_script_smoke!(smoke_test_poly_ext_sympy, "test_poly_ext");
bin_script_smoke!(smoke_test_modular_sympy, "test_modular");
bin_script_smoke!(smoke_test_factor_sympy, "test_factor");
bin_script_smoke!(smoke_test_groebner_sympy, "test_groebner");

// Compile-time guard: smoke tests cover the same scripts as ALL_SYMPTY_SCRIPTS.
const _: () = assert!(ALL_SYMPTY_SCRIPTS.len() == 6);
const _: () = assert!(PHASE1_SCRIPTS.len() == 1);
const _: () = assert!(PHASE2_SCRIPTS.len() == 5);
const _: () = assert!(PHASE3_SCRIPTS.len() == 4);

fn is_known_sympy_gap(line: &str) -> bool {
    matches!(line, "roots(x^3-1,x)")
}

fn is_phase3_skip(_line: &str) -> bool {
    false
}

fn is_known_testcas_gap(line: &str) -> bool {
    matches!(
        line,
        "integrate(1/(1-x^4),x)"
            | "normal(-2+-3/-2*2)"
            | "normal(-6+-3/-2*4)"
    ) || line.parse::<i64>().is_ok()
}
