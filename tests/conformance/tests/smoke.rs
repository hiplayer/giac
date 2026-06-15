//! Golden regression harness for giac-rs (Phase 1 + SymPy).

use giac_conformance::{run_script, sympy_verify_script, upstream_root, verify_sympy, run_line, PHASE1_SCRIPTS};
use giac_core::Context;
use giac_parse::parse_program;
use std::fs;

#[test]
fn test_cas_basic_matches_giac() -> Result<(), String> {
    let script = upstream_root().join("bin/test_cas_basic");
    let got = run_script(&script)?;
    let want = vec![
        "sqrt(5)".to_string(),
        "15".to_string(),
        "-3-4*i".to_string(),
    ];
    assert_eq!(got, want, "bin/test_cas_basic");
    let lines = giac_conformance::script_lines(&script)?;
    for (line, out) in lines.iter().zip(got.iter()) {
        verify_sympy(line, out)?;
    }
    Ok(())
}

#[test]
fn phase1_scripts_sympy() -> Result<(), String> {
    for name in PHASE1_SCRIPTS {
        let results = sympy_verify_script(name)?;
        for r in results {
            assert!(r.ok, "SymPy: `{}` -> `{}`", r.line, r.output);
        }
    }
    Ok(())
}

#[test]
fn harness_phase0_bin_scripts_parseable() {
    for name in PHASE1_SCRIPTS {
        let path = upstream_root().join("bin").join(name);
        let input = fs::read_to_string(&path).unwrap();
        let ctx = Context::xcas_default();
        assert!(
            parse_program(&input, &ctx).is_ok(),
            "parse failed for {name}"
        );
    }
}

#[test]
fn phase1_complex_gcd_sympy() -> Result<(), String> {
    let line = "gcd(1999,2001)";
    let got = run_line(line)?;
    verify_sympy(line, &got)?;
    Ok(())
}
