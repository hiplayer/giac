//! Phase 2 triple validation: giac-rs + Giac reference + SymPy (third-party CAS).

use giac_conformance::{
    phase2_format_diff, phase2_giac_gap, phase2_sympy_gap, run_giac, run_line, sympy_equiv,
    triple_assert_sympy_rs, triple_check_script_filtered, triple_check_script_line,
    triple_note_format_diffs, giac_binary, verify_sympy,
};

/// Lines from bin/test_poly — giac-rs must pass SymPy; cross-check Giac.
#[test]
fn test_poly_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_poly", |_| false)?;
    assert_eq!(results.len(), 5);
    for r in &results {
        assert!(r.sympy_rs_ok, "giac-rs failed SymPy: {} -> {}", r.line, r.giac_rs);
        assert!(
            r.rs_giac_equiv || phase2_format_diff(&r.line),
            "giac-rs `{}` = `{}` differs from giac `{}`",
            r.line,
            r.giac_rs,
            r.giac
        );
    }
    Ok(())
}

/// Extended polynomial ops — SymPy property checks on giac-rs; format diffs vs giac allowed.
#[test]
fn test_poly_ext_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_poly_ext", |_| false)?;
    assert_eq!(results.len(), 7);
    triple_assert_sympy_rs(&results, phase2_sympy_gap)?;
    triple_note_format_diffs(&results, phase2_format_diff);
    Ok(())
}

macro_rules! test_modular_line_triple {
    ($fn_name:ident, $idx:literal) => {
        #[test]
        fn $fn_name() -> Result<(), String> {
            let r = triple_check_script_line("test_modular", $idx, |_| false)?;
            if !r.sympy_rs_ok && !phase2_sympy_gap(&r.line) {
                return Err(format!(
                    "giac-rs failed SymPy on {}: {}",
                    r.line, r.giac_rs
                ));
            }
            if !r.sympy_giac_ok && !phase2_giac_gap(&r.line) {
                eprintln!(
                    "note: giac failed SymPy on {}: {}",
                    r.line, r.giac
                );
            }
            Ok(())
        }
    };
}

test_modular_line_triple!(test_modular_line_0_triple, 0);
test_modular_line_triple!(test_modular_line_1_triple, 1);
test_modular_line_triple!(test_modular_line_2_triple, 2);
test_modular_line_triple!(test_modular_line_3_triple, 3);
test_modular_line_triple!(test_modular_line_4_triple, 4);
test_modular_line_triple!(test_modular_line_5_triple, 5);
test_modular_line_triple!(test_modular_line_6_triple, 6);
test_modular_line_triple!(test_modular_line_7_triple, 7);

#[test]
fn test_factor_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_factor", |_| false)?;
    assert_eq!(results.len(), 3);
    for r in &results {
        assert!(r.sympy_rs_ok, "giac-rs SymPy: {} -> {}", r.line, r.giac_rs);
        if !r.rs_giac_equiv && !phase2_format_diff(&r.line) {
            eprintln!("note: {} giac-rs={} giac={}", r.line, r.giac_rs, r.giac);
        }
    }
    Ok(())
}

/// greduce: Giac + SymPy agree; giac-rs MVP may not reduce (known gap).
#[test]
fn test_groebner_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_groebner", |_| false)?;
    assert_eq!(results.len(), 2);
    for r in &results {
        assert!(
            r.sympy_rs_ok,
            "giac-rs greduce failed SymPy on {}: {}",
            r.line,
            r.giac_rs
        );
        assert!(
            r.sympy_giac_ok,
            "Giac should match SymPy on {}: {}",
            r.line,
            r.giac
        );
    }
    Ok(())
}

#[test]
fn factor_x4_minus_1_all_three() -> Result<(), String> {
    let line = "factor(x^4-1)";
    let rs = run_line(line)?;
    let giac = run_giac(line)?;
    verify_sympy(line, &rs)?;
    verify_sympy(line, &giac)?;
    sympy_equiv(&rs, &giac)?;
    assert_eq!(rs, "(x-1)*(x+1)*(x^2+1)");
    Ok(())
}

#[test]
fn gcd_mod_13_triple() -> Result<(), String> {
    let line = "gcd((2*x^2+5) % 13,(5*x^2+2*x-3) % 13)";
    let rs = run_line(line)?;
    verify_sympy(line, &rs)?;
    let giac = run_giac(line)?;
    verify_sympy(line, &giac)?;
    assert_eq!(rs, "7*x+1 mod 13");
    Ok(())
}

#[test]
fn giac_binary_available() {
    let giac = giac_binary();
    assert!(
        giac.exists(),
        "Giac reference binary required at {} for triple tests",
        giac.display()
    );
}
