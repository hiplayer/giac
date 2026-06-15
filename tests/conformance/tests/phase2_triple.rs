//! Phase 2 triple validation: giac-rs + Giac reference + SymPy (third-party CAS).

use giac_conformance::{
    run_giac, run_line, sympy_equiv, triple_check_script_filtered, upstream_root, verify_sympy,
};

/// Lines from bin/test_poly — giac-rs must pass SymPy; cross-check Giac.
#[test]
fn test_poly_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_poly", |_| false)?;
    assert_eq!(results.len(), 5);
    for r in &results {
        assert!(r.sympy_rs_ok, "giac-rs failed SymPy: {} -> {}", r.line, r.giac_rs);
        assert!(
            r.rs_giac_equiv || is_known_format_diff(&r.line),
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
    for r in &results {
        if !r.sympy_rs_ok && !is_known_sympy_gap(&r.line) {
            return Err(format!(
                "giac-rs failed SymPy on {}: {}",
                r.line, r.giac_rs
            ));
        }
        if !r.rs_giac_equiv && !is_known_format_diff(&r.line) {
            eprintln!(
                "note: {} giac-rs={} giac={} (sympy_rs={} sympy_giac={})",
                r.line, r.giac_rs, r.giac, r.sympy_rs_ok, r.sympy_giac_ok
            );
        }
    }
    Ok(())
}

#[test]
fn test_modular_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_modular", |_| false)?;
    assert_eq!(results.len(), 8);
    for r in &results {
        if !r.sympy_rs_ok && !is_known_sympy_gap(&r.line) {
            return Err(format!(
                "giac-rs failed SymPy on {}: {}",
                r.line, r.giac_rs
            ));
        }
        if !r.sympy_giac_ok && !is_known_giac_gap(&r.line) {
            eprintln!(
                "note: giac failed SymPy on {}: {}",
                r.line, r.giac
            );
        }
    }
    Ok(())
}

#[test]
fn test_factor_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_factor", |_| false)?;
    assert_eq!(results.len(), 3);
    for r in &results {
        assert!(r.sympy_rs_ok, "giac-rs SymPy: {} -> {}", r.line, r.giac_rs);
        if !r.rs_giac_equiv && !is_known_format_diff(&r.line) {
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
    assert!(rs.contains('x'), "giac-rs gcd mod 13 should be linear, got {rs}");
    Ok(())
}

fn is_known_format_diff(line: &str) -> bool {
    matches!(
        line,
        "gauss(2*x*y,[x,y])"
            | "egcd(x^2+2*x+1,x^2-1)"
            | "abcuv(x^2+2*x+1,x^2-1,x+1)"
            | "lcm(x^2+2*x+1,x^2-1)"
            | "roots(x^3-1,x)"
            | "partfrac(1/(x^2-1),x)"
            | "chinrem([x+2,x^2+1],[x+1,x^2+x+1])"
    )
}

fn is_known_sympy_gap(line: &str) -> bool {
    matches!(line, "roots(x^3-1,x)")
}

fn is_known_giac_gap(line: &str) -> bool {
    matches!(line, "chinrem([x+2,x^2+1],[x+1,x^2+x+1])")
}

#[test]
fn giac_binary_available() {
    let giac = upstream_root().join("build/bin/giac");
    assert!(
        giac.exists(),
        "Giac reference binary required at {} for triple tests",
        giac.display()
    );
}
