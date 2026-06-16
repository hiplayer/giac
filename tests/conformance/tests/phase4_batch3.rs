//! Phase 4 Batch 3 acceptance (GIAC-207, 208, 212, 213, 215).

use giac_conformance::{run_line, script_lines, sympy_verify_lines, upstream_root, verify_sympy};

#[test]
fn batch3_giac212_integrate_partfrac() -> Result<(), String> {
    for line in [
        "integrate(1/(x^2-1),x)",
        "integrate(1/(1-x^2),x)",
        "integrate(1/(x^3+1),x)",
        "integrate(1/(x^4-1),x)",
    ] {
        let got = run_line(line)?;
        verify_sympy(line, &got)?;
    }
    Ok(())
}

#[test]
fn batch3_giac213_integrate_exp_trig() -> Result<(), String> {
    for line in [
        "integrate(exp(x)*sin(x),x)",
        "integrate(exp(x)*cos(x),x)",
    ] {
        let got = run_line(line)?;
        verify_sympy(line, &got)?;
    }
    Ok(())
}

#[test]
fn batch3_giac215_limit_basics() -> Result<(), String> {
    for line in [
        "limit(sin(x)/x,x,0)",
        "limit((1+1/x)^x,x,+infinity)",
        "limit((1-cos(x))/x^2,x,0)",
    ] {
        let got = run_line(line)?;
        verify_sympy(line, &got)?;
    }
    Ok(())
}

#[test]
fn batch3_giac215_test_limit_bin() -> Result<(), String> {
    let path = upstream_root().join("bin/test_limit");
    let lines = script_lines(&path)?;
    let outputs = giac_conformance::run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
    assert!(failures.is_empty(), "test_limit failures: {failures:?}");
    Ok(())
}

#[test]
fn batch3_giac206_sturm_odd_multiplicity() -> Result<(), String> {
    let got = run_line("sturm((x^3+1)^2)")?;
    verify_sympy("sturm((x^3+1)^2)", &got)?;
    Ok(())
}

#[test]
fn batch3_giac207_realroot() -> Result<(), String> {
    let got = run_line("realroot(x^4-1)")?;
    assert!(got.contains("-1"), "expected -1 in {got}");
    assert!(got.contains("1"), "expected 1 in {got}");
    Ok(())
}

#[test]
fn batch3_giac208_solve_sin() -> Result<(), String> {
    let got = run_line("solve(sin(x)=0,x)")?;
    assert!(got.contains("0"), "expected 0 in {got}");
    Ok(())
}
