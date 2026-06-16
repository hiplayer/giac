//! Phase 4 Batch 4 acceptance (GIAC-216, 217, 218).

use giac_conformance::{run_line, script_lines, sympy_verify_lines, upstream_root, verify_sympy};

#[test]
fn batch4_giac216_series_taylor() -> Result<(), String> {
    for line in [
        "taylor(sin(x),x=0,5)",
        "taylor(exp(x),x=0,4)",
        "series(exp(x),x,0,4)",
    ] {
        let got = run_line(line)?;
        verify_sympy(line, &got)?;
    }
    Ok(())
}

#[test]
fn batch4_giac216_test_series_bin() -> Result<(), String> {
    let path = upstream_root().join("bin/test_series");
    let lines = script_lines(&path)?;
    let outputs = giac_conformance::run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
    assert!(failures.is_empty(), "test_series failures: {failures:?}");
    Ok(())
}

#[test]
fn batch4_giac217_risch() -> Result<(), String> {
    let line = "risch(exp(x)*cos(x),x)";
    let got = run_line(line)?;
    verify_sympy("integrate(exp(x)*cos(x),x)", &got)?;
    Ok(())
}

#[test]
fn batch4_giac218_desolve_basics() -> Result<(), String> {
    for line in [
        "desolve(y''+y=0,y(x))",
        "desolve(y'=x*y,y(x))",
        "desolve(y''-3*y'+2*y=0,y(x))",
        "desolve(y''-2*y'+y=0,y(x))",
        "desolve(y'+y=x,y(x))",
        "desolve(y''+4*y=sin(x),y(x))",
    ] {
        let got = run_line(line)?;
        verify_sympy(line, &got)?;
    }
    Ok(())
}

#[test]
fn batch4_giac218_test_desolve_bins() -> Result<(), String> {
    for name in ["test_desolve", "test_desolve_ext"] {
        let path = upstream_root().join("bin").join(name);
        let lines = script_lines(&path)?;
        let outputs = giac_conformance::run_lines(&lines)?;
        let results = sympy_verify_lines(&lines, &outputs)?;
        let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
        assert!(failures.is_empty(), "{name} failures: {failures:?}");
    }
    Ok(())
}
