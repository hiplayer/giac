//! GIAC-114: `bin/test_solve` parse + end-to-end verification.

use giac_conformance::{run_line, run_lines, script_lines, sympy_verify_lines, upstream_root, verify_sympy};
use giac_calculus::xcas_default;
use giac_parse::parse_program;

#[test]
fn test_solve_parse_full_file() -> Result<(), String> {
    let path = upstream_root().join("bin/test_solve");
    let input = std::fs::read_to_string(&path)
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    parse_program(&input, &xcas_default()).map_err(|e| format!("parse test_solve: {e}"))?;
    Ok(())
}

#[test]
fn test_solve_quadratic_double_root() -> Result<(), String> {
    let got = run_line("solve(x^2-2*x+1=0,x)")?;
    verify_sympy("solve(x^2-2*x+1=0,x)", &got)
}

#[test]
fn test_solve_at_least_one_line() -> Result<(), String> {
    let path = upstream_root().join("bin/test_solve");
    let lines: Vec<String> = script_lines(&path)?;
    let outputs = run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let pass = results.iter().filter(|r| r.ok).count();
    assert!(
        pass >= 1,
        "test_solve: need >=1 line passing, got {pass}/{}; failures: {:?}",
        lines.len(),
        results
            .iter()
            .filter(|r| !r.ok)
            .map(|r| (&r.line, &r.output))
            .collect::<Vec<_>>()
    );
    Ok(())
}
