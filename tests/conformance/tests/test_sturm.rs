//! GIAC-206: `bin/test_sturm*` SymPy / property verification.

use giac_conformance::{run_lines, script_lines, sympy_verify_lines, upstream_root};

#[test]
fn test_sturm_sympy() -> Result<(), String> {
    let path = upstream_root().join("bin/test_sturm");
    let lines = script_lines(&path)?;
    let outputs = run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
    assert!(failures.is_empty(), "test_sturm failures: {failures:?}");
    Ok(())
}

#[test]
fn test_sturm_ext_sturm_lines_sympy() -> Result<(), String> {
    let path = upstream_root().join("bin/test_sturm_ext");
    let lines: Vec<String> = script_lines(&path)?
        .into_iter()
        .filter(|l| l.starts_with("sturm("))
        .collect();
    let outputs = run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
    assert!(failures.is_empty(), "test_sturm_ext sturm failures: {failures:?}");
    Ok(())
}
