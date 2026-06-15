//! GIAC-112: `bin/test_trig` parse + SymPy / equivalence verification.

use giac_conformance::{run_lines, script_lines, sympy_verify_lines, upstream_root};
use giac_calculus::xcas_default;
use giac_parse::parse_program;

const MAX_LINES: usize = 10;
const MIN_PASS: usize = 8;

#[test]
fn test_trig_parse_full_file() -> Result<(), String> {
    let path = upstream_root().join("bin/test_trig");
    let input = std::fs::read_to_string(&path)
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    parse_program(&input, &xcas_default()).map_err(|e| format!("parse test_trig: {e}"))?;
    Ok(())
}

#[test]
fn test_trig_sympy_subset() -> Result<(), String> {
    let path = upstream_root().join("bin/test_trig");
    let lines: Vec<String> = script_lines(&path)?
        .into_iter()
        .take(MAX_LINES)
        .collect();
    let outputs = run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let pass = results.iter().filter(|r| r.ok).count();
    let required = MIN_PASS.min(lines.len());
    let failures: Vec<_> = results
        .iter()
        .filter(|r| !r.ok)
        .map(|r| (r.line.clone(), r.output.clone()))
        .collect();
    assert!(
        pass >= required,
        "test_trig: {pass}/{} SymPy-verified, need {required}; failures: {failures:?}",
        lines.len()
    );
    Ok(())
}
