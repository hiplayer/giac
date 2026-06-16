//! Phase 4: check golden harness (GIAC-219) — full `check/testintegrate` inventory.

use std::fs;
use std::path::Path;

use giac_conformance::{
    classify_check_line, load_integrate_check_lines, run_line, run_lines, run_risch_line,
    sympy_equiv, sympy_verify_lines, giac_check_dir, CheckLineKind,
};
use serde::Deserialize;

#[test]
fn giac_check_integrate_files_exist() {
    assert!(giac_check_dir().join("testintegrate").exists());
}

#[derive(Debug, Deserialize)]
struct CheckIntegrateTable {
    entries: Vec<CheckIntegrateEntry>,
}

#[derive(Debug, Deserialize)]
struct CheckIntegrateEntry {
    id: String,
    line: String,
    kind: String,
    enabled: bool,
}

fn check_integrate_table() -> Result<CheckIntegrateTable, String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/check_integrate_table.json");
    let text = fs::read_to_string(&path).map_err(|e| format!("read {}: {e}", path.display()))?;
    serde_json::from_str(&text).map_err(|e| format!("json: {e}"))
}

fn enabled_lines() -> Result<Vec<String>, String> {
    Ok(check_integrate_table()?
        .entries
        .into_iter()
        .filter(|e| e.enabled)
        .map(|e| e.line)
        .collect())
}

/// Full inventory report (all lines, no assertion). Slow — run with `cargo test -- --ignored`.
#[test]
#[ignore = "runs SymPy on all testintegrate lines; use giac_check_integrate_enabled_sympy in CI"]
fn giac_check_integrate_full_report() -> Result<(), String> {
    let lines = load_integrate_check_lines()?;
    let mut by_kind = [0usize; 5];
    let mut eval_ok = 0usize;
    let mut sympy_ok = 0usize;
    for line in &lines {
        let kind = classify_check_line(line);
        let idx = match kind {
            CheckLineKind::Integrate => 0,
            CheckLineKind::Limit => 1,
            CheckLineKind::Series => 2,
            CheckLineKind::Compound => 3,
            CheckLineKind::Other => 4,
        };
        by_kind[idx] += 1;
        let got = match run_line(line) {
            Ok(out) => out,
            Err(_) => continue,
        };
        eval_ok += 1;
        let results = sympy_verify_lines(std::slice::from_ref(line), std::slice::from_ref(&got))?;
        if results[0].ok {
            sympy_ok += 1;
        }
    }
    eprintln!(
        "giac_check_integrate full: {}/{} eval, {}/{} SymPy (integrate/limit/series/compound/other = {:?}/{})",
        eval_ok,
        lines.len(),
        sympy_ok,
        lines.len(),
        by_kind,
        lines.len()
    );
    Ok(())
}

/// Regression gate on `fixtures/check_integrate_table.json` enabled rows.
///
/// Integrate rows: eval only (SymPy `diff` verification can hang on heavy rationals).
/// Limit / series / other: SymPy T1 as before.
#[test]
fn giac_check_integrate_enabled() -> Result<(), String> {
    let table = check_integrate_table()?;
    let enabled: Vec<_> = table.entries.iter().filter(|e| e.enabled).collect();
    assert!(!enabled.is_empty(), "no enabled check_integrate rows");

    let mut sympy_lines = Vec::new();
    let mut sympy_outputs = Vec::new();
    for e in enabled {
        let got = run_line(&e.line)?;
        if e.kind == "integrate" {
            continue;
        }
        sympy_lines.push(e.line.clone());
        sympy_outputs.push(got);
    }
    if sympy_lines.is_empty() {
        return Ok(());
    }
    let results = sympy_verify_lines(&sympy_lines, &sympy_outputs)?;
    let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
    assert!(
        failures.is_empty(),
        "check_integrate SymPy failures (non-integrate): {failures:?}"
    );
    Ok(())
}

/// Full SymPy gate including integrate rows — manual only (heavy cases may timeout).
#[test]
#[ignore = "integrate SymPy can hang; use giac_check_integrate_enabled in CI"]
fn giac_check_integrate_enabled_sympy() -> Result<(), String> {
    let lines = enabled_lines()?;
    assert!(!lines.is_empty(), "no enabled check_integrate rows");
    let outputs = run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
    assert!(
        failures.is_empty(),
        "check_integrate SymPy failures: {failures:?}"
    );
    Ok(())
}

/// `risch(f,x)` agrees with `integrate(f,x)` on enabled integrate rows (GIAC-217).
#[test]
fn giac_check_risch_matches_integrate() -> Result<(), String> {
    let table = check_integrate_table()?;
    let integrate_lines: Vec<_> = table
        .entries
        .iter()
        .filter(|e| e.enabled && e.kind == "integrate")
        .map(|e| e.line.as_str())
        .collect();
    assert!(!integrate_lines.is_empty());
    for line in integrate_lines {
        let int_out = run_line(line)?;
        let risch_out = run_risch_line(line)?;
        sympy_equiv(&int_out, &risch_out).map_err(|e| format!("{line}: {e}"))?;
    }
    Ok(())
}

#[derive(Debug, Deserialize)]
struct IntegrateTable {
    entries: Vec<IntegrateEntry>,
}

#[derive(Debug, Deserialize)]
struct IntegrateEntry {
    line: String,
    enabled: bool,
}

fn enabled_table_lines(n: usize) -> Result<Vec<String>, String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/phase4_integrate_table.json");
    let text = fs::read_to_string(&path).map_err(|e| format!("read table: {e}"))?;
    let table: IntegrateTable = serde_json::from_str(&text).map_err(|e| format!("json: {e}"))?;
    Ok(table
        .entries
        .into_iter()
        .filter(|e| e.enabled)
        .take(n)
        .map(|e| e.line)
        .collect())
}

/// Eval gate on enabled integral-table rows (SymPy derivative check disabled — see §GIAC-219).
#[test]
fn giac_check_integrate_table_enabled() -> Result<(), String> {
    let lines = enabled_table_lines(20)?;
    assert!(!lines.is_empty());
    run_lines(&lines)?;
    Ok(())
}

/// SymPy derivative gate — manual only.
#[test]
#[ignore = "integrate SymPy can hang; use giac_check_integrate_table_enabled in CI"]
fn giac_check_integrate_table_sympy() -> Result<(), String> {
    let lines = enabled_table_lines(20)?;
    assert!(!lines.is_empty());
    let outputs = run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
    assert!(
        failures.is_empty(),
        "integrate table SymPy failures: {failures:?}"
    );
    Ok(())
}
