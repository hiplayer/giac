//! Phase 4: check golden pilot harness (GIAC-219).

use std::fs;
use std::path::Path;

use giac_conformance::{run_line, sympy_verify_lines, upstream_root};
use serde::Deserialize;

#[test]
fn giac_check_integrate_files_exist() {
    let root = upstream_root();
    assert!(root.join("giac/giac-1.5.0/check/testintegrate").exists());
}

/// Pilot: first 20 `check/testintegrate` lines (report only — most need future `risch`).
#[test]
fn giac_check_integrate_pilot_report() -> Result<(), String> {
    let path = upstream_root().join("giac/giac-1.5.0/check/testintegrate");
    let text = fs::read_to_string(&path).map_err(|e| format!("read: {e}"))?;
    let lines: Vec<String> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with("cas_setup"))
        .take(20)
        .map(|l| l.replace("**", "^").trim_end_matches(';').to_string())
        .collect();
    let mut eval_ok = 0usize;
    let mut sympy_ok = 0usize;
    for line in &lines {
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
    eprintln!("giac_check_integrate pilot: {sympy_ok}/{eval_ok} SymPy ({}/20 lines)", lines.len());
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

/// SymPy derivative gate on enabled integral-table rows (proxy for check golden T1).
#[test]
fn giac_check_integrate_table_sympy() -> Result<(), String> {
    let lines = enabled_table_lines(20)?;
    assert!(!lines.is_empty());
    let outputs = giac_conformance::run_lines(&lines)?;
    let results = sympy_verify_lines(&lines, &outputs)?;
    let failures: Vec<_> = results.iter().filter(|r| !r.ok).collect();
    assert!(
        failures.is_empty(),
        "integrate table SymPy failures: {failures:?}"
    );
    Ok(())
}
