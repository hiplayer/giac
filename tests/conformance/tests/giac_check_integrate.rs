//! Phase 4: check golden harness (GIAC-219) — full `check/testintegrate` inventory.

use std::fs;
use std::path::Path;

use giac_conformance::{
    assert_check_integrate_non_integrate_sympy, assert_check_integrate_risch,
    classify_check_line, load_check_integrate_table, load_integrate_check_lines, run_line,
    run_lines, sympy_verify_lines, giac_check_dir, CheckLineKind,
};

#[test]
fn giac_check_integrate_files_exist() {
    assert!(giac_check_dir().join("testintegrate").exists());
}

fn enabled_lines() -> Result<Vec<String>, String> {
    Ok(load_check_integrate_table()?
        .entries
        .into_iter()
        .filter(|e| e.enabled)
        .map(|e| e.line)
        .collect())
}

/// Full inventory report (all lines, no assertion). Slow — run with `cargo test -- --ignored`.
#[test]
#[ignore = "runs SymPy on all testintegrate lines; use per-entry tests in CI"]
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

macro_rules! check_integrate_non_integrate {
    ($fn_name:ident, $id:literal) => {
        #[test]
        fn $fn_name() -> Result<(), String> {
            assert_check_integrate_non_integrate_sympy($id)
        }
    };
}

check_integrate_non_integrate!(giac_check_integrate_ck_int_55, "CK-INT-55");
check_integrate_non_integrate!(giac_check_integrate_ck_int_56, "CK-INT-56");
check_integrate_non_integrate!(giac_check_integrate_ck_int_57, "CK-INT-57");
check_integrate_non_integrate!(giac_check_integrate_ck_int_58, "CK-INT-58");
check_integrate_non_integrate!(giac_check_integrate_ck_int_59, "CK-INT-59");
check_integrate_non_integrate!(giac_check_integrate_ck_int_60, "CK-INT-60");
check_integrate_non_integrate!(giac_check_integrate_ck_int_61, "CK-INT-61");
check_integrate_non_integrate!(giac_check_integrate_ck_int_62, "CK-INT-62");
check_integrate_non_integrate!(giac_check_integrate_ck_int_63, "CK-INT-63");
check_integrate_non_integrate!(giac_check_integrate_ck_int_64, "CK-INT-64");
check_integrate_non_integrate!(giac_check_integrate_ck_int_65, "CK-INT-65");
check_integrate_non_integrate!(giac_check_integrate_ck_int_66, "CK-INT-66");

macro_rules! check_integrate_risch {
    ($fn_name:ident, $id:literal) => {
        #[test]
        fn $fn_name() -> Result<(), String> {
            assert_check_integrate_risch($id)
        }
    };
}

check_integrate_risch!(giac_check_risch_ck_int_02, "CK-INT-02");
check_integrate_risch!(giac_check_risch_ck_int_03, "CK-INT-03");
check_integrate_risch!(giac_check_risch_ck_int_04, "CK-INT-04");
check_integrate_risch!(giac_check_risch_ck_int_05, "CK-INT-05");
check_integrate_risch!(giac_check_risch_ck_int_06, "CK-INT-06");
check_integrate_risch!(giac_check_risch_ck_int_07, "CK-INT-07");
check_integrate_risch!(giac_check_risch_ck_int_08, "CK-INT-08");
check_integrate_risch!(giac_check_risch_ck_int_09, "CK-INT-09");
check_integrate_risch!(giac_check_risch_ck_int_11, "CK-INT-11");
check_integrate_risch!(giac_check_risch_ck_int_12, "CK-INT-12");
check_integrate_risch!(giac_check_risch_ck_int_13, "CK-INT-13");
check_integrate_risch!(giac_check_risch_ck_int_14, "CK-INT-14");
check_integrate_risch!(giac_check_risch_ck_int_18, "CK-INT-18");
check_integrate_risch!(giac_check_risch_ck_int_19, "CK-INT-19");
check_integrate_risch!(giac_check_risch_ck_int_20, "CK-INT-20");
check_integrate_risch!(giac_check_risch_ck_int_21, "CK-INT-21");
check_integrate_risch!(giac_check_risch_ck_int_22, "CK-INT-22");
check_integrate_risch!(giac_check_risch_ck_int_28, "CK-INT-28");
check_integrate_risch!(giac_check_risch_ck_int_29, "CK-INT-29");
check_integrate_risch!(giac_check_risch_ck_int_30, "CK-INT-30");
check_integrate_risch!(giac_check_risch_ck_int_32, "CK-INT-32");
check_integrate_risch!(giac_check_risch_ck_int_37, "CK-INT-37");
check_integrate_risch!(giac_check_risch_ck_int_43, "CK-INT-43");

/// Full SymPy gate including integrate rows — manual only (heavy cases may timeout).
#[test]
#[ignore = "integrate SymPy can hang; use per-entry tests in CI"]
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

#[derive(Debug, serde::Deserialize)]
struct IntegrateTable {
    entries: Vec<IntegrateEntry>,
}

#[derive(Debug, serde::Deserialize)]
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
