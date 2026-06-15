//! Phase 4 standard integral table — SymPy derivative verification (§2.7).

use std::fs;
use std::path::PathBuf;

use giac_conformance::{run_line, verify_sympy};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct IntegrateTable {
    entries: Vec<IntegrateEntry>,
}

#[derive(Debug, Deserialize)]
struct IntegrateEntry {
    id: String,
    line: String,
    enabled: bool,
}

fn fixture_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/phase4_integrate_table.json")
}

use std::path::Path;

fn load_table() -> IntegrateTable {
    let text = fs::read_to_string(fixture_path()).expect("phase4_integrate_table.json");
    serde_json::from_str(&text).expect("valid integrate table JSON")
}

#[test]
fn phase4_integrate_table_enabled_sympy() -> Result<(), String> {
    let table = load_table();
    let enabled: Vec<_> = table.entries.iter().filter(|e| e.enabled).collect();
    assert!(!enabled.is_empty(), "no enabled integral table entries");

    let mut failures = Vec::new();
    for entry in enabled {
        let got = match run_line(&entry.line) {
            Ok(out) => out,
            Err(e) => {
                failures.push((entry.id.clone(), entry.line.clone(), e));
                continue;
            }
        };
        if let Err(e) = verify_sympy(&entry.line, &got) {
            failures.push((entry.id.clone(), got, e));
        }
    }
    assert!(
        failures.is_empty(),
        "integral table SymPy failures: {:?}",
        failures
    );
    Ok(())
}

#[test]
fn phase4_integrate_table_json_has_unique_ids() {
    let table = load_table();
    let mut ids = std::collections::HashSet::new();
    for e in &table.entries {
        assert!(ids.insert(e.id.clone()), "duplicate id {}", e.id);
    }
}
