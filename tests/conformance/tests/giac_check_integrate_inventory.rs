//! One-off inventory: disabled check_integrate rows that already pass SymPy.

use std::fs;
use std::path::Path;

use giac_conformance::{run_line, sympy_verify_lines};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct Table {
    entries: Vec<Entry>,
}

#[derive(Debug, Deserialize)]
struct Entry {
    id: String,
    line: String,
    enabled: bool,
}

#[test]
#[ignore = "inventory helper; run manually to find rows to enable"]
fn giac_check_integrate_disabled_sympy_green() -> Result<(), String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/check_integrate_table.json");
    let text = fs::read_to_string(&path).map_err(|e| e.to_string())?;
    let table: Table = serde_json::from_str(&text).map_err(|e| e.to_string())?;
    for e in table.entries.iter().filter(|e| !e.enabled) {
        let Ok(got) = run_line(&e.line) else {
            eprintln!("{} eval fail: {}", e.id, e.line);
            continue;
        };
        let results = sympy_verify_lines(&[e.line.clone()], &[got])?;
        if results[0].ok {
            eprintln!("GREEN disabled: {} {}", e.id, e.line);
        }
    }
    Ok(())
}
