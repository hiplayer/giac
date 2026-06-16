//! One-off inventory: disabled check_integrate rows that already pass SymPy.

use std::fs;
use std::io::Write;
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
    kind: String,
    enabled: bool,
}

#[derive(Debug, serde::Serialize)]
struct InventoryRow {
    id: String,
    line: String,
    kind: String,
    status: &'static str,
}

#[test]
#[ignore = "inventory helper: eval only, no SymPy"]
fn giac_check_integrate_disabled_eval_only() -> Result<(), String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/check_integrate_table.json");
    let text = fs::read_to_string(&path).map_err(|e| e.to_string())?;
    let table: Table = serde_json::from_str(&text).map_err(|e| e.to_string())?;
    let mut ok = 0usize;
    for e in table.entries.iter().filter(|e| !e.enabled) {
        match run_line(&e.line) {
            Ok(_) => {
                eprintln!("EVAL OK: {} {}", e.id, e.line);
                ok += 1;
            }
            Err(err) => eprintln!("EVAL FAIL: {} {} ({err})", e.id, e.line),
        }
    }
    eprintln!("eval inventory: {}/{} ok", ok, table.entries.iter().filter(|e| !e.enabled).count());
    Ok(())
}

#[test]
#[ignore = "inventory helper; run manually to find rows to enable"]
fn giac_check_integrate_disabled_sympy_green() -> Result<(), String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/check_integrate_table.json");
    let text = fs::read_to_string(&path).map_err(|e| e.to_string())?;
    let table: Table = serde_json::from_str(&text).map_err(|e| e.to_string())?;
    let out_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/ck_int_inventory.jsonl");
    let mut out = fs::File::create(&out_path).map_err(|e| e.to_string())?;
    let mut rows = Vec::new();
    for e in table.entries.iter().filter(|e| !e.enabled) {
        let status = match run_line(&e.line) {
            Ok(got) => {
                let results = sympy_verify_lines(&[e.line.clone()], &[got])?;
                if results[0].ok {
                    "green"
                } else {
                    "sympy_fail"
                }
            }
            Err(_) => "eval_fail",
        };
        if status == "green" {
            eprintln!("GREEN disabled: {} {}", e.id, e.line);
        } else if status == "eval_fail" {
            eprintln!("{} eval fail: {}", e.id, e.line);
        }
        let row = InventoryRow {
            id: e.id.clone(),
            line: e.line.clone(),
            kind: e.kind.clone(),
            status,
        };
        serde_json::to_writer(&mut out, &row).map_err(|e| e.to_string())?;
        writeln!(out).map_err(|e| e.to_string())?;
        out.flush().ok();
        rows.push(row);
    }
    let green: Vec<_> = rows.iter().filter(|r| r.status == "green").collect();
    eprintln!(
        "inventory: {} disabled, {} green -> {}",
        rows.len(),
        green.len(),
        out_path.display()
    );
    Ok(())
}
