//! Phase 4 Batch 1 acceptance (GIAC-201–204, GIAC-210).

use std::fs;
use std::path::Path;

use giac_conformance::{run_line, sympy_script, upstream_root, verify_sympy};
use giac_ode::xcas_default;
use giac_parse::parse_program;
use serde::Deserialize;
use std::process::Command;

#[derive(Debug, Deserialize)]
struct IntegrateTable {
    entries: Vec<IntegrateEntry>,
}

#[derive(Debug, Deserialize)]
struct IntegrateEntry {
    id: String,
    line: String,
    enabled: bool,
    giac_issue: Option<String>,
}

fn integrate_table() -> IntegrateTable {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/phase4_integrate_table.json");
    let text = fs::read_to_string(path).expect("phase4_integrate_table.json");
    serde_json::from_str(&text).expect("valid integrate table JSON")
}

fn sympy_supported(line: &str) -> Result<(), String> {
    let out = Command::new("python3")
        .arg(sympy_script())
        .arg("supported")
        .arg(line)
        .output()
        .map_err(|e| format!("spawn sympy supported: {e}"))?;
    if out.status.success() {
        Ok(())
    } else {
        Err(String::from_utf8_lossy(&out.stderr).into_owned())
    }
}

#[test]
fn batch1_giac201_sympy_supported_limit_series() -> Result<(), String> {
    for line in [
        "limit(sin(x)/x,x,0)",
        "series(exp(x),x,0,4)",
        "limit((1+1/x)^x,x,+infinity)",
    ] {
        sympy_supported(line)?;
    }
    Ok(())
}

#[test]
fn batch1_giac204_parse_infinity_and_taylor_x0() -> Result<(), String> {
    let ctx = xcas_default();
    parse_program("limit((1+1/x)^x,x,+infinity);", &ctx)
        .map_err(|e| format!("parse +infinity: {e}"))?;
    parse_program("taylor(sin(x),x=0,5);", &ctx).map_err(|e| format!("parse x=0: {e}"))?;
    parse_program("series(exp(x),x,0,4);", &ctx).map_err(|e| format!("parse series: {e}"))?;
    Ok(())
}

#[test]
fn batch1_giac202_desolve_bin_scripts_parse() -> Result<(), String> {
    let ctx = xcas_default();
    for name in ["test_desolve", "test_desolve_ext"] {
        let path = upstream_root().join("bin").join(name);
        let input = std::fs::read_to_string(&path)
            .map_err(|e| format!("read {}: {e}", path.display()))?;
        parse_program(&input, &ctx).map_err(|e| format!("parse {name}: {e}"))?;
    }
    Ok(())
}

#[test]
fn batch1_giac203_integrate_bins_parse() -> Result<(), String> {
    let ctx = xcas_default();
    for name in ["test_integrate", "test_integrate_more", "test_integrate_ext"] {
        let path = upstream_root().join("bin").join(name);
        let input = std::fs::read_to_string(&path)
            .map_err(|e| format!("read {}: {e}", path.display()))?;
        parse_program(&input, &ctx).map_err(|e| format!("parse {name}: {e}"))?;
    }
    Ok(())
}

#[test]
fn batch1_giac210_table_entries_sympy() -> Result<(), String> {
    let table = integrate_table();
    let g210: Vec<_> = table
        .entries
        .iter()
        .filter(|e| e.enabled && e.giac_issue.as_deref() == Some("GIAC-210"))
        .collect();
    assert!(
        g210.len() >= 4,
        "expected at least 4 enabled GIAC-210 integral table rows, got {}",
        g210.len()
    );
    let mut failures = Vec::new();
    for entry in g210 {
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
        "GIAC-210 integral table SymPy failures: {:?}",
        failures
    );
    Ok(())
}

#[test]
fn batch1_giac210_test_integrate_line2() -> Result<(), String> {
    let got = run_line("integrate(1/(x^2+1),x)")?;
    verify_sympy("integrate(1/(x^2+1),x)", &got)?;
    Ok(())
}
