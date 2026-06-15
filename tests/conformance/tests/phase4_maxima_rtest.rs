//! Maxima rtest-derived conformance draft (§4.x external corpus).

use std::fs;
use std::path::{Path, PathBuf};

use giac_conformance::{run_line, verify_sympy};
use giac_ode::xcas_default;
use giac_parse::parse_program;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct MaximaRtestFixture {
    version: u32,
    entries: Vec<MaximaRtestEntry>,
}

#[derive(Debug, Deserialize)]
struct MaximaRtestEntry {
    id: String,
    line: String,
    api: String,
    #[allow(dead_code)]
    verify: String,
    enabled: bool,
    #[serde(default)]
    #[allow(dead_code)]
    maxima_input: Option<String>,
    #[serde(default)]
    #[allow(dead_code)]
    maxima_expected: Option<String>,
    #[serde(default)]
    #[allow(dead_code)]
    source_file: Option<String>,
    #[serde(default)]
    #[allow(dead_code)]
    source_line: Option<u32>,
}

fn fixture_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/phase4_maxima_rtest.json")
}

fn load_fixture() -> MaximaRtestFixture {
    let text = fs::read_to_string(fixture_path()).expect("phase4_maxima_rtest.json");
    serde_json::from_str(&text).expect("valid maxima rtest JSON")
}

#[test]
fn phase4_maxima_rtest_json_valid() {
    let fixture = load_fixture();
    assert_eq!(fixture.version, 1);
    assert!(!fixture.entries.is_empty(), "maxima rtest fixture is empty");
    let mut ids = std::collections::HashSet::new();
    for e in &fixture.entries {
        assert!(ids.insert(e.id.clone()), "duplicate id {}", e.id);
        assert!(
            matches!(e.api.as_str(), "integrate" | "limit" | "solve" | "diff"),
            "unknown api {} on {}",
            e.api,
            e.id
        );
    }
}

#[test]
fn phase4_maxima_rtest_all_parse() -> Result<(), String> {
    let fixture = load_fixture();
    let ctx = xcas_default();
    let mut failures = Vec::new();
    for entry in &fixture.entries {
        let src = format!("{};", entry.line);
        if let Err(e) = parse_program(&src, &ctx) {
            failures.push((entry.id.clone(), entry.line.clone(), e.to_string()));
        }
    }
    assert!(
        failures.is_empty(),
        "maxima rtest parse failures: {:?}",
        failures
    );
    Ok(())
}

#[test]
fn phase4_maxima_rtest_enabled_sympy() -> Result<(), String> {
    let fixture = load_fixture();
    let enabled: Vec<_> = fixture.entries.iter().filter(|e| e.enabled).collect();
    if enabled.is_empty() {
        return Ok(());
    }

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
        "maxima rtest SymPy failures: {:?}",
        failures
    );
    Ok(())
}
