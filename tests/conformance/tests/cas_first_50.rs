//! Golden regression: `check/testcas` first 50 lines vs `cas.out.norm`.

use std::fs;
use std::path::{Path, PathBuf};

use giac_core::{exec_stmt, format_expr, Context, StmtResult};
use giac_parse::parse_program;

fn upstream_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .expect("upstream giac root")
}

fn load_lines(path: &Path) -> Vec<String> {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    text.lines()
        .map(|l| l.trim_end_matches(',').trim().to_string())
        .filter(|l| !l.is_empty())
        .collect()
}

fn run_script_lines(lines: &[&str]) -> Vec<String> {
    let input = lines
        .iter()
        .map(|l| {
            if l.ends_with(';') {
                l.to_string()
            } else {
                format!("{l};")
            }
        })
        .collect::<String>();
    let mut ctx = Context::xcas_default();
    let stmts = parse_program(&input, &ctx).unwrap_or_else(|e| panic!("parse: {e}"));
    let mut out = Vec::new();
    for stmt in &stmts {
        match exec_stmt(stmt, &mut ctx).unwrap() {
            StmtResult::Value(v) | StmtResult::Assign { value: v, .. } => {
                out.push(format_expr(v.as_ref()));
            }
            StmtResult::NoValue => {}
        }
    }
    out
}

#[test]
fn cas_tst_first_25_batch() {
    let root = upstream_root();
    let inputs = load_lines(&root.join("giac/giac-1.5.0/check/testcas"));
    let expected = load_lines(&root.join("giac/giac-1.5.0/check/cas.out.norm"));
    let n = 25;
    let input_slice: Vec<&str> = inputs.iter().take(n).map(String::as_str).collect();
    let got = run_script_lines(&input_slice);
    assert_eq!(got.len(), n);
    let passed = (0..n).filter(|i| got.get(*i) == expected.get(*i)).count();
    eprintln!("passed {passed}/{n}");
    assert!(passed >= 23);
}

#[test]
fn cas_tst_first_30_batch() {
    let root = upstream_root();
    let inputs = load_lines(&root.join("giac/giac-1.5.0/check/testcas"));
    let expected = load_lines(&root.join("giac/giac-1.5.0/check/cas.out.norm"));
    let n = 30;
    let input_slice: Vec<&str> = inputs.iter().take(n).map(String::as_str).collect();
    let got = run_script_lines(&input_slice);
    assert_eq!(got.len(), n);
    let passed = (0..n).filter(|i| got.get(*i) == expected.get(*i)).count();
    eprintln!("passed {passed}/{n}");
    assert!(passed >= 25);
}

#[test]
fn cas_tst_first_50_lines() {
    let root = upstream_root();
    let inputs = load_lines(&root.join("giac/giac-1.5.0/check/testcas"));
    let expected = load_lines(&root.join("giac/giac-1.5.0/check/cas.out.norm"));
    let n = 50.min(inputs.len()).min(expected.len());
    let input_slice: Vec<&str> = inputs.iter().take(n).map(String::as_str).collect();
    let got = run_script_lines(&input_slice);
    assert_eq!(got.len(), n, "expected {n} results, got {}", got.len());
    let passed = (0..n).filter(|i| got.get(*i) == expected.get(*i)).count();
    eprintln!("passed {passed}/{n}");
    assert!(passed >= 25, "cas.tst first {n}: {passed}/{n} golden matches (need >= 25)");
}

#[test]
fn cas_tst_first_20_lines_progress() {
    let root = upstream_root();
    let inputs = load_lines(&root.join("giac/giac-1.5.0/check/testcas"));
    let expected = load_lines(&root.join("giac/giac-1.5.0/check/cas.out.norm"));
    let n = 20.min(inputs.len()).min(expected.len());
    let input_slice: Vec<&str> = inputs.iter().take(n).map(String::as_str).collect();
    let got = run_script_lines(&input_slice);
    let passed = (0..n).filter(|i| got.get(*i) == expected.get(*i)).count();
    assert!(
        passed >= 15,
        "expected at least 15/20 cas lines passing, got {passed}/20\nfailures: {:?}",
        (0..n)
            .filter(|i| got.get(*i) != expected.get(*i))
            .map(|i| (i + 1, &inputs[i], got.get(i), &expected[i]))
            .collect::<Vec<_>>()
    );
}
