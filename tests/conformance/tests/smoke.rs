//! Golden regression harness for giac-rs.

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

fn run_script(path: &Path) -> Vec<String> {
    let input = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut ctx = Context::xcas_default();
    let stmts = parse_program(&input, &ctx).unwrap_or_else(|e| {
        panic!("parse {}: {e}", path.display())
    });
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
fn test_cas_basic_matches_giac() {
    let script = upstream_root().join("bin/test_cas_basic");
    let got = run_script(&script);
    let want = vec![
        "sqrt(5)".to_string(),
        "15".to_string(),
        "-3-4*i".to_string(),
    ];
    assert_eq!(got, want, "bin/test_cas_basic");
}

#[test]
fn harness_phase0_bin_scripts_parseable() {
    let scripts = ["test_cas_basic"];
    for name in scripts {
        let path = upstream_root().join("bin").join(name);
        let input = fs::read_to_string(&path).unwrap();
        let ctx = Context::xcas_default();
        assert!(
            parse_program(&input, &ctx).is_ok(),
            "parse failed for {name}"
        );
    }
}
