//! Shared harness for giac-rs conformance tests (giac reference + SymPy third-party).

mod triple_skip;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use giac_core::{assert_equiv, exec_stmt, format_expr, Stmt, StmtResult};
use giac_linalg::xcas_default;
use giac_parse::parse_program;

pub use triple_skip::{
    phase2_format_diff, phase2_giac_gap, phase2_sympy_gap, phase3_format_diff,
    phase3_numerical_decomp, phase3_skip, phase3_sympy_gap,
};

pub fn upstream_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .expect("upstream giac root")
}

pub fn giac_binary() -> PathBuf {
    upstream_root().join("build/bin/giac")
}

pub fn sympy_script() -> PathBuf {
    let primary = Path::new(env!("CARGO_MANIFEST_DIR")).join("scripts/sympy_verify.py");
    if primary.exists() {
        return primary;
    }
    Path::new(env!("CARGO_MANIFEST_DIR")).join("scripts/phase2_sympy.py")
}

/// Phase 1 bin scripts (basic CAS).
pub const PHASE1_SCRIPTS: &[&str] = &["test_cas_basic"];

/// Phase 2 bin scripts (polynomial ring).
pub const PHASE2_SCRIPTS: &[&str] = &[
    "test_poly",
    "test_poly_ext",
    "test_modular",
    "test_factor",
    "test_groebner",
];

/// Phase 3 bin scripts (linear algebra).
pub const PHASE3_SCRIPTS: &[&str] = &[
    "test_linalg",
    "test_linalg_ext",
    "test_linalg_decomp",
    "test_gauss_ext",
];

/// Upstream `check/testfactor` + `check/factor.out` (giac_check_factor).
pub fn factor_check_paths() -> (PathBuf, PathBuf) {
    let root = upstream_root().join("giac/giac-1.5.0/check");
    (root.join("testfactor"), root.join("factor.out"))
}

/// Load factor regression inputs and golden outputs (skips `cas_setup` line).
pub fn load_factor_check_lines() -> Result<(Vec<String>, Vec<String>), String> {
    let (input_path, golden_path) = factor_check_paths();
    let inputs = script_lines(&input_path)?;
    let golden_text =
        fs::read_to_string(&golden_path).map_err(|e| format!("read {}: {e}", golden_path.display()))?;
    let golden: Vec<String> = golden_text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(String::from)
        .collect();
    let inputs: Vec<String> = inputs.into_iter().skip(1).collect();
    let golden: Vec<String> = golden.into_iter().skip(1).collect();
    if inputs.len() != golden.len() {
        return Err(format!(
            "testfactor {} lines vs factor.out {} lines",
            inputs.len(),
            golden.len()
        ));
    }
    Ok((inputs, golden))
}

pub fn run_factor_check() -> Result<Vec<String>, String> {
    let (inputs, _) = load_factor_check_lines()?;
    run_lines(&inputs)
}

pub const ALL_SYMPTY_SCRIPTS: &[&str] = &[
    "test_cas_basic",
    "test_poly",
    "test_poly_ext",
    "test_modular",
    "test_factor",
    "test_groebner",
];

pub fn run_script(path: &Path) -> Result<Vec<String>, String> {
    let input = fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
    let mut ctx = xcas_default();
    let stmts = parse_program(&input, &ctx).map_err(|e| format!("parse {}: {e}", path.display()))?;
    let mut out = Vec::new();
    for (idx, stmt) in stmts.iter().enumerate() {
        match exec_stmt(stmt, &mut ctx)
            .map_err(|e| format!("eval stmt {idx} in {}: {e}", path.display()))?
        {
            StmtResult::Value(v) | StmtResult::Assign { value: v, .. } => {
                out.push(format_expr(v.as_ref()));
            }
            StmtResult::NoValue => {}
        }
    }
    Ok(out)
}

pub fn run_line(line: &str) -> Result<String, String> {
    let mut ctx = xcas_default();
    let stmts = parse_program(&format!("{line};"), &ctx)
        .map_err(|e| format!("parse {line}: {e}"))?;
    let stmt = stmts.first().ok_or_else(|| format!("empty {line}"))?;
    match exec_stmt(stmt, &mut ctx).map_err(|e| format!("eval {line}: {e}"))? {
        StmtResult::Value(v) | StmtResult::Assign { value: v, .. } => {
            Ok(format_expr(v.as_ref()))
        }
        StmtResult::NoValue => Ok(String::new()),
    }
}

pub fn script_lines(path: &Path) -> Result<Vec<String>, String> {
    let input = fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
    Ok(input
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with("//"))
        .map(|l| l.trim_end_matches(';').to_string())
        .collect())
}

pub fn run_giac(line: &str) -> Result<String, String> {
    let giac = giac_binary();
    if !giac.exists() {
        return Err(format!("giac binary not found at {}", giac.display()));
    }
    let mut child = Command::new(&giac)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("spawn giac: {e}"))?;
    use std::io::Write;
    if let Some(mut stdin) = child.stdin.take() {
        writeln!(stdin, "{line};").map_err(|e| format!("write giac stdin: {e}"))?;
    }
    let finished = child.wait_with_output().map_err(|e| format!("wait giac: {e}"))?;
    parse_giac_output(&finished)
}

fn parse_giac_output(output: &Output) -> Result<String, String> {
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("giac failed: {stderr}"));
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    let line = stdout
        .lines()
        .map(str::trim)
        .find(|l| !l.is_empty())
        .unwrap_or("")
        .to_string();
    Ok(line)
}

pub fn verify_sympy(line: &str, output: &str) -> Result<(), String> {
    let script = sympy_script();
    if !script.exists() {
        return Err(format!("sympy script not found at {}", script.display()));
    }
    let finished = Command::new("python3")
        .arg(&script)
        .arg("verify")
        .arg(line)
        .arg(output)
        .output()
        .map_err(|e| format!("spawn python3: {e}"))?;
    if finished.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&finished.stderr);
        Err(format!("SymPy verify failed for `{line}` output `{output}`: {stderr}"))
    }
}

pub fn sympy_equiv(a: &str, b: &str) -> Result<(), String> {
    if outputs_assert_equiv(a, b)? {
        return Ok(());
    }
    let script = sympy_script();
    let finished = Command::new("python3")
        .arg(&script)
        .arg("equiv")
        .arg(a)
        .arg(b)
        .output()
        .map_err(|e| format!("spawn python3: {e}"))?;
    if finished.status.success() {
        Ok(())
    } else {
        Err(format!("not equivalent: `{a}` vs `{b}`"))
    }
}

/// Parse a giac-style output string as a single expression.
pub fn parse_output_expr(s: &str) -> Result<std::sync::Arc<giac_core::Expr>, String> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return Err("empty output".to_string());
    }
    let input = if trimmed.ends_with(';') {
        trimmed.to_string()
    } else {
        format!("{trimmed};")
    };
    let ctx = xcas_default();
    let stmts = parse_program(&input, &ctx).map_err(|e| format!("parse `{trimmed}`: {e}"))?;
    match stmts.first() {
        Some(Stmt::ExprStmt(e)) => Ok(std::sync::Arc::clone(e)),
        Some(Stmt::Assign(_, e)) => Ok(std::sync::Arc::clone(e)),
        _ => Err(format!("expected expression output, got stmt for `{trimmed}`")),
    }
}

/// Mathematical equivalence via `normal(a-b)==0` (conformance-testing.md §3).
pub fn outputs_assert_equiv(a: &str, b: &str) -> Result<bool, String> {
    if a == b {
        return Ok(true);
    }
    let ctx = xcas_default();
    let ea = parse_output_expr(a)?;
    let eb = parse_output_expr(b)?;
    assert_equiv(ea.as_ref(), eb.as_ref(), &ctx).map_err(|e| e.to_string())
}

/// Triple-check one giac line: giac-rs vs SymPy, giac reference, cross-equivalence.
pub fn triple_check(line: &str) -> Result<TripleResult, String> {
    let rs = run_line(line)?;
    let sympy_rs_ok = verify_sympy(line, &rs).is_ok();

    let giac_out = run_giac(line)?;
    let sympy_giac_ok = verify_sympy(line, &giac_out).is_ok();

    let rs_giac_equiv = sympy_equiv(&rs, &giac_out).is_ok();

    Ok(TripleResult {
        line: line.to_string(),
        giac_rs: rs,
        giac: giac_out,
        sympy_rs_ok,
        sympy_giac_ok,
        rs_giac_equiv,
    })
}

#[derive(Debug)]
pub struct TripleResult {
    pub line: String,
    pub giac_rs: String,
    pub giac: String,
    pub sympy_rs_ok: bool,
    pub sympy_giac_ok: bool,
    pub rs_giac_equiv: bool,
}

pub fn triple_check_script(name: &str) -> Result<Vec<TripleResult>, String> {
    triple_check_script_filtered(name, |_| false)
}

/// Triple-check script lines, skipping those for which `skip(line)` is true.
pub fn triple_check_script_filtered(
    name: &str,
    skip: impl Fn(&str) -> bool,
) -> Result<Vec<TripleResult>, String> {
    let path = upstream_root().join("bin").join(name);
    let lines = script_lines(&path)?;
    lines
        .iter()
        .filter(|l| !skip(l))
        .map(|l| triple_check(l))
        .collect()
}

/// Require SymPy verification on giac-rs for every triple result (with optional gaps).
pub fn triple_assert_sympy_rs(
    results: &[TripleResult],
    known_gap: impl Fn(&str) -> bool,
) -> Result<(), String> {
    for r in results {
        if !r.sympy_rs_ok && !known_gap(&r.line) {
            return Err(format!(
                "giac-rs failed SymPy on {}: {}",
                r.line, r.giac_rs
            ));
        }
    }
    Ok(())
}

/// Log giac-rs vs giac format differences that are not SymPy failures.
pub fn triple_note_format_diffs(results: &[TripleResult], format_diff: impl Fn(&str) -> bool) {
    for r in results {
        if !r.rs_giac_equiv && !format_diff(&r.line) {
            eprintln!(
                "note: {} giac-rs={} giac={} (sympy_rs={} sympy_giac={})",
                r.line, r.giac_rs, r.giac, r.sympy_rs_ok, r.sympy_giac_ok
            );
        }
    }
}

/// Outcome of comparing giac-rs output to golden (conformance-testing.md §3.5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckOutcome {
    LiteralMatch,
    EquivMatch,
    Fail,
}

/// Compare outputs: literal first, then `assert_equiv`.
pub fn check_output_equiv(got: &str, expected: &str) -> Result<CheckOutcome, String> {
    if got == expected {
        return Ok(CheckOutcome::LiteralMatch);
    }
    if outputs_assert_equiv(got, expected)? {
        return Ok(CheckOutcome::EquivMatch);
    }
    Ok(CheckOutcome::Fail)
}

/// Run giac-rs on a script and SymPy-verify every output line.
pub fn sympy_verify_script(name: &str) -> Result<Vec<SympyResult>, String> {
    let path = upstream_root().join("bin").join(name);
    let lines = script_lines(&path)?;
    let outputs = run_script(&path)?;
    if lines.len() != outputs.len() {
        return Err(format!(
            "script {name}: {} lines vs {} outputs",
            lines.len(),
            outputs.len()
        ));
    }
    lines
        .iter()
        .zip(outputs.iter())
        .map(|(line, out)| sympy_verify_line(line, out))
        .collect()
}

/// SymPy-verify giac-rs output for a single input line.
pub fn sympy_verify_line(line: &str, output: &str) -> Result<SympyResult, String> {
    Ok(SympyResult {
        line: line.to_string(),
        output: output.to_string(),
        ok: verify_sympy(line, output).is_ok(),
    })
}

/// SymPy-verify a batch of input lines (e.g. testcas subset).
pub fn sympy_verify_lines(lines: &[String], outputs: &[String]) -> Result<Vec<SympyResult>, String> {
    if lines.len() != outputs.len() {
        return Err(format!(
            "sympy_verify_lines: {} inputs vs {} outputs",
            lines.len(),
            outputs.len()
        ));
    }
    lines
        .iter()
        .zip(outputs.iter())
        .map(|(l, o)| sympy_verify_line(l, o))
        .collect()
}

#[derive(Debug)]
pub struct SympyResult {
    pub line: String,
    pub output: String,
    pub ok: bool,
}

pub fn load_testcas_lines(n: usize) -> Result<(Vec<String>, Vec<String>), String> {
    let root = upstream_root();
    let inputs = script_lines(&root.join("giac/giac-1.5.0/check/testcas"))?;
    let expected_text = fs::read_to_string(root.join("giac/giac-1.5.0/check/cas.out.norm"))
        .map_err(|e| format!("read cas.out.norm: {e}"))?;
    let expected: Vec<String> = expected_text
        .lines()
        .map(|l| l.trim_end_matches(',').trim())
        .filter(|l| !l.is_empty())
        .map(String::from)
        .collect();
    let n = n.min(inputs.len()).min(expected.len());
    Ok((
        inputs.into_iter().take(n).collect(),
        expected.into_iter().take(n).collect(),
    ))
}

pub fn run_lines(lines: &[String]) -> Result<Vec<String>, String> {
    let mut ctx = xcas_default();
    let mut out = Vec::new();
    for (idx, line) in lines.iter().enumerate() {
        let stmts = parse_program(&format!("{line};"), &ctx)
            .map_err(|e| format!("parse line {idx} `{line}`: {e}"))?;
        let stmt = stmts.first().ok_or_else(|| format!("empty line {idx}"))?;
        match exec_stmt(stmt, &mut ctx).map_err(|e| format!("eval line {idx} `{line}`: {e}"))? {
            StmtResult::Value(v) | StmtResult::Assign { value: v, .. } => {
                out.push(format_expr(v.as_ref()));
            }
            StmtResult::NoValue => out.push(String::new()),
        }
    }
    Ok(out)
}
