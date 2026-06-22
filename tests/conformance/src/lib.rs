//! Shared harness for giac-rs conformance tests (giac reference + SymPy third-party).

mod triple_skip;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::thread;
use std::time::{Duration, Instant};

use giac_core::{exec_stmt, exec_stmts, format_expr, Context, Stmt, StmtResult};
use giac_simplify::assert_equiv;
use giac_ode::xcas_default;
use giac_parse::{parse_compound_line, parse_program};
use serde::Deserialize;

pub use triple_skip::{
    phase2_format_diff, phase2_giac_gap, phase2_sympy_gap, phase3_format_diff,
    phase3_numerical_decomp, phase3_skip, phase3_sympy_gap, phase4_skip, trig_format_diff,
};

/// Subprocess (SymPy `python3`) wall-clock cap; override with `GIAC_CHECK_TIMEOUT_SECS`.
///
/// Rust eval hangs are capped by [cargo-nextest](https://nexte.st/) per-test `slow-timeout`
/// (see `giac-rs/.config/nextest.toml`); this applies only to external processes.
pub const DEFAULT_CHECK_TIMEOUT_SECS: u64 = 10;

/// Timeout for `python3` SymPy verify/equiv subprocesses.
pub fn subprocess_timeout() -> Duration {
    std::env::var("GIAC_CHECK_TIMEOUT_SECS")
        .ok()
        .and_then(|s| s.parse().ok())
        .map(Duration::from_secs)
        .unwrap_or_else(|| Duration::from_secs(DEFAULT_CHECK_TIMEOUT_SECS))
}

fn command_with_timeout(
    mut cmd: Command,
    timeout: Duration,
    label: &str,
) -> Result<Output, String> {
    let mut child = cmd
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("spawn {label}: {e}"))?;
    let deadline = Instant::now() + timeout;
    loop {
        match child.try_wait() {
            Ok(Some(_)) => {
                return child
                    .wait_with_output()
                    .map_err(|e| format!("wait {label}: {e}"));
            }
            Ok(None) => {
                if Instant::now() >= deadline {
                    let _ = child.kill();
                    let _ = child.wait();
                    return Err(format!("timeout ({timeout:?}) on {label}"));
                }
                thread::sleep(Duration::from_millis(50));
            }
            Err(e) => return Err(format!("try_wait {label}: {e}")),
        }
    }
}

struct EvalLineResult {
    output: String,
    ctx: Context,
}

fn eval_line_in_ctx(line: &str, mut ctx: Context) -> Result<EvalLineResult, String> {
    let stmts = if line.contains(',') && (line.starts_with("assume(") || line.contains("purge(")) {
        parse_compound_line(line, &ctx).map_err(|e| format!("parse `{line}`: {e}"))?
    } else {
        parse_program(&format!("{line};"), &ctx)
            .map_err(|e| format!("parse `{line}`: {e}"))?
    };
    if stmts.is_empty() {
        return Err(format!("empty `{line}`"));
    }
    let result = exec_stmts(&stmts, &mut ctx).map_err(|e| format!("eval `{line}`: {e}"))?;
    let output = match result {
        StmtResult::Value(v) | StmtResult::Assign { value: v, .. } => format_expr(v.as_ref()),
        StmtResult::NoValue => String::new(),
    };
    Ok(EvalLineResult { output, ctx })
}

pub fn upstream_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .expect("upstream giac root")
}

/// Upstream source tree for `check/` golden scripts (aligned with CMake `GIAC_VERSION_DIR`).
pub const GIAC_VERSION_DIR: &str = "giac-2.0.0";

/// `giac/<GIAC_VERSION_DIR>/check`
pub fn giac_check_dir() -> PathBuf {
    upstream_root().join("giac").join(GIAC_VERSION_DIR).join("check")
}

pub fn giac_binary() -> PathBuf {
    if let Ok(path) = std::env::var("GIAC_BINARY") {
        return PathBuf::from(path);
    }
    let root = upstream_root();
    let build_20 = root.join("build-2.0/bin/giac");
    if build_20.is_file() {
        return build_20;
    }
    root.join("build/bin/giac")
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

/// Phase 4 bin scripts (solve / calculus / ODE).
pub const PHASE4_SCRIPTS: &[&str] = &[
    "test_solve",
    "test_solve_ext",
    "test_diff",
    "test_integrate",
    "test_integrate_ext",
    "test_integrate_more",
    "test_limit",
    "test_series",
    "test_sturm",
    "test_sturm_ext",
    "test_desolve",
    "test_desolve_ext",
    "test_partfrac_ext",
];

/// Upstream `check/testfactor` + `check/factor.out` (giac_check_factor).
pub fn factor_check_paths() -> (PathBuf, PathBuf) {
    let root = giac_check_dir();
    (root.join("testfactor"), root.join("factor.out"))
}

/// Upstream `check/testintegrate` (GIAC-219).
pub fn integrate_check_path() -> PathBuf {
    giac_check_dir().join("testintegrate")
}

/// Upstream `check/testlimit` (GIAC-219).
pub fn limit_check_path() -> PathBuf {
    giac_check_dir().join("testlimit")
}

/// Kind of line in upstream check golden scripts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CheckLineKind {
    Integrate,
    Limit,
    Series,
    Compound,
    Other,
}

/// Normalize a raw upstream check line (`**` → `^`, trim `;`).
pub fn normalize_check_line(raw: &str) -> String {
    raw.trim()
        .replace("**", "^")
        .trim_end_matches(';')
        .to_string()
}

/// Classify a normalized check line.
pub fn classify_check_line(line: &str) -> CheckLineKind {
    if line.contains(',') && (line.starts_with("assume(") || line.contains("purge(")) {
        return CheckLineKind::Compound;
    }
    if line.starts_with("integrate(") || line.starts_with("int(") {
        return CheckLineKind::Integrate;
    }
    if line.starts_with("limit(") {
        return CheckLineKind::Limit;
    }
    if line.starts_with("series(") || line.starts_with("taylor(") {
        return CheckLineKind::Series;
    }
    CheckLineKind::Other
}

/// Load all executable lines from `check/testintegrate` (skips `cas_setup`).
pub fn load_integrate_check_lines() -> Result<Vec<String>, String> {
    let path = integrate_check_path();
    let inputs = script_lines(&path)?;
    Ok(inputs
        .into_iter()
        .filter(|l| !l.starts_with("cas_setup"))
        .map(|l| normalize_check_line(&l))
        .collect())
}

/// Load all executable lines from `check/testlimit` (skips numeric-only setup lines).
pub fn load_limit_check_lines() -> Result<Vec<String>, String> {
    let path = limit_check_path();
    let inputs = script_lines(&path)?;
    Ok(inputs
        .into_iter()
        .filter(|l| l.starts_with("limit("))
        .map(|l| normalize_check_line(&l))
        .collect())
}

/// Run `risch(f,x)` for the same integrand as `integrate(f,x)`.
pub fn run_risch_line(integrate_line: &str) -> Result<String, String> {
    let inner = integrate_line
        .strip_prefix("integrate(")
        .or_else(|| integrate_line.strip_prefix("int("))
        .and_then(|s| s.strip_suffix(')'))
        .ok_or_else(|| format!("not an integrate line: {integrate_line}"))?;
    run_line(&format!("risch({inner})"))
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
    eval_line_in_ctx(line, xcas_default()).map(|r| r.output)
}

pub fn verify_sympy_with_timeout(
    line: &str,
    output: &str,
    timeout: Duration,
) -> Result<(), String> {
    let script = sympy_script();
    if !script.exists() {
        return Err(format!("sympy script not found at {}", script.display()));
    }
    let label = format!("SymPy verify `{line}`");
    let mut cmd = Command::new("python3");
    cmd.arg(&script).arg("verify").arg(line).arg(output);
    let finished = command_with_timeout(cmd, timeout, &label)?;
    if finished.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&finished.stderr);
        Err(format!(
            "SymPy verify failed for `{line}` output `{output}`: {stderr}"
        ))
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
    verify_sympy_with_timeout(line, output, subprocess_timeout())
}

pub fn sympy_equiv(a: &str, b: &str) -> Result<(), String> {
    if outputs_assert_equiv(a, b)? {
        return Ok(());
    }
    let script = sympy_script();
    let timeout = subprocess_timeout();
    let mut cmd = Command::new("python3");
    cmd.arg(&script).arg("equiv").arg(a).arg(b);
    let finished = command_with_timeout(cmd, timeout, "SymPy equiv")?;
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
    sympy_verify_line_with_timeout(line, output, subprocess_timeout())
}

pub fn sympy_verify_line_with_timeout(
    line: &str,
    output: &str,
    timeout: Duration,
) -> Result<SympyResult, String> {
    Ok(SympyResult {
        line: line.to_string(),
        output: output.to_string(),
        ok: verify_sympy_with_timeout(line, output, timeout).is_ok(),
    })
}

/// SymPy-verify a batch of input lines (e.g. testcas subset).
pub fn sympy_verify_lines(lines: &[String], outputs: &[String]) -> Result<Vec<SympyResult>, String> {
    sympy_verify_lines_with_timeout(lines, outputs, subprocess_timeout())
}

pub fn sympy_verify_lines_with_timeout(
    lines: &[String],
    outputs: &[String],
    timeout: Duration,
) -> Result<Vec<SympyResult>, String> {
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
        .map(|(l, o)| sympy_verify_line_with_timeout(l, o, timeout))
        .collect()
}

#[derive(Debug)]
pub struct SympyResult {
    pub line: String,
    pub output: String,
    pub ok: bool,
}

pub fn load_testcas_lines(n: usize) -> Result<(Vec<String>, Vec<String>), String> {
    let check = giac_check_dir();
    let inputs = script_lines(&check.join("testcas"))?;
    let expected_text = fs::read_to_string(check.join("cas.out.norm"))
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
    for line in lines {
        let result = eval_line_in_ctx(line, ctx.clone())?;
        out.push(result.output);
        ctx = result.ctx;
    }
    Ok(out)
}

/// Assert SymPy verification on every line of a bin script.
pub fn assert_sympy_script(name: &str, allow_gap: fn(&str) -> bool) -> Result<(), String> {
    let results = sympy_verify_script(name)?;
    for r in &results {
        if r.ok || allow_gap(&r.line) {
            continue;
        }
        return Err(format!(
            "SymPy failed on {name}: `{}` -> `{}`",
            r.line, r.output
        ));
    }
    Ok(())
}

/// SymPy-verify each non-skipped line of a bin script (phase 3 style).
pub fn assert_sympy_script_lines(
    name: &str,
    skip: fn(&str) -> bool,
    allow_gap: fn(&str) -> bool,
) -> Result<(), String> {
    let path = upstream_root().join("bin").join(name);
    for line in script_lines(&path)? {
        if skip(&line) {
            continue;
        }
        let got = run_line(&line)?;
        let r = sympy_verify_line(&line, &got)?;
        if !r.ok && !allow_gap(&line) {
            return Err(format!("SymPy failed on {name}: `{line}` -> `{got}`"));
        }
    }
    Ok(())
}

/// SymPy-verify `testcas` lines in `[start, end)`; failures unless `allow_gap`.
pub fn assert_testcas_sympy_range(
    start: usize,
    end: usize,
    allow_gap: fn(&str) -> bool,
    min_verified: usize,
) -> Result<(), String> {
    let (inputs, _) = load_testcas_lines(end)?;
    let inputs: Vec<String> = inputs.into_iter().skip(start).take(end - start).collect();
    let outputs = run_lines(&inputs)?;
    let results = sympy_verify_lines(&inputs, &outputs)?;
    let mut ok = 0usize;
    let mut skipped = 0usize;
    let mut failed = Vec::new();
    for r in &results {
        if r.ok {
            ok += 1;
        } else if allow_gap(&r.line) {
            skipped += 1;
        } else {
            failed.push((r.line.clone(), r.output.clone()));
        }
    }
    assert!(
        failed.is_empty(),
        "testcas [{start},{end}) SymPy failures: {:?}",
        failed.iter().take(5).collect::<Vec<_>>()
    );
    assert!(
        ok + skipped >= min_verified,
        "testcas [{start},{end}): need >= {min_verified} ok+gap, got {ok} ok + {skipped} gaps"
    );
    Ok(())
}

/// SymPy-verify one upstream factor check line by index (skips `cas_setup`).
pub fn assert_factor_line_sympy(index: usize) -> Result<(), String> {
    let (inputs, _) = load_factor_check_lines()?;
    let line = inputs.get(index).ok_or_else(|| {
        format!(
            "factor line index {index} out of range ({} lines)",
            inputs.len()
        )
    })?;
    let timeout = subprocess_timeout();
    let outputs = run_lines(std::slice::from_ref(line))?;
    let results =
        sympy_verify_lines_with_timeout(std::slice::from_ref(line), &outputs, timeout)?;
    assert!(
        results[0].ok,
        "SymPy failed on factor check `{}` -> `{}`",
        results[0].line, results[0].output
    );
    Ok(())
}

/// Triple-check one line from a bin script by line index.
pub fn triple_check_script_line(
    name: &str,
    line_index: usize,
    skip: fn(&str) -> bool,
) -> Result<TripleResult, String> {
    let path = upstream_root().join("bin").join(name);
    let lines = script_lines(&path)?;
    let line = lines.get(line_index).ok_or_else(|| {
        format!(
            "script {name} line index {line_index} out of range ({} lines)",
            lines.len()
        )
    })?;
    if skip(line) {
        return Err(format!("script {name} line {line_index} skipped: `{line}`"));
    }
    triple_check(line)
}

#[derive(Debug, Deserialize)]
pub struct CheckIntegrateTable {
    pub entries: Vec<CheckIntegrateEntry>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CheckIntegrateEntry {
    pub id: String,
    pub line: String,
    pub kind: String,
    pub enabled: bool,
}

pub fn load_check_integrate_table() -> Result<CheckIntegrateTable, String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/check_integrate_table.json");
    let text = fs::read_to_string(&path).map_err(|e| format!("read {}: {e}", path.display()))?;
    serde_json::from_str(&text).map_err(|e| format!("json: {e}"))
}

pub fn check_integrate_entry(id: &str) -> Result<CheckIntegrateEntry, String> {
    load_check_integrate_table()?
        .entries
        .into_iter()
        .find(|e| e.id == id)
        .ok_or_else(|| format!("check_integrate entry {id} not found"))
}

/// SymPy gate for one enabled non-integrate check_integrate row (limit/series).
pub fn assert_check_integrate_non_integrate_sympy(id: &str) -> Result<(), String> {
    let entry = check_integrate_entry(id)?;
    assert!(entry.enabled, "{id} is not enabled");
    assert_ne!(
        entry.kind, "integrate",
        "{id} is integrate; use assert_check_integrate_risch"
    );
    let timeout = subprocess_timeout();
    let got = run_line(&entry.line)?;
    let results = sympy_verify_lines_with_timeout(
        std::slice::from_ref(&entry.line),
        std::slice::from_ref(&got),
        timeout,
    )?;
    assert!(
        results[0].ok,
        "SymPy failed on {id} `{}` -> `{}`",
        entry.line, got
    );
    Ok(())
}

/// `risch(f,x)` agrees with `integrate(f,x)` for one enabled integrate row.
pub fn assert_check_integrate_risch(id: &str) -> Result<(), String> {
    let entry = check_integrate_entry(id)?;
    assert!(entry.enabled, "{id} is not enabled");
    assert_eq!(entry.kind, "integrate", "{id} is not integrate");
    let int_out = run_line(&entry.line)?;
    let risch_out = run_risch_line(&entry.line)?;
    sympy_equiv(&int_out, &risch_out).map_err(|e| format!("{id} `{}`: {e}", entry.line))
}

#[cfg(test)]
mod timeout_tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn command_with_timeout_kills_slow_subprocess() {
        let mut cmd = Command::new("sleep");
        cmd.arg("5");
        let err = command_with_timeout(cmd, Duration::from_millis(200), "sleep")
            .unwrap_err();
        assert!(err.contains("timeout"), "{err}");
    }

    #[test]
    fn subprocess_timeout_default_is_ten_seconds() {
        assert_eq!(subprocess_timeout(), Duration::from_secs(10));
    }
}
