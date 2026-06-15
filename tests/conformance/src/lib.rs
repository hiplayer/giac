//! Shared harness for giac-rs conformance tests (giac reference + SymPy third-party).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use giac_core::{exec_stmt, format_expr, Context, StmtResult};
use giac_parse::parse_program;

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
    Path::new(env!("CARGO_MANIFEST_DIR")).join("scripts/phase2_sympy.py")
}

pub fn run_script(path: &Path) -> Result<Vec<String>, String> {
    let input = fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
    let mut ctx = Context::xcas_default();
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
    let mut ctx = Context::xcas_default();
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
    let path = upstream_root().join("bin").join(name);
    let lines = script_lines(&path)?;
    lines.iter().map(|l| triple_check(l)).collect()
}
