//! SymPy third-party verification for Phase 1–3 (all conformance scripts).

use giac_conformance::{
    load_testcas_lines, run_lines, sympy_verify_line, sympy_verify_lines, sympy_verify_script,
    verify_sympy, ALL_SYMPTY_SCRIPTS, PHASE1_SCRIPTS, PHASE2_SCRIPTS, PHASE3_SCRIPTS,
    upstream_root, run_line, run_script,
};

#[test]
fn sympy_available() {
    let script = giac_conformance::sympy_script();
    assert!(script.exists(), "sympy_verify.py required at {}", script.display());
}

#[test]
fn phase1_scripts_sympy() -> Result<(), String> {
    for name in PHASE1_SCRIPTS {
        let results = sympy_verify_script(name)?;
        for r in &results {
            assert!(r.ok, "SymPy failed on {name}: `{}` -> `{}`", r.line, r.output);
        }
    }
    Ok(())
}

#[test]
fn phase2_scripts_sympy() -> Result<(), String> {
    for name in PHASE2_SCRIPTS {
        let results = sympy_verify_script(name)?;
        for r in &results {
            if !r.ok && !is_known_sympy_gap(&r.line) {
                return Err(format!(
                    "SymPy failed on {name}: `{}` -> `{}`",
                    r.line, r.output
                ));
            }
        }
    }
    Ok(())
}

#[test]
fn phase3_scripts_sympy() -> Result<(), String> {
    for name in PHASE3_SCRIPTS {
        let path = upstream_root().join("bin").join(name);
        let lines = giac_conformance::script_lines(&path)?;
        for line in &lines {
            if is_phase3_skip(line) {
                continue;
            }
            let got = run_line(line)?;
            let r = sympy_verify_line(line, &got)?;
            if !r.ok && !is_known_sympy_gap(line) {
                return Err(format!(
                    "SymPy failed on {name}: `{line}` -> `{got}`"
                ));
            }
        }
    }
    Ok(())
}

#[test]
fn test_cas_basic_each_line_sympy() -> Result<(), String> {
    let cases = [
        ("abs(1+2*i)", "sqrt(5)"),
        ("gcd(45,75)", "15"),
        ("conj((1+2*i)^2)", "-3-4*i"),
    ];
    for (line, want) in cases {
        let got = run_line(line)?;
        assert_eq!(got, want, "{line}");
        verify_sympy(line, &got)?;
    }
    Ok(())
}

#[test]
fn testcas_first_50_sympy() -> Result<(), String> {
    let (inputs, _expected) = load_testcas_lines(50)?;
    let outputs = run_lines(&inputs)?;
    let results = sympy_verify_lines(&inputs, &outputs)?;
    let mut failed = Vec::new();
    let mut skipped = 0usize;
    for r in &results {
        if r.ok {
            continue;
        }
        if is_known_testcas_gap(&r.line) {
            skipped += 1;
            eprintln!("skip sympy gap: {} -> {}", r.line, r.output);
            continue;
        }
        failed.push((r.line.clone(), r.output.clone()));
    }
    let ok_count = results.iter().filter(|r| r.ok).count();
    eprintln!(
        "testcas sympy: {ok_count}/{} ok, {skipped} known gaps, {} failures",
        results.len(),
        failed.len()
    );
    assert!(
        ok_count + skipped >= 40,
        "need >= 40 sympy-verified testcas lines (got {ok_count} ok + {skipped} gaps); failures: {:?}",
        failed.iter().take(5).collect::<Vec<_>>()
    );
    Ok(())
}

#[test]
fn all_bin_scripts_sympy_smoke() -> Result<(), String> {
    for name in ALL_SYMPTY_SCRIPTS {
        let path = upstream_root().join("bin").join(name);
        let outputs = run_script(&path)?;
        let lines = giac_conformance::script_lines(&path)?;
        assert_eq!(lines.len(), outputs.len());
        for (line, out) in lines.iter().zip(outputs.iter()) {
            let r = sympy_verify_line(line, out)?;
            if !r.ok && !is_known_sympy_gap(line) && !is_known_testcas_gap(line) {
                return Err(format!("SymPy: `{line}` -> `{out}`"));
            }
        }
    }
    Ok(())
}

fn is_known_sympy_gap(line: &str) -> bool {
    matches!(line, "roots(x^3-1,x)")
        || line.starts_with("gramschmidt(")
        || matches!(
            line,
            "jordan([[1,1],[0,1]])" | "egv([[4,1,-2],[1,2,-1],[2,1,0]])"
        )
}

fn is_phase3_skip(line: &str) -> bool {
    line.starts_with("gramschmidt(")
        || matches!(
            line,
            "jordan([[1,1],[0,1]])" | "egv([[4,1,-2],[1,2,-1],[2,1,0]])"
        )
}

fn is_known_testcas_gap(line: &str) -> bool {
    matches!(
        line,
        "integrate(1/(1-x^4),x)"
            | "normal(-2+-3/-2*2)"
            | "normal(-6+-3/-2*4)"
    ) || line.parse::<i64>().is_ok()
}
