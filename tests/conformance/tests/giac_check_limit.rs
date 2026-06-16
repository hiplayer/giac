//! Phase 4: check golden pilot for `check/testlimit` (GIAC-219).

use giac_conformance::{giac_check_dir, load_limit_check_lines, run_line, sympy_verify_lines};

#[test]
fn giac_check_limit_files_exist() {
    assert!(giac_check_dir().join("testlimit").exists());
}

/// Full `testlimit` inventory report (no assertion). Slow — run with `cargo test -- --ignored`.
#[test]
#[ignore = "runs SymPy on all testlimit lines"]
fn giac_check_limit_full_report() -> Result<(), String> {
    let lines = load_limit_check_lines()?;
    let mut eval_ok = 0usize;
    let mut sympy_ok = 0usize;
    for line in &lines {
        let got = match run_line(line) {
            Ok(out) => out,
            Err(_) => continue,
        };
        eval_ok += 1;
        let results = sympy_verify_lines(std::slice::from_ref(line), std::slice::from_ref(&got))?;
        if results[0].ok {
            sympy_ok += 1;
        }
    }
    eprintln!(
        "giac_check_limit full: {eval_ok}/{} eval, {sympy_ok}/{} SymPy",
        lines.len(),
        lines.len()
    );
    Ok(())
}
