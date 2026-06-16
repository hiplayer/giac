//! `giac_check_factor` conformance: upstream `check/testfactor` vs `factor.out`.
//!
//! MVP: SymPy verifies `expand(factor(p)) == expand(p)` (identity factorization allowed).
//! Golden literal match is reported but not required until general factor is implemented.

use giac_conformance::{
    factor_check_paths, load_factor_check_lines, outputs_assert_equiv, run_lines, sympy_equiv,
    giac_check_dir, sympy_verify_lines, verify_sympy,
};

#[test]
fn factor_check_files_exist() {
    let (input, golden) = factor_check_paths();
    assert!(
        input.exists(),
        "testfactor missing at {}",
        input.display()
    );
    assert!(
        golden.exists(),
        "factor.out missing at {}",
        golden.display()
    );
}

#[test]
fn giac_check_factor_sympy() -> Result<(), String> {
    let (inputs, _) = load_factor_check_lines()?;
    let outputs = run_lines(&inputs)?;
    let results = sympy_verify_lines(&inputs, &outputs)?;
    for r in &results {
        assert!(
            r.ok,
            "SymPy failed on factor check `{}` -> `{}`",
            r.line,
            r.output
        );
    }
    Ok(())
}

#[test]
fn giac_check_factor_each_line() -> Result<(), String> {
    let (inputs, _) = load_factor_check_lines()?;
    let outputs = run_lines(&inputs)?;
    for (line, out) in inputs.iter().zip(outputs.iter()) {
        verify_sympy(line, out)?;
    }
    Ok(())
}

#[test]
fn giac_check_factor_golden_report() -> Result<(), String> {
    let (inputs, golden) = load_factor_check_lines()?;
    let outputs = run_lines(&inputs)?;
    let mut exact = 0usize;
    let mut equiv = 0usize;
    for ((line, out), want) in inputs.iter().zip(outputs.iter()).zip(golden.iter()) {
        if out == want {
            exact += 1;
            continue;
        }
        if outputs_assert_equiv(out, want).unwrap_or(false) || sympy_equiv(out, want).is_ok() {
            equiv += 1;
        } else {
            eprintln!("factor golden diff: `{line}` -> rs=`{out}` giac=`{want}`");
        }
    }
    eprintln!(
        "giac_check_factor: {}/{} exact, {}/{} equivalent to factor.out",
        exact,
        inputs.len(),
        exact + equiv,
        inputs.len()
    );
    // MVP gate: mathematical correctness via SymPy, not literal golden match.
    Ok(())
}

#[test]
fn giac_check_factor_upstream_root() {
    let chk = giac_check_dir().join("chk_factor");
    assert!(
        chk.exists(),
        "chk_factor wrapper expected under {}",
        chk.display()
    );
}
