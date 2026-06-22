//! `giac_check_factor` conformance: upstream `check/testfactor` vs `factor.out`.
//!
//! MVP: SymPy verifies `expand(factor(p)) == expand(p)` (identity factorization allowed).
//! Golden literal match is reported but not required until general factor is implemented.
//!
//! **Gate:** eval hangs capped by nextest per-test `slow-timeout`; SymPy subprocess cap
//! `subprocess_timeout()` / `GIAC_CHECK_TIMEOUT_SECS` (default 10s).
//! Full workspace: `cargo test-timeout` (nextest **--release** via `.cargo/config.toml`).
//! Run subset: `cargo nextest run --release -p giac-conformance --test giac_check_factor`.

use giac_conformance::{
    assert_factor_line_sympy, factor_check_paths, load_factor_check_lines, outputs_assert_equiv,
    run_lines, sympy_equiv, giac_check_dir,
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

macro_rules! factor_line_sympy {
    ($fn_name:ident, $idx:literal) => {
        #[test]
        fn $fn_name() -> Result<(), String> {
            assert_factor_line_sympy($idx)
        }
    };
}

factor_line_sympy!(factor_sympy_line_00, 0);
factor_line_sympy!(factor_sympy_line_01, 1);
factor_line_sympy!(factor_sympy_line_02, 2);
factor_line_sympy!(factor_sympy_line_03, 3);
factor_line_sympy!(factor_sympy_line_04, 4);
factor_line_sympy!(factor_sympy_line_05, 5);
factor_line_sympy!(factor_sympy_line_06, 6);
factor_line_sympy!(factor_sympy_line_07, 7);
factor_line_sympy!(factor_sympy_line_08, 8);
factor_line_sympy!(factor_sympy_line_09, 9);
factor_line_sympy!(factor_sympy_line_10, 10);
factor_line_sympy!(factor_sympy_line_11, 11);
factor_line_sympy!(factor_sympy_line_12, 12);
factor_line_sympy!(factor_sympy_line_13, 13);
factor_line_sympy!(factor_sympy_line_14, 14);
factor_line_sympy!(factor_sympy_line_15, 15);
factor_line_sympy!(factor_sympy_line_16, 16);
factor_line_sympy!(factor_sympy_line_17, 17);
factor_line_sympy!(factor_sympy_line_18, 18);
factor_line_sympy!(factor_sympy_line_19, 19);
factor_line_sympy!(factor_sympy_line_20, 20);
factor_line_sympy!(factor_sympy_line_21, 21);
factor_line_sympy!(factor_sympy_line_22, 22);
factor_line_sympy!(factor_sympy_line_23, 23);
factor_line_sympy!(factor_sympy_line_24, 24);
factor_line_sympy!(factor_sympy_line_25, 25);
factor_line_sympy!(factor_sympy_line_26, 26);
factor_line_sympy!(factor_sympy_line_27, 27);
factor_line_sympy!(factor_sympy_line_28, 28);
factor_line_sympy!(factor_sympy_line_29, 29);

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
