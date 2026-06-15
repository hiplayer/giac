//! Phase 4 triple-check: giac-rs + upstream giac + SymPy.

use giac_conformance::{
    phase4_skip, triple_assert_sympy_rs, triple_check_script_filtered, triple_note_format_diffs,
};

#[test]
fn test_solve_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_solve", phase4_skip)?;
    triple_assert_sympy_rs(&results, |_| false)?;
    Ok(())
}

#[test]
fn test_diff_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_diff", phase4_skip)?;
    triple_assert_sympy_rs(&results, |_| false)?;
    Ok(())
}

#[test]
fn test_integrate_table_enabled_triple() -> Result<(), String> {
    let results = triple_check_script_filtered("test_integrate", phase4_skip)?;
    triple_note_format_diffs(&results, |_| true);
    triple_assert_sympy_rs(&results, phase4_skip)?;
    Ok(())
}
