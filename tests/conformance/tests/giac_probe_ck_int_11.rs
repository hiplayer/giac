//! Manual SymPy probe — not run in CI.

use giac_conformance::{run_line, sympy_verify_lines};

#[test]
#[ignore = "manual SymPy integrate probe"]
fn probe_ck_int_11_sin2x_over_cos2x() {
    let line = "integrate((sin(2*x)+1)/(cos(2*x)),x)";
    let got = run_line(line).expect("integrate");
    let ok = sympy_verify_lines(&[line.to_string()], &[got.clone()])
        .map(|r| r[0].ok)
        .unwrap_or(false);
    assert!(ok, "got {got}");
}

#[test]
fn probe_ck_int_11_eval() {
    let line = "integrate((sin(2*x)+1)/(cos(2*x)),x)";
    let got = run_line(line).expect("integrate");
    assert_eq!(got, "ln(abs(1-sin(2*x)))*-1/2");
}
