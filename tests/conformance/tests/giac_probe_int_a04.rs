//! Manual SymPy probe — not run in CI.

use giac_conformance::{run_line, sympy_verify_lines};

#[test]
#[ignore = "manual SymPy integrate probe"]
fn probe_int_a04_x_over_x_squared_plus_one() {
    let line = "integrate(x/(x^2+1),x)";
    let got = run_line(line).expect("eval");
    eprintln!("got={got}");
    let ok = sympy_verify_lines(&[line.to_string()], &[got.clone()])
        .map(|r| r[0].ok)
        .unwrap_or(false);
    assert!(ok, "got {got}");
}
