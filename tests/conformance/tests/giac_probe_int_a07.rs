use giac_conformance::{run_line, sympy_verify_lines};

#[test]
fn probe_int_a07_one_over_one_plus_x_fourth() {
    let line = "integrate(1/(1+x^4),x)";
    let got = run_line(line).expect("eval");
    let ok = sympy_verify_lines(&[line.to_string()], &[got.clone()])
        .map(|r| r[0].ok)
        .unwrap_or(false);
    assert!(ok, "line={line} got={got}");
}
