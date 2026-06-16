use giac_conformance::{run_line, sympy_verify_lines};

fn assert_sympy(line: &str) {
    let got = run_line(line).expect("eval");
    let ok = sympy_verify_lines(&[line.to_string()], &[got.clone()])
        .map(|r| r[0].ok)
        .unwrap_or(false);
    assert!(ok, "line={line} got={got}");
}

#[test]
fn probe_integrate_x_over_x_squared_plus_one_squared() {
    assert_sympy("integrate(x/(x^2+1)^2,x)");
}

#[test]
fn probe_integrate_one_over_x_fourth_plus_one_squared() {
    assert_sympy("integrate(1/(x^4+1)^2,x)");
}

#[test]
fn probe_integrate_one_over_x_fourth_plus_one() {
    assert_sympy("integrate(1/(x^4+1),x)");
}
