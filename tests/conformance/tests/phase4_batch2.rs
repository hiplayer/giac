//! Phase 4 Batch 2 acceptance (GIAC-205, 206, 209, 211, 214).

use giac_conformance::{run_line, verify_sympy};

#[test]
fn batch2_giac205_solve_quadratic_rootof() -> Result<(), String> {
    let got = run_line("solve(t^2-2=0,t)")?;
    verify_sympy("solve(t^2-2=0,t)", &got)?;
    Ok(())
}

#[test]
fn batch2_giac206_sturmab() -> Result<(), String> {
    let got = run_line("sturmab(x^2*(x^3+2),x,-2,0)")?;
    verify_sympy("sturmab(x^2*(x^3+2),x,-2,0)", &got)?;
    Ok(())
}

#[test]
fn batch2_giac209_fsolve() -> Result<(), String> {
    let got = run_line("fsolve(x^2-2,x)")?;
    verify_sympy("fsolve(x^2-2,x)", &got)?;
    Ok(())
}

#[test]
fn batch2_giac211_integrate_sin_squared() -> Result<(), String> {
    let got = run_line("integrate(sin(x)^2,x)")?;
    verify_sympy("integrate(sin(x)^2,x)", &got)?;
    Ok(())
}

#[test]
fn batch2_giac211_integrate_tan_via_simplify() -> Result<(), String> {
    let got = run_line("simplify(int(tan(x),x))")?;
    verify_sympy("integrate(tan(x),x)", &got)?;
    Ok(())
}

#[test]
fn batch2_giac214_derive_multivariate() -> Result<(), String> {
    let line = "derive(2*x^2*y-x*z^3,[x,y,z])";
    let got = run_line(line)?;
    verify_sympy(line, &got)?;
    Ok(())
}
