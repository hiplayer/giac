//! Phase 3 linear algebra conformance tests.

use giac_conformance::{run_line, upstream_root, verify_sympy, PHASE3_SCRIPTS};

#[test]
fn phase3_scripts_parseable() {
    use giac_core::Context;
    use giac_parse::parse_program;
    use std::fs;

    let root = upstream_root().join("bin");
    for name in PHASE3_SCRIPTS {
        let path = root.join(name);
        let input = fs::read_to_string(&path).unwrap();
        let ctx = Context::xcas_default();
        parse_program(&input, &ctx).unwrap_or_else(|e| panic!("parse {name}: {e}"));
    }
}

#[test]
fn test_linalg_core() -> Result<(), String> {
    let cases = [
        ("[[1,2],[3,4]]^2", "[[7,10],[15,22]]"),
        ("det([[1,2],[3,4]])", "-2"),
        ("tran([[1,2],[3,4]])", "matrix[[1,3],[2,4]]"),
        ("trace([[1,2],[3,4]])", "5"),
    ];
    for (line, want) in cases {
        let got = run_line(line)?;
        assert_eq!(got, want, "{line}");
    }
    Ok(())
}

#[test]
fn test_linalg_rref_linsolve_charpoly() -> Result<(), String> {
    let rref = run_line("rref([[1,2,3],[4,5,6]])")?;
    assert!(rref.contains("1,0") && rref.contains("0,1"), "rref got {rref}");

    let ls = run_line("linsolve([2*x+y=3,x-y=1],[x,y])")?;
    assert!(ls.contains("4/3") || ls.contains("1/3"), "linsolve got {ls}");

    let cp = run_line("charpoly([[1,2],[3,4]],x)")?;
    assert!(
        cp.contains("x^2") && cp.contains("-5") && cp.contains("-2"),
        "charpoly got {cp}"
    );
    Ok(())
}

#[test]
fn test_linalg_ext_partial() -> Result<(), String> {
    let ker = run_line("ker([[1,2],[3,6]])")?;
    assert!(ker.contains("2") && ker.contains("1"), "ker got {ker}");

    let j = run_line("jordan([[1,1],[0,1]])")?;
    assert!(j.contains("[[1,0],[0,1]]"), "jordan got {j}");

    let egv = run_line("egv([[4,1,-2],[1,2,-1],[2,1,0]])")?;
    assert!(
        egv.contains("Not diagonalizable"),
        "egv got {egv}"
    );
    Ok(())
}

#[test]
fn test_linalg_decomp_numeric() -> Result<(), String> {
    let lu = run_line("lu([[3,5],[4,5]])")?;
    assert!(
        lu.contains("matrix[[") && lu.matches(',').count() > 3,
        "lu got {lu}"
    );

    let qr = run_line("qr([[3,5],[4,5]])")?;
    assert!(qr.contains("matrix[["), "qr got {qr}");

    let svd = run_line("svd([[1,2],[3,4]])")?;
    assert!(svd.contains("matrix[[") && svd.contains(','), "svd got {svd}");
    Ok(())
}

#[test]
fn test_linalg_gramschmidt() -> Result<(), String> {
    let s = run_line("gramschmidt([1,1+x],(p,q)->integrate(p*q,x,-1,1))")?;
    assert!(
        s.contains("sqrt"),
        "gramschmidt should be orthonormal with sqrt, got {s}"
    );
    verify_sympy(
        "gramschmidt([1,1+x],(p,q)->integrate(p*q,x,-1,1))",
        &s,
    )?;
    Ok(())
}
