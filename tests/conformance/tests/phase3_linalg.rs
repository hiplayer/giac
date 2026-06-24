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
    assert_eq!(rref, "[[1,0,-1],[0,1,2]]");

    let ls = run_line("linsolve([2*x+y=3,x-y=1],[x,y])")?;
    assert_eq!(ls, "[1^-1*4/3,1^-1*1/3]");

    let cp = run_line("charpoly([[1,2],[3,4]],x)")?;
    assert_eq!(cp, "x^2-5*x-2");
    Ok(())
}

#[test]
fn test_linalg_ext_partial() -> Result<(), String> {
    let ker = run_line("ker([[1,2],[3,6]])")?;
    assert_eq!(ker, "[[-1*2,1]]");

    let j = run_line("jordan([[1,1],[0,1]])")?;
    assert_eq!(j, "[[1,1],[0,1]],matrix[[1,0],[0,1]]");

    let egv = run_line("egv([[4,1,-2],[1,2,-1],[2,1,0]])")?;
    assert_eq!(egv, "\"Not diagonalizable at eigenvalue 2\"");
    Ok(())
}

#[test]
fn test_linalg_decomp_numeric() -> Result<(), String> {
    let lu = run_line("lu([[3,5],[4,5]])")?;
    assert_eq!(
        lu,
        "[1,0],matrix[[1,0],[3/4,1]],matrix[[4,5],[0,5/4]]"
    );

    let qr = run_line("qr([[3,5],[4,5]])")?;
    assert_eq!(
        qr,
        "matrix[[3/5,4/5],[4/5,-3/5]],matrix[[5,7],[0,1]]"
    );

    let svd = run_line("svd([[1,2],[3,4]])")?;
    assert_eq!(
        svd,
        "matrix[[0.4045535848,-0.9145142957],[0.9145142957,0.4045535848]],[5.464985704,0.3659661906],matrix[[0.5760484368,0.8174155605],[0.8174155605,-0.5760484368]]"
    );
    Ok(())
}

#[test]
fn test_linalg_gramschmidt() -> Result<(), String> {
    let s = run_line("gramschmidt([1,1+x],(p,q)->integrate(p*q,x,-1,1))")?;
    assert_eq!(s, "[sqrt(2)^-1,x*sqrt(2/3)^-1]");
    verify_sympy(
        "gramschmidt([1,1+x],(p,q)->integrate(p*q,x,-1,1))",
        &s,
    )?;
    Ok(())
}
