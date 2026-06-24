//! Phase 2 polynomial ring conformance tests against giac reference output.

use giac_conformance::{run_line, run_script, upstream_root, verify_sympy};
use giac_core::Context;
use giac_parse::parse_program;
use std::fs;

#[test]
fn phase2_scripts_parseable() {
    let root = upstream_root().join("bin");
    for name in [
        "test_poly",
        "test_poly_ext",
        "test_modular",
        "test_factor",
        "test_groebner",
    ] {
        let path = root.join(name);
        let input = fs::read_to_string(&path).unwrap();
        let ctx = Context::xcas_default();
        parse_program(&input, &ctx).unwrap_or_else(|e| panic!("parse {name}: {e}"));
    }
}

#[test]
fn test_poly_matches_giac() -> Result<(), String> {
    let path = upstream_root().join("bin/test_poly");
    let lines = giac_conformance::script_lines(&path)?;
    let got = run_script(&path)?;
    let want = ["x-1", "x^2+x+1", "0", "6", "2*x*y"];
    assert_eq!(got, want.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    for (line, out) in lines.iter().zip(got.iter()) {
        verify_sympy(line, out)?;
    }
    Ok(())
}

#[test]
fn test_poly_ext_core() -> Result<(), String> {
    let cases = [
        ("horner(x^4+2*x^3-3*x^2+x-2,1)", "-1"),
        ("resultant(x^2-1,x^3-1,x)", "0"),
        ("simp2(x^3-1,x^2-1)", "[x^2+x+1,x+1]"),
    ];
    for (line, want) in cases {
        let got = run_line(line)?;
        assert_eq!(got, want, "{line}");
        verify_sympy(line, &got)?;
    }
    Ok(())
}

#[test]
fn test_modular_core() -> Result<(), String> {
    let cases = [
        ("modp(x^2+1,3)", "x^2+1"),
        ("smod(17,5)", "2"),
        ("irem(17,5)", "2"),
        ("factor(x^4-1) mod 2", "x^4+1 mod 2"),
    ];
    for (line, want) in cases {
        let got = run_line(line)?;
        assert_eq!(got, want, "{line}");
        verify_sympy(line, &got)?;
    }
    let chin = run_line("chinrem([x+2,x^2+1],[x+1,x^2+x+1])")?;
    verify_sympy("chinrem([x+2,x^2+1],[x+1,x^2+x+1])", &chin)?;

    let line = "normal(((2*x+1) % 13)^5)";
    let got = run_line(line)?;
    assert_eq!(
        got,
        "(6 % 13)*x^5+(2 % 13)*x^4+(2 % 13)*x^3+(1 % 13)*x^2+(10 % 13)*x+(1 % 13)"
    );
    verify_sympy(line, &got)?;
    Ok(())
}

#[test]
fn test_factor_core() -> Result<(), String> {
    let line = "factor(x^4-1)";
    let got = run_line(line)?;
    assert_eq!(got, "(x-1)*(x+1)*(x^2+1)");
    verify_sympy(line, &got)?;

    let line = "normal((x+1)^3)";
    let got = run_line(line)?;
    assert_eq!(got, "x^3+3*x^2+3*x+1");
    verify_sympy(line, &got)?;
    Ok(())
}

#[test]
fn test_groebner_runs() -> Result<(), String> {
    let got = run_script(&upstream_root().join("bin/test_groebner"))?;
    assert_eq!(got.len(), 2);
    assert_eq!(got[0], "1/2*y^2-1");
    assert_eq!(got[1], "2*y^2-1");
    for (line, out) in [
        (
            "greduce(x*y-1,[x^2-y^2,2*x*y-y^2,y^3],[x,y,z])",
            got[0].as_str(),
        ),
        (
            "greduce(x^2+y^2-1,[x^2-y^2,2*x*y-1],[x,y])",
            got[1].as_str(),
        ),
    ] {
        verify_sympy(line, out)?;
    }
    Ok(())
}

#[test]
fn gcd_mod_13_linear() -> Result<(), String> {
    let line = "gcd((2*x^2+5) % 13,(5*x^2+2*x-3) % 13)";
    let got = run_line(line)?;
    assert_eq!(got, "7*x+1 mod 13");
    verify_sympy(line, &got)?;
    Ok(())
}
