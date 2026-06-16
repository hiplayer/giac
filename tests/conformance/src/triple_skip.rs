//! Centralized skip / known-gap predicates for triple conformance tests.

/// Phase 2: output may differ from giac while SymPy still passes.
pub fn phase2_format_diff(line: &str) -> bool {
    matches!(
        line,
        "gauss(2*x*y,[x,y])"
            | "egcd(x^2+2*x+1,x^2-1)"
            | "abcuv(x^2+2*x+1,x^2-1,x+1)"
            | "lcm(x^2+2*x+1,x^2-1)"
            | "roots(x^3-1,x)"
            | "partfrac(1/(x^2-1),x)"
            | "chinrem([x+2,x^2+1],[x+1,x^2+x+1])"
    )
}

pub fn phase2_sympy_gap(line: &str) -> bool {
    matches!(line, "roots(x^3-1,x)")
}

pub fn phase2_giac_gap(line: &str) -> bool {
    matches!(line, "chinrem([x+2,x^2+1],[x+1,x^2+x+1])")
}

/// Phase 3: matrix / linalg lines where giac-rs and giac print differently.
pub fn phase3_format_diff(line: &str) -> bool {
    matches!(
        line,
        "[[1,2],[3,4]]^2"
            | "tran([[1,2],[3,4]])"
            | "ker([[1,2],[3,6]])"
            | "image([[1,2],[3,6]])"
            | "pcar([[4,1,-2],[1,2,-1],[2,1,0]])"
    )
}

pub fn phase3_sympy_gap(_line: &str) -> bool {
    false
}

pub fn phase3_skip(_line: &str) -> bool {
    false
}

/// Trig simplification: giac-rs and giac may print different but equivalent forms.
pub fn trig_format_diff(line: &str) -> bool {
    matches!(line, "texpand(cos(3*x))")
}

pub fn phase3_numerical_decomp(line: &str) -> bool {
    line.starts_with("lu(") || line.starts_with("qr(") || line.starts_with("svd(")
}

/// Phase 4: skip lines whose eval is not implemented yet (parse-only or hard gaps).
pub fn phase4_skip(line: &str) -> bool {
    let line = line.trim().trim_end_matches(';');
    line.starts_with("proot(")
        || matches!(line, "proot(x^3-2)")
        || line.starts_with("partfrac(")
}
