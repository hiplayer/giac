#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod lexer;
mod parser;

pub use lexer::Lexer;
pub use parser::{parse_compound_line, parse_program, ParseError};

use giac_core::Context;

/// Parse a giac script into statements with default Xcas context semantics.
pub fn parse_script(input: &str) -> Result<Vec<giac_core::Stmt>, ParseError> {
    parse_program(input, &Context::xcas_default())
}

#[cfg(test)]
mod integration {
    use giac_core::{eval, format_expr, Context};

    use super::parse_script;

    #[test]
    fn parse_and_eval_cas_basic_lines() {
        let lines = [
            "abs(1+2*i);",
            "gcd(45,75);",
            "conj((1+2*i)^2);",
        ];
        let expected = ["sqrt(5)", "15", "-3-4*i"];
        let ctx = Context::xcas_default();

        for (src, want) in lines.iter().zip(expected.iter()) {
            let stmts = parse_script(src).expect(src);
            assert_eq!(stmts.len(), 1);
            let value = match &stmts[0] {
                giac_core::Stmt::ExprStmt(e) => eval(e, &ctx).unwrap(),
                _ => panic!("expected expr stmt"),
            };
            assert_eq!(format_expr(value.as_ref()), *want, "for {src}");
        }
    }
}
