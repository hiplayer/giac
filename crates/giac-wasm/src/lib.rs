#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

use giac_calculus::xcas_default;
use giac_core::{exec_stmt, format_expr, StmtResult};
use giac_parse::parse_program;
use wasm_bindgen::prelude::*;

/// Evaluate a giac script fragment and return the formatted result.
///
/// Parses `input` as a program (statements separated by `;`), executes with the
/// default Xcas context (linear algebra, solve, calculus plugins), and joins
/// expression results with commas. Returns an empty string when there is no value.
pub fn eval_to_string(input: &str) -> Result<String, String> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Ok(String::new());
    }
    let script = if trimmed.ends_with(';') {
        trimmed.to_string()
    } else {
        format!("{trimmed};")
    };
    let mut ctx = xcas_default();
    let stmts = parse_program(&script, &ctx).map_err(|e| format!("parse error: {e}"))?;
    let mut parts = Vec::new();
    for stmt in &stmts {
        match exec_stmt(stmt, &mut ctx).map_err(|e| e.to_string())? {
            StmtResult::Value(v) | StmtResult::Assign { value: v, .. } => {
                parts.push(format_expr(v.as_ref()));
            }
            StmtResult::NoValue => {}
        }
    }
    Ok(parts.join(","))
}

/// WASM export: evaluate input; errors become `"error: …"` strings.
#[wasm_bindgen(js_name = evalToString)]
pub fn eval_to_string_wasm(input: &str) -> String {
    eval_to_string(input).unwrap_or_else(|e| format!("error: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eval_one_plus_two() {
        assert_eq!(eval_to_string("1+2").unwrap(), "3");
    }

    #[test]
    fn eval_empty_input() {
        assert_eq!(eval_to_string("").unwrap(), "");
    }

    #[test]
    fn eval_parse_error() {
        assert!(eval_to_string("1+").is_err());
    }
}
