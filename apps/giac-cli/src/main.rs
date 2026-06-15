use std::env;
use std::fs;
use std::io::{self, Read};
use std::path::PathBuf;

use anyhow::{Context as AnyhowCtx, Result};
use giac_core::{exec_stmt, format_expr, StmtResult};
use giac_calculus::xcas_default;
use giac_parse::parse_program;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<()> {
    let script = env::args().nth(1).map(PathBuf::from);
    let input = read_input(script.as_deref())?;
    let mut ctx = xcas_default();
    let stmts = parse_program(&input, &ctx).map_err(|e| anyhow::anyhow!("parse error: {e}"))?;

    let mut first = true;
    for stmt in &stmts {
        match exec_stmt(stmt, &mut ctx)? {
            StmtResult::Value(v) => {
                if !first {
                    print!(",");
                }
                print!("{}", format_expr(v.as_ref()));
                first = false;
            }
            StmtResult::Assign { value, .. } => {
                if !first {
                    print!(",");
                }
                print!("{}", format_expr(value.as_ref()));
                first = false;
            }
            StmtResult::NoValue => {}
        }
    }
    if !first {
        println!();
    }
    Ok(())
}

#[cfg(target_arch = "wasm32")]
fn main() {
    panic!("giac-cli is not available on wasm32; use giac-wasm instead");
}

#[cfg(not(target_arch = "wasm32"))]
fn read_input(path: Option<&std::path::Path>) -> Result<String> {
    match path {
        Some(p) => fs::read_to_string(p)
            .with_context(|| format!("failed to read {}", p.display())),
        None => {
            let mut buf = String::new();
            io::stdin().read_to_string(&mut buf)?;
            Ok(buf)
        }
    }
}
