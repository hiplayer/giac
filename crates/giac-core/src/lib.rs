#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod algebra;
mod context;
mod display;
mod error;
mod eval;
mod eval_poly;
#[cfg(test)]
mod eval_poly_tests;
mod expr;
mod ident;
mod integrate;
pub mod limits;
mod linalg;
mod matrix;
mod num_util;
mod simplify;
mod stmt;

pub use context::{Assumption, Context};
pub use display::format_expr;
pub use error::EvalError;
pub use eval::eval;
pub use expr::{Expr, ExprArc, FuncKind, RelOp};
pub use ident::Ident;
pub use limits::MAX_POLY_EXPONENT;
pub use algebra::{assert_equiv, expand, factor, is_zero, normal, ratnormal, sub};
pub use simplify::simplify;
pub use stmt::{exec_stmt, Stmt, StmtResult};
