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
pub mod float_format;
mod ident;
mod integrate;
mod diff;
pub mod limits;
mod linalg_plugin;
mod solve_plugin;
mod matrix;
mod num_util;
mod simplify;
mod stmt;

pub use context::{Assumption, Context};
pub use display::format_expr;
pub use error::EvalError;
pub use eval::{eval, eval_subst_map};
pub use expr::{Expr, ExprArc, FuncKind, RelOp};
pub use float_format::format_float;
pub use ident::Ident;
pub use limits::MAX_POLY_EXPONENT;
pub use linalg_plugin::LinalgPlugin;
pub use solve_plugin::SolvePlugin;
pub use algebra::poly::{expr_to_poly, poly_to_expr};
pub use num_util::{bigint_to_i64, reduce_rational_pair};
pub use algebra::{assert_equiv, expand, factor, is_zero, normal, ratnormal, sub};
pub use simplify::simplify;
pub use stmt::{exec_stmt, Stmt, StmtResult};
