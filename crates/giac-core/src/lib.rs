#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod algebra;
mod algebra_plugin;
mod context;
mod display;
mod error;
mod eval;
mod eval_poly;
#[cfg(test)]
mod eval_poly_tests;
mod expr;
pub mod float_format;
mod calculus_plugin;
mod ident;
pub mod limits;
mod linalg_plugin;
mod ode_plugin;
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
pub use algebra_plugin::AlgebraPlugin;
pub use calculus_plugin::CalculusPlugin;
pub use linalg_plugin::LinalgPlugin;
pub use ode_plugin::OdePlugin;
pub use solve_plugin::SolvePlugin;
pub use algebra::alg_ext::{
    algext_square_roots, common_ext, contains_algext, try_as_algext_data, AlgExtData,
    try_rootof_to_algext,
};
pub use algebra::poly::{
    expr_to_poly, poly_mod_to_expr, poly_to_expr, ratio_to_expr, vars_from_expr,
};
pub use num_util::{
    bigint_to_i64, bigint_to_nonneg_u32, bigint_to_u32_abs, reduce_rational_pair,
};
pub use simplify::simplify;
pub use stmt::{exec_stmt, Stmt, StmtResult};
