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

pub use context::{Assumption, Context, RelationAssumption, expr_mentions_ident};
pub use stmt::{exec_stmt, exec_stmts, Stmt, StmtResult};
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
    algext_sqrt_branches, algext_square_roots, common_ext, contains_algext, try_as_algext_data,
    AlgExtData,
    try_rootof_to_algext,
};
pub use algebra::alg_ext_c::{canonicalize_to_algext_c, AlgExtCData};
pub use algebra::ext_tower::{CommonFieldPair, ExtensionField, ExtensionTower, FieldEmbedding};
pub use algebra::poly::{
    expr_contains_alg_coeff, expr_to_poly, poly_alg_from_expr, poly_mod_to_expr, poly_to_expr,
    ratio_to_expr, vars_from_expr,
};
pub use algebra::poly_conv::{
    ERR_ALG_EXT_COEFF, ERR_ALG_EXT_C_COEFF, ERR_POLY_ALG_NO_ALG_COEFF, ERR_POLY_ALG_UNIMPL,
    ERR_ROOTOF_COEFF,
};
pub use num_util::{
    bigint_to_i64, bigint_to_nonneg_u32, bigint_to_u32_abs, reduce_rational_pair,
};
pub use simplify::simplify;
