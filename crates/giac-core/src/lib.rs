#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod algebra;
mod algebra_plugin;
mod context;
mod expr_shape;
mod display;
mod error;
mod eval;
mod eval_poly;
mod expr_rational;
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
pub use expr_shape::{
    is_cos_of_var, is_ln_of_var, is_sin_of_var, is_sin_of_var_expr, is_var, is_var_expr, var_to_expr,
};
pub use ident::{ident_from_expr, Ident};
pub use limits::MAX_POLY_EXPONENT;
pub use algebra_plugin::AlgebraPlugin;
pub use calculus_plugin::CalculusPlugin;
pub use linalg_plugin::LinalgPlugin;
pub use ode_plugin::OdePlugin;
pub use solve_plugin::SolvePlugin;
pub use algebra::alg_ext::{
    algext_sqrt_branches, algext_square_roots, common_ext, contains_algext, quadratic_rootof_branches,
    rootof_from_minpoly, try_as_algext_data,
    AlgExtData,
    try_rootof_to_algext,
};
pub use algebra::alg_ext_c::{canonicalize_to_algext_c, AlgExtCData};
pub use algebra::ext_tower::{AlignedElements, CommonFieldPair, ExtensionField, ExtensionTower, FieldEmbedding};
pub use algebra::poly::{
    algext_poly_to_expr, expr_contains_alg_coeff, expr_to_poly, poly_alg_from_expr,
    poly_algext_from_poly, poly_mod_to_expr, poly_to_expr, univariate_poly_to_poly1_expr, ratio_to_expr, vars_from_expr,
    PolyAlgExt,
};
pub use algebra::poly_alg_coeff::AlgExtCPolyCoeff;
pub use algebra::poly_alg_ops::{
    align_algext_polys, div_rem_wrt_algext, ensure_common_field_for_polys, infer_ambient_field,
    monic_wrt_algext, normalize_algext_poly, quo_exact_wrt_algext,
};
pub use algebra::poly_roots::{poly_algext_roots, poly_algext_roots_for_ctx};
pub use algebra::poly_conv::{
    ERR_ALG_EXT_COEFF, ERR_ALG_EXT_C_COEFF, ERR_POLY_ALG_NO_ALG_COEFF, ERR_ROOTOF_COEFF,
};
pub use expr_rational::{expr_to_rational_polys, expr_to_rational_polys_eval};
pub use num_util::{
    bigint_to_i64, bigint_to_nonneg_u32, bigint_to_u32_abs, integer_nth_root, integer_sqrt,
    reduce_rational_pair,
};
pub use simplify::simplify;
