#![deny(unsafe_code)]

mod context;
mod display;
mod error;
mod eval;
mod expr;
mod ident;
mod simplify;
mod stmt;

pub use context::{Assumption, Context};
pub use display::format_expr;
pub use error::EvalError;
pub use eval::eval;
pub use expr::{Expr, ExprArc, FuncKind, RelOp};
pub use ident::Ident;
pub use simplify::simplify;
pub use stmt::{exec_stmt, Stmt, StmtResult};
