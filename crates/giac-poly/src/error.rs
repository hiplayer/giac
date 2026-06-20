//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
pub use giac_error::EvalError;

pub type PolyResult<T> = Result<T, EvalError>;
