//! Re-exports for solve functions implemented in dedicated modules.

//!
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-solve-api-stability.md`.
//!
pub use crate::fsolve::eval_fsolve;
pub use crate::sturm::{eval_sturm, eval_sturmab};
