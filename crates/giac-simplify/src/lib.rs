#![deny(unsafe_code)]
//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-simplify-api-stability.md`.
//!
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

//! Symbolic simplification: expand, normal, ratnormal, factor, equivalence.
//!
//! # API stability
//!
//! 三层分类与完整清单见 [`.doc/giac-simplify-api-stability.md`](../../.doc/giac-simplify-api-stability.md)
//! 与 [`.doc/algorithm-expr-api.md`](../../.doc/algorithm-expr-api.md)。
//! 上游缺口索引：[`.doc/issues/GIAC-simplify-poly-upstream-gaps.md`](../../.doc/issues/GIAC-simplify-poly-upstream-gaps.md) §1。
//!
//! | Tier | Public API | Scope |
//! |------|------------|-------|
//! | **Stable** | `normal`, `ratnormal`, `expand_polynomial`, `expand_with_policy`, `assert_equiv`, `is_zero`, `sub`, `ifactor`, `install_simplify`, `xcas_default` | Polynomial / rational / integer paths with documented `Result` semantics |
//! | **Stable (bounded)** | `factor` | Univariate and structural factorization; multivariate Hensel gaps live in `giac-poly` |
//! | **Stable (legacy default)** | `expand` | Same as `expand_with_policy(..., Full)`; prefer `expand_polynomial` when `exp` shapes must be preserved |
//! | **Partial** | `texpand`, `lin`, `halftan` | Rule-table subsets of upstream `usual.cc` / `lin.cc`; general input → `NotImplemented` |
//!
//! **Temporary (private):** `ratnormal_algext` (shim); `canonical_radical` / `inv_sqrt_to_mul` (equiv drift only).

mod expand;
mod equiv;
mod ifactor;
mod factor;
mod ratnormal;
mod trig;
mod plugin;

pub mod test_verify;

pub use expand::{expand, expand_polynomial, expand_with_policy, normal, ExpandPolicy};
pub use equiv::{assert_equiv, is_zero, sub};
pub use factor::factor;
pub use ifactor::ifactor;
pub use ratnormal::ratnormal;
pub use trig::{halftan, lin, texpand};
pub use plugin::{install_simplify, xcas_default, DefaultAlgebraPlugin};
