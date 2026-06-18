#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

//! Symbolic simplification: expand, normal, ratnormal, factor, equivalence.
//!
//! # API stability
//!
//! | Tier | Public API | Scope |
//! |------|------------|-------|
//! | **Stable** | `normal`, `ratnormal`, `expand_polynomial`, `expand_with_policy`, `assert_equiv`, `is_zero`, `sub`, `ifactor`, `install_simplify`, `xcas_default` | Polynomial / rational / integer paths with documented `Result` semantics |
//! | **Stable (bounded)** | `factor` | Univariate and structural factorization; multivariate Hensel gaps live in `giac-poly` |
//! | **Stable (legacy default)** | `expand` | Same as `expand_with_policy(..., Full)`; prefer `expand_polynomial` when `exp` shapes must be preserved |
//! | **Partial** | `texpand`, `lin`, `halftan` | Rule-table subsets of upstream `usual.cc` / `lin.cc`; general input → `NotImplemented` |
//!
//! Internal helpers are **pipeline-private** unless noted `canonical_*` (equiv drift only).
//!
//! # Known technical debt
//!
//! - No crate-wide `canonical_*` naming layer; AST drift handled ad hoc (`equiv::canonical_radical` only).
//! - `ratnormal` on `AlgExt` delegates to `eval`, not `ext_reduce` (GIAC-algext-adoption A-03).
//! - Many upstream `usual.cc` / `subst.cc` builtins unported (`tlin`, `trig2exp`, `reorder`, `simplify`, …).

mod expand;
mod equiv;
mod ifactor;
mod factor;
mod ratnormal;
mod trig;
mod plugin;

pub use expand::{expand, expand_polynomial, expand_with_policy, normal, ExpandPolicy};
pub use equiv::{assert_equiv, is_zero, sub};
pub use factor::factor;
pub use ifactor::ifactor;
pub use ratnormal::ratnormal;
pub use trig::{halftan, lin, texpand};
pub use plugin::{install_simplify, xcas_default, DefaultAlgebraPlugin};
