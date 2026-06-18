#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

//! Symbolic calculus: diff, integrate, limit, series, Risch subset.
//!
//! # API stability
//!
//! 三层分类与注释格式见 [`.doc/giac-calculus-api-stability.md`](../../.doc/giac-calculus-api-stability.md)
//! 与通用 [`.doc/algorithm-expr-api.md`](../../.doc/algorithm-expr-api.md)。
//!
//! | 标记 | 含义 |
//! |------|------|
//! | `/// **Stable**` | 可跨模块调用，有契约与单测 |
//! | `/// **Partial**` | 窄路径或待收敛；注释含退役条件 |
//! | `/// **Pipeline**` | crate 内编排，非规范形保证 |
//! | `// **Pipeline private**` | 模块内私有辅助 |
//!
//! | Tier | 入口 |
//! |------|------|
//! | **Stable** | `diff`, `integrate`, `eval_limit`, `eval_series`, `depends_on_var`, `limit_engine` 文档 §3 清单 |
//! | **Partial** | `eval_risch`, `try_integrate_*`, `first_order_exp_vanishing_epsilon`, `limit_at_plus_infinity_fallback` |
//! | **Pipeline** | `limit_preprocess_*`, `mrv_lead_term_*`, `remove_lnexp` |
//!
//! `limit_engine` 专项契约：[limit-engine-expr-api.md](../../.doc/limit-engine-expr-api.md)、
//! [exp-diff-expr-api.md](../../.doc/exp-diff-expr-api.md)。

mod diff;
mod eval_diff;
mod eval_integrate;
mod expr_util;
mod integrate;
mod integrate_heuristics;
mod partfrac_integrate;
mod plugin;
mod stubs;
mod limit;
mod limit_engine;
mod series;
mod risch;

pub use diff::diff;
pub use eval_diff::eval_diff;
pub use eval_integrate::eval_integrate;
pub use integrate::integrate;
pub use plugin::{install_calculus, xcas_default, DefaultCalculusPlugin};
pub use limit::eval_limit;
pub use series::eval_series;
pub use risch::{eval_risch, hermite_reduce, pow2expln, risch_tower, rlvarx, rothstein_trager_integrate, HermiteTerm, RischTowerError};
