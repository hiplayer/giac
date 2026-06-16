#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

mod diff;
mod eval_diff;
mod eval_integrate;
mod integrate;
mod integrate_heuristics;
mod partfrac_integrate;
mod plugin;
mod stubs;
mod limit;
mod series;
mod risch;

pub use diff::diff;
pub use eval_diff::eval_diff;
pub use eval_integrate::eval_integrate;
pub use integrate::integrate;
pub use plugin::{install_calculus, xcas_default, DefaultCalculusPlugin};
pub use limit::eval_limit;
pub use series::eval_series;
pub use risch::{eval_risch, risch_tower, rlvarx, RischTowerError};
