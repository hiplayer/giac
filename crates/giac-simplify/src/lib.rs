//! Re-exports symbolic simplification from `giac-core` (facade crate per migration plan).

#![deny(unsafe_code)]
#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

pub use giac_core::{expand, factor, normal, ratnormal};
