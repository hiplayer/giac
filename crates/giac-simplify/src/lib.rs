//! Re-exports symbolic simplification from `giac-core` (facade crate per migration plan).

#![deny(unsafe_code)]

pub use giac_core::{expand, factor, normal, ratnormal};
