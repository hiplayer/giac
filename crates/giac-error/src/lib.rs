#![deny(unsafe_code)]

use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
pub enum EvalError {
    #[error("too few arguments for {0}")]
    TooFewArgs(&'static str),
    #[error("too many arguments for {0}")]
    TooManyArgs(&'static str),
    #[error("division by zero")]
    DivisionByZero,
    #[error("type error: {0}")]
    TypeError(&'static str),
    #[error("unknown function: {0}")]
    UnknownFunction(String),
    #[error("unknown variable: {0}")]
    UnknownVariable(String),
    #[error("not implemented: {0}")]
    NotImplemented(&'static str),
}
