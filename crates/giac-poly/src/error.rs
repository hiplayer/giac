use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
pub enum PolyError {
    #[error("division by zero")]
    DivisionByZero,
    #[error("type error: {0}")]
    TypeError(&'static str),
    #[error("not implemented: {0}")]
    NotImplemented(&'static str),
}

pub type PolyResult<T> = Result<T, PolyError>;
