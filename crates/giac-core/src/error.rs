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

impl From<giac_poly::PolyError> for EvalError {
    fn from(e: giac_poly::PolyError) -> Self {
        match e {
            giac_poly::PolyError::DivisionByZero => EvalError::DivisionByZero,
            giac_poly::PolyError::TypeError(m) => EvalError::TypeError(m),
            giac_poly::PolyError::NotImplemented(m) => EvalError::NotImplemented(m),
        }
    }
}

/// Legacy mapping used by solve/calculus (division-by-zero → `TypeError` message).
pub fn poly_error_compat(e: giac_poly::PolyError) -> EvalError {
    match e {
        giac_poly::PolyError::DivisionByZero => EvalError::TypeError("division by zero"),
        other => EvalError::from(other),
    }
}
