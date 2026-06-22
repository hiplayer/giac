use std::fmt;
use std::sync::Arc;

use crate::{EvalError, Expr};

/// A symbol name (variable or function).
#[derive(Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Ident(Arc<str>);

impl Ident {
    pub fn new(s: impl AsRef<str>) -> Self {
        Self(Arc::from(s.as_ref()))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn is_imaginary_unit(&self) -> bool {
        self.as_str() == "i" || self.as_str() == "ii"
    }
}

impl fmt::Display for Ident {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for Ident {
    fn from(s: &str) -> Self {
        Self::new(s)
    }
}

/// Parse a single variable symbol from an expression leaf.
pub fn ident_from_expr(e: &Expr) -> Result<Ident, EvalError> {
    match e {
        Expr::Symbol(id) => Ok(id.clone()),
        _ => Err(EvalError::TypeError("variable name expected")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ident_equality_and_display() {
        let a = Ident::new("x");
        let b = Ident::new("x");
        assert_eq!(a, b);
        assert_eq!(a.to_string(), "x");
    }

    #[test]
    fn imaginary_unit_detection() {
        assert!(Ident::new("i").is_imaginary_unit());
        assert!(Ident::new("ii").is_imaginary_unit());
        assert!(!Ident::new("x").is_imaginary_unit());
    }

    #[test]
    fn ident_from_str() {
        let id: Ident = "y".into();
        assert_eq!(id.as_str(), "y");
    }

    #[test]
    fn ident_from_expr_symbol() {
        let id = ident_from_expr(Expr::sym("z").as_ref()).unwrap();
        assert_eq!(id.as_str(), "z");
        assert!(ident_from_expr(Expr::int(1).as_ref()).is_err());
    }
}
