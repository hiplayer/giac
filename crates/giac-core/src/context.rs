use std::collections::HashMap;
use std::sync::Arc;

use crate::expr::ExprArc;
use crate::ident::Ident;
use crate::linalg_plugin::LinalgPlugin;
use crate::solve_plugin::SolvePlugin;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Assumption {
    Integer(Ident),
    Real(Ident),
    Positive(Ident),
}

#[derive(Clone)]
pub struct Context {
    pub vars: HashMap<Ident, ExprArc>,
    pub assumptions: Vec<Assumption>,
    pub complex_mode: bool,
    pub epsilon: f64,
    pub series_order: u32,
    pub float_digits: u32,
    pub(crate) linalg_plugin: Option<Arc<dyn LinalgPlugin>>,
    pub(crate) solve_plugin: Option<Arc<dyn SolvePlugin>>,
}

impl Default for Context {
    fn default() -> Self {
        Self {
            vars: HashMap::new(),
            assumptions: Vec::new(),
            complex_mode: true,
            epsilon: 1e-10,
            series_order: 6,
            float_digits: 12,
            linalg_plugin: None,
            solve_plugin: None,
        }
    }
}

impl Context {
    pub fn new() -> Self {
        Self::default()
    }

    /// Default Xcas script setup (`xcas_mode(0)`).
    pub fn xcas_default() -> Self {
        Self::default()
    }

    pub fn set(&mut self, name: Ident, value: ExprArc) {
        self.vars.insert(name, value);
    }

    pub fn get(&self, name: &Ident) -> Option<&ExprArc> {
        self.vars.get(name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Expr;

    #[test]
    fn variable_binding_roundtrip() {
        let mut ctx = Context::new();
        let x = Ident::new("x");
        ctx.set(x.clone(), Expr::int(42));
        assert_eq!(ctx.get(&x), Some(&Expr::int(42)));
    }
}
