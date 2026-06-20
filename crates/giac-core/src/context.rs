use std::collections::HashMap;
use std::sync::Arc;

use crate::expr::ExprArc;
use crate::ident::Ident;
use crate::algebra_plugin::AlgebraPlugin;
use crate::calculus_plugin::CalculusPlugin;
use crate::linalg_plugin::LinalgPlugin;
use crate::ode_plugin::OdePlugin;
use crate::solve_plugin::SolvePlugin;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Assumption {
    Integer(Ident),
    Real(Ident),
    Positive(Ident),
}

/// Symbolic relation stored by `assume` (e.g. `t>2`).
#[derive(Clone, Debug, PartialEq)]
pub struct RelationAssumption {
    pub relation: ExprArc,
}

#[derive(Clone)]
pub struct Context {
    pub vars: HashMap<Ident, ExprArc>,
    pub assumptions: Vec<Assumption>,
    /// Stack of relations from `assume(...)` (GIAC-204b).
    pub relation_assumptions: Vec<RelationAssumption>,
    pub complex_mode: bool,
    pub epsilon: f64,
    pub series_order: u32,
    pub float_digits: u32,
    pub(crate) linalg_plugin: Option<Arc<dyn LinalgPlugin>>,
    pub(crate) solve_plugin: Option<Arc<dyn SolvePlugin>>,
    pub(crate) algebra_plugin: Option<Arc<dyn AlgebraPlugin>>,
    pub(crate) calculus_plugin: Option<Arc<dyn CalculusPlugin>>,
    pub(crate) ode_plugin: Option<Arc<dyn OdePlugin>>,
    pub with_sqrt: bool,
}

impl Default for Context {
    fn default() -> Self {
        Self {
            vars: HashMap::new(),
            assumptions: Vec::new(),
            relation_assumptions: Vec::new(),
            complex_mode: true,
            epsilon: 1e-10,
            series_order: 6,
            float_digits: 12,
            linalg_plugin: None,
            solve_plugin: None,
            algebra_plugin: None,
            calculus_plugin: None,
            ode_plugin: None,
            with_sqrt: true,
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

    /// Push a relation from `assume(relation)`.
    pub fn push_relation_assumption(&mut self, relation: ExprArc) {
        self.relation_assumptions
            .push(RelationAssumption { relation });
    }

    /// Remove assumptions involving `name` (upstream `purgenoassume`).
    pub fn purge_relation_assumptions_for(&mut self, name: &Ident) {
        self.relation_assumptions.retain(|a| {
            !expr_mentions_ident(a.relation.as_ref(), name)
        });
        self.assumptions.retain(|a| match a {
            Assumption::Integer(id)
            | Assumption::Real(id)
            | Assumption::Positive(id) => id != name,
        });
        self.vars.remove(name);
    }
}

/// True if `e` contains symbol `name`.
pub fn expr_mentions_ident(e: &crate::expr::Expr, name: &Ident) -> bool {
    use crate::expr::Expr;
    match e {
        Expr::Symbol(id) => id == name,
        Expr::Add(terms) | Expr::Mul(terms) => {
            terms.iter().any(|t| expr_mentions_ident(t.as_ref(), name))
        }
        Expr::Pow(base, exp) => {
            expr_mentions_ident(base.as_ref(), name)
                || expr_mentions_ident(exp.as_ref(), name)
        }
        Expr::Frac(n, d) | Expr::Mod(n, d) => {
            expr_mentions_ident(n.as_ref(), name) || expr_mentions_ident(d.as_ref(), name)
        }
        Expr::Complex(re, im) => {
            expr_mentions_ident(re.as_ref(), name) || expr_mentions_ident(im.as_ref(), name)
        }
        Expr::Relation(_, lhs, rhs) => {
            expr_mentions_ident(lhs.as_ref(), name) || expr_mentions_ident(rhs.as_ref(), name)
        }
        Expr::Func(_, args) => args.iter().any(|a| expr_mentions_ident(a.as_ref(), name)),
        Expr::List(items) | Expr::Seq(items) => items
            .iter()
            .any(|a| expr_mentions_ident(a.as_ref(), name)),
        Expr::Matrix(items) => items
            .iter()
            .flatten()
            .any(|a| expr_mentions_ident(a.as_ref(), name)),
        Expr::GiacMatrix(rows) => rows
            .iter()
            .flatten()
            .any(|a| expr_mentions_ident(a.as_ref(), name)),
        Expr::Int(_) | Expr::Rat(_) | Expr::AlgExt(_) | Expr::AlgExtC(_) | Expr::Str(_) | Expr::Undefined => false,
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
