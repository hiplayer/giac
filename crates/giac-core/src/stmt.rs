use std::sync::Arc;

use crate::context::Context;
use crate::error::EvalError;
use crate::eval::eval;
use crate::expr::ExprArc;
use crate::ident::Ident;

/// A program statement.
#[derive(Clone, Debug, PartialEq)]
pub enum Stmt {
    ExprStmt(ExprArc),
    Assign(Ident, ExprArc),
}

/// Result of executing one statement.
#[derive(Clone, Debug, PartialEq)]
pub enum StmtResult {
    Value(ExprArc),
    Assign { name: Ident, value: ExprArc },
    NoValue,
}

/// Execute a statement, updating context for assignments.
pub fn exec_stmt(stmt: &Stmt, ctx: &mut Context) -> Result<StmtResult, EvalError> {
    match stmt {
        Stmt::ExprStmt(expr) => {
            let value = eval(expr, ctx)?;
            Ok(StmtResult::Value(value))
        }
        Stmt::Assign(name, expr) => {
            let value = eval(expr, ctx)?;
            ctx.set(name.clone(), Arc::clone(&value));
            Ok(StmtResult::Assign {
                name: name.clone(),
                value,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Expr;

    #[test]
    fn assignment_updates_context() {
        let mut ctx = Context::new();
        let stmt = Stmt::Assign(Ident::new("x"), Expr::int(5));
        exec_stmt(&stmt, &mut ctx).unwrap();
        assert_eq!(ctx.get(&Ident::new("x")), Some(&Expr::int(5)));
    }
}
