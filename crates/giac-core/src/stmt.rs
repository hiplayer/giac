use std::sync::Arc;

use crate::context::Context;
use crate::error::EvalError;
use crate::eval::eval;
use crate::expr::{Expr, ExprArc, FuncKind};
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

/// Execute a statement, updating context for assignments and assume/purge.
pub fn exec_stmt(stmt: &Stmt, ctx: &mut Context) -> Result<StmtResult, EvalError> {
    match stmt {
        Stmt::ExprStmt(expr) => {
            if let Some(value) = try_exec_special_expr(expr, ctx)? {
                return Ok(StmtResult::Value(value));
            }
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

/// Execute a comma-separated script (GIAC-204b: `assume,...,purge`).
pub fn exec_stmts(stmts: &[Stmt], ctx: &mut Context) -> Result<StmtResult, EvalError> {
    let mut last = StmtResult::NoValue;
    for stmt in stmts {
        last = exec_stmt(stmt, ctx)?;
    }
    Ok(last)
}

fn try_exec_special_expr(expr: &ExprArc, ctx: &mut Context) -> Result<Option<ExprArc>, EvalError> {
    match expr.as_ref() {
        Expr::Func(FuncKind::Assume, args) => Ok(Some(eval_assume(args, ctx)?)),
        Expr::Func(FuncKind::Purge, args) => Ok(Some(eval_purge(args, ctx)?)),
        _ => Ok(None),
    }
}

fn eval_assume(args: &[ExprArc], ctx: &mut Context) -> Result<ExprArc, EvalError> {
    if args.len() != 1 {
        return Err(EvalError::TooFewArgs("assume"));
    }
    let relation = eval(args[0].as_ref(), ctx)?;
    ctx.push_relation_assumption(Arc::clone(&relation));
    Ok(relation)
}

fn eval_purge(args: &[ExprArc], ctx: &mut Context) -> Result<ExprArc, EvalError> {
    if args.is_empty() {
        return Err(EvalError::TooFewArgs("purge"));
    }
    let sym = eval(args[0].as_ref(), ctx)?;
    let name = match sym.as_ref() {
        Expr::Symbol(id) => id.clone(),
        _ => return Err(EvalError::TypeError("variable name expected")),
    };
    ctx.purge_relation_assumptions_for(&name);
    Ok(Expr::int(0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Expr;

    #[test]
    fn assume_purge_roundtrip() {
        use crate::expr::RelOp;
        let mut ctx = Context::new();
        let rel = Arc::new(Expr::Relation(RelOp::Gt, Expr::sym("t"), Expr::int(2)));
        let stmt = Stmt::ExprStmt(Expr::func(FuncKind::Assume, vec![rel]));
        exec_stmt(&stmt, &mut ctx).unwrap();
        assert_eq!(ctx.relation_assumptions.len(), 1);
        let purge = Stmt::ExprStmt(Expr::func(FuncKind::Purge, vec![Expr::sym("t")]));
        exec_stmt(&purge, &mut ctx).unwrap();
        assert!(ctx.relation_assumptions.is_empty());
    }

    #[test]
    fn assignment_updates_context() {
        let mut ctx = Context::new();
        let stmt = Stmt::Assign(Ident::new("x"), Expr::int(5));
        exec_stmt(&stmt, &mut ctx).unwrap();
        assert_eq!(ctx.get(&Ident::new("x")), Some(&Expr::int(5)));
    }
}
