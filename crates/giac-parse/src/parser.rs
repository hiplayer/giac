use std::sync::Arc;

use giac_core::{Context, Expr, ExprArc, FuncKind, Ident, RelOp, Stmt};
use thiserror::Error;

use crate::lexer::{Lexer, Token};

#[derive(Debug, Error, PartialEq)]
pub enum ParseError {
    #[error("unexpected token: {0}")]
    Unexpected(&'static str),
    #[error("unexpected end of input")]
    Eof,
    #[error("lexer error")]
    Lexer,
}

pub fn parse_program(input: &str, ctx: &Context) -> Result<Vec<Stmt>, ParseError> {
    let mut p = Parser::new(input, ctx);
    let mut stmts = Vec::new();
    while !p.at_end() {
        stmts.push(p.parse_stmt()?);
    }
    Ok(stmts)
}

struct Parser<'a, 'ctx> {
    lexer: Lexer<'a>,
    current: Option<Result<Token<'a>, ()>>,
    ctx: &'ctx Context,
}

impl<'a, 'ctx> Parser<'a, 'ctx> {
    fn new(input: &'a str, ctx: &'ctx Context) -> Self {
        let mut lexer = Lexer::new(input);
        let current = lexer.next_token();
        Self {
            lexer,
            current,
            ctx,
        }
    }

    fn at_end(&self) -> bool {
        self.current.is_none()
    }

    fn bump(&mut self) -> Result<Token<'a>, ParseError> {
        match self.current.take() {
            Some(Ok(tok)) => {
                self.current = self.lexer.next_token();
                Ok(tok)
            }
            Some(Err(())) => Err(ParseError::Lexer),
            None => Err(ParseError::Eof),
        }
    }

    fn peek(&self) -> Option<&Token<'a>> {
        match &self.current {
            Some(Ok(tok)) => Some(tok),
            _ => None,
        }
    }

    fn expect(&mut self, want: Token<'a>) -> Result<(), ParseError> {
        match self.bump() {
            Ok(got) if std::mem::discriminant(&got) == std::mem::discriminant(&want) => Ok(()),
            Ok(_) => Err(ParseError::Unexpected("token mismatch")),
            Err(e) => Err(e),
        }
    }

    fn parse_stmt(&mut self) -> Result<Stmt, ParseError> {
        if let Some(Token::Ident(name)) = self.peek().cloned() {
            self.bump()?;
            if matches!(self.peek(), Some(Token::Assign)) {
                self.bump()?;
                let expr = self.parse_expr()?;
                self.expect_semi()?;
                return Ok(Stmt::Assign(Ident::new(name), expr));
            }
            let mut expr = Expr::sym(name);
            expr = self.finish_expr_from(expr)?;
            self.expect_semi()?;
            return Ok(Stmt::ExprStmt(expr));
        }
        let expr = self.parse_expr()?;
        self.expect_semi()?;
        Ok(Stmt::ExprStmt(expr))
    }

    /// Continue parsing infix/postfix after a leading identifier was consumed.
    fn finish_expr_from(&mut self, mut expr: ExprArc) -> Result<ExprArc, ParseError> {
        loop {
            match self.peek() {
                Some(Token::LParen) => {
                    self.bump()?;
                    let args = self.parse_arg_list()?;
                    self.expect(Token::RParen)?;
                    if let Expr::Symbol(id) = expr.as_ref() {
                        let kind = lookup_func(id.as_str())
                            .ok_or(ParseError::Unexpected("unknown function"))?;
                        expr = Expr::func(kind, args);
                    } else {
                        return Err(ParseError::Unexpected("call on non-function"));
                    }
                }
                Some(Token::Caret | Token::StarStar) => {
                    self.bump()?;
                    let exp = self.parse_unary()?;
                    expr = Expr::pow(expr, exp);
                }
                Some(Token::Star) => {
                    self.bump()?;
                    let rhs = self.parse_unary()?;
                    expr = Expr::mul(vec![expr, rhs]);
                }
                Some(Token::Slash) => {
                    self.bump()?;
                    let rhs = self.parse_unary()?;
                    expr = Expr::mul(vec![expr, Expr::pow(rhs, Expr::int(-1))]);
                }
                Some(Token::Plus | Token::Minus) => {
                    return self.parse_add_from(expr);
                }
                _ => break,
            }
        }
        Ok(expr)
    }

    fn parse_add_from(&mut self, first: ExprArc) -> Result<ExprArc, ParseError> {
        let mut terms = vec![first];
        while matches!(self.peek(), Some(Token::Plus | Token::Minus)) {
            let neg = matches!(self.peek(), Some(Token::Minus));
            self.bump()?;
            let mut term = self.parse_mul()?;
            if neg {
                term = Expr::mul(vec![Expr::int(-1), term]);
            }
            terms.push(term);
        }
        Ok(Expr::add(terms))
    }

    fn expect_semi(&mut self) -> Result<(), ParseError> {
        match self.peek() {
            Some(Token::Semi) => {
                self.bump()?;
                Ok(())
            }
            None => Ok(()),
            _ => Err(ParseError::Unexpected("expected ';'")),
        }
    }

    fn parse_expr(&mut self) -> Result<ExprArc, ParseError> {
        self.parse_relation()
    }

    fn parse_relation(&mut self) -> Result<ExprArc, ParseError> {
        let mut lhs = self.parse_add()?;
        while let Some(op) = self.peek_relation_op() {
            self.bump()?;
            let rhs = self.parse_add()?;
            lhs = Arc::new(Expr::Relation(op, lhs, rhs));
        }
        Ok(lhs)
    }

    fn peek_relation_op(&self) -> Option<RelOp> {
        match self.peek()? {
            Token::EqEq => Some(RelOp::Eq),
            Token::Ne => Some(RelOp::Ne),
            Token::Lt => Some(RelOp::Lt),
            Token::Le => Some(RelOp::Le),
            Token::Gt => Some(RelOp::Gt),
            Token::Ge => Some(RelOp::Ge),
            Token::Eq => Some(RelOp::Eq),
            _ => None,
        }
    }

    fn parse_add(&mut self) -> Result<ExprArc, ParseError> {
        let mut terms = vec![self.parse_mul()?];
        while matches!(self.peek(), Some(Token::Plus | Token::Minus)) {
            let neg = matches!(self.peek(), Some(Token::Minus));
            self.bump()?;
            let mut term = self.parse_mul()?;
            if neg {
                term = Expr::mul(vec![Expr::int(-1), term]);
            }
            terms.push(term);
        }
        Ok(Expr::add(terms))
    }

    fn parse_mul(&mut self) -> Result<ExprArc, ParseError> {
        let mut factors = vec![self.parse_pow()?];
        loop {
            match self.peek() {
                Some(Token::Star) => {
                    self.bump()?;
                    factors.push(self.parse_pow()?);
                }
                Some(Token::Slash) => {
                    self.bump()?;
                    let rhs = self.parse_pow()?;
                    factors.push(Expr::pow(rhs, Expr::int(-1)));
                }
                Some(Token::Mod) => {
                    self.bump()?;
                    let rhs = self.parse_pow()?;
                    let lhs = Expr::mul(factors);
                    return Ok(Arc::new(Expr::Mod(lhs, rhs)));
                }
                // implicit multiplication: number followed by ident or '('
                Some(Token::Ident(_) | Token::LParen) if is_numeric_factor(factors.last().unwrap())
                => {
                    factors.push(self.parse_pow()?);
                }
                _ => break,
            }
        }
        Ok(Expr::mul(factors))
    }

    fn parse_pow(&mut self) -> Result<ExprArc, ParseError> {
        let mut base = self.parse_unary()?;
        while matches!(self.peek(), Some(Token::Caret | Token::StarStar)) {
            self.bump()?;
            let exp = self.parse_unary()?;
            base = Expr::pow(base, exp);
        }
        Ok(base)
    }

    fn parse_unary(&mut self) -> Result<ExprArc, ParseError> {
        if matches!(self.peek(), Some(Token::Minus)) {
            self.bump()?;
            let inner = self.parse_unary()?;
            return Ok(Expr::mul(vec![Expr::int(-1), inner]));
        }
        if matches!(self.peek(), Some(Token::Plus)) {
            self.bump()?;
            return self.parse_unary();
        }
        self.parse_postfix()
    }

    fn parse_postfix(&mut self) -> Result<ExprArc, ParseError> {
        let mut expr = self.parse_atom()?;
        loop {
            match self.peek() {
                Some(Token::LParen) => {
                    self.bump()?;
                    let args = self.parse_arg_list()?;
                    self.expect(Token::RParen)?;
                    if let Expr::Symbol(id) = expr.as_ref() {
                        let kind = lookup_func(id.as_str()).ok_or(ParseError::Unexpected(
                            "unknown function",
                        ))?;
                        expr = Expr::func(kind, args);
                    } else {
                        return Err(ParseError::Unexpected("call on non-function"));
                    }
                }
                Some(Token::LBracket) => {
                    self.bump()?;
                    let idx = self.parse_expr()?;
                    self.expect(Token::RBracket)?;
                    expr = Expr::func(FuncKind::Sign, vec![expr, idx]); // placeholder: subscript later
                }
                _ => break,
            }
        }
        Ok(expr)
    }

    fn parse_atom(&mut self) -> Result<ExprArc, ParseError> {
        match self.bump()? {
            Token::Number(s) => Ok(parse_number(&s)),
            Token::Ident(name) => {
                if name == "pi" {
                    Ok(Expr::sym("pi"))
                } else if self.ctx.complex_mode && (name == "i" || name == "ii") {
                    Ok(Expr::sym("i"))
                } else {
                    Ok(Expr::sym(name))
                }
            }
            Token::LParen => {
                let e = self.parse_expr()?;
                self.expect(Token::RParen)?;
                Ok(e)
            }
            Token::LMat => self.parse_matrix(),
            Token::LBracket => {
                let items = self.parse_seq_items()?;
                self.expect(Token::RBracket)?;
                Ok(Arc::new(Expr::Seq(items)))
            }
            other => Err(ParseError::Unexpected(match other {
                Token::Semi => "';'",
                Token::Assign => "':='",
                _ => "atom",
            })),
        }
    }

    fn parse_matrix(&mut self) -> Result<ExprArc, ParseError> {
        let mut rows = Vec::new();
        if !matches!(self.peek(), Some(Token::RMat)) {
            loop {
                rows.push(self.parse_seq_items()?);
                if matches!(self.peek(), Some(Token::Comma)) {
                    self.bump()?;
                } else {
                    break;
                }
            }
        }
        self.expect(Token::RMat)?;
        Ok(Arc::new(Expr::Matrix(rows)))
    }

    fn parse_arg_list(&mut self) -> Result<Vec<ExprArc>, ParseError> {
        if matches!(self.peek(), Some(Token::RParen)) {
            return Ok(Vec::new());
        }
        let mut args = vec![self.parse_expr()?];
        while matches!(self.peek(), Some(Token::Comma)) {
            self.bump()?;
            args.push(self.parse_expr()?);
        }
        Ok(args)
    }

    fn parse_seq_items(&mut self) -> Result<Vec<ExprArc>, ParseError> {
        let mut items = vec![self.parse_expr()?];
        while matches!(self.peek(), Some(Token::Comma)) {
            self.bump()?;
            items.push(self.parse_expr()?);
        }
        Ok(items)
    }
}

fn is_numeric_factor(expr: &ExprArc) -> bool {
    matches!(expr.as_ref(), Expr::Int(_) | Expr::Rat(_))
}

fn parse_number(s: &str) -> ExprArc {
    if let Some((a, b)) = s.split_once('.') {
        if b.parse::<u64>().is_ok() && a.parse::<i64>().is_ok() {
            // keep as rational approximation for now — eval path uses integers primarily
            let num: i64 = format!("{}{}", a, b).parse().unwrap_or(0);
            let den = 10_i64.pow(b.len() as u32);
            return Expr::rat(num, den);
        }
    }
    if let Ok(n) = s.parse::<i64>() {
        Expr::int(n)
    } else {
        Expr::int(0)
    }
}

fn lookup_func(name: &str) -> Option<FuncKind> {
    match name {
        "abs" => Some(FuncKind::Abs),
        "gcd" => Some(FuncKind::Gcd),
        "conj" => Some(FuncKind::Conj),
        "sqrt" => Some(FuncKind::Sqrt),
        "sin" => Some(FuncKind::Sin),
        "cos" => Some(FuncKind::Cos),
        "exp" => Some(FuncKind::Exp),
        "ln" => Some(FuncKind::Ln),
        "re" => Some(FuncKind::Re),
        "im" => Some(FuncKind::Im),
        "arg" => Some(FuncKind::Arg),
        "sign" => Some(FuncKind::Sign),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use giac_core::Context;

    use super::*;

    #[test]
    fn parse_add_mul_power() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("1+2*i^2;", &ctx).unwrap();
        assert_eq!(stmts.len(), 1);
    }

    #[test]
    fn parse_function_call() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("gcd(12,18);", &ctx).unwrap();
        match &stmts[0] {
            Stmt::ExprStmt(e) => match e.as_ref() {
                Expr::Func(FuncKind::Gcd, args) => assert_eq!(args.len(), 2),
                other => panic!("unexpected {other:?}"),
            },
            other => panic!("unexpected {other:?}"),
        }
    }

    #[test]
    fn parse_assignment() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("x:=1+2;", &ctx).unwrap();
        match &stmts[0] {
            Stmt::Assign(id, _) => assert_eq!(id.as_str(), "x"),
            _ => panic!("expected assign"),
        }
    }
}
