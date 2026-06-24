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

/// Parse a comma-separated Xcas script line (`assume(...),integrate(...),purge(...)`).
pub fn parse_compound_line(input: &str, ctx: &Context) -> Result<Vec<Stmt>, ParseError> {
    let trimmed = input.trim().trim_end_matches(';');
    let parts = split_top_level_commas(trimmed);
    let mut stmts = Vec::new();
    for part in parts {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let src = format!("{part};");
        let mut p = Parser::new(&src, ctx);
        stmts.push(p.parse_stmt()?);
    }
    Ok(stmts)
}

fn split_top_level_commas(input: &str) -> Vec<String> {
    let mut parts = Vec::new();
    let mut start = 0usize;
    let mut depth = 0i32;
    for (i, ch) in input.char_indices() {
        match ch {
            '(' | '[' => depth += 1,
            ')' | ']' => depth = depth.saturating_sub(1),
            ',' if depth == 0 => {
                parts.push(input[start..i].to_string());
                start = i + ch.len_utf8();
            }
            _ => {}
        }
    }
    parts.push(input[start..].to_string());
    parts
}

struct Parser<'a, 'ctx> {
    tokens: Vec<Token<'a>>,
    pos: usize,
    ctx: &'ctx Context,
}

impl<'a, 'ctx> Parser<'a, 'ctx> {
    fn new(input: &'a str, ctx: &'ctx Context) -> Self {
        let mut lexer = Lexer::new(input);
        let mut tokens = Vec::new();
        while let Some(tok) = lexer.next_token() {
            match tok {
                Ok(t) => tokens.push(t),
                Err(()) => tokens.push(Token::Error),
            }
        }
        Self {
            tokens,
            pos: 0,
            ctx,
        }
    }

    fn at_end(&self) -> bool {
        self.pos >= self.tokens.len()
    }

    fn bump(&mut self) -> Result<Token<'a>, ParseError> {
        if self.pos >= self.tokens.len() {
            return Err(ParseError::Eof);
        }
        let tok = self.tokens[self.pos].clone();
        self.pos += 1;
        if tok == Token::Error {
            return Err(ParseError::Lexer);
        }
        Ok(tok)
    }

    fn peek(&self) -> Option<&Token<'a>> {
        self.tokens.get(self.pos)
    }

    fn save(&self) -> usize {
        self.pos
    }

    fn restore(&mut self, pos: usize) {
        self.pos = pos;
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

    fn try_parse_lambda(&mut self) -> Result<ExprArc, ParseError> {
        let mut params = Vec::new();
        loop {
            match self.bump()? {
                Token::Ident(name) => params.push(Expr::sym(name)),
                _ => return Err(ParseError::Unexpected("lambda parameter")),
            }
            match self.peek() {
                Some(Token::Comma) => {
                    self.bump()?;
                }
                Some(Token::RParen) => {
                    self.bump()?;
                    break;
                }
                _ => return Err(ParseError::Unexpected("lambda parameter list")),
            }
        }
        self.expect(Token::Arrow)?;
        let body = self.parse_expr()?;
        Ok(Expr::func(
            FuncKind::Lambda,
            vec![Arc::new(Expr::List(params)), body],
        ))
    }

    /// Continue parsing infix/postfix after a leading identifier was consumed.
    fn finish_expr_from(&mut self, mut expr: ExprArc) -> Result<ExprArc, ParseError> {
        loop {
            match self.peek() {
                Some(Token::LParen) => {
                    self.bump()?;
                    let args = self.parse_arg_list()?;
                    self.expect(Token::RParen)?;
                    expr = self.finish_call(expr, args)?;
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
                Some(Token::Mod) => {
                    self.bump()?;
                    let rhs = self.parse_unary()?;
                    expr = Arc::new(Expr::Mod(expr, rhs));
                }
                Some(Token::Prime) => {
                    let mut order = 0i64;
                    while matches!(self.peek(), Some(Token::Prime)) {
                        self.bump()?;
                        order += 1;
                    }
                    expr = Expr::func(FuncKind::Prime, vec![expr, Expr::int(order)]);
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
                Some(Token::Ident(_) | Token::LParen | Token::Number(_))
                    if factors
                        .last()
                        .is_some_and(imp_mult_after) =>
                {
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
            let mut order = 0i64;
            while matches!(self.peek(), Some(Token::Prime)) {
                self.bump()?;
                order += 1;
            }
            if order > 0 {
                expr = Expr::func(FuncKind::Prime, vec![expr, Expr::int(order)]);
            }

            match self.peek() {
                Some(Token::LParen) => {
                    self.bump()?;
                    let args = self.parse_arg_list()?;
                    self.expect(Token::RParen)?;
                    expr = self.finish_call(expr, args)?;
                }
                Some(Token::LBracket) => {
                    self.bump()?;
                    if let Expr::Symbol(id) = expr.as_ref() {
                        if id.as_str() == "poly1" {
                            let items = self.parse_seq_items()?;
                            self.expect(Token::RBracket)?;
                            expr = Expr::func(FuncKind::Poly1, vec![Arc::new(Expr::Seq(items))]);
                            continue;
                        }
                    }
                    let idx = self.parse_expr()?;
                    self.expect(Token::RBracket)?;
                    expr = Expr::func(FuncKind::Sign, vec![expr, idx]);
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
                let saved = self.save();
                // `(` already consumed; try `(p,q)->body` before treating as grouped expr.
                match self.try_parse_lambda() {
                    Ok(lambda) => return Ok(lambda),
                    Err(_) => self.restore(saved),
                }
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
        loop {
            let row = if matches!(self.peek(), Some(Token::LBracket)) {
                self.bump()?;
                let items = self.parse_seq_items()?;
                match self.peek() {
                    Some(Token::RMat) => {
                        rows.push(items);
                        self.bump()?;
                        return Ok(Arc::new(Expr::Matrix(rows)));
                    }
                    Some(Token::RBracket) => {
                        self.bump()?;
                        items
                    }
                    _ => return Err(ParseError::Unexpected("expected ']' or ']]'")),
                }
            } else {
                let items = self.parse_seq_items()?;
                if matches!(self.peek(), Some(Token::RBracket)) {
                    self.bump()?;
                }
                items
            };
            rows.push(row);
            match self.peek() {
                Some(Token::Comma) => {
                    self.bump()?;
                }
                Some(Token::RMat) => {
                    self.bump()?;
                    break;
                }
                _ => return Err(ParseError::Unexpected("expected ',' or ']]'")),
            }
        }
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

    fn finish_call(&self, callee: ExprArc, args: Vec<ExprArc>) -> Result<ExprArc, ParseError> {
        if let Expr::Symbol(id) = callee.as_ref() {
            if let Some(kind) = lookup_func(id.as_str()) {
                return Ok(Expr::func(kind, args));
            }
        }
        let arg = if args.len() == 1 {
            args.into_iter()
                .next()
                .ok_or(ParseError::Unexpected("empty call"))?
        } else {
            Arc::new(Expr::Seq(args))
        };
        Ok(Expr::func(FuncKind::Apply, vec![callee, arg]))
    }
}

fn imp_mult_after(factor: &ExprArc) -> bool {
    matches!(
        factor.as_ref(),
        Expr::Int(_)
            | Expr::Rat(_)
            | Expr::Symbol(_)
            | Expr::Pow(_, _)
            | Expr::Func(_, _)
            | Expr::Add(_)
            | Expr::Mul(_)
            | Expr::Matrix(_)
            | Expr::GiacMatrix(_)
    )
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
        "atan" => Some(FuncKind::Atan),
        "tan" => Some(FuncKind::Tan),
        "exp" => Some(FuncKind::Exp),
        "ln" => Some(FuncKind::Ln),
        "re" => Some(FuncKind::Re),
        "im" => Some(FuncKind::Im),
        "arg" => Some(FuncKind::Arg),
        "sign" => Some(FuncKind::Sign),
        "normal" => Some(FuncKind::Normal),
        "ratnormal" => Some(FuncKind::Ratnormal),
        "expand" => Some(FuncKind::Expand),
        "texpand" => Some(FuncKind::Texpand),
        "tlin" => Some(FuncKind::Tlin),
        "halftan" => Some(FuncKind::Halftan),
        "lin" => Some(FuncKind::Lin),
        "factor" => Some(FuncKind::Factor),
        "ifactor" => Some(FuncKind::Ifactor),
        "quo" => Some(FuncKind::Quo),
        "rem" => Some(FuncKind::Rem),
        "content" => Some(FuncKind::Content),
        "gauss" => Some(FuncKind::Gauss),
        "egcd" => Some(FuncKind::Egcd),
        "abcuv" => Some(FuncKind::Abcuv),
        "simp2" => Some(FuncKind::Simp2),
        "lcm" => Some(FuncKind::Lcm),
        "horner" => Some(FuncKind::Horner),
        "resultant" => Some(FuncKind::Resultant),
        "roots" => Some(FuncKind::Roots),
        "modp" => Some(FuncKind::Modp),
        "smod" => Some(FuncKind::Smod),
        "irem" => Some(FuncKind::Irem),
        "chinrem" => Some(FuncKind::Chinrem),
        "partfrac" => Some(FuncKind::Partfrac),
        "greduce" => Some(FuncKind::Greduce),
        "rref" => Some(FuncKind::Rref),
        "integrate" | "int" => Some(FuncKind::Integrate),
        "diff" => Some(FuncKind::Diff),
        "derive" => Some(FuncKind::Derive),
        "solve" => Some(FuncKind::Solve),
        "fsolve" => Some(FuncKind::Fsolve),
        "sturm" => Some(FuncKind::Sturm),
        "sturmab" => Some(FuncKind::Sturmab),
        "realroot" => Some(FuncKind::Realroot),
        "limit" => Some(FuncKind::Limit),
        "series" => Some(FuncKind::Series),
        "taylor" => Some(FuncKind::Taylor),
        "desolve" => Some(FuncKind::Desolve),
        "risch" => Some(FuncKind::Risch),
        "proot" => Some(FuncKind::Proot),
        "simplify" => Some(FuncKind::Simplify),
        "idn" => Some(FuncKind::Idn),
        "inv" => Some(FuncKind::Inv),
        "det" => Some(FuncKind::Det),
        "tran" => Some(FuncKind::Tran),
        "ker" => Some(FuncKind::Ker),
        "image" => Some(FuncKind::Image),
        "pcar" => Some(FuncKind::Pcar),
        "charpoly" => Some(FuncKind::Charpoly),
        "linsolve" => Some(FuncKind::Linsolve),
        "jordan" => Some(FuncKind::Jordan),
        "egv" => Some(FuncKind::Egv),
        "lu" => Some(FuncKind::Lu),
        "qr" => Some(FuncKind::Qr),
        "svd" => Some(FuncKind::Svd),
        "gramschmidt" => Some(FuncKind::Gramschmidt),
        "trace" => Some(FuncKind::Trace),
        "subst" => Some(FuncKind::Subst),
        "assume" => Some(FuncKind::Assume),
        "purge" => Some(FuncKind::Purge),
        "froot" => Some(FuncKind::Froot),
        "froots" => Some(FuncKind::Froots),
        "rootof" => Some(FuncKind::RootOf),
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
    fn parse_matrix_inv() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("inv([[1,2],[3,4]]);", &ctx).unwrap();
        assert_eq!(stmts.len(), 1);
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

    #[test]
    fn parse_relations_and_modulo() {
        let ctx = Context::xcas_default();
        for src in ["1==2;", "1!=2;", "1<2;", "1>=2;", "7 mod 3;"] {
            let stmts = parse_program(src, &ctx).unwrap();
            assert_eq!(stmts.len(), 1, "failed on {src}");
        }
        let stmts = parse_program("7 mod 3;", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Mod(_, _))
        ));
        let rel = parse_program("(x==1);", &ctx).unwrap();
        assert!(matches!(
            &rel[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Relation(RelOp::Eq, _, _))
        ));
        let fm = parse_program("factor(x^4-1) mod 2;", &ctx).unwrap();
        assert!(matches!(
            &fm[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Mod(_, _))
        ));
    }

    #[test]
    fn parse_implicit_mult_and_decimal() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("2x; 1.5;", &ctx).unwrap();
        match &stmts[0] {
            Stmt::ExprStmt(e) => match e.as_ref() {
                Expr::Mul(f) => assert_eq!(f.len(), 2),
                _ => panic!("expected implicit mult"),
            },
            _ => panic!("expected expr"),
        }
        assert!(matches!(&stmts[1], Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Rat(_))));
    }

    #[test]
    fn parse_seq_poly1_and_star_star() {
        let ctx = Context::xcas_default();
        let seq = parse_program("[1,2,3];", &ctx).unwrap();
        assert!(matches!(
            &seq[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Seq(_))
        ));
        let poly = parse_program("a:=poly1[1,0];", &ctx).unwrap();
        assert!(matches!(
            &poly[0],
            Stmt::Assign(_, e) if matches!(e.as_ref(), Expr::Func(FuncKind::Poly1, _))
        ));
        let pow = parse_program("2**3;", &ctx).unwrap();
        let ev = giac_core::eval(
            match &pow[0] {
                Stmt::ExprStmt(e) => e,
                _ => panic!("expected expr"),
            },
            &ctx,
        )
        .unwrap();
        assert_eq!(ev, Expr::int(8));
    }

    #[test]
    fn parse_matrix_bracket_form() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("[[1,2],[3,4]];", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Matrix(_))
        ));
    }

    #[test]
    fn parse_unary_plus_and_finish_expr() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("f:=x/y+1;", &ctx).unwrap();
        match &stmts[0] {
            Stmt::Assign(_, e) => assert!(matches!(e.as_ref(), Expr::Add(_))),
            _ => panic!("expected assign"),
        }
        let stmts = parse_program("+3;", &ctx).unwrap();
        assert_eq!(
            giac_core::eval(
                match &stmts[0] {
                    Stmt::ExprStmt(e) => e,
                    _ => panic!("expected expr"),
                },
                &ctx,
            )
            .unwrap(),
            Expr::int(3)
        );
    }

    #[test]
    fn parse_unknown_function_becomes_apply() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("unknown(1);", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(
                e.as_ref(),
                Expr::Func(FuncKind::Apply, args)
                    if args.len() == 2
                        && matches!(args[0].as_ref(), Expr::Symbol(id) if id.as_str() == "unknown")
            )
        ));
        let stmts = parse_program("(1)(2);", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Func(FuncKind::Apply, _))
        ));
    }

    #[test]
    fn parse_lexer_error() {
        let ctx = Context::xcas_default();
        assert_eq!(parse_program("@;", &ctx), Err(ParseError::Lexer));
    }

    #[test]
    fn parse_finish_expr_from_ident_statements() {
        let ctx = Context::xcas_default();
        for (src, check) in [
            ("x^2;", "Pow"),
            ("x*y;", "Mul"),
            ("x/y;", "Mul"),
            ("x+y-1;", "Add"),
            ("gcd(12,18);", "Func"),
        ] {
            let stmts = parse_program(src, &ctx).unwrap();
            let e = match &stmts[0] {
                Stmt::ExprStmt(e) => e.as_ref(),
                other => panic!("{other:?}"),
            };
            match (check, e) {
                ("Pow", Expr::Pow(_, _)) => {}
                ("Mul", Expr::Mul(_)) => {}
                ("Add", Expr::Add(_)) => {}
                ("Func", Expr::Func { .. }) => {}
                _ => panic!("{src} -> {e:?}"),
            }
        }
    }

    #[test]
    fn parse_single_eq_relation_and_no_semi() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("1=2;", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Relation(RelOp::Eq, _, _))
        ));
        let stmts = parse_program("1+2", &ctx).unwrap();
        assert_eq!(stmts.len(), 1);
    }

    #[test]
    fn parse_implicit_mult_paren() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("2*(x+1);", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Mul(_))
        ));
    }

    #[test]
    fn parse_index_via_sign() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("a:=sign(x,2);", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::Assign(_, e) if matches!(e.as_ref(), Expr::Func(FuncKind::Sign, _))
        ));
    }

    #[test]
    fn parse_matrix_row_bracket_form() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("[[1],[2,3]];", &ctx).unwrap();
        match &stmts[0] {
            Stmt::ExprStmt(e) => match e.as_ref() {
                Expr::Matrix(rows) => {
                    assert_eq!(rows.len(), 2);
                    assert_eq!(rows[0].len(), 1);
                    assert_eq!(rows[1].len(), 2);
                }
                other => panic!("{other:?}"),
            },
            _ => panic!("expected expr"),
        }
    }

    #[test]
    fn parse_empty_call_and_identifiers() {
        let ctx = Context::xcas_default();
        let stmts = parse_program("gcd();", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Func(FuncKind::Gcd, args) if args.is_empty())
        ));
        let pi = parse_program("pi;", &ctx).unwrap();
        assert!(matches!(
            &pi[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Symbol(id) if id.as_str() == "pi")
        ));
        let ii = parse_program("a:=ii;", &ctx).unwrap();
        assert!(matches!(
            &ii[0],
            Stmt::Assign(_, e) if matches!(e.as_ref(), Expr::Symbol(id) if id.as_str() == "i")
        ));
    }

    #[test]
    fn parse_le_gt_and_invalid_number() {
        let ctx = Context::xcas_default();
        for src in ["1<=2;", "3>2;"] {
            parse_program(src, &ctx).expect(src);
        }
        let stmts = parse_program("9223372036854775808;", &ctx).unwrap();
        assert_eq!(
            giac_core::eval(
                match &stmts[0] {
                    Stmt::ExprStmt(e) => e,
                    _ => panic!("expected expr"),
                },
                &ctx,
            )
            .unwrap(),
            Expr::int(0)
        );
    }

    #[test]
    fn parse_eof_and_token_mismatch() {
        let ctx = Context::xcas_default();
        assert_eq!(parse_program("1+", &ctx), Err(ParseError::Eof));
        assert!(parse_program("1+);", &ctx).is_err());
    }

    #[test]
    fn parse_phase4_infinity_and_series() {
        let ctx = Context::xcas_default();
        parse_program("limit((1+1/x)^x,x,+infinity);", &ctx).unwrap();
        parse_program("series(exp(x),x,0,4);", &ctx).unwrap();
        let stmts = parse_program("taylor(sin(x),x=0,5);", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Func(FuncKind::Taylor, _))
        ));
    }

    #[test]
    fn parse_compound_assume_integrate_purge() {
        let ctx = Context::xcas_default();
        let stmts = parse_compound_line(
            "assume(t>2),integrate(x,x,2,t),purge(t)",
            &ctx,
        )
        .unwrap();
        assert_eq!(stmts.len(), 3);
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Func(FuncKind::Assume, _))
        ));
        assert!(matches!(
            &stmts[2],
            Stmt::ExprStmt(e) if matches!(e.as_ref(), Expr::Func(FuncKind::Purge, _))
        ));
    }

    #[test]
    fn parse_phase4_desolve_and_prime() {
        let ctx = Context::xcas_default();
        parse_program("desolve(y''+y=0,y(x));", &ctx).unwrap();
        parse_program("limit(sin(x)/x,x,0);", &ctx).unwrap();
        parse_program("taylor(sin(x),x=0,5);", &ctx).unwrap();
        let stmts = parse_program("y'+y;", &ctx).unwrap();
        assert!(matches!(
            &stmts[0],
            Stmt::ExprStmt(e) if matches!(
                e.as_ref(),
                Expr::Add(terms) if terms.len() == 2
                    && matches!(terms[0].as_ref(), Expr::Func(FuncKind::Prime, _))
            )
        ));
    }
}
