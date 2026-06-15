use logos::Logos;

#[derive(Logos, Debug, PartialEq, Clone)]
#[logos(skip r"[ \t\r\n]+")]
pub enum Token<'a> {
    #[token(";")]
    Semi,

    #[token(":=")]
    Assign,

    #[token("+")]
    Plus,
    #[token("-")]
    Minus,
    #[token("*")]
    Star,
    #[token("/")]
    Slash,
    #[token("^")]
    Caret,
    #[token("**")]
    StarStar,

    #[token("mod", priority = 2)]
    #[token("%")]
    Mod,

    #[token("(")]
    LParen,
    #[token(")")]
    RParen,
    #[token("[")]
    LBracket,
    #[token("]")]
    RBracket,
    #[token("[[")]
    LMat,
    #[token("]]")]
    RMat,
    #[token(",")]
    Comma,

    #[token("==")]
    EqEq,
    #[token("!=")]
    Ne,
    #[token("<=")]
    Le,
    #[token(">=")]
    Ge,
    #[token("<")]
    Lt,
    #[token(">")]
    Gt,
    #[token("=")]
    Eq,

    #[regex(r"[0-9]+(\.[0-9]+)?", |lex| lex.slice().to_string())]
    Number(String),

    #[regex(r"[a-zA-Z_][a-zA-Z0-9_]*", |lex| lex.slice())]
    Ident(&'a str),

    Error,
}

pub struct Lexer<'a> {
    inner: logos::Lexer<'a, Token<'a>>,
}

impl<'a> Lexer<'a> {
    pub fn new(input: &'a str) -> Self {
        Self {
            inner: Token::lexer(input),
        }
    }

    pub fn next_token(&mut self) -> Option<Result<Token<'a>, ()>> {
        match self.inner.next() {
            None => None,
            Some(Ok(Token::Error)) => Some(Err(())),
            Some(Ok(tok)) => Some(Ok(tok)),
            Some(Err(())) => Some(Err(())),
        }
    }

    pub fn span(&self) -> std::ops::Range<usize> {
        self.inner.span()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lex_all(input: &str) -> Vec<Token<'_>> {
        let mut lexer = Lexer::new(input);
        let mut out = Vec::new();
        while let Some(tok) = lexer.next_token() {
            out.push(tok.unwrap());
        }
        out
    }

    #[test]
    fn lex_basic_expression() {
        let toks = lex_all("1+2*i;");
        assert_eq!(
            toks,
            vec![
                Token::Number("1".into()),
                Token::Plus,
                Token::Number("2".into()),
                Token::Star,
                Token::Ident("i"),
                Token::Semi,
            ]
        );
    }

    #[test]
    fn lex_matrix_literal() {
        let toks = lex_all("[[1,2],[3,4]]");
        eprintln!("{toks:?}");
        assert!(matches!(toks.first(), Some(Token::LMat)));
    }

    #[test]
    fn lex_assign_and_power() {
        let toks = lex_all("x:=a^2;");
        assert!(matches!(toks.as_slice(), [
            Token::Ident("x"),
            Token::Assign,
            Token::Ident("a"),
            Token::Caret,
            Token::Number(_),
            Token::Semi,
        ]));
    }

    #[test]
    fn lex_relation_and_mod_tokens() {
        let toks = lex_all("a==b!=c mod 3");
        assert!(toks.contains(&Token::EqEq));
        assert!(toks.contains(&Token::Ne));
        assert!(toks.contains(&Token::Mod));
    }

    #[test]
    fn lex_error_on_invalid_char() {
        let mut lexer = Lexer::new("@");
        assert!(matches!(lexer.next_token(), Some(Err(()))));
    }

    #[test]
    fn lexer_span_tracks_position() {
        let mut lexer = Lexer::new("abc");
        let _ = lexer.next_token();
        assert!(lexer.span().end > lexer.span().start);
    }
}
