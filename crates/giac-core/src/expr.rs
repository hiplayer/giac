use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Zero;

use crate::algebra::alg_ext::AlgExtData;
use crate::ident::Ident;

pub type ExprArc = Arc<Expr>;

/// Top-level symbolic expression (user-visible CAS value).
#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    Int(BigInt),
    Rat(Ratio<BigInt>),
    Frac(ExprArc, ExprArc),
    Mod(ExprArc, ExprArc),
    Complex(ExprArc, ExprArc),
    Symbol(Ident),
    Add(Vec<ExprArc>),
    Mul(Vec<ExprArc>),
    Pow(ExprArc, ExprArc),
    Func(FuncKind, Vec<ExprArc>),
    Seq(Vec<ExprArc>),
    List(Vec<ExprArc>),
    /// `[[...]]` literal / computational result
    Matrix(Vec<Vec<ExprArc>>),
    /// `matrix[[...]]` display (idn, tran, …)
    GiacMatrix(Vec<Vec<ExprArc>>),
    Relation(RelOp, ExprArc, ExprArc),
    /// Algebraic extension element (upstream `_EXT` / `rootof` value).
    AlgExt(Arc<AlgExtData>),
    Str(String),
    Undefined,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RelOp {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FuncKind {
    Abs,
    Gcd,
    Conj,
    Sqrt,
    Sin,
    Cos,
    Atan,
    Exp,
    Ln,
    Re,
    Im,
    Arg,
    Sign,
    Normal,
    Ratnormal,
    Expand,
    Factor,
    Quo,
    Rem,
    Content,
    Gauss,
    Egcd,
    Abcuv,
    Simp2,
    Lcm,
    Horner,
    Resultant,
    Roots,
    Modp,
    Smod,
    Irem,
    Chinrem,
    Partfrac,
    Greduce,
    Rref,
    Integrate,
    Int,
    Idn,
    Inv,
    Det,
    Tran,
    Ker,
    Image,
    Pcar,
    Charpoly,
    Linsolve,
    Jordan,
    Egv,
    Lu,
    Qr,
    Svd,
    Gramschmidt,
    Trace,
    Lambda,
    Subst,
    RootOf,
    Poly1,
    Tan,
    Texpand,
    Tlin,
    Halftan,
    Lin,
    Diff,
    Derive,
    Solve,
    Fsolve,
    Sturm,
    Sturmab,
    Realroot,
    Limit,
    Series,
    Taylor,
    Desolve,
    Risch,
    Proot,
    Simplify,
    Ifactor,
    Assume,
    Purge,
    Froot,
    Froots,
    /// Unknown identifier applied to arguments, e.g. `y(x)` in ODEs.
    Apply,
    /// Prime notation: `y'`, `y''` → args `[base, order]`.
    Prime,
}

impl Expr {
    pub fn int(n: i64) -> ExprArc {
        Arc::new(Expr::Int(BigInt::from(n)))
    }

    pub fn rat(n: i64, d: i64) -> ExprArc {
        Arc::new(Expr::Rat(Ratio::new(
            BigInt::from(n),
            BigInt::from(d),
        )))
    }

    pub fn sym(name: &str) -> ExprArc {
        Arc::new(Expr::Symbol(Ident::new(name)))
    }

    pub fn add(terms: Vec<ExprArc>) -> ExprArc {
        if terms.len() == 1 {
            return Arc::clone(&terms[0]);
        }
        Arc::new(Expr::Add(terms))
    }

    pub fn mul(factors: Vec<ExprArc>) -> ExprArc {
        if factors.len() == 1 {
            return Arc::clone(&factors[0]);
        }
        Arc::new(Expr::Mul(factors))
    }

    pub fn pow(base: ExprArc, exp: ExprArc) -> ExprArc {
        Arc::new(Expr::Pow(base, exp))
    }

    pub fn func(kind: FuncKind, args: Vec<ExprArc>) -> ExprArc {
        Arc::new(Expr::Func(kind, args))
    }

    pub fn is_zero(&self) -> bool {
        match self {
            Expr::Int(n) => n.is_zero(),
            Expr::Rat(r) => r.is_zero(),
            Expr::Complex(re, im) => re.is_zero() && im.is_zero(),
            Expr::AlgExt(a) => a.is_zero(),
            _ => false,
        }
    }

    pub fn is_one(&self) -> bool {
        matches!(
            self,
            Expr::Int(n) if n == &BigInt::from(1)
        ) || matches!(
            self,
            Expr::Rat(r) if *r == Ratio::from_integer(BigInt::from(1))
        ) || matches!(self, Expr::AlgExt(a) if a.is_one())
    }

    pub fn alg_ext(data: AlgExtData) -> ExprArc {
        Arc::new(Expr::AlgExt(Arc::new(data)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructors_and_predicates() {
        assert!(Expr::int(0).is_zero());
        assert!(!Expr::int(1).is_zero());
        assert!(Expr::int(1).is_one());
        assert_eq!(Expr::add(vec![Expr::int(1)]), Expr::int(1));
        assert_eq!(Expr::mul(vec![Expr::sym("x")]), Expr::sym("x"));
    }

    #[test]
    fn rat_and_complex_predicates() {
        assert!(Expr::rat(0, 1).is_zero());
        assert!(!Expr::rat(3, 2).is_one());
        assert!(Expr::rat(2, 2).is_one());
        let zero_c = Expr::Complex(Expr::int(0), Expr::int(0));
        assert!(zero_c.is_zero());
        assert!(!Expr::Complex(Expr::int(1), Expr::int(0)).is_zero());
    }
}
