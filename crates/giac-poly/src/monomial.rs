//! **API inventory:** inline `/// **Tier**` / `// **Tier**` on every function;
//! full module index in `.doc/giac-poly-api-stability.md`.
//!
//!
use std::collections::BTreeMap;
use std::sync::Arc;

/// Variable name in a multivariate polynomial.
pub type Var = Arc<str>;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Monomial(BTreeMap<Var, u64>);

impl Monomial {
    /// **Stable** — Poly one
    pub fn one() -> Self {
        Self(BTreeMap::new())
    }

    /// **Stable** — Poly univariate generator
    pub fn var(name: impl Into<Var>) -> Self {
        let mut m = BTreeMap::new();
        m.insert(name.into(), 1);
        Self(m)
    }

    /// **Stable** — total degree
    pub fn degree(&self) -> u64 {
        self.0.values().sum()
    }

    /// **Stable** — `Poly::exp_of`
    pub fn exp_of(&self, var: &Var) -> u64 {
        self.0.get(var).copied().unwrap_or(0)
    }

    /// **Stable** — `is_const`
    pub fn is_const(&self) -> bool {
        self.0.is_empty()
    }

    /// **Stable** — `iter`
    pub fn iter(&self) -> impl Iterator<Item = (&Var, u64)> + '_ {
        self.0.iter().map(|(v, e)| (v, *e))
    }

    /// **Stable** — Poly multiplication
    pub fn mul(&self, other: &Self) -> Self {
        let mut out = self.0.clone();
        for (v, e) in &other.0 {
            *out.entry(v.clone()).or_insert(0) += e;
        }
        Self(out)
    }

    /// **Stable** — exact division if remainder zero
    pub fn div_exact(&self, other: &Self) -> Option<Self> {
        let mut out = self.0.clone();
        for (v, e) in &other.0 {
            let entry = out.get_mut(v)?;
            if *entry < *e {
                return None;
            }
            *entry -= e;
            if *entry == 0 {
                out.remove(v);
            }
        }
        Some(Self(out))
    }

    /// **Stable** — `is_dividing`
    pub fn is_dividing(&self, other: &Self) -> bool {
        other.0.iter().all(|(v, e)| self.0.get(v).copied().unwrap_or(0) >= *e)
    }

    /// Lex comparison: `order[0]` is the most significant variable.
    /// **Stable** — `cmp_lex`
    pub fn cmp_lex(&self, other: &Self, order: &[Var]) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        for v in order {
            match self.exp_of(v).cmp(&other.exp_of(v)) {
                Ordering::Equal => {}
                o => return o,
            }
        }
        Ordering::Equal
    }

    /// True if `self` divides `other` (i.e. `other/self` is a monomial).
    /// **Stable** — `divides`
    pub fn divides(&self, other: &Self) -> bool {
        self.0
            .iter()
            .all(|(v, e)| other.0.get(v).copied().unwrap_or(0) >= *e)
    }

    /// **Stable** — Poly lcm
    pub fn lcm(&self, other: &Self) -> Self {
        let mut out = self.0.clone();
        for (v, e) in &other.0 {
            let cur = out.entry(v.clone()).or_insert(0);
            *cur = (*cur).max(*e);
        }
        Self(out)
    }
}
