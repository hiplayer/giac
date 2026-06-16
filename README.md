# giac-rs

Rust reimplementation of the [Giac](https://www-fourier.univ-grenoble-alpes.fr/~parisse/giac/) computer algebra system (headless CAS subset).

## Status

**Phase 0–3** complete for the current MVP scope; **Phase 4** (solve / calculus) in progress — engineering scaffolding landed (see below).

| Phase | Scope | Crates |
|-------|--------|--------|
| **0–1** | Expr, parse, eval, simplify, basic CAS | `giac-core`, `giac-parse`, `giac-simplify` |
| **2** | Polynomial ring, gcd/factor, modular, Groebner `greduce` | `giac-poly`, `giac-groebner` |
| **3** | Symbolic + numeric linear algebra | `giac-linalg` |
| **4** (WIP) | Solve, calculus, ODE | `giac-solve`, `giac-calculus`, `giac-ode` |

### Phase 3 APIs (via `giac-cli` / eval)

**Symbolic** (`giac-linalg`): `rref`, `det`, `inv`, `tran`, `ker`, `image`, `pcar`, `charpoly`, `linsolve`, `trace`, `gauss`, `egv`, `jordan`, `gramschmidt`

**Numeric** (`nalgebra` backend): `lu`, `qr`, `svd`

Matrix builtins require the linalg plugin — CLI and conformance call `giac_linalg::install_linalg` automatically. For embedded use:

```rust
use giac_core::Context;
use giac_linalg::install_linalg;

let mut ctx = Context::xcas_default();
install_linalg(&mut ctx);
```

### Earlier phases (summary)

- **Algebra:** `expand`, `normal`, `ratnormal`, `factor`, `gcd`, `integrate` (rule subset), `subst`
- **Poly / modular:** `quo`, `rem`, `content`, `egcd`, `abcuv`, `roots`, `chinrem`, `modp`, `smod`, `irem`, `greduce`, …
- **Complex / misc:** `arg`, `re`, `im`, `sign`, `abs`, `sqrt`, …

Known output differences vs upstream Giac are tracked in [known-divergences.md](../.doc/known-divergences.md) (DIV-070–073 for Phase 3).

## Build, test & lint

Requires **Rust 1.75+** (dependencies pinned for older toolchains).

**Clippy** is required before merge (see [supplement §7](../.doc/rust-migration-supplement.md#7-工程门禁)). Install: `rustup component add clippy`, or on Ubuntu/Debian system rustc: `sudo apt install rust-clippy` (must match `rustc --version`; do **not** `cargo install clippy`).

```bash
cd giac-rs
cargo test --workspace      # required
cargo ci-clippy             # required (-D warnings)
echo 'sqrt(5)' | cargo run -q --bin giac-cli
```

### WebAssembly (`giac-wasm`)

See [crates/giac-wasm/README.md](crates/giac-wasm/README.md). Minimal build:

```bash
rustup target add wasm32-unknown-unknown
cargo build --release --target wasm32-unknown-unknown -p giac-wasm
cargo test -p giac-wasm    # host test: eval_to_string("1+2") == "3"
```

Optional coverage (see [supplement §7](../.doc/rust-migration-supplement.md#7-工程门禁)):

```bash
cargo install cargo-tarpaulin   # once
cargo tarpaulin -p giac-core --out Stdout    # target ≥70%
cargo tarpaulin -p giac-linalg --out Stdout  # Phase 3 linalg
```

Triple validation (giac-rs + upstream Giac + SymPy) needs a **2.0** reference binary. By default conformance uses `../build-2.0/bin/giac` when present, else `../build/bin/giac`. Override with `GIAC_BINARY`. Golden scripts come from `giac/giac-2.0.0/check/` (`GIAC_VERSION_DIR` in `tests/conformance/src/lib.rs`).

## Conformance tests

Integration tests live in `tests/conformance/`. Run subsets with:

```bash
# Phase 3 — linear algebra
cargo test -p giac-conformance --test phase3_linalg
cargo test -p giac-conformance --test phase3_triple

# Phase 4 — solve / calculus (scaffolding + integral table)
cargo test -p giac-conformance --test phase4_parse
cargo test -p giac-conformance --test phase4_integrate_table
cargo test -p giac-conformance --test phase4_triple
cargo test -p giac-conformance --test test_diff
cargo test -p giac-conformance --test test_solve

# Phase 2 — polynomials / Groebner
cargo test -p giac-conformance --test phase2_poly
cargo test -p giac-conformance --test phase2_triple

# Cross-phase SymPy smoke (bin scripts)
cargo test -p giac-conformance --test sympy_all

# Phase 0–1 smoke + golden
cargo test -p giac-conformance --test smoke
cargo test -p giac-conformance --test cas_first_50
cargo test -p giac-conformance --test giac_check_factor
```

Upstream script coverage: `bin/test_linalg`, `test_linalg_ext`, `test_linalg_decomp`, `test_gauss_ext`, `test_poly`, `test_factor`, `test_groebner`, `test_cas_basic`, …

## Layout

```
giac-rs/
├── crates/
│   ├── giac-core/        # Expr, Context, eval, display, LinalgPlugin trait
│   ├── giac-parse/       # lexer + recursive-descent parser
│   ├── giac-simplify/    # expand / normal / ratnormal / factor facade
│   ├── giac-poly/        # polynomial ring, gcd, factor, modular
│   ├── giac-groebner/    # greduce
│   ├── giac-linalg/      # symbolic + numeric linear algebra
│   ├── giac-calculus/    # integrate, diff (CalculusPlugin)
│   ├── giac-ode/         # desolve (OdePlugin)
│   ├── giac-solve/       # solve, linsolve (SolvePlugin)
│   └── giac-wasm/        # WASM eval_to_string export
├── apps/giac-cli/        # script runner (native only; not wasm32)
└── tests/conformance/    # golden + triple checks vs upstream bin/
```

## Documentation

- [Rust migration plan](../.doc/rust-migration-plan.md) — crate graph, phases, type design
- [Known divergences](../.doc/known-divergences.md) — Rust vs Giac golden diffs and verification strategy
- [Phase 4 issues](../.doc/phase4-issues.md) — GIAC-201+ backlog and integral table
- Error-handling rules: [supplement §6](../.doc/rust-migration-supplement.md#6-错误处理与-unwrap-规范)
