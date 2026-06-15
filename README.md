# giac-rs

Rust reimplementation of the [Giac](https://www-fourier.univ-grenoble-alpes.fr/~parisse/giac/) computer algebra system (headless CAS subset).

## Status

**Phase 0** (skeleton): workspace, `giac-core`, `giac-parse`, `giac-cli`, conformance harness.

`bin/test_cas_basic` (3 lines) passes golden equivalence.

## Build & test

Requires **Rust 1.75+** (dependencies pinned for older toolchains).

```bash
cd giac-rs
cargo test --workspace
cargo run -p giac-cli -- ../bin/test_cas_basic
```

Expected output: `sqrt(5),15,-3-4*i`

## Layout

```
giac-rs/
├── crates/giac-core/     # Expr, Context, eval, display
├── crates/giac-parse/    # lexer + recursive-descent parser
├── apps/giac-cli/        # script runner
└── tests/conformance/    # golden tests vs upstream bin/
```

See [`.doc/rust-migration-plan.md`](../.doc/rust-migration-plan.md) for the full roadmap.
