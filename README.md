# giac-rs

Rust reimplementation of the [Giac](https://www-fourier.univ-grenoble-alpes.fr/~parisse/giac/) computer algebra system (headless CAS subset).

## Status

**Phase 0–1** (in progress):

- `giac-core` + `giac-simplify` (facade): `expand`, `normal`, `ratnormal`, `factor`
- Extended builtins: `arg`, `re`, `im`, `sign`, `gcd` (int + poly), `integrate`, `idn`, `inv`, `det`, `tran`, `ker`, `image`, `pcar`, `subst`
- `giac-parse`: implicit multiply, `poly1[...]`, `rootof`, matrix `[[...]]` grammar
- Conformance: `test_cas_basic`, `cas.tst` lines 1–25 (21/25 golden), lines 1–20 (≥15/25)

## Build, test & lint

Requires **Rust 1.75+** (dependencies pinned for older toolchains).

**Clippy** is **required** before merge (see [supplement §7](../.doc/rust-migration-supplement.md#7-工程门禁)). Install: `rustup component add clippy`, or on Ubuntu/Debian system rustc: `sudo apt install rust-clippy` (must match `rustc --version`; do **not** `cargo install clippy`).

```bash
cd giac-rs
cargo test --workspace      # required
cargo ci-clippy             # required (-D warnings)
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

See [`.doc/rust-migration-plan.md`](../.doc/rust-migration-plan.md) for the full roadmap. Error-handling rules: [supplement §6](../.doc/rust-migration-supplement.md#6-错误处理与-unwrap-规范).
