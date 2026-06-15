# giac-wasm

WebAssembly bindings for giac-rs: evaluate giac-style expressions in the browser or Node via [`wasm-bindgen`](https://rustwasm.github.io/wasm-bindgen/).

## Prerequisites

```bash
rustup target add wasm32-unknown-unknown
```

Optional (for bundling): [wasm-pack](https://rustwasm.github.io/wasm-pack/installer/) and a JavaScript bundler.

## Build

From the `giac-rs` workspace root:

```bash
cargo build --release --target wasm32-unknown-unknown -p giac-wasm
```

The artifact is `target/wasm32-unknown-unknown/release/giac_wasm.wasm` (with `cdylib` crate type).

### wasm-pack (browser / npm)

```bash
cd crates/giac-wasm
wasm-pack build --target web --release
```

This produces a `pkg/` directory with `giac_wasm.js` and `.wasm`. Call the exported function:

```javascript
import init, { evalToString } from "./pkg/giac_wasm.js";
await init();
console.log(evalToString("1+2")); // "3"
```

## API

| Rust | WASM (JS) | Description |
|------|-----------|-------------|
| `eval_to_string(input) -> Result<String, String>` | `evalToString(input) -> string` | Parse and evaluate; JS returns `"error: …"` on failure |

Uses `giac_ode::xcas_default()` (linear algebra, solve, calculus, and ODE plugin stubs).

## Test

Host unit tests (including `1+2` → `3`):

```bash
cargo test -p giac-wasm
```

Run the same check in the browser after `wasm-pack build` (see above).

## Notes

- Pure Rust dependency chain — no C FFI in giac-wasm or its transitive deps.
- `giac-cli` is excluded on `wasm32` (`cfg(not(target_arch = "wasm32"))`); use this crate instead.
