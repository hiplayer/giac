#!/usr/bin/env bash
# Workspace clippy + substring-golden lint (conformance-testing.md §3).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
cargo clippy --workspace --all-targets -- -D warnings
./scripts/lint-substring-golden.sh
