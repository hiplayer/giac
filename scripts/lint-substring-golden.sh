#!/usr/bin/env bash
# Reject substring golden asserts (assert! / assert_eq! + .contains).
# Prefer assert_equiv or exact assert_eq!(format_expr(...), "...") — conformance-testing.md §3.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

violations="$(rg 'assert(_eq)?!\([^;]*\.contains\([\"'\'']' crates tests apps \
  --glob '*.rs' \
  --glob '!**/target/**' \
  || true)"

if [[ -n "$violations" ]]; then
  echo "lint-substring-golden: substring .contains() in assertions is disallowed." >&2
  echo "$violations" >&2
  echo "Use assert_equiv or exact assert_eq!; see .doc/conformance-testing.md §3." >&2
  exit 1
fi
