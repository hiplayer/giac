#!/usr/bin/env bash
# Run all workspace tests with a per-test timeout (default 10s).
#
# Preferred: cargo-nextest (fast, reads .config/nextest.toml).
# Fallback: GNU `timeout` per test (slow; no extra install).
#
# Usage:
#   ./scripts/test-with-timeout.sh [cargo nextest args...]
#   TEST_TIMEOUT_SECS=10 ./scripts/test-with-timeout.sh
#
# On timeout, re-run the failing test with:
#   RUST_BACKTRACE=1 cargo test -p PKG 'TEST_FILTER' -- --exact --nocapture

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TIMEOUT_SECS="${TEST_TIMEOUT_SECS:-10}"
export NEXTEST_PROFILE="${NEXTEST_PROFILE:-default}"

if command -v cargo-nextest >/dev/null 2>&1 || cargo nextest --version >/dev/null 2>&1; then
  echo "==> cargo nextest --release (per-test timeout: ${TIMEOUT_SECS}s via .config/nextest.toml)"
  if [[ -n "${TEST_PACKAGES:-}" ]]; then
    args=(--release)
    for pkg in ${TEST_PACKAGES}; do
      args+=(-p "$pkg")
    done
    exec cargo nextest run "${args[@]}" "$@"
  fi
  exec cargo nextest run --workspace --release "$@"
fi

echo "warning: cargo-nextest not found; using GNU timeout fallback (${TIMEOUT_SECS}s per test, slower)." >&2
echo "  Install (rustc 1.75): cargo install cargo-nextest --locked --version 0.9.85" >&2
echo >&2

if ! command -v timeout >/dev/null 2>&1; then
  echo "error: need cargo-nextest or GNU coreutils 'timeout'" >&2
  exit 1
fi

echo "==> building test binaries (release, no run)"
if [[ -n "${TEST_PACKAGES:-}" ]]; then
  for pkg in ${TEST_PACKAGES}; do
    cargo test --release -p "$pkg" --no-run >/dev/null
  done
else
  cargo test --workspace --release --no-run >/dev/null
fi

mapfile -t PACKAGES < <(
  cargo metadata --format-version=1 --no-deps 2>/dev/null \
    | python3 -c "
import json, sys
data = json.load(sys.stdin)
ws = set(data.get('workspace_members', []))
names = [p['name'] for p in data.get('packages', []) if p['id'] in ws]
filter_env = __import__('os').environ.get('TEST_PACKAGES', '').split()
if filter_env:
    names = [n for n in names if n in filter_env]
for n in sorted(names):
    print(n)
"
)

if [[ ${#PACKAGES[@]} -eq 0 ]]; then
  echo "error: no packages to test (check TEST_PACKAGES filter)" >&2
  exit 1
fi

TIMEOUT_HITS=()
FAILURES=0
RAN=0

run_one() {
  local pkg="$1"
  local kind="$2"   # test | doc
  local test_id="$3"
  local -a cmd

  RAN=$((RAN + 1))
  if [[ "$kind" == doc ]]; then
    cmd=(cargo test --release -p "$pkg" --doc "$test_id" -- --exact --nocapture)
  else
    cmd=(cargo test --release -p "$pkg" "$test_id" -- --exact --nocapture)
  fi

  printf '  [%ss] %s (%s) %s\n' "$TIMEOUT_SECS" "$pkg" "$kind" "$test_id"
  set +e
  timeout --kill-after=2s "${TIMEOUT_SECS}s" "${cmd[@]}" >/dev/null 2>&1
  local ec=$?
  set -e

  if [[ $ec -eq 0 ]]; then
    return 0
  fi
  if [[ $ec -eq 124 || $ec -eq 137 ]]; then
    echo "  *** TIMEOUT after ${TIMEOUT_SECS}s: ${pkg} :: ${test_id} ***" >&2
    if [[ "$kind" == doc ]]; then
      echo "  Debug: RUST_BACKTRACE=1 cargo test -p ${pkg} --doc '${test_id}' -- --exact --nocapture" >&2
    else
      echo "  Debug: RUST_BACKTRACE=1 cargo test -p ${pkg} '${test_id}' -- --exact --nocapture" >&2
    fi
    TIMEOUT_HITS+=("${pkg} :: ${test_id}")
    return 1
  fi
  echo "  *** FAILED (exit ${ec}): ${pkg} :: ${test_id} ***" >&2
  if [[ "$kind" == doc ]]; then
    echo "  Debug: RUST_BACKTRACE=1 cargo test -p ${pkg} --doc '${test_id}' -- --exact --nocapture" >&2
  else
    echo "  Debug: RUST_BACKTRACE=1 cargo test -p ${pkg} '${test_id}' -- --exact --nocapture" >&2
  fi
  FAILURES=$((FAILURES + 1))
  return 1
}

list_tests() {
  local pkg="$1"
  cargo test -p "$pkg" -- --list 2>/dev/null | awk '/: test$/ {print $1}'
}

list_doc_tests() {
  local pkg="$1"
  cargo test -p "$pkg" --doc -- --list 2>/dev/null | awk '/: test$/ {print $1}'
}

for pkg in "${PACKAGES[@]}"; do
  echo "=== ${pkg} ==="
  while IFS= read -r t; do
    [[ -z "$t" ]] && continue
    run_one "$pkg" test "$t" || true
  done < <(list_tests "$pkg" || true)

  while IFS= read -r t; do
    [[ -z "$t" ]] && continue
    run_one "$pkg" doc "$t" || true
  done < <(list_doc_tests "$pkg" || true)
done

echo
echo "Ran ${RAN} test(s)."

if ((${#TIMEOUT_HITS[@]} > 0)); then
  echo "Timed out (${#TIMEOUT_HITS[@]}):"
  for hit in "${TIMEOUT_HITS[@]}"; do
    echo "  - ${hit}"
  done
fi

if ((${#TIMEOUT_HITS[@]} > 0 || FAILURES > 0)); then
  exit 1
fi

echo "All tests passed within ${TIMEOUT_SECS}s each."
