#!/usr/bin/env bash
# Remove giac-2.0.0 compile outputs from the git index (for filter-branch --index-filter).
set -euo pipefail
while IFS= read -r f; do
  case "$f" in
    giac/giac-2.0.0/src/.libs/*|giac/giac-2.0.0/src/icas|giac/giac-2.0.0/src/*.o)
      git rm -rf --cached --ignore-unmatch "$f" >/dev/null 2>&1 || true
      ;;
  esac
done < <(git ls-files 'giac/giac-2.0.0/src')
