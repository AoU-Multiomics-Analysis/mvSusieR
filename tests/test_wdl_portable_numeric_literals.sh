#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

scientific_literals="$({
  rg --pcre2 -n --glob '*.wdl' \
    '(?<![A-Za-z0-9_])[0-9]+(?:\.[0-9]+)?[eE][+-]?[0-9]+' \
    workflows || true
})"

if [[ -n "$scientific_literals" ]]; then
  echo "WDL files contain scientific-notation literals that Dockstore cannot parse:" >&2
  echo "$scientific_literals" >&2
  exit 1
fi

echo "WDL numeric literals are portable"
