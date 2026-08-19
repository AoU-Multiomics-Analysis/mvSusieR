#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

dockerfile="envs/trans-window-mvsusie.Dockerfile"
workflow=".github/workflows/trans-window-mvsusie-image.yml"

test -s "$dockerfile"
test -s "$workflow"
rg -q '^FROM rocker/r-ver:4[.]4[.]1$' "$dockerfile"
rg -q 'stephenslab/mvsusieR' "$dockerfile"
rg -q 'stephenslab/susieR' "$dockerfile"
rg -q 'mvsusieR.*0[.]3[.]0' "$workflow"
rg -q 'susieR.*0[.]15[.]54' "$workflow"

for script in \
  trans_window_io.R \
  trans_window_preprocess.R \
  trans_window_model.R \
  trans_window_cli.R \
  run_window_mvsusie.R \
  summarize_window.R \
  merge_window_outputs.R; do
  rg -q "scripts/${script}" "$dockerfile"
  rg -q "scripts/${script}" "$workflow"
done

if rg -q 'prepare_window[.]R' "$dockerfile" "$workflow"; then
  echo "The mvSuSiE model image must not include prepare_window.R." >&2
  exit 1
fi

if rg -q 'workflows/trans_window_mvsusie[.]wdl' "$workflow"; then
  echo "The model image must not rebuild on WDL-only changes." >&2
  exit 1
fi

actionlint "$workflow"
echo "mvSuSiE model container definition passed"
