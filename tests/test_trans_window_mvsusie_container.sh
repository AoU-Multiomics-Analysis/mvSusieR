#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

dockerfile="envs/trans-window-mvsusie.Dockerfile"
workflow=".github/workflows/trans-window-mvsusie-image.yml"

test -s "$dockerfile"
test -s "$workflow"
rg -q '^FROM rocker/r-ver:4[.]4[.]1$' "$dockerfile"
rg -q '65f3586a865fb6748cb4f9df50510ac577706348' "$dockerfile"
rg -q 'ebd1133953005fa70c6b338727b5fe9222e2a1c2' "$dockerfile"
rg -q 'mashr' "$dockerfile"
rg -q 'ggplot2' "$dockerfile"
rg -q 'ripgrep' "$dockerfile"
rg -q 'install_github' "$dockerfile"
for package in dplyr purrr readr R.utils stringr tibble; do
  rg -q "    ${package}" "$dockerfile"
done
rg -q 'RemoteSha.*65f3586a865fb6748cb4f9df50510ac577706348' "$workflow"
rg -q 'RemoteSha.*ebd1133953005fa70c6b338727b5fe9222e2a1c2' "$workflow"
rg -q 'mvsusie_plot' "$workflow"

for script in \
  trans_window_io.R \
  trans_window_logging.R \
  trans_window_preprocess.R \
  trans_window_model.R \
  trans_window_prior.R \
  trans_window_cli.R \
  fit_window.R \
  run_window_mvsusie.R \
  summarize_window.R \
  merge_window_outputs.R \
  plot_window_mvsusie.R; do
  rg -q "scripts/${script}" "$dockerfile"
  rg -q "scripts/${script}" "$workflow"
done

rg -q 'compute_marginal_bhat_shat_matrix' "$workflow"

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
