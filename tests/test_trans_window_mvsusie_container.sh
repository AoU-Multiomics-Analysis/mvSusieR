#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

dockerfile="envs/trans-window-mvsusie.Dockerfile"
environment="envs/trans-window-mvsusie.environment.yml"
workflow=".github/workflows/trans-window-mvsusie-image.yml"

test -s "$dockerfile"
test -s "$environment"
test -s "$workflow"
rg -q '^FROM mambaorg/micromamba:2[.]3[.]3$' "$dockerfile"
rg -q 'micromamba config set channel_priority strict' "$dockerfile"
rg -q 'micromamba install --yes --name base --override-channels --strict-channel-priority --file /tmp/environment[.]yml' "$dockerfile"
rg -q 'ENV MAMBA_DOCKERFILE_ACTIVATE=1' "$dockerfile"
rg -q 'ENV PATH=/opt/conda/bin:' "$dockerfile"
if rg -q 'apt-get|install_github|rocker/r-ver' "$dockerfile"; then
  echo "The mvSuSiE image must install its runtime with micromamba." >&2
  exit 1
fi

for channel in dnachun conda-forge bioconda; do
  rg -Fq "  - ${channel}" "$environment"
done
for package in \
  'r-base=4.4' \
  'r-mvsusier=0.3.0' \
  'r-susier>=0.15' \
  r-mashr \
  r-data.table \
  r-dplyr \
  r-ggplot2 \
  r-optparse \
  r-purrr \
  r-readr \
  r-r.utils \
  r-stringr \
  r-tibble \
  ripgrep; do
  rg -Fq "  - ${package}" "$environment"
done
rg -q 'packageVersion[(]"mvsusieR"[)] >= "0.3.0"' "$workflow"
rg -q 'packageVersion[(]"susieR"[)] >= "0.15.0"' "$workflow"
rg -q 'mvsusie_plot' "$workflow"

for script in \
  trans_window_io.R \
  trans_window_logging.R \
  trans_window_preprocess.R \
  trans_window_model.R \
  trans_window_prior.R \
  trans_window_cli.R \
  prepare_window.R \
  fit_window.R \
  run_window_mvsusie.R \
  summarize_window.R \
  merge_window_outputs.R \
  plot_window_mvsusie.R; do
  rg -q "scripts/${script}" "$dockerfile"
  rg -q "scripts/${script}" "$workflow"
done

rg -q 'compute_marginal_bhat_shat_matrix' "$workflow"
rg -q 'tests/test_trans_window_wdl_smoke[.]sh' "$workflow"
rg -q 'python3 -m pip install miniwdl' "$workflow"

rg -q 'workflows/trans_window_mvsusie[.]wdl' "$workflow"

actionlint "$workflow"
echo "mvSuSiE model container definition passed"
