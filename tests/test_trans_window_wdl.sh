#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

miniwdl check workflows/trans_window_mvsusie.wdl

for token in \
  PrepareWindowData \
  RunMvSusie \
  SummarizeMvSusie \
  MergeWindowOutputs \
  variant_pip \
  credible_sets \
  component_effects \
  window_qc; do
  rg -q "$token" workflows/trans_window_mvsusie.wdl
done

echo "Trans-window WDL validation passed"
