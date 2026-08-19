#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

miniwdl check workflows/trans_window_mvsusie.wdl
miniwdl check workflows/prepare_trans_window.wdl
rg -q '^version 1[.]0$' workflows/trans_window_mvsusie.wdl
rg -q 'covariate_modalities' workflows/trans_window_mvsusie.wdl
rg -q 'sep="[,]" covariate_files' workflows/trans_window_mvsusie.wdl
rg -q 'File phenotype_data' workflows/trans_window_mvsusie.wdl
test "$(rg -c 'disks: "local-disk 500 SSD"' workflows/trans_window_mvsusie.wdl)" -eq 4
test "$(rg -c 'memory: "16 GiB"' workflows/trans_window_mvsusie.wdl)" -eq 4

if rg -q 'Array\[File\] phenotype_files' workflows/trans_window_mvsusie.wdl; then
  echo "The mvSuSiE workflow should consume one combined phenotype file." >&2
  exit 1
fi

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
bash tests/test_prepare_trans_window_wdl.sh
bash tests/test_prepare_trans_window_containers.sh
