#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

miniwdl check workflows/prepare_trans_window.wdl

for declaration in \
  'File expression_phenotypes' \
  'File splicing_phenotypes' \
  'File protein_phenotypes' \
  'File target_phenotypes' \
  'Int top_n_expression = 25' \
  'Int top_n_splicing = 25' \
  'Int top_n_protein = 15'; do
  rg -q "$declaration" workflows/prepare_trans_window.wdl
done

for removed in \
  phenotype_files \
  phenotype_modalities \
  extract_cis_window_phenotypes \
  top_n_trans_phenotypes; do
  if rg -q "$removed" workflows/prepare_trans_window.wdl; then
    echo "Removed prepare input remains: $removed" >&2
    exit 1
  fi
done

for cli_flag in \
  expression-phenotypes \
  splicing-phenotypes \
  protein-phenotypes \
  target-phenotypes \
  top-n-expression \
  top-n-splicing \
  top-n-protein; do
  rg -q -- "--$cli_flag" workflows/prepare_trans_window.wdl
done

rg -q 'call PrepareWindowGenotypes' workflows/prepare_trans_window.wdl
rg -q 'call PrepareWindowPhenotypes' workflows/prepare_trans_window.wdl
test "$(rg -c 'disks: "local-disk 500 SSD"' workflows/prepare_trans_window.wdl)" -eq 2
test "$(rg -c 'memory: "16 GiB"' workflows/prepare_trans_window.wdl)" -eq 2

if rg -q 'scatter[[:space:]]*\(' workflows/prepare_trans_window.wdl; then
  echo "The single-window preparation workflow must not scatter." >&2
  exit 1
fi

echo "Joint preparation WDL validation passed"
