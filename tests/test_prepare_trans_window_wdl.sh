#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

miniwdl check workflows/prepare_trans_window.wdl
rg -q '^version 1[.]0$' workflows/prepare_trans_window.wdl
rg -q 'sep="[,]" phenotype_files' workflows/prepare_trans_window.wdl
rg -q 'sep="[,]" phenotype_modalities' workflows/prepare_trans_window.wdl

if rg -q '^[[:space:]]*String phenotype_(files|modalities)_csv[[:space:]]*=' workflows/prepare_trans_window.wdl; then
  echo "The prepare workflow must not use sep() as a standalone WDL expression." >&2
  exit 1
fi

for token in \
  PrepareTransWindow \
  PrepareWindowGenotypes \
  PrepareWindowPhenotypes \
  genome_dosage_tbi \
  trans_window_associations \
  phenotype_modalities \
  extract_cis_window_phenotypes \
  top_n_trans_phenotypes \
  top-n-trans-phenotypes \
  window_manifest \
  phenotype_data \
  tabix \
  gzip \
  prepare_trans_window.R; do
  rg -q "$token" workflows/prepare_trans_window.wdl
done

rg -q 'call PrepareWindowGenotypes' workflows/prepare_trans_window.wdl
rg -q 'call PrepareWindowPhenotypes' workflows/prepare_trans_window.wdl
rg -q 'top_n_trans_phenotypes = 25' workflows/prepare_trans_window.wdl
rg -q 'trans_window_associations = trans_window_associations' workflows/prepare_trans_window.wdl
rg -q 'distinct coordinates|seen\[row\]' workflows/prepare_trans_window.wdl
test "$(rg -c 'disks: "local-disk 500 SSD"' workflows/prepare_trans_window.wdl)" -eq 2
test "$(rg -c 'memory: "16 GiB"' workflows/prepare_trans_window.wdl)" -eq 2

if rg -q 'scatter[[:space:]]*\(' workflows/prepare_trans_window.wdl; then
  echo "The single-window preparation workflow must not scatter." >&2
  exit 1
fi

if rg -q 'windows_tsv|--windows' workflows/prepare_trans_window.wdl; then
  echo "The prepare workflow should use the long trans-window association file for coordinates." >&2
  exit 1
fi

echo "Single-window preparation WDL validation passed"
