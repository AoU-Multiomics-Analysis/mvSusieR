#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

miniwdl check workflows/prepare_trans_window.wdl
rg -q '^version 1[.]0$' workflows/prepare_trans_window.wdl

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
  tabix \
  prepare_trans_window.R; do
  rg -q "$token" workflows/prepare_trans_window.wdl
done

rg -q 'call PrepareWindowGenotypes' workflows/prepare_trans_window.wdl
rg -q 'call PrepareWindowPhenotypes' workflows/prepare_trans_window.wdl
rg -q 'top_n_trans_phenotypes = 25' workflows/prepare_trans_window.wdl

if rg -q 'scatter[[:space:]]*\(' workflows/prepare_trans_window.wdl; then
  echo "The single-window preparation workflow must not scatter." >&2
  exit 1
fi

echo "Single-window preparation WDL validation passed"
