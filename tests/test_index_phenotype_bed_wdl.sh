#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

wdl="workflows/index_phenotype_bed.wdl"
miniwdl check "$wdl"

rg -Fq 'workflow IndexPhenotypeBed' "$wdl"
rg -Fq 'File phenotype_bed' "$wdl"
rg -Fq 'String modality' "$wdl"
rg -Fq 'Int threads = 4' "$wdl"
for output_name in indexed_phenotypes phenotype_tbi phenotype_lookup index_qc; do
  rg -q "File ${output_name}[[:space:]]*=" "$wdl"
done
rg -Fq '/opt/mvsusie/scripts/index_phenotype_bed.sh' "$wdl"
test "$(rg -Fc "date '+%Y-%m-%d %H:%M:%S'" "$wdl")" -ge 2
rg -Fq 'cpu: 4' "$wdl"
rg -Fq 'memory: "16 GiB"' "$wdl"
rg -Fq 'disks: "local-disk 500 SSD"' "$wdl"

echo "Phenotype index WDL validation passed"
