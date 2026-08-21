#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

miniwdl check workflows/transqtl_bam_filtering.wdl

for token in \
  TransQTLBamFiltering \
  FilterTransQTLBam \
  input_bam \
  input_bai \
  low_mappability_bed \
  genes_gtf \
  sample_id \
  strandedness \
  excluded_read_names \
  filtered_bam \
  filtered_bai \
  rnaseqc \
  gene_reads.gct \
  gene_tpm.gct \
  NH \
  'samtools view' \
  'samtools sort' \
  'samtools index' \
  'samtools view.*-L'; do
  rg -q "$token" workflows/transqtl_bam_filtering.wdl
done

rg -q 'outFilterMultimapNmax|mappability_threshold' workflows/transqtl_bam_filtering.wdl
rg -q 'low_mappability_bed' workflows/transqtl_bam_filtering.wdl
rg -q 'legacy' workflows/transqtl_bam_filtering.wdl

echo "TransQTLBamFiltering WDL validation passed"
