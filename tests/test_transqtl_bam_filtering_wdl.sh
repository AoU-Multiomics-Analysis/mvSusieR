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
  filter_metrics \
  rnaseqc \
  gene_reads.gct \
  gene_tpm.gct \
  NH \
  'samtools view' \
  'samtools sort' \
  'samtools index' \
  'samtools view.*-L' \
  'log()' \
  'Starting TransQTLBamFiltering' \
  'Preparing low-mappability BED' \
  'Identifying nonunique templates' \
  'Identifying low-mappability templates' \
  'Writing filtered BAM' \
  'Sorting and indexing filtered BAM' \
  'Running RNA-SeQC2' \
  'Computing filtering summaries' \
  'Completed TransQTLBamFiltering'; do
  rg -q "$token" workflows/transqtl_bam_filtering.wdl
done

rg -q 'outFilterMultimapNmax|mappability_threshold' workflows/transqtl_bam_filtering.wdl
rg -q 'low_mappability_bed' workflows/transqtl_bam_filtering.wdl
rg -q 'legacy' workflows/transqtl_bam_filtering.wdl
for metric in \
  input_alignment_records \
  output_alignment_records \
  removed_alignment_records \
  total_reads_before_filtering \
  total_reads_after_filtering \
  reads_removed_total \
  nonunique_filter_alignment_records \
  low_mappability_filter_alignment_records \
  both_filters_alignment_records \
  input_read_templates \
  output_read_templates \
  excluded_read_templates; do
  rg -q "$metric" workflows/transqtl_bam_filtering.wdl
done

echo "TransQTLBamFiltering WDL validation passed"
