#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

miniwdl check workflows/transqtl_bam_filtering.wdl

for token in \
  TransQTLBamFiltering \
  FilterTransQTLBam \
  write_filter_metrics \
  preemptible_tries \
  memory \
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
  'if samtools index' \
  'Filtered BAM was already coordinate sorted; skipped resorting' \
  'Filtered BAM was not indexable as coordinate sorted; sorting and retrying' \
  'log()' \
  'Starting TransQTLBamFiltering' \
  'Preparing low-mappability BED' \
  'Identifying nonunique templates' \
  'Identifying low-mappability templates' \
  'Writing filtered BAM' \
  'Validating filtered BAM order and creating index' \
  'Running RNA-SeQC2' \
  'Computing filtering summaries' \
  'Completed TransQTLBamFiltering'; do
  rg -q "$token" workflows/transqtl_bam_filtering.wdl
done

rg -q 'Boolean write_filter_metrics = true' workflows/transqtl_bam_filtering.wdl
rg -q 'write_filter_metrics = write_filter_metrics' workflows/transqtl_bam_filtering.wdl
rg -q 'Int preemptible_tries = 1' workflows/transqtl_bam_filtering.wdl
rg -q 'preemptible_tries = preemptible_tries' workflows/transqtl_bam_filtering.wdl
rg -q 'preemptible: preemptible_tries' workflows/transqtl_bam_filtering.wdl
rg -q 'String memory = "16 GiB"' workflows/transqtl_bam_filtering.wdl
rg -q 'memory = memory' workflows/transqtl_bam_filtering.wdl
rg -q 'memory: memory' workflows/transqtl_bam_filtering.wdl
rg -q 'File\? filter_metrics' workflows/transqtl_bam_filtering.wdl
rg -q 'if \[ "~\{write_filter_metrics\}" = "true" \]' workflows/transqtl_bam_filtering.wdl
rg -q 'filter_metrics = FilterTransQTLBam.filter_metrics' workflows/transqtl_bam_filtering.wdl

for output_name in \
  'TransQTLFiltered.bam' \
  'TransQTLFiltered.bam.bai' \
  'TransQTLFiltered.excluded_read_names.txt' \
  'TransQTLFiltered.filter_metrics.tsv' \
  'TransQTLFiltered.flagstat.txt' \
  'TransQTLFiltered.gene_reads.gct' \
  'TransQTLFiltered.gene_tpm.gct' \
  'TransQTLFiltered.rnaseqc_metrics.tsv' \
  'TransQTLFiltered.rnaseqc_coverage.tsv'; do
  rg -q "$output_name" workflows/transqtl_bam_filtering.wdl
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
