#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

miniwdl check workflows/transqtl_mappability.wdl

for token in \
  TransQTLMappability \
  GenerateTransQTLMappability \
  reference_fasta \
  kmer_length \
  max_mismatches \
  mappability_threshold \
  preemptible_tries \
  'genmap index' \
  'genmap map' \
  'bedtools merge' \
  'if [ -s low_mappability.candidate.bed ]' \
  low_mappability_bed \
  mappability_bedgraph \
  mappability_metadata \
  'log()' \
  'Mappability outputs verified'; do
  rg -Fq "$token" workflows/transqtl_mappability.wdl
done

rg -q 'threshold > 0 && threshold <= 1' workflows/transqtl_mappability.wdl
rg -q 'preemptible: preemptible_tries' workflows/transqtl_mappability.wdl

echo "TransQTLMappability WDL validation passed"
