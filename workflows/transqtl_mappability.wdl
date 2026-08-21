version 1.0

workflow TransQTLMappability {
  input {
    File reference_fasta
    Int kmer_length = 146
    Int max_mismatches = 2
    Float mappability_threshold = 1.0
    String output_prefix = "GRCh38.K146.m2"
    Int threads = 16
    Int preemptible_tries = 3
  }

  call GenerateTransQTLMappability {
    input:
      reference_fasta = reference_fasta,
      kmer_length = kmer_length,
      max_mismatches = max_mismatches,
      mappability_threshold = mappability_threshold,
      output_prefix = output_prefix,
      threads = threads,
      preemptible_tries = preemptible_tries
  }

  output {
    File low_mappability_bed = GenerateTransQTLMappability.low_mappability_bed
    File mappability_bedgraph = GenerateTransQTLMappability.mappability_bedgraph
    File mappability_metadata = GenerateTransQTLMappability.mappability_metadata
  }
}

task GenerateTransQTLMappability {
  input {
    File reference_fasta
    Int kmer_length
    Int max_mismatches
    Float mappability_threshold
    String output_prefix
    Int threads
    Int preemptible_tries
  }

  command <<<
    set -euo pipefail

    log() {
      printf '[%s] [TransQTLMappability] %s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2
    }

    log "Starting TransQTLMappability with k=~{kmer_length}, e=~{max_mismatches}, threshold=~{mappability_threshold}, threads=~{threads}"

    if [ ~{kmer_length} -lt 1 ]; then
      echo "kmer_length must be at least 1" >&2
      exit 1
    fi
    if [ ~{max_mismatches} -lt 0 ]; then
      echo "max_mismatches must be nonnegative" >&2
      exit 1
    fi
    if ! awk -v threshold='~{mappability_threshold}' \
      'BEGIN { exit !(threshold > 0 && threshold <= 1) }'; then
      echo "mappability_threshold must be greater than 0 and at most 1" >&2
      exit 1
    fi
    if [ ~{threads} -lt 1 ]; then
      echo "threads must be at least 1" >&2
      exit 1
    fi
    test -s ~{reference_fasta}

    log "Copying reference FASTA and creating GenMap index"
    cp ~{reference_fasta} reference.fasta
    genmap index \
      -F reference.fasta \
      -I genmap_index

    log "Computing GenMap mappability"
    genmap map \
      -K ~{kmer_length} \
      -E ~{max_mismatches} \
      -I genmap_index \
      -O genmap_output \
      -t \
      -w \
      -bg

    mappability_source="$(find . -type f -name '*.bedGraph' -print -quit)"
    if [ -z "${mappability_source}" ]; then
      echo "GenMap did not produce a BEDGraph" >&2
      exit 1
    fi
    cp "${mappability_source}" "~{output_prefix}.mappability.bedGraph"
    test -s "~{output_prefix}.mappability.bedGraph"
    log "Mappability BEDGraph created"

    log "Applying mappability threshold and merging intervals"
    awk -v threshold='~{mappability_threshold}' \
      'BEGIN { OFS = "\t" }
       /^#/ { next }
       NF >= 4 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && ($4 + 0) < threshold {
         print $1, $2, $3
       }' \
      "~{output_prefix}.mappability.bedGraph" \
      > low_mappability.candidate.bed
    : > "~{output_prefix}.low_mappability.bed"
    if [ -s low_mappability.candidate.bed ]; then
      bedtools sort -i low_mappability.candidate.bed \
        | bedtools merge -i - \
        > "~{output_prefix}.low_mappability.bed"
    fi
    low_mappability_intervals="$(wc -l < "~{output_prefix}.low_mappability.bed" | tr -d ' ')"
    log "Low-mappability BED created with ${low_mappability_intervals} merged intervals"

    reference_size_bytes="$(wc -c < reference.fasta | tr -d ' ')"
    reference_sha256="$(sha256sum reference.fasta | awk '{ print $1 }')"
    {
      printf 'metric\tvalue\n'
      printf 'kmer_length\t%s\n' '~{kmer_length}'
      printf 'max_mismatches\t%s\n' '~{max_mismatches}'
      printf 'mappability_threshold\t%s\n' '~{mappability_threshold}'
      printf 'reference_fasta\t%s\n' '~{reference_fasta}'
      printf 'reference_size_bytes\t%s\n' "${reference_size_bytes}"
      printf 'reference_sha256\t%s\n' "${reference_sha256}"
      printf 'low_mappability_intervals\t%s\n' "${low_mappability_intervals}"
    } > "~{output_prefix}.mappability.metadata.tsv"

    test -s "~{output_prefix}.mappability.metadata.tsv"
    log "Mappability outputs verified"
    log "Completed TransQTLMappability"
  >>>

  output {
    File low_mappability_bed = output_prefix + ".low_mappability.bed"
    File mappability_bedgraph = output_prefix + ".mappability.bedGraph"
    File mappability_metadata = output_prefix + ".mappability.metadata.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/transqtl-mappability:latest"
    cpu: threads
    memory: "64 GiB"
    disks: "local-disk 200 SSD"
    preemptible: preemptible_tries
  }
}
