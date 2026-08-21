version 1.0

workflow TransQTLBamFiltering {
  input {
    File input_bam
    File input_bai
    File low_mappability_bed
    File genes_gtf
    String sample_id
    String output_prefix = "transqtl"
    Float mappability_threshold = 1.0
    String strandedness = "rf"
    Boolean legacy = false
    Int threads = 4
  }

  call FilterTransQTLBam {
    input:
      input_bam = input_bam,
      input_bai = input_bai,
      low_mappability_bed = low_mappability_bed,
      genes_gtf = genes_gtf,
      sample_id = sample_id,
      output_prefix = output_prefix,
      mappability_threshold = mappability_threshold,
      strandedness = strandedness,
      legacy = legacy,
      threads = threads
  }

  output {
    File filtered_bam = FilterTransQTLBam.filtered_bam
    File filtered_bai = FilterTransQTLBam.filtered_bai
    File excluded_read_names = FilterTransQTLBam.excluded_read_names
    File filter_metrics = FilterTransQTLBam.filter_metrics
    File filtered_flagstat = FilterTransQTLBam.filtered_flagstat
    File gene_reads = FilterTransQTLBam.gene_reads
    File gene_tpm = FilterTransQTLBam.gene_tpm
    File rnaseqc_metrics = FilterTransQTLBam.rnaseqc_metrics
    File rnaseqc_coverage = FilterTransQTLBam.rnaseqc_coverage
  }
}

task FilterTransQTLBam {
  input {
    File input_bam
    File input_bai
    File low_mappability_bed
    File genes_gtf
    String sample_id
    String output_prefix
    Float mappability_threshold
    String strandedness
    Boolean legacy
    Int threads
  }

  command <<<
    set -euo pipefail

    log() {
      printf '[%s] [TransQTLBamFiltering] %s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2
    }

    log "Starting TransQTLBamFiltering for sample ~{sample_id} with ~{threads} threads"

    if [ ~{threads} -lt 1 ]; then
      echo "threads must be at least 1" >&2
      exit 1
    fi

    log "Preparing low-mappability BED and linking input BAM"
    mkdir -p bam_input
    ln -s ~{input_bam} bam_input/input.bam
    ln -s ~{input_bai} bam_input/input.bam.bai

    # The BED input is expected to be precomputed from an ENCODE 36-mer
    # mappability track using the documented score threshold and mismatch
    # definition. Normalize and sort it for samtools -L.
    awk 'BEGIN { OFS = "\t" }
      /^#/ { next }
      NF < 3 { next }
      $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 < $2 {
        print "Invalid low-mappability BED interval" > "/dev/stderr"
        exit 1
      }
      { print }
    ' ~{low_mappability_bed} > low_mappability.normalized.bed
    LC_ALL=C sort -k1,1V -k2,2n -k3,3n \
      low_mappability.normalized.bed > low_mappability.sorted.bed
    log "Low-mappability BED prepared"

    # Exclude a template if any of its alignments is multi-mapped or lacks
    # the NH tag. NH:i:1 is the unique-alignment criterion emitted by STAR.
    log "Identifying nonunique templates from NH tags"
    samtools view -@ ~{threads} bam_input/input.bam \
      | awk -F '\t' '
        {
          nh = ""
          for (i = 12; i <= NF; i++) {
            if ($i ~ /^NH:i:/) {
              nh = substr($i, 6)
              break
            }
          }
          if (nh == "" || nh != 1) print $1
        }
      ' \
      | LC_ALL=C sort -u > nonunique.templates.txt

    # Exclude a template if any aligned record overlaps the low-mappability
    # mask. samtools -L uses reference-coordinate overlap, so this applies to
    # the actual alignment rather than to the read sequence alone.
    log "Identifying low-mappability templates from reference overlap"
    samtools view -@ ~{threads} -L low_mappability.sorted.bed bam_input/input.bam \
      | cut -f1 \
      | LC_ALL=C sort -u > low_mappability.templates.txt

    LC_ALL=C sort -u nonunique.templates.txt low_mappability.templates.txt \
      > ~{output_prefix}.excluded_read_names.txt
    log "Template exclusion list created"

    # Remove every SAM record belonging to an excluded template, preserving
    # the header and therefore the read-group metadata.
    log "Writing filtered BAM"
    samtools view -@ ~{threads} -h bam_input/input.bam \
      | awk -F '\t' -v bad_names='~{output_prefix}.excluded_read_names.txt' '
        BEGIN {
          while ((getline name < bad_names) > 0) excluded[name] = 1
          close(bad_names)
        }
        /^@/ { print; next }
        !($1 in excluded) { print }
      ' \
      | samtools view -@ ~{threads} -b -o ~{output_prefix}.unsorted.bam -

    log "Sorting and indexing filtered BAM"
    samtools sort -@ ~{threads} \
      -o ~{output_prefix}.filtered.bam \
      ~{output_prefix}.unsorted.bam
    samtools index -@ ~{threads} \
      ~{output_prefix}.filtered.bam \
      ~{output_prefix}.filtered.bam.bai
    samtools quickcheck ~{output_prefix}.filtered.bam
    samtools flagstat -@ ~{threads} \
      ~{output_prefix}.filtered.bam \
      > ~{output_prefix}.filtered.flagstat.txt
    log "Filtered BAM sorted, indexed, and validated"

    strandedness_value="~{strandedness}"
    case "${strandedness_value}" in
      RF|rf|FR|fr) ;;
      *)
        echo "strandedness must be one of RF, rf, FR, or fr" >&2
        exit 1
        ;;
    esac

    legacy_flag=""
    if [ "~{legacy}" = "true" ]; then
      legacy_flag="--legacy"
    fi

    mkdir -p rnaseqc_output
    log "Running RNA-SeQC2"
    rnaseqc \
      ~{genes_gtf} \
      ~{output_prefix}.filtered.bam \
      rnaseqc_output \
      --sample ~{sample_id} \
      --stranded "${strandedness_value}" \
      --coverage \
      ${legacy_flag}

    test -s rnaseqc_output/~{sample_id}.gene_reads.gct
    test -s rnaseqc_output/~{sample_id}.gene_tpm.gct
    test -s rnaseqc_output/~{sample_id}.metrics.tsv
    test -s rnaseqc_output/~{sample_id}.coverage.tsv
    log "RNA-SeQC2 outputs verified"

    log "Computing filtering summaries"
    input_records="$(samtools view -@ ~{threads} -c bam_input/input.bam)"
    output_records="$(samtools view -@ ~{threads} -c ~{output_prefix}.filtered.bam)"
    removed_records="$((input_records - output_records))"
    samtools view -@ ~{threads} bam_input/input.bam \
      | cut -f1 \
      | LC_ALL=C sort -u > input.templates.txt
    samtools view -@ ~{threads} ~{output_prefix}.filtered.bam \
      | cut -f1 \
      | LC_ALL=C sort -u > output.templates.txt

    input_templates="$(wc -l < input.templates.txt | tr -d ' ')"
    output_templates="$(wc -l < output.templates.txt | tr -d ' ')"
    excluded_templates="$(wc -l < ~{output_prefix}.excluded_read_names.txt | tr -d ' ')"
    nonunique_templates="$(wc -l < nonunique.templates.txt | tr -d ' ')"
    low_mappability_templates="$(wc -l < low_mappability.templates.txt | tr -d ' ')"
    both_filter_templates="$(comm -12 nonunique.templates.txt low_mappability.templates.txt | wc -l | tr -d ' ')"
    nonunique_only_templates="$(comm -23 nonunique.templates.txt low_mappability.templates.txt | wc -l | tr -d ' ')"
    low_mappability_only_templates="$(comm -13 nonunique.templates.txt low_mappability.templates.txt | wc -l | tr -d ' ')"

    # Count alignment records removed because their template was affected by
    # each filter. A record can be counted in both filters when its template
    # has both an NH failure and a low-mappability overlap.
    samtools view -@ ~{threads} bam_input/input.bam \
      | awk -F '\t' \
        -v nonunique_file='nonunique.templates.txt' \
        -v low_mappability_file='low_mappability.templates.txt' '
        BEGIN {
          while ((getline name < nonunique_file) > 0) nonunique[name] = 1
          close(nonunique_file)
          while ((getline name < low_mappability_file) > 0) low_mappability[name] = 1
          close(low_mappability_file)
        }
        {
          has_nonunique = ($1 in nonunique)
          has_low_mappability = ($1 in low_mappability)
          if (has_nonunique) nonunique_records++
          if (has_low_mappability) low_mappability_records++
          if (has_nonunique && has_low_mappability) both_records++
          if (has_nonunique && !has_low_mappability) nonunique_only_records++
          if (!has_nonunique && has_low_mappability) low_mappability_only_records++
        }
        END {
          printf "nonunique_filter_alignment_records\t%d\n", nonunique_records + 0
          printf "low_mappability_filter_alignment_records\t%d\n", low_mappability_records + 0
          printf "both_filters_alignment_records\t%d\n", both_records + 0
          printf "nonunique_only_alignment_records\t%d\n", nonunique_only_records + 0
          printf "low_mappability_only_alignment_records\t%d\n", low_mappability_only_records + 0
        }
      ' > record_filter_metrics.tsv

    {
      printf 'metric\tvalue\n'
      printf 'total_reads_before_filtering\t%s\n' "${input_records}"
      printf 'total_reads_after_filtering\t%s\n' "${output_records}"
      printf 'reads_removed_total\t%s\n' "${removed_records}"
      printf 'input_alignment_records\t%s\n' "${input_records}"
      printf 'output_alignment_records\t%s\n' "${output_records}"
      printf 'removed_alignment_records\t%s\n' "${removed_records}"
      printf 'input_read_templates\t%s\n' "${input_templates}"
      printf 'output_read_templates\t%s\n' "${output_templates}"
      printf 'excluded_read_templates\t%s\n' "${excluded_templates}"
      printf 'nonunique_filter_read_templates\t%s\n' "${nonunique_templates}"
      printf 'low_mappability_filter_read_templates\t%s\n' "${low_mappability_templates}"
      printf 'both_filters_read_templates\t%s\n' "${both_filter_templates}"
      printf 'nonunique_only_read_templates\t%s\n' "${nonunique_only_templates}"
      printf 'low_mappability_only_read_templates\t%s\n' "${low_mappability_only_templates}"
      cat record_filter_metrics.tsv
      printf 'mappability_threshold\t%s\n' '~{mappability_threshold}'
      printf 'nh_policy\tNH:i:1 required; missing or non-1 excluded\n'
      printf 'low_mappability_policy\tany reference overlap excludes the whole template\n'
      printf 'low_mappability_bed\t%s\n' '~{low_mappability_bed}'
    } > ~{output_prefix}.filter_metrics.tsv

    log "Completed TransQTLBamFiltering: alignment_records_before=${input_records} alignment_records_after=${output_records} reads_removed=${removed_records} nonunique_templates=${nonunique_templates} low_mappability_templates=${low_mappability_templates}"

    rm -f ~{output_prefix}.unsorted.bam
  >>>

  output {
    File filtered_bam = output_prefix + ".filtered.bam"
    File filtered_bai = output_prefix + ".filtered.bam.bai"
    File excluded_read_names = output_prefix + ".excluded_read_names.txt"
    File filter_metrics = output_prefix + ".filter_metrics.tsv"
    File filtered_flagstat = output_prefix + ".filtered.flagstat.txt"
    File gene_reads = "rnaseqc_output/" + sample_id + ".gene_reads.gct"
    File gene_tpm = "rnaseqc_output/" + sample_id + ".gene_tpm.gct"
    File rnaseqc_metrics = "rnaseqc_output/" + sample_id + ".metrics.tsv"
    File rnaseqc_coverage = "rnaseqc_output/" + sample_id + ".coverage.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/transqtl-bam-filtering:latest"
    cpu: threads
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}
