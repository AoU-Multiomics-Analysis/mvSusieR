version 1.0

workflow IndexPhenotypeBed {
  input {
    File phenotype_bed
    String modality
    Int threads = 4
  }

  call BuildPhenotypeIndex {
    input:
      phenotype_bed = phenotype_bed,
      modality = modality,
      threads = threads
  }

  output {
    File indexed_phenotypes = BuildPhenotypeIndex.indexed_phenotypes
    File phenotype_tbi = BuildPhenotypeIndex.phenotype_tbi
    File phenotype_lookup = BuildPhenotypeIndex.phenotype_lookup
    File index_qc = BuildPhenotypeIndex.index_qc
  }
}

task BuildPhenotypeIndex {
  input {
    File phenotype_bed
    String modality
    Int threads
  }

  command <<<
    set -euo pipefail

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting phenotype indexing for ~{modality}."
    bash /opt/mvsusie/scripts/index_phenotype_bed.sh \
      '~{phenotype_bed}' \
      '~{modality}' \
      output \
      '~{threads}'
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Validating phenotype index outputs."
    test -s output/phenotypes.bed.gz
    test -s output/phenotypes.bed.gz.tbi
    test -s output/phenotype_lookup.tsv.gz
    test -s output/index_qc.tsv
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Phenotype indexing complete for ~{modality}."
  >>>

  output {
    File indexed_phenotypes = "output/phenotypes.bed.gz"
    File phenotype_tbi = "output/phenotypes.bed.gz.tbi"
    File phenotype_lookup = "output/phenotype_lookup.tsv.gz"
    File index_qc = "output/index_qc.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-prepare-window-phenotypes:latest"
    cpu: 4
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}
