version 1.0

workflow PrepareTransWindow {
  input {
    File windows_tsv
    String window_id
    File genome_dosage
    File genome_dosage_tbi
    File trans_window_associations
    Array[File] phenotype_files
    Array[String] phenotype_modalities
    Boolean extract_cis_window_phenotypes = true
    Int top_n_trans_phenotypes = 25
  }

  String phenotype_files_csv = sep(",", phenotype_files)
  String phenotype_modalities_csv = sep(",", phenotype_modalities)

  call PrepareWindowGenotypes {
    input:
      windows_tsv = windows_tsv,
      window_id = window_id,
      genome_dosage = genome_dosage,
      genome_dosage_tbi = genome_dosage_tbi
  }

  call PrepareWindowPhenotypes {
    input:
      windows_tsv = windows_tsv,
      window_id = window_id,
      trans_window_associations = trans_window_associations,
      phenotype_files = phenotype_files,
      phenotype_files_csv = phenotype_files_csv,
      phenotype_modalities_csv = phenotype_modalities_csv,
      extract_cis_window_phenotypes = extract_cis_window_phenotypes,
      top_n_trans_phenotypes = top_n_trans_phenotypes
  }

  output {
    File window_dosage = PrepareWindowGenotypes.window_dosage
    File window_manifest = PrepareWindowGenotypes.window_manifest
    File window_phenotypes = PrepareWindowPhenotypes.window_phenotypes
    File phenotype_data = PrepareWindowPhenotypes.phenotype_data
    File window_qc = PrepareWindowPhenotypes.window_qc
  }
}

task PrepareWindowGenotypes {
  input {
    File windows_tsv
    String window_id
    File genome_dosage
    File genome_dosage_tbi
  }

  command <<<
    set -euo pipefail

    mkdir -p output

    dosage_name="$(basename ~{genome_dosage})"
    ln -sf ~{genome_dosage_tbi} "${dosage_name}.tbi"

    window_row="$(awk -F '\t' -v requested_id='~{window_id}' '
      NR == 1 {
        for (i = 1; i <= NF; i++) column[$i] = i
        next
      }
      $(column["window_id"]) == requested_id {
        print $(column["chrom"]) "\t" $(column["start"]) "\t" $(column["end"])
        matches++
      }
      END {
        if (matches != 1) exit 1
      }
    ' ~{windows_tsv})"

    IFS=$'\t' read -r window_chrom window_start window_end <<< "${window_row}"
    tabix -H "${dosage_name}" > output/window_dosage.tsv
    tabix "${dosage_name}" \
      "${window_chrom}:$((window_start + 1))-${window_end}" \
      >> output/window_dosage.tsv
    test "$(wc -l < output/window_dosage.tsv)" -gt 1

    {
      printf 'window_id\tchrom\tstart\tend\tdosage_file\n'
      printf '%s\t%s\t%s\t%s\t%s\n' \
        '~{window_id}' \
        "${window_chrom}" \
        "${window_start}" \
        "${window_end}" \
        'window_dosage.tsv'
    } > output/window_manifest.tsv
  >>>

  output {
    File window_dosage = "output/window_dosage.tsv"
    File window_manifest = "output/window_manifest.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-prepare-window-genotypes:latest"
    cpu: 2
    memory: "2 GiB"
  }
}

task PrepareWindowPhenotypes {
  input {
    File windows_tsv
    String window_id
    File trans_window_associations
    Array[File] phenotype_files
    String phenotype_files_csv
    String phenotype_modalities_csv
    Boolean extract_cis_window_phenotypes
    Int top_n_trans_phenotypes
  }

  command <<<
    set -euo pipefail
    test ~{length(phenotype_files)} -gt 0

    Rscript /opt/mvsusie/scripts/prepare_trans_window.R \
      --windows ~{windows_tsv} \
      --window-id ~{window_id} \
      --trans-associations ~{trans_window_associations} \
      --phenotype-files "~{phenotype_files_csv}" \
      --phenotype-modalities "~{phenotype_modalities_csv}" \
      --extract-cis-window-phenotypes ~{extract_cis_window_phenotypes} \
      --top-n-trans-phenotypes ~{top_n_trans_phenotypes} \
      --output-dir output
  >>>

  output {
    File window_phenotypes = "output/window_phenotypes.tsv"
    File phenotype_data = "output/window_phenotypes.bed.gz"
    File window_qc = "output/window_qc.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-prepare-window-phenotypes:latest"
    cpu: 2
    memory: "16 GiB"
  }
}
