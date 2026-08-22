version 1.0

workflow TransWindowMvSusie {
  input {
    File windows_tsv
    File window_phenotypes_tsv
    File phenotype_data
    Array[File] covariate_files
    Array[String] covariate_modalities = ["shared"]
    File? keep_samples
    Int L = 10
    Int max_iter = 100
    Float tol = 1e-4
    Float coverage = 0.95
    Float min_abs_corr = 0.5
    Float min_genotype_variance = 1e-8
    Float min_phenotype_variance = 1e-8
    Int n_thread = 1
    String prior_method = "canonical"
    Int mashr_n_pca = 5
    Int? mashr_seed
  }

  Array[Array[String]] window_rows = read_tsv(windows_tsv)

  scatter (window_index in range(length(window_rows) - 1)) {
    Array[String] window = window_rows[window_index + 1]
    File dosage = window[4]

    call RunMvSusie {
      input:
        windows_tsv = windows_tsv,
        window_phenotypes_tsv = window_phenotypes_tsv,
        window_id = window[0],
        dosage = dosage,
        phenotype_data = phenotype_data,
        covariate_files = covariate_files,
        covariate_modalities = covariate_modalities,
        keep_samples = keep_samples,
        min_genotype_variance = min_genotype_variance,
        min_phenotype_variance = min_phenotype_variance,
        L = L,
        max_iter = max_iter,
        tol = tol,
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        n_thread = n_thread,
        prior_method = prior_method,
        mashr_n_pca = mashr_n_pca,
        mashr_seed = mashr_seed
    }

    call SummarizeMvSusie {
      input:
        prepared_window = RunMvSusie.prepared_window,
        mvsusie_fit = RunMvSusie.mvsusie_fit
    }
  }

  call MergeWindowOutputs {
    input:
      variant_pips = SummarizeMvSusie.variant_pip,
      credible_sets = SummarizeMvSusie.credible_sets,
      component_effects = SummarizeMvSusie.component_effects,
      window_qc = SummarizeMvSusie.window_qc
  }

  output {
    Array[File] prepared_windows = RunMvSusie.prepared_window
    Array[File] mvsusie_fits = RunMvSusie.mvsusie_fit
    Array[File] variant_pip = SummarizeMvSusie.variant_pip
    Array[File] credible_sets = SummarizeMvSusie.credible_sets
    Array[File] component_effects = SummarizeMvSusie.component_effects
    Array[File] window_qc = SummarizeMvSusie.window_qc
    File merged_variant_pip = MergeWindowOutputs.merged_variant_pip
    File merged_credible_sets = MergeWindowOutputs.merged_credible_sets
    File merged_component_effects = MergeWindowOutputs.merged_component_effects
    File merged_window_qc = MergeWindowOutputs.merged_window_qc
  }
}

task RunMvSusie {
  input {
    File windows_tsv
    File window_phenotypes_tsv
    String window_id
    File dosage
    File phenotype_data
    Array[File] covariate_files
    Array[String] covariate_modalities
    File? keep_samples
    Float min_genotype_variance
    Float min_phenotype_variance
    Int L
    Int max_iter
    Float tol
    Float coverage
    Float min_abs_corr
    Int n_thread
    String prior_method
    Int mashr_n_pca
    Int? mashr_seed
  }

  command <<<
    set -euo pipefail

    Rscript /opt/mvsusie/scripts/run_window_mvsusie.R \
      --windows ~{windows_tsv} \
      --window-phenotypes ~{window_phenotypes_tsv} \
      --window-id ~{window_id} \
      --dosage ~{dosage} \
      --phenotype-files ~{phenotype_data} \
      --covariate-files "~{sep="," covariate_files}" \
      --covariate-modalities "~{sep="," covariate_modalities}" \
      ~{if defined(keep_samples) then "--keep-samples " + select_first([keep_samples]) else ""} \
      --min-genotype-variance ~{min_genotype_variance} \
      --min-phenotype-variance ~{min_phenotype_variance} \
      --L ~{L} \
      --max-iter ~{max_iter} \
      --tol ~{tol} \
      --coverage ~{coverage} \
      --min-abs-corr ~{min_abs_corr} \
      --n-thread ~{n_thread} \
      --prior-method ~{prior_method} \
      --mashr-n-pca ~{mashr_n_pca} \
      ~{if defined(mashr_seed) then "--mashr-seed " + select_first([mashr_seed]) else ""} \
      --prepared-output prepared_window.rds \
      --fit-output mvsusie_fit.rds
  >>>

  output {
    File prepared_window = "prepared_window.rds"
    File mvsusie_fit = "mvsusie_fit.rds"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-trans-window-mvsusie:latest"
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}

task SummarizeMvSusie {
  input {
    File prepared_window
    File mvsusie_fit
  }

  command <<<
    set -euo pipefail
    mkdir -p window_outputs

    Rscript /opt/mvsusie/scripts/summarize_window.R \
      --prepared ~{prepared_window} \
      --fit ~{mvsusie_fit} \
      --output-dir window_outputs
  >>>

  output {
    File variant_pip = "window_outputs/variant_pip.tsv.gz"
    File credible_sets = "window_outputs/credible_sets.tsv.gz"
    File component_effects = "window_outputs/component_effects.tsv.gz"
    File window_qc = "window_outputs/window_qc.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-trans-window-mvsusie:latest"
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}

task MergeWindowOutputs {
  input {
    Array[File] variant_pips
    Array[File] credible_sets
    Array[File] component_effects
    Array[File] window_qc
  }

  command <<<
    set -euo pipefail
    mkdir -p merged

    Rscript /opt/mvsusie/scripts/merge_window_outputs.R \
      --variant-pips "~{sep="," variant_pips}" \
      --credible-sets "~{sep="," credible_sets}" \
      --component-effects "~{sep="," component_effects}" \
      --window-qc "~{sep="," window_qc}" \
      --output-dir merged
  >>>

  output {
    File merged_variant_pip = "merged/variant_pip.tsv.gz"
    File merged_credible_sets = "merged/credible_sets.tsv.gz"
    File merged_component_effects = "merged/component_effects.tsv.gz"
    File merged_window_qc = "merged/window_qc.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-trans-window-mvsusie:latest"
    cpu: 1
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}
