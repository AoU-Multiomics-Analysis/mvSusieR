version 1.0

workflow TransWindowMvSusie {
  input {
    File windows_tsv
    File window_phenotypes_tsv
    File phenotype_data
    File expression_covariates
    File splicing_covariates
    File protein_covariates
    File? keep_samples
    Int start_L = 10
    Int step_L = 5
    Int max_L = 40
    Float greedy_lbf_cutoff = 1.0
    Int max_iter = 100
    Float tol = 1e-4
    Float coverage = 0.95
    Float min_abs_corr = 0.5
    Float min_genotype_variance = 1e-8
    Float min_phenotype_variance = 1e-8
    Int n_thread = 1
    Int mashr_n_pca = 5
    Float mashr_strong_lfsr = 0.05
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
        expression_covariates = expression_covariates,
        splicing_covariates = splicing_covariates,
        protein_covariates = protein_covariates,
        keep_samples = keep_samples,
        min_genotype_variance = min_genotype_variance,
        min_phenotype_variance = min_phenotype_variance,
        start_L = start_L,
        step_L = step_L,
        max_L = max_L,
        greedy_lbf_cutoff = greedy_lbf_cutoff,
        max_iter = max_iter,
        tol = tol,
        coverage = coverage,
        min_abs_corr = min_abs_corr,
        n_thread = n_thread,
        mashr_n_pca = mashr_n_pca,
        mashr_strong_lfsr = mashr_strong_lfsr,
        mashr_seed = mashr_seed
    }

    call SummarizeMvSusie {
      input:
        prepared_window = RunMvSusie.prepared_window,
        mvsusie_fit = RunMvSusie.mvsusie_fit
    }

    call PlotMvSusie {
      input:
        prepared_window = RunMvSusie.prepared_window,
        mvsusie_fit = RunMvSusie.mvsusie_fit
    }
  }

  call MergeWindowOutputs {
    input:
      variant_pips = SummarizeMvSusie.variant_pip,
      credible_sets = SummarizeMvSusie.credible_sets,
      credible_set_members = SummarizeMvSusie.credible_set_members,
      component_feature_support = SummarizeMvSusie.component_feature_support,
      window_qc = SummarizeMvSusie.window_qc
  }

  output {
    Array[File] prepared_windows = RunMvSusie.prepared_window
    Array[File] mvsusie_fits = RunMvSusie.mvsusie_fit
    Array[File] mashr_training = RunMvSusie.mashr_training
    Array[File] greedy_L_history = RunMvSusie.greedy_L_history
    Array[File] covariate_provenance = RunMvSusie.covariate_provenance
    Array[File] run_stdout = RunMvSusie.run_stdout
    Array[File] run_stderr = RunMvSusie.run_stderr
    Array[File] session_info = RunMvSusie.session_info
    Array[File] variant_pip = SummarizeMvSusie.variant_pip
    Array[File] credible_sets = SummarizeMvSusie.credible_sets
    Array[File] credible_set_members = SummarizeMvSusie.credible_set_members
    Array[File] component_feature_support = SummarizeMvSusie.component_feature_support
    Array[File] window_qc = SummarizeMvSusie.window_qc
    Array[File] effect_plot_png = PlotMvSusie.effect_plot_png
    Array[File] effect_plot_pdf = PlotMvSusie.effect_plot_pdf
    Array[File] effect_plot_rds = PlotMvSusie.effect_plot_rds
    File merged_variant_pip = MergeWindowOutputs.merged_variant_pip
    File merged_credible_sets = MergeWindowOutputs.merged_credible_sets
    File merged_credible_set_members = MergeWindowOutputs.merged_credible_set_members
    File merged_component_feature_support = MergeWindowOutputs.merged_component_feature_support
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
    File expression_covariates
    File splicing_covariates
    File protein_covariates
    File? keep_samples
    Float min_genotype_variance
    Float min_phenotype_variance
    Int start_L
    Int step_L
    Int max_L
    Float greedy_lbf_cutoff
    Int max_iter
    Float tol
    Float coverage
    Float min_abs_corr
    Int n_thread
    Int mashr_n_pca
    Float mashr_strong_lfsr
    Int? mashr_seed
  }

  command <<<
    set -euo pipefail

    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }

    log "Starting joint mvSuSiE for ~{window_id}: start L=~{start_L}, step L=~{step_L}, maximum L=~{max_L}."
    log "Resolved controls: cutoff=~{greedy_lbf_cutoff}, PCA=~{mashr_n_pca}, coverage=~{coverage}, minimum correlation=~{min_abs_corr}."
    Rscript /opt/mvsusie/scripts/run_window_mvsusie.R \
      --windows ~{windows_tsv} \
      --window-phenotypes ~{window_phenotypes_tsv} \
      --window-id ~{window_id} \
      --dosage ~{dosage} \
      --phenotype-files ~{phenotype_data} \
      --expression-covariates ~{expression_covariates} \
      --splicing-covariates ~{splicing_covariates} \
      --protein-covariates ~{protein_covariates} \
      ~{if defined(keep_samples) then "--keep-samples " + select_first([keep_samples]) else ""} \
      --min-genotype-variance ~{min_genotype_variance} \
      --min-phenotype-variance ~{min_phenotype_variance} \
      --start-L ~{start_L} \
      --step-L ~{step_L} \
      --max-L ~{max_L} \
      --greedy-lbf-cutoff ~{greedy_lbf_cutoff} \
      --max-iter ~{max_iter} \
      --tol ~{tol} \
      --coverage ~{coverage} \
      --min-abs-corr ~{min_abs_corr} \
      --n-thread ~{n_thread} \
      --mashr-n-pca ~{mashr_n_pca} \
      --mashr-strong-lfsr ~{mashr_strong_lfsr} \
      ~{if defined(mashr_seed) then "--mashr-seed " + select_first([mashr_seed]) else ""} \
      --covariate-provenance-output covariate_provenance.tsv.gz \
      --mashr-output mashr_training.rds \
      --greedy-history-output greedy_L_history.tsv \
      --prepared-output prepared_window.rds \
      --fit-output mvsusie_fit.rds
    log "Writing the R session information."
    Rscript -e 'writeLines(capture.output(sessionInfo()), "session_info.txt")'
    log "Verifying the joint model outputs."
    test -s prepared_window.rds
    test -s mvsusie_fit.rds
    test -s mashr_training.rds
    test -s greedy_L_history.tsv
    test -s covariate_provenance.tsv.gz
    test -s session_info.txt
    log "Completed joint mvSuSiE for ~{window_id}."
  >>>

  output {
    File prepared_window = "prepared_window.rds"
    File mvsusie_fit = "mvsusie_fit.rds"
    File mashr_training = "mashr_training.rds"
    File greedy_L_history = "greedy_L_history.tsv"
    File covariate_provenance = "covariate_provenance.tsv.gz"
    File run_stdout = stdout()
    File run_stderr = stderr()
    File session_info = "session_info.txt"
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
    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }
    log "Starting the joint mvSuSiE summary."
    mkdir -p window_outputs
    Rscript /opt/mvsusie/scripts/summarize_window.R \
      --prepared ~{prepared_window} \
      --fit ~{mvsusie_fit} \
      --output-dir window_outputs
    log "Verifying the joint mvSuSiE summary outputs."
    test -s window_outputs/variant_pip.tsv.gz
    test -e window_outputs/credible_sets.tsv.gz
    test -e window_outputs/credible_set_members.tsv.gz
    test -s window_outputs/component_feature_support.tsv.gz
    test -s window_outputs/window_qc.tsv
    log "Completed the joint mvSuSiE summary."
  >>>

  output {
    File variant_pip = "window_outputs/variant_pip.tsv.gz"
    File credible_sets = "window_outputs/credible_sets.tsv.gz"
    File credible_set_members = "window_outputs/credible_set_members.tsv.gz"
    File component_feature_support = "window_outputs/component_feature_support.tsv.gz"
    File window_qc = "window_outputs/window_qc.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-trans-window-mvsusie:latest"
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}

task PlotMvSusie {
  input {
    File prepared_window
    File mvsusie_fit
  }

  command <<<
    set -euo pipefail
    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }
    log "Starting the mvSuSiE credible-set-by-feature plot."
    Rscript /opt/mvsusie/scripts/plot_window_mvsusie.R \
      --prepared ~{prepared_window} \
      --fit ~{mvsusie_fit} \
      --png effect_plot.png \
      --pdf effect_plot.pdf \
      --plot-rds effect_plot.rds
    log "Verifying the mvSuSiE plot outputs."
    test -s effect_plot.png
    test -s effect_plot.pdf
    test -s effect_plot.rds
    log "Completed the mvSuSiE credible-set-by-feature plot."
  >>>

  output {
    File effect_plot_png = "effect_plot.png"
    File effect_plot_pdf = "effect_plot.pdf"
    File effect_plot_rds = "effect_plot.rds"
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
    Array[File] credible_set_members
    Array[File] component_feature_support
    Array[File] window_qc
  }

  command <<<
    set -euo pipefail
    log() {
      printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2
    }
    log "Starting the merge of joint mvSuSiE summaries."
    mkdir -p merged
    Rscript /opt/mvsusie/scripts/merge_window_outputs.R \
      --variant-pips "~{sep="," variant_pips}" \
      --credible-sets "~{sep="," credible_sets}" \
      --credible-set-members "~{sep="," credible_set_members}" \
      --component-feature-support "~{sep="," component_feature_support}" \
      --window-qc "~{sep="," window_qc}" \
      --output-dir merged
    log "Verifying the merged joint mvSuSiE summaries."
    test -s merged/variant_pip.tsv.gz
    test -e merged/credible_sets.tsv.gz
    test -e merged/credible_set_members.tsv.gz
    test -s merged/component_feature_support.tsv.gz
    test -s merged/window_qc.tsv
    log "Completed the merge of joint mvSuSiE summaries."
  >>>

  output {
    File merged_variant_pip = "merged/variant_pip.tsv.gz"
    File merged_credible_sets = "merged/credible_sets.tsv.gz"
    File merged_credible_set_members = "merged/credible_set_members.tsv.gz"
    File merged_component_feature_support = "merged/component_feature_support.tsv.gz"
    File merged_window_qc = "merged/window_qc.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-trans-window-mvsusie:latest"
    cpu: 1
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}
