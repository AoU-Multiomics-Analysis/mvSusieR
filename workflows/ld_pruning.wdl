version 1.0

workflow LDPruning {
  input {
    File input_pgen
    File input_pvar
    File input_psam
    String output_prefix = "ld_pruning"
    File? keep_samples
    File? exclude_variants
    Float maf = 0.05
    Float geno = 0.01
    Int ld_window_kb = 1000
    Int ld_step_variants = 50
    Float ld_r2 = 0.1
  }

  call LDPruningTask {
    input:
      input_pgen = input_pgen,
      input_pvar = input_pvar,
      input_psam = input_psam,
      output_prefix = output_prefix,
      keep_samples = keep_samples,
      exclude_variants = exclude_variants,
      maf = maf,
      geno = geno,
      ld_window_kb = ld_window_kb,
      ld_step_variants = ld_step_variants,
      ld_r2 = ld_r2
  }

  output {
    File pruned_variants = LDPruningTask.pruned_variants
    File excluded_variants = LDPruningTask.excluded_variants
    File pruned_pgen = LDPruningTask.pruned_pgen
    File pruned_pvar = LDPruningTask.pruned_pvar
    File pruned_psam = LDPruningTask.pruned_psam
    File retained_variant_ids = LDPruningTask.retained_variant_ids
  }
}

task LDPruningTask {
  input {
    File input_pgen
    File input_pvar
    File input_psam
    String output_prefix
    File? keep_samples
    File? exclude_variants
    Float maf
    Float geno
    Int ld_window_kb
    Int ld_step_variants
    Float ld_r2
  }

  command <<<
    set -euo pipefail

    cp ~{input_pgen} input.pgen
    cp ~{input_pvar} input.pvar
    cp ~{input_psam} input.psam

    plink2 \
      --pfile input \
      --snps-only just-acgt \
      --max-alleles 2 \
      --maf ~{maf} \
      --geno ~{geno} \
      ~{if defined(keep_samples) then "--keep " + select_first([keep_samples]) else ""} \
      ~{if defined(exclude_variants) then "--exclude " + select_first([exclude_variants]) else ""} \
      --indep-pairwise ~{ld_window_kb}kb ~{ld_step_variants} ~{ld_r2} \
      --out ~{output_prefix}

    test -s ~{output_prefix}.prune.in

    plink2 \
      --pfile input \
      ~{if defined(keep_samples) then "--keep " + select_first([keep_samples]) else ""} \
      --extract ~{output_prefix}.prune.in \
      --make-pgen \
      --write-snplist \
      --out ~{output_prefix}.pruned
  >>>

  output {
    File pruned_variants = output_prefix + ".prune.in"
    File excluded_variants = output_prefix + ".prune.out"
    File pruned_pgen = output_prefix + ".pruned.pgen"
    File pruned_pvar = output_prefix + ".pruned.pvar"
    File pruned_psam = output_prefix + ".pruned.psam"
    File retained_variant_ids = output_prefix + ".pruned.snplist"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-ld-pruning:latest"
    cpu: 2
    memory: "64 GiB"
  }
}
