# Joint RNA–Protein mvSuSiE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build one production workflow that jointly preprocesses expression, splicing, and protein outcomes and fits a correctly scaled mashr-informed mvSuSiE model with a 10, 15, 20, ... greedy schedule and cutoff 1.0.

**Architecture:** Keep the existing two-stage WDL structure: `PrepareTransWindow` selects locus genotypes and modality-specific outcomes, and `TransWindowMvSusie` aligns samples, residualizes the joint data, learns the mashr prior, fits mvSuSiE, and summarizes the final fit. Replace flexible modality arrays with explicit expression, splicing, and protein interfaces. Keep the fit object from only the final greedy round while writing a compact history for every round.

**Tech Stack:** WDL 1.0, R 4.4.1, data.table, dplyr, purrr, readr, optparse, mashr, mvsusieR, susieR, ggplot2, GitHub Actions, miniwdl, actionlint.

**Spec:** `docs/superpowers/specs/2026-08-24-joint-rna-protein-mvsusie-design.md`

## Global Constraints

- The only accepted modalities are `expression`, `splicing`, and `protein`; all three are required.
- Remove `isoform_usage` from production interfaces, examples, fixtures, and tests.
- Default trans counts are expression 25, splicing 25, and protein 15; each count remains independently configurable.
- Target-feature inputs may contain expression and splicing only.
- Phenotypes use modality-specific covariates; the genotype matrix uses the aligned union of every modality covariate.
- Same-name covariates with different values remain distinct and receive modality prefixes.
- Use all retained locus SNPs for the final mash mixture and strong rows at lfsr 0.05 for PCA covariance learning.
- Call `mashr::cov_pca` with `npc = 5`; do not use canonical covariances or extreme deconvolution.
- Preserve fitted mashr weights and the corrected covariance scale in mvSuSiE.
- Fix the residual covariance at the initial joint-outcome covariance.
- Greedy defaults are start 10, step 5, maximum 40, and minimum component log-BF cutoff 1.0.
- Save only the final mvSuSiE fit. Write one compact history row for each greedy round.
- Generate effect plots only through `mvsusieR::mvsusie_plot(..., conditional_effect = TRUE)`.
- Pin mvsusieR commit `ebd1133953005fa70c6b338727b5fe9222e2a1c2` and susieR commit `65f3586a865fb6748cb4f9df50510ac577706348`.
- Add logging to every WDL command and each long R preprocessing or model stage.
- Do not build Docker locally. GitHub Actions owns the image build and smoke test.

---

### Task 1: Joint phenotype selection and preparation interface

**Files:**
- Modify: `scripts/prepare_trans_window.R`
- Modify: `workflows/prepare_trans_window.wdl`
- Modify: `tests/fixtures/trans_window/generate_prepare_fixture.R`
- Modify: `tests/test_prepare_trans_window.R`
- Modify: `tests/test_prepare_trans_window_wdl.sh`

**Interfaces:**
- Consumes: trans associations with `window_id`, coordinates, `modality`, `molecular_trait_id`, and `p_value`; three BED-like phenotype files; target TSV with `window_id`, `modality`, and `phenotype_id`.
- Produces: `select_top_trans_phenotypes(trans_associations, top_n_by_modality)`, `read_target_phenotypes(path)`, and a combined manifest with `window_id`, `outcome_key`, `phenotype_id`, `modality`, and `phenotype_file`.

- [ ] **Step 1: Extend the preparation fixture with all required modalities and targets**

Update `generate_prepare_fixture.R` so it writes 27 expression, 27 splicing, and 17 protein trans-association rows. Write the corresponding BED-like phenotype rows, plus three target-only phenotype rows. Add `target_phenotypes.tsv` with one expression and two splicing targets. Give the three target IDs no trans-association row so their inclusion is independent of the top-N selection. Use distinct outcome IDs and sample values. The default selection must retain 25 trans expression, 25 trans splicing, 15 trans protein, and all target rows after de-duplication.

```r
modalities <- c("expression", "splicing", "protein")
trans_counts <- c(expression = 27L, splicing = 27L, protein = 17L)
targets <- tibble::tribble(
  ~window_id, ~modality, ~phenotype_id,
  "w1", "expression", "expr_target",
  "w1", "splicing", "splice_target_1",
  "w1", "splicing", "splice_target_2"
)
```

- [ ] **Step 2: Write failing joint-selection tests**

In `test_prepare_trans_window.R`, assert exact required modalities, independent default counts, target inclusion, and unique outcome keys. Add one focused case where a target is already in the top trans set and assert that it occurs once in the manifest. Add failures for a missing modality, `isoform_usage`, a missing target, and duplicate target rows.

```r
expected_counts <- c(expression = 26L, splicing = 27L, protein = 15L)
actual_counts <- table(factor(
  result$manifest$modality,
  levels = names(expected_counts)
))
stopifnot(identical(as.integer(actual_counts), unname(expected_counts)))
stopifnot(identical(
  result$manifest$outcome_key,
  paste(result$manifest$modality, result$manifest$phenotype_id, sep = "::")
))
```

Run: `Rscript tests/test_prepare_trans_window.R`

Expected: FAIL because the current selector accepts one shared `top_n`, has no target TSV, and has no `outcome_key`.

- [ ] **Step 3: Implement explicit joint selection**

Change the selector signature and validate exact named counts:

```r
required_joint_modalities <- function() c("expression", "splicing", "protein")

validate_joint_modalities <- function(modalities, label) {
  expected <- required_joint_modalities()
  actual <- sort(unique(as.character(modalities)))
  if (!identical(actual, sort(expected))) {
    stop(label, " must contain exactly: ", paste(expected, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

select_top_trans_phenotypes <- function(trans_associations, top_n_by_modality) {
  expected <- required_joint_modalities()
  validate_joint_modalities(trans_associations$modality, "Trans associations")
  if (
    !identical(sort(names(top_n_by_modality)), sort(expected)) ||
    any(!is.finite(top_n_by_modality)) ||
    any(top_n_by_modality < 1L) ||
    any(top_n_by_modality != as.integer(top_n_by_modality))
  ) {
    stop("top_n_by_modality must contain one positive integer per required modality.",
         call. = FALSE)
  }
  eligible <- trans_associations |>
    dplyr::group_by(.data$modality, .data$molecular_trait_id) |>
    dplyr::summarise(min_pval = min(.data$p_value), .groups = "drop") |>
    dplyr::arrange(.data$modality, .data$min_pval, .data$molecular_trait_id)
  available <- table(eligible$modality)
  if (any(available[expected] < top_n_by_modality[expected])) {
    stop("Each modality must contain at least its requested number of trans phenotypes.",
         call. = FALSE)
  }
  eligible |>
    dplyr::group_by(.data$modality) |>
    dplyr::group_modify(function(.x, .y) {
      dplyr::slice_head(.x, n = unname(top_n_by_modality[[.y$modality[[1L]]]]))
    }) |>
    dplyr::ungroup()
}
```

Read target rows, validate expression/splicing only, require every target exactly once, and form `outcome_key = paste(modality, phenotype_id, sep = "::")`.

- [ ] **Step 4: Replace WDL modality arrays with explicit phenotype inputs**

Change `PrepareTransWindow` and `PrepareWindowPhenotypes` inputs to:

```wdl
File expression_phenotypes
File splicing_phenotypes
File protein_phenotypes
File target_phenotypes
Int top_n_expression = 25
Int top_n_splicing = 25
Int top_n_protein = 15
```

Pass explicit CLI flags and add log lines before selection and after output validation. Remove `extract_cis_window_phenotypes`, `phenotype_files`, `phenotype_modalities`, and `top_n_trans_phenotypes`.

- [ ] **Step 5: Update WDL tests and run task tests**

Add exact `rg` checks for the seven new inputs and checks that the removed inputs are absent.

Run:

```bash
Rscript tests/test_prepare_trans_window.R
bash tests/test_prepare_trans_window_wdl.sh
miniwdl check workflows/prepare_trans_window.wdl
```

Expected: PASS with explicit counts 25, 25, and 15.

- [ ] **Step 6: Commit Task 1**

```bash
git add scripts/prepare_trans_window.R workflows/prepare_trans_window.wdl tests/fixtures/trans_window/generate_prepare_fixture.R tests/test_prepare_trans_window.R tests/test_prepare_trans_window_wdl.sh
git commit -m "feat: select joint RNA protein phenotypes"
```

---

### Task 2: Joint phenotype and covariate readers

**Files:**
- Modify: `scripts/trans_window_io.R`
- Modify: `scripts/prepare_window.R`
- Modify: `tests/fixtures/trans_window/generate_reader_fixture.R`
- Modify: `tests/fixtures/trans_window/generate_model_fixture.R`
- Modify: `tests/test_trans_window_r.R`
- Modify: `tests/test_trans_window_cli.R`

**Interfaces:**
- Consumes: manifest from Task 1 and explicit expression, splicing, and protein covariate files.
- Produces: `read_joint_covariates(expression_path, splicing_path, protein_path)` returning a named three-element list and `read_window_phenotypes(...)` returning unique `outcome_key` columns.

- [ ] **Step 1: Write failing reader fixtures and tests**

Remove isoform fixture generation. Add a protein phenotype table with the same four metadata columns and a protein covariate file. Give all three covariate files a `PC1` column with different aligned values and a `SHARED` column with identical values. Shuffle sample columns and covariate rows independently. Give the model fixture at least two outcomes per modality and at least 12 usable SNPs so it can exercise five PCA inputs and a first greedy round at `L = 10`.

Assert:

```r
stopifnot(identical(sort(unique(phenotype_data$modalities)), c("expression", "protein", "splicing")))
stopifnot(identical(names(covariates), c("expression", "splicing", "protein")))
stopifnot(any(grepl("^protein::", phenotype_data$phenotype_ids)))
```

Add negative tests for RNA-only input, `isoform_usage`, duplicate sample IDs, duplicate outcome keys, and a covariate file with no exact shared sample intersection.

Run: `bash tests/test_trans_window_r.sh`

Expected: FAIL because the existing reader rejects protein and accepts isoform usage.

- [ ] **Step 2: Restrict manifest and phenotype layouts**

Use one constant in `trans_window_io.R`:

```r
required_joint_modalities <- function() c("expression", "splicing", "protein")
```

Make `validate_phenotype_manifest` require exactly those modalities. Treat all three as four-metadata-column BED-like files. Allow a combined file to contain exactly those three layouts. Use `outcome_key` for matrix column names while preserving `phenotype_id` in metadata.

- [ ] **Step 3: Implement explicit joint covariate reading**

Add:

```r
read_joint_covariates <- function(expression_path, splicing_path, protein_path) {
  paths <- c(
    expression = expression_path,
    splicing = splicing_path,
    protein = protein_path
  )
  result <- lapply(paths, read_covariate_file)
  if (any(vapply(result, function(x) anyDuplicated(rownames(x)) > 0L, logical(1L)))) {
    stop("Joint covariate files must contain unique sample IDs.", call. = FALSE)
  }
  result
}
```

Remove `read_covariate_matrices`, shared-modality expansion, and isoform branches from the production call path.

- [ ] **Step 4: Replace preparation CLI arrays with explicit joint flags**

In `prepare_window.R`, replace comma-separated covariate inputs with:

```text
--expression-covariates
--splicing-covariates
--protein-covariates
```

Call `read_joint_covariates` and log each matrix dimension.

- [ ] **Step 5: Run reader and CLI tests**

Run:

```bash
Rscript tests/test_trans_window_cli.R
bash tests/test_trans_window_r.sh
```

Expected: PASS with exact joint modality validation and reordered samples.

- [ ] **Step 6: Commit Task 2**

```bash
git add scripts/trans_window_io.R scripts/prepare_window.R tests/fixtures/trans_window/generate_reader_fixture.R tests/fixtures/trans_window/generate_model_fixture.R tests/test_trans_window_r.R tests/test_trans_window_cli.R
git commit -m "feat: read required joint modalities"
```

---

### Task 3: Correct joint preprocessing and covariate provenance

**Files:**
- Modify: `scripts/trans_window_preprocess.R`
- Modify: `scripts/prepare_window.R`
- Modify: `scripts/run_window_mvsusie.R`
- Modify: `tests/test_trans_window_r.R`
- Modify: `tests/test_trans_window_r.sh`

**Interfaces:**
- Consumes: aligned joint phenotype data, locus dosage, and the named covariate list from Task 2.
- Produces: `prepare_joint_window_data(...)` and `covariate_provenance` with source modality, original/final names, ranks, and MD5 checksums.

- [ ] **Step 1: Write failing preprocessing regression tests**

Add tests that place one non-finite value in one phenotype row and one covariate row. Assert that only affected samples are removed. Add a zero-variance phenotype between retained columns and assert that modality labels remain aligned after filtering.

Verify orthogonality:

```r
for (modality in c("expression", "splicing", "protein")) {
  idx <- which(prepared$phenotype_metadata$modality == modality)
  C <- cbind(covariates[[modality]][prepared$samples, , drop = FALSE], intercept = 1)
  stopifnot(max(abs(crossprod(C, prepared$Y[, idx, drop = FALSE]))) < 1e-6)
}
G <- cbind(prepared$genotype_covariates, intercept = 1)
stopifnot(max(abs(crossprod(G, prepared$X))) < 1e-6)
stopifnot(max(abs(apply(prepared$Y, 2L, sd) - 1)) < 1e-10)
```

Run the focused R test and confirm failure at the current undefined `C` branch or modality-index mismatch.

- [ ] **Step 2: Align samples before covariate comparison**

Replace `make_genotype_covariates(modality_covariates)` with an interface that aligns samples before it compares columns:

```r
make_genotype_covariates <- function(modality_covariates, samples) {
  aligned <- lapply(modality_covariates, function(x) x[samples, , drop = FALSE])
  candidates <- Map(function(x, modality) {
    colnames(x) <- paste(modality, colnames(x), sep = "::")
    x
  }, aligned, names(aligned))
  collapse_aligned_covariates(candidates, source_modalities = names(aligned))
}
```

Implement `collapse_aligned_covariates` to compare numeric columns after the exact sample reorder. Preserve same-name columns when their values differ. De-duplicate columns only when their aligned numeric values are identical. Return `list(matrix = ..., provenance = ...)`, with one provenance row per source column and fields `source_modality`, `original_name`, `final_name`, `deduplicated_to`, `rank_before`, `rank_after`, and `md5`. Use `tools::md5sum` on a temporary uncompressed RDS serialization of each aligned numeric column to create stable checksums without a new dependency.

- [ ] **Step 3: Replace the incomplete-sample branch**

Build a row-wise mask:

```r
finite_by_row <- function(x) apply(is.finite(x), 1L, all)
complete <- finite_by_row(X_raw) & finite_by_row(Y_raw)
for (matrix in aligned_covariates) complete <- complete & finite_by_row(matrix)
if (!all(complete)) {
  pipeline_log(sprintf("Removing %d incomplete joint samples.", sum(!complete)))
  samples <- samples[complete]
  X_raw <- X_raw[complete, , drop = FALSE]
  Y_raw <- Y_raw[complete, , drop = FALSE]
  aligned_covariates <- lapply(aligned_covariates, function(x) x[complete, , drop = FALSE])
}
```

Delete the undefined `C <- C[complete, ]` statement.

- [ ] **Step 4: Keep phenotype metadata aligned through both filters**

Filter `phenotype_modalities`, `phenotype_ids`, and metadata immediately after each phenotype mask. Require at least one retained outcome per required modality after residualization.

Rename the public entry point to:

```r
prepare_joint_window_data(
  window,
  phenotype_data,
  dosage,
  covariates_by_modality,
  keep_samples = NULL,
  min_genotype_variance = 1e-8,
  min_phenotype_variance = 1e-8
)
```

Return `genotype_covariates` and `covariate_provenance` in the prepared bundle.

- [ ] **Step 5: Add stage logging and save provenance**

Log sample intersection, input reorder status, incomplete removal, rank-INT, each modality residualization, genotype union construction, final scaling, variant filtering, and retained dimensions. Write `covariate_provenance.tsv.gz` in the entry points. Compute MD5 checksums for every input file in `prepare_window.R` and `run_window_mvsusie.R`, then store them as `prepared$input_checksums` before saving the prepared bundle.

Assert in the test that the saved residualized genotype columns are not scaled to unit variance. Pass `standardize = TRUE` only in the mvSuSiE call so genotype standardization occurs once, inside mvSuSiE.

- [ ] **Step 6: Run preprocessing tests**

Run:

```bash
bash tests/test_trans_window_r.sh
```

Expected: PASS with different expression/splicing/protein `PC1` columns preserved and the identical `SHARED` column de-duplicated.

- [ ] **Step 7: Commit Task 3**

```bash
git add scripts/trans_window_preprocess.R scripts/prepare_window.R scripts/run_window_mvsusie.R tests/test_trans_window_r.R tests/test_trans_window_r.sh
git commit -m "fix: harmonize joint preprocessing"
```

---

### Task 4: Joint mashr prior and separate greedy schedule controls

**Files:**
- Modify: `scripts/trans_window_prior.R`
- Modify: `scripts/trans_window_model.R`
- Modify: `scripts/fit_window.R`
- Modify: `tests/test_mashr_prior.R`
- Modify: `tests/test_trans_window_r.R`
- Modify: `tests/test_trans_window_r.sh`

**Interfaces:**
- Consumes: prepared joint bundle from Task 3.
- Produces: `learn_joint_mashr_prior`, `fit_mvsusie_greedy_schedule`, final fit, mashr training bundle, and greedy history.

- [ ] **Step 1: Write failing joint-prior default tests**

Assert these individual fields from `make_model_config()` so unrelated metadata does not make the test brittle:

```r
config <- make_model_config()
stopifnot(
  identical(config$start_L, 10L),
  identical(config$step_L, 5L),
  identical(config$max_L, 40L),
  identical(config$greedy_lbf_cutoff, 1),
  identical(config$mashr_n_pca, 5L),
  identical(config$mashr_strong_lfsr, 0.05),
  identical(config$coverage, 0.95),
  identical(config$min_abs_corr, 0.5)
)
stopifnot(!"estimate_residual_variance" %in% names(config))
```

Assert that canonical, ED, prior-method, and residual-estimation switches are absent from the production config and CLI. Retain the existing numeric round-trip test for `prepare_mashr_prior_for_mvsusie`.

Run:

```bash
Rscript tests/test_mashr_prior.R
bash tests/test_trans_window_r.sh
```

Expected: FAIL because current defaults are canonical, ED on, residual estimation on, and cutoff 0.1.

- [ ] **Step 2: Make mashr PCA-only and joint-only**

Replace `learn_mashr_prior(..., use_extreme_deconvolution)` with:

```r
learn_joint_mashr_prior <- function(
  Bhat,
  Shat,
  n_pca = 5L,
  seed = NULL,
  strong_lfsr = 0.05,
  cov_pca_fun = mashr::cov_pca
) {
  mash_data <- make_mashr_data(Bhat, Shat)
  selection <- select_mashr_covariance_rows(
    mash_data,
    n_pca = n_pca,
    lfsr_threshold = strong_lfsr
  )
  pca_covariances <- cov_pca_fun(
    mash_data,
    npc = n_pca,
    subset = selection$rows
  )
  fit <- mashr::mash(
    data = mash_data,
    Ulist = pca_covariances,
    usepointmass = TRUE,
    outputlevel = 0,
    verbose = TRUE
  )
  build_joint_mashr_bundle(fit, selection, pca_covariances, seed, nrow(Bhat))
}
```

Keep `select_mashr_covariance_rows`: run `mashr::mash_1by1` on all rows, select rows with minimum lfsr at or below 0.05, and fall back to the five smallest minimum-lfsr rows when fewer than five pass. Implement `build_joint_mashr_bundle` to store the final fitted mash object, fitted weights, covariance matrices, selected row IDs, requested and returned PCA covariance counts, seed, and all-SNP/strong-row training counts. Test the exact `cov_pca_fun(..., npc = 5)` call with a stub so the test does not depend on PCA rank in the fixture.

Remove calls to `cov_ed` and canonical-prior creation from the production model path. Keep the blocked all-SNP matrix association code, its numerical fallback behavior, and one progress log per SNP block. Remove `write_marginal_association_table` and the `marginal_output` config/CLI path. Keep `Bhat` and `Shat` only in `mashr_training_bundle.rds`.

- [ ] **Step 3: Implement an independently tested greedy scheduler**

Add this interface:

```r
fit_mvsusie_greedy_schedule <- function(
  X,
  Y,
  prior,
  start_L = 10L,
  step_L = 5L,
  max_L = 40L,
  greedy_lbf_cutoff = 1,
  fit_fun = mvsusieR::mvsusie,
  ...
)
```

For each `requested_L`, call `fit_fun` with `model_init = previous_fit`, except in round one where it is `NULL`. Compute history fields from the returned fit. Stop when `min(fit$lbf) < greedy_lbf_cutoff` or at `max_L`. Keep only `previous_fit`, the current fit, and the history in memory.

Use a stub `fit_fun` test that records `L` and `model_init`, returns controlled LBFs `rep(2, 10)` and `c(rep(2, 14), 0.8)`, and verifies requests `10, 15`, warm start in round two, action `saturated`, and no list of intermediate fits in the result.

- [ ] **Step 4: Wire the joint fit with fixed covariance and weights**

In `fit_window_mvsusie`, always:

```r
residual_variance <- stats::cov(prepared$Y)
prepared_prior <- prepare_mashr_prior_for_mvsusie(
  mashr_training$raw_prior,
  prepared$Y
)
scheduled <- fit_mvsusie_greedy_schedule(
  X = prepared$X,
  Y = prepared$Y,
  prior = prepared_prior,
  start_L = config$start_L,
  step_L = config$step_L,
  max_L = config$max_L,
  greedy_lbf_cutoff = config$greedy_lbf_cutoff,
  residual_variance = residual_variance,
  estimate_residual_variance = FALSE,
  estimate_prior_variance = FALSE,
  estimate_prior_mixture_weights = FALSE,
  verbose = TRUE
)
```

Keep the scale conversion inside `prepare_mashr_prior_for_mvsusie`: calculate `outcome_se[r] = sd(Y[, r]) / sqrt(sum(is.finite(Y[, r])))`, then divide each covariance entry `U[r, s]` by `outcome_se[r] * outcome_se[s]`. Test the conversion and inverse round trip for unequal expression, splicing, and protein scales.

Return `fit`, `greedy_history`, `mashr_training`, and metadata. The mashr training object contains `Bhat`, `Shat`, raw prior, converted prior scale, outcome standard-error vector, raw and converted covariance ranges, fitted weights, selected row IDs, and full training counts.

- [ ] **Step 5: Replace CLI controls**

Remove `--L`, `--L-greedy`, `--prior-method`, `--mashr-skip-ed`, and `--fix-residual-variance`. Add:

```text
--start-L 10
--step-L 5
--max-L 40
--greedy-lbf-cutoff 1
--mashr-n-pca 5
--mashr-strong-lfsr 0.05
--mashr-seed
```

Validate `start_L <= max_L`, positive integer step, and a schedule that reaches a final value not greater than `max_L`.

Keep numeric overrides for `coverage` and `min_abs_corr`, with defaults 0.95 and 0.5. These values affect credible-set reporting only; they do not affect the greedy stopping rule.

- [ ] **Step 6: Run prior and model tests**

Run:

```bash
Rscript tests/test_mashr_prior.R
bash tests/test_trans_window_r.sh
```

Expected: PASS with PCA-only metadata, fixed weights/covariance, corrected scale, and history-only greedy output.

Add negative assertions for non-finite mash covariances, dimension mismatch between every covariance and `ncol(Y)`, non-finite mvSuSiE results, and a non-converged greedy round.

- [ ] **Step 7: Commit Task 4**

```bash
git add scripts/trans_window_prior.R scripts/trans_window_model.R scripts/fit_window.R tests/test_mashr_prior.R tests/test_trans_window_r.R tests/test_trans_window_r.sh
git commit -m "feat: fit joint mashr mvSuSiE schedule"
```

---

### Task 5: Final summaries and mvSuSiE API plot

**Files:**
- Modify: `scripts/trans_window_model.R`
- Modify: `scripts/summarize_window.R`
- Modify: `scripts/merge_window_outputs.R`
- Create: `scripts/plot_window_mvsusie.R`
- Modify: `tests/test_trans_window_r.sh`
- Create: `tests/test_plot_window_mvsusie.R`

**Interfaces:**
- Consumes: final fit bundle, prepared joint bundle, and greedy history from Task 4.
- Produces: PIPs, CS summaries/members, component-feature support, QC, final fit RDS, history TSV, and API PNG/PDF.

- [ ] **Step 1: Write failing summary-schema tests**

Extend the integration test to require these files and columns:

```r
required_support <- c(
  "component", "outcome_key", "modality", "phenotype_id",
  "single_effect_lfsr", "outcome_lbf"
)
required_history <- c(
  "round", "requested_L", "fitted_L", "niter", "minimum_lbf",
  "credible_set_count", "supported_component_count", "maximum_alpha", "action"
)
```

Create a unit test where a fake fit has only `mu2_diag`, then only legacy `mu2`, and assert that both produce identical uncertainty values.

Run: `Rscript tests/test_plot_window_mvsusie.R`

Expected: FAIL because the plot script and fallback helper do not exist.

- [ ] **Step 2: Correct the posterior second-moment fallback**

Use:

```r
get_fit_mu2 <- function(fit) {
  if (!is.null(fit$mu2_diag)) return(fit$mu2_diag)
  if (!is.null(fit$mu2)) return(fit$mu2)
  stop("The mvSuSiE fit contains neither mu2_diag nor mu2.", call. = FALSE)
}
```

Restrict large component-effect exports to credible-set members or sentinels. Do not expand all components × all variants × all outcomes.

- [ ] **Step 3: Implement final support and credible-set tables**

Write one row per component/outcome in `component_feature_support.tsv.gz`. Join `outcome_key` to modality and original phenotype ID from prepared metadata. Use `fit$single_effect_lfsr` and `fit$lbf_outcome`. Keep pairwise CS member rows in `credible_set_members.tsv.gz` and one aggregate row per set in `credible_sets.tsv.gz`. Keep coverage 0.95 and minimum absolute CS correlation 0.5 in the QC and CS tables.

- [ ] **Step 4: Implement API-only plotting**

Create `plot_window_mvsusie.R` with:

```r
plot_result <- mvsusieR::mvsusie_plot(
  fit = bundle$fit,
  chr = chromosome,
  pos = positions_mb,
  markers = bundle$fit$variable_names,
  outcomes = display_labels,
  lfsr_cutoff = 0.05,
  sentinel_only = FALSE,
  add_cs = TRUE,
  conditional_effect = TRUE,
  sort_by_cs = TRUE
)
stopifnot(!is.null(plot_result$effect_plot))
ggplot2::ggsave(png_path, plot_result$effect_plot, width = 15, height = 10, dpi = 300)
ggplot2::ggsave(pdf_path, plot_result$effect_plot, width = 15, height = 10)
```

Save the returned API result RDS for traceability, but do not calculate a replacement effect matrix.

- [ ] **Step 5: Update merge logic and run output tests**

Make merge inputs include CS summaries, CS members, component support, and QC. PIP and support tables must retain `window_id`; outcome tables must retain `modality` and `phenotype_id`.

Run:

```bash
Rscript tests/test_plot_window_mvsusie.R
bash tests/test_trans_window_r.sh
```

Expected: PASS and nonempty PNG/PDF files.

- [ ] **Step 6: Commit Task 5**

```bash
git add scripts/trans_window_model.R scripts/summarize_window.R scripts/merge_window_outputs.R scripts/plot_window_mvsusie.R tests/test_trans_window_r.sh tests/test_plot_window_mvsusie.R
git commit -m "feat: summarize and plot joint mvSuSiE fits"
```

---

### Task 6: Joint-only model WDL and container contract

**Files:**
- Modify: `scripts/run_window_mvsusie.R`
- Modify: `workflows/trans_window_mvsusie.wdl`
- Modify: `tests/test_trans_window_wdl.sh`
- Modify: `tests/test_trans_window_mvsusie_container.sh`
- Modify: `envs/trans-window-mvsusie.Dockerfile`
- Modify: `.github/workflows/trans-window-mvsusie-image.yml`

**Interfaces:**
- Consumes: prepared phenotype output from Task 1, locus dosage, and explicit three-modality covariates.
- Produces: all final artifacts from Tasks 3–5 through WDL outputs.

- [ ] **Step 1: Write failing WDL contract tests**

Assert the workflow contains explicit covariate files and defaults:

```bash
rg -q 'File expression_covariates' workflows/trans_window_mvsusie.wdl
rg -q 'File splicing_covariates' workflows/trans_window_mvsusie.wdl
rg -q 'File protein_covariates' workflows/trans_window_mvsusie.wdl
rg -q 'Int start_L = 10' workflows/trans_window_mvsusie.wdl
rg -q 'Int step_L = 5' workflows/trans_window_mvsusie.wdl
rg -q 'Int max_L = 40' workflows/trans_window_mvsusie.wdl
rg -q 'Float greedy_lbf_cutoff = 1.0' workflows/trans_window_mvsusie.wdl
```

Assert the old modality arrays, canonical prior, ED, isoform, and estimated-residual switches are absent.

Assert the workflow keeps `Float coverage = 0.95`, `Float min_abs_corr = 0.5`, and the existing CPU, memory, disk, and retry runtime inputs.

Run: `bash tests/test_trans_window_wdl.sh`

Expected: FAIL on the current flexible contract.

- [ ] **Step 2: Replace model WDL inputs and commands**

Pass explicit covariate files and the new numeric controls to `run_window_mvsusie.R`. Add WDL log lines for task start, all resolved defaults, model completion, summary completion, plot completion, and each output verification. Return prepared bundle, mashr bundle, final fit, history, summaries, provenance, logs, session information, and plots.

- [ ] **Step 3: Pin exact dependency commits**

Replace master tarballs with immutable URLs:

```dockerfile
https://github.com/stephenslab/susieR/archive/65f3586a865fb6748cb4f9df50510ac577706348.tar.gz
https://github.com/stephenslab/mvsusieR/archive/ebd1133953005fa70c6b338727b5fe9222e2a1c2.tar.gz
```

Ensure the model image contains `ggplot2` and copies `plot_window_mvsusie.R`.

- [ ] **Step 4: Extend the GitHub Actions smoke test**

Without running Docker locally, update the workflow’s container test to verify:

```r
stopifnot(all(c("start_L", "step_L", "max_L") %in% names(formals(make_model_config))))
stopifnot("protein" %in% required_joint_modalities())
stopifnot(all(file.exists(file.path("/opt/mvsusie/scripts", required_scripts))))
```

Run a small in-container joint fixture through preparation, mashr, one saturated greedy round, summary, and API plot creation.

- [ ] **Step 5: Run local WDL/static checks**

Run:

```bash
miniwdl check workflows/trans_window_mvsusie.wdl
bash tests/test_trans_window_wdl.sh
bash tests/test_trans_window_mvsusie_container.sh
actionlint .github/workflows/trans-window-mvsusie-image.yml
```

Expected: PASS. Do not run `docker build`.

- [ ] **Step 6: Commit Task 6**

```bash
git add scripts/run_window_mvsusie.R workflows/trans_window_mvsusie.wdl tests/test_trans_window_wdl.sh tests/test_trans_window_mvsusie_container.sh envs/trans-window-mvsusie.Dockerfile .github/workflows/trans-window-mvsusie-image.yml
git commit -m "feat: expose joint mvSuSiE WDL"
```

---

### Task 7: Documentation, complete verification, and GitHub handoff

**Files:**
- Modify: `README.md`
- Modify: `docs/trans-window-fine-mapping.md`
- Verify: `.dockstore.yml`
- Modify: `tests/test_trans_window_r.sh`
- Modify: `tests/test_trans_window_wdl.sh`
- Create: `tests/test_joint_workflow_requirements.sh`

**Interfaces:**
- Consumes: completed joint preparation and model workflows.
- Produces: user documentation, Dockstore metadata, verified branch, and pull-request-ready changes.

- [ ] **Step 1: Update documentation and examples**

Document the exact three modalities, explicit phenotype/covariate inputs, target TSV schema, per-modality top counts, sample/covariate residualization, mashr training scopes, corrected prior scale, fixed residual covariance, greedy defaults, final-only fit behavior, output paths, and API plot scale. Remove RNA-only and isoform examples.

- [ ] **Step 2: Update Dockstore metadata**

Keep the existing workflow names and descriptor paths. Validate that every `.dockstore.yml` descriptor path exists after the interface change. This repository has no Dockstore test parameter files, so do not add an unused parameter-file reference.

- [ ] **Step 3: Run the complete non-Docker suite**

Run:

```bash
bash tests/test_trans_window_r.sh
bash tests/test_trans_window_wdl.sh
Rscript tests/test_plot_window_mvsusie.R
Rscript tools/lint_r.R
miniwdl check workflows/prepare_trans_window.wdl
miniwdl check workflows/trans_window_mvsusie.wdl
git diff --check
```

Expected: all commands exit zero. Record any skipped real-data test and its required environment variable.

- [ ] **Step 4: Verify requirements against the approved spec**

Create `tests/test_joint_workflow_requirements.sh`. Use exact `rg` assertions and focused R assertions to cover: required modalities, default selection counts, target handling, separate phenotype covariates, union genotype covariates, finite-sample masking, PCA-only mashr, corrected scale conversion, fixed weights and residual covariance, the 10/5/40/1 greedy defaults, final-only fit persistence, API plotting, both dependency SHAs, WDL logging, and the GitHub Actions container smoke test. Run:

```bash
bash tests/test_joint_workflow_requirements.sh
```

Expected: PASS with each approved design heading represented by at least one assertion.

- [ ] **Step 5: Commit documentation and verification updates**

```bash
git add README.md docs/trans-window-fine-mapping.md tests/test_trans_window_r.sh tests/test_trans_window_wdl.sh tests/test_joint_workflow_requirements.sh
git commit -m "docs: document joint RNA protein workflow"
```

- [ ] **Step 6: Push and open a pull request after user approval**

Use the branch `codex/joint-rna-protein-workflow`. The PR summary must list the breaking joint-only interface, protein preprocessing, statistical defaults, preprocessing fixes, dependency pins, and tests. Monitor R lint, WDL validation, and the GitHub container smoke test until completion.
