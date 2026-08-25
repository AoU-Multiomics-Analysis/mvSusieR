# Single-window Staged mvSuSiE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run one window per workflow invocation, optionally bypass preprocessing with one prepared RDS, and align phenotype assays to their shared samples.

**Architecture:** `PrepareWindowPhenotypes` will align contributing assays before it writes the combined phenotype file. `TransWindowMvSusie` will validate one of two input modes, conditionally run a 16 GiB preparation task, and always run a separate 8 GiB fitting task from one prepared RDS. Summary and plotting remain per-window scalar operations.

**Tech Stack:** WDL 1.0, R 4.4.1, data.table, readr, dplyr, mashr, mvsusieR, MiniWDL, Bash, GitHub Actions

**Spec:** `docs/superpowers/specs/2026-08-24-single-window-staged-mvsusie-design.md`

## Global Constraints

- One external workflow invocation processes one genomic window.
- `prepared_window` is one optional `File?`, not an array.
- Preprocessing requests 16 GiB; fitting requests 8 GiB.
- One prepared RDS is written per window; greedy-L rounds do not write separate fit RDS files.
- Phenotype intersection uses only modalities that contribute outcomes.
- Sample IDs use the downstream normalization rule.
- The mashr prior, greedy-L algorithm, statistical transformations, summary scale, and mvsusieR plotting API do not change.
- WDL commands include timestamped logging.
- Local smoke tests do not build containers; GitHub Actions builds and tests images.

---

### Task 1: Align phenotype samples across contributing assays

**Files:**
- Modify: `tests/fixtures/trans_window/generate_prepare_fixture.R`
- Modify: `tests/test_prepare_trans_window.R`
- Modify: `scripts/prepare_trans_window.R`

**Interfaces:**
- Produces: `align_prepare_phenotype_samples(selected_tables)` returning `list(tables, qc, shared_samples)`.
- Produces QC columns: `n_input_samples`, `n_shared_samples`, and `n_samples_removed`.
- Preserves: `write_prepare_phenotype_subset(selected_tables, output_dir)` as the file writer for already aligned tables.

- [ ] **Step 1: Write fixtures with different sample membership and order**

Change `make_phenotypes` to accept a named sample vector and generate these headers:

```r
expression_samples <- c("X1001", "1002", "1003")
splicing_samples <- c("1003", "X1001", "1004")
protein_samples <- c("1005", "1003", "1001")
```

Write deterministic values from the phenotype row index plus the sample
position so that the test can verify reordering, not only column names.

- [ ] **Step 2: Write failing alignment tests**

In `tests/test_prepare_trans_window.R`, assert:

```r
combined <- read_tsv(result$phenotype_data, show_col_types = FALSE)
stopifnot(identical(names(combined)[-(1:4)], c("1001", "1003")))
stopifnot(identical(as.integer(qc$n_input_samples), c(3L, 3L, 3L)))
stopifnot(identical(as.integer(qc$n_shared_samples), c(2L, 2L, 2L)))
stopifnot(identical(as.integer(qc$n_samples_removed), c(1L, 1L, 1L)))
```

Also hand-check one expression, splicing, and protein row against literal
fixture values. Add a fixture with disjoint contributing sample sets and
assert an error matching `no shared phenotype samples`. Capture messages and
assert that they report the input, shared, and removed counts.

- [ ] **Step 3: Run the test and verify the intended failure**

Run:

```bash
Rscript tests/test_prepare_trans_window.R
```

Expected: FAIL because the current writer requires identical columns and does
not write the three sample-count QC fields.

- [ ] **Step 4: Implement sample normalization and alignment**

Add:

```r
normalize_prepare_sample_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  sub("^X(?=[0-9])", "", ids, perl = TRUE)
}

align_prepare_phenotype_samples <- function(selected_tables) {
  sample_columns <- lapply(selected_tables, function(data) {
    names(data)[!startsWith(names(data), ".")][-(1:4)]
  })
  normalized <- lapply(sample_columns, normalize_prepare_sample_ids)
  if (any(vapply(normalized, anyDuplicated, integer(1L)) > 0L)) {
    stop("Sample ID normalization creates duplicates within a modality.", call. = FALSE)
  }
  shared <- Reduce(intersect, normalized)
  shared <- normalized[[1L]][normalized[[1L]] %in% shared]
  if (!length(shared)) {
    stop("Contributing modalities have no shared phenotype samples.", call. = FALSE)
  }
  aligned <- Map(function(data, source_columns, ids) {
    metadata_columns <- names(data)[seq_len(4L)]
    internal_columns <- names(data)[startsWith(names(data), ".")]
    selected_columns <- source_columns[match(shared, ids)]
    output <- data[, c(metadata_columns, selected_columns, internal_columns)]
    names(output)[seq.int(5L, 4L + length(shared))] <- shared
    output
  }, selected_tables, sample_columns, normalized)
  qc <- tibble(
    modality = names(selected_tables),
    n_input_samples = lengths(normalized),
    n_shared_samples = length(shared),
    n_samples_removed = lengths(normalized) - length(shared)
  )
  list(tables = aligned, qc = qc, shared_samples = shared)
}
```

Call it after empty modalities are removed and before
`write_prepare_phenotype_subset`. Join its QC values into `window_qc.tsv` and
emit one `prepare_log()` message per contributing modality.

- [ ] **Step 5: Run preparation tests**

Run:

```bash
Rscript tests/test_prepare_trans_window.R
bash tests/test_prepare_trans_window_wdl.sh
Rscript tools/lint_r.R
```

Expected: all pass and R lint reports `R lint ok`.

- [ ] **Step 6: Commit**

```bash
git add tests/fixtures/trans_window/generate_prepare_fixture.R \
  tests/test_prepare_trans_window.R scripts/prepare_trans_window.R
git commit -m "fix: intersect phenotype preparation samples"
```

---

### Task 2: Validate and fit one prepared window

**Files:**
- Modify: `scripts/trans_window_io.R`
- Modify: `scripts/fit_window.R`
- Modify: `tests/test_trans_window_r.R`
- Modify: `tests/test_trans_window_r.sh`

**Interfaces:**
- Produces: `validate_prepared_window(prepared, expected_window_id)` returning the validated object invisibly.
- `fit_window.R` consumes `--prepared`, `--window-id`, and `--covariate-provenance-output`.
- `fit_window.R` continues to produce the final fit, mashr training bundle, and one greedy-L history.

- [ ] **Step 1: Write failing prepared-bundle validation tests**

Add literal test bundles to `tests/test_trans_window_r.R` and assert:

```r
expect_error_matching(
  validate_prepared_window(list(X = matrix(1)), "w1"),
  "missing required fields"
)
expect_error_matching(
  validate_prepared_window(valid_prepared, "wrong_window"),
  "does not match requested window"
)
stopifnot(inherits(validate_prepared_window(valid_prepared, "w1"), "list"))
```

Update the entrypoint fixture command in `tests/test_trans_window_r.sh` to pass
`--window-id w1` and `--covariate-provenance-output`. Assert that the
provenance file is nonempty.

- [ ] **Step 2: Run tests and verify the intended failure**

Run:

```bash
bash tests/test_trans_window_r.sh
```

Expected: FAIL because `validate_prepared_window` and the new fit CLI options
do not exist.

- [ ] **Step 3: Implement prepared-bundle validation**

Add this contract to `scripts/trans_window_io.R`:

```r
validate_prepared_window <- function(prepared, expected_window_id) {
  required <- c(
    "window", "X", "Y", "variant_metadata", "phenotype_metadata",
    "samples", "covariate_provenance", "covariate_rank",
    "phenotype_covariate_rank", "qc"
  )
  missing <- setdiff(required, names(prepared))
  if (length(missing)) {
    stop("Prepared window is missing required fields: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  actual_window_id <- as.character(prepared$window$window_id)
  if (length(actual_window_id) != 1L || actual_window_id != expected_window_id) {
    stop("Prepared window ID does not match requested window: ",
         expected_window_id, call. = FALSE)
  }
  if (!is.matrix(prepared$X) || !is.matrix(prepared$Y) ||
      nrow(prepared$X) != nrow(prepared$Y)) {
    stop("Prepared X and Y must be aligned matrices.", call. = FALSE)
  }
  prepared
}
```

- [ ] **Step 4: Update the fit entry point**

Source `trans_window_io.R`. Add `--window-id` and
`--covariate-provenance-output`. Immediately after `readRDS`, call:

```r
prepared <- validate_prepared_window(
  prepared,
  require_cli_arg(args, "window_id")
)
```

Write `prepared$covariate_provenance` with `data.table::fwrite`, verify the
file is nonempty, and log its path. Do not copy `prepared$X` or `prepared$Y`.

- [ ] **Step 5: Run model entrypoint tests**

Run:

```bash
bash tests/test_trans_window_r.sh
Rscript tools/lint_r.R
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add scripts/trans_window_io.R scripts/fit_window.R \
  tests/test_trans_window_r.R tests/test_trans_window_r.sh
git commit -m "feat: validate prepared window fits"
```

---

### Task 3: Convert the model WDL to a staged single-window workflow

**Files:**
- Modify: `workflows/trans_window_mvsusie.wdl`
- Modify: `tests/test_trans_window_wdl.sh`
- Modify: `tests/test_joint_workflow_requirements.sh`

**Interfaces:**
- Workflow input: `String window_id` and `File? prepared_window`.
- Workflow input: `String docker_image` with the current GHCR image as its default.
- Raw-mode optional inputs: `File? window_manifest`, `File? window_phenotypes_tsv`, `File? dosage`, `File? phenotype_data`, and three `File?` covariates.
- Workflow output: scalar `File prepared_window` plus scalar fit, summary, and plot files.
- Task output: `ValidateMvSusieInputs.run_preparation` controls the conditional preparation call.

- [ ] **Step 1: Write failing WDL contract tests**

Replace array/scatter assertions with:

```bash
rg -Fq 'File? prepared_window' "$workflow"
rg -Fq 'task ValidateMvSusieInputs' "$workflow"
rg -Fq 'task PrepareMvSusieInput' "$workflow"
rg -Fq 'task FitMvSusie' "$workflow"
rg -Fq 'memory: "8 GiB"' "$workflow"
rg -Fq 'memory: "16 GiB"' "$workflow"
if rg -q 'scatter[[:space:]]*[(]|Array\[File\] prepared_windows' "$workflow"; then
  exit 1
fi
```

Assert that all workflow result outputs are scalar `File` declarations and
that the merge call is absent.

- [ ] **Step 2: Run WDL tests and verify the intended failure**

Run:

```bash
bash tests/test_trans_window_wdl.sh
```

Expected: FAIL because the current workflow scatters and has one combined
16 GiB `RunMvSusie` task.

- [ ] **Step 3: Implement `ValidateMvSusieInputs`**

Give the task Boolean inputs for `has_prepared_window` and each required raw
file. Its logged Bash command must write `false` to `run_preparation.txt` in
prepared mode. In raw mode, collect missing input labels in a Bash array,
print them, exit nonzero when the array is nonempty, and otherwise write
`true`. Declare:

```wdl
output {
  Boolean run_preparation = read_boolean("run_preparation.txt")
}
```

Use 1 GiB memory and one CPU because the task receives no files.

- [ ] **Step 4: Implement conditional preparation**

Define `PrepareMvSusieInput` with the existing raw inputs and call
`prepare_window.R`. Request 16 GiB and verify `prepared_window.rds` before the
task exits.

At workflow level:

```wdl
if (ValidateMvSusieInputs.run_preparation) {
  call PrepareMvSusieInput {
    input:
      window_id = window_id,
      window_manifest = select_first([window_manifest]),
      window_phenotypes_tsv = select_first([window_phenotypes_tsv]),
      dosage = select_first([dosage]),
      phenotype_data = select_first([phenotype_data]),
      expression_covariates = select_first([expression_covariates]),
      splicing_covariates = select_first([splicing_covariates]),
      protein_covariates = select_first([protein_covariates])
  }
}

File resolved_prepared_window = select_first([
  prepared_window,
  PrepareMvSusieInput.prepared_window
])
```

Pass the optional keep-samples input and preprocessing thresholds through the
conditional task.

- [ ] **Step 5: Implement the 8 GiB fitting task and scalar consumers**

`FitMvSusie` must call `fit_window.R` with the resolved prepared RDS,
`--window-id`, model controls, `--covariate-provenance-output`, mashr output,
greedy history, and final fit output. Preserve verbose output. Request 8 GiB.
Pass the workflow `docker_image` input to every task and use it in each
runtime block. This permits GitHub Actions to run the WDL against the image
built in the same job while production keeps the GHCR default.

Call `SummarizeMvSusie` and `PlotMvSusie` directly with scalar files. Remove
the workflow scatter, merge call, and merged outputs. Keep all existing
per-window summary and plot outputs as scalar `File` values.

- [ ] **Step 6: Validate WDL behavior**

Run:

```bash
miniwdl check workflows/trans_window_mvsusie.wdl
bash tests/test_trans_window_wdl.sh
bash tests/test_joint_workflow_requirements.sh
```

Expected: all pass with no scatter and no array output.

- [ ] **Step 7: Commit**

```bash
git add workflows/trans_window_mvsusie.wdl \
  tests/test_trans_window_wdl.sh tests/test_joint_workflow_requirements.sh
git commit -m "feat: stage single-window mvsusie workflow"
```

---

### Task 4: Package and document both execution modes

**Files:**
- Modify: `envs/trans-window-mvsusie.Dockerfile`
- Modify: `.github/workflows/trans-window-mvsusie-image.yml`
- Modify: `tests/test_trans_window_mvsusie_container.sh`
- Modify: `tests/fixtures/trans_window/generate_model_fixture.R`
- Modify: `docs/trans-window-fine-mapping.md`
- Modify: `README.md`

**Interfaces:**
- Container includes both `prepare_window.R` and `fit_window.R`.
- GitHub smoke tests execute the WDL in raw-input and prepared-input modes.
- Documentation shows raw-input and prepared-input WDL examples.

- [ ] **Step 1: Write failing packaging checks**

Require `prepare_window.R` in both the Dockerfile copy list and GitHub script
existence list. Remove the old assertion that forbids it. Require the workflow
documentation to contain `prepared_window`, `16 GiB`, and `8 GiB`.

- [ ] **Step 2: Run packaging tests and verify the intended failure**

Run:

```bash
bash tests/test_trans_window_mvsusie_container.sh
```

Expected: FAIL because the model image currently excludes `prepare_window.R`.

- [ ] **Step 3: Update the container and GitHub smoke test**

Add `scripts/prepare_window.R` to the image. Add it to the GitHub existence
check. Update `generate_model_fixture.R` to also write one combined
`model_phenotypes.tsv` and make its phenotype manifest refer to that file.
Ensure `tests/test_trans_window_r.sh` runs the staged prepare and fit entry
points and repeats the fit with the saved prepared RDS.

Keep container building inside GitHub Actions. Install MiniWDL in the image
workflow job after the local image is built. Generate the small model fixture
and run `TransWindowMvSusie` once with raw inputs and
`docker_image=trans-window-mvsusie-ci:latest`. Save its output JSON, extract
the scalar prepared RDS path with `jq`, and run the WDL again with only
`window_id`, `prepared_window`, model controls, and the local image tag. Check
that both invocations produce nonempty fit, summary, and plot outputs. Do not
build the image locally.

- [ ] **Step 4: Update documentation**

Document:

```text
Raw mode: raw window inputs -> PrepareMvSusieInput (16 GiB) -> FitMvSusie (8 GiB)
Prepared mode: prepared_window.rds --------------------------> FitMvSusie (8 GiB)
```

State that the workflow is one window per invocation, phenotype preparation
intersects contributing assay samples, and cross-window merge is an outer
aggregation step.

- [ ] **Step 5: Run packaging and documentation checks**

Run:

```bash
bash tests/test_trans_window_mvsusie_container.sh
bash tests/test_joint_workflow_requirements.sh
git diff --check
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add envs/trans-window-mvsusie.Dockerfile \
  .github/workflows/trans-window-mvsusie-image.yml \
  tests/test_trans_window_mvsusie_container.sh \
  tests/fixtures/trans_window/generate_model_fixture.R \
  docs/trans-window-fine-mapping.md README.md
git commit -m "docs: expose staged mvsusie execution modes"
```

---

### Task 5: Verify and publish

**Files:**
- Verify all files changed in Tasks 1 through 4.

**Interfaces:**
- Produces: a fast-forward commit series suitable for `origin/main`.
- Produces: GitHub Actions evidence for WDL, R lint, and image smoke tests.

- [ ] **Step 1: Run the complete local verification suite**

Run:

```bash
Rscript tools/lint_r.R
bash tests/test_trans_window_r.sh
bash tests/test_trans_window_wdl.sh
bash tests/test_joint_workflow_requirements.sh
git diff --check
git status --short
```

Expected: all tests pass, lint reports `R lint ok`, diff check is empty, and
status contains only the expected committed branch state.

- [ ] **Step 2: Review the final commit range**

Run:

```bash
git log --oneline origin/main..HEAD
git diff --stat origin/main...HEAD
```

Expected: the design, phenotype alignment, prepared fit, staged WDL, and
documentation commits only.

- [ ] **Step 3: Push to main**

Run:

```bash
git push origin HEAD:main
```

Expected: a fast-forward update of `main`.

- [ ] **Step 4: Verify GitHub Actions**

Run:

```bash
gh run list --branch main --limit 10 \
  --json name,status,conclusion,url,headSha
```

Wait for R lint, WDL validation, the prepare-window phenotype image, and the
trans-window mvSuSiE image smoke test. Report any failure with its workflow
URL and failing step. Do not claim completion until these relevant runs pass.
