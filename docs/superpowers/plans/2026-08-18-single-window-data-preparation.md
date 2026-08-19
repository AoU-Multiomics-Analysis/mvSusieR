# Single-Window mvSuSiE Data Preparation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a non-scattered WDL workflow and tidyverse R entrypoint that prepare one mvSuSiE window bundle per invocation from a full windows table, tabix-indexed genome-wide dosage, gzipped phenotype BED files, and one window’s trans associations.

**Architecture:** The caller supplies exactly one `window_id`; an independent genotype task uses `tabix` to write the dosage subset and one-row dosage manifest, while an independent phenotype task resolves the window, reads each phenotype BED fully into memory, retains matching trans phenotypes and optionally all coordinate-overlapping cis phenotypes, and writes one combined phenotype BED plus a modality key and QC table. Neither task scatters over the windows table.

**Tech Stack:** WDL 1.1, R, tidyverse (`dplyr`, `purrr`, `readr`, `stringr`, `tibble`), `optparse`, `tabix`.

**Spec:** `docs/superpowers/specs/2026-08-18-trans-window-mv-finemapping-design.md`

## Global Constraints

- One workflow invocation handles one `window_id`; no WDL scatter is used.
- Window coordinates are 0-based, half-open intervals.
- Phenotype files are gzip-compressed tab-delimited BED-like files and are not assumed to be indexed.
- Expression and splicing phenotype IDs are in column 4; the first three columns are chromosome, start, and end.
- `extract_cis_window_phenotypes` defaults to `true`.
- Trans association inputs contain `phenotype_id` and `modality` and are already restricted to the requested window.
- Output phenotype manifests retain the existing `window_id`, `phenotype_id`, `modality`, and `phenotype_file` interface used by `prepare_window.R`.

### Task 1: Add the failing single-window extraction tests

**Files:**
- Create: `tests/test_prepare_trans_window.R`
- Create: `tests/fixtures/trans_window/generate_prepare_fixture.R`

**Interfaces:**
- Consumes: the planned `prepare_trans_window_data()` function and its standard BED assumptions.
- Produces: executable tests for trans-only selection, cis overlap selection, compressed input handling, combined phenotype output, and one-window manifest output.

- [ ] **Step 1: Write the failing test**

  Generate a two-window fixture with gzipped expression and splicing files, a per-window trans association file, and a dosage subset file. Assert that `prepare_trans_window_data()` for `w1` writes only `w1` outputs, keeps trans IDs even when outside the cis coordinates, adds coordinate-overlapping cis IDs when enabled, excludes non-overlapping rows, and writes a manifest pointing to the subset files.

- [ ] **Step 2: Run the test to verify it fails**

  Run `Rscript tests/test_prepare_trans_window.R` from the repository root. It should fail because `scripts/prepare_trans_window.R` does not yet exist.

### Task 2: Implement the tidyverse extractor

**Files:**
- Create: `scripts/prepare_trans_window.R`

**Interfaces:**
- `read_prepare_window_manifest(path)` returns a validated tibble with one row per `window_id`.
- `prepare_trans_window_data(windows, window_id, trans_associations, phenotype_inputs, output_dir, extract_cis_window_phenotypes = TRUE)` returns the combined phenotype BED, modality key, and QC paths for the phenotype bundle.
- CLI accepts `--windows`, `--window-id`, `--trans-associations`, `--phenotype-files`, `--phenotype-modalities`, `--output-dir`, and `--extract-cis-window-phenotypes`.

- [ ] **Step 1: Implement validation and selection helpers**
- [ ] **Step 2: Implement full-table phenotype reading and trans/cis filtering**
- [ ] **Step 3: Implement subset BED, manifest, dosage, and QC writers**
- [ ] **Step 4: Add the optparse CLI entrypoint**
- [ ] **Step 5: Run the focused test and make it pass**

### Task 3: Add the non-scattered WDL workflow

**Files:**
- Create: `workflows/prepare_trans_window.wdl`
- Create: `tests/test_prepare_trans_window_wdl.sh`

**Interfaces:**
- Workflow `PrepareTransWindow` accepts one `windows_tsv`, one `window_id`, one `genome_dosage` plus its tabix index, one `trans_window_associations`, and parallel `phenotype_files`/`phenotype_modalities` arrays.
- Genotype outputs are `window_dosage` and `window_manifest`; phenotype outputs are `window_phenotypes`, `phenotype_data`, and `window_qc` for that one window.

- [ ] **Step 1: Define the single-call workflow with independent genotype and phenotype calls**
- [ ] **Step 2: Extract the requested genomic interval with tabix in `PrepareWindowGenotypes`**
- [ ] **Step 3: Invoke the R phenotype extractor once with the requested window ID in `PrepareWindowPhenotypes`**
- [ ] **Step 4: Validate the WDL with MiniWDL**

### Task 4: Integrate focused regression coverage and verify

**Files:**
- Modify: `tests/test_trans_window_r.sh`
- Modify: `tests/test_trans_window_wdl.sh`

- [ ] **Step 1: Add the new R and WDL tests to the existing harnesses**
- [ ] **Step 2: Run the focused R test, existing R harness, R lint, and WDL validation**
- [ ] **Step 3: Review the diff and confirm no user test data is staged**
