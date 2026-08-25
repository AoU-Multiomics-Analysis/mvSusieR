# Indexed Phenotype Access Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add reusable BGZF/tabix phenotype indexes and use them to retrieve selected expression, splicing, and protein rows without reading each complete wide phenotype matrix for every window.

**Architecture:** A one-time WDL workflow converts one BED-like phenotype file into a coordinate-sorted BGZF file, a tabix index, a phenotype-coordinate lookup, and QC. `PrepareTransWindow` accepts optional index and lookup companions for each modality. The R preparation entrypoint uses one batched tabix query when both companions are present and otherwise uses the existing full-scan reader.

**Tech Stack:** WDL 1.0, Bash, GNU coreutils and sort, BGZF, HTSlib tabix, R 4.4, dplyr, purrr, readr, tibble, MiniWDL, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-08-24-indexed-phenotype-access-design.md`

## Global Constraints

- Use phenotype-file coordinates. Do not add a GTF dependency.
- Keep the existing top-N, target selection, sample intersection, manifest, and combined phenotype output behavior.
- Keep full-scan access when neither indexed companion is present.
- Stop when only one of a modality's index and lookup inputs is present.
- Make one batched tabix query per indexed modality and filter exact phenotype IDs after the query.
- Keep timestamped logging in every WDL command.
- Do not build Docker images locally. Build and test images in GitHub Actions.
- Preserve existing WDL input JSON compatibility.

## File Structure

- `scripts/index_phenotype_bed.sh`: stream validation, coordinate sorting, BGZF compression, tabix indexing, lookup generation, and index QC.
- `workflows/index_phenotype_bed.wdl`: one-file reusable index workflow.
- `scripts/prepare_trans_window.R`: full-scan and indexed row access, shared selection logic, timing, and QC.
- `workflows/prepare_trans_window.wdl`: optional indexed inputs and CLI wiring.
- `envs/prepare-window-phenotypes.Dockerfile`: add bgzip/tabix and install the index script.
- `.github/workflows/prepare-window-phenotypes-image.yml`: run index and equivalence tests inside the image.
- `tests/fixtures/trans_window/generate_index_fixture.R`: small unsorted wide phenotype fixture with overlapping intervals.
- `tests/test_index_phenotype_bed.sh`: index artifact and query behavior tests.
- `tests/test_prepare_trans_window_indexed.R`: indexed/full-scan equivalence and mixed-mode tests.
- `tests/test_prepare_trans_window.R`: unit tests for lookup validation and indexed-input pair validation.
- `tests/test_prepare_trans_window_wdl.sh`: WDL interface and logging validation.
- `.dockstore.yml`, `README.md`, and `docs/trans-window-fine-mapping.md`: workflow registration and user instructions.

---

### Task 1: Reusable phenotype index artifacts

**Files:**
- Create: `scripts/index_phenotype_bed.sh`
- Create: `tests/fixtures/trans_window/generate_index_fixture.R`
- Create: `tests/test_index_phenotype_bed.sh`
- Modify: `envs/prepare-window-phenotypes.Dockerfile`
- Modify: `.github/workflows/prepare-window-phenotypes-image.yml`

**Interfaces:**
- Consumes: `index_phenotype_bed.sh INPUT MODALITY OUTPUT_DIR THREADS` where `INPUT` is gzip/BGZF or plain text.
- Produces: `OUTPUT_DIR/phenotypes.bed.gz`, `OUTPUT_DIR/phenotypes.bed.gz.tbi`, `OUTPUT_DIR/phenotype_lookup.tsv.gz`, and `OUTPUT_DIR/index_qc.tsv`.

- [ ] **Step 1: Write the failing definition test**

Create `tests/test_index_phenotype_bed.sh`. Before the implementation exists,
the test must fail at `test -x scripts/index_phenotype_bed.sh`. After that
check, it must generate the fixture, run the script, and assert:

```bash
test -s "$output_dir/phenotypes.bed.gz"
test -s "$output_dir/phenotypes.bed.gz.tbi"
test -s "$output_dir/phenotype_lookup.tsv.gz"
test -s "$output_dir/index_qc.tsv"
test "$(tabix "$output_dir/phenotypes.bed.gz" chr1:100-300 | wc -l)" -ge 2
gzip -cd "$output_dir/phenotype_lookup.tsv.gz" | head -n 1 | \
  grep -Fx $'phenotype_id\tchrom\tstart\tend'
```

The fixture generator must write an unsorted file with five phenotype rows,
three sample columns, two overlapping intervals, and metadata header names
that are not canonical.

- [ ] **Step 2: Run the definition test and verify the expected failure**

Run:

```bash
bash tests/test_index_phenotype_bed.sh
```

Expected: FAIL because `scripts/index_phenotype_bed.sh` does not exist.

- [ ] **Step 3: Implement the streaming index script**

The script must:

```bash
#!/usr/bin/env bash
set -euo pipefail

input="$1"
modality="$2"
output_dir="$3"
threads="$4"
```

It must select `gzip -cd -- "$input"` for `.gz` or `.bgz` input and `cat --
"$input"` for plain text. It must write one temporary normalized TSV. The
first four header values must become `#chrom`, `start`, `end`, and
`phenotype_id`; all sample headers must remain unchanged.

Use `awk -F '\t'` to validate a fixed field count, numeric nonnegative start,
end greater than start, nonempty unique phenotype IDs, and at least one data
row. Write the lookup as:

```text
phenotype_id  chrom  start  end
```

Sort data records with:

```bash
LC_ALL=C sort -t $'\t' -k1,1V -k2,2n -k3,3n
```

Prepend the normalized header, run `bgzip -@ "$threads"`, and create the
index with `tabix -p bed`. Query the first lookup interval and require at least
one returned record. Write QC columns for modality, phenotype rows, sample
columns, input bytes, BGZF bytes, lookup bytes, and elapsed seconds. Emit a
timestamped log before and after each major operation.

- [ ] **Step 4: Add tabix and the script to the phenotype image**

Add `tabix` to the existing `apt-get install` list in
`envs/prepare-window-phenotypes.Dockerfile`. Copy
`scripts/index_phenotype_bed.sh` to `/opt/mvsusie/scripts/` and make it
executable.

Add an in-image action step that runs:

```bash
docker run --rm \
  -v "$PWD:/workspace" \
  -w /workspace \
  prepare-window-phenotypes-ci:latest \
  bash tests/test_index_phenotype_bed.sh
```

- [ ] **Step 5: Run local structural checks**

Run:

```bash
bash -n scripts/index_phenotype_bed.sh
bash -n tests/test_index_phenotype_bed.sh
actionlint .github/workflows/prepare-window-phenotypes-image.yml
git diff --check
```

Expected: PASS. Do not build the image locally. The real bgzip/tabix behavior
will run in GitHub Actions.

- [ ] **Step 6: Commit Task 1**

```bash
git add scripts/index_phenotype_bed.sh \
  tests/fixtures/trans_window/generate_index_fixture.R \
  tests/test_index_phenotype_bed.sh \
  envs/prepare-window-phenotypes.Dockerfile \
  .github/workflows/prepare-window-phenotypes-image.yml
git commit -m "feat: build reusable phenotype indexes"
```

---

### Task 2: Indexed phenotype reader and output equivalence

**Files:**
- Modify: `scripts/prepare_trans_window.R`
- Modify: `tests/test_prepare_trans_window.R`
- Create: `tests/test_prepare_trans_window_indexed.R`
- Modify: `.github/workflows/prepare-window-phenotypes-image.yml`

**Interfaces:**
- Produces: `validate_indexed_input_pair(index_path, lookup_path, modality)`.
- Produces: `read_phenotype_lookup(path, modality)` returning a tibble with `phenotype_id`, `chrom`, `start`, and `end`.
- Produces: `read_prepare_phenotype_table_indexed(path, index_path, lookup_path, modality, requested_ids)` returning the same internal table shape as `read_prepare_phenotype_table()`.
- Produces: per-modality access QC fields `access_method`, `lookup_seconds`, `query_seconds`, `parse_seconds`, `query_rows`, and `modality_seconds`.

- [ ] **Step 1: Write failing lookup and pair-validation tests**

In `tests/test_prepare_trans_window.R`, add assertions that:

```r
expect_error_matching(
  validate_indexed_input_pair("file.tbi", NULL, "expression"),
  "both.*index.*lookup"
)
expect_error_matching(
  validate_indexed_input_pair(NULL, "lookup.tsv.gz", "expression"),
  "both.*index.*lookup"
)
stopifnot(identical(
  validate_indexed_input_pair(NULL, NULL, "expression"),
  "full_scan"
))
```

Write a lookup with duplicate IDs and require `read_phenotype_lookup()` to
reject it. Write a valid lookup and require canonical types and columns.

- [ ] **Step 2: Run the unit test and verify the expected failure**

Run:

```bash
Rscript tests/test_prepare_trans_window.R
```

Expected: FAIL because `validate_indexed_input_pair()` is not defined.

- [ ] **Step 3: Implement input-pair and lookup validation**

Implement:

```r
validate_indexed_input_pair <- function(index_path, lookup_path, modality) {
  has_index <- !is.null(index_path) && nzchar(index_path)
  has_lookup <- !is.null(lookup_path) && nzchar(lookup_path)
  if (xor(has_index, has_lookup)) {
    stop(modality, " requires both phenotype index and lookup inputs.", call. = FALSE)
  }
  if (!has_index) return("full_scan")
  if (!file.exists(index_path) || !file.exists(lookup_path)) {
    stop(modality, " indexed phenotype inputs do not exist.", call. = FALSE)
  }
  "tabix"
}
```

`read_phenotype_lookup()` must use `read_tsv()` with character defaults,
require the four canonical columns, convert coordinates to integer, reject
invalid coordinates, reject empty or duplicate IDs, and return only those
four columns.

- [ ] **Step 4: Write the failing indexed/full-scan equivalence test**

Create `tests/test_prepare_trans_window_indexed.R`. It must use the Task 1
fixture generator and index script to create indexed expression, splicing, and
protein files. Run `prepare_trans_window_data()` once with full-scan inputs and
once with all six indexed companions. Compare:

```r
stopifnot(identical(full_manifest, indexed_manifest))
stopifnot(identical(full_phenotypes, indexed_phenotypes))
stopifnot(identical(
  full_qc[preexisting_qc_columns],
  indexed_qc[preexisting_qc_columns]
))
stopifnot(all(indexed_qc$access_method == "tabix"))
```

Add a mixed-mode run with indexed expression, full-scan splicing, and indexed
protein. Require `c("tabix", "full_scan", "tabix")` in modality order. The
fixture must include overlapping coordinates; require no duplicate
`outcome_key` values.

- [ ] **Step 5: Run the equivalence test and verify the expected failure**

Run inside an environment with tabix:

```bash
Rscript tests/test_prepare_trans_window_indexed.R
```

Expected: FAIL because `prepare_trans_window_data()` does not accept indexed
companions.

- [ ] **Step 6: Implement one batched tabix query per modality**

Add optional indexed arguments to `prepare_trans_window_data()`. The indexed
reader must:

```r
lookup <- read_phenotype_lookup(lookup_path, modality)
matched <- lookup[match(requested_ids, lookup$phenotype_id), ]
```

Write `matched[c("chrom", "start", "end")]` to a temporary file with a `.bed`
suffix. Create a temporary symlink for the BGZF file and a second symlink at
`paste0(local_bgzf, ".tbi")`. Run one command:

```r
system2(
  "tabix",
  c("-R", shQuote(region_path), shQuote(local_bgzf)),
  stdout = query_path,
  stderr = error_path
)
```

Read the original header with `readLines(gzfile(path), n = 1L)`. Remove a
leading `#` from the first header field for the internal table. Parse only the
query result, derive `.phenotype_id` from column four, retain exact requested
IDs, reject duplicate IDs, and restore `requested_ids` order. Reuse
`select_joint_phenotype_rows()` for target/trans error semantics and outcome
key generation.

Return access timings with the selected table. Store the original phenotype
row count from the lookup on the indexed path and from the loaded table on the
full-scan path. Remove the current second full-file read during QC generation.

- [ ] **Step 7: Add extraction logging and QC fields**

For each modality, log the access method, requested IDs, lookup seconds, query
seconds, parse seconds, query rows, retained rows, and total seconds. Add the
six access fields from the Interfaces section to `window_qc.tsv`. Preserve all
existing QC values.

- [ ] **Step 8: Run R tests**

Run:

```bash
Rscript tools/lint_r.R
Rscript tests/test_prepare_trans_window.R
Rscript tests/test_prepare_trans_window_indexed.R
git diff --check
```

Expected: PASS in the tabix-enabled GitHub image. Locally, run the first two
commands and use the image workflow for the real tabix integration test.

- [ ] **Step 9: Commit Task 2**

```bash
git add scripts/prepare_trans_window.R \
  tests/test_prepare_trans_window.R \
  tests/test_prepare_trans_window_indexed.R \
  .github/workflows/prepare-window-phenotypes-image.yml
git commit -m "feat: query selected phenotypes with tabix"
```

---

### Task 3: WDL interfaces for indexing and optional fast access

**Files:**
- Create: `workflows/index_phenotype_bed.wdl`
- Modify: `workflows/prepare_trans_window.wdl`
- Modify: `tests/test_prepare_trans_window_wdl.sh`
- Modify: `tests/test_joint_workflow_requirements.sh`
- Create: `tests/test_index_phenotype_bed_wdl.sh`
- Modify: `.dockstore.yml`

**Interfaces:**
- `IndexPhenotypeBed` consumes `File phenotype_bed`, `String modality`, and `Int threads = 4`.
- `IndexPhenotypeBed` produces `indexed_phenotypes`, `phenotype_tbi`, `phenotype_lookup`, and `index_qc`.
- `PrepareTransWindow` adds optional index and lookup files for expression, splicing, and protein.

- [ ] **Step 1: Write failing WDL interface tests**

Make `tests/test_index_phenotype_bed_wdl.sh` require the new WDL file, workflow
name, four outputs, the index script command, timestamped log messages, and
runtime values. Extend `tests/test_prepare_trans_window_wdl.sh` to require:

```text
File? expression_phenotypes_tbi
File? expression_phenotype_lookup
File? splicing_phenotypes_tbi
File? splicing_phenotype_lookup
File? protein_phenotypes_tbi
File? protein_phenotype_lookup
```

Require all six matching CLI flags in the phenotype task. Extend the joint
workflow requirements test with the same declarations.

- [ ] **Step 2: Run WDL tests and verify the expected failure**

Run:

```bash
bash tests/test_index_phenotype_bed_wdl.sh
bash tests/test_prepare_trans_window_wdl.sh
```

Expected: FAIL because the index WDL and optional inputs do not exist.

- [ ] **Step 3: Create the one-file index workflow**

Create `workflows/index_phenotype_bed.wdl` with one task. Its command must run:

```bash
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting phenotype indexing for ~{modality}."
bash /opt/mvsusie/scripts/index_phenotype_bed.sh \
  '~{phenotype_bed}' '~{modality}' output '~{threads}'
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Phenotype indexing complete for ~{modality}."
```

Use the prepare-window phenotype image, 4 CPU, 16 GiB memory, and 500 GiB SSD.
Expose all four scalar outputs.

- [ ] **Step 4: Wire optional indexed inputs into PrepareTransWindow**

Add the six optional workflow and task inputs. Forward them to the R CLI with
WDL optional interpolation. Each optional file expression must use
`default=""` and must emit its flag only when the value is defined. Add a
timestamped message that states whether each modality requested indexed or
full-scan access.

Add these CLI options in `main()`:

```r
--expression-phenotypes-tbi
--expression-phenotype-lookup
--splicing-phenotypes-tbi
--splicing-phenotype-lookup
--protein-phenotypes-tbi
--protein-phenotype-lookup
```

Pass `NULL` for an omitted option and validate every pair before phenotype
loading begins.

- [ ] **Step 5: Register the index workflow**

Add this `.dockstore.yml` entry:

```yaml
  - name: index-phenotype-bed
    subclass: WDL
    primaryDescriptorPath: /workflows/index_phenotype_bed.wdl
```

- [ ] **Step 6: Run WDL and shell validation**

Run:

```bash
miniwdl check workflows/index_phenotype_bed.wdl
miniwdl check workflows/prepare_trans_window.wdl
bash tests/test_index_phenotype_bed_wdl.sh
bash tests/test_prepare_trans_window_wdl.sh
bash tests/test_joint_workflow_requirements.sh
git diff --check
```

Expected: PASS.

- [ ] **Step 7: Commit Task 3**

```bash
git add workflows/index_phenotype_bed.wdl \
  workflows/prepare_trans_window.wdl \
  scripts/prepare_trans_window.R \
  tests/test_index_phenotype_bed_wdl.sh \
  tests/test_prepare_trans_window_wdl.sh \
  tests/test_joint_workflow_requirements.sh \
  .dockstore.yml
git commit -m "feat: add indexed phenotype WDL inputs"
```

---

### Task 4: Documentation and end-to-end verification

**Files:**
- Modify: `README.md`
- Modify: `docs/trans-window-fine-mapping.md`
- Modify: `.github/workflows/prepare-window-phenotypes-image.yml`

**Interfaces:**
- Documents the index-once/run-many workflow and the optional fallback.
- Produces a published `latest` phenotype preparation image tested with real bgzip/tabix extraction.

- [ ] **Step 1: Extend the in-image smoke test**

Make the prepare-window image action run, in order:

```bash
bash tests/test_index_phenotype_bed.sh
Rscript tests/test_prepare_trans_window.R
Rscript tests/test_prepare_trans_window_indexed.R
```

The indexed equivalence test must print per-modality timing logs. The action
must fail before image publication if any test fails.

- [ ] **Step 2: Document index creation and fast-path inputs**

Update `README.md` with the index workflow link. Update
`docs/trans-window-fine-mapping.md` with:

- one `IndexPhenotypeBed` launch per source phenotype file;
- the four index outputs;
- the six optional `PrepareTransWindow` inputs;
- the rule that both companions for a modality must be present;
- mixed indexed/full-scan mode;
- the absence of a GTF dependency;
- timing fields in `window_qc.tsv` and task logs.

- [ ] **Step 3: Run the complete local verification set**

Run:

```bash
Rscript tools/lint_r.R
Rscript tests/test_prepare_trans_window.R
miniwdl check workflows/index_phenotype_bed.wdl
miniwdl check workflows/prepare_trans_window.wdl
bash tests/test_index_phenotype_bed_wdl.sh
bash tests/test_prepare_trans_window_wdl.sh
bash tests/test_joint_workflow_requirements.sh
bash tests/test_trans_window_mvsusie_container.sh
actionlint .github/workflows/prepare-window-phenotypes-image.yml
git diff --check
```

Expected: PASS. Do not run a local Docker build.

- [ ] **Step 4: Commit documentation and CI changes**

```bash
git add README.md docs/trans-window-fine-mapping.md \
  .github/workflows/prepare-window-phenotypes-image.yml
git commit -m "docs: describe indexed phenotype preparation"
```

- [ ] **Step 5: Push and verify GitHub Actions**

Push the branch to `main`. Monitor the prepare-window phenotype image action.
Require all of these steps to pass before reporting completion:

- image build;
- package and binary checks;
- phenotype index test;
- full-scan preparation test;
- indexed/full-scan equivalence test;
- image publication.

Record the successful run URL and the final commit in the handoff.
