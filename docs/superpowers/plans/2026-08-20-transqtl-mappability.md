# TransQTL mappability generation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reproducible GenMap-based WDL workflow, micromamba image, CI smoke test, and Dockstore registration for producing read-length-specific low-mappability BED files.

**Architecture:** `TransQTLMappability` will run one containerized task that indexes a supplied reference FASTA, computes GenMap `(k,e)` mappability, thresholds and merges low-scoring intervals, and emits a BED plus provenance metadata. GitHub Actions will validate the WDL, build the image, run a tiny-reference smoke test, and publish the image to GHCR on non-PR events.

**Tech Stack:** WDL 1.0, miniwdl, Cromwell-compatible shell commands, Docker, micromamba, GenMap, bedtools, GitHub Actions, Ruby YAML.

**Spec:** `docs/superpowers/specs/2026-08-20-transqtl-mappability-design.md`

## Global Constraints

- Use GenMap for `(k,e)` mappability with defaults `k=146`, `e=2`, and threshold `<1.0`.
- Require the caller to provide the reference FASTA; do not download or silently select a GRCh38 contig set.
- Use `mambaorg/micromamba:2.3.3`, strict `conda-forge` then `bioconda`, and expose `/opt/conda/bin` directly on `PATH`.
- Preserve the repository’s existing unrelated working-tree changes; stage only files created for this feature.
- Keep generated index files internal to the task; publish the filtered BED, full BEDGraph, and metadata as workflow outputs.

---

### Task 1: Add failing validation tests

**Files:**
- Create: `tests/test_transqtl_mappability_wdl.sh`
- Create: `tests/test_transqtl_mappability_metadata.sh`
- Create: `tests/test_transqtl_mappability_micromamba.sh`
- Create: `tests/test_transqtl_mappability_container.sh`

**Interfaces:**
- Tests consume the future files `workflows/transqtl_mappability.wdl`, `envs/transqtl-mappability.Dockerfile`, `envs/transqtl-mappability.environment.yml`, and `.dockstore.yml`.
- Tests produce clear failures until those files and their expected interfaces exist.

- [ ] **Step 1: Write the WDL static test**

Create a Bash test that runs `miniwdl check workflows/transqtl_mappability.wdl` and checks for `TransQTLMappability`, `GenerateTransQTLMappability`, `reference_fasta`, `kmer_length`, `max_mismatches`, `mappability_threshold`, `preemptible_tries`, `genmap index`, `genmap map`, `bedtools merge`, `low_mappability_bed`, `mappability_bedgraph`, `mappability_metadata`, `log()`, and the output-validation log message.

- [ ] **Step 2: Write the Dockstore metadata test**

Create a Ruby YAML check that finds workflow name `transqtl-mappability` and requires primary descriptor path `/workflows/transqtl_mappability.wdl`.

- [ ] **Step 3: Write the micromamba environment test**

Create a Bash test that requires `conda-forge`, `bioconda`, `genmap`, and `bedtools` in the environment file; requires the micromamba base image, activation flag, strict channel priority, and direct Conda PATH in the Dockerfile; and rejects `apt-get`.

- [ ] **Step 4: Write the container smoke-test contract**

Create a Bash test that builds `envs/transqtl-mappability.Dockerfile`, verifies `genmap --version`, `bedtools --version`, and `command -v` under `/bin/bash`, then runs GenMap on a small synthetic FASTA and asserts that the generated BEDGraph is nonempty and BEDTools can sort and merge an interval from it.

- [ ] **Step 5: Run the tests to verify the red state**

Run `bash tests/test_transqtl_mappability_wdl.sh`, `bash tests/test_transqtl_mappability_metadata.sh`, and `bash tests/test_transqtl_mappability_micromamba.sh`. Expected result: the tests fail because the new workflow, image, environment, and Dockstore entry do not yet exist.

### Task 2: Implement the WDL workflow

**Files:**
- Create: `workflows/transqtl_mappability.wdl`
- Test: `tests/test_transqtl_mappability_wdl.sh`

**Interfaces:**
- Consumes `reference_fasta`, `kmer_length`, `max_mismatches`, `mappability_threshold`, `output_prefix`, `threads`, and `preemptible_tries`.
- Produces `low_mappability_bed`, `mappability_bedgraph`, and `mappability_metadata`.

- [ ] **Step 1: Add workflow and task inputs**

Define WDL 1.0 workflow `TransQTLMappability` calling task `GenerateTransQTLMappability`. Use defaults `kmer_length = 146`, `max_mismatches = 2`, `mappability_threshold = 1.0`, `output_prefix = "GRCh38.K146.m2"`, `threads = 16`, and `preemptible_tries = 3`.

- [ ] **Step 2: Add the GenMap command block**

In the task command, use `set -euo pipefail`, timestamped logging, and run `genmap index -F reference.fasta -I genmap_index` followed by `genmap map -K ~{kmer_length} -E ~{max_mismatches} -I genmap_index -O genmap_output -t -w -bg`. Locate the generated BEDGraph from `genmap_output`, fail if none exists, and copy it to `~{output_prefix}.mappability.bedGraph`.

- [ ] **Step 3: Add thresholding and interval normalization**

Use `awk` to retain BEDGraph records whose fourth column is below `~{mappability_threshold}`, then run `bedtools sort` and `bedtools merge` to emit the three-column file `~{output_prefix}.low_mappability.bed`. Require the full BEDGraph to be nonempty, allow the low-mappability BED to be empty, and count the resulting intervals in the log.

- [ ] **Step 4: Add provenance metadata and outputs**

Write a tab-delimited metadata file containing `kmer_length`, `max_mismatches`, `mappability_threshold`, `reference_fasta`, `reference_size_bytes`, and `reference_sha256`. Expose the BED, BEDGraph, and metadata files in the task and workflow output blocks.

- [ ] **Step 5: Add validation and runtime settings**

Reject invalid values where `kmer_length < 1`, `max_mismatches < 0`, `mappability_threshold <= 0`, `mappability_threshold > 1`, or `threads < 1`. Set runtime Docker to `ghcr.io/aou-multiomics-analysis/transqtl-mappability:latest`, CPU to `threads`, memory to `64 GiB`, disks to `local-disk 200 SSD`, and preemptible to `preemptible_tries`.

- [ ] **Step 6: Run the WDL test to verify green**

Run `bash tests/test_transqtl_mappability_wdl.sh`. Expected result: miniwdl parses the workflow and all required interface tokens are present.

### Task 3: Implement the micromamba environment and container

**Files:**
- Create: `envs/transqtl-mappability.environment.yml`
- Create: `envs/transqtl-mappability.Dockerfile`
- Test: `tests/test_transqtl_mappability_micromamba.sh`
- Test: `tests/test_transqtl_mappability_container.sh`

**Interfaces:**
- Produces an image containing `genmap`, `bedtools`, and `micromamba` with all tools available without relying on an entrypoint activation shell.

- [ ] **Step 1: Define the Conda environment**

Create a base environment with channels in this order: `conda-forge`, then `bioconda`; dependencies are `genmap` and `bedtools`.

- [ ] **Step 2: Define the Dockerfile**

Use `FROM mambaorg/micromamba:2.3.3`, copy the environment file, set `MAMBA_DOCKERFILE_ACTIVATE=1`, set the direct Conda PATH, configure strict channel priority, install the base environment, clean caches, and run both tool version checks during the build.

- [ ] **Step 3: Run static environment validation**

Run `bash tests/test_transqtl_mappability_micromamba.sh`. Expected result: all package, channel, and PATH checks pass.

- [ ] **Step 4: Build and run the container smoke test**

Run `bash tests/test_transqtl_mappability_container.sh`. Expected result: the image builds, both tools are discoverable from the default and `/bin/bash` entrypoints, and a tiny reference produces a nonempty low-mappability BED.

### Task 4: Add GitHub Actions and Dockstore registration

**Files:**
- Create: `.github/workflows/transqtl-mappability-image.yml`
- Modify: `.dockstore.yml`
- Test: `tests/test_transqtl_mappability_metadata.sh`

**Interfaces:**
- CI builds the image from `envs/transqtl-mappability.Dockerfile`, runs WDL/static/container smoke checks, and publishes `ghcr.io/${{ github.repository_owner }}/transqtl-mappability` on non-PR events.
- Dockstore exposes `/workflows/transqtl_mappability.wdl` as `transqtl-mappability`.

- [ ] **Step 1: Register the Dockstore workflow**

Append a workflow entry with `name: transqtl-mappability`, `subclass: WDL`, and `primaryDescriptorPath: /workflows/transqtl_mappability.wdl`.

- [ ] **Step 2: Add CI triggers and build configuration**

Create a workflow triggered by manual dispatch, pushes to `main`, and pull requests targeting `main`, limited to the new WDL, Dockerfile, environment, tests, workflow file, and Dockstore metadata. Use Docker Buildx and build the image with `load: true` for PR smoke tests.

- [ ] **Step 3: Add validation and smoke-test steps**

Run the static tests, install miniwdl with pipx, run `miniwdl check`, build the image, verify GenMap and bedtools versions, and run the synthetic reference smoke test inside the image.

- [ ] **Step 4: Add GHCR publication**

Grant `contents: read` and `packages: write`, log into GHCR only for non-PR events, generate branch/SHA/latest tags with `docker/metadata-action`, and publish with `docker/build-push-action`.

- [ ] **Step 5: Validate Dockstore metadata**

Run `bash tests/test_transqtl_mappability_metadata.sh`. Expected result: the registered name and descriptor path match the new workflow.

### Task 5: Full verification and handoff

**Files:**
- Test: all new `tests/test_transqtl_mappability_*.sh` files

- [ ] **Step 1: Run all new tests**

Run the four new tests: WDL, metadata, micromamba, and container smoke test. Each must exit with status 0.

- [ ] **Step 2: Re-run existing TransQTL tests**

Run `bash tests/test_transqtl_bam_filtering_wdl.sh`, `bash tests/test_transqtl_bam_filtering_metadata.sh`, and `bash tests/test_transqtl_bam_filtering_micromamba.sh`.

- [ ] **Step 3: Inspect the final diff and status**

Run `git diff --check`, `git status --short`, and `git diff --stat`. Confirm that only the specification, plan, new mappability files, and Dockstore entry are present in the feature branch; do not stage unrelated files from the original checkout.

- [ ] **Step 4: Report verification evidence**

Report the exact test commands and their observed exit status, identify the generated workflow and image paths, and state that no Docker image build or GitHub push was performed unless separately requested.
