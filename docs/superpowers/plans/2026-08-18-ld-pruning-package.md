# LD-pruning package Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reproducible PLINK 2 LD-pruning WDL, its GHCR container, CI publishing, and Dockstore workflow metadata.

**Architecture:** A WDL task receives explicit PGEN, PVAR, and PSAM files, stages them to a fixed PLINK prefix, produces `--indep-pairwise` IDs, and extracts them as a compact PGEN set. A pinned PLINK 2 container executes the task. Existing repository workflows validate the WDL and build/publish the specialized image.

**Tech Stack:** WDL 1.1, PLINK 2, Docker, GitHub Actions, GitHub Container Registry, Dockstore.

**Spec:** `docs/superpowers/specs/2026-08-18-ld-pruning-package-design.md`

## Global Constraints

- Inputs are explicit PGEN, PVAR, and PSAM files; callers never supply a shared prefix.
- Defaults are MAF 0.05, genotype missingness 0.01, and 1000 kb / 50 variants / r2 0.1 pruning.
- The OCI source label is `https://github.com/evin-padhi/mvSusieR`.
- Images publish only on `main` or manual dispatch, never on pull requests.
- The GHCR image is `ghcr.io/evin-padhi/mvsusier-ld-pruning`.

---

### Task 1: Implement and validate the WDL descriptor

**Files:**
- Create: `workflows/ld_pruning.wdl`
- Create: `tests/test_ld_pruning_wdl.sh`

**Interfaces:**
- Consumes: `File input_pgen`, `File input_pvar`, `File input_psam`, optional keep/exclude files, and pruning parameters.
- Produces: `ld_pruning.prune.in`, `ld_pruning.prune.out`, `ld_pruned.pgen`, `ld_pruned.pvar`, `ld_pruned.psam`, and `ld_pruned.snplist`.

- [ ] **Step 1: Write the failing test**

```bash
#!/usr/bin/env bash
set -euo pipefail
miniwdl check workflows/ld_pruning.wdl
```

The test catches a missing or malformed WDL descriptor; it runs the same MiniWDL command as CI.

- [ ] **Step 2: Verify the test fails**

Run: `bash tests/test_ld_pruning_wdl.sh`

Expected: failure because `workflows/ld_pruning.wdl` does not exist.

- [ ] **Step 3: Write the minimal WDL**

Create `workflow LDPruning` and task `LDPruningTask`. Stage the three inputs as `input.pgen`, `input.pvar`, and `input.psam`; run PLINK filters and `--indep-pairwise`; fail on an empty prune set; use `--extract`, `--make-pgen`, and `--write-snplist`; and set the GHCR image in the runtime.

- [ ] **Step 4: Verify the test passes**

Run: `bash tests/test_ld_pruning_wdl.sh`

Expected: MiniWDL exits 0.

### Task 2: Implement and smoke-test the container

**Files:**
- Create: `envs/ld-pruning.Dockerfile`
- Create: `tests/test_ld_pruning_container.sh`

**Interfaces:**
- Consumes: Docker and `envs/ld-pruning.Dockerfile`.
- Produces: an image tagged `ld-pruning-test` with `plink2` as its entrypoint.

- [ ] **Step 1: Write the failing test**

```bash
#!/usr/bin/env bash
set -euo pipefail
docker build --tag ld-pruning-test --file envs/ld-pruning.Dockerfile .
docker run --rm ld-pruning-test plink2 --version
```

The test catches an image that cannot expose a working PLINK 2 executable.

- [ ] **Step 2: Verify the test fails**

Run: `bash tests/test_ld_pruning_container.sh`

Expected: failure because `envs/ld-pruning.Dockerfile` does not exist.

- [ ] **Step 3: Write the minimal Dockerfile**

Use a pinned Ubuntu base and pinned PLINK 2 release. Install `curl`, `ca-certificates`, and `unzip`; install PLINK at `/usr/local/bin/plink2`; set source, title, description, and license OCI labels; and set `CMD ["plink2"]` so a WDL backend can still execute its shell command.

- [ ] **Step 4: Verify the test passes**

Run: `bash tests/test_ld_pruning_container.sh`

Expected: Docker builds successfully and PLINK reports its version.

### Task 3: Configure CI and Dockstore

**Files:**
- Create: `tests/test_ld_pruning_metadata.sh`
- Modify: `.github/workflows/docker-image.yml`
- Modify: `.dockstore.yml`

**Interfaces:**
- Consumes: `envs/ld-pruning.Dockerfile`, `workflows/ld_pruning.wdl`.
- Produces: PR image validation, main-branch/dispatch GHCR publication, and a Dockstore workflow declaration.

- [ ] **Step 1: Write the failing metadata test**

```bash
#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PY'
import yaml
with open('.dockstore.yml') as handle:
    config = yaml.safe_load(handle)
assert config['workflows'][0]['name'] == 'ld-pruning'
assert config['workflows'][0]['primaryDescriptorPath'] == '/workflows/ld_pruning.wdl'
PY
```

The test catches Dockstore metadata that still points to the template placeholder.

- [ ] **Step 2: Verify the test fails**

Run: `bash tests/test_ld_pruning_metadata.sh`

Expected: failure because the existing Dockstore descriptor uses a placeholder path.

- [ ] **Step 3: Write minimal workflow and metadata changes**

Update `docker-image.yml` to build `envs/ld-pruning.Dockerfile`, tag `ghcr.io/evin-padhi/mvsusier-ld-pruning`, run `plink2 --version`, and condition `push` on `main` or manual dispatch. Replace `.dockstore.yml` with an `ld-pruning` WDL declaration whose descriptor path is `/workflows/ld_pruning.wdl`.

- [ ] **Step 4: Verify the tests pass**

Run: `bash tests/test_ld_pruning_metadata.sh && bash tests/test_ld_pruning_wdl.sh`

Expected: both commands exit 0.

### Task 4: Document and verify the complete package

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes: the WDL input interface and package image name.
- Produces: a runnable WDL input example and output description.

- [ ] **Step 1: Document usage**

Add an input JSON example that supplies PGEN, PVAR, and PSAM independently, lists default pruning parameters, and describes the compact PGEN and prune-list outputs.

- [ ] **Step 2: Run the full verification set**

Run: `bash tests/test_ld_pruning_wdl.sh && bash tests/test_ld_pruning_metadata.sh && bash tests/test_ld_pruning_container.sh`

Expected: every test exits 0.
