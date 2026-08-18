# LD-pruning package design

## Goal

Add a containerized, Dockstore-discoverable PLINK 2 LD-pruning workflow to `AoU-Multiomics-Analysis/mvSusieR` for making reproducible genome-wide sentinel variant sets.

## Workflow interface

`workflows/ld_pruning.wdl` exposes three required, independent genotype files: `input_pgen`, `input_pvar`, and `input_psam`. It accepts optional sample and variant lists (`keep_samples`, `exclude_variants`) plus QTL-oriented pruning parameters: `maf` (0.05), `geno` (0.01), a 1000-kb window, a 50-variant step, and `r2` (0.1).

The task stages the three inputs under a fixed local PLINK prefix, runs `plink2 --indep-pairwise`, verifies that `prune.in` is nonempty, and extracts the retained IDs to a compact PGEN/PVAR/PSAM dataset. It returns both PLINK's prune lists and the compact dataset.

## Container

`envs/ld-pruning.Dockerfile` supplies a pinned PLINK 2 release on Ubuntu and sets standard OCI labels. The source label points at this repository so GitHub links the published image to `mvSusieR`. The image name is `ghcr.io/aou-multiomics-analysis/mvsusier-ld-pruning`.

## Automation

The existing Docker workflow becomes the LD-pruning image workflow. It builds on pull requests for validation and publishes only on pushes to `main` or an explicitly selected `workflow_dispatch` run. The existing WDL-validation workflow continues to run `miniwdl check` on the new WDL.

## Dockstore

`.dockstore.yml` declares `ld-pruning` as a WDL workflow, with `/workflows/ld_pruning.wdl` as its primary descriptor. The task runtime uses the GHCR image above.

## Verification

Run `miniwdl check workflows/ld_pruning.wdl` locally when available. Validate the container file by building it and running `plink2 --version`. The image workflow must not publish images on pull requests.
