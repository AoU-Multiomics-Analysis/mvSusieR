# mvSusieR workflows

<!-- workflow-badges:start -->
[![LD-pruning container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/docker-image.yml)
[![Prepare-window genotypes container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/prepare-window-genotypes-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/prepare-window-genotypes-image.yml)
[![Prepare-window phenotypes container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/prepare-window-phenotypes-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/prepare-window-phenotypes-image.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/r-lint.yml)
[![Trans-window mvSuSiE container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/trans-window-mvsusie-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/trans-window-mvsusie-image.yml)
[![TransQTL BAM filtering container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/transqtl-bam-filtering-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/transqtl-bam-filtering-image.yml)
[![TransQTL mappability container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/transqtl-mappability-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/transqtl-mappability-image.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->

This repository contains containerized WDL workflows for preparing molecular
and genotype data, filtering RNA-seq alignments, and performing multivariate
fine-mapping of trans-QTL signals.

## Pipeline overview

### TransQTL mappability generation

[`TransQTLMappability`](workflows/transqtl_mappability.wdl) computes a
read-length-specific GenMap mappability track from a reference FASTA and
produces a low-mappability BED. Its purpose is to provide a reproducible,
reference-compatible mask that can be used to remove reads aligned to regions
where mapping is less reliable. See the [mappability guide](docs/transqtl-mappability.md).

### TransQTL BAM filtering and RNA-SeQC2

[`TransQTLBamFiltering`](workflows/transqtl_bam_filtering.wdl) removes whole
read templates with non-unique alignments or overlap with the low-mappability
BED, then requantifies the retained BAM with RNA-SeQC2. Its purpose is to
reduce mapping-related artifacts before generating gene-level expression
measurements while keeping filtering and quantification in one task. See the
[BAM filtering guide](docs/transqtl-bam-filtering.md).

### Trans-window data preparation

[`PrepareTransWindow`](workflows/prepare_trans_window.wdl) extracts the dosage
records and phenotype features needed for one trans window. Its purpose is to
turn genome-wide dosage, association, and molecular-phenotype inputs into
small, window-specific files that can be processed efficiently downstream.

### Trans-window multivariate fine-mapping

[`TransWindowMvSusie`](workflows/trans_window_mvsusie.wdl) fits mvSuSiE across
trans windows using multiple molecular phenotypes jointly, then summarizes and
merges the results. Its purpose is to identify shared genetic signals and
quantify variant contributions with PIPs, credible sets, component effects, and
window-level QC. See the [trans-window fine-mapping guide](docs/trans-window-fine-mapping.md).

### LD pruning

[`LDPruning`](workflows/ld_pruning.wdl) applies PLINK 2 quality and
linkage-disequilibrium filters and creates a pruned variant set. Its purpose is
to reduce redundant correlated variants before analyses that require a more
independent set of markers. See the [LD-pruning guide](docs/ld-pruning.md).

## Using the workflows

All workflows are WDL 1.0 and can be run with MiniWDL or a Cromwell-compatible
engine. Workflow registrations and descriptor paths are listed in
[`.dockstore.yml`](.dockstore.yml), and container images are published to
GitHub Container Registry by the corresponding GitHub Actions workflows.

The [workflows directory](workflows/) contains the executable WDL files. The
[docs directory](docs/) contains detailed input, output, and implementation
notes for individual pipelines.
