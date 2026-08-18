# mvSusieR workflows

<!-- workflow-badges:start -->
[![LD-pruning container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/docker-image.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->


## LD-pruning workflow

`workflows/ld_pruning.wdl` creates an LD-pruned sentinel-variant set with
PLINK 2 from explicit PGEN, PVAR, and PSAM inputs.

The container is published as `ghcr.io/evin-padhi/mvsusier-ld-pruning` after a
successful push to `main`.

See [the LD-pruning guide](docs/ld-pruning.md) for inputs, a MiniWDL example,
and outputs.
