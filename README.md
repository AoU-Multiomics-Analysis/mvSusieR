# mvSusieR workflows

<!-- workflow-badges:start -->
[![LD-pruning container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/docker-image.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->


## Trans-window multivariate fine-mapping

`workflows/trans_window_mvsusie.wdl` runs the canonical-prior `mvsusieR` model independently across trans windows. The v1 workflow learns no mashr covariance and accepts raw window-level matrices plus manifests. Genome-wide association summaries can be added later for mashr covariance learning.

### Inputs

The WDL input JSON supplies these fields:

```json
{
  "TransWindowMvSusie.windows_tsv": "windows.tsv",
  "TransWindowMvSusie.window_phenotypes_tsv": "window_phenotypes.tsv",
  "TransWindowMvSusie.phenotype_files": ["expression.tsv", "splicing.tsv", "isoform_usage.tsv"],
  "TransWindowMvSusie.covariate_files": ["covariates.tsv"],
  "TransWindowMvSusie.keep_samples": null,
  "TransWindowMvSusie.L": 10,
  "TransWindowMvSusie.max_iter": 100,
  "TransWindowMvSusie.tol": 0.0001,
  "TransWindowMvSusie.coverage": 0.95,
  "TransWindowMvSusie.min_abs_corr": 0.5,
  "TransWindowMvSusie.min_genotype_variance": 1e-8,
  "TransWindowMvSusie.min_phenotype_variance": 1e-8,
  "TransWindowMvSusie.n_thread": 1
}
```

`windows.tsv` is tab-delimited with one header row and columns `window_id`, `chrom`, `start`, `end`, and `dosage_file`. Coordinates are 0-based, half-open. Each dosage file is window-specific and wide: `CHROM`, `POS`, `REF`, `ALT`, followed by one column per sample.

`window_phenotypes.tsv` has columns `window_id`, `phenotype_id`, `modality`, and `phenotype_file`. Supported modalities are `expression`, `splicing`, and `isoform_usage`. Expression and splicing files use the established feature-by-sample layout with the phenotype ID in column 4; isoform-usage files use the phenotype ID in column 1. `phenotype_files` must contain every file referenced by the manifest.

Covariate files are covariate-by-sample tables: the first column is the covariate ID and the remaining columns are sample IDs. `keep_samples` is optional and contains one sample ID per line.

The R entrypoints use the `optparse` package and expose the same long-form flags shown in the WDL commands, including `--window-phenotypes`, `--max-iter`, and `--output-dir`.

Run locally with:

```bash
miniwdl run workflows/trans_window_mvsusie.wdl -i trans_window.inputs.json
```

The WDL declares a future container runtime image, but this repository slice does not build or publish that image yet.

### Model defaults and behavior

The canonical mixture prior is fixed in v1. Residual covariance uses the `mvsusieR` default behavior (`residual_variance = NULL`); the fit estimates residual variance and does not estimate prior variance. Other defaults are `L = 10`, `max_iter = 100`, `tol = 1e-4`, `coverage = 0.95`, `min_abs_corr = 0.5`, genotype and phenotype variance thresholds of `1e-8`, and one thread. Samples are intersected across inputs, phenotypes are rank-inverse-normal transformed, and covariates are residualized before fitting. A failed or non-converged window fails the overall workflow.

### Outputs

The workflow emits per-window prepared data, fit RDS files, variant PIP tables, credible-set tables, component-effect tables, and QC tables. It also emits merged `variant_pip.tsv.gz`, `credible_sets.tsv.gz`, `component_effects.tsv.gz`, and `window_qc.tsv` files across all windows.

## LD-pruning workflow

`workflows/ld_pruning.wdl` creates an LD-pruned sentinel-variant set with
PLINK 2 from explicit PGEN, PVAR, and PSAM inputs.

The container is published as `ghcr.io/evin-padhi/mvsusier-ld-pruning` after a
successful push to `main`.

See [the LD-pruning guide](docs/ld-pruning.md) for inputs, a MiniWDL example,
and outputs.
