# mvSusieR workflows

<!-- workflow-badges:start -->
[![LD-pruning container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/docker-image.yml)
[![Prepare-window genotypes container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/prepare-window-genotypes-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/prepare-window-genotypes-image.yml)
[![Prepare-window phenotypes container](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/prepare-window-phenotypes-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/mvSusieR/actions/workflows/prepare-window-phenotypes-image.yml)
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
  "TransWindowMvSusie.phenotype_data": "window_phenotypes.bed.gz",
  "TransWindowMvSusie.covariate_files": ["covariates.tsv"],
  "TransWindowMvSusie.covariate_modalities": ["shared"],
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

`window_phenotypes.tsv` has columns `window_id`, `phenotype_id`, `modality`, and `phenotype_file`. The `phenotype_data` input should normally be the single combined `window_phenotypes.bed.gz` produced by `PrepareWindowPhenotypes`. Supported modalities are `expression`, `splicing`, and `isoform_usage`. Expression and splicing use the established feature-by-sample layout with the phenotype ID in column 4; isoform-usage files use the phenotype ID in column 1.

Covariate files are covariate-by-sample tables: the first column is the covariate ID and the remaining columns are sample IDs. `covariate_modalities` labels each covariate file as `shared`, `expression`, `splicing`, or `isoform_usage`; the arrays must have matching lengths. Shared covariates are applied to every modality, while modality-specific files are applied only to that modality. `keep_samples` is optional and contains one sample ID per line.

The R entrypoints use the `optparse` package and expose the same long-form flags shown in the WDL commands, including `--window-phenotypes`, `--max-iter`, and `--output-dir`.

Run locally with:

```bash
miniwdl run workflows/trans_window_mvsusie.wdl -i trans_window.inputs.json
```

The model image is built and published to GHCR by the `trans-window-mvsusie-image` GitHub Actions workflow. Each scattered `RunMvSusie` task uses `scripts/run_window_mvsusie.R` to prepare one window and fit mvSuSiE in the same job, writing both the prepared-window RDS and the fit RDS. The summary and merge tasks consume those outputs downstream.

### Model defaults and behavior

The canonical mixture prior is fixed in v1. Residual covariance uses the `mvsusieR` default behavior (`residual_variance = NULL`); the fit estimates residual variance and does not estimate prior variance. Other defaults are `L = 10`, `max_iter = 100`, `tol = 1e-4`, `coverage = 0.95`, `min_abs_corr = 0.5`, genotype and phenotype variance thresholds of `1e-8`, one thread, `16 GiB` memory, and a `500 GB` local disk. Samples are intersected across inputs, phenotypes are rank-inverse-normal transformed, and covariates are residualized before fitting. The genotype matrix is residualized against the union of all nuisance covariates, while each phenotype is residualized only against the covariates assigned to its modality; expression is not regressed out of splicing. A failed or non-converged window fails the overall workflow.

### Outputs

The workflow emits per-window prepared data, fit RDS files, variant PIP tables, credible-set tables, component-effect tables, and QC tables. It also emits merged `variant_pip.tsv.gz`, `credible_sets.tsv.gz`, `component_effects.tsv.gz`, and `window_qc.tsv` files across all windows.

### Per-window input preparation

`workflows/prepare_trans_window.wdl` prepares one window per workflow invocation using two independent tasks. The caller supplies one `window_id`, the complete `windows.tsv`, a tabix-indexed genome-wide dosage file plus its `.tbi` index, that window's TensorQTL association file, and matching arrays of phenotype files and modality labels. The workflow does not scatter over the windows table; launch one invocation per window from the outer scheduler. `PrepareWindowGenotypes` and `PrepareWindowPhenotypes` can therefore be scheduled independently when desired.

`PrepareWindowGenotypes` extracts the requested dosage interval with `tabix` and writes the one-row dosage manifest. `PrepareWindowPhenotypes` runs `scripts/prepare_trans_window.R`; for now, each gzipped phenotype BED is read fully into R and filtered in memory. Expression and splicing files use the first three columns for chromosome/start/end and column 4 for the phenotype ID. The default `top_n_trans_phenotypes` is 25 per modality: trans phenotypes are ranked within each modality by their minimum TensorQTL p-value across variants, and the top N from each modality are retained. All cis phenotypes overlapping the selected window are retained independently; set `extract_cis_window_phenotypes` to `false` for trans-only preparation.

The outputs are `window_dosage.tsv`, a one-row `window_manifest.tsv` compatible with the dosage-manifest reader, one combined `window_phenotypes.bed.gz` file, `window_phenotypes.tsv` as the phenotype/modality key consumed by `run_window_mvsusie.R`, and `window_qc.tsv`. The combined phenotype file is intentionally not tabix-indexed in this first implementation; indexing or a one-time in-memory/columnar cache can be added later if full-file scans become the bottleneck.

## LD-pruning workflow

`workflows/ld_pruning.wdl` creates an LD-pruned sentinel-variant set with
PLINK 2 from explicit PGEN, PVAR, and PSAM inputs.

The container is published as `ghcr.io/aou-multiomics-analysis/mvsusier-ld-pruning` after a
successful push to `main`.

See [the LD-pruning guide](docs/ld-pruning.md) for inputs, a MiniWDL example,
and outputs.
