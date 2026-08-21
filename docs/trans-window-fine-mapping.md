# Trans-window preparation and multivariate fine-mapping

The trans-window workflows split a large trans-QTL analysis into manageable
units. `PrepareTransWindow` creates the inputs for one window, and
`TransWindowMvSusie` fits and summarizes the multivariate fine-mapping model
across all listed windows.

## Data preparation

[`workflows/prepare_trans_window.wdl`](../workflows/prepare_trans_window.wdl)
has two independent responsibilities:

- `PrepareWindowGenotypes` uses the window association manifest and a
  tabix-indexed genome-wide dosage file to extract the dosage interval and
  write a one-row window manifest.
- `PrepareWindowPhenotypes` selects cis phenotypes overlapping the window and
  the top trans phenotypes by association strength from expression, splicing,
  or isoform-usage inputs. It writes a combined phenotype file, a phenotype
  manifest, and QC output.

The workflow is invoked once per window. Its main outputs are
`window_dosage.tsv`, `window_manifest.tsv`, `window_phenotypes.bed.gz`,
`window_phenotypes.tsv`, and `window_qc.tsv`.

## Multivariate fine-mapping

[`workflows/trans_window_mvsusie.wdl`](../workflows/trans_window_mvsusie.wdl)
scatters one `RunMvSusie` task per window. Each task aligns the dosage and
phenotype samples, applies the configured covariates, fits mvSuSiE jointly
across the molecular phenotypes, and writes prepared data and a fit object.
Downstream tasks create per-window summaries and merge them across windows.

The workflow produces variant PIP tables, credible sets, component-effect
tables, and window-level QC, along with merged versions of those outputs. The
current WDL uses the canonical mvSuSiE prior.

## Inputs and execution

The preparation workflow consumes a window ID, a tabix-indexed dosage file, a
window association table, and phenotype files with matching modality labels.
The fine-mapping workflow consumes the resulting window table, phenotype data,
covariate files, and one dosage file per window.

Both workflows use container images published to GHCR. See the WDL files for
the complete input schema and defaults, and use MiniWDL or a Cromwell-compatible
engine to run them.
