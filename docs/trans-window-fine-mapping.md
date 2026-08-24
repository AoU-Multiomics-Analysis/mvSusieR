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
default uses the canonical mvSuSiE prior. Set prior_method to mashr to learn a
window-specific prior. This option calculates marginal effects for all retained
SNPs in the window. It obtains PCA covariance inputs from SNPs that pass the
configured lfsr threshold. It can refine those inputs with extreme
deconvolution. The mashr mixture fit uses all retained SNPs. The pipeline
preserves the effect-covariance scale from mashr when it supplies the prior to
mvSuSiE.

By default, mvSuSiE fits the fixed number of components in `L`. To use greedy
component selection, set `L_greedy` to a positive step size. In this mode, `L`
is the maximum number of components. The model starts with
`min(L_greedy, L)` components and adds `L_greedy` components after each round.
Each larger round uses the preceding fitted model as `model_init`; it does not
restart the model from an empty state. The iteration counter starts again at
one for each round.
It stops when the minimum component log Bayes factor is less than
`greedy_lbf_cutoff`, or when it reaches `L`. The default cutoff is `0.1`.
Leave `L_greedy` unset to keep fixed-L behavior. The verbose task log records
each greedy round. Window QC records the maximum, step, cutoff, and final
number of components.

## Inputs and execution

The preparation workflow consumes a window ID, a tabix-indexed dosage file, a
window association table, and phenotype files with matching modality labels.
The fine-mapping workflow consumes the resulting window table, phenotype data,
covariate files, and one dosage file per window.

Both workflows use container images published to GHCR. See the WDL files for
the complete input schema and defaults, and use MiniWDL or a Cromwell-compatible
engine to run them.
