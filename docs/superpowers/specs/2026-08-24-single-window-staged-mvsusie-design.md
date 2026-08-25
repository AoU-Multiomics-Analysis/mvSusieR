# Single-window staged mvSuSiE workflow design

## Objective

Change the joint mvSuSiE workflow so that one workflow invocation processes
one genomic window. Separate short, memory-intensive preprocessing from the
long model fit. Permit an existing prepared-window RDS file to bypass
preprocessing.

Also change phenotype preparation so that expression, splicing, and protein
inputs can have different sample sets. The output phenotype file must contain
the shared samples in one consistent order.

## Scope

This change covers:

- the `TransWindowMvSusie` WDL interface and task structure;
- the preprocessing and fitting R entry points;
- phenotype sample alignment in `PrepareWindowPhenotypes`;
- preparation, fitting, WDL, and container smoke tests;
- workflow documentation.

This change does not alter the statistical transformations, mashr prior,
greedy-L schedule, mvSuSiE settings, summary calculations, or plotting API.
It does not create one fit RDS per greedy-L round.

## Single-window workflow interface

`TransWindowMvSusie` will not scatter. The external scheduler will launch one
workflow invocation for each window.

The workflow will accept:

- one `String window_id`;
- one optional `File? prepared_window`;
- optional raw-preparation inputs: the window manifest, phenotype manifest,
  dosage file, combined phenotype data, three covariate files, and optional
  keep-samples file;
- the existing preprocessing thresholds and model controls.

The workflow will have two input modes.

### Raw-input mode

When `prepared_window` is absent, every raw-preparation input except the
keep-samples file is required. A validation task will receive only Boolean
`defined(...)` values. It will stop with a clear message that lists missing
raw inputs. It will not localize the large input files.

After validation, `PrepareMvSusieInput` will run `prepare_window.R` and write
one `prepared_window.rds`.

### Prepared-input mode

When `prepared_window` is present, the workflow will skip
`PrepareMvSusieInput`. Extra raw inputs, if supplied, will not be used. The
prepared RDS will go directly to `FitMvSusie`.

The fit entry point will validate that the RDS contains the required window,
genotype, phenotype, metadata, sample, QC, and covariate-provenance fields. It
will also verify that the RDS window ID equals the requested `window_id`.

## Task and memory boundaries

The workflow will use these tasks:

1. `ValidateMvSusieInputs`: small Boolean-only validation task.
2. Conditional `PrepareMvSusieInput`: 16 GiB memory and 2 CPUs.
3. `FitMvSusie`: 8 GiB memory and 2 CPUs.
4. `SummarizeMvSusie`: unchanged statistical behavior.
5. `PlotMvSusie`: unchanged and continues to use the mvSuSiE plotting API.

The preprocessing task will end after it writes the prepared RDS. Thus, raw
tables, transposed matrices, and residualization copies cannot remain in the
long-running fit process.

`FitMvSusie` will load only the prepared RDS. It will run the all-SNP
association calculations, PCA-only mashr training, and every greedy-L round
in the same process. The 8 GiB request gives headroom over the observed
approximately 3.7 GiB R allocation for an IKZF1 window with 8,895 samples and
23,393 variants.

## Outputs

The single-window workflow will return scalar files:

- `prepared_window`;
- `mvsusie_fit`;
- `mashr_training`;
- `greedy_L_history`;
- `covariate_provenance`;
- fit standard output, standard error, and session information;
- variant PIP, credible-set, credible-set-member, feature-support, and QC
  tables;
- PNG, PDF, and RDS effect plots.

The fit task will always write `covariate_provenance` from the prepared RDS.
This keeps the output contract identical in raw-input and prepared-input
modes.

The per-window workflow will not merge results across windows. The existing
merge script remains available for an outer aggregation stage.

## Phenotype sample intersection

`PrepareWindowPhenotypes` will align samples after it selects outcomes and
before it combines modality tables.

For each modality that contributes at least one outcome, the preparation code
will:

1. Treat the first four columns as phenotype metadata and the remaining
   columns as samples.
2. Normalize sample IDs with the downstream rule: trim surrounding spaces and
   remove an R-added leading `X` when it precedes a digit.
3. Stop if normalization creates duplicate sample IDs within a modality.
4. Calculate the intersection across contributing modalities.
5. Preserve the normalized sample order from the first contributing modality
   in expression, splicing, protein order.
6. Subset and reorder every contributing table to this common order.
7. Write normalized sample IDs in the combined phenotype file.

A modality with no selected outcome will not reduce the sample intersection.
The code will stop if contributing modalities have no shared samples.

The log will report, for each contributing modality:

- the input sample count;
- the shared sample count;
- the number of samples removed.

The window QC table will add these same counts so the alignment is recorded in
a machine-readable output.

Downstream joint preprocessing will continue to intersect the prepared
phenotype samples with genotype and covariate samples. The phenotype
preparation step will not inspect genotype or covariate inputs.

## Logging and failures

Every WDL command will retain timestamped progress messages. New messages will
identify the input mode, preprocessing start and completion, prepared-RDS
validation, fit start and completion, and requested memory boundary.

Failures will identify:

- missing raw inputs when no prepared RDS is supplied;
- a prepared RDS with missing fields;
- a prepared RDS for the wrong window;
- duplicate normalized sample IDs;
- zero shared phenotype samples;
- missing or empty task outputs.

## Tests

Tests will cover:

- phenotype inputs with different sample membership and order;
- normalization of R-added sample-name prefixes;
- exact sample-value alignment after reordering;
- exclusion of a modality with no selected outcomes from the intersection;
- failure when contributing modalities have no shared samples;
- sample-count fields in logs and window QC;
- WDL single-window scalar inputs and outputs;
- absence of the model-workflow scatter;
- conditional preprocessing when `prepared_window` is absent;
- direct fitting when `prepared_window` is present;
- clear failure for incomplete raw-input mode;
- prepared-RDS schema and window-ID validation;
- 16 GiB preprocessing and 8 GiB fit runtime requests;
- unchanged mashr, greedy-L, summary, and plotting behavior.

GitHub Actions will build the model image and run both workflow modes with
small fixtures. Local smoke tests will validate R and WDL behavior without
building a container.

## Compatibility

This is an intentional WDL interface change. Existing callers that pass a
multirow window manifest to `TransWindowMvSusie` must instead launch one job
per window. Array outputs become scalar outputs. `PrepareTransWindow` remains
a single-window workflow and its primary output names remain unchanged.

Existing prepared-window RDS files can be reused if they pass the new schema
and window-ID checks.
