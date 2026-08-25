# Joint RNA–protein mvSuSiE workflow design

## Goal

Replace the old trans-window path with one production workflow that jointly fine-maps available expression, splicing, and protein outcomes. The workflow must apply one preprocessing contract, learn one data-driven mashr prior from the full locus genotype matrix, and run one joint mvSuSiE model with a reproducible greedy effect-count schedule.

## Scope

The workflow supports these three modalities:

- `expression`
- `splicing`
- `protein`

All three source phenotype files and all three covariate files are required as global workflow inputs. An individual window can contain any nonempty subset of the supported modalities. The workflow does not support `isoform_usage`. It rejects unsupported or duplicated modalities and rejects a window only when no usable outcomes remain.

The workflow supports one locus window per model task. It keeps every usable SNP in that window. It selects trans outcomes separately by modality and adds explicit target-gene expression and splicing outcomes.

## Production interfaces

### Preparation workflow

`workflows/prepare_trans_window.wdl` will replace the modality arrays with explicit inputs:

- `File expression_phenotypes`
- `File splicing_phenotypes`
- `File protein_phenotypes`
- `File target_phenotypes`
- `Int top_n_expression = 25`
- `Int top_n_splicing = 25`
- `Int top_n_protein = 15`

The three phenotype files use the same BED-like layout: chromosome, zero-based start, half-open end, phenotype ID, then one column per sample. `target_phenotypes` is a TSV with `window_id`, `modality`, and `phenotype_id`. Target rows may use only `expression` or `splicing`. Every requested target must exist exactly once in its source phenotype file.

The trans-association table can contain any subset of expression, splicing, and protein rows for the requested window. The preparation task ranks the minimum p-value per phenotype within each modality, selects up to the requested number, adds the explicit target rows, and removes duplicate phenotype selections. A modality can contribute zero outcomes. Preparation fails only when the complete window has no usable outcomes.

The output phenotype manifest stores both the original `phenotype_id` and a unique `outcome_key` formed as `modality::phenotype_id`. Internal matrices and model outputs use `outcome_key`. User-facing tables retain the modality and original phenotype ID.

### Model workflow

`workflows/trans_window_mvsusie.wdl` will require explicit covariate inputs:

- `File expression_covariates`
- `File splicing_covariates`
- `File protein_covariates`

The model workflow will not accept generic modality arrays. Its production defaults are:

- mashr prior;
- five requested PCA factors, limited by the available outcomes and selected SNP rows, with requested and used counts recorded;
- strong-row threshold lfsr of 0.05;
- extreme deconvolution disabled;
- canonical covariance matrices disabled;
- residual covariance initialized with `cov(Y)` and then estimated;
- prior scale and mixture weights estimated from the mashr initialization;
- greedy start `L = 10`;
- greedy step `5`;
- maximum `L = 40`;
- greedy log-Bayes-factor cutoff `1.0`;
- credible-set coverage `0.95`;
- minimum absolute credible-set correlation `0.5`;
- verbose model output enabled.

The WDL and command-line interfaces keep explicit overrides for numeric settings. They do not provide alternative modality or prior modes.

## Joint preprocessing contract

### Sample alignment

The preprocessing step normalizes and validates sample IDs before matrix construction. It requires unique sample IDs in the genotype, phenotype, and covariate inputs. It forms one exact intersection across:

- the locus genotype matrix;
- expression outcomes;
- splicing outcomes;
- protein outcomes;
- expression covariates;
- splicing covariates;
- protein covariates.

The genotype file defines the final row order. Every other matrix is reordered by exact sample ID. The workflow logs input counts, intersection counts, reordered inputs, and dropped sample counts.

For incomplete rows, the workflow computes one row-wise finite-value mask across the selected genotype, phenotype, and covariate matrices. It removes only the affected samples, logs the number removed, and fails if no samples remain. It does not use a scalar completeness flag or an undefined intermediate matrix.

### Covariate adjustment

Each phenotype modality uses only its respective covariates plus an intercept:

- expression outcomes use expression covariates;
- splicing outcomes use splicing covariates;
- protein outcomes use protein covariates.

The genotype covariate matrix is the unique union of all three covariate matrices plus an intercept. If two modalities use the same column name with different values, the workflow preserves both columns and prefixes them with the modality, for example `expression::PC1`, `splicing::PC1`, and `protein::PC1`. If same-name columns have identical aligned values, the union keeps one copy.

The workflow validates alignment before it compares same-name covariates. It records the source modality, original column name, final column name, matrix rank, and a checksum for every covariate column.

### Transformations

The workflow applies these steps on the final shared samples:

1. Remove raw phenotype columns with non-finite or zero variance.
2. Apply rank-based inverse-normal transformation separately to each retained outcome.
3. Residualize each transformed outcome against its modality-specific covariates and an intercept.
4. Remove residualized outcomes with non-finite or zero variance.
5. Center and scale each retained residualized outcome to unit sample variance.
6. Residualize the genotype matrix once against the union covariate model and an intercept.
7. Remove residualized SNP columns with non-finite or zero variance.

The modality vector is filtered with the phenotype matrix at every phenotype-removal step. The prepared bundle records input and retained dimensions, per-modality outcome counts, transformation labels, phenotype covariate ranks, genotype covariate rank, sample IDs, variant IDs, outcome keys, and input checksums.

mvSuSiE performs its own genotype-column standardization. The saved residualized genotype matrix is not standardized a second time during preprocessing.

## Joint mashr prior

The workflow computes marginal `Bhat` and `Shat` matrices for all retained SNPs and all retained joint outcomes. It uses blocked matrix operations and emits progress messages for every block. It stores these matrices in the mashr training RDS and does not write a flattened SNP-by-outcome association table.

Prior learning uses this sequence:

1. Run the one-by-one mash model on all SNP rows.
2. Select rows with minimum lfsr at or below 0.05 for covariance learning.
3. If fewer than five rows pass, use the five rows with the smallest minimum lfsr values.
4. Call `mashr::cov_pca` with `npc = 5` on the selected rows and record the number of covariance matrices it returns.
5. Fit the mash mixture weights on all SNP rows.
6. Convert the fitted mixture to an mvSuSiE mixture prior with no canonical matrices and no extreme-deconvolution step.
7. Use the fitted mash mixture weights to initialize mvSuSiE.

The workflow passes the raw mash covariance matrices to mvSuSiE. It does not
apply the outcome-scale conversion that is used for a fixed prior. The training
bundle records the raw prior, raw covariance range, covariance-training rows,
mixture-training row count, fitted weights, and random seed.

## Joint mvSuSiE model

The workflow initializes the residual covariance with the sample covariance of
the prepared joint outcome matrix. It lets mvSuSiE update the residual
covariance, raw mashr prior scale, and mixture weights. It prints every mvSuSiE
iteration.

The greedy scheduler is a pipeline function with separate `start_L`, `step_L`, and `max_L` arguments. It runs `L = 10, 15, 20, ...` and warm-starts each round from the preceding fit. A round is saturated when its minimum component log Bayes factor is below `1.0`. The scheduler stops at the first saturated round or at `L = 40`.

Only the final fit is saved. Intermediate fits remain in memory for warm starts and are discarded after the next round starts. The scheduler writes one history row per round with:

- round number;
- requested and fitted `L`;
- iteration count;
- minimum component log Bayes factor;
- reported credible-set count;
- number of components with at least one feature lfsr below 0.05;
- maximum component-level variant probability;
- stopping action.

Credible-set purity, feature lfsr, and per-outcome log Bayes factors are diagnostics. They do not change the greedy stopping rule.

## Outputs

The model task returns:

- `prepared_joint_window.rds`;
- `mashr_training_bundle.rds` containing `Bhat` and `Shat`;
- `mvsusie_fit_bundle.rds` containing only the final fit;
- `greedy_L_history.tsv`;
- `variant_pip.tsv.gz`;
- `credible_sets.tsv.gz`;
- `credible_set_members.tsv.gz`;
- `component_feature_support.tsv.gz` with lfsr and per-outcome log BF;
- `window_qc.tsv`;
- `covariate_provenance.tsv.gz`;
- `session_info.txt`;
- complete preparation and model logs;
- credible-set-by-feature PNG and PDF plots.

The plot task must call `mvsusieR::mvsusie_plot` with `conditional_effect = TRUE`, `add_cs = TRUE`, and the final model outcome labels. It must not reconstruct effects in a custom plotting layer.

Summary code accepts either `fit$mu2_diag` or the legacy `fit$mu2` field. It fails if neither field exists.

## Errors and diagnostics

The workflow fails with a specific message when:

- a required modality or covariate file is absent;
- `isoform_usage` or another modality is present;
- phenotype or sample IDs are duplicated;
- a target phenotype is absent or duplicated;
- sample alignment cannot be completed;
- no complete shared samples remain;
- no outcomes remain for a required modality;
- no usable SNPs remain;
- mash covariance dimensions do not match the joint outcome count;
- a mash or mvSuSiE result contains non-finite values;
- an mvSuSiE round does not converge;
- plotting does not return an effect plot.

Every WDL command logs task start, resolved input counts, each preprocessing stage, each mashr stage, each greedy round, saved output paths, and task completion.

## Dependency reproducibility

The model container pins these upstream commits:

- `mvsusieR`: `ebd1133953005fa70c6b338727b5fe9222e2a1c2`
- `susieR`: `65f3586a865fb6748cb4f9df50510ac577706348`

The container smoke test verifies the installed commit metadata or an equivalent immutable source reference, required function arguments, and required fit fields. Local development does not build the Docker image.

## Test and automation strategy

Tests use one synthetic joint fixture with expression, splicing, and protein outcomes. They also test windows that omit each modality in turn. The fixture includes:

- different sample order in each input;
- different `PC1` values for all three modalities;
- one identical covariate column shared across modalities;
- incomplete samples in phenotype and covariate inputs;
- one zero-variance phenotype;
- enough associations to test expression 25, splicing 25, and protein 15 selection;
- explicit expression and splicing target outcomes;
- a small genotype matrix with known covariate relationships.

The R tests verify exact sample intersection and order, outcome keys,
per-modality selection counts, correct phenotype residualization, genotype
residualization against the full covariate union, phenotype unit variance,
variant filtering, the raw mashr prior handoff, enabled prior and residual
updates, greedy warm starts, the `1.0` stopping threshold, final-only fit
storage, summary tables, and mvSuSiE API plot creation.

Negative tests verify rejection of an empty outcome set, `isoform_usage`, duplicate IDs, missing targets, mismatched covariate samples, and invalid greedy settings.

WDL tests validate explicit joint inputs, per-modality top counts, production defaults, command logging, outputs, and runtime settings. GitHub Actions runs R lint, R integration tests, WDL validation, and the model-container smoke test. No local Docker build is required.

## Migration

The existing production preparation and model workflow names remain stable. Documentation and examples show the required expression, splicing, and protein source inputs. A window can use one, two, or three of these modalities. Isoform-usage inputs fail validation.

The joint workflow supersedes the standalone exploratory scripts used for the IKZF1 and CREB5 protein comparisons. Those analysis outputs remain outside the repository and are not added to Git.
