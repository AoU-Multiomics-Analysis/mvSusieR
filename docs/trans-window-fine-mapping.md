# Joint trans-window fine-mapping

The pipeline has two workflows. `PrepareTransWindow` creates the locus dosage
and joint phenotype inputs. `TransWindowMvSusie` preprocesses expression,
splicing, and protein data and fits one mvSuSiE model for each locus.

## Prepare the window data

[`workflows/prepare_trans_window.wdl`](../workflows/prepare_trans_window.wdl)
has two tasks:

- `PrepareWindowGenotypes` extracts every dosage record in the requested locus
  from a tabix-indexed genome-wide dosage file.
- `PrepareWindowPhenotypes` selects the top trans features for each modality
  and adds the requested target features.

The default trans-feature counts are 25 expression, 25 splicing, and 15
protein. You can set each count independently. The target-feature table can
contain expression and splicing features. Its required columns are:

```text
window_id  modality  phenotype_id
```

The workflow requires separate BED-like phenotype files for expression,
splicing, and protein. Each file has four metadata columns followed by one
column for each sample. The workflow writes one combined phenotype file and a
manifest with these columns:

```text
window_id  outcome_key  phenotype_id  modality  phenotype_file
```

`outcome_key` has the form `modality::phenotype_id`. This key prevents name
collisions between modalities. A window can retain any nonempty subset of the
three supported modalities. A top-N value is a maximum, not a minimum. Thus,
a modality can contribute fewer than N outcomes or no outcomes in a window.
The workflow rejects a window only when it has no usable outcomes from any
modality.

## Preprocess the joint data

[`workflows/trans_window_mvsusie.wdl`](../workflows/trans_window_mvsusie.wdl)
requires one combined phenotype file and three explicit covariate files:

- `expression_covariates`
- `splicing_covariates`
- `protein_covariates`

The workflow keeps samples present in the genotype, phenotype, and all three
covariate inputs. All three source phenotype and covariate files remain global
workflow inputs, even when one modality has no outcomes in a given window. It
removes samples with non-finite values. It applies a rank-based inverse-normal
transform to each outcome. It then regresses each present outcome on the
covariates for its own modality and scales the residual to unit variance.

The genotype matrix uses the aligned union of all expression, splicing, and
protein covariates. If two covariates have the same name and values, the
workflow keeps one copy. If they have the same name but different values, the
workflow keeps both and adds modality prefixes. The workflow writes this
decision and an MD5 checksum for each source column to
`covariate_provenance.tsv.gz`.

## Learn the mashr prior

The workflow calculates `Bhat` and `Shat` for every retained SNP and every
joint outcome with matrix operations. It uses a one-by-one mashr fit to select
strong SNP rows at lfsr 0.05. If fewer than five rows pass, it uses the five
rows with the smallest lfsr values. It requests five PCA covariance matrices,
but limits the count to the number of outcomes and selected SNP rows. For one
outcome, it uses the equivalent one-dimensional covariance because PCA cannot
add extra covariance directions. It then fits the mash mixture on all retained
SNPs.

The workflow supplies only the PCA covariance matrices to mashr. It does not
add another covariance family. It keeps the fitted mash mixture weights for
mvSuSiE.

Before mvSuSiE uses the mash prior, the workflow divides each covariance
element by the product of the corresponding outcome standard-error scales.
This conversion cancels mvSuSiE's automatic outcome-scale conversion and
preserves the mash effect covariance. It does not standardize the prior twice.

## Fit mvSuSiE

The workflow initializes the residual covariance with the covariance of the
prepared outcome matrix and keeps it fixed. It also keeps the mash prior and
mixture weights fixed. mvSuSiE runs with verbose output.

The greedy schedule starts at `L = 10`, increases in steps of 5, and stops at
`L = 40`. It stops earlier when the minimum component log Bayes factor is less
than 1. Each larger round uses the preceding fit as `model_init`. The iteration
counter starts at one in each round, but the model does not restart from an
empty state. The workflow saves only the final fit and writes one compact row
per round to `greedy_L_history.tsv`.

## Outputs

For each window, the workflow writes:

- the prepared joint-data bundle;
- the mashr training bundle with `Bhat`, `Shat`, covariance inputs, and scale
  information;
- the final mvSuSiE fit;
- the greedy-L history;
- covariate provenance, run logs, and R session information;
- variant PIPs;
- one summary row per credible set;
- the member variants for each credible set;
- component support by feature, including lfsr and outcome log Bayes factor;
- window QC;
- PNG, PDF, and RDS plot outputs.

The workflow also merges the tabular summaries across windows. It does not
write the full component-by-variant-by-feature posterior tensor.

The effect plot comes from `mvsusieR::mvsusie_plot` with
`conditional_effect = TRUE`. The pipeline does not transform the effect values
after this API call. Thus, the figure keeps the mvSuSiE effect scale.

## Reproducibility

The model image pins susieR commit
`65f3586a865fb6748cb4f9df50510ac577706348` and mvsusieR commit
`ebd1133953005fa70c6b338727b5fe9222e2a1c2`. GitHub Actions builds the image
and runs a small joint preparation, prior, model, summary, and plot test. Local
smoke tests do not build the image.

Both workflows use WDL 1.0. You can run them with MiniWDL or a
Cromwell-compatible engine. Dockstore keeps the existing workflow names and
descriptor paths.
