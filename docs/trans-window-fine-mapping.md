# Joint trans-window fine-mapping

The pipeline has three workflows. `IndexPhenotypeBed` creates reusable
phenotype indexes. `PrepareTransWindow` creates the locus dosage and joint
phenotype inputs. `TransWindowMvSusie` processes one locus per workflow
launch. It preprocesses expression, splicing, and protein data and fits one
joint mvSuSiE model.

## Index phenotype files once

Run
[`workflows/index_phenotype_bed.wdl`](../workflows/index_phenotype_bed.wdl)
once for each source expression, splicing, and protein phenotype file. The
input is a BED-like table with chromosome, zero-based start, half-open end,
and phenotype ID in the first four columns. Sample columns follow these
metadata columns. The input metadata names do not have to use canonical names.

Each launch writes four outputs:

- `indexed_phenotypes`: a coordinate-sorted BGZF phenotype file;
- `phenotype_tbi`: its tabix index;
- `phenotype_lookup`: a compressed table with phenotype ID and coordinates;
- `index_qc`: row counts, sample counts, file sizes, and elapsed time.

The index workflow uses the coordinates in the phenotype file. It does not
use a GTF. Thus, the same workflow supports genes, splicing events, and
proteins. Keep the four outputs for reuse across all fine-mapping windows.

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
window_id  outcome_key  phenotype_id  modality  phenotype_file  p_value
```

`outcome_key` has the form `modality::phenotype_id`. This key prevents name
collisions between modalities. `p_value` is the minimum window association P
value for a selected trans feature. A target-only feature has a missing
`p_value`. This one manifest is an input to both mvSuSiE and checkpointed
univariate SuSiE. mvSuSiE uses all rows. Univariate SuSiE uses rows with a
finite `p_value`.

A window can retain any nonempty subset of the three supported modalities. A
top-N value is a maximum, not a minimum. Thus, a modality can contribute fewer
than N outcomes or no outcomes in a window. The workflow rejects a window only
when it has no usable outcomes from any modality.

The phenotype task uses the intersection of samples in the nonempty assay
files. It reports the input, retained, and removed sample counts for each
contributing assay. A missing modality does not reduce the sample set.

Each preparation output starts with the `window_id`. For example, window
`chr7_50000000_52000000` writes:

```text
chr7_50000000_52000000.window_dosage.tsv
chr7_50000000_52000000.window_manifest.tsv
chr7_50000000_52000000.window_phenotypes.tsv
chr7_50000000_52000000.window_phenotypes.bed.gz
chr7_50000000_52000000.window_qc.tsv
```

The `window_id` can contain letters, numbers, periods, underscores, and
hyphens. The workflow rejects other characters before it writes output files.

### Use indexed phenotype access

`PrepareTransWindow` has two optional companion inputs for each modality:

```text
expression_phenotypes_tbi       expression_phenotype_lookup
splicing_phenotypes_tbi         splicing_phenotype_lookup
protein_phenotypes_tbi          protein_phenotype_lookup
```

For indexed access, set the main modality phenotype input to the corresponding
`indexed_phenotypes` output. Set its TBI input to `phenotype_tbi` and its
lookup input to `phenotype_lookup`. Both companions must be present. The task
stops with an error if only one companion is present.

If both companions are absent, the task reads the complete phenotype file.
This rule preserves existing input JSON files. It also permits mixed access.
For example, expression and protein can use tabix while splicing uses a full
scan.

The task makes one batched tabix query per indexed, contributing modality. It
then applies an exact phenotype-ID filter because nearby or overlapping
intervals can add records to a coordinate query. Full-scan and tabix access
use the same feature selection, sample intersection, output writer, and error
checks.

The window-prefixed `window_qc.tsv` reports these access fields for each
modality:

```text
access_method  lookup_seconds  query_seconds  parse_seconds  query_rows  modality_seconds
```

The task log also reports the requested, queried, and retained phenotype
counts. Use these values to compare tabix and full-scan performance on the
source data.

## Preprocess the joint data

[`workflows/trans_window_mvsusie.wdl`](../workflows/trans_window_mvsusie.wdl)
has two input modes. In raw mode, it requires one combined phenotype file and
three explicit covariate files:

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

In prepared mode, set `prepared_window` to an existing prepared RDS. The
workflow skips preprocessing and starts with model fitting. The raw genotype,
phenotype, and covariate inputs are optional in this mode. The workflow checks
the RDS structure and the `window_id` before it fits the model.

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
add another covariance family. It uses the fitted mash mixture weights to
initialize mvSuSiE. It also preserves mashr's fitted point-mass weight as the
initial mvSuSiE null weight.

The workflow passes the raw fitted mash covariance matrices to mvSuSiE. This
handoff is required when mvSuSiE estimates the prior scale. The workflow does
not apply the fixed-prior outcome-scale conversion.

## Fit mvSuSiE

The workflow initializes the residual covariance with the covariance of the
prepared outcome matrix. mvSuSiE then updates the residual covariance, the raw
mash prior scale, the mash mixture weights, and the null weight. mvSuSiE runs
with verbose output.

The greedy schedule starts at `L = 10`, increases in steps of 5, and stops at
`L = 40`. It stops earlier when the minimum component log Bayes factor is less
than 1. Each larger round uses the preceding fit as `model_init`. The iteration
counter starts at one in each round, but the model does not restart from an
empty state. The workflow saves only the final fit and writes one compact row
per round to `greedy_L_history.tsv`.

## Outputs

For one window, the workflow writes scalar outputs:

- the prepared joint-data bundle;
- the mashr training bundle with `Bhat`, `Shat`, covariance inputs, and the raw
  mash prior;
- the final mvSuSiE fit;
- the greedy-L history;
- covariate provenance, run logs, and R session information;
- variant PIPs;
- one summary row per credible set;
- the member variants for each credible set;
- component support by feature, including lfsr and outcome log Bayes factor;
- window QC;
- PNG, PDF, and RDS plot outputs.

Each file starts with the `window_id`. For window `w1`, the file basenames are:

```text
w1.prepared_window.rds
w1.mvsusie_fit.rds
w1.mashr_training.rds
w1.greedy_L_history.tsv
w1.covariate_provenance.tsv.gz
w1.variant_pip.tsv.gz
w1.credible_sets.tsv.gz
w1.credible_set_members.tsv.gz
w1.component_feature_support.tsv.gz
w1.window_qc.tsv
w1.effect_plot.png
w1.effect_plot.pdf
w1.effect_plot.rds
w1.run.stdout.log
w1.run.stderr.log
w1.session_info.txt
```

Raw mode and prepared mode use the same output names. The workflow output
variable names do not change.

The workflow does not merge results across windows. Run
`scripts/merge_window_outputs.R` after all window jobs finish. The model does
not write the full component-by-variant-by-feature posterior tensor.

The effect plot comes from `mvsusieR::mvsusie_plot` with
`conditional_effect = TRUE`. The pipeline does not transform the effect values
after this API call. Thus, the figure keeps the mvSuSiE effect scale.

## Reproducibility

The raw preparation task requests 16 GiB of memory. The long-running fit task
requests 8 GiB. Thus, a scheduler does not need to hold preparation memory for
the full model run. The summary and plot tasks remain separate.

The model image uses micromamba with strict channel priority. The `dnachun`
channel supplies mvSuSiER 0.3.0 and susieR 0.15 or newer. The image build
checks both installed versions before it copies the pipeline scripts. GitHub
Actions builds the image and runs small raw-input and prepared-window WDL
smoke tests. Local smoke tests do not build the image.

All workflows use WDL 1.0. You can run them with MiniWDL or a
Cromwell-compatible engine. Dockstore keeps the existing workflow names and
descriptor paths.
