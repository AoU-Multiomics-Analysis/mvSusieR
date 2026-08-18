# LD-pruning workflow

`workflows/ld_pruning.wdl` produces an LD-pruned sentinel-variant set with
PLINK 2. It accepts the three components of a PLINK 2 data set explicitly, so
callers do not need to construct or localize a shared filename prefix.

The workflow first applies the requested variant and optional sample filters,
then runs `plink2 --indep-pairwise`. It verifies that at least one variant is
retained and extracts the retained IDs into a compact PGEN/PVAR/PSAM dataset.

## Container

The task runs in `ghcr.io/evin-padhi/mvsusier-ld-pruning`. The GitHub Action
builds and smoke-tests this image for pull requests, and publishes it after a
successful push to `main`.

## Inputs

| Input | Required | Default | Description |
| --- | --- | --- | --- |
| `input_pgen` | Yes | — | PLINK 2 genotype file. |
| `input_pvar` | Yes | — | Matching PLINK 2 variant file. |
| `input_psam` | Yes | — | Matching PLINK 2 sample file. |
| `keep_samples` | No | — | PLINK two-column sample list passed to `--keep`. |
| `exclude_variants` | No | — | Variant-ID list passed to `--exclude`. |
| `maf` | No | `0.05` | Minimum minor-allele frequency. |
| `geno` | No | `0.01` | Maximum variant missingness. |
| `ld_window_kb` | No | `1000` | LD-pruning window size in kilobases. |
| `ld_step_variants` | No | `50` | Number of variants advanced per pruning step. |
| `ld_r2` | No | `0.1` | Pairwise LD threshold. |

## Example input

Save the following as `ld_pruning.inputs.json` and replace the example URIs
with your PGEN, PVAR, and PSAM paths.

```json
{
  "LDPruning.input_pgen": "gs://example-bucket/genotypes.pgen",
  "LDPruning.input_pvar": "gs://example-bucket/genotypes.pvar",
  "LDPruning.input_psam": "gs://example-bucket/genotypes.psam"
}
```

Run it locally with MiniWDL:

```bash
miniwdl run workflows/ld_pruning.wdl -i ld_pruning.inputs.json
```

## Outputs

- `pruned_variants`: retained variant IDs (`ld_pruning.prune.in`).
- `excluded_variants`: IDs removed during LD pruning (`ld_pruning.prune.out`).
- `pruned_pgen`, `pruned_pvar`, `pruned_psam`: compact PLINK 2 dataset for
  the retained variants.
- `retained_variant_ids`: retained IDs written by `plink2 --write-snplist`.
