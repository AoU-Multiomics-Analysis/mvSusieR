# mvSusieR workflows

<!-- workflow-badges:start -->
[![LD-pruning container](https://github.com/evin-padhi/mvSusieR/actions/workflows/docker-image.yml/badge.svg)](https://github.com/evin-padhi/mvSusieR/actions/workflows/docker-image.yml)
[![R lint](https://github.com/evin-padhi/mvSusieR/actions/workflows/r-lint.yml/badge.svg)](https://github.com/evin-padhi/mvSusieR/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/evin-padhi/mvSusieR/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/evin-padhi/mvSusieR/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/evin-padhi/mvSusieR/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/evin-padhi/mvSusieR/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->


## LD-pruning workflow

`workflows/ld_pruning.wdl` creates an LD-pruned sentinel-variant set with
PLINK 2. It accepts the three components of a PLINK 2 data set explicitly, so
the caller does not need to rely on a shared input-file prefix.

The container is published as `ghcr.io/evin-padhi/mvsusier-ld-pruning` after a
successful push to `main`.

### Inputs

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

Example MiniWDL input JSON:

```json
{
  "LDPruning.input_pgen": "gs://example-bucket/genotypes.pgen",
  "LDPruning.input_pvar": "gs://example-bucket/genotypes.pvar",
  "LDPruning.input_psam": "gs://example-bucket/genotypes.psam"
}
```

```bash
miniwdl run workflows/ld_pruning.wdl -i ld_pruning.inputs.json
```

### Outputs

- `pruned_variants`: retained variant IDs (`ld_pruning.prune.in`).
- `excluded_variants`: IDs removed during LD pruning (`ld_pruning.prune.out`).
- `pruned_pgen`, `pruned_pvar`, `pruned_psam`: compact PLINK 2 dataset for
  the retained variants.
- `retained_variant_ids`: the retained IDs written by `plink2 --write-snplist`.
