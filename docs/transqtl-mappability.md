# TransQTL mappability generation

`workflows/transqtl_mappability.wdl` creates a low-mappability BED for use by
the TransQTL BAM-filtering workflow. It runs GenMap on a reference FASTA,
selects positions below the configured mappability threshold, and merges the
resulting intervals with BEDTools.

## Important compatibility rule

The reference FASTA must have the same genome build and contig naming as the
BAMs that will later be filtered. The workflow does not download a FASTA or
choose a GRCh38 contig set for you. Decide whether the reference contains only
primary chromosomes or also alternate/unplaced contigs before running it, and
use the same choice for downstream alignment files.

The generated BED represents the configured GenMap `(k,e)` definition. It is
not automatically the ENCODE 36-mer track: use `kmer_length = 36` and a
matching reference if that exact definition is required. The default `k=146`
is intended for 146-base reads. The default `max_mismatches = 2` follows the
two-mismatch criterion discussed for Trans-PCO-style filtering.

## Inputs and outputs

Required input:

- `reference_fasta`: reference FASTA used to build the GenMap index.

Optional inputs and defaults:

| Input | Default | Meaning |
| --- | ---: | --- |
| `kmer_length` | `146` | k-mer/read length used for mappability |
| `max_mismatches` | `2` | mismatches allowed when assessing uniqueness |
| `mappability_threshold` | `1.0` | emit intervals with score `<` this value |
| `output_prefix` | `GRCh38.K146.m2` | output filename prefix |
| `threads` | `16` | task CPU allocation |
| `preemptible_tries` | `3` | retry count for preemptible execution |

The workflow emits:

- `<prefix>.low_mappability.bed`: three-column, sorted/merged BED for
  `samtools view -L`.
- `<prefix>.mappability.bedGraph`: full GenMap mappability output.
- `<prefix>.mappability.metadata.tsv`: parameters, reference size, reference
  path, reference SHA-256, and interval count.

An empty low-mappability BED is valid if no positions fall below the threshold;
the full BEDGraph must still be nonempty.

## Example MiniWDL launch

Create an inputs JSON such as:

```json
{
  "TransQTLMappability.reference_fasta": "GRCh38.primary.fa",
  "TransQTLMappability.kmer_length": 146,
  "TransQTLMappability.max_mismatches": 2,
  "TransQTLMappability.mappability_threshold": 1.0,
  "TransQTLMappability.output_prefix": "GRCh38.K146.m2",
  "TransQTLMappability.threads": 16,
  "TransQTLMappability.preemptible_tries": 3
}
```

Run it with:

```bash
miniwdl run workflows/transqtl_mappability.wdl -i transqtl_mappability.inputs.json
```

The container is published as
`ghcr.io/aou-multiomics-analysis/transqtl-mappability` by
`.github/workflows/transqtl-mappability-image.yml`. It uses micromamba with
`genmap` and `bedtools` installed from `conda-forge` and `bioconda`.

## Resource expectations

GenMap indexes the entire supplied reference and can require substantial RAM
and temporary/local disk. The WDL requests 64 GiB memory and a 200 GB local
SSD by default. A human genome run is expected to be much more expensive than
the CI smoke test, which only checks the tool interface on a tiny synthetic
reference.
