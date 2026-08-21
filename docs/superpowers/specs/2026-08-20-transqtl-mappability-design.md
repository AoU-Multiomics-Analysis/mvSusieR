# TransQTL mappability generation

## Goal

Add a reproducible WDL workflow that generates a read-length-specific whole-genome low-mappability BED for TransQTL BAM filtering, using GenMap with configurable k-mer length, mismatch allowance, and score threshold.

## Scope

The feature adds four pieces:

1. A WDL workflow named `TransQTLMappability`.
2. A micromamba container containing GenMap and BED-processing utilities.
3. GitHub Actions validation, container smoke testing, and GHCR publication.
4. Dockstore registration and repository tests for the workflow and metadata.

The workflow will calculate genomic mappability from a user-supplied reference FASTA. It will not download or select a reference automatically, because the reference contig set and naming must match the BAMs that will later be filtered.

## Design

### Mappability calculation

GenMap will build an index from the input FASTA and compute `(k,e)` mappability:

- `kmer_length`: default `146`, representing the individual read length.
- `max_mismatches`: default `2`.
- `mappability_threshold`: default `1.0`; intervals with a score below this value are emitted as low mappability.

The workflow will run GenMap with its text, WIG, and BEDGraph outputs enabled. It will normalize the generated BEDGraph, select records with score `< mappability_threshold`, sort and merge adjacent/overlapping intervals, and write a three-column BED suitable for `samtools view -L`. The full BEDGraph must be nonempty; an empty low-mappability BED is valid for a reference with no intervals below the threshold.

### WDL interface

Workflow inputs:

- `File reference_fasta`
- `Int kmer_length = 146`
- `Int max_mismatches = 2`
- `Float mappability_threshold = 1.0`
- `String output_prefix = "GRCh38.K146.m2"`
- `Int threads = 16`
- `Int preemptible_tries = 3`

Workflow outputs:

- `File low_mappability_bed`
- `File mappability_bedgraph`
- `File mappability_metadata`

The task will validate positive k-mer length, nonnegative mismatch count, threshold in a valid range, and at least one thread. It will log the major index, map, conversion, and output-validation stages.

### Container

The Docker image will use `mambaorg/micromamba:2.3.3`, strict channel priority, and the existing repository convention:

```text
conda-forge
bioconda
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
```

The environment will install `genmap` and `bedtools` through the Conda channels. The image build will verify both executables.

### Continuous integration and Dockstore

The new GitHub Actions workflow will:

- run on changes to the mappability Dockerfile, environment, WDL, tests, or workflow itself;
- run `miniwdl check` on the WDL;
- build the image for a local smoke test;
- run GenMap and BED validation on a small synthetic reference;
- publish the image to `ghcr.io/${{ github.repository_owner }}/transqtl-mappability` on non-PR events.

`.dockstore.yml` will register `/workflows/transqtl_mappability.wdl` under the name `transqtl-mappability`.

## Testing

Static tests will verify:

- the WDL parses with miniwdl;
- the expected workflow/task names, parameters, GenMap commands, outputs, logging, and runtime settings are present;
- the environment uses conda-forge and bioconda and contains GenMap and bedtools;
- Dockstore points to the correct WDL;
- the container smoke test can produce a nonempty GenMap BEDGraph and process an interval with BEDTools from a tiny reference.

## Non-goals

- Reproducing the exact byte-level ENCODE/GEM K36 track.
- Modeling splice-junction mappability for RNA-seq reads.
- Automatically downloading a GRCh38 FASTA or choosing primary/ALT contigs.
- Uploading the generated BED to a public release location.
