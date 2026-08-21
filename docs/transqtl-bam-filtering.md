# TransQTL BAM filtering and RNA-SeQC2

`workflows/transqtl_bam_filtering.wdl` filters aligned RNA-seq reads at the
read-template level and then requantifies the retained BAM with RNA-SeQC2 in
the same task. Keeping both operations together avoids relocating the filtered
BAM between jobs.

## Filtering semantics

The task performs these steps:

1. Normalize and sort the supplied low-mappability BED.
2. Collect template names where any alignment has a missing or non-unique
   `NH` tag. Only `NH:i:1` passes this criterion.
3. Collect template names where any alignment overlaps the BED in reference
   coordinates, using `samtools view -L`.
4. Remove every alignment record belonging to the union of those template
   names. Thus, a paired-end template is excluded as a whole when either mate
   triggers a filter.
5. Validate/index the filtered BAM. Coordinate-sorted inputs are retained as
   they are; `samtools sort` is used only if indexing shows that resorting is
   necessary.
6. Run RNA-SeQC2 on the retained BAM.

The low-mappability filter therefore acts on aligned records and their
reference-coordinate overlap. It does not remove genes based on a gene-level
mappability score, and it does not compare the read sequence to the BED.

## Inputs

Required inputs are:

- `input_bam` and `input_bai`;
- `low_mappability_bed` generated for the same reference/contig set;
- `genes_gtf`; and
- `sample_id`.

Other defaults:

- `strandedness = "rf"`; use `fr` when appropriate for the library.
- `legacy = false`; set `true` only when the RNA-SeQC 1.1.9-compatible
  behavior is specifically required.
- `threads = 4`.
- `mappability_threshold = 1.0`, recorded in the metrics as provenance.
- `write_filter_metrics = true`; set `false` to skip the filtering summary
  calculations and leave the optional `filter_metrics` output undefined.
- `preemptible_tries = 1`; increase this to allow retries on preemptible/spot
  instances when supported by the execution backend.

## Outputs

The workflow emits:

- coordinate-validated filtered BAM and BAI;
- excluded template names;
- filtering metrics and filtered BAM `flagstat` output; and
- RNA-SeQC2 gene-level read counts (`gene_reads.gct`), TPMs (`gene_tpm.gct`),
  metrics, and coverage tables.

The filtering metrics file is emitted only when `write_filter_metrics` is
`true`.

The filtering metrics distinguish alignment records from read templates and
report totals before/after filtering. Filter-specific counts can overlap:
templates or records counted in both the non-unique and low-mappability
categories are not counted twice in the total removed count.

## Example launch

```bash
miniwdl run workflows/transqtl_bam_filtering.wdl \
  -i transqtl_bam_filtering.inputs.json
```

The BAM-filtering image is published as
`ghcr.io/aou-multiomics-analysis/transqtl-bam-filtering` by
`.github/workflows/transqtl-bam-filtering-image.yml`.
