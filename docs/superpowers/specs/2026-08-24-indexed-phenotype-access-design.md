# Indexed phenotype access design

## Objective

Reduce the phenotype-loading time in `PrepareTransWindow`. The current reader
decompresses and parses each complete wide phenotype file for every window.
The new path must retrieve only the selected phenotype rows while it preserves
the current selection, sample-intersection, manifest, and QC behavior.

## Scope

This change applies to expression, splicing, and protein phenotype inputs. It
adds reusable phenotype indexing and an optional indexed access path to
`PrepareTransWindow`.

This change does not alter:

- top-N selection by modality;
- target phenotype selection;
- the sample intersection across contributing modalities;
- phenotype transformations or covariate regression;
- mashr or mvSuSiE fitting;
- the combined phenotype file and manifest formats.

The design does not use a GTF. It uses the coordinates in each phenotype file.
Thus, it supports gene expression, splicing events, and proteins through one
interface.

## Reusable phenotype index workflow

Add a workflow that prepares one phenotype file per invocation. The workflow
accepts:

- a phenotype BED-like file;
- a modality label for logging and QC.

The input must have four metadata columns followed by sample columns. The four
metadata columns are chromosome, zero-based start, half-open end, and phenotype
ID. Their input names can differ.

The workflow produces:

- a coordinate-sorted BGZF phenotype file;
- a tabix index for the BGZF file;
- a lookup table with `phenotype_id`, `chrom`, `start`, and `end`;
- an index QC table.

The workflow runs once for each source phenotype file. Its outputs can be used
by all fine-mapping windows.

The index workflow must use streaming commands where possible. It must not load
the full phenotype matrix into R. It must log the input file, modality, row
count, sample count, coordinate validation, BGZF compression, tabix indexing,
lookup generation, elapsed time, and output sizes.

The workflow must reject:

- fewer than five columns;
- empty or duplicate phenotype IDs;
- invalid coordinates;
- records that are not in coordinate order after the indexing preparation;
- an empty data file;
- an index that cannot answer a test query.

## PrepareTransWindow interface

Keep the existing phenotype file inputs:

- `expression_phenotypes`;
- `splicing_phenotypes`;
- `protein_phenotypes`.

Add two optional companion inputs for each modality:

- the tabix index;
- the phenotype coordinate lookup table.

The indexed path is active for a modality when both companion inputs are
present. The full-scan path is active when neither companion input is present.
The workflow must stop with a clear error when only one companion input is
present.

This optional interface preserves existing input JSON files. A caller can move
one modality at a time to indexed access. A window can use indexed access for
one modality and full-scan access for another modality.

## Indexed phenotype selection

Top-N and target selection occurs before phenotype matrix extraction. For each
modality on the indexed path, the preparation code must:

1. Read the small coordinate lookup table.
2. Match all requested phenotype IDs.
3. Stop if a requested target or selected trans phenotype is absent.
4. Write one temporary BED region list with the matched coordinates.
5. Make one batched `tabix -R` query against the localized BGZF file.
6. Parse only the returned rows and the original header.
7. Filter by exact phenotype ID because an interval query can return nearby or
   overlapping phenotype records.
8. Remove duplicate query results and restore the requested phenotype order.
9. Stop if any requested phenotype is missing or duplicated after extraction.

The task must localize each `.tbi` file next to its BGZF file with the filename
that tabix expects.

The region list must use a `.bed` suffix so tabix applies zero-based,
half-open BED coordinates. The code must not convert these coordinates through
a GTF.

## Full-scan fallback

The existing full-scan reader remains available for compatibility. It must
produce the same selected rows as the indexed reader. The fallback must log
that it is reading the complete phenotype file so that slow runs are clear in
the task log.

No workflow can silently use the fallback when indexed inputs were supplied but
are invalid.

## Sample alignment and output

After row extraction, both access paths use the same downstream code. That code
must:

- normalize sample IDs;
- compute the intersection across nonempty modalities;
- retain the order from the first contributing modality;
- report input, shared, and removed sample counts;
- normalize the first four output metadata names;
- write the current combined phenotype BED, manifest, and QC files.

The phenotype QC table must add the access method and extraction timing for each
modality. Existing QC columns must remain unchanged.

## Logging

The preparation task must log these values for each modality:

- access method: `tabix` or `full_scan`;
- requested phenotype count;
- lookup time;
- tabix query time, when applicable;
- parse time;
- returned row count before exact-ID filtering;
- retained row count after exact-ID filtering;
- input, shared, and removed sample counts;
- total modality time.

The index workflow and the preparation WDL commands must include timestamped
logging messages.

## Containers

The phenotype indexing and preparation images must include `bgzip` and `tabix`.
The image definitions must be tested in GitHub Actions. Local smoke tests must
not build the images.

The prepare-window phenotype image test must run both the indexed and full-scan
paths inside the image.

## Tests

The automated tests must cover:

- one-time index creation from a fixture phenotype file;
- successful tabix queries for scattered phenotype IDs;
- overlapping query intervals without duplicate output outcomes;
- different input metadata header names;
- different sample sets and sample orders across modalities;
- exact output equivalence between indexed and full-scan paths;
- preservation of requested phenotype order;
- one modality on the indexed path and another on the fallback path;
- a missing target phenotype;
- a missing selected trans phenotype;
- duplicate phenotype IDs in the lookup;
- an incomplete index and lookup input pair;
- no shared samples;
- WDL validation;
- timestamped WDL logging;
- GitHub Actions image build and in-container smoke tests.

The equivalence test must compare the combined phenotype data, manifest, and
all pre-existing QC values. Timing and access-method fields can differ.

## Expected performance

The one-time index workflow reads each full phenotype file. Each window job then
reads only the header, the coordinate lookup table, and the selected phenotype
rows. For the default joint analysis, this is approximately 15 to 25 rows per
modality instead of every row in each genome-wide phenotype file.

Performance depends on storage latency and BGZF block layout. The logs provide
the measurements needed to compare the indexed and fallback paths on the real
data.
