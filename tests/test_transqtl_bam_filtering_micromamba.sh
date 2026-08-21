#!/usr/bin/env bash
set -euo pipefail

dockerfile="envs/transqtl-bam-filtering.Dockerfile"
environment="envs/transqtl-bam-filtering.environment.yml"

test -f "$environment"
grep -Fx '  - conda-forge' "$environment"
grep -Fx '  - bioconda' "$environment"
grep -Eq '^  - samtools(=|$)' "$environment"
grep -Fx '  - rna-seqc=2.4.2' "$environment"

grep -Eq '^FROM mambaorg/micromamba:' "$dockerfile"
grep -Fx 'ENV MAMBA_DOCKERFILE_ACTIVATE=1' "$dockerfile"
grep -F 'micromamba install' "$dockerfile"
grep -F 'micromamba config set channel_priority strict' "$dockerfile"
grep -Fx 'ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin' "$dockerfile"
! grep -F 'apt-get' "$dockerfile"
! grep -F 'https://github.com/getzlab/rnaseqc/releases' "$dockerfile"
