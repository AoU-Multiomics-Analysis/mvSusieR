#!/usr/bin/env bash
set -euo pipefail

docker build \
  --tag transqtl-bam-filtering-test \
  --file envs/transqtl-bam-filtering.Dockerfile \
  .

docker run --rm transqtl-bam-filtering-test samtools --version
docker run --rm transqtl-bam-filtering-test rnaseqc --version
docker run --rm transqtl-bam-filtering-test micromamba --version
docker run --rm transqtl-bam-filtering-test micromamba list --name base

# WDL engines and Docker-to-Apptainer execution can bypass the micromamba
# entrypoint. The base environment must therefore remain on PATH directly.
docker run --rm --entrypoint /bin/bash transqtl-bam-filtering-test \
  -c 'command -v samtools && command -v rnaseqc'
