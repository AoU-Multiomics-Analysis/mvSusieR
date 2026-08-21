#!/usr/bin/env bash
set -euo pipefail

docker build \
  --tag transqtl-bam-filtering-test \
  --file envs/transqtl-bam-filtering.Dockerfile \
  .

docker run --rm transqtl-bam-filtering-test samtools --version
docker run --rm transqtl-bam-filtering-test rnaseqc --version
