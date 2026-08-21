#!/usr/bin/env bash
set -euo pipefail

dockerfile="envs/transqtl-mappability.Dockerfile"
environment="envs/transqtl-mappability.environment.yml"

test -f "$environment"
grep -Fx '  - conda-forge' "$environment"
grep -Fx '  - bioconda' "$environment"
grep -Eq '^  - genmap(=|$)' "$environment"
grep -Eq '^  - bedtools(=|$)' "$environment"

grep -Eq '^FROM mambaorg/micromamba:' "$dockerfile"
grep -Fx 'ENV MAMBA_DOCKERFILE_ACTIVATE=1' "$dockerfile"
grep -F 'micromamba install' "$dockerfile"
grep -F 'micromamba config set channel_priority strict' "$dockerfile"
grep -Fx 'ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin' "$dockerfile"
! grep -F 'apt-get' "$dockerfile"

echo "TransQTLMappability micromamba validation passed"
