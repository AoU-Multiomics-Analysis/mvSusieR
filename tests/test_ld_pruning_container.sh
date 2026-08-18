#!/usr/bin/env bash
set -euo pipefail

docker build --tag ld-pruning-test --file envs/ld-pruning.Dockerfile .
docker run --rm ld-pruning-test plink2 --version
