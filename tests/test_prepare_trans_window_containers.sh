#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

test -s envs/prepare-window-genotypes.Dockerfile
test -s envs/prepare-window-phenotypes.Dockerfile
test -s .github/workflows/prepare-window-genotypes-image.yml
test -s .github/workflows/prepare-window-phenotypes-image.yml

rg -q '^FROM ubuntu:24[.]04$' envs/prepare-window-genotypes.Dockerfile
rg -q 'apt-get.*install|apt-get install' envs/prepare-window-genotypes.Dockerfile
rg -q 'tabix' envs/prepare-window-genotypes.Dockerfile

rg -q '^FROM rocker/r-ver:' envs/prepare-window-phenotypes.Dockerfile
rg -q 'optparse' envs/prepare-window-phenotypes.Dockerfile
for package in dplyr purrr readr tibble; do
  rg -q "^[[:space:]]+${package}([[:space:]]|$)" envs/prepare-window-phenotypes.Dockerfile
done
rg -q 'COPY scripts/prepare_trans_window[.]R scripts/trans_window_cli[.]R /opt/mvsusie/scripts/' envs/prepare-window-phenotypes.Dockerfile

rg -q 'name: prepare-trans-window' .dockstore.yml
rg -q 'primaryDescriptorPath: /workflows/prepare_trans_window[.]wdl' .dockstore.yml

rg -q 'prepare-window-genotypes' workflows/prepare_trans_window.wdl
rg -q 'prepare-window-phenotypes' workflows/prepare_trans_window.wdl

for workflow in \
  .github/workflows/prepare-window-genotypes-image.yml \
  .github/workflows/prepare-window-phenotypes-image.yml; do
  rg -q 'docker/build-push-action@v7' "$workflow"
  rg -q 'docker/metadata-action@v6' "$workflow"
  rg -q 'packages: write' "$workflow"
  rg -q 'push: true' "$workflow"
done

rg -q 'envs/prepare-window-genotypes[.]Dockerfile' .github/workflows/prepare-window-genotypes-image.yml
rg -q 'envs/prepare-window-phenotypes[.]Dockerfile' .github/workflows/prepare-window-phenotypes-image.yml

if rg -q 'workflows/prepare_trans_window[.]wdl|prepare-window-(genotypes|phenotypes)-image[.]yml' \
  .github/workflows/prepare-window-genotypes-image.yml \
  .github/workflows/prepare-window-phenotypes-image.yml; then
  echo "Container rebuild triggers must not include workflow files." >&2
  exit 1
fi

echo "Preparation container definitions passed"
