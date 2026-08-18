#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

test -s envs/prepare-window-genotypes.Dockerfile
test -s envs/prepare-window-phenotypes.Dockerfile

rg -q '^FROM ubuntu:24[.]04$' envs/prepare-window-genotypes.Dockerfile
rg -q 'apt-get.*install|apt-get install' envs/prepare-window-genotypes.Dockerfile
rg -q 'tabix' envs/prepare-window-genotypes.Dockerfile

rg -q '^FROM rocker/r-ver:' envs/prepare-window-phenotypes.Dockerfile
rg -q 'optparse' envs/prepare-window-phenotypes.Dockerfile
for package in dplyr purrr readr stringr tibble; do
  rg -q "^[[:space:]]+${package}([[:space:]]|$)" envs/prepare-window-phenotypes.Dockerfile
done
rg -q 'COPY scripts/prepare_trans_window[.]R scripts/trans_window_cli[.]R /opt/mvsusie/scripts/' envs/prepare-window-phenotypes.Dockerfile

rg -q 'name: prepare-trans-window' .dockstore.yml
rg -q 'primaryDescriptorPath: /workflows/prepare_trans_window[.]wdl' .dockstore.yml

rg -q 'prepare-window-genotypes' workflows/prepare_trans_window.wdl
rg -q 'prepare-window-phenotypes' workflows/prepare_trans_window.wdl

echo "Preparation container definitions passed"
