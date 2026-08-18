#!/usr/bin/env bash
set -euo pipefail

ruby -ryaml -e '
  config = YAML.load_file(".dockstore.yml")
  workflow = config.fetch("workflows").fetch(0)
  abort "unexpected workflow name" unless workflow["name"] == "ld-pruning"
  abort "unexpected WDL path" unless workflow["primaryDescriptorPath"] == "/workflows/ld_pruning.wdl"
'
