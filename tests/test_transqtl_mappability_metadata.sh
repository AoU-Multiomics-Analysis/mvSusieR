#!/usr/bin/env bash
set -euo pipefail

ruby -ryaml -e '
  config = YAML.load_file(".dockstore.yml")
  workflow = config.fetch("workflows").find { |item| item["name"] == "transqtl-mappability" }
  abort "missing workflow" unless workflow
  abort "unexpected WDL path" unless workflow["primaryDescriptorPath"] == "/workflows/transqtl_mappability.wdl"
'

echo "TransQTLMappability Dockstore metadata validation passed"
