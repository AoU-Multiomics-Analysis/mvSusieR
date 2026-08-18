#!/usr/bin/env bash
set -euo pipefail

ruby -ryaml -e '
  config = YAML.load_file(".dockstore.yml")
  workflow = config.fetch("workflows").find { |item| item["name"] == "trans-window-mvsusie" }
  abort "missing workflow" unless workflow
  abort "unexpected WDL path" unless workflow["primaryDescriptorPath"] == "/workflows/trans_window_mvsusie.wdl"
'
