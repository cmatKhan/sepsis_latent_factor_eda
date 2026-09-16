#!/usr/bin/env bash
# Re-runs R/create_slurm_bundle.R (every method present in each dataset's
# `methods:` block, no --only filter) for every dataset config in
# config/, skipping the config/dataset_metadata.example.yml template.
# Continues past a failing dataset (e.g. one still missing its
# preprocessing_script) rather than aborting the whole run, and prints a
# pass/fail summary at the end.
#
# Usage:
#   ./run_setup_all.sh
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

ok=()
failed=()

for cfg in config/*_config.yml; do
  [ "$(basename "$cfg")" = "dataset_metadata.example.yml" ] && continue

  echo "=================================================================="
  echo "== $cfg"
  echo "=================================================================="
  if Rscript R/create_slurm_bundle.R "$cfg"; then
    ok+=("$cfg")
  else
    echo "!! setup failed for $cfg -- continuing with the next dataset" >&2
    failed+=("$cfg")
  fi
  echo
done

echo "=================================================================="
echo "done: ${#ok[@]} succeeded, ${#failed[@]} failed"
[ ${#ok[@]} -gt 0 ] && printf '  ok:     %s\n' "${ok[@]}"
[ ${#failed[@]} -gt 0 ] && printf '  failed: %s\n' "${failed[@]}"

# nonzero exit if anything failed, so this is CI/script-friendly too
[ ${#failed[@]} -eq 0 ]
