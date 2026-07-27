#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${repo_dir}/.cache}"

exec snakemake \
  --directory "${repo_dir}" \
  --cores "${SNAKEMAKE_CORES:-32}" \
  --use-conda \
  --conda-prefix "${repo_dir}/.conda/rules" \
  --printshellcmds \
  "$@"
