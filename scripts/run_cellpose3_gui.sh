#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CELLPOSE_PYTHON="${ROOT_DIR}/.venv-cellpose3/bin/python"

if [[ ! -x "${CELLPOSE_PYTHON}" ]]; then
  echo "missing Cellpose 3 interpreter: ${CELLPOSE_PYTHON}" >&2
  exit 1
fi

exec "${CELLPOSE_PYTHON}" -m cellpose "$@"
