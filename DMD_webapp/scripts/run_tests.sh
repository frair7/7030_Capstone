#!/usr/bin/env bash
# Run the DMD Explorer test suite from DMD_webapp/
set -euo pipefail
cd "$(dirname "$0")/.."
python -m pytest -q "$@"
