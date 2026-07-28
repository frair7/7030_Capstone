#!/usr/bin/env bash
# Run the DMD Explorer test suite from cursor_ai/
set -euo pipefail
cd "$(dirname "$0")/.."
python -m pytest -q "$@"
