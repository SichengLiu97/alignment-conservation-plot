#!/bin/bash
# Create a virtual environment and install the two dependencies.
#
#   ./setup.sh
#
# Override the interpreter if needed:
#   PYTHON=/usr/bin/python3 ./setup.sh
#
# If your network requires an HTTP proxy:
#   http_proxy=http://127.0.0.1:8080 https_proxy=http://127.0.0.1:8080 ./setup.sh

set -euo pipefail
cd "$(dirname "$0")"

if [ -z "${PYTHON:-}" ]; then
  for candidate in /opt/homebrew/bin/python3.14 /opt/homebrew/bin/python3 /usr/local/bin/python3 /usr/bin/python3; do
    if [ -x "$candidate" ]; then
      PYTHON="$candidate"
      break
    fi
  done
fi

if [ -z "${PYTHON:-}" ]; then
  echo "No python3 found. Set PYTHON=/path/to/python3 and retry." >&2
  exit 1
fi

echo "Using interpreter: $PYTHON ($("$PYTHON" --version 2>&1))"
"$PYTHON" -m venv .venv
./.venv/bin/python -m pip install --upgrade pip >/dev/null
./.venv/bin/python -m pip install -r requirements.txt

cat <<'EOF'

Done.

  CLI : ./.venv/bin/python align_conservation_plot.py -i ALIGNMENT -o figure.svg --pdf
  GUI : ./.venv/bin/python align_conservation_plot.py --gui
        (or double-click run_gui.command)

Matching the interpreter used by this script:
  /opt/homebrew/bin/python3.14  ->  Tk 9.1 (recommended on macOS)
EOF
