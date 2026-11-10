#!/bin/bash
# Double-click launcher for the Tkinter GUI (macOS).
cd "$(dirname "$0")"
if [ ! -x .venv/bin/python ]; then
  echo "Missing .venv - run ./setup.sh first."
  read -r -p "Press Enter to close..."
  exit 1
fi
exec ./.venv/bin/python alignment_conservation_plot.py --gui
