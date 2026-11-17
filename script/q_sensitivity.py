"""Compatibility entry point; use python -m dcs for the full CLI."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dcs.cli import main

if __name__ == "__main__":
    main(["scan", "q"] + sys.argv[1:])
