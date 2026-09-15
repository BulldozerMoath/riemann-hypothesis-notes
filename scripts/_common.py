"""
Shared argument parsing and reporting helpers for the command-line scripts.

Every script prints the same standing disclaimer and labels its output
``EXPERIMENTAL``, so a reader who runs one in isolation cannot come away
thinking a certified result was produced.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from critical_line.constants import DISCLAIMER_LONG

#: Default working precision for the scripts.
#:
#: The library defaults to 50 digits; the scripts default lower because they do
#: many more evaluations and 30 digits is ample for everything they report.
#: Raise it with ``--dps`` when you want the extra digits.
SCRIPT_DPS = 30

#: Repository root, resolved from this file's location.
REPO_ROOT = Path(__file__).resolve().parent.parent

#: Where figures are written.
OUTPUT_DIR = REPO_ROOT / "outputs"


def base_parser(description: str) -> argparse.ArgumentParser:
    """An argument parser carrying the options every script accepts."""
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--dps", type=int, default=SCRIPT_DPS,
        help="decimal digits of working precision (library default is 50)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=OUTPUT_DIR,
        help="directory for generated figures",
    )
    return parser


def banner(title: str) -> None:
    """Print a script header including the standing disclaimer."""
    print("=" * 74)
    print(title)
    print("=" * 74)
    print(DISCLAIMER_LONG)
    print("=" * 74)
    print()


def section(title: str) -> None:
    """Print a section heading."""
    print()
    print("-" * 74)
    print(title)
    print("-" * 74)


def experimental(message: str) -> None:
    """Print a line explicitly labelled as experimental output."""
    print(f"  [EXPERIMENTAL] {message}")


def saved(path: Path) -> None:
    """Report a written figure."""
    print(f"  [SAVED] {path}")
