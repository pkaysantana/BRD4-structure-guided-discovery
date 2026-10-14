"""Command-line interface.

Reproduce the full Phase 1 baseline from scratch:

    python -m brd4_analysis.cli all

Individual steps:

    python -m brd4_analysis.cli fetch     # download + fingerprint raw data
    python -m brd4_analysis.cli analyse   # parse, descriptors, contacts, figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

from brd4_analysis.constants import DEFAULT_DATA_DIR, DEFAULT_RESULTS_DIR
from brd4_analysis.fetch import fetch_all
from brd4_analysis.pipeline import run_analysis


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="brd4_analysis",
        description="Reproducible experimental baseline for BRD4-BD1 "
        "ligand complexes (PDB 4MEN, 4MEO; Vidler et al. 2013).",
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("fetch", help="download raw data and write provenance.json")
    sub.add_parser("analyse", help="run the analysis over data/raw")
    sub.add_parser("all", help="fetch, then analyse")
    args = parser.parse_args(argv)

    if args.command in {"fetch", "all"}:
        prov = fetch_all(args.data_dir)
        print(f"Wrote {prov}")
    if args.command in {"analyse", "all"}:
        out = run_analysis(args.data_dir, args.results_dir)
        print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
