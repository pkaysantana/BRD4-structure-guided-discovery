"""Command-line interface.

Reproduce the full Phase 1 baseline:

    python -m brd4_analysis.cli all

Individual steps:

    python -m brd4_analysis.cli fetch [--force]  # download missing inputs, verify checksums
    python -m brd4_analysis.cli analyse          # parse, descriptors, contacts, figures
    python -m brd4_analysis.cli report           # re-render numbers in docs/REPORT.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

from brd4_analysis.constants import DEFAULT_DATA_DIR, DEFAULT_RESULTS_DIR
from brd4_analysis.fetch import fetch_all
from brd4_analysis.pipeline import run_analysis
from brd4_analysis.report import DEFAULT_REPORT, update_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="brd4_analysis",
        description="Reproducible experimental baseline for BRD4-BD1 "
        "ligand complexes (PDB 4MEN, 4MEO; Vidler et al. 2013).",
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, help_text in [
        ("fetch", "download missing inputs, verify existing ones against provenance.json"),
        ("all", "fetch, analyse, then report"),
    ]:
        p = sub.add_parser(name, help=help_text)
        p.add_argument(
            "--force", action="store_true",
            help="re-download every input even if checksums match (overwrites data/)",
        )
    sub.add_parser("analyse", help="run the analysis over data/")
    sub.add_parser("report", help="re-render generated number blocks in the report")
    args = parser.parse_args(argv)

    if args.command in {"fetch", "all"}:
        print(f"Wrote {fetch_all(args.data_dir, force=args.force)}")
    if args.command in {"analyse", "all"}:
        print(f"Wrote {run_analysis(args.data_dir, args.results_dir)}")
    if args.command in {"report", "all"}:
        print(f"Wrote {update_report(args.report, args.data_dir, args.results_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
