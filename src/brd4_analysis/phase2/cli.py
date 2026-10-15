"""Offline Phase 2 CLI. There is no execution or download subcommand."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from brd4_analysis.constants import PROJECT_ROOT
from brd4_analysis.phase2.audit import write_audit
from brd4_analysis.phase2.runtime import capture_runtime
from brd4_analysis.phase2.workflow import analyse_run, prepare, validate

DEFAULT_CONFIG = PROJECT_ROOT / "configs/phase2/4men_pilot.json"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ["validate", "prepare"]:
        child = sub.add_parser(name)
        child.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    audit = sub.add_parser("audit", help="reproduce selected research values and the 32SC duplicate check")
    audit.add_argument("--output", type=Path, default=PROJECT_ROOT / "data/phase2/review")
    for name in ["capture-runtime", "analyse"]:
        child = sub.add_parser(name)
        child.add_argument("run", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            _, _, report = validate(args.config)
            print(json.dumps(report, indent=2, allow_nan=False))
        elif args.command == "prepare":
            print(prepare(args.config))
        elif args.command == "audit":
            print(write_audit(args.output))
        elif args.command == "capture-runtime":
            print(capture_runtime(args.run))
        elif args.command == "analyse":
            print(analyse_run(args.run))
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f"{args.command} stopped: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
