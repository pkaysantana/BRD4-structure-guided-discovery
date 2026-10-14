"""Download and fingerprint the primary experimental data.

Retrieves PDB coordinate files, Chemical Component Dictionary (CCD)
ideal-coordinate SDFs for the bound ligands, and RCSB entry metadata.
Every file stored under ``data/`` is fingerprinted (SHA-256) and recorded
in ``data/provenance/provenance.json`` together with its source URL and
UTC retrieval timestamp.

Original files are never silently replaced: if a file already exists and
its checksum matches the manifest, it is kept (and its original retrieval
record preserved). A checksum mismatch is an error unless ``force=True``.
Paths in the manifest are POSIX paths relative to the data directory, so
the manifest is portable across machines.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

from brd4_analysis.constants import (
    CCD_IDEAL_SDF_URL,
    DEFAULT_DATA_DIR,
    ENTRY_METADATA_URL,
    LIGAND_COMPONENTS,
    PDB_DOWNLOAD_URL,
    PDB_IDS,
)

HTTP_TIMEOUT_S = 30
USER_AGENT = "brd4-structure-guided-discovery/0.1.0 (reproducible research baseline)"
PROVENANCE_RELPATH = Path("provenance") / "provenance.json"
SCHEMA_VERSION = 2


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_pdb_dates(path: Path) -> dict:
    """Extract deposition date and latest REVDAT record from a PDB file."""
    header_dep = None
    revdat_entries: list[tuple[int, str]] = []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("HEADER"):
            # cols 51-59: deposition date DD-MMM-YY
            header_dep = line[50:59].strip() or None
        elif line.startswith("REVDAT"):
            try:
                modnum = int(line[7:10])
            except ValueError:
                continue
            revdat_entries.append((modnum, line[13:22].strip()))
    latest = max(revdat_entries) if revdat_entries else None
    return {
        "header_deposition_date": header_dep,
        "revdat_count": len(revdat_entries),
        "latest_revdat_num": latest[0] if latest else None,
        "latest_revdat_date": latest[1] if latest else None,
    }


def parse_entry_metadata(path: Path) -> dict:
    try:
        meta = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {"metadata_parse_error": True}
    return {
        "title": meta.get("struct", {}).get("title"),
        "deposition_date": meta.get("pdbx_database_status", {}).get(
            "recvd_initial_deposition_date"
        ),
        "experimental_method": (meta.get("exptl") or [{}])[0].get("method"),
        "resolution_angstrom": (meta.get("rcsb_entry_info") or {}).get("resolution_combined"),
    }


def input_targets() -> list[tuple[str, str, str, Path]]:
    """(pdb_id, record_kind, url, path relative to data dir) for every input."""
    targets = []
    for pdb_id in PDB_IDS:
        comp_id = LIGAND_COMPONENTS[pdb_id]
        targets += [
            (pdb_id, "pdb_file", PDB_DOWNLOAD_URL.format(pdb_id=pdb_id), Path("raw") / f"{pdb_id}.pdb"),
            (pdb_id, "entry_metadata", ENTRY_METADATA_URL.format(pdb_id=pdb_id), Path("raw") / f"{pdb_id}_entry.json"),
            (pdb_id, "ligand_ccd_ideal_sdf", CCD_IDEAL_SDF_URL.format(comp_id=comp_id), Path("ccd") / f"{comp_id}_ideal.sdf"),
        ]
    return targets


def _derived_fields(kind: str, path: Path) -> dict:
    if kind == "pdb_file":
        return parse_pdb_dates(path)
    if kind == "entry_metadata":
        return parse_entry_metadata(path)
    return {}


def _download(url: str, dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(url, timeout=HTTP_TIMEOUT_S, headers={"User-Agent": USER_AGENT})
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_provenance(data_dir: Path) -> dict:
    path = Path(data_dir) / PROVENANCE_RELPATH
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def fetch_all(data_dir: Path = DEFAULT_DATA_DIR, force: bool = False) -> Path:
    """Fetch missing inputs, verify existing ones, and write provenance.json.

    Existing files whose SHA-256 matches the manifest are kept together with
    their original retrieval record. ``force=True`` re-downloads everything.
    The manifest is only rewritten if its content changes.
    """
    data_dir = Path(data_dir)
    previous = load_provenance(data_dir)
    prev_entries = previous.get("entries", {})
    entries: dict[str, dict] = {}
    downloaded = False

    for pdb_id, kind, url, rel in input_targets():
        dest = data_dir / rel
        prev = prev_entries.get(pdb_id, {}).get(kind)
        if not force and dest.exists() and prev is not None:
            if sha256_file(dest) != prev["sha256"]:
                raise RuntimeError(
                    f"{rel.as_posix()} does not match its recorded SHA-256; "
                    "refusing to overwrite. Re-run with --force to re-download."
                )
            retrieved = prev["retrieved_utc"]
        else:
            retrieved = _download(url, dest)
            downloaded = True
        entries.setdefault(pdb_id, {})[kind] = {
            "path": rel.as_posix(),
            "url": url,
            "retrieved_utc": retrieved,
            "sha256": sha256_file(dest),
            "bytes": dest.stat().st_size,
            **_derived_fields(kind, dest),
        }

    if downloaded:
        environment = {
            "package": "brd4-analysis",
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        }
    else:
        # schema v1 called this block "generator"
        environment = previous.get("retrieval_environment") or previous.get("generator")
    provenance = {
        "schema_version": SCHEMA_VERSION,
        "paths_relative_to": "data/",
        "retrieval_environment": environment,
        "entries": entries,
    }
    out = data_dir / PROVENANCE_RELPATH
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(provenance, indent=2) + "\n").encode("utf-8")
    # write_bytes: deterministic LF bytes regardless of platform newline rules
    if not out.exists() or out.read_bytes() != payload:
        out.write_bytes(payload)
    return out
