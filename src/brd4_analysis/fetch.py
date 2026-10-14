"""Download and fingerprint the primary experimental data.

Retrieves PDB coordinate files, Chemical Component Dictionary (CCD)
ideal-coordinate SDFs for the bound ligands, and RCSB entry metadata.
Every byte stored under ``data/`` is fingerprinted (SHA-256) and recorded
in ``data/provenance/provenance.json`` together with its source URL and
UTC retrieval timestamp.
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


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(url, timeout=HTTP_TIMEOUT_S, headers={"User-Agent": USER_AGENT})
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return {
        "path": str(dest),
        "url": url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": sha256_file(dest),
        "bytes": dest.stat().st_size,
    }


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


def fetch_pdb(pdb_id: str, data_dir: Path) -> dict:
    record = _download(PDB_DOWNLOAD_URL.format(pdb_id=pdb_id), data_dir / "raw" / f"{pdb_id}.pdb")
    record.update(parse_pdb_dates(data_dir / "raw" / f"{pdb_id}.pdb"))
    return record


def fetch_ccd_sdf(comp_id: str, data_dir: Path) -> dict:
    return _download(
        CCD_IDEAL_SDF_URL.format(comp_id=comp_id), data_dir / "ccd" / f"{comp_id}_ideal.sdf"
    )


def fetch_entry_metadata(pdb_id: str, data_dir: Path) -> dict:
    dest = data_dir / "raw" / f"{pdb_id}_entry.json"
    record = _download(ENTRY_METADATA_URL.format(pdb_id=pdb_id), dest)
    try:
        meta = json.loads(dest.read_text())
        record["title"] = meta.get("struct", {}).get("title")
        record["deposition_date"] = (
            meta.get("pdbx_database_status", {}).get("recvd_initial_deposition_date")
        )
        record["experimental_method"] = (meta.get("exptl") or [{}])[0].get("method")
        record["resolution_angstrom"] = (meta.get("rcsb_entry_info") or {}).get(
            "resolution_combined"
        )
    except (json.JSONDecodeError, KeyError):
        record["metadata_parse_error"] = True
    return record


def fetch_all(data_dir: Path = DEFAULT_DATA_DIR) -> Path:
    """Fetch all primary data and write provenance.json."""
    data_dir = Path(data_dir)
    records: dict[str, dict] = {}
    for pdb_id in PDB_IDS:
        records[pdb_id] = {"pdb_file": fetch_pdb(pdb_id, data_dir)}
        records[pdb_id]["entry_metadata"] = fetch_entry_metadata(pdb_id, data_dir)
        comp_id = LIGAND_COMPONENTS[pdb_id]
        records[pdb_id]["ligand_ccd_ideal_sdf"] = fetch_ccd_sdf(comp_id, data_dir)

    provenance = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": {
            "package": "brd4-analysis",
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "entries": records,
    }
    out = data_dir / "provenance" / "provenance.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(provenance, indent=2) + "\n")
    return out
