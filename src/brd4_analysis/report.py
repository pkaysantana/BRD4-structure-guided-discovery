"""Inject computed values into docs/REPORT.md.

The narrative in REPORT.md is hand-written; every number lives inside
generated blocks delimited by

    <!-- BEGIN GENERATED: name -->
    ...
    <!-- END GENERATED: name -->

which are re-rendered from ``results/summary.json`` and
``data/provenance/provenance.json``. A test asserts the committed report
equals a fresh render, so reported values cannot drift from the outputs.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Callable

from brd4_analysis.constants import (
    CONTACT_THRESHOLDS,
    DEFAULT_DATA_DIR,
    DEFAULT_RESULTS_DIR,
    PDB_IDS,
    PROJECT_ROOT,
)
from brd4_analysis.fetch import load_provenance

DEFAULT_REPORT = PROJECT_ROOT / "docs" / "REPORT.md"
BLOCK_RE = re.compile(
    r"(<!-- BEGIN GENERATED: (?P<name>[\w-]+) -->\n)(?P<body>.*?)(<!-- END GENERATED: (?P=name) -->)",
    re.DOTALL,
)


def _table(header: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines) + "\n"


def _provenance(summary: dict, prov: dict) -> str:
    rows = []
    for pid in PDB_IDS:
        e = prov["entries"][pid]
        meta, pdb = e["entry_metadata"], e["pdb_file"]
        res = ", ".join(f"{r:.2f}" for r in meta.get("resolution_angstrom") or [])
        rows.append([
            pid, summary[pid]["ligand_component"], meta.get("experimental_method"), res,
            (meta.get("deposition_date") or "")[:10],
            f"#{pdb['latest_revdat_num']}, {pdb['latest_revdat_date']}",
            pdb["retrieved_utc"][:10], f"`{pdb['sha256']}`",
        ])
    return _table(
        ["Entry", "Ligand (CCD)", "Method", "Resolution (Å)", "Deposited",
         "Latest REVDAT", "Retrieved", "PDB file SHA-256"],
        rows,
    ) + f"\nInput files recorded in the manifest: {sum(len(e) for e in prov['entries'].values())}.\n"


def _descriptors(summary: dict, _prov: dict) -> str:
    spec: list[tuple[str, str, Callable[[Any], str]]] = [
        ("Formula", "ccd_formula", str),
        ("Isomeric SMILES", "isomeric_smiles", lambda v: f"`{v}`"),
        ("InChIKey", "inchikey", str),
        ("MW (Da)", "molecular_weight", lambda v: f"{v:.2f}"),
        ("cLogP (Crippen)", "clogp", lambda v: f"{v:.2f}"),
        ("HBD (RDKit CalcNumHBD)", "hbd", str),
        ("HBA (RDKit CalcNumHBA)", "hba", str),
        ("HBD, Lipinski N-H + O-H", "hbd_lipinski_nhoh", str),
        ("HBA, Lipinski N + O", "hba_lipinski_no", str),
        ("TPSA (Å²)", "tpsa", lambda v: f"{v:.2f}"),
        ("Rotatable bonds", "rotatable_bonds", str),
        ("Heavy atoms", "heavy_atom_count", str),
        ("Aromatic rings", "aromatic_ring_count", str),
        ("Fraction Csp3", "fraction_csp3", lambda v: f"{v:.2f}"),
        ("Formal charge", "formal_charge", str),
    ]
    header = ["Descriptor"] + [
        f"{summary[p]['ligand_component']} ({p})" for p in PDB_IDS
    ]
    rows = [[label] + [fmt(summary[p]["descriptors"][key]) for p in PDB_IDS] for label, key, fmt in spec]
    return _table(header, rows)


def _similarity(summary: dict, _prov: dict) -> str:
    s = summary["ligand_similarity"]
    b = s["bit_counts"]
    return (
        f"Tanimoto (4MEN vs 4MEO) = **{s['morgan_ecfp4_tanimoto_4MEN_vs_4MEO']:.4f}** "
        f"({b['common_bits']} common / {b['union_bits']} union on-bits; "
        f"{b['on_bits_a']} and {b['on_bits_b']} on-bits per ligand). "
        f"Fingerprint: {s['fingerprint']}.\n"
    )


def _thresholds(summary: dict, _prov: dict) -> str:
    rows = [
        [f"≤ {t} Å"] + [summary[p]["residues_per_threshold"][str(t)] for p in PDB_IDS]
        for t in CONTACT_THRESHOLDS
    ]
    return _table(["Threshold"] + [f"{p} residues" for p in PDB_IDS], rows)


def _closest(summary: dict, _prov: dict) -> str:
    rows = [
        [p, c["protein_atom"], c["ligand_atom"], f"{c['distance_A']:.3f}"]
        for p in PDB_IDS
        for c in summary[p]["closest_atom_pairs"]
    ]
    return _table(["Entry", "Protein atom", "Ligand atom", "Distance (Å)"], rows)


def _model_quality(summary: dict, _prov: dict) -> str:
    rows = []
    for p in PDB_IDS:
        q, n = summary[p]["ligand_model_quality"], summary[p]["nonprotein_neighbours"]
        het = ", ".join(f"{k} {v:.2f}" for k, v in n["het_group_min_distance_A"].items())
        rows.append([
            p, ", ".join(f"{o:.2f}" for o in q["occupancies"]),
            f"{q['b_factor_min']:.1f} / {q['b_factor_mean']:.1f} / {q['b_factor_max']:.1f}",
            n["waters_within_3.5A"], het,
        ])
    return _table(
        ["Entry", "Ligand occupancy", "Ligand B min / mean / max (Å²)",
         "Waters ≤ 3.5 Å of ligand", "Other HET groups, min distance (Å)"],
        rows,
    )


def _vdw(summary: dict, _prov: dict) -> str:
    rows = []
    for p in PDB_IDS:
        scr = summary[p]["vdw_overlap_screen"]
        rows += [
            [p, r["protein_atom"], r["ligand_atom"], f"{r['distance_A']:.3f}",
             f"{r['bondi_sum_A']:.2f}", f"{r['overlap_A']:.3f}"]
            for r in scr["pairs"]
        ]
    thr = summary[PDB_IDS[0]]["vdw_overlap_screen"]["min_overlap_A"]
    note = f"Screen: carbon-containing heavy-atom pairs with Bondi overlap ≥ {thr} Å."
    empty = [p for p in PDB_IDS if not summary[p]["vdw_overlap_screen"]["pairs"]]
    if empty:
        note += f" No pairs flagged in: {', '.join(empty)}."
    if not rows:
        return note + "\n"
    return _table(
        ["Entry", "Protein atom", "Ligand atom", "Distance (Å)", "Bondi sum (Å)", "Overlap (Å)"],
        rows,
    ) + "\n" + note + "\n"


RENDERERS: dict[str, Callable[[dict, dict], str]] = {
    "provenance": _provenance,
    "descriptors": _descriptors,
    "similarity": _similarity,
    "thresholds": _thresholds,
    "closest-pairs": _closest,
    "model-quality": _model_quality,
    "vdw-overlaps": _vdw,
}


def render_report(text: str, summary: dict, prov: dict) -> str:
    text = text.replace("\r\n", "\n")
    seen = set()

    def repl(m: re.Match) -> str:
        name = m.group("name")
        if name not in RENDERERS:
            raise KeyError(f"Unknown generated block: {name}")
        seen.add(name)
        return m.group(1) + RENDERERS[name](summary, prov) + m.group(4)

    out = BLOCK_RE.sub(repl, text)
    missing = set(RENDERERS) - seen
    if missing:
        raise KeyError(f"Report is missing generated blocks: {sorted(missing)}")
    return out


def update_report(
    report_path: Path = DEFAULT_REPORT,
    data_dir: Path = DEFAULT_DATA_DIR,
    results_dir: Path = DEFAULT_RESULTS_DIR,
) -> Path:
    report_path = Path(report_path)
    summary = json.loads((Path(results_dir) / "summary.json").read_text(encoding="utf-8"))
    rendered = render_report(
        report_path.read_text(encoding="utf-8"), summary, load_provenance(data_dir)
    )
    report_path.write_bytes(rendered.encode("utf-8"))
    return report_path
