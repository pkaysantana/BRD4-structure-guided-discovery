"""End-to-end Phase 1 analysis pipeline.

Reads the raw files under ``data/`` (populated by ``fetch``) and writes:
  * ``results/tables/ligand_descriptors.csv`` — comparison table
  * ``results/tables/contacts_summary_<PDB>.csv``
  * ``results/tables/contact_pairs_<PDB>.csv`` (all pairs <= 5.0 Å)
  * ``results/tables/vdw_overlaps_<PDB>.csv`` (geometric screen)
  * ``results/tables/threshold_counts.csv``
  * ``results/figures/*.png``
  * ``results/summary.json`` — every headline number used in the report.
    Deterministic: depends only on the input files and pinned software.
  * ``results/run_metadata.json`` — variable run information (timestamp,
    platform); kept separate so summary.json is byte-reproducible.
"""

from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import RDLogger

from brd4_analysis.constants import (
    CONTACT_THRESHOLDS,
    DEFAULT_DATA_DIR,
    DEFAULT_RESULTS_DIR,
    LIGAND_COMPONENTS,
    PDB_IDS,
)
from brd4_analysis.contacts import (
    VDW_OVERLAP_SCREEN_A,
    contact_atom_pairs,
    contacts_by_residue,
    distance_matrix,
    threshold_counts,
    vdw_overlaps,
)
from brd4_analysis.descriptors import (
    compute_descriptors,
    load_ccd_mol,
    tanimoto_morgan,
    verify_pdb_ligand_identity,
)
from brd4_analysis.fetch import input_targets, sha256_file
from brd4_analysis.parse import (
    ParsedEntry,
    extract_het_record_block,
    is_heavy_atom,
    ligand_element_counts,
    load_entry,
)
from brd4_analysis.visualize import plot_residue_contacts, save_ligand_2d

RDLogger.DisableLog("rdApp.warning")

DESCRIPTOR_COLUMNS = [
    "component_id", "ccd_formula", "canonical_smiles", "isomeric_smiles",
    "inchikey", "molecular_weight", "exact_mw", "clogp", "hbd", "hba",
    "hbd_lipinski_nhoh", "hba_lipinski_no",
    "tpsa", "rotatable_bonds", "heavy_atom_count", "ring_count",
    "aromatic_ring_count", "fraction_csp3", "formal_charge",
]
SOFTWARE_PACKAGES = ("rdkit", "biopython", "numpy", "pandas", "matplotlib")
WATER_SHELL_A = 3.5


def ligand_model_quality(entry: ParsedEntry) -> dict[str, Any]:
    """Occupancy and B-factor summary of the deposited ligand atoms."""
    b = np.array([a.bfactor for a in entry.ligand_atoms], dtype=float)
    return {
        "occupancies": sorted({float(a.occupancy) for a in entry.ligand_atoms}),
        "b_factor_min": round(float(b.min()), 2),
        "b_factor_mean": round(float(b.mean()), 2),
        "b_factor_max": round(float(b.max()), 2),
        "b_factor_by_atom": {
            a.get_id(): round(float(a.bfactor), 2) for a in entry.ligand_atoms
        },
    }


def nonprotein_neighbours(entry: ParsedEntry) -> dict[str, Any]:
    """Ligand proximity to water and non-ligand HET groups (excluded from
    the protein contact tables, reported here so the exclusion is visible)."""
    _, lig = entry.atom_coord_arrays()
    het = {}
    for r in entry.solvent_residues + entry.other_het_residues:
        coords = np.array([a.coord for a in r if is_heavy_atom(a)], dtype=float)
        het[f"{r.get_resname().strip()}{r.id[1]}"] = round(
            float(distance_matrix(coords, lig).min()), 3
        )
    water_o = np.array(
        [a.coord for r in entry.water_residues for a in r if is_heavy_atom(a)], dtype=float
    )
    n_water = int((distance_matrix(water_o, lig).min(axis=1) <= WATER_SHELL_A).sum()) if len(water_o) else 0
    return {
        "het_group_min_distance_A": het,
        f"waters_within_{WATER_SHELL_A}A": n_water,
    }


def analyse_entry(
    pdb_id: str, data_dir: Path, results_dir: Path
) -> tuple[dict[str, Any], Any]:
    comp_id = LIGAND_COMPONENTS[pdb_id]
    pdb_path = data_dir / "raw" / f"{pdb_id}.pdb"
    ccd_path = data_dir / "ccd" / f"{comp_id}_ideal.sdf"

    entry = load_entry(pdb_path)
    ccd_mol = load_ccd_mol(ccd_path)

    identity = verify_pdb_ligand_identity(
        extract_het_record_block(pdb_path, {comp_id}),
        ligand_element_counts(entry),
        ccd_mol,
    )

    contacts = contacts_by_residue(entry, CONTACT_THRESHOLDS)
    pairs = contact_atom_pairs(entry, max(CONTACT_THRESHOLDS))
    overlaps = vdw_overlaps(entry)

    tables = results_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    contacts.to_csv(tables / f"contacts_summary_{pdb_id}.csv", index=False)
    pairs.to_csv(tables / f"contact_pairs_{pdb_id}.csv", index=False)
    overlaps.to_csv(tables / f"vdw_overlaps_{pdb_id}.csv", index=False)

    save_ligand_2d(ccd_mol, comp_id, pdb_id, results_dir / "figures" / f"ligand_{pdb_id}.png")
    plot_residue_contacts(
        contacts, pdb_id, results_dir / "figures" / f"contacts_{pdb_id}.png"
    )

    summary = {
        "pdb_id": pdb_id,
        "ligand_component": comp_id,
        "n_models": entry.n_models,
        "n_polymer_residues": len(entry.polymer_residues),
        "n_ligand_residues": len(entry.ligand_residues),
        "n_ligand_heavy_atoms": len(entry.ligand_atoms),
        "n_water_residues": len(entry.water_residues),
        "solvent_other_het": sorted(
            {r.get_resname().strip() for r in entry.solvent_residues + entry.other_het_residues}
        ),
        "identity_check": identity,
        "ligand_model_quality": ligand_model_quality(entry),
        "descriptors": compute_descriptors(ccd_mol, comp_id),
        "contact_thresholds_A": list(CONTACT_THRESHOLDS),
        "residues_per_threshold": {
            str(t): int(contacts[f"within_{t}A"].sum()) for t in CONTACT_THRESHOLDS
        },
        "atom_pairs_within_max": len(pairs),
        "closest_residues": contacts.nsmallest(5, "min_distance_A")[
            ["residue", "min_distance_A"]
        ].to_dict(orient="records"),
        "closest_atom_pairs": pairs.sort_values(
            ["distance_A", "protein_atom", "ligand_atom"]
        ).head(6)[["protein_atom", "ligand_atom", "distance_A"]].to_dict(orient="records"),
        "vdw_overlap_screen": {
            "min_overlap_A": VDW_OVERLAP_SCREEN_A,
            "pairs": overlaps[
                ["protein_atom", "ligand_atom", "distance_A", "bondi_sum_A", "overlap_A"]
            ].to_dict(orient="records"),
        },
        "nonprotein_neighbours": nonprotein_neighbours(entry),
    }
    return summary, ccd_mol


def run_analysis(
    data_dir: Path = DEFAULT_DATA_DIR, results_dir: Path = DEFAULT_RESULTS_DIR
) -> Path:
    data_dir, results_dir = Path(data_dir), Path(results_dir)
    summaries: dict[str, Any] = {
        "inputs_sha256": {
            rel.as_posix(): sha256_file(data_dir / rel) for _, _, _, rel in input_targets()
        },
        "software": {p: version(p) for p in SOFTWARE_PACKAGES},
    }
    mols = {}
    descriptor_rows = []

    for pdb_id in PDB_IDS:
        summary, mol = analyse_entry(pdb_id, data_dir, results_dir)
        summaries[pdb_id] = summary
        mols[pdb_id] = mol
        descriptor_rows.append(
            {k: v for k, v in summary["descriptors"].items() if k in DESCRIPTOR_COLUMNS}
            | {"pdb_id": pdb_id}
        )

    df = pd.DataFrame(descriptor_rows).set_index("pdb_id")
    df.to_csv(results_dir / "tables" / "ligand_descriptors.csv")

    counts_frames = []
    for pid in PDB_IDS:
        c = threshold_counts(
            pd.read_csv(results_dir / "tables" / f"contacts_summary_{pid}.csv")
        )
        c.insert(0, "pdb_id", pid)
        counts_frames.append(c)
    pd.concat(counts_frames).to_csv(
        results_dir / "tables" / "threshold_counts.csv", index=False
    )

    sim = tanimoto_morgan(mols["4MEN"], mols["4MEO"])
    summaries["ligand_similarity"] = {
        "morgan_ecfp4_tanimoto_4MEN_vs_4MEO": sim["tanimoto"],
        "fingerprint": "Morgan radius 2, 2048 bits, heavy-atom graph (explicit H removed)",
        "bit_counts": {k: v for k, v in sim.items() if k != "tanimoto"},
    }
    summaries["units"] = {"distance": "angstrom", "mass": "dalton", "tpsa": "angstrom^2"}

    out = results_dir / "summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes((json.dumps(summaries, indent=2) + "\n").encode("utf-8"))

    run_meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    (results_dir / "run_metadata.json").write_bytes(
        (json.dumps(run_meta, indent=2) + "\n").encode("utf-8")
    )
    return out
