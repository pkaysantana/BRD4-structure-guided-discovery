"""End-to-end Phase 1 analysis pipeline.

Reads the raw files under ``data/`` (populated by ``fetch``) and writes:
  * ``results/tables/ligand_descriptors.csv`` — comparison table
  * ``results/tables/contacts_summary_<PDB>.csv``
  * ``results/tables/contact_pairs_<PDB>.csv`` (all pairs <= 5.0 Å)
  * ``results/tables/threshold_counts.csv``
  * ``results/figures/*.png``
  * ``results/summary.json`` — every headline number used in the report
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

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
    contact_atom_pairs,
    contacts_by_residue,
    threshold_counts,
)
from brd4_analysis.descriptors import (
    compute_descriptors,
    load_ccd_mol,
    tanimoto_morgan,
    verify_pdb_ligand_identity,
)
from brd4_analysis.parse import (
    extract_het_record_block,
    ligand_element_counts,
    load_entry,
)
from brd4_analysis.visualize import plot_residue_contacts, save_ligand_2d

RDLogger.DisableLog("rdApp.warning")

DESCRIPTOR_COLUMNS = [
    "component_id", "ccd_formula", "canonical_smiles", "isomeric_smiles",
    "inchikey", "molecular_weight", "exact_mw", "clogp", "hbd", "hba",
    "tpsa", "rotatable_bonds", "heavy_atom_count", "ring_count",
    "aromatic_ring_count", "fraction_csp3", "formal_charge",
]


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

    tables = results_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    contacts.to_csv(tables / f"contacts_summary_{pdb_id}.csv", index=False)
    pairs.to_csv(tables / f"contact_pairs_{pdb_id}.csv", index=False)

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
        "descriptors": compute_descriptors(ccd_mol, comp_id),
        "contact_thresholds_A": list(CONTACT_THRESHOLDS),
        "residues_per_threshold": {
            str(t): int(contacts[f"within_{t}A"].sum()) for t in CONTACT_THRESHOLDS
        },
        "atom_pairs_within_max": len(pairs),
        "closest_residues": contacts.nsmallest(5, "min_distance_A")[
            ["residue", "min_distance_A"]
        ].to_dict(orient="records"),
    }
    return summary, ccd_mol


def run_analysis(
    data_dir: Path = DEFAULT_DATA_DIR, results_dir: Path = DEFAULT_RESULTS_DIR
) -> Path:
    data_dir, results_dir = Path(data_dir), Path(results_dir)
    summaries: dict[str, Any] = {}
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

    summaries["ligand_similarity"] = {
        "morgan_ecfp4_tanimoto_4MEN_vs_4MEO": tanimoto_morgan(
            mols["4MEN"], mols["4MEO"]
        )
    }
    summaries["generated_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    summaries["units"] = {"distance": "angstrom", "mass": "dalton", "tpsa": "angstrom^2"}

    out = results_dir / "summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summaries, indent=2) + "\n")
    return out
