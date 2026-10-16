"""Reproduce selected research measurements and test the proposed later-release case."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser

from brd4_analysis.constants import PROJECT_ROOT
from brd4_analysis.descriptors import load_ccd_mol, tanimoto_morgan
from brd4_analysis.parse import is_heavy_atom, load_entry
from brd4_analysis.phase2.config import SOURCE_DIR, read_json, sha256, write_json
from brd4_analysis.phase2.metrics import protein_alignment, rmsd
from brd4_analysis.phase2.workflow import verify_sources

SUPPLIED = [
    "docs/research/_claude_science_01_evidence_brief_vidler2013.md",
    "docs/research/02_4MEO_ligand_identity_assessment_claude_science_.md",
    "docs/protocols/03_boltz2_protocol_4MEN_claude_science.md",
    "docs/research/04_interview_questions_claude_science.md",
    "docs/research/05_assessment_and_uncertainties_claude_science_.md",
    "data/derived/_claude_science_measurements_4MEN_4MEO.csv",
]


def later_case_review() -> dict:
    verify_sources()
    parser = MMCIFParser(QUIET=True)
    models = {
        pid: parser.get_structure(pid, str(SOURCE_DIR / f"{pid}.cif"))[0]
        for pid in ["32SC", "9FWX"]
    }
    entries = {pid: read_json(SOURCE_DIR / f"{pid}_entry.json") for pid in models}
    coords = {}
    for pid, model in models.items():
        coords[pid] = {
            (c.id, r.id, a.name): np.asarray(a.coord, dtype=float)
            for c in model for r in c for a in r
            if is_heavy_atom(a) and (r.id[0] == " " or r.resname == "A1IGX")
        }
    if coords["32SC"].keys() != coords["9FWX"].keys():
        raise ValueError("Later-case atom correspondence changed; re-review required")
    ca = sorted(key for key in coords["32SC"] if key[1][0] == " " and key[2] == "CA")
    ligand = sorted(key for key in coords["32SC"] if key[1][0] == "H_A1IGX")
    a = np.array([coords["32SC"][key] for key in ca])
    b = np.array([coords["9FWX"][key] for key in ca])
    rotation, translation, ca_rmsd = protein_alignment(a, b)
    ligand_rmsd = rmsd(
        np.array([coords["32SC"][key] for key in ligand]) @ rotation + translation,
        np.array([coords["9FWX"][key] for key in ligand]),
    )
    exactly_equal = all(np.array_equal(coords["32SC"][key], coords["9FWX"][key]) for key in coords["32SC"])
    similarity = tanimoto_morgan(
        load_ccd_mol(PROJECT_ROOT / "data/ccd/25K_ideal.sdf"),
        load_ccd_mol(SOURCE_DIR / "A1IGX_ideal.sdf"),
    )
    return {
        "entries": {
            pid: {
                "title": entry["struct"]["title"],
                "initial_release_date": entry["rcsb_accession_info"]["initial_release_date"],
                "resolution_A": entry["rcsb_entry_info"]["resolution_combined"],
                "sha256_mmcif": sha256(SOURCE_DIR / f"{pid}.cif"),
            } for pid, entry in entries.items()
        },
        "protein_CA_atoms_compared": len(ca), "ligand_heavy_atoms_compared": len(ligand),
        "all_protein_and_ligand_heavy_atoms_compared": len(coords["32SC"]),
        "all_compared_coordinates_exactly_equal": exactly_equal,
        "protein_CA_rmsd_A": ca_rmsd, "same_named_ligand_rmsd_after_protein_fit_A": ligand_rmsd,
        "ligand_component": "A1IGX", "graph_search": "A1IGX_exact_graph_search.json",
        "25K_vs_A1IGX_morgan": similarity,
        "similarity_scope": "pairwise only, NOT maximum similarity to the pre-cutoff BET corpus",
        "decision": "not_an_independent_2026_case_use_earliest_public_complex_date_2024-08-07",
        "structural_cutoff_comparison": "2024-08-07 is after the report's 2023-06-01 PDB cutoff; checkpoint membership remains unverified",
        "citation_conflict": {
            "32SC_RCSB_and_mmcif_DOI": entries["32SC"]["rcsb_primary_citation"]["pdbx_database_id_DOI"],
            "9FWX_RCSB_DOI": entries["9FWX"]["rcsb_primary_citation"]["pdbx_database_id_DOI"],
            "32SC_PDBe_DOI": read_json(SOURCE_DIR / "32SC_pdbe_publications.json")["32sc"][0]["doi"],
        },
    }


def measurement_review() -> list[dict]:
    provenance = read_json(PROJECT_ROOT / "data/provenance/provenance.json")
    for entry in provenance["entries"].values():
        for item in entry.values():
            if sha256(PROJECT_ROOT / "data" / item["path"]) != item["sha256"]:
                raise ValueError("Phase 1 input no longer matches its provenance")
    entries = {pid: load_entry(PROJECT_ROOT / f"data/raw/{pid}.pdb", pid) for pid in ["4MEN", "4MEO"]}
    ligand = {pid: {a.name: a for a in e.ligand_atoms} for pid, e in entries.items()}

    def distance(pid, atom, residue_name, residue_number, other_atom):
        e = entries[pid]
        matches = [
            r for r in e.polymer_residues + e.water_residues
            if r.resname == residue_name and r.id[1] == residue_number
        ]
        if len(matches) != 1:
            raise ValueError("Ambiguous contact identity")
        return float(np.linalg.norm(ligand[pid][atom].coord.astype(float) - matches[0][other_atom].coord))

    values = {
        "Resolution (A)": {
            pid: read_json(PROJECT_ROOT / f"data/raw/{pid}_entry.json")["rcsb_entry_info"]["resolution_combined"][0]
            for pid in entries
        },
        "Ligand mean B (A^2)": {
            pid: float(np.mean([a.bfactor for a in atoms.values()])) for pid, atoms in ligand.items()
        },
        "Ligand occupancy": {
            pid: sorted({a.occupancy for a in atoms.values()}) for pid, atoms in ligand.items()
        },
        "Deposited formula": {"4MEN": "C15H17N5", "4MEO": "C18H16N2O2"},
        "Key KAc-mimetic H-bond": {
            "4MEN": distance("4MEN", "N2", "ASN", 140, "ND2"),
            "4MEO": distance("4MEO", "O", "ASN", 140, "ND2"),
        },
        "Conserved-water contact": {
            "4MEN": distance("4MEN", "N1", "HOH", 316, "O"),
            "4MEO": distance("4MEO", "O", "HOH", 327, "O"),
        },
        "Additional polar contact": {"4MEN": None, "4MEO": distance("4MEO", "O1", "ASP", 88, "N")},
        "N1-O distance (deposited)": {
            "4MEN": None, "4MEO": float(np.linalg.norm(
                ligand["4MEO"]["N1"].coord.astype(float) - ligand["4MEO"]["O"].coord
            )),
        },
        "B(N1) / B(O) deposited": {
            "4MEN": None, "4MEO": [ligand["4MEO"][name].bfactor for name in ["N1", "O"]],
        },
    }
    from rdkit.Chem import rdMolDescriptors
    for pid, comp in [("4MEN", "25K"), ("4MEO", "25V")]:
        values["Deposited formula"][pid] = rdMolDescriptors.CalcMolFormula(
            load_ccd_mol(PROJECT_ROOT / f"data/ccd/{comp}_ideal.sdf")
        )
    source = PROJECT_ROOT / SUPPLIED[-1]
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    reviewed = []
    for row in rows:
        metric = row["metric"]
        reviewed.append({
            "supplied": row,
            "status": (
                "deposited_metadata_rechecked" if metric == "Resolution (A)"
                else "ccd_graph_calculation_reproduced" if metric == "Deposited formula"
                else "locally_reproduced_from_deposited_model" if metric in values
                else "external_or_author_computation_not_reproduced"
            ),
            "local_values": values.get(metric),
            "interpretation": (
                "Distances are geometry; no H-bond assignment or experimental sample identity established."
                if "contact" in metric.lower() or "bond" in metric.lower() else
                "Model parameters and template chemistry do not establish the chemical species in the crystal."
            ),
        })
    return reviewed


def write_audit(output: Path) -> Path:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    measurements = measurement_review()
    write_json(output / "research_inventory.json", {
        "review_date": "2026-10-15",
        "supplied_files": [{"path": p, "sha256": sha256(PROJECT_ROOT / p)} for p in SUPPLIED],
        "missing_deliverables": [],
        "review_status": "all six supplied files reviewed; original source files preserved without endorsement",
        "downloads_intake": [
            {
                "original_filename": "03_boltz2_protocol_4MEN.md",
                "original_folder": "user Downloads",
                "archived_path": "docs/protocols/03_boltz2_protocol_4MEN_claude_science.md",
            },
            {
                "original_filename": "04_interview_questions.md",
                "original_folder": "user Downloads",
                "archived_path": "docs/research/04_interview_questions_claude_science.md",
            },
        ],
    })
    write_json(output / "measurements_reconciliation.json", measurements)
    write_json(output / "32SC_assessment.json", later_case_review())
    return output
