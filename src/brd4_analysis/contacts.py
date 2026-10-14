"""Heavy-atom distance contacts between the bound ligand and protein.

Purely geometric: a protein residue is "in contact" at threshold t iff the
minimum heavy-atom distance between any of its atoms and any ligand atom
is <= t. No inference of interaction chemistry (hydrogen bonds, hydrophobic
contacts, etc.) is made — that is a Phase 2+ question.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

from brd4_analysis.constants import CONTACT_THRESHOLDS
from brd4_analysis.parse import ParsedEntry, atom_element, is_heavy_atom


def residue_label(residue: Residue) -> str:
    """e.g. 'A:ASN140' — chain, component ID, sequence number."""
    chain_id = residue.get_parent().id
    return f"{chain_id}:{residue.get_resname().strip()}{residue.id[1]}{residue.id[2].strip()}"


def distance_matrix(a_coords: np.ndarray, b_coords: np.ndarray) -> np.ndarray:
    """Pairwise Euclidean distances between two (N,3) coordinate arrays."""
    a = np.asarray(a_coords, dtype=float)
    b = np.asarray(b_coords, dtype=float)
    if a.ndim != 2 or a.shape[1] != 3 or b.ndim != 2 or b.shape[1] != 3:
        raise ValueError("Coordinate arrays must have shape (N, 3)")
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=-1))


def contacts_by_residue(
    entry: ParsedEntry, thresholds: tuple[float, ...] = CONTACT_THRESHOLDS
) -> pd.DataFrame:
    """One row per polymer residue; minimum ligand heavy-atom distance and
    per-threshold contact flags. All distances in Angstrom."""
    rows = []
    lig_coords = np.array([a.coord for a in entry.ligand_atoms], dtype=float)
    for residue in entry.polymer_residues:
        heavy = [a for a in residue if is_heavy_atom(a)]
        coords = np.array([a.coord for a in heavy], dtype=float)
        d = distance_matrix(coords, lig_coords).min() if len(heavy) else np.inf
        row = {
            "pdb_id": entry.pdb_id,
            "residue": residue_label(residue),
            "chain": residue.get_parent().id,
            "resname": residue.get_resname().strip(),
            "resseq": residue.id[1],
            "min_distance_A": round(float(d), 3),
        }
        for t in thresholds:
            row[f"within_{t}A"] = bool(d <= t)
        rows.append(row)
    return pd.DataFrame(rows)


def contact_atom_pairs(entry: ParsedEntry, threshold: float) -> pd.DataFrame:
    """Every ligand-atom/protein-atom pair within ``threshold`` Angstrom."""
    rows = []
    lig_atoms = entry.ligand_atoms
    lig_coords = np.array([a.coord for a in lig_atoms], dtype=float)
    prot_atoms: list[tuple[Residue, Atom]] = [
        (r, a) for r in entry.polymer_residues for a in r if is_heavy_atom(a)
    ]
    prot_coords = np.array([a.coord for _, a in prot_atoms], dtype=float)
    if len(prot_atoms) == 0 or len(lig_atoms) == 0:
        return pd.DataFrame(
            columns=[
                "pdb_id", "residue", "protein_atom", "protein_element",
                "ligand_atom", "ligand_element", "distance_A",
            ]
        )
    dist = distance_matrix(prot_coords, lig_coords)
    for i, j in zip(*np.nonzero(dist <= threshold)):
        res, pa = prot_atoms[i]
        la = lig_atoms[j]
        rows.append(
            {
                "pdb_id": entry.pdb_id,
                "residue": residue_label(res),
                "protein_atom": f"{residue_label(res)}:{pa.name.strip()}",
                "protein_element": atom_element(pa),
                "ligand_atom": f"{entry.ligand_comp_id}:{la.name.strip()}",
                "ligand_element": atom_element(la),
                "distance_A": round(float(dist[i, j]), 3),
            }
        )
    return pd.DataFrame(rows)


def threshold_counts(contacts_df: pd.DataFrame) -> pd.DataFrame:
    """Number of distinct contacting residues at each threshold."""
    rows = []
    for col in contacts_df.columns:
        if col.startswith("within_"):
            t = col[len("within_") : -1]
            rows.append({"threshold_A": float(t), "n_residues": int(contacts_df[col].sum())})
    return pd.DataFrame(rows)
