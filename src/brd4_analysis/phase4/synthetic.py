"""Write synthetic Boltz-format mmCIF predictions for adversarial tests.

This reproduces the strict adapter contract of ``phase2.inputs.read_prediction``:
chain A polymer ATOM records keyed by ``label_seq_id`` (Boltz input index),
one chain B HETATM ligand, full occupancy, no alternate conformers. Rows are
written reversed so the adapter cannot depend on file order.

These files are synthetic scoring-harness controls, never model outputs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from brd4_analysis.phase2.inputs import atom_names
from brd4_analysis.phase2.metrics import coordinates, graph_mappings
from brd4_analysis.phase4.perturb import require_proper_rotation, rotation_matrix

ATOM_SITE_FIELDS = [
    "group_PDB", "type_symbol", "label_atom_id", "label_comp_id",
    "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
    "label_alt_id", "occupancy", "pdbx_PDB_model_num",
]

# A fixed, prespecified rigid displacement applied to synthetic structures so
# the writer exercises the harness's protein-alignment path instead of
# handing it already-aligned coordinates.
SYNTHETIC_ROTATION = rotation_matrix([0.3, -1.0, 0.6], 37.0)
SYNTHETIC_TRANSLATION = np.array([11.25, -7.5, 3.75])


def reference_protein_atoms(ref: dict) -> dict:
    """Reference protein coordinates as ``label -> atom name -> xyz``."""
    return {
        label: {name: atom.coord.astype(float) for name, atom in atoms.items()}
        for label, atoms in ref["protein_by_residue"].items()
    }


def prediction_atom_names(path: Path, ref: dict) -> dict:
    """Recover chain-A heavy-atom names in ``read_prediction`` coordinate order.

    The registered reader strips Atom objects, so a faithful serialization of a
    transformed prediction needs the source mmCIF to name residue-level extras
    or reference-ambiguous atoms without reordering them.
    """
    data = MMCIF2Dict(str(path))
    input_by_label = {item["input_index_1"]: item["reference_label"]
                      for item in ref["mapping"]}
    out = {label: [] for label in input_by_label.values()}
    rows = zip(
        data["_atom_site.group_PDB"],
        data["_atom_site.label_asym_id"],
        data["_atom_site.label_seq_id"],
        data["_atom_site.label_atom_id"],
        data["_atom_site.type_symbol"],
    )
    for group, chain, seq, atom_name, element in rows:
        if group != "ATOM" or chain != "A" or element.upper() == "H":
            continue
        label = input_by_label.get(int(seq))
        if label is None:
            raise ValueError(f"Predicted residue A:{seq} is not in the reference map")
        out[label].append(atom_name)
    for label in input_by_label.values():
        if not out[label]:
            raise ValueError(f"No predicted heavy atoms found for {label}")
    return out


def prediction_protein_atoms(path: Path, ref: dict, pred: dict) -> dict:
    """Name a transformed ``read_prediction`` coordinate block faithfully."""
    names_by_residue = prediction_atom_names(path, ref)
    out = {}
    for item in ref["mapping"]:
        label = item["reference_label"]
        xyz = coordinates(pred["protein_by_residue"][label])
        names = names_by_residue[label]
        if len(xyz) != len(names):
            raise ValueError(f"Predicted atom count changed for {label}")
        out[label] = dict(zip(names, xyz))
    return out


def write_prediction(
    path: Path,
    ref: dict,
    ligand_xyz_reference_frame,
    *,
    protein_atoms=None,
    rotation=SYNTHETIC_ROTATION,
    translation=SYNTHETIC_TRANSLATION,
    omit: tuple[int, str] | None = None,
) -> Path:
    """Write a complete synthetic prediction.

    ``ligand_xyz_reference_frame`` and optional ``protein_atoms`` are
    expressed in the reference frame. By default protein heavy atoms are
    replayed from the reference; callers may instead pass a mapping of
    residue label to ``atom name -> coordinates`` to preserve an actual
    predicted protein. Both ligand and protein are then moved by the same
    proper rigid transform.
    """
    rotation = require_proper_rotation(rotation)
    translation = np.asarray(translation, dtype=float)
    if translation.shape != (3,) or not np.isfinite(translation).all():
        raise ValueError("Synthetic translation must be a finite 3-vector")
    ligand = coordinates(ligand_xyz_reference_frame)
    if len(ligand) != ref["mol"].GetNumAtoms():
        raise ValueError("Ligand coordinate block does not match the reference graph")
    protein_atoms = reference_protein_atoms(ref) if protein_atoms is None else protein_atoms
    if set(protein_atoms) != set(ref["protein_by_residue"]):
        raise ValueError("Protein atom mapping must cover every reference residue")
    rows = []
    for item in ref["mapping"]:
        label = item["reference_label"]
        for name, xyz in protein_atoms[label].items():
            if (item["input_index_1"], name) == omit:
                continue
            point = coordinates([xyz])[0] @ rotation + translation
            element = ref["protein_by_residue"][label].get(name)
            element = element.element if element is not None else _atom_element(name)
            rows.append(["ATOM", element, name, item["resname"], "A",
                         item["input_index_1"], *map(float, point), ".", 1, 1])
    mapping = graph_mappings(ref["mol"], ref["predicted_mol"])[0]
    pred_names = atom_names(ref["predicted_mol"])
    for ref_index, pred_index in enumerate(mapping):
        point = ligand[ref_index] @ rotation + translation
        atom = ref["predicted_mol"].GetAtomWithIdx(pred_index)
        rows.append(["HETATM", atom.GetSymbol(), pred_names[pred_index], "LIG1", "B",
                     ".", *map(float, point), ".", 1, 1])
    text = "data_phase4_synthetic_control\nloop_\n"
    text += "".join(f"_atom_site.{field}\n" for field in ATOM_SITE_FIELDS)
    text += "".join(" ".join(map(str, row)) + "\n" for row in reversed(rows)) + "#\n"
    path = Path(path)
    path.write_text(text, encoding="utf-8")
    return path


def _atom_element(name: str) -> str:
    element = name.strip().lstrip("0123456789")[:2]
    if element and element[0].isalpha() and (len(element) == 1 or element[1].islower()):
        return element
    return name.strip().lstrip("0123456789")[:1]
