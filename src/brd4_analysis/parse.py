"""PDB parsing and entity classification via Biopython.

Separates each entry into:
  * polymer protein atoms (standard amino-acid residues),
  * the bound ligand (expected component ID per entry, see constants),
  * solvent / crystallisation additives and water (reported, not analysed
    as ligand contacts).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue

from brd4_analysis.constants import LIGAND_COMPONENTS, SOLVENT_COMPONENT_IDS


def atom_element(atom: Atom) -> str:
    """Return the element symbol, falling back to the atom name if needed."""
    el = (atom.element or "").strip()
    if el:
        return el.upper()
    name = "".join(ch for ch in atom.name if ch.isalpha())
    return name[:1].upper() if name else "?"


def is_heavy_atom(atom: Atom) -> bool:
    return atom_element(atom) not in {"H", "D"}


@dataclass
class ParsedEntry:
    pdb_id: str
    ligand_comp_id: str
    n_models: int
    polymer_residues: list[Residue] = field(default_factory=list)
    ligand_residues: list[Residue] = field(default_factory=list)
    solvent_residues: list[Residue] = field(default_factory=list)
    water_residues: list[Residue] = field(default_factory=list)
    other_het_residues: list[Residue] = field(default_factory=list)

    @property
    def protein_atoms(self) -> list[Atom]:
        return [a for r in self.polymer_residues for a in r if is_heavy_atom(a)]

    @property
    def ligand_atoms(self) -> list[Atom]:
        return [a for r in self.ligand_residues for a in r if is_heavy_atom(a)]

    def atom_coord_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        """(protein_coords, ligand_coords) as (N,3) float arrays, Angstrom."""
        prot = np.array([a.coord for a in self.protein_atoms], dtype=float)
        lig = np.array([a.coord for a in self.ligand_atoms], dtype=float)
        return prot, lig


def classify_model(model: Model, pdb_id: str, ligand_comp_id: str) -> ParsedEntry:
    entry = ParsedEntry(pdb_id=pdb_id, ligand_comp_id=ligand_comp_id, n_models=1)
    for chain in model:
        for residue in chain:
            hetfield, resseq, icode = residue.id
            resname = residue.get_resname().strip()
            if hetfield == " ":
                entry.polymer_residues.append(residue)
            elif hetfield == "W" or resname == "HOH":
                entry.water_residues.append(residue)
            elif resname == ligand_comp_id:
                entry.ligand_residues.append(residue)
            elif resname in SOLVENT_COMPONENT_IDS:
                entry.solvent_residues.append(residue)
            else:
                entry.other_het_residues.append(residue)
    return entry


def load_entry(pdb_path: Path, pdb_id: str | None = None) -> ParsedEntry:
    """Parse a PDB file; analyse model 0 only (X-ray entries are single-model)."""
    pdb_path = Path(pdb_path)
    pid = (pdb_id or pdb_path.stem).upper()
    ligand_comp_id = LIGAND_COMPONENTS.get(pid)
    structure = PDBParser(QUIET=True).get_structure(pid, str(pdb_path))
    models = list(structure.get_models())
    if not models:
        raise ValueError(f"No models parsed from {pdb_path}")
    entry = classify_model(models[0], pid, ligand_comp_id or "")
    entry.n_models = len(models)
    if ligand_comp_id is not None and not entry.ligand_residues:
        raise ValueError(
            f"Expected ligand component {ligand_comp_id} not found in {pdb_path}"
        )
    return entry


def ligand_element_counts(entry: ParsedEntry) -> Counter:
    """Heavy-atom element multiset of the extracted ligand — for identity checks."""
    return Counter(atom_element(a) for a in entry.ligand_atoms)


def extract_het_record_block(pdb_path: Path, resnames: set[str]) -> str:
    """Return the raw HETATM lines for the given component IDs plus all
    CONECT records, as a single PDB-format block (for RDKit parsing)."""
    atoms, conect = [], []
    for line in Path(pdb_path).read_text(errors="replace").splitlines():
        if line.startswith("HETATM") and line[17:20].strip() in resnames:
            atoms.append(line)
        elif line.startswith("CONECT"):
            conect.append(line)
    return "\n".join(atoms + conect) + "\n"
