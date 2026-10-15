"""Parsing and entity-classification tests against committed raw data."""

import numpy as np

from brd4_analysis.parse import load_entry


def test_4men_ligand_is_25k(pdb_path_4men):
    entry = load_entry(pdb_path_4men)
    assert {r.get_resname().strip() for r in entry.ligand_residues} == {"25K"}
    assert len(entry.ligand_residues) >= 1
    # X-ray entry: ligand atoms are heavy atoms only
    assert len(entry.ligand_atoms) > 0


def test_4meo_ligand_is_25v(pdb_path_4meo):
    entry = load_entry(pdb_path_4meo)
    assert {r.get_resname().strip() for r in entry.ligand_residues} == {"25V"}


def test_protein_excludes_water_and_solvent(pdb_path_4men):
    entry = load_entry(pdb_path_4men)
    prot_resnames = {r.get_resname().strip() for r in entry.polymer_residues}
    assert "HOH" not in prot_resnames
    assert "EDO" not in prot_resnames
    assert "25K" not in prot_resnames
    assert len(entry.water_residues) > 0
    assert len(entry.solvent_residues) > 0  # EDO present in 4MEN


def test_coordinate_units_sane(pdb_path_4men):
    """PDB coordinates are Angstrom: protein-ligand distances must be on a
    molecular scale (not nm, not microns)."""
    entry = load_entry(pdb_path_4men)
    prot, lig = entry.atom_coord_arrays()
    d = np.sqrt(((prot[:, None, :] - lig[None, :, :]) ** 2).sum(-1)).min()
    assert 1.0 < d < 10.0  # bound ligand must sit inside the binding site
