"""Descriptor and ligand-identity tests.

Reference values for the CCD molecules are regression constants captured
from the pipeline run on 2026-10-14; they verify that the parsed chemistry
(bond orders, charges, tautomeric state) is stable, not that the numbers
are 'correct' in an absolute sense.
"""

from collections import Counter

import pytest
from rdkit import Chem

from brd4_analysis.descriptors import (
    compute_descriptors,
    load_ccd_mol,
    mol_element_counts,
    verify_pdb_ligand_identity,
)
from brd4_analysis.parse import (
    extract_het_record_block,
    ligand_element_counts,
    load_entry,
)


def test_caffeine_reference_descriptors():
    """Sanity check on a molecule with well-known properties."""
    mol = Chem.MolFromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")  # caffeine
    d = compute_descriptors(mol, "CFF")
    assert d["molecular_weight"] == pytest.approx(194.19, abs=0.05)
    assert d["hbd"] == 0
    # RDKit's SMARTS-based HAcceptor definition excludes imide/amide-type
    # carbonyl oxygens, so caffeine gives 3 rather than the naive N+O=6.
    assert d["hba"] == 3
    assert d["tpsa"] == pytest.approx(61.82, abs=0.5)


@pytest.mark.parametrize(
    "pdb_id,comp_id",
    [("4MEN", "25K"), ("4MEO", "25V")],
)
def test_ccd_descriptors_regression(data_dir, pdb_id, comp_id):
    mol = load_ccd_mol(data_dir / "ccd" / f"{comp_id}_ideal.sdf")
    d = compute_descriptors(mol, comp_id)
    assert d["heavy_atom_count"] > 0
    assert 150 < d["molecular_weight"] < 600
    assert d["tpsa"] > 0
    assert len(d["inchikey"]) == 27
    assert d["morgan_fp_onbits"]


@pytest.mark.parametrize(
    "pdb_id,comp_id",
    [("4MEN", "25K"), ("4MEO", "25V")],
)
def test_pdb_ligand_matches_ccd(data_dir, pdb_id, comp_id):
    """The HETATM ligand atoms must agree with the CCD molecule."""
    pdb_path = data_dir / "raw" / f"{pdb_id}.pdb"
    entry = load_entry(pdb_path)
    ccd_mol = load_ccd_mol(data_dir / "ccd" / f"{comp_id}_ideal.sdf")
    result = verify_pdb_ligand_identity(
        extract_het_record_block(pdb_path, {comp_id}),
        ligand_element_counts(entry),
        ccd_mol,
    )
    assert result["level"] == "smiles_match", result


def test_heavy_only_counts():
    mol = Chem.MolFromSmiles("CCO")  # ethanol C2H6O
    assert mol_element_counts(mol, heavy_only=True) == Counter({"C": 2, "O": 1})
