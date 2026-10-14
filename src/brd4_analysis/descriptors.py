"""Ligand physicochemical descriptors and chemical identity verification.

Descriptors are computed from the Chemical Component Dictionary (CCD)
ideal-coordinate structure, which carries authoritative bond orders and
formal charges. The ligand atoms extracted from the experimental PDB file
(X-ray: heavy atoms only, no explicit hydrogens) are then checked against
the CCD molecule, primarily by full-graph template matching and secondly
by heavy-element multiset comparison.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem,
    Crippen,
    Descriptors,
    Lipinski,
    MACCSkeys,
    rdFingerprintGenerator,
    rdMolDescriptors,
)


def load_ccd_mol(sdf_path: Path) -> Chem.Mol:
    # removeHs=True: heavy-atom graph only. This matches the X-ray PDB
    # representation (no explicit H records) so the CCD molecule can be
    # used directly as a template for bond-order assignment.
    suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=True)
    mols = [m for m in suppl if m is not None]
    if len(mols) != 1:
        raise ValueError(f"Expected 1 molecule in {sdf_path}, parsed {len(mols)}")
    return mols[0]


def mol_element_counts(mol: Chem.Mol, heavy_only: bool = True) -> Counter:
    return Counter(
        a.GetSymbol()
        for a in mol.GetAtoms()
        if not heavy_only or a.GetAtomicNum() > 1
    )


def morgan_fingerprint(mol: Chem.Mol, radius: int = 2, nbits: int = 2048):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    return gen.GetFingerprint(mol)


def maccs_fingerprint(mol: Chem.Mol):
    return MACCSkeys.GenMACCSKeys(mol)


def compute_descriptors(mol: Chem.Mol, comp_id: str) -> dict[str, Any]:
    """Descriptor set required for the baseline comparison table."""
    morgan = morgan_fingerprint(mol)
    maccs = maccs_fingerprint(mol)
    return {
        "component_id": comp_id,
        "ccd_formula": rdMolDescriptors.CalcMolFormula(mol),
        "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=False),
        "isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "inchikey": Chem.MolToInchiKey(mol),
        "molecular_weight": Descriptors.MolWt(mol),
        "exact_mw": Descriptors.ExactMolWt(mol),
        "clogp": Crippen.MolLogP(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "heavy_atom_count": Lipinski.HeavyAtomCount(mol),
        "ring_count": Lipinski.RingCount(mol),
        "aromatic_ring_count": Lipinski.NumAromaticRings(mol),
        "fraction_csp3": Lipinski.FractionCSP3(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "morgan_fp_radius": 2,
        "morgan_fp_nbits": 2048,
        "morgan_fp_onbits": sorted(morgan.GetOnBits()),
        "maccs_fp_onbits": sorted(maccs.GetOnBits()),
    }


def verify_pdb_ligand_identity(
    pdb_het_block: str,
    pdb_element_counts: Counter,
    ccd_mol: Chem.Mol,
) -> dict[str, Any]:
    """Check the PDB-extracted ligand atoms against the CCD molecule.

    Returns a dict with ``level``:
      * ``smiles_match``  — template-assigned PDB ligand gives the same
        canonical isomeric SMILES as the CCD molecule (strongest check);
      * ``formula_match`` — heavy-element multisets agree but template
        assignment failed (weaker; recorded as reduced confidence);
      * ``mismatch``      — element counts disagree.
    """
    ccd_counts = mol_element_counts(ccd_mol, heavy_only=True)
    result: dict[str, Any] = {
        "pdb_heavy_element_counts": dict(pdb_element_counts),
        "ccd_heavy_element_counts": dict(ccd_counts),
        "formula_match": pdb_element_counts == ccd_counts,
    }

    pdb_mol = Chem.MolFromPDBBlock(pdb_het_block, sanitize=False, removeHs=False)
    smiles_match = False
    if pdb_mol is not None and pdb_mol.GetNumAtoms() > 0:
        try:
            assigned = AllChem.AssignBondOrdersFromTemplate(ccd_mol, pdb_mol)
            smiles_match = (
                assigned is not None
                and Chem.MolToSmiles(assigned, isomericSmiles=True)
                == Chem.MolToSmiles(ccd_mol, isomericSmiles=True)
            )
        except Exception:
            smiles_match = False
    result["smiles_match"] = smiles_match

    if smiles_match:
        result["level"] = "smiles_match"
    elif result["formula_match"]:
        result["level"] = "formula_match"
    else:
        result["level"] = "mismatch"
    return result


def tanimoto_morgan(mol_a: Chem.Mol, mol_b: Chem.Mol) -> float:
    """ECFP4 (Morgan r=2, 2048-bit) Tanimoto between two ligands."""
    return DataStructs.TanimotoSimilarity(morgan_fingerprint(mol_a), morgan_fingerprint(mol_b))
