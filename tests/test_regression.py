"""Exact regression tests for the Phase 1 baseline.

Expected values are pinned constants from the reviewed 2026-10-14 run with
the environment in requirements.txt. Several are re-derived independently
here (from raw coordinates or SMILES) instead of being read back from
results/, so a pipeline bug cannot silently confirm itself.
"""

import json
from pathlib import Path

import pytest
from Bio.PDB import PDBParser
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator, rdMolDescriptors

from brd4_analysis.contacts import contacts_by_residue
from brd4_analysis.descriptors import compute_descriptors, load_ccd_mol, tanimoto_morgan
from brd4_analysis.fetch import load_provenance, sha256_file
from brd4_analysis.parse import extract_het_record_block, load_entry
from brd4_analysis.report import render_report

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"

EXPECTED_SHA256 = {
    "raw/4MEN.pdb": "a1b21bdc30365fb0ca905a6f9f00cc5955e3b382e6b726c02e4a2162be55b6a0",
    "raw/4MEN_entry.json": "0ab358f1c41ae1008405ac3613aab027d7a625cbfe3e44a1a3073ddf8c7cc29b",
    "ccd/25K_ideal.sdf": "2e978f1ec28868bcb8595e0001dd7964e8dc196a6f8a236a6c0dc6e3e48eb81b",
    "raw/4MEO.pdb": "16f7df00cb65a5658c5e3bbb97458327ce73a27993b7304f035753215fa1b3de",
    "raw/4MEO_entry.json": "6c47f4fb65484bbd3d87d835096a5ae82d31f5fb45533fd2813af4cd71e18cd0",
    "ccd/25V_ideal.sdf": "18482eef7c8936325b4ce3455b219af3eb679c503f36145ba6f623f0eec31ef0",
}

LIGANDS = {
    "4MEN": {
        "comp_id": "25K",
        "formula": "C15H17N5",
        "smiles": "Cc1ccc(CN(C)c2cc(C)nc3ncnn23)cc1",
        "inchikey": "JAEVUMCKPBRDAJ-UHFFFAOYSA-N",
        "molecular_weight": 267.336, "clogp": 2.37754, "tpsa": 46.32,
        "hbd": 0, "hba": 4, "hbd_lipinski_nhoh": 0, "hba_lipinski_no": 5,
        "rotatable_bonds": 3, "heavy_atom_count": 20, "aromatic_ring_count": 3,
        "formal_charge": 0,
    },
    "4MEO": {
        "comp_id": "25V",
        "formula": "C18H16N2O2",
        "smiles": "CC(=O)Nc1cccc(-c2cc(C)[n+]([O-])c3ccccc23)c1",
        "inchikey": "YSGOTUAMXVUKOD-UHFFFAOYSA-N",
        "molecular_weight": 292.338, "clogp": 3.40702, "tpsa": 56.04,
        "hbd": 1, "hba": 2, "hbd_lipinski_nhoh": 1, "hba_lipinski_no": 4,
        "rotatable_bonds": 2, "heavy_atom_count": 22, "aromatic_ring_count": 3,
        "formal_charge": 0,
    },
}

RESIDUE_COUNTS = {"4MEN": {3.5: 5, 4.0: 11, 5.0: 13}, "4MEO": {3.5: 5, 4.0: 11, 5.0: 14}}
RESIDUES_WITHIN_3_5 = {
    "4MEN": {"ASN140", "PHE83", "TRP81", "LEU92", "PRO82"},
    "4MEO": {"ASP88", "ASN140", "LEU92", "TRP81", "VAL87"},
}
# (protein residue, protein atom, ligand atom) -> distance in Angstrom
KEY_DISTANCES = {
    "4MEN": {(140, "ND2", "N2"): 2.844, (140, "OD1", "N2"): 3.427},
    "4MEO": {(140, "ND2", "O"): 2.881, (88, "N", "O1"): 2.670, (92, "CD1", "O1"): 2.910},
}


def _summary() -> dict:
    return json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))


# ---------- input integrity ----------

def test_input_checksums_pinned(data_dir):
    assert {k: sha256_file(data_dir / k) for k in EXPECTED_SHA256} == EXPECTED_SHA256


def test_provenance_manifest_portable_and_consistent(data_dir):
    prov = load_provenance(data_dir)
    records = [r for e in prov["entries"].values() for r in e.values()]
    assert len(records) == 6
    for r in records:
        assert not Path(r["path"]).is_absolute() and "\\" not in r["path"]
        assert EXPECTED_SHA256[r["path"]] == r["sha256"]


# ---------- ligand identity ----------

@pytest.mark.parametrize("pdb_id", LIGANDS)
def test_ligand_identity_exact(data_dir, pdb_id):
    exp = LIGANDS[pdb_id]
    entry = load_entry(data_dir / "raw" / f"{pdb_id}.pdb")
    assert {r.get_resname().strip() for r in entry.ligand_residues} == {exp["comp_id"]}
    mol = load_ccd_mol(data_dir / "ccd" / f"{exp['comp_id']}_ideal.sdf")
    assert rdMolDescriptors.CalcMolFormula(mol) == exp["formula"]
    assert Chem.MolToSmiles(mol) == exp["smiles"]
    assert Chem.MolToInchiKey(mol) == exp["inchikey"]
    assert Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False) == []


def test_25v_oxygen_roles(data_dir):
    """O1 is the acetamide carbonyl O; O is the N-oxide O (guards against
    the atom-role swap found in review)."""
    pdb = data_dir / "raw" / "4MEO.pdb"
    block = extract_het_record_block(pdb, {"25V"})
    raw = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
    mol = AllChem.AssignBondOrdersFromTemplate(load_ccd_mol(data_dir / "ccd" / "25V_ideal.sdf"), raw)
    by_name = {a.GetPDBResidueInfo().GetName().strip(): a for a in mol.GetAtoms()}

    def neighbours(name):
        a = by_name[name]
        return {(n.GetPDBResidueInfo().GetName().strip(), str(mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondType())) for n in a.GetNeighbors()}

    assert neighbours("O1") == {("C1", "DOUBLE")}
    assert neighbours("O") == {("N1", "SINGLE")}
    assert by_name["O"].GetFormalCharge() == -1 and by_name["N1"].GetFormalCharge() == 1
    lig = [r for r in PDBParser(QUIET=True).get_structure("x", str(pdb))[0]["A"] if r.get_resname() == "25V"][0]
    assert lig["O1"] - lig["C1"] == pytest.approx(1.227, abs=1e-3)
    assert lig["O"] - lig["N1"] == pytest.approx(1.358, abs=1e-3)


# ---------- descriptors ----------

@pytest.mark.parametrize("pdb_id", LIGANDS)
def test_reference_descriptors_exact(data_dir, pdb_id):
    exp = LIGANDS[pdb_id]
    d = compute_descriptors(load_ccd_mol(data_dir / "ccd" / f"{exp['comp_id']}_ideal.sdf"), exp["comp_id"])
    for key in ("molecular_weight", "clogp", "tpsa"):
        assert d[key] == pytest.approx(exp[key], abs=5e-4), key
    for key in ("hbd", "hba", "hbd_lipinski_nhoh", "hba_lipinski_no", "rotatable_bonds",
                "heavy_atom_count", "aromatic_ring_count", "formal_charge"):
        assert d[key] == exp[key], key


def test_hba_definition_pinned():
    """CalcNumHBA (behind Lipinski.NumHAcceptors) counts carbonyl O and
    pyridine-type n but, from RDKit 2025.09, not N-substituted aromatic n.
    Caffeine: 2 C=O + imidazole N3 = 3. Lipinski N+O = 6. If this fails after
    an RDKit upgrade, re-check the documented HBA definition."""
    caffeine = Chem.MolFromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    assert rdMolDescriptors.CalcNumHBA(caffeine) == 3
    assert rdMolDescriptors.CalcNumLipinskiHBA(caffeine) == 6


# ---------- fingerprint similarity ----------

def test_tanimoto_exact(data_dir):
    sim = tanimoto_morgan(load_ccd_mol(data_dir / "ccd" / "25K_ideal.sdf"),
                          load_ccd_mol(data_dir / "ccd" / "25V_ideal.sdf"))
    assert (sim["on_bits_a"], sim["on_bits_b"], sim["common_bits"], sim["union_bits"]) == (38, 40, 6, 72)
    assert sim["tanimoto"] == 6 / 72
    assert _summary()["ligand_similarity"]["morgan_ecfp4_tanimoto_4MEN_vs_4MEO"] == 6 / 72


def test_tanimoto_independent_from_smiles():
    g = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    a, b = (g.GetFingerprint(Chem.MolFromSmiles(LIGANDS[p]["smiles"])) for p in LIGANDS)
    assert DataStructs.TanimotoSimilarity(a, b) == 6 / 72


# ---------- contact geometry ----------

@pytest.mark.parametrize("pdb_id", LIGANDS)
def test_key_distances_from_raw_coordinates(data_dir, pdb_id):
    chain = PDBParser(QUIET=True).get_structure(pdb_id, str(data_dir / "raw" / f"{pdb_id}.pdb"))[0]["A"]
    lig = [r for r in chain if r.get_resname() == LIGANDS[pdb_id]["comp_id"]][0]
    for (resseq, patom, latom), expected in KEY_DISTANCES[pdb_id].items():
        assert chain[resseq][patom] - lig[latom] == pytest.approx(expected, abs=5e-4)


@pytest.mark.parametrize("pdb_id", LIGANDS)
def test_residue_contact_counts_exact(data_dir, pdb_id):
    contacts = contacts_by_residue(load_entry(data_dir / "raw" / f"{pdb_id}.pdb"))
    for t, n in RESIDUE_COUNTS[pdb_id].items():
        assert int(contacts[f"within_{t}A"].sum()) == n, t
    close = {f"{r.resname}{r.resseq}" for r in contacts[contacts["within_3.5A"]].itertuples()}
    assert close == RESIDUES_WITHIN_3_5[pdb_id]


def test_vdw_overlap_screen_flags_only_4meo_leu92():
    s = _summary()
    assert s["4MEN"]["vdw_overlap_screen"]["pairs"] == []
    flagged = {(p["protein_atom"], p["ligand_atom"]) for p in s["4MEO"]["vdw_overlap_screen"]["pairs"]}
    assert flagged == {("A:LEU92:CD1", "25V:O1"), ("A:LEU92:CD1", "25V:C3")}


# ---------- report / summary consistency ----------

def test_summary_has_no_run_variable_fields():
    text = (RESULTS / "summary.json").read_text(encoding="utf-8")
    assert "generated_utc" not in text and "retrieved_utc" not in text


def test_summary_inputs_match_pinned_checksums():
    assert _summary()["inputs_sha256"] == EXPECTED_SHA256


def test_report_numbers_match_results(data_dir):
    report = (PROJECT_ROOT / "docs" / "REPORT.md").read_text(encoding="utf-8").replace("\r\n", "\n")
    assert render_report(report, _summary(), load_provenance(data_dir)) == report
    assert "**0.0833**" in report
