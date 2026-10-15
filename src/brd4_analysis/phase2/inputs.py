"""Construct and chemical-graph verification; strict Boltz mmCIF adapter."""

from __future__ import annotations

from importlib.metadata import version
from pathlib import Path

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem

from brd4_analysis.contacts import residue_label
from brd4_analysis.descriptors import compute_descriptors, load_ccd_mol
from brd4_analysis.parse import atom_element, is_heavy_atom, load_entry
from brd4_analysis.phase2.config import RDKIT_VERSION, read_json, source_path
from brd4_analysis.phase2.metrics import coordinates, graph_mappings


def ccd_graph(path: Path) -> Chem.Mol:
    """Read named heavy atoms and ideal coordinates from the pinned CCD mmCIF."""
    data = MMCIF2Dict(str(path))
    rw = Chem.RWMol()
    indices, ideal = {}, []
    atom_rows = zip(
        data["_chem_comp_atom.atom_id"], data["_chem_comp_atom.type_symbol"],
        data["_chem_comp_atom.charge"], data["_chem_comp_atom.pdbx_model_Cartn_x_ideal"],
        data["_chem_comp_atom.pdbx_model_Cartn_y_ideal"],
        data["_chem_comp_atom.pdbx_model_Cartn_z_ideal"],
    )
    for name, element, charge, x, y, z in atom_rows:
        if element in {"H", "D"}:
            continue
        if name in indices:
            raise ValueError("Duplicate CCD atom name")
        atom = Chem.Atom(element)
        atom.SetFormalCharge(int(charge))
        atom.SetProp("name", name)
        indices[name] = rw.AddAtom(atom)
        ideal.append([float(x), float(y), float(z)])
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for a, b, order in zip(
        data["_chem_comp_bond.atom_id_1"], data["_chem_comp_bond.atom_id_2"],
        data["_chem_comp_bond.value_order"],
    ):
        if a in indices and b in indices:
            rw.AddBond(indices[a], indices[b], orders[order])
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(len(ideal))
    for i, point in enumerate(coordinates(ideal)):
        conf.SetAtomPosition(i, point)
    mol.AddConformer(conf)
    return mol


def boltz_ligand(smiles: str) -> Chem.Mol:
    """Reproduce v2.2.1 schema.py naming, without embedding or importing Boltz.

    The version gate matters: canonical ranks are not a cross-version atom ID.
    Runtime preflight must reproduce this mapping in the prediction environment.
    """
    if version("rdkit") != RDKIT_VERSION:
        raise ValueError(f"Atom naming requires RDKit {RDKIT_VERSION}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Expected one valid, connected ligand SMILES")
    mol = Chem.AddHs(mol)
    ranks = Chem.CanonicalRankAtoms(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    for atom, rank in zip(mol.GetAtoms(), ranks):
        atom.SetProp("name", atom.GetSymbol().upper() + str(rank + 1))
    return Chem.RemoveHs(mol)


def atom_names(mol):
    return [a.GetProp("name") for a in mol.GetAtoms()]


def reference_data(config: dict) -> dict:
    refs = {key: source_path(spec) for key, spec in config["reference"].items()}
    metadata = read_json(refs["polymer_metadata"])
    if metadata["entity_poly"]["pdbx_seq_one_letter_code_can"] != config["protein"]["sequence"]:
        raise ValueError("RCSB entity sequence disagrees with the preregistered construct")
    alignments = metadata["rcsb_polymer_entity_align"]
    if not any(
        a["reference_database_accession"] == "O60885"
        and a["aligned_regions"] == [{"entity_beg_seq_id": 3, "length": 125, "ref_beg_seq_id": 44}]
        for a in alignments
    ):
        raise ValueError("Unreviewed UniProt/construct mapping")
    raw = refs["pdb"].read_text(encoding="ascii")
    tags = [line for line in raw.splitlines() if line.startswith("SEQADV")]
    if not any("SER A   42" in line and "EXPRESSION TAG" in line for line in tags) or not any(
        "MET A   43" in line and "EXPRESSION TAG" in line for line in tags
    ):
        raise ValueError("SM cloning-remnant annotation not found")
    entry = load_entry(refs["pdb"], "4MEN")
    if entry.n_models != 1 or len(entry.ligand_residues) != 1:
        raise ValueError("Ambiguous reference model or ligand instance")
    residues = entry.polymer_residues
    if len(residues) != 127 or any(r.get_parent().id != "A" for r in residues):
        raise ValueError("Reference must contain one complete chain A")
    if [(r.id[1], r.id[2]) for r in residues] != [(i, " ") for i in range(42, 169)]:
        raise ValueError("Unexpected author numbering, insertion codes, gap or residue order")
    sequence = "".join(protein_letters_3to1[r.resname] for r in residues)
    if sequence != config["protein"]["sequence"]:
        raise ValueError("Observed sequence and construct sequence differ")
    if any("CA" not in r for r in residues):
        raise ValueError("Missing reference CA")
    mapping = [
        {
            "input_index_1": i, "reference_auth_chain": "A",
            "reference_auth_resseq": r.id[1], "reference_icode": r.id[2].strip(),
            "reference_label": residue_label(r), "resname": r.resname,
            "one_letter": protein_letters_3to1[r.resname],
            "uniprot_accession": "O60885" if i > 2 else None,
            "uniprot_position": r.id[1] if i > 2 else None,
            "fit_CA": i > 2,
        }
        for i, r in enumerate(residues, 1)
    ]
    mol = ccd_graph(refs["ccd_cif"])
    sdf = load_ccd_mol(refs["ccd_sdf"])
    predicted_mol = boltz_ligand(config["ligand"]["smiles"])
    graph_mappings(mol, sdf)
    graph_mappings(mol, predicted_mol)
    descriptors = compute_descriptors(predicted_mol, "25K")
    if (
        descriptors["inchikey"] != config["ligand"]["expected_inchikey"]
        or descriptors["ccd_formula"] != "C15H17N5"
        or descriptors["heavy_atom_count"] != 20
        or descriptors["formal_charge"] != 0
        or Chem.FindMolChiralCenters(predicted_mol, includeUnassigned=True)
    ):
        raise ValueError("SMILES chemistry is not the preregistered neutral, achiral 25K graph")
    ligand_atoms = {a.name: a for a in entry.ligand_atoms}
    if set(ligand_atoms) != set(atom_names(mol)):
        raise ValueError("Reference ligand has missing, duplicated or extra heavy atoms")
    if any(atom_element(ligand_atoms[a.GetProp("name")]) != a.GetSymbol().upper() for a in mol.GetAtoms()):
        raise ValueError("Reference CCD atom-name/element mismatch")
    protein_by_residue = {
        residue_label(r): {a.name: a for a in r if is_heavy_atom(a)} for r in residues
    }
    return {
        "mapping": mapping, "mol": mol, "predicted_mol": predicted_mol,
        "ligand_xyz": coordinates([ligand_atoms[name].coord for name in atom_names(mol)]),
        "ideal_xyz": coordinates(mol.GetConformer().GetPositions()),
        "ca": coordinates([r["CA"].coord for r in residues]),
        "protein_by_residue": protein_by_residue,
        "descriptors": descriptors,
    }


def read_prediction(path: Path, reference: dict) -> dict:
    """Require one complete 127-residue protein and one 20-heavy-atom ligand.

    Polymer mapping uses label_seq_id (Boltz input index), never an inferred
    constant offset or a proximity-based bond assignment.
    """
    data = MMCIF2Dict(str(path))
    fields = [
        "group_PDB", "type_symbol", "label_atom_id", "label_comp_id",
        "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
        "label_alt_id", "occupancy", "pdbx_PDB_model_num",
    ]
    try:
        columns = [data[f"_atom_site.{field}"] for field in fields]
    except KeyError as exc:
        raise ValueError(f"Prediction lacks required mmCIF column: {exc}") from exc
    if len({len(c) for c in columns}) != 1:
        raise ValueError("Inconsistent mmCIF atom column lengths")
    if len(set(columns[-1])) != 1:
        raise ValueError("Multiple prediction models are not supported")
    protein, ligand, ligand_components = {}, {}, set()
    for row in zip(*columns):
        group, element, name, resname, chain, seqid, x, y, z, alt, occupancy, _ = row
        element = element.upper()
        if element in {"H", "D"}:
            continue
        if alt not in {".", "?"} or float(occupancy) != 1:
            raise ValueError("Predicted alternate conformers or partial occupancy are unsupported")
        point = coordinates([[float(x), float(y), float(z)]])[0]
        if chain == "A" and group == "ATOM":
            index = int(seqid)
            if not 1 <= index <= 127 or resname != reference["mapping"][index - 1]["resname"]:
                raise ValueError("Prediction protein residue mapping or identity mismatch")
            atoms = protein.setdefault(index, {})
            if name in atoms:
                raise ValueError("Duplicate protein atom")
            atoms[name] = (point, element)
        elif chain == "B" and group == "HETATM":
            ligand_components.add(resname)
            if name in ligand:
                raise ValueError("Duplicate ligand atom")
            ligand[name] = (point, element)
        else:
            raise ValueError("Unexpected chain, solvent, additive or atom record in prediction")
    if set(protein) != set(range(1, 128)) or len(ligand_components) != 1:
        raise ValueError("Missing protein residue or ambiguous ligand instance")
    model_mol = reference["predicted_mol"]
    if set(ligand) != set(atom_names(model_mol)):
        raise ValueError("Missing/extra ligand atoms or unverified Boltz atom-name mapping")
    if any(ligand[a.GetProp("name")][1] != a.GetSymbol().upper() for a in model_mol.GetAtoms()):
        raise ValueError("Prediction ligand atom-name/element mismatch")
    by_residue, elements, all_xyz, ca = {}, [], [], []
    for item in reference["mapping"]:
        index, label = item["input_index_1"], item["reference_label"]
        atoms = protein[index]
        expected = reference["protein_by_residue"][label]
        extras = set(atoms) - set(expected)
        if not set(expected) <= set(atoms) or (extras and not (index == 127 and extras == {"OXT"})):
            raise ValueError(f"Missing/extra protein heavy atoms at input residue {index}")
        if any(atoms[name][1] != atom_element(atom) for name, atom in expected.items()):
            raise ValueError("Protein atom-name/element mismatch")
        points = coordinates([xyz for xyz, _ in atoms.values()])
        by_residue[label] = points
        all_xyz.extend(points)
        elements.extend(element for _, element in atoms.values())
        ca.append(atoms["CA"][0])
    return {
        "ca": coordinates(ca), "protein_by_residue": by_residue,
        "protein_xyz": coordinates(all_xyz), "protein_elements": elements,
        "ligand_xyz": coordinates([ligand[name][0] for name in atom_names(model_mol)]),
    }
