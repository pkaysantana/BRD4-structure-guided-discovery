"""Small offline controls for errors that would invalidate the GPU pilot."""

from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from brd4_analysis.constants import PROJECT_ROOT
from brd4_analysis.phase2.audit import later_case_review
from brd4_analysis.phase2.config import read_json, sha256, source_path, validate_a3m, validate_config
from brd4_analysis.phase2.inputs import atom_names, reference_data
from brd4_analysis.phase2.metrics import (
    contact_residues, contact_scores, geometry_checks, graph_mappings,
    protein_alignment, rmsd, symmetry_rmsd,
)
from brd4_analysis.phase2.workflow import (
    planned_command, prepare, score_prediction, validate, verify_run,
)

CONFIG = PROJECT_ROOT / "configs/phase2/4men_smoke.json"


@pytest.fixture(scope="module")
def config():
    return read_json(CONFIG)


@pytest.fixture(scope="module")
def reference(config):
    return reference_data(validate_config(config))


def synthetic_prediction(path, reference, ligand_shift=(0, 0, 0), omit=None, extra_atoms=()):
    """Known reference coordinates moved rigidly, in Boltz's actual naming scheme.

    This is a synthetic adapter control, never represented as a model prediction.
    Atom rows are reversed to ensure the adapter does not depend on file order.
    Optional extra coordinates are synthetic coverage/contact-accounting probes,
    not reconstructed experimental side chains.
    """
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    translation = np.array([12.3, -8.2, 4.7])
    fields = [
        "group_PDB", "type_symbol", "label_atom_id", "label_comp_id",
        "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
        "label_alt_id", "occupancy", "pdbx_PDB_model_num",
    ]
    rows = []
    for item in reference["mapping"]:
        for name, atom in reference["protein_by_residue"][item["reference_label"]].items():
            if (item["input_index_1"], name) == omit:
                continue
            point = atom.coord.astype(float) @ rotation + translation
            rows.append(["ATOM", atom.element, name, item["resname"], "A",
                         item["input_index_1"], *point, ".", 1, 1])
    for index, name, element, original_point in extra_atoms:
        point = np.asarray(original_point) @ rotation + translation
        rows.append(["ATOM", element, name, reference["mapping"][index - 1]["resname"],
                     "A", index, *point, ".", 1, 1])
    mapping = graph_mappings(reference["mol"], reference["predicted_mol"])[0]
    pred_names = atom_names(reference["predicted_mol"])
    for ref_index, pred_index in enumerate(mapping):
        point = reference["ligand_xyz"][ref_index] @ rotation + translation + ligand_shift
        atom = reference["predicted_mol"].GetAtomWithIdx(pred_index)
        rows.append(["HETATM", atom.GetSymbol(), pred_names[pred_index], "LIG1", "B",
                     ".", *point, ".", 1, 1])
    text = "data_synthetic_control\nloop_\n"
    text += "".join(f"_atom_site.{field}\n" for field in fields)
    text += "".join(" ".join(map(str, row)) + "\n" for row in reversed(rows)) + "#\n"
    path.write_text(text, encoding="utf-8")
    return path


def test_construct_and_chemistry_are_independent_of_smiles_spelling(config, reference):
    alternate = deepcopy(config)
    alternate["ligand"]["smiles"] = "Cc1ccc(cc1)CN(C)c2cc(nc3n2ncn3)C"
    other = reference_data(validate_config(alternate))
    assert other["descriptors"]["inchikey"] == reference["descriptors"]["inchikey"]
    assert reference["mapping"][0]["uniprot_position"] is None
    assert reference["mapping"][2]["uniprot_position"] == 44
    assert reference["mapping"][-1]["uniprot_position"] == 168
    for number, residue in [(81, "TRP"), (88, "ASP"), (97, "TYR"), (140, "ASN")]:
        assert reference["mapping"][number - 42]["resname"] == residue
    assert reference["descriptors"]["ccd_formula"] == "C15H17N5"
    assert reference["descriptors"]["hba"] == 4


@pytest.mark.parametrize("section,key,value", [
    ("boltz", "diffusion_samples", 5),
    ("boltz", "seed", True),
    ("boltz", "sampling_steps", -1),
    ("boltz", "recycling_steps", 10),
    ("boltz", "affinity", True),
    ("boltz", "templates", [{"cif": "4MEN.cif"}]),
    ("boltz", "version", "latest"),
    ("msa", "mode", "server"),
    ("protein", "sequence", "AAAA"),
    ("metrics", "min_contact_recall", float("nan")),
    ("metrics", "min_contact_precision", 1.1),
    ("resources", "wall_time_limit_minutes", 60),
])
def test_unreviewed_configuration_rejected(config, section, key, value):
    bad = deepcopy(config)
    bad[section][key] = value
    with pytest.raises(ValueError):
        validate_config(bad)


def test_unknown_keys_and_reference_tampering_rejected(config):
    bad = deepcopy(config)
    bad["boltz"]["use_msa_server"] = True
    with pytest.raises(ValueError):
        validate_config(bad)
    bad = deepcopy(config)
    bad["reference"]["pdb"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="changed source"):
        validate_config(bad)
    bad["reference"]["pdb"]["path"] = "../outside.pdb"
    with pytest.raises(ValueError, match="escapes"):
        validate_config(bad)


def test_a3m_query_mapping_and_insertion_handling(tmp_path):
    path = tmp_path / "test.a3m"
    path.write_text(">query\nACD\n>homolog\nAcC-\n")
    assert validate_a3m(path, "ACD") == 2
    path.write_text(">query\nCD\n")
    with pytest.raises(ValueError, match="query"):
        validate_a3m(path, "ACD")
    path.write_text(">query\nACD\n>bad\nACDE\n")
    with pytest.raises(ValueError, match="length"):
        validate_a3m(path, "ACD")


def test_kabsch_rotation_translation_and_reflection():
    fixed = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
    moving = fixed @ np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]) + 5
    rotation, translation, error = protein_alignment(moving, fixed)
    assert np.linalg.det(rotation) == pytest.approx(1)
    assert error < 1e-12
    assert rmsd(moving @ rotation + translation, fixed) < 1e-12
    reflected = fixed * [-1, 1, 1]
    assert protein_alignment(reflected, fixed)[2] > 0.1


@pytest.mark.parametrize("points", [
    [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
    [[0, 0, 0], [1, 0, 0]],
    [[0, 0, 0], [1, float("nan"), 0], [0, 1, 0]],
])
def test_invalid_alignment_is_not_silently_scored(points):
    with pytest.raises(ValueError):
        protein_alignment(points, points)


def test_symmetry_without_ligand_refitting():
    mol = Chem.MolFromSmiles("CC")
    fixed = np.array([[0., 0., 0.], [1.5, 0., 0.]])
    assert symmetry_rmsd(mol, mol, fixed, fixed[::-1])["rmsd_A"] == 0
    translated = symmetry_rmsd(mol, mol, fixed, fixed + [0, 4, 0])
    assert translated["rmsd_A"] == 4  # a ligand-only fit would falsely return zero
    with pytest.raises(ValueError, match="mismatch"):
        symmetry_rmsd(mol, Chem.MolFromSmiles("CN"), fixed, fixed)
    with pytest.raises(ValueError, match="cap"):
        graph_mappings(mol, mol, cap=1)


def test_contacts_penalize_missing_and_extra_residues_and_empty_predictions():
    score = contact_scores({"A:ASN140", "A:TRP81"}, {"A:ASN140", "A:ASP88", "A:LEU92"})
    assert score["precision"] == pytest.approx(1 / 3)
    assert score["recall"] == 0.5
    assert (score["false_positive"], score["false_negative"]) == (2, 1)
    assert contact_scores({"A:ASN140"}, set())["precision"] == 0
    assert contact_scores({"A:ASN140"}, set())["recall"] == 0
    assert contact_scores(set(), set())["recall"] is None
    assert contact_scores(set(), set())["precision"] is None
    protein = {"boundary": [[4, 0, 0]], "outside": [[4.00001, 0, 0]]}
    assert contact_residues(protein, [[0, 0, 0]], 4.0) == {"boundary"}


def test_geometry_rejects_distortion_and_severe_protein_clash(config, reference):
    ideal = reference["ideal_xyz"]
    far_protein = np.array([[100., 100., 100.]])
    args = (reference["mol"], ideal, ideal, far_protein, ["C"], config["metrics"]["geometry"])
    assert geometry_checks(*args)["valid"]
    damaged = ideal.copy()
    damaged[0] += [5, 0, 0]
    result = geometry_checks(reference["mol"], damaged, ideal, far_protein, ["C"], config["metrics"]["geometry"])
    assert not result["valid"] and "bond_lengths" in result["failed_checks"]
    result = geometry_checks(reference["mol"], ideal, ideal, ideal[:1], ["C"], config["metrics"]["geometry"])
    assert not result["valid"] and "protein_ligand_clash" in result["failed_checks"]
    ring = next(r for r in reference["mol"].GetRingInfo().AtomRings() if len(r) == 6)
    damaged = ideal.copy()
    damaged[ring[0]] += [0, 0, 2]
    result = geometry_checks(reference["mol"], damaged, ideal, far_protein, ["C"], config["metrics"]["geometry"])
    assert "aromatic_planarity" in result["failed_checks"]


def test_complete_adapter_scores_synthetic_rigid_transform_and_displaced_pose(tmp_path, config, reference):
    correct = synthetic_prediction(tmp_path / "synthetic_correct.cif", reference)
    good = score_prediction(config, correct)
    assert good["ligand_pose"]["rmsd_A"] < 1e-10
    assert good["protein_ca_rmsd_A"] < 1e-10
    assert good["contacts"]["4.0"]["precision"] == 1
    assert good["contacts"]["4.0"]["recall"] == 1
    assert good["harness_success"]
    wrong = synthetic_prediction(tmp_path / "synthetic_displaced.cif", reference, [50, 0, 0])
    bad = score_prediction(config, wrong)
    assert bad["ligand_pose"]["rmsd_A"] > 40
    assert bad["contacts"]["4.0"]["recall"] == 0
    assert not bad["harness_success"]


def test_terminal_oxt_can_be_absent_from_boltz_canonical_output(tmp_path, config, reference):
    # Boltz 2.2.1's canonical GLU atom set excludes OXT; 4MEN contains it.
    oxt = reference["protein_by_residue"]["A:GLU168"]["OXT"].coord.astype(float)
    assert np.linalg.norm(reference["ligand_xyz"] - oxt, axis=1).min() > 8
    path = synthetic_prediction(tmp_path / "canonical_terminal_atoms.cif", reference, omit=(127, "OXT"))
    score = score_prediction(config, path)
    assert score["harness_success"]
    assert score["ligand_pose"]["rmsd_A"] < 1e-10
    assert score["contacts"]["4.0"]["precision"] == 1
    assert score["contacts"]["4.0"]["recall"] == 1
    assert score["atom_coverage"]["optional_terminal_atoms_missing"] == ["A:GLU168:OXT"]


def test_canonical_atoms_missing_from_reference_are_reported_without_hiding_contacts(tmp_path, config, reference):
    # All six atoms exist in Boltz's canonical residue sets but lack coordinates
    # in 4MEN. Artificial 3.9 A contacts test accounting, not side-chain physics.
    probe = reference["ligand_xyz"][reference["ligand_xyz"][:, 0].argmin()] - [3.9, 0, 0]
    extras = [
        (50, "CD", "C", probe), (50, "CE", "C", probe), (50, "NZ", "N", probe),
        (82, "CD", "C", probe), (82, "NE2", "N", probe), (82, "OE1", "O", probe),
    ]
    path = synthetic_prediction(tmp_path / "full_canonical.cif", reference, omit=(127, "OXT"), extra_atoms=extras)
    score = score_prediction(config, path)
    assert score["ligand_pose"]["rmsd_A"] < 1e-10
    assert set(score["atom_coverage"]["predicted_atoms_not_observed_in_reference"]) == {
        "A:LYS91:CD", "A:LYS91:CE", "A:LYS91:NZ", "A:GLN123:CD", "A:GLN123:NE2", "A:GLN123:OE1",
    }
    assert score["contacts"]["4.0"]["false_positive"] == 2  # Original primary definition retained.
    assert score["contacts"]["4.0"]["precision"] == pytest.approx(11 / 13)
    assert score["contacts_reference_observed_atoms"]["4.0"]["precision"] == 1
    assert score["contacts_reference_observed_atoms"]["4.0"]["recall"] == 1
    assert score["geometry"]["valid"]


@pytest.mark.parametrize("name,element", [("ZZ", "C"), ("NZ", "C")])
def test_additional_atoms_require_a_valid_residue_template_and_element(tmp_path, config, reference, name, element):
    path = synthetic_prediction(
        tmp_path / "invalid_extra.cif", reference,
        extra_atoms=[(50, name, element, reference["ligand_xyz"][0])],
    )
    with pytest.raises(ValueError, match="[Pp]rotein"):
        score_prediction(config, path)


@pytest.mark.parametrize("missing", [(99, "CA"), (99, "ND2"), (127, "O")])
def test_adapter_rejects_missing_ca_or_contact_atom(tmp_path, config, reference, missing):
    path = synthetic_prediction(tmp_path / "incomplete.cif", reference, omit=missing)
    with pytest.raises(ValueError, match="Missing/extra"):
        score_prediction(config, path)


def test_preparation_is_offline_immutable_and_contains_no_reference_template(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Offline preparation attempted network or process execution")

    monkeypatch.setattr("requests.sessions.Session.request", forbidden)
    monkeypatch.setattr("subprocess.run", forbidden)
    run = prepare(CONFIG, tmp_path)
    frozen = verify_run(run)
    document = read_json(run / "inputs/4men.yaml")
    assert set(document) == {"version", "sequences"}
    assert document["sequences"][0]["protein"]["msa"] == "empty"
    assert set(document["sequences"][1]["ligand"]) == {"id", "smiles"}
    assert "--use_msa_server" not in planned_command(frozen)
    assert not (run / "outputs").exists()
    with pytest.raises(FileExistsError):
        prepare(CONFIG, tmp_path)
    (run / "inputs/4men.yaml").write_text("{}")
    with pytest.raises(ValueError, match="modified"):
        verify_run(run)
    with pytest.raises(ValueError, match="must not overwrite"):
        prepare(CONFIG, PROJECT_ROOT / "results")


def test_reference_validation_and_duplicate_later_case():
    _, _, report = validate(CONFIG)
    assert report["reference_geometry_screen"]["valid"]
    assert {key: len(value) for key, value in report["reference_contacts"].items()} == {
        "3.5": 5, "4.0": 11, "5.0": 13,
    }
    assert report["interface"]["runtime_cli_executed"] is False
    later = later_case_review()
    assert later["all_compared_coordinates_exactly_equal"]
    assert later["protein_CA_atoms_compared"] == 127
    assert later["ligand_heavy_atoms_compared"] == 19
    assert later["entries"]["9FWX"]["initial_release_date"].startswith("2024-08-07")


def test_pinned_primary_is_preparable_but_unauthorised_while_smoke_is_separate():
    """The committed primary MSA is pinned; readiness never grants execution."""
    primary = PROJECT_ROOT / "configs/phase2/4men_primary.json"
    _, _, plan = validate(primary)
    assert plan["experiment_role"] == "primary_structural_evaluation"
    assert plan["status"] == "offline_validated_not_executed"
    assert plan["ready_for_preparation"]
    assert plan["execution_authorized"] is False
    assert plan["proposed_input"]["sequences"][0]["protein"]["msa"] == "inputs/protein.a3m"
    assert "--use_msa_server" not in plan["planned_command"]
    _, _, smoke = validate(CONFIG)
    assert smoke["experiment_role"] == "smoke_test"
    assert smoke["ready_for_preparation"]
    assert smoke["execution_authorized"] is False


@pytest.mark.parametrize("unset", ["checksum_and_provenance", "checksum", "provenance"])
def test_unpinned_primary_msa_is_blocked_from_strict_preparation(tmp_path, unset):
    """A synthetic copy of the primary design with its MSA pin removed."""
    pending = read_json(PROJECT_ROOT / "configs/phase2/4men_primary.json")
    pending["experiment_id"] = "4men_msa_pending_fixture"
    if unset != "provenance":
        pending["msa"]["file"]["sha256"] = None
    if unset != "checksum":
        pending["msa"]["provenance"] = dict.fromkeys(pending["msa"]["provenance"])
    path = tmp_path / "pending.json"
    path.write_text(json.dumps(pending))
    _, _, plan = validate(path, allow_pending_msa=True)
    assert plan["status"] == "design_validated_MSA_pending_not_executable"
    assert not plan["ready_for_preparation"]
    assert plan["execution_authorized"] is False
    with pytest.raises(ValueError, match="MSA|SHA-256"):
        validate(path)
    output = tmp_path / "experiments"
    with pytest.raises(ValueError, match="MSA|SHA-256"):
        prepare(path, output)
    assert not output.exists()


@pytest.fixture
def local_msa_fixture(tmp_path, monkeypatch):
    """An explicit synthetic alignment overlay; never a production MSA."""
    primary = read_json(PROJECT_ROOT / "configs/phase2/4men_primary.json")
    primary["experiment_id"] = "4men_msa_synthetic_fixture"
    relative = primary["msa"]["file"]["path"]
    msa = tmp_path / relative
    msa.parent.mkdir(parents=True)
    query = primary["protein"]["sequence"]
    altered = query[:2] + "Q" + query[3:]
    msa.write_bytes(f">query\n{query}\n>synthetic_test_only\n{altered}\n".encode())
    primary["msa"]["file"]["sha256"] = sha256(msa)
    primary["msa"]["provenance"] = {
        "generator": "unit test", "database": "synthetic fixture only",
        "database_release": "not applicable", "retrieved_utc": "2026-10-15T00:00:00Z",
        "server": "none; synthetic fixture",
    }

    def lookup(spec, root=PROJECT_ROOT, **kwargs):
        return source_path(spec, tmp_path if spec["path"] == relative else root, **kwargs)

    monkeypatch.setattr("brd4_analysis.phase2.config.source_path", lookup)
    monkeypatch.setattr("brd4_analysis.phase2.workflow.source_path", lookup)
    return primary, msa


def test_cached_msa_is_copied_hashed_and_never_refreshed(tmp_path, local_msa_fixture, monkeypatch):
    primary, msa = local_msa_fixture
    config_path = tmp_path / "primary.json"
    config_path.write_text(json.dumps(primary), encoding="utf-8")

    def forbidden(*args, **kwargs):
        raise AssertionError("Local-MSA preparation attempted network or process execution")

    monkeypatch.setattr("requests.sessions.Session.request", forbidden)
    monkeypatch.setattr("subprocess.run", forbidden)
    run = prepare(config_path, tmp_path / "runs")
    assert read_json(run / "preregistration.json")["experiment_role"] == "primary_structural_evaluation"
    assert (run / "inputs/protein.a3m").read_bytes() == msa.read_bytes()
    assert "--use_msa_server" not in read_json(run / "command.json")["argv"]
    assert verify_run(run)["msa"]["file"]["sha256"] == sha256(msa)
    (run / "inputs/protein.a3m").write_bytes(b"changed alignment")
    with pytest.raises(ValueError, match="modified"):
        verify_run(run)


def test_query_only_alignment_cannot_masquerade_as_primary(local_msa_fixture):
    primary, msa = local_msa_fixture
    msa.write_bytes(f">query\n{primary['protein']['sequence']}\n".encode())
    primary["msa"]["file"]["sha256"] = sha256(msa)
    with pytest.raises(ValueError, match="at least one alignment hit"):
        validate_config(primary)


@pytest.mark.parametrize("padding", ["-", "X", "insertion-"])
def test_non_hit_rows_cannot_make_primary_preparable(tmp_path, local_msa_fixture, padding):
    primary, msa = local_msa_fixture
    query = primary["protein"]["sequence"]
    msa.write_bytes(f">query\n{query}\n>synthetic_padding_only\n{padding * len(query)}\n".encode())
    primary["msa"]["file"]["sha256"] = sha256(msa)
    config_path = tmp_path / "padding.json"
    config_path.write_text(json.dumps(primary), encoding="utf-8")
    with pytest.raises(ValueError, match="at least one alignment hit"):
        prepare(config_path, tmp_path / "not_prepared")
    assert not (tmp_path / "not_prepared").exists()


@pytest.mark.parametrize("row", ["AC", "ACDE", "AC1", "A.D"])
def test_malformed_alignment_rows_remain_rejected(tmp_path, row):
    msa = tmp_path / "malformed.a3m"
    msa.write_bytes(f">query\nACD\n>synthetic\n{row}\n".encode())
    with pytest.raises(ValueError, match="A3M"):
        validate_a3m(msa, "ACD")


def test_informative_alignment_preserves_lowercase_insertions_and_padding(tmp_path):
    msa = tmp_path / "valid.a3m"
    msa.write_bytes(b">query\nACD\n>padding\n---\n>synthetic_informative\nAcd-D\n")
    assert validate_a3m(msa, "ACD") == 3
    with pytest.raises(ValueError, match="at least one alignment hit"):
        validate_a3m(msa, "ACD", max_records=2)


def test_plan_mode_does_not_ignore_a_changed_pinned_alignment(local_msa_fixture):
    primary, msa = local_msa_fixture
    msa.write_bytes(msa.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="changed source"):
        validate_config(primary, allow_pending_msa=True)
