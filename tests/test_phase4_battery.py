"""Phase 4 adversarial battery tests.

Pins the in-memory evaluator to the committed Phase 3 metrics, reproduces
the auditor's 12%-expansion false acceptance through both the replica and
the registered ``score_prediction`` path, and checks every battery case's
expected verdict so the pass/fail matrix is itself regression-tested.
"""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from brd4_analysis.constants import PROJECT_ROOT
from brd4_analysis.phase2.config import read_json
from brd4_analysis.phase2.workflow import score_prediction
from brd4_analysis.phase4 import battery as battery_module
from brd4_analysis.phase4 import cli, perturb
from brd4_analysis.phase4.battery import CASES, run_battery
from brd4_analysis.phase4.evaluate import (
    ANALYSIS_METRICS,
    PREDICTION_CIF,
    aligned_ligand_ref_order,
    alignment,
    evaluate,
    load_context,
    replace_all_coordinates,
)
from brd4_analysis.phase4.synthetic import prediction_protein_atoms, write_prediction


@pytest.fixture(scope="module")
def ctx():
    return load_context()


@pytest.fixture(scope="module")
def battery(ctx):
    return run_battery(ctx)


@pytest.fixture(scope="module")
def case_map(battery):
    return {c["name"]: c for c in battery["cases"]}


def write_faithful_prediction(path, ctx, pred):
    """Serialize a prediction dictionary without replacing its protein by
    deposited coordinates. Atom names come from the source mmCIF because
    read_prediction() exposes residue coordinate arrays."""
    rotation, translation, _ = alignment(pred, ctx["ref"])
    atoms = {
        label: {
            name: coordinates @ rotation + translation
            for name, coordinates in residue.items()
        }
        for label, residue in prediction_protein_atoms(
            PREDICTION_CIF, ctx["ref"], pred).items()
    }
    return write_prediction(
        path, ctx["ref"], aligned_ligand_ref_order(pred, ctx["ref"]),
        protein_atoms=atoms,
    )


def test_parity_with_committed_phase3_metrics(ctx):
    """evaluate() on the real prediction must reproduce the frozen run's
    committed metrics within the tight scientific tolerance used for portable
    artifact comparison — this is what licenses the battery to stand in for
    score_prediction on synthetic variants."""
    committed = read_json(ANALYSIS_METRICS)
    result = evaluate(ctx["config"], ctx["ref"], ctx["pred"])
    tolerance = {"rel": 1e-12, "abs": 1e-12}
    assert result["ligand_pose"]["rmsd_A"] == pytest.approx(
        committed["ligand_pose"]["rmsd_A"], **tolerance)
    assert result["protein_ca_rmsd_A"] == pytest.approx(
        committed["protein_ca_rmsd_A"], **tolerance)
    assert result["acceptance_gates"] == committed["acceptance_gates"]
    assert result["harness_success"] == committed["harness_success"]
    for cutoff in ("3.5", "4.0", "5.0"):
        got, exp = result["contacts"][cutoff], committed["contacts"][cutoff]
        assert got["predicted_residues"] == exp["predicted_residues"]
        assert got["precision"] == pytest.approx(exp["precision"], **tolerance)
        assert got["recall"] == pytest.approx(exp["recall"], **tolerance)
        assert (got["true_positive"], got["false_positive"], got["false_negative"]) == (
            exp["true_positive"], exp["false_positive"], exp["false_negative"])
    for key in ("max_bond_length_error_A", "max_bond_angle_error_deg",
                "max_aromatic_plane_deviation_A", "min_nonbonded_vdw_ratio",
                "severe_protein_ligand_overlap_pairs", "valid"):
        expected = committed["geometry"][key]
        if isinstance(expected, float):
            assert result["geometry"][key] == pytest.approx(expected, **tolerance)
        else:
            assert result["geometry"][key] == expected


def test_parity_with_registered_score_prediction(ctx):
    """The replica and the real scorer agree field-for-field on the same
    committed prediction file."""
    registered = score_prediction(ctx["config"], PREDICTION_CIF)
    replica = evaluate(ctx["config"], ctx["ref"], ctx["pred"])
    for key in ("ligand_pose", "protein_ca_rmsd_A", "pocket_ca_rmsd_in_global_frame_A",
                "contacts", "contacts_reference_observed_atoms", "geometry",
                "acceptance_gates", "harness_success"):
        assert replica[key] == registered[key]


def test_registered_gates_pass_on_real_prediction(case_map):
    case = case_map["control_real_prediction"]
    assert case["harness_success"]
    assert case["metrics"]["ligand_pose_rmsd_A"] == pytest.approx(0.1404848997959765)
    assert case["metrics"]["protein_ca_rmsd_A"] == pytest.approx(0.19136727025432587)
    assert case["metrics"]["contacts"]["4.0"]["precision"] == 1.0
    assert case["metrics"]["contacts"]["4.0"]["recall"] == 1.0
    assert case["metrics"]["contacts"]["3.5"]["precision"] == pytest.approx(2 / 3)
    assert case["metrics"]["contacts"]["3.5"]["recall"] == 0.4
    assert case["metrics"]["contacts"]["5.0"]["precision"] == 1.0


def test_the_12_percent_expansion_false_acceptance(case_map):
    """The prior audit's finding, reproduced: a uniform 12% inflation of the
    committed prediction's ligand passes every registered gate."""
    case = case_map["centroid_scale_1.12"]
    assert case["harness_success"]  # false acceptance
    assert all(case["gates"].values())
    assert case["metrics"]["ligand_pose_rmsd_A"] == pytest.approx(0.3972135, abs=1e-5)
    # Internal evidence that the accepted structure is chemically impossible:
    # every bond is ~12% too long, yet inside the absolute 0.25 A gate.
    assert case["metrics"]["geometry"]["max_bond_length_error_A"] == pytest.approx(0.205, abs=0.01)
    assert case["metrics"]["geometry"]["max_bond_angle_error_deg"] < 4.0
    # And the proposed checks do catch it.
    for proposal in ("bond_signed_bias", "bond_relative_error",
                     "similarity_scale", "internal_distance_rmsd"):
        assert not case["candidate_checks"][proposal]["pass"]


def test_end_to_end_variants_through_registered_scorer(ctx, tmp_path):
    """Serialize faithful predicted complexes and score them with the
    unmodified registered harness, including the audited ligand expansion."""
    case_by_name = {case["name"]: case for case in CASES}

    expanded_pred = case_by_name["centroid_scale_1.12"]["apply"](ctx)
    path = write_faithful_prediction(
        tmp_path / "faithful_expanded_112.cif", ctx, expanded_pred)
    result = score_prediction(ctx["config"], path)
    assert result["harness_success"]  # accepted synthetic variant
    assert result["ligand_pose"]["rmsd_A"] == pytest.approx(0.3972135, abs=1e-6)
    assert result["protein_ca_rmsd_A"] == pytest.approx(0.1913673, abs=1e-6)
    assert result["geometry"]["valid"]

    opened = case_by_name["angle_bend_C5-N-C7_+20deg"]["apply"](ctx)
    opened_result = score_prediction(
        ctx["config"],
        write_faithful_prediction(tmp_path / "angle_opening_20.cif", ctx, opened),
    )
    assert opened_result["harness_success"]
    assert opened_result["ligand_pose"]["rmsd_A"] == pytest.approx(0.5525381, abs=1e-6)
    assert opened_result["geometry"]["max_bond_angle_error_deg"] == pytest.approx(
        23.2018766, abs=1e-6)

    expanded_protein = case_by_name["protein_scale_1.25"]["apply"](ctx)
    protein_result = score_prediction(
        ctx["config"],
        write_faithful_prediction(
            tmp_path / "protein_scale_125.cif", ctx, expanded_protein),
    )
    assert not protein_result["harness_success"]
    assert protein_result["acceptance_gates"]["protein_ca_rmsd"] is False
    assert protein_result["protein_ca_rmsd_A"] == pytest.approx(3.6284726, abs=1e-6)

    deposited_expanded = write_prediction(
        tmp_path / "expanded_deposited.cif", ctx["ref"],
        perturb.centroid_scale(ctx["ref"]["ligand_xyz"], 1.12),
    )
    dep_result = score_prediction(ctx["config"], deposited_expanded)
    assert not dep_result["harness_success"]  # deposited bonds already deviate
    assert "bond_lengths" in dep_result["geometry"]["failed_checks"]


def test_symmetric_atom_relabelling_is_pose_invariant(case_map):
    """Heavy-atom correspondence is symmetry-aware: permuting the ligand
    through its non-identity graph automorphism must not change the pose
    RMSD the harness reports."""
    real = case_map["control_real_prediction"]
    swap = case_map["control_symmetry_swap"]
    assert swap["harness_success"]
    assert swap["metrics"]["ligand_pose_rmsd_A"] == pytest.approx(
        real["metrics"]["ligand_pose_rmsd_A"])


def test_controls_pass_gates_and_all_proposals(case_map):
    for name in ("control_real_prediction", "control_reference_replay",
                 "control_complex_rigid", "control_symmetry_swap"):
        case = case_map[name]
        assert case["harness_success"], name
        assert all(c["pass"] for c in case["candidate_checks"].values()), name


def test_rigid_motion_is_invariant_internal_distortion_is_not(case_map):
    """Rigid whole-complex motion keeps every metric; internal ligand
    distortion keeps none of the internal statistics."""
    rigid = case_map["control_complex_rigid"]
    real = case_map["control_real_prediction"]
    assert rigid["metrics"]["ligand_pose_rmsd_A"] == pytest.approx(
        real["metrics"]["ligand_pose_rmsd_A"], abs=1e-9)
    distorted = case_map["centroid_scale_1.12"]
    assert distorted["metrics"]["geometry"]["max_bond_length_error_A"] > 3 * \
        real["metrics"]["geometry"]["max_bond_length_error_A"]


# --- the post-audit pass/fail matrix ----------------------------------------
# Every non-control case is pinned to an expected verdict so a future change
# in gates, thresholds or perturbation code is loud, not silent. Acceptance
# is not treated as an empirical false-positive label.

EXPECTED_ACCEPTED_VARIANTS = {
    "centroid_scale_1.02", "centroid_scale_1.05", "centroid_scale_1.08",
    "centroid_scale_1.10", "centroid_scale_1.12", "centroid_scale_1.15",
    "ligand_translate_0.25A", "ligand_translate_0.50A", "ligand_translate_1.00A",
    "ligand_rotate_15deg", "bond_stretch_C5-N_+0.20A",
    "angle_bend_C5-N-C7_+20deg",
    "protein_scale_1.05", "jitter_sigma_0.05",
}

EXPECTED_REJECTION = {
    "centroid_scale_0.85", "centroid_scale_0.90", "centroid_scale_1.20",
    "centroid_scale_1.30",
    "ligand_translate_1.50A", "ligand_translate_2.00A",
    "ligand_translate_3.00A", "ligand_translate_4.00A",
    "clash_translate_1.0A", "clash_translate_2.0A", "clash_translate_3.0A",
    "ligand_rotate_45deg", "ligand_rotate_90deg", "ligand_rotate_180deg",
    "torsion_C5-N_90deg", "torsion_C5-N_180deg",
    "torsion_N-C7_90deg", "torsion_N-C7_180deg",
    "bond_stretch_C5-N_-0.30A", "bond_stretch_C5-N_+0.30A",
    "angle_bend_C5-N-C7_-35deg", "angle_bend_C5-N-C7_-20deg",
    "angle_bend_C5-N-C7_+35deg",
    "reflect_ligand", "anisotropic_scale_z_1.25", "collapse_0.70",
    "protein_scale_1.25", "jitter_sigma_0.20", "jitter_sigma_0.50",
}


def test_battery_covers_every_case_exactly_once(case_map):
    names = set(case_map)
    controls = {n for n in names if n.startswith("control_")}
    assert controls == {"control_real_prediction", "control_reference_replay",
                        "control_complex_rigid", "control_symmetry_swap"}
    assert names == controls | EXPECTED_ACCEPTED_VARIANTS | EXPECTED_REJECTION
    assert len(names) == len(CASES)


def test_accepted_variant_set_is_exact(case_map):
    passing = {
        name for name, case in case_map.items()
        if case["harness_success"] and not name.startswith("control_")
    }
    assert passing == EXPECTED_ACCEPTED_VARIANTS


def test_deliberately_defective_cases_fail_specific_gates(case_map):
    """Controls on the failure side: defects trip the gate designed for them."""
    assert case_map["collapse_0.70"]["gates"]["geometry"] is False
    assert "intraligand_clash" in \
        case_map["collapse_0.70"]["metrics"]["geometry"]["failed_checks"]
    assert case_map["bond_stretch_C5-N_+0.30A"]["gates"]["geometry"] is False
    assert "bond_lengths" in \
        case_map["bond_stretch_C5-N_+0.30A"]["metrics"]["geometry"]["failed_checks"]
    for clash in ("clash_translate_1.0A", "clash_translate_2.0A", "clash_translate_3.0A"):
        assert case_map[clash]["gates"]["geometry"] is False
        assert "protein_ligand_clash" in \
            case_map[clash]["metrics"]["geometry"]["failed_checks"]
    assert case_map["ligand_translate_4.00A"]["gates"]["ligand_rmsd"] is False
    assert case_map["ligand_translate_4.00A"]["gates"]["contact_precision"] is False
    assert case_map["ligand_translate_4.00A"]["gates"]["contact_recall"] is False


def test_reflection_is_rejected_only_via_pose_not_internal_geometry(case_map):
    """The mirror image keeps perfect internal geometry; it is caught by the
    pose/contact consequences, not by any chirality-aware check — 25K is
    achiral, so nothing internal *can* see the flip."""
    case = case_map["reflect_ligand"]
    assert not case["harness_success"]
    assert case["gates"]["ligand_rmsd"] is False
    internal = case["metrics"]["geometry"]
    assert "bond_lengths" not in internal["failed_checks"]
    assert "bond_angles" not in internal["failed_checks"]
    assert "aromatic_planarity" not in internal["failed_checks"]


def test_protein_side_inflation_is_accepted_by_registered_gates(case_map):
    """5% protein inflation: CA gate granularity documented — the protein
    fit absorbs a correlated 5% coordinate error."""
    case = case_map["protein_scale_1.05"]
    assert case["harness_success"]
    assert 0.5 < case["metrics"]["protein_ca_rmsd_A"] < 2.0


def test_uncorrelated_noise_vs_correlated_drift(case_map):
    """The designed distinction: sigma=0.05 jitter passes the *bias* check
    (noise averages out) while 5% uniform inflation fails it (correlated)."""
    jitter = case_map["jitter_sigma_0.05"]
    scale = case_map["centroid_scale_1.05"]
    assert jitter["candidate_checks"]["bond_signed_bias"]["pass"]
    assert not scale["candidate_checks"]["bond_signed_bias"]["pass"]


def test_contacts_evaluated_at_all_three_cutoffs(battery):
    for case in battery["cases"]:
        assert set(case["metrics"]["contacts"]) == {"3.5", "4.0", "5.0"}


def test_battery_is_deterministic(ctx):
    a = json.dumps(run_battery(ctx), sort_keys=True)
    b = json.dumps(run_battery(ctx), sort_keys=True)
    assert a == b


def test_no_registered_threshold_was_altered(ctx, battery):
    """The battery must report the *frozen* thresholds, not Phase 4 edits."""
    committed = read_json(
        PROJECT_ROOT / "experiments" / "phase2" / "4men_msa_s20261019_v1" / "config.json")
    assert battery["registered_thresholds"] == committed["metrics"] == ctx["config"]["metrics"]
    baseline = read_json(PROJECT_ROOT / "configs" / "phase2" / "4men_primary.json")
    assert ctx["config"]["metrics"] == baseline["metrics"]


def _assert_scientific_artifact_equal(generated: dict, stored: dict) -> None:
    """Compare deterministic scientific fields while retaining truthful
    environment provenance. Floating-point fields use a tight numerical
    tolerance so ordinary platform-dependent BLAS/library rounding cannot
    produce a false scientific regression."""
    assert set(generated) == set(stored)
    for artifact in (generated, stored):
        assert set(artifact["environment"]) == {"python", "platform", "packages"}
        assert set(artifact["environment"]["packages"]) == {"biopython", "numpy", "rdkit"}

    def assert_node(actual, expected, path="$"):
        if isinstance(expected, float):
            assert actual == pytest.approx(expected, rel=1e-12, abs=1e-12), path
        elif isinstance(expected, dict):
            assert set(actual) == set(expected), path
            for key, value in expected.items():
                assert_node(actual[key], value, f"{path}.{key}")
        elif isinstance(expected, list):
            assert len(actual) == len(expected), path
            for index, value in enumerate(expected):
                assert_node(actual[index], value, f"{path}[{index}]")
        else:
            assert actual == expected, path

    for key in stored:
        if key != "environment":
            assert_node(generated[key], stored[key], f"$.{key}")


def test_battery_records_post_audit_provenance(battery):
    assert battery["schema_version"] == 2
    assert battery["battery_revision"] == "post_audit_v1"
    assert "post-audit" in battery["design_status"]
    assert battery["random_seeds"] == {"gaussian_jitter": 20261018}
    assert battery["case_parameters"]["jitter_sigma_0.05"] == {
        "sigma_A": 0.05, "seed": 20261018}
    assert set(battery["phase4_source_sha256_LF"]) == {
        "__init__.py", "battery.py", "checks.py", "cli.py",
        "evaluate.py", "perturb.py", "synthetic.py"}
    assert set(battery["environment"]["packages"]) == {"biopython", "numpy", "rdkit"}


def _markdown_rows(document: str) -> list[dict]:
    rows = []
    header = None
    for line in document.splitlines():
        if not line.startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if cells[0] == "case":
            header = cells
        elif header is not None and cells[0] != "---":
            rows.append(dict(zip(header, cells)))
    return rows


def _assert_findings_document_matches_result(document: str, result: dict) -> None:
    rows = _markdown_rows(document)
    gate_rows = [row for row in rows if "all gates" in row]
    value_rows = [row for row in rows if "all gates" not in row]
    expected_names = [f"`{case['name']}`" for case in result["cases"]]
    assert [row["case"] for row in gate_rows] == expected_names
    assert [row["case"] for row in value_rows] == expected_names

    for row, case in zip(gate_rows, result["cases"], strict=True):
        assert row["family"] == case["family"]
        for gate, passed in case["gates"].items():
            assert row[gate] == ("PASS" if passed else "FAIL")
        assert row["all gates"] == ("PASS" if case["harness_success"] else "FAIL")
        proposals = set() if row["failed proposals"] == "-" else set(
            row["failed proposals"].split(", "))
        assert proposals == {
            name for name, check in case["candidate_checks"].items()
            if not check["pass"]
        }

    for row, case in zip(value_rows, result["cases"], strict=True):
        assert set(row) - {"case"} == set(case["candidate_checks"])
        for name, check in case["candidate_checks"].items():
            expected = f"{check['value']:.4g}" + (
                "" if check["pass"] else " **FAIL**")
            assert row[name] == expected


def test_generated_result_matches_versioned_artifact(battery):
    artifact_dir = PROJECT_ROOT / "experiments" / "phase4" / "adversarial_screen_v1"
    stored = read_json(artifact_dir / "results.json")
    _assert_scientific_artifact_equal(battery, stored)
    _assert_findings_document_matches_result(
        cli.findings_markdown(battery), battery)
    _assert_findings_document_matches_result(
        (artifact_dir / "FINDINGS.md").read_text(encoding="utf-8"), stored)


def test_scientific_artifact_comparison_is_environment_portable(ctx, monkeypatch):
    """Windows/Linux Python micro-version and kernel strings are provenance,
    not scientific output."""
    monkeypatch.setattr(battery_module.sys, "version", "3.12.10 |test build|")
    monkeypatch.setattr(battery_module.platform, "platform",
                        lambda: "Windows-2022Server-x86_64-test")
    generated = run_battery(ctx)
    assert generated["environment"]["python"].startswith("3.12.10")
    assert generated["environment"]["platform"].startswith("Windows-2022")
    stored = read_json(
        PROJECT_ROOT / "experiments" / "phase4" / "adversarial_screen_v1" / "results.json")
    _assert_scientific_artifact_equal(generated, stored)


def test_generated_artifact_manifest_hashes():
    artifact_dir = PROJECT_ROOT / "experiments" / "phase4" / "adversarial_screen_v1"
    manifest = read_json(artifact_dir / "manifest.json")
    assert manifest["battery_revision"] == "post_audit_v1"
    for name, digest in manifest["output_artifacts_sha256"].items():
        from brd4_analysis.phase2.config import sha256
        assert sha256(artifact_dir / name) == digest


def test_rotation_contracts_reject_nonorthonormal_scale(ctx, tmp_path):
    matrix = np.diag([2.0, 0.5, 1.0])
    with pytest.raises(ValueError, match="orthonormal"):
        replace_all_coordinates(ctx["pred"], matrix, np.zeros(3))
    with pytest.raises(ValueError, match="orthonormal"):
        write_prediction(
            tmp_path / "bad_rotation.cif", ctx["ref"], ctx["ref"]["ligand_xyz"],
            rotation=matrix)


def test_protein_ca_gate_has_a_failing_control(case_map):
    case = case_map["protein_scale_1.25"]
    assert case["gates"]["protein_ca_rmsd"] is False
    assert case["metrics"]["protein_ca_rmsd_A"] == pytest.approx(3.6284726)


def test_recommended_exploratory_combination_has_documented_misses(case_map):
    recommended = (
        "similarity_scale", "internal_distance_rmsd",
        "bond_signed_bias", "pose_rmsd_tightened")
    for name in ("bond_stretch_C5-N_+0.20A", "protein_scale_1.05"):
        assert case_map[name]["harness_success"]
        assert all(case_map[name]["candidate_checks"][check]["pass"]
                   for check in recommended)


def _fake_battery_result(success=True):
    return {
        "battery_revision": "post_audit_v1",
        "design_status": "post-audit test result",
        "provenance": {},
        "phase4_source_sha256_LF": {},
        "cases": [{
            "name": "control_real_prediction",
            "family": "control",
            "harness_success": success,
            "gates": {"ligand_rmsd": success},
            "metrics": {
                "ligand_pose_rmsd_A": 0.0,
                "contacts": {"4.0": {"precision": 0.0, "recall": 0.0}},
            },
            "candidate_checks": {},
        }],
    }


def _symlink(link, target, *, target_is_directory=False):
    try:
        link.symlink_to(target, target_is_directory=target_is_directory)
    except (NotImplementedError, OSError) as exc:
        pytest.skip(f"symlinks unavailable in this environment: {exc}")


def _directory_symlink(link, target):
    _symlink(link, target, target_is_directory=True)


def test_battery_output_path_is_constrained(tmp_path, monkeypatch):
    allowed = tmp_path / "experiments" / "phase4"
    monkeypatch.setattr(cli, "PHASE4_OUTPUT_ROOT", allowed.resolve())
    assert cli.battery_output_dir(allowed / "screen") == allowed.resolve() / "screen"
    with pytest.raises(ValueError, match="child of"):
        cli.battery_output_dir(tmp_path / "experiments" / "phase2" / "run")
    link_target = tmp_path / "outside"
    link_target.mkdir()
    (allowed / "linked").mkdir(parents=True)
    (allowed / "linked").rmdir()
    _directory_symlink(allowed / "linked", link_target)
    with pytest.raises(ValueError, match="symlink"):
        cli.battery_output_dir(allowed / "linked")


def test_cli_refuses_symlinked_artifact(tmp_path, monkeypatch):
    out_dir = tmp_path / "battery"
    out_dir.mkdir()
    target = tmp_path / "outside_results.json"
    target.write_text("sentinel", encoding="utf-8")
    _symlink(out_dir / "results.json", target)
    monkeypatch.setattr(cli, "PHASE4_OUTPUT_ROOT", tmp_path.resolve())
    with pytest.raises(ValueError, match="symlinked artifact"):
        cli.cmd_battery(out_dir)
    assert target.read_text(encoding="utf-8") == "sentinel"


@pytest.mark.parametrize("artifact_name", cli.ARTIFACT_NAMES)
def test_cli_hard_linked_artifact_replaces_entry_not_linked_inode(
    tmp_path, monkeypatch, artifact_name
):
    out_dir = tmp_path / "battery"
    out_dir.mkdir()
    sentinel = tmp_path / f"outside_{artifact_name}"
    sentinel.write_text("sentinel", encoding="utf-8")
    os.link(sentinel, out_dir / artifact_name)
    monkeypatch.setattr(cli, "load_context", lambda: {})
    monkeypatch.setattr(cli, "run_battery", lambda ctx: _fake_battery_result())
    monkeypatch.setattr(cli, "PHASE4_OUTPUT_ROOT", tmp_path.resolve())

    assert cli.cmd_battery(out_dir) == 0
    assert sentinel.read_text(encoding="utf-8") == "sentinel"
    assert (out_dir / artifact_name).read_text(encoding="utf-8") != "sentinel"
    assert (out_dir / artifact_name).stat().st_nlink == 1


def test_cli_rejects_symlinked_phase4_output_root(tmp_path, monkeypatch):
    real_root = tmp_path / "real_phase4"
    real_root.mkdir()
    linked_root = tmp_path / "phase4"
    _directory_symlink(linked_root, real_root)
    monkeypatch.setattr(cli, "PHASE4_OUTPUT_ROOT", linked_root)
    with pytest.raises(ValueError, match="output root.*symlink"):
        cli.battery_output_dir(linked_root / "battery")


def test_cli_public_path_uses_descriptor_writer_on_posix(tmp_path, monkeypatch):
    if os.name != "posix":
        pytest.skip("descriptor-relative writer selection is POSIX-specific")
    assert cli._POSIX_DIRFD

    out_dir = tmp_path / "phase4" / "battery"
    out_dir.mkdir(parents=True)
    opened, renamed = [], []
    real_open, real_rename = cli.os.open, cli.os.rename

    def open_spy(path, flags, mode=0o777, *, dir_fd=None):
        fd = real_open(path, flags, mode, dir_fd=dir_fd)
        if dir_fd is not None:
            opened.append((str(path), flags, dir_fd))
        return fd

    def rename_spy(src, dst, *, src_dir_fd=None, dst_dir_fd=None):
        renamed.append((str(src), str(dst), src_dir_fd, dst_dir_fd))
        return real_rename(
            src, dst, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd)

    monkeypatch.setattr(cli.os, "open", open_spy)
    monkeypatch.setattr(cli.os, "rename", rename_spy)
    monkeypatch.setattr(
        cli.tempfile, "NamedTemporaryFile",
        lambda *args, **kwargs: pytest.fail("portable writer used"))
    monkeypatch.setattr(cli, "load_context", lambda: {})
    monkeypatch.setattr(cli, "run_battery", lambda ctx: _fake_battery_result())
    monkeypatch.setattr(cli, "PHASE4_OUTPUT_ROOT", tmp_path.resolve() / "phase4")

    assert cli.cmd_battery(out_dir) == 0
    temporary_opens = [
        op for op in opened
        if op[1] & cli.os.O_CREAT and Path(op[0]).name.startswith(".")
    ]
    artifact_renames = [op for op in renamed if op[1] in cli.ARTIFACT_NAMES]
    descriptor_ids = {op[2] for op in temporary_opens}
    assert len(temporary_opens) == len(cli.ARTIFACT_NAMES)
    assert {op[1] for op in artifact_renames} == set(cli.ARTIFACT_NAMES)
    assert len(descriptor_ids) == 1
    descriptor = descriptor_ids.pop()
    assert all(
        op[2] == descriptor and op[3] == descriptor for op in artifact_renames)


@pytest.mark.parametrize("artifact_name", cli.ARTIFACT_NAMES)
def test_cli_pinned_output_directory_does_not_follow_pathname_swap(
    tmp_path, monkeypatch, artifact_name
):
    out_dir = tmp_path / "phase4" / "battery"
    out_dir.mkdir(parents=True)
    external = tmp_path / "external"
    external.mkdir()
    sentinel = external / artifact_name
    sentinel.write_text("sentinel", encoding="utf-8")

    def replace_destination(ctx):
        out_dir.rename(out_dir.with_name("kept"))
        _directory_symlink(out_dir, external)
        return _fake_battery_result()

    monkeypatch.setattr(cli, "load_context", lambda: {})
    monkeypatch.setattr(cli, "run_battery", replace_destination)
    monkeypatch.setattr(
        cli, "PHASE4_OUTPUT_ROOT", tmp_path.resolve() / "phase4")
    if cli._POSIX_DIRFD:
        assert cli.cmd_battery(out_dir) == 0
        for name in cli.ARTIFACT_NAMES:
            assert (out_dir.with_name("kept") / name).is_file()
    else:
        with pytest.raises(ValueError, match="unsafe"):
            cli.cmd_battery(out_dir)
    assert sentinel.read_text(encoding="utf-8") == "sentinel"
    for name in cli.ARTIFACT_NAMES:
        if name != artifact_name:
            assert not (external / name).exists()


def test_cli_safely_fails_if_pinned_output_directory_is_removed(
    tmp_path, monkeypatch
):
    out_dir = tmp_path / "battery"
    out_dir.mkdir()
    external = tmp_path / "external"
    external.mkdir()
    (external / "sentinel").write_text("sentinel", encoding="utf-8")

    def remove_destination(ctx):
        out_dir.rmdir()
        _directory_symlink(out_dir, external)
        return _fake_battery_result()

    monkeypatch.setattr(cli, "load_context", lambda: {})
    monkeypatch.setattr(cli, "run_battery", remove_destination)
    monkeypatch.setattr(cli, "PHASE4_OUTPUT_ROOT", tmp_path.resolve())
    with pytest.raises((FileNotFoundError, ValueError)):
        cli.cmd_battery(out_dir)
    assert (external / "sentinel").read_text(encoding="utf-8") == "sentinel"
    for name in cli.ARTIFACT_NAMES:
        assert not (external / name).exists()


def test_cli_returns_nonzero_when_required_control_fails(tmp_path, monkeypatch):
    result = _fake_battery_result(success=False)
    monkeypatch.setattr(cli, "load_context", lambda: {})
    monkeypatch.setattr(cli, "run_battery", lambda ctx: result)
    monkeypatch.setattr(cli, "PHASE4_OUTPUT_ROOT", tmp_path.resolve())
    assert cli.cmd_battery(tmp_path / "battery") == 1
    assert (tmp_path / "battery" / "results.json").is_file()


def test_cli_rejects_phase3_output_before_running(ctx, monkeypatch):
    monkeypatch.setattr(
        cli, "load_context", lambda: pytest.fail("battery ran before path validation"))
    with pytest.raises(ValueError, match="child of"):
        cli.cmd_battery(PROJECT_ROOT / "experiments" / "phase2" / "4men_msa_s20261019_v1")
