"""Prepare immutable experiment bundles and analyse one externally approved run.

There is deliberately no command to install, download weights, invoke Boltz,
submit a job, or iterate over seeds.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import io
import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np

from brd4_analysis.constants import PROJECT_ROOT
from brd4_analysis.parse import atom_element
from brd4_analysis.phase2.config import (
    SOURCE_DIR, experiment_role, msa_is_pinned, read_json, sha256, source_path,
    validate_config, write_json,
)
from brd4_analysis.phase2.inputs import atom_names, read_prediction, reference_data
from brd4_analysis.phase2.metrics import (
    contact_residues, contact_scores, geometry_checks, protein_alignment, rmsd,
    symmetry_rmsd,
)

OUTPUT_ROOT = PROJECT_ROOT / "experiments" / "phase2"


def analysis_source_hashes() -> dict[str, str]:
    """Hash new code and reused Phase 1 utilities, normalizing source EOL to LF.

    Data and prepared input hashes remain byte-exact. Source normalization avoids
    invalidating a bundle solely because Git checks out Python with Windows EOL.
    """
    package = Path(__file__).parent.parent
    paths = [package / name for name in [
        "__init__.py", "constants.py", "contacts.py", "descriptors.py", "parse.py",
    ]] + sorted((package / "phase2").glob("*.py"))
    return {
        path.relative_to(package).as_posix(): hashlib.sha256(
            path.read_text(encoding="utf-8").encode("utf-8")
        ).hexdigest() for path in paths
    }


def verify_sources() -> dict:
    manifest = read_json(SOURCE_DIR / "manifest.json")
    for item in manifest["files"]:
        path = SOURCE_DIR / item["path"]
        if path.parent != SOURCE_DIR or sha256(path) != item["sha256"]:
            raise ValueError(f"Changed source snapshot: {item['path']}")
    return manifest


def planned_command(config: dict) -> list[str]:
    b = config["boltz"]
    argv = [
        "boltz", "predict", "inputs/4men.yaml", "--out_dir", "outputs",
        "--cache", b["cache_dir"], "--checkpoint", b["cache_dir"] + "/boltz2_conf.ckpt",
        "--model", "boltz2", "--accelerator", "gpu", "--devices", "1",
        "--output_format", "mmcif", "--num_workers", "0", "--preprocessing-threads", "1",
    ]
    for name in ["seed", "recycling_steps", "sampling_steps", "diffusion_samples",
                 "max_parallel_samples", "step_scale", "max_msa_seqs"]:
        argv += ["--" + name, str(b[name])]
    argv.append("--no_kernels")
    return argv


def verify_interface(config: dict) -> dict:
    """Inspect the exact wheel's Click decorators; do not import boltz.main."""
    verify_sources()
    source = SOURCE_DIR / "wheel_boltz_main.py.txt"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    predict = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "predict")
    options = {}
    for decorator in predict.decorator_list:
        if isinstance(decorator, ast.Call) and isinstance(decorator.func, ast.Attribute) and decorator.func.attr == "option":
            name = ast.literal_eval(decorator.args[0])
            flag = next((ast.literal_eval(k.value) for k in decorator.keywords if k.arg == "is_flag"), False)
            options[name] = flag
    argv = planned_command(config)
    for arg in argv:
        if arg.startswith("--") and arg not in options:
            raise ValueError(f"Option absent from the pinned wheel: {arg}")
    pypi = read_json(SOURCE_DIR / "boltz_pypi_2.2.1.json")
    if pypi["info"]["version"] != config["boltz"]["version"]:
        raise ValueError("Wheel version mismatch")
    try:
        installed = version("boltz")
    except PackageNotFoundError:
        installed = None
    wheel = next(f for f in pypi["urls"] if f["filename"].endswith(".whl"))
    return {
        "method": "static_inspection_of_checksum_verified_PyPI_wheel",
        "version": pypi["info"]["version"], "wheel_sha256": wheel["digests"]["sha256"],
        "locally_installed_boltz": installed, "runtime_cli_executed": False,
        "predict_options_verified": [arg for arg in argv if arg.startswith("--")],
        "upstream_git_tag_commit": "cb04aeccdd480fd4db707f0bbafde538397fa2ac",
    }


def model_provenance() -> dict:
    data = read_json(SOURCE_DIR / "boltz_model_metadata.json")
    return {
        "repository": "boltz-community/boltz-2", "revision": data["sha"],
        "files": [
            {
                "filename": item["rfilename"], "sha256": item["lfs"]["sha256"],
                "size_bytes": item["size"],
                "url": f"https://huggingface.co/boltz-community/boltz-2/resolve/{data['sha']}/{item['rfilename']}",
            }
            for item in data["siblings"] if "lfs" in item
        ],
        "downloaded": False,
        "note": "v2.2.1 downloads mols.tar AND both checkpoints, even without an affinity request",
    }


def reference_arrays(ref):
    by_residue = {
        label: np.asarray([a.coord for a in atoms.values()], dtype=float)
        for label, atoms in ref["protein_by_residue"].items()
    }
    atoms = [a for residue in ref["protein_by_residue"].values() for a in residue.values()]
    return by_residue, np.asarray([a.coord for a in atoms], dtype=float), [atom_element(a) for a in atoms]


def input_document(config: dict) -> dict:
    msa_value = "empty" if config["msa"]["mode"] == "single_sequence" else "inputs/protein.a3m"
    return {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": config["protein"]["sequence"], "msa": msa_value}},
            {"ligand": {"id": "B", "smiles": config["ligand"]["smiles"]}},
        ],
    }


def validate(config_path: Path, *, allow_pending_msa: bool = False) -> tuple[dict, dict, dict]:
    config = validate_config(read_json(config_path), allow_pending_msa=allow_pending_msa)
    ref = reference_data(config)
    by_residue, xyz, elements = reference_arrays(ref)
    report = {
        "status": "offline_validated_not_executed" if msa_is_pinned(config) else "design_validated_MSA_pending_not_executable",
        "experiment_role": experiment_role(config),
        "ready_for_preparation": msa_is_pinned(config),
        "execution_authorized": False,
        "msa": config["msa"],
        "proposed_input": input_document(config),
        "planned_command": planned_command(config),
        "interface": verify_interface(config),
        "sequence_length": len(ref["mapping"]), "fit_CA_count": 125,
        "construct": "SM expression remnant + UniProt O60885 residues 44-168",
        "ligand": ref["descriptors"],
        "reference_contacts": {
            str(cutoff): sorted(contact_residues(by_residue, ref["ligand_xyz"], cutoff))
            for cutoff in config["metrics"]["contact_cutoffs_A"]
        },
        "reference_geometry_screen": geometry_checks(
            ref["mol"], ref["ligand_xyz"], ref["ideal_xyz"],
            xyz, elements, config["metrics"]["geometry"],
        ),
        "model": model_provenance(),
    }
    return config, ref, report


def prepare(config_path: Path, output_root: Path = OUTPUT_ROOT) -> Path:
    config, ref, report = validate(config_path)
    output_root = Path(output_root).resolve()
    protected = ["data", "results", "docs", "src", "tests"]
    if output_root == PROJECT_ROOT or any(
        output_root.is_relative_to(PROJECT_ROOT / name) for name in protected
    ):
        raise ValueError("Experiment output must not overwrite a Phase 1 or source directory")
    run = output_root / config["experiment_id"]
    run.mkdir(parents=True, exist_ok=False)
    inputs = run / "inputs"
    inputs.mkdir()
    if config["msa"]["mode"] == "local_a3m":
        shutil.copyfile(source_path(config["msa"]["file"]), inputs / "protein.a3m")
    # JSON is a YAML subset understood by Boltz's yaml.safe_load, avoiding a new dependency.
    write_json(inputs / "4men.yaml", input_document(config))
    write_json(run / "config.json", config)
    write_json(run / "validation.json", report)
    write_json(run / "command.json", {"argv": planned_command(config), "working_directory": "this experiment directory", "executed": False})
    write_json(run / "ligand_atom_names.json", {
        "rdkit_version": version("rdkit"),
        "smiles_input_order": atom_names(ref["predicted_mol"]),
        "ccd_reference_order": atom_names(ref["mol"]),
    })
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(ref["mapping"][0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(ref["mapping"])
    (run / "residue_mapping.csv").write_bytes(buffer.getvalue().encode())
    write_json(run / "preregistration.json", {
        "status": "awaiting_user_review_no_prediction_generated",
        "experiment_role": experiment_role(config),
        "experiment_id": config["experiment_id"], "metrics": config["metrics"],
        "fit_residues": [m["reference_label"] for m in ref["mapping"] if m["fit_CA"]],
        "pocket_residues": sorted(contact_residues(
            reference_arrays(ref)[0], ref["ligand_xyz"], config["metrics"]["pocket_cutoff_A"]
        )),
        "primary_pose": "one seed, one sample, model_0; no reference-guided selection",
        "rmsd": "proper Kabsch fit on all 125 matched protein CAs, then minimum ligand heavy-atom RMSD over exact graph isomorphisms; no ligand fit",
        "contact_unit": "unique mapped polymer residue; minimum heavy-atom distance <= cutoff; no waters, additives, or inferred H bonds",
        "protein_atom_coverage": "all reference heavy atoms required except final-residue OXT; additional atoms must match a reference residue-type template and are reported",
        "contact_reference_coverage": "primary contacts retain all predicted protein heavy atoms; a separate diagnostic uses only atoms observed in the reference",
        "success": "all declared thresholds AND geometry screen pass; interpret within the recorded experiment role",
        "confidence": "diagnostic only; no calibrated probability or acceptance threshold",
        "msa": config["msa"],
        "model": model_provenance(),
    })
    artifacts = {
        p.relative_to(run).as_posix(): sha256(p)
        for p in sorted(run.rglob("*")) if p.is_file()
    }
    write_json(run / "manifest.json", {
        "schema_version": 1, "files": artifacts, "analysis_source_sha256_LF": analysis_source_hashes(),
        "source_manifest_sha256": sha256(SOURCE_DIR / "manifest.json"),
        "reference": config["reference"],
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "packages": {p: version(p) for p in ["rdkit", "numpy", "biopython"]},
        "execution_authorized": False,
    })
    return run


def verify_run(run: Path) -> dict:
    run = Path(run).resolve()
    manifest = read_json(run / "manifest.json")
    for rel, digest in manifest["files"].items():
        path = (run / rel).resolve()
        if not path.is_relative_to(run) or sha256(path) != digest:
            raise ValueError(f"Prepared experiment modified: {rel}")
    if analysis_source_hashes() != manifest["analysis_source_sha256_LF"]:
        raise ValueError("Analysis code changed since preregistration; review a new experiment")
    if sha256(SOURCE_DIR / "manifest.json") != manifest["source_manifest_sha256"]:
        raise ValueError("Evidence manifest changed since preparation")
    verify_sources()
    return validate_config(read_json(run / "config.json"))


def score_prediction(config: dict, path: Path) -> dict:
    ref = reference_data(config)
    pred = read_prediction(path, ref)
    selected = [i for i, row in enumerate(ref["mapping"]) if row["fit_CA"]]
    rotation, translation, protein_rmsd = protein_alignment(pred["ca"][selected], ref["ca"][selected])
    pose = symmetry_rmsd(
        ref["mol"], ref["predicted_mol"], ref["ligand_xyz"],
        pred["ligand_xyz"] @ rotation + translation,
    )
    reordered_ligand = pred["ligand_xyz"][pose["reference_to_prediction_atom_indices"]]
    by_residue, _, _ = reference_arrays(ref)
    cutoffs = config["metrics"]["contact_cutoffs_A"]
    contacts = {
        str(c): contact_scores(
            contact_residues(by_residue, ref["ligand_xyz"], c),
            contact_residues(pred["protein_by_residue"], pred["ligand_xyz"], c),
        )
        for c in cutoffs
    }
    matched_contacts = {
        str(c): contact_scores(
            contact_residues(by_residue, ref["ligand_xyz"], c),
            contact_residues(pred["reference_matched_protein_by_residue"], pred["ligand_xyz"], c),
        )
        for c in cutoffs
    }
    geometry = geometry_checks(
        ref["mol"], reordered_ligand, ref["ideal_xyz"],
        pred["protein_xyz"], pred["protein_elements"], config["metrics"]["geometry"],
    )
    pocket = contact_residues(by_residue, ref["ligand_xyz"], config["metrics"]["pocket_cutoff_A"])
    pocket_indices = [i for i, row in enumerate(ref["mapping"]) if row["reference_label"] in pocket]
    primary = contacts[str(config["metrics"]["primary_contact_cutoff_A"])]
    gates = {
        "ligand_rmsd": pose["rmsd_A"] <= config["metrics"]["max_ligand_rmsd_A"],
        "protein_ca_rmsd": protein_rmsd <= config["metrics"]["max_protein_ca_rmsd_A"],
        "contact_precision": primary["precision"] is not None and primary["precision"] >= config["metrics"]["min_contact_precision"],
        "contact_recall": primary["recall"] is not None and primary["recall"] >= config["metrics"]["min_contact_recall"],
        "geometry": geometry["valid"],
    }
    return {
        "experiment_role": experiment_role(config),
        "atom_coverage": pred["atom_coverage"],
        "ligand_pose": pose, "protein_ca_rmsd_A": protein_rmsd,
        "protein_fit_CA_count": len(selected),
        "pocket_ca_rmsd_in_global_frame_A": rmsd(
            pred["ca"][pocket_indices] @ rotation + translation, ref["ca"][pocket_indices]
        ),
        "contacts": contacts, "contacts_reference_observed_atoms": matched_contacts,
        "geometry": geometry,
        "acceptance_gates": gates, "harness_success": all(gates.values()),
        "interpretation": "Single-case pipeline control; not independent evidence of generalisation, affinity, or sample identity.",
    }


def analyse_run(run: Path) -> Path:
    """Only evaluates the sole planned output; does not search for a better pose."""
    from brd4_analysis.phase2.runtime import load_runtime_manifest

    run = Path(run).resolve()
    config = verify_run(run)
    analysis_packages = {p: version(p) for p in ["rdkit", "numpy", "biopython"]}
    if analysis_packages != read_json(run / "manifest.json")["packages"]:
        raise ValueError("Use the frozen Phase 1 analysis environment, separate from the Boltz runtime")
    load_runtime_manifest(run, config)
    folder = run / "outputs" / "boltz_results_4men" / "predictions" / "4men"
    prediction = folder / "4men_model_0.cif"
    if sorted(folder.glob("4men_model_*.cif")) != [prediction]:
        raise ValueError("Expected exactly the single preregistered model_0")
    confidence_path = folder / "confidence_4men_model_0.json"
    confidence = read_json(confidence_path)
    for key in ["confidence_score", "ptm", "iptm", "ligand_iptm", "complex_plddt"]:
        value = confidence.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
            raise ValueError(f"Missing or invalid confidence field: {key}")
    result = score_prediction(config, prediction)
    result.update({
        "confidence_diagnostics": confidence,
        "prediction_sha256": sha256(prediction), "confidence_sha256": sha256(confidence_path),
        "runtime_manifest_sha256": sha256(run / "runtime.json"),
        "prepared_manifest_sha256": sha256(run / "manifest.json"),
        "analysis_packages": analysis_packages,
        "analysis_python": sys.version,
        "analysis_platform": platform.platform(),
    })
    output = run / "analysis"
    output.mkdir(exist_ok=False)
    write_json(output / "metrics.json", result)
    return output / "metrics.json"
