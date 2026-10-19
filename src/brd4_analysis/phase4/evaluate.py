"""Evaluate the registered Phase 3 acceptance gates on in-memory poses.

``evaluate`` reproduces ``phase2.workflow.score_prediction`` exactly, except
that the "prediction" is a coordinate dictionary rather than an mmCIF file.
That lets the adversarial battery interrogate the *scoring* layer without
rewriting structures to disk. ``test_parity_with_registered_scorer`` pins
this replica to the committed ``analysis/metrics.json`` and to
``score_prediction`` on a real prediction file, so drift between the two
implementations would fail the suite.

Coordinate frames: the committed prediction stores coordinates in the
predicted (model) frame. ``aligned_ligand`` maps them into the reference
frame using the same protein-CA Kabsch fit the harness performs, so
perturbations can be specified in interpretable pocket coordinates and
then mapped back with ``to_prediction_frame``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from brd4_analysis.constants import PROJECT_ROOT
from brd4_analysis.phase2.config import read_json, sha256, validate_config
from brd4_analysis.phase2.inputs import read_prediction, reference_data
from brd4_analysis.phase2.metrics import (
    contact_residues,
    contact_scores,
    geometry_checks,
    graph_mappings,
    protein_alignment,
    rmsd,
    symmetry_rmsd,
)
from brd4_analysis.phase2.workflow import reference_arrays
from brd4_analysis.phase4.perturb import require_proper_rotation

PRIMARY_CONFIG = PROJECT_ROOT / "configs" / "phase2" / "4men_primary.json"
PRIMARY_RUN = PROJECT_ROOT / "experiments" / "phase2" / "4men_msa_s20261019_v1"
PREDICTION_CIF = (
    PRIMARY_RUN / "outputs" / "boltz_results_4men" / "predictions" / "4men" / "4men_model_0.cif"
)
ANALYSIS_METRICS = PRIMARY_RUN / "analysis" / "metrics.json"


def load_context(
    config_path: Path = PRIMARY_CONFIG,
    prediction_path: Path = PREDICTION_CIF,
) -> dict:
    """Load the preregistered config, the 4MEN reference and the committed
    Boltz prediction once, so a battery can evaluate many variants."""
    config = validate_config(read_json(config_path))
    ref = reference_data(config)
    pred = read_prediction(prediction_path, ref)
    return {
        "config": config,
        "ref": ref,
        "pred": pred,
        "provenance": {
            "config_sha256": sha256(config_path),
            "prediction_sha256": sha256(prediction_path),
            "analysis_metrics_sha256": sha256(ANALYSIS_METRICS),
        },
    }


def alignment(pred: dict, ref: dict) -> tuple[np.ndarray, np.ndarray, float]:
    """The harness's protein-CA Kabsch fit: predicted frame -> reference."""
    selected = [i for i, row in enumerate(ref["mapping"]) if row["fit_CA"]]
    return protein_alignment(pred["ca"][selected], ref["ca"][selected])


def aligned_ligand(pred: dict, ref: dict) -> np.ndarray:
    """Predicted ligand heavy atoms in the reference protein frame."""
    rotation, translation, _ = alignment(pred, ref)
    return pred["ligand_xyz"] @ rotation + translation


def to_prediction_frame(pred: dict, ref: dict, aligned_xyz) -> np.ndarray:
    """Inverse of :func:`aligned_ligand` (rotation is orthogonal)."""
    rotation, translation, _ = alignment(pred, ref)
    return (np.asarray(aligned_xyz, dtype=float) - translation) @ rotation.T


def aligned_ligand_ref_order(pred: dict, ref: dict) -> np.ndarray:
    """Predicted ligand in the reference frame AND reference atom order,
    using the canonical graph correspondence (first sorted isomorphism)."""
    mapping = list(graph_mappings(ref["mol"], ref["predicted_mol"])[0])
    return aligned_ligand(pred, ref)[mapping]


def replace_ligand(pred: dict, ligand_xyz) -> dict:
    """Copy ``pred`` with a new ligand coordinate block (prediction frame)."""
    new = dict(pred)
    new["ligand_xyz"] = np.asarray(ligand_xyz, dtype=float)
    return new


def replace_all_coordinates(pred: dict, rotation, translation) -> dict:
    """Rigidly move the ENTIRE predicted complex (protein and ligand) by the
    same proper transform — the legitimate invariance control."""
    rotation = require_proper_rotation(rotation)
    translation = np.asarray(translation, dtype=float)
    if translation.shape != (3,) or not np.isfinite(translation).all():
        raise ValueError("Whole-complex translation must be a finite 3-vector")
    new = dict(pred)
    new["ca"] = pred["ca"] @ rotation + translation
    new["ligand_xyz"] = pred["ligand_xyz"] @ rotation + translation
    new["protein_xyz"] = pred["protein_xyz"] @ rotation + translation
    new["protein_by_residue"] = {
        label: xyz @ rotation + translation for label, xyz in pred["protein_by_residue"].items()
    }
    new["reference_matched_protein_by_residue"] = {
        label: xyz @ rotation + translation
        for label, xyz in pred["reference_matched_protein_by_residue"].items()
    }
    return new


def reference_replay(ref: dict) -> dict:
    """A "prediction" that replays the deposited reference coordinates:
    must pass every gate (it is the structure the gates were tuned on).
    The ligand block is expressed in predicted (canonical-rank) atom order,
    matching the contract of ``read_prediction``."""
    by_residue, all_xyz, elements = reference_arrays(ref)
    mapping = graph_mappings(ref["mol"], ref["predicted_mol"])[0]
    ligand = np.empty_like(ref["ligand_xyz"])
    ligand[list(mapping)] = ref["ligand_xyz"]
    return {
        "ca": ref["ca"],
        "protein_by_residue": by_residue,
        "reference_matched_protein_by_residue": by_residue,
        "protein_xyz": all_xyz,
        "protein_elements": elements,
        "ligand_xyz": ligand,
        "atom_coverage": {"note": "synthetic replay of the deposited reference"},
    }


def evaluate(config: dict, ref: dict, pred: dict) -> dict:
    """The five registered gates, identical in logic to ``score_prediction``.

    Returns the same metric field names so results can be compared
    value-for-value against the committed analysis/metrics.json.
    """
    metrics = config["metrics"]
    rotation, translation, protein_rmsd = alignment(pred, ref)
    pose = symmetry_rmsd(
        ref["mol"], ref["predicted_mol"], ref["ligand_xyz"],
        pred["ligand_xyz"] @ rotation + translation,
    )
    reordered_ligand = pred["ligand_xyz"][pose["reference_to_prediction_atom_indices"]]
    by_residue, _, _ = reference_arrays(ref)
    cutoffs = metrics["contact_cutoffs_A"]
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
        pred["protein_xyz"], pred["protein_elements"], metrics["geometry"],
    )
    pocket = contact_residues(by_residue, ref["ligand_xyz"], metrics["pocket_cutoff_A"])
    pocket_indices = [
        i for i, row in enumerate(ref["mapping"]) if row["reference_label"] in pocket
    ]
    primary = contacts[str(metrics["primary_contact_cutoff_A"])]
    gates = {
        "ligand_rmsd": pose["rmsd_A"] <= metrics["max_ligand_rmsd_A"],
        "protein_ca_rmsd": protein_rmsd <= metrics["max_protein_ca_rmsd_A"],
        "contact_precision": primary["precision"] is not None
        and primary["precision"] >= metrics["min_contact_precision"],
        "contact_recall": primary["recall"] is not None
        and primary["recall"] >= metrics["min_contact_recall"],
        "geometry": geometry["valid"],
    }
    return {
        "ligand_pose": pose,
        "protein_ca_rmsd_A": protein_rmsd,
        "protein_fit_CA_count": sum(1 for r in ref["mapping"] if r["fit_CA"]),
        "pocket_ca_rmsd_in_global_frame_A": rmsd(
            pred["ca"][pocket_indices] @ rotation + translation,
            ref["ca"][pocket_indices],
        ),
        "contacts": contacts,
        "contacts_reference_observed_atoms": matched_contacts,
        "geometry": geometry,
        "acceptance_gates": gates,
        "harness_success": all(gates.values()),
    }
