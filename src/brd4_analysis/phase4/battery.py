"""The Phase 4 adversarial battery: fixed ligand/coordinate transformations
scored by the frozen Phase 3 harness.

Every case declares the physical assumption it challenges. The post-audit
battery revision is exploratory and was corrected after earlier outcomes were
reviewed. Cases are grouped
as:

* ``control``: correct or legitimately invariant inputs that MUST pass;
* ``internal_distortion``: chemically impossible coordinate manipulations
  that the geometry screen should catch;
* ``pose_displacement``: ligand-only rigid motions that should move the
  pose RMSD and eventually break contact gates;
* ``symmetry``: automorphism and reflection handling;
* ``protein_side``: protein-only distortions.

Results are computed by ``evaluate()`` — the exact registered gate logic —
plus ``candidate_checks()`` proposals. This is a screen of the harness's
discriminating power on prespecified synthetic cases, not a population-level
estimate of predictive accuracy or failure rate.
"""

from __future__ import annotations

import hashlib
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np

from brd4_analysis.phase2.metrics import graph_mappings
from brd4_analysis.phase4 import perturb
from brd4_analysis.phase4.checks import candidate_checks
from brd4_analysis.phase4.evaluate import (
    aligned_ligand,
    evaluate,
    reference_replay,
    replace_all_coordinates,
    replace_ligand,
    to_prediction_frame,
)

JITTER_SEED = 20261018
BATTERY_REVISION = "post_audit_v1"

ACCEPTANCE_INTERPRETATIONS = {
    "control": "Required control; acceptance is expected and not an adversarial outcome.",
    "internal_distortion": "Acceptance means the registered gates tolerate this fixed internal coordinate change; it is not an empirical model error rate.",
    "pose_displacement": "Acceptance means the registered pose/contact tolerances admit this rigid displacement; it does not by itself establish chemical invalidity.",
    "symmetry": "Acceptance concerns molecular symmetry handling and is not an empirical error rate.",
    "protein_side": "Acceptance means the registered gates tolerate this protein-coordinate distortion; it does not by itself establish ligand chemical invalidity.",
}


def phase4_source_hashes() -> dict[str, str]:
    """LF-normalized hashes of the Phase 4 analysis modules."""
    package = Path(__file__).parent
    return {
        path.name: hashlib.sha256(path.read_text(encoding="utf-8").encode("utf-8")).hexdigest()
        for path in sorted(package.glob("*.py"))
    }


def _ligand_case(fn):
    """Wrap a ligand transform into a pred-dict transform.

    ``fn`` receives the predicted ligand in the REFERENCE frame and in
    REFERENCE atom order, so graph-aware perturbations (which resolve
    targets through ref["mol"]) address the intended atoms. The result is
    converted back to predicted atom order and the predicted frame.
    """
    def apply(ctx):
        pred, ref = ctx["pred"], ctx["ref"]
        mapping = list(graph_mappings(ref["mol"], ref["predicted_mol"])[0])
        aligned = aligned_ligand(pred, ref)          # predicted atom order
        ref_order = aligned[mapping]                 # reference atom order
        moved = fn(ref_order, ctx)
        new_aligned = np.empty_like(aligned)
        new_aligned[mapping] = moved
        return replace_ligand(pred, to_prediction_frame(pred, ref, new_aligned))
    return apply


def _translate_toward_nearest_protein(aligned, ctx, step_A):
    """Push the whole ligand toward its nearest protein atom (aligned frame).

    Direction: ligand centroid -> nearest reference protein atom. This is
    the prespecified steric-clash probe."""
    ref = ctx["ref"]
    centroid = aligned.mean(axis=0)
    protein = np.asarray(
        [a.coord for res in ref["protein_by_residue"].values() for a in res.values()],
        dtype=float,
    )
    nearest = protein[np.argmin(np.linalg.norm(protein - centroid, axis=1))]
    direction = nearest - centroid
    direction /= np.linalg.norm(direction)
    return aligned + step_A * direction


def _automorphism_swap(ref_order, ctx):
    """Permute ligand coordinates by the molecule's non-identity graph
    automorphism (para-substituted tolyl ring mirror). Internal distances
    are exactly preserved and a correct symmetry-aware RMSD must be
    unchanged, not merely small."""
    return ref_order[perturb.automorphism_permutation(ctx["ref"]["mol"])]


# Prespecified 25K distortion targets (CCD atom names).
LINKER_BOND = ("C5", "N")
BICYCLE_BOND = ("N", "C7")
ANGLE_TRIPLE = ("C5", "N", "C7")
REFLECTION_NORMAL = [0.0, 0.0, 1.0]
ANISOTROPIC_AXIS = [0.0, 0.0, 1.0]
TRANSLATE_DIRECTION = [1.0, 0.0, 0.0]
ROTATION_AXIS = [0.0, 0.0, 1.0]

CASES = [
    # --- controls: must pass every gate --------------------------------------
    {"name": "control_real_prediction", "family": "control",
     "assumption": "The unmodified committed Boltz prediction passes all gates "
                   "(baseline pin for the whole battery).",
     "parameters": {},
     "apply": lambda ctx: ctx["pred"]},
    {"name": "control_reference_replay", "family": "control",
     "assumption": "Replaying deposited coordinates must pass: the gates were "
                   "tuned on this structure.",
     "parameters": {},
     "apply": lambda ctx: reference_replay(ctx["ref"])},
    {"name": "control_complex_rigid", "family": "control",
     "assumption": "A rigid motion of the whole complex changes no measurable "
                   "quantity; any failure here is a harness bug.",
     "parameters": {"rotation_axis": [0.4, -0.7, 0.2],
                    "rotation_degrees": 51.0,
                    "translation_A": [9.0, 13.5, -6.0]},
     "apply": lambda ctx: replace_all_coordinates(
         ctx["pred"], perturb.rotation_matrix([0.4, -0.7, 0.2], 51.0),
         np.array([9.0, 13.5, -6.0]))},
    {"name": "control_symmetry_swap", "family": "symmetry",
     "assumption": "The tolyl ring's mirror automorphism permutes equivalent "
                   "atoms; a symmetry-aware RMSD must remain ~0.",
     "parameters": {"automorphism": "non_identity"},
     "apply": _ligand_case(_automorphism_swap)},

    # --- the prior auditor's finding ------------------------------------------
    {"name": "centroid_scale_1.12", "family": "internal_distortion",
     "assumption": "Uniform 12% inflation: chemically implausible for a static "
                   "model (all bond lengths +12%) yet within the registered "
                   "bond, angle, contact and pose tolerances.",
     "parameters": {"factor": 1.12},
     "apply": _ligand_case(lambda a, ctx: perturb.centroid_scale(a, 1.12))},
]

# Parameter sweeps — the same transform at prespecified magnitudes.
for factor in [0.85, 0.90, 1.02, 1.05, 1.08, 1.10, 1.15, 1.20, 1.30]:
    CASES.append({
        "name": f"centroid_scale_{factor:.2f}", "family": "internal_distortion",
        "assumption": "Uniform scaling about the ligand centroid: probes the "
                      "absolute bond-length gate and pose-RMSD gate granularity.",
        "parameters": {"factor": factor},
        "apply": _ligand_case(
            lambda a, ctx, f=factor: perturb.centroid_scale(a, f)),
    })
for step in [0.25, 0.50, 1.00, 1.50, 2.00, 3.00, 4.00]:
    CASES.append({
        "name": f"ligand_translate_{step:.2f}A", "family": "pose_displacement",
        "assumption": "Rigid ligand displacement: internal geometry is unchanged, "
                      "so only pose RMSD and the quantised contact set respond.",
        "parameters": {"translation_vector_A": [step, 0.0, 0.0]},
        "apply": _ligand_case(
            lambda a, ctx, s=step: perturb.translate(
                a, np.asarray(TRANSLATE_DIRECTION) * s)),
    })
for step in [1.0, 2.0, 3.0]:
    CASES.append({
        "name": f"clash_translate_{step:.1f}A", "family": "internal_distortion",
        "assumption": "Translation of the whole ligand toward its nearest "
                      "protein atom: should trigger severe protein-ligand "
                      "overlap once inside van der Waals contact.",
        "parameters": {"step_A": step,
                       "direction": "ligand_centroid_to_nearest_reference_protein_atom"},
        "apply": _ligand_case(
            lambda a, ctx, s=step: _translate_toward_nearest_protein(a, ctx, s)),
    })
for degrees in [15.0, 45.0, 90.0, 180.0]:
    CASES.append({
        "name": f"ligand_rotate_{degrees:.0f}deg", "family": "pose_displacement",
        "assumption": "Rigid rotation of the ligand about its own centroid: "
                      "unchanged internal geometry, large pose error.",
        "parameters": {"axis": ROTATION_AXIS, "degrees": degrees},
        "apply": _ligand_case(
            lambda a, ctx, d=degrees: perturb.rotate_about(
                a, range(len(a)), a.mean(axis=0), ROTATION_AXIS, d)),
    })
for bond_name, bond in [("C5-N", LINKER_BOND), ("N-C7", BICYCLE_BOND)]:
    for degrees in [90.0, 180.0]:
        CASES.append({
            "name": f"torsion_{bond_name}_{degrees:.0f}deg", "family": "internal_distortion",
            "assumption": "Dihedral rotation about bond %s: preserves every "
                          "bond length and valence angle; only nonbonded "
                          "contacts and pose respond." % bond_name,
            "parameters": {"bond": list(bond), "degrees": degrees},
            "apply": _ligand_case(
                lambda a, ctx, b=bond, d=degrees: perturb.torsion_rotate(
                    ctx["ref"]["mol"], a, b[0], b[1], d)),
        })
for delta in [-0.30, 0.20, 0.30]:
    CASES.append({
        "name": f"bond_stretch_C5-N_{delta:+.2f}A", "family": "internal_distortion",
        "assumption": "Single-bond length error injected on C5-N: directly "
                      "probes the 0.25 A absolute bond gate granularity.",
        "parameters": {"bond": ["C5", "N"], "delta_A": delta},
        "apply": _ligand_case(
            lambda a, ctx, d=delta: perturb.bond_stretch(
                ctx["ref"]["mol"], a, "C5", "N", d)),
    })
for degrees in [-35.0, -20.0, 20.0, 35.0]:
    CASES.append({
        "name": f"angle_bend_C5-N-C7_{degrees:+.0f}deg", "family": "internal_distortion",
        "assumption": "Signed valence-angle change at the linker nitrogen: "
                      "probes the 30 deg absolute angle gate granularity.",
        "parameters": {"angle_atoms": list(ANGLE_TRIPLE),
                       "signed_angle_change_deg": degrees},
        "apply": _ligand_case(
            lambda a, ctx, d=degrees: perturb.angle_bend(
                ctx["ref"]["mol"], a, *ANGLE_TRIPLE, d)),
    })
CASES += [
    {"name": "reflect_ligand", "family": "symmetry",
     "assumption": "Mirror image (det=-1) of the achiral ligand: internal "
                   "distances are exactly preserved, so no internal check can "
                   "respond; tests whether pose RMSD or graph chirality catch it.",
     "parameters": {"reflection_normal": REFLECTION_NORMAL},
     "apply": _ligand_case(
         lambda a, ctx: perturb.reflect_plane(a, REFLECTION_NORMAL))},
    {"name": "anisotropic_scale_z_1.25", "family": "internal_distortion",
     "assumption": "25% stretch along one axis: distorts both bond lengths "
                   "and angles nonuniformly.",
     "parameters": {"axis": ANISOTROPIC_AXIS, "factor": 1.25},
     "apply": _ligand_case(
         lambda a, ctx: perturb.anisotropic_scale(a, ANISOTROPIC_AXIS, 1.25))},
    {"name": "collapse_0.70", "family": "internal_distortion",
     "assumption": "Uniform 30% compression: should trip intraligand clash "
                   "and bond-length checks.",
     "parameters": {"factor": 0.70},
     "apply": _ligand_case(lambda a, ctx: perturb.centroid_scale(a, 0.70))},
    {"name": "protein_scale_1.05", "family": "protein_side",
     "assumption": "5% uniform inflation of the predicted protein about its "
                   "own centroid: tests whether the CA-RMSD gate and contact "
                   "distances absorb a systematic protein-frame error.",
     "parameters": {"factor": 1.05},
     "apply": lambda ctx: _scale_protein(ctx, 1.05)},
    {"name": "protein_scale_1.25", "family": "protein_side",
     "assumption": "25% uniform inflation of the predicted protein must exceed "
                   "the registered protein-CA RMSD gate.",
     "parameters": {"factor": 1.25},
     "apply": lambda ctx: _scale_protein(ctx, 1.25)},
]
for sigma in [0.05, 0.20, 0.50]:
    CASES.append({
        "name": f"jitter_sigma_{sigma:.2f}", "family": "internal_distortion",
        "assumption": "Fixed-seed isotropic coordinate noise: uncorrelated "
                      "error, the control counterpart of correlated inflation.",
        "parameters": {"sigma_A": sigma, "seed": JITTER_SEED},
        "apply": _ligand_case(
            lambda a, ctx, s=sigma: perturb.gaussian_jitter(a, s, JITTER_SEED)),
    })


def _scale_protein(ctx, factor):
    pred = dict(ctx["pred"])
    ref_ca = pred["ca"]
    centroid = ref_ca.mean(axis=0)
    pred["ca"] = centroid + factor * (ref_ca - centroid)
    pred["protein_xyz"] = centroid + factor * (pred["protein_xyz"] - centroid)
    pred["protein_by_residue"] = {
        label: centroid + factor * (xyz - centroid)
        for label, xyz in pred["protein_by_residue"].items()
    }
    pred["reference_matched_protein_by_residue"] = {
        label: centroid + factor * (xyz - centroid)
        for label, xyz in pred["reference_matched_protein_by_residue"].items()
    }
    return pred


def run_case(ctx, case: dict) -> dict:
    pred = case["apply"](ctx)
    result = evaluate(ctx["config"], ctx["ref"], pred)
    aligned = aligned_ligand(pred, ctx["ref"])
    # Candidate checks are evaluated over ALL graph automorphisms (best value
    # wins) inside candidate_checks; we pass the ligand in predicted atom
    # order plus the mapping list is reconstructed there. The registered
    # geometry gate keeps the harness's own pose-optimal reordering — that
    # behaviour is part of what is under test.
    result["aligned_pred_order"] = aligned
    result["pred_ligand_xyz"] = pred["ligand_xyz"]
    checks = candidate_checks(ctx, pred, result)
    contacts = {
        c: {
            "precision": s["precision"], "recall": s["recall"],
            "true_positive": s["true_positive"],
            "false_positive": s["false_positive"],
            "false_negative": s["false_negative"],
        }
        for c, s in result["contacts"].items()
    }
    return {
        "name": case["name"],
        "family": case["family"],
        "assumption_challenged": case["assumption"],
        "parameters": case.get("parameters", {}),
        "acceptance_interpretation": ACCEPTANCE_INTERPRETATIONS[case["family"]],
        "gates": result["acceptance_gates"],
        "harness_success": result["harness_success"],
        "metrics": {
            "ligand_pose_rmsd_A": result["ligand_pose"]["rmsd_A"],
            "protein_ca_rmsd_A": result["protein_ca_rmsd_A"],
            "pocket_ca_rmsd_A": result["pocket_ca_rmsd_in_global_frame_A"],
            "contacts": contacts,
            "geometry": {
                "valid": result["geometry"]["valid"],
                "failed_checks": result["geometry"]["failed_checks"],
                "max_bond_length_error_A": result["geometry"]["max_bond_length_error_A"],
                "max_bond_angle_error_deg": result["geometry"]["max_bond_angle_error_deg"],
                "max_aromatic_plane_deviation_A": result["geometry"]["max_aromatic_plane_deviation_A"],
                "min_nonbonded_vdw_ratio": result["geometry"]["min_nonbonded_vdw_ratio"],
                "severe_protein_ligand_overlap_pairs": result["geometry"]["severe_protein_ligand_overlap_pairs"],
                "max_protein_ligand_overlap_A": result["geometry"]["max_protein_ligand_overlap_A"],
            },
        },
        "candidate_checks": checks,
    }


def run_battery(ctx) -> dict:
    """Run every registered case. Deterministic given the pinned inputs."""
    return {
        "schema_version": 2,
        "battery_revision": BATTERY_REVISION,
        "design_status": "post-audit corrected exploratory battery; case outcomes were reviewed during remediation",
        "purpose": "Adversarial screen of the frozen Phase 3 acceptance gates; "
                   "synthetic controls, not population-level accuracy estimates.",
        "provenance": ctx["provenance"],
        "phase4_source_sha256_LF": phase4_source_hashes(),
        "case_parameters": {
            case["name"]: case.get("parameters", {}) for case in CASES
        },
        "random_seeds": {"gaussian_jitter": JITTER_SEED},
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": {
                name: version(name) for name in ("biopython", "numpy", "rdkit")
            },
        },
        "registered_thresholds": ctx["config"]["metrics"],
        "cases": [run_case(ctx, case) for case in CASES],
    }


def summary_table(results: dict) -> list[dict]:
    """Compact per-case verdicts for the findings document."""
    rows = []
    for case in results["cases"]:
        rows.append({
            "name": case["name"],
            "family": case["family"],
            "harness_success": case["harness_success"],
            "failed_gates": [g for g, ok in case["gates"].items() if not ok],
            "failed_candidate_checks": [
                c for c, r in case["candidate_checks"].items() if not r["pass"]
            ],
            "ligand_pose_rmsd_A": case["metrics"]["ligand_pose_rmsd_A"],
            "contact_precision_4A": case["metrics"]["contacts"]["4.0"]["precision"],
            "contact_recall_4A": case["metrics"]["contacts"]["4.0"]["recall"],
        })
    return rows
