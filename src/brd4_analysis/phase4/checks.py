"""Candidate geometry checks, proposed but NOT registered.

Each check here is evaluated alongside the five frozen gates in the
adversarial battery. Thresholds are exploratory proposals: they were
chosen so that the *unmodified* committed prediction passes and the
fixed exploratory adversarial distortions fail. None of these values is a
registered acceptance criterion, and adopting any of them would require
a separate preregistered evaluation on more than one case.

Two reference frames matter:

* checks against the pinned CCD *ideal* conformer are conformation-blind:
  the supplied 25K coordinates are a specific three-dimensional reference
  geometry, not the deposited bound conformation. Bond lengths and valence
  angles remain informative because they are conformation-invariant.
* checks against the *deposited* ligand coordinates measure internal
  agreement with the experimental structure; they conflate pose and
  conformation fidelity, which is acceptable for reconstruction
  validation but must be declared.
"""

from __future__ import annotations

import numpy as np

from brd4_analysis.contacts import distance_matrix
from brd4_analysis.phase2.metrics import coordinates, graph_mappings


def bond_error_stats(mol, xyz, reference_xyz) -> dict:
    """Per-bond signed and relative length errors against a reference
    geometry (same atom order). Signed bias is the key observable: a
    uniform scaling inflates *every* bond in the same direction, a
    correlated pattern that random coordinate perturbations are less likely
    to reproduce exactly. The aggregate bias is diagnostic, not definitive;
    independent noise can also change bond lengths and its signed mean."""
    xyz, ref = coordinates(xyz), coordinates(reference_xyz)
    signed, relative = [], []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ideal = np.linalg.norm(ref[i] - ref[j])
        err = np.linalg.norm(xyz[i] - xyz[j]) - ideal
        signed.append(float(err))
        relative.append(float(abs(err) / ideal))
    return {
        "mean_signed_A": float(np.mean(signed)),
        "max_abs_A": float(np.max(np.abs(signed))),
        "max_relative": float(np.max(relative)),
    }


def similarity_scale(xyz, target_xyz) -> float:
    """Optimal uniform scale s mapping ``xyz`` onto ``target_xyz`` under a
    proper rotation (Umeyama scale factor, same atom order). Values below 1
    mean ``xyz`` is larger than the target and must be shrunk to fit."""
    x, y = coordinates(xyz), coordinates(target_xyz)
    xc, yc = x - x.mean(axis=0), y - y.mean(axis=0)
    u, _, vt = np.linalg.svd(xc.T @ yc)
    correction = np.eye(3)
    correction[-1, -1] = 1 if np.linalg.det(u @ vt) > 0 else -1
    fitted = xc @ (u @ correction @ vt)
    denom = float(np.sum(fitted * fitted))
    if denom < 1e-12:
        raise ValueError("Degenerate point cloud")
    return float(np.sum(fitted * yc) / denom)


def internal_distance_rmsd(xyz, target_xyz) -> float:
    """RMSD between the full pairwise heavy-atom distance matrices of two
    conformers (same atom order). Sensitive to scaling, conformational
    change and local distortion; insensitive to pose."""
    x, y = coordinates(xyz), coordinates(target_xyz)
    tri = np.triu_indices(len(x), 1)
    diff = distance_matrix(x, x)[tri] - distance_matrix(y, y)[tri]
    return float(np.sqrt(np.mean(diff * diff)))


def contact_distance_profile(by_residue, ref_ligand_xyz, pred_by_residue,
                             pred_ligand_xyz, residues) -> np.ndarray:
    """Per-residue (predicted minus reference) minimum heavy-atom contact
    distance for the given residue set. Detects ligand shifts that leave
    the contact *set* unchanged while systematically changing distances."""
    deltas = []
    for label in residues:
        d0 = distance_matrix(coordinates(by_residue[label]), ref_ligand_xyz).min()
        d1 = distance_matrix(coordinates(pred_by_residue[label]), pred_ligand_xyz).min()
        deltas.append(float(d1 - d0))
    return np.asarray(deltas)


# --- Candidate gate registry -------------------------------------------------
# Each entry: a statistic, a proposed limit, and a comparison direction.
# "le" means pass when value <= limit; "abs_le" means pass when |value| <= limit.

def _le(value, limit):
    return bool(value <= limit)


def _abs_le(value, limit):
    return bool(abs(value) <= limit)


def candidate_checks(ctx, pred: dict, evaluated: dict) -> dict:
    """Run every proposed check for one battery case.

    ``ctx`` is the load_context() bundle; ``pred`` is the evaluated variant's
    prediction dict (its protein may differ from ctx["pred"]); ``evaluated``
    is the dict returned by evaluate(), augmented by the caller with
    ``aligned_pred_order`` (predicted ligand in the reference frame, still in
    predicted atom order) and ``pred_ligand_xyz``.

    For a symmetric ligand neither the canonical nor the pose-optimal atom
    correspondence is universally correct: internal-geometry statistics are
    therefore reported as the *best* (smallest absolute) value over all
    exact graph automorphisms — "agreement with the reference up to
    symmetry". The best-over-mappings convention cannot mask a real
    internal distortion, because no automorphism rescues stretched bonds or
    a scaled distance matrix.
    """
    mol = ctx["ref"]["mol"]
    ideal = ctx["ref"]["ideal_xyz"]
    deposited = ctx["ref"]["ligand_xyz"]
    aligned = evaluated["aligned_pred_order"]
    mappings = graph_mappings(mol, ctx["ref"]["predicted_mol"])

    def best(stat_fn):
        return min(stat_fn(aligned[list(m)]) for m in mappings)

    ideal_bias = best(lambda lig: abs(bond_error_stats(mol, lig, ideal)["mean_signed_A"]))
    ideal_rel = best(lambda lig: bond_error_stats(mol, lig, ideal)["max_relative"])
    dep_bias = best(lambda lig: abs(bond_error_stats(mol, lig, deposited)["mean_signed_A"]))
    log_scale_dep = best(lambda lig: abs(np.log(similarity_scale(lig, deposited))))
    dist_rmsd_dep = best(lambda lig: internal_distance_rmsd(lig, deposited))

    primary_cutoff = ctx["config"]["metrics"]["primary_contact_cutoff_A"]
    primary_residues = evaluated["contacts"][str(primary_cutoff)]["reference_residues"]
    by_residue = {
        label: np.asarray([a.coord for a in atoms.values()], dtype=float)
        for label, atoms in ctx["ref"]["protein_by_residue"].items()
    }
    deltas = contact_distance_profile(
        by_residue, deposited,
        pred["protein_by_residue"], evaluated["pred_ligand_xyz"],
        primary_residues,
    )

    stats = {
        "abs_bond_signed_bias_vs_ccd_A": ideal_bias,
        "max_bond_relative_error_vs_ccd": ideal_rel,
        "abs_bond_signed_bias_vs_deposited_A": dep_bias,
        "abs_log_similarity_scale_vs_deposited": log_scale_dep,
        "internal_distance_rmsd_vs_deposited_A": dist_rmsd_dep,
        "contact_distance_mean_shift_A": float(np.mean(deltas)),
        "contact_distance_max_abs_shift_A": float(np.max(np.abs(deltas))),
        "ligand_pose_rmsd_A": evaluated["ligand_pose"]["rmsd_A"],
        "pocket_ca_rmsd_A": evaluated["pocket_ca_rmsd_in_global_frame_A"],
    }
    proposals = [
        ("bond_signed_bias", "abs_bond_signed_bias_vs_ccd_A", 0.04, _le,
         "reference-free when the predicted graph has a trusted CCD template; "
         "otherwise retrospective",
         "correlated bond-length drift; uniform inflation/deflation signature. "
         "Deposited 25K itself carries +0.030 A bias vs CCD ideal, so the "
         "exploratory limit must exceed that; catches ~4%+ uniform inflation"),
        ("bond_relative_error", "max_bond_relative_error_vs_ccd", 0.10, _le,
         "reference-free when the predicted graph has a trusted CCD template; "
         "otherwise retrospective",
         "per-bond relative tolerance. Deposited max is 8.6% vs CCD ideal, "
         "so the exploratory limit sits just above it; catches ~9%+ uniform inflation"),
        ("similarity_scale", "abs_log_similarity_scale_vs_deposited", 0.03, _le,
         "reference-dependent retrospective reconstruction check",
         "|ln s| of the optimal scale onto the deposited ligand; direct scaling detector"),
        ("internal_distance_rmsd", "internal_distance_rmsd_vs_deposited_A", 0.15, _le,
         "reference-dependent retrospective reconstruction check",
         "full pairwise-distance agreement with the deposited conformer; can reject "
         "a legitimate alternate conformer"),
        ("contact_distance_bias", "contact_distance_mean_shift_A", 0.10, _abs_le,
         "reference-dependent retrospective pose/contact check",
         "systematic shift of minimum contact distances within an unchanged contact set"),
        ("pose_rmsd_tightened", "ligand_pose_rmsd_A", 0.5, _le,
         "reference-dependent retrospective reconstruction check",
         "registered metric with a tighter (non-registered) threshold"),
        ("pocket_ca_rmsd", "pocket_ca_rmsd_A", 1.0, _le,
         "reference-dependent retrospective protein-fit check",
         "local fit quality of pocket CAs in the global frame"),
    ]
    out = {}
    for name, key, limit, cmp_, applicability, rationale in proposals:
        value = float(stats[key])
        out[name] = {
            "value": value,
            "statistic": key,
            "proposed_limit": limit,
            "pass": cmp_(value, limit),
            "applicability": applicability,
            "rationale": rationale,
        }
    return out
