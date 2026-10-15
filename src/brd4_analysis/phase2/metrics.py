"""Protein-frame pose metrics and deliberately limited heavy-atom geometry checks."""

from __future__ import annotations

from itertools import combinations

import numpy as np
from rdkit import Chem

from brd4_analysis.contacts import BONDI_RADII_A, distance_matrix


def coordinates(values, *, minimum: int = 1) -> np.ndarray:
    points = np.asarray(values, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < minimum:
        raise ValueError(f"Expected at least {minimum} coordinates with shape (N, 3)")
    if not np.isfinite(points).all():
        raise ValueError("Non-finite coordinates")
    return points


def rmsd(a, b) -> float:
    a, b = coordinates(a), coordinates(b)
    if a.shape != b.shape:
        raise ValueError("RMSD requires a complete one-to-one atom correspondence")
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def protein_alignment(moving, fixed) -> tuple[np.ndarray, np.ndarray, float]:
    """Unweighted Kabsch fit, row vectors: aligned = moving @ rotation + translation.

    Reflections are forbidden. Three non-collinear points are the minimum.
    Callers must supply the preregistered residue correspondence, not a best subset.
    """
    moving, fixed = coordinates(moving, minimum=3), coordinates(fixed, minimum=3)
    if moving.shape != fixed.shape:
        raise ValueError("Protein alignment has incomplete residue correspondence")
    cm, cf = moving.mean(axis=0), fixed.mean(axis=0)
    x, y = moving - cm, fixed - cf
    if min(np.linalg.matrix_rank(x), np.linalg.matrix_rank(y)) < 2:
        raise ValueError("Protein alignment is underdetermined (collinear points)")
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(3)
    correction[-1, -1] = 1 if np.linalg.det(u @ vt) > 0 else -1
    rotation = u @ correction @ vt
    translation = cf - cm @ rotation
    return rotation, translation, rmsd(moving @ rotation + translation, fixed)


def graph_mappings(reference: Chem.Mol, predicted: Chem.Mol, cap: int = 4096):
    """Full graph isomorphisms; each tuple maps reference index -> predicted index."""
    for mol in (reference, predicted):
        if mol is None or any(a.GetAtomicNum() <= 1 for a in mol.GetAtoms()):
            raise ValueError("A sanitized heavy-atom molecular graph is required")
    if Chem.MolToSmiles(reference, isomericSmiles=True) != Chem.MolToSmiles(
        predicted, isomericSmiles=True
    ):
        raise ValueError("Ligand graph, charge, isotope or stereochemistry mismatch")
    matches = predicted.GetSubstructMatches(
        reference, uniquify=False, useChirality=True, maxMatches=cap + 1
    )
    if not matches or len(matches) > cap:
        raise ValueError("No complete graph mapping or symmetry enumeration cap exceeded")
    return sorted(set(matches))


def symmetry_rmsd(reference_mol, predicted_mol, reference_xyz, aligned_predicted_xyz):
    """Minimum ligand RMSD over exact graph automorphisms, WITHOUT a ligand fit.

    The prediction must already be in the reference protein frame. Never call
    GetBestRMS/AlignMol here: fitting the ligand would erase pose-placement error.
    """
    reference_xyz = coordinates(reference_xyz)
    predicted_xyz = coordinates(aligned_predicted_xyz)
    if len(reference_xyz) != reference_mol.GetNumAtoms() or len(predicted_xyz) != predicted_mol.GetNumAtoms():
        raise ValueError("Missing or extra ligand coordinates")
    mappings = graph_mappings(reference_mol, predicted_mol)
    scores = [rmsd(reference_xyz, predicted_xyz[list(m)]) for m in mappings]
    best = int(np.argmin(scores))  # mappings sorted: deterministic tie handling
    return {
        "rmsd_A": scores[best],
        "reference_to_prediction_atom_indices": list(mappings[best]),
        "graph_isomorphisms_evaluated": len(mappings),
    }


def contact_residues(protein_by_residue: dict, ligand_xyz, cutoff_A: float) -> set:
    if not np.isfinite(cutoff_A) or cutoff_A <= 0:
        raise ValueError("Contact cutoff must be positive and finite")
    ligand_xyz = coordinates(ligand_xyz)
    return {
        key for key, xyz in protein_by_residue.items()
        if distance_matrix(coordinates(xyz), ligand_xyz).min() <= cutoff_A
    }


def contact_scores(reference: set, predicted: set) -> dict:
    """Residue-level precision/recall. Empty denominators are never perfect scores."""
    common = reference & predicted
    precision = len(common) / len(predicted) if predicted else (0.0 if reference else None)
    recall = len(common) / len(reference) if reference else None
    return {
        "true_positive": len(common),
        "false_positive": len(predicted - reference),
        "false_negative": len(reference - predicted),
        "precision": precision,
        "recall": recall,
        "reference_residues": sorted(reference),
        "predicted_residues": sorted(predicted),
        "scorable": bool(reference),
    }


def geometry_checks(mol, xyz, ideal_xyz, protein_xyz, protein_elements, limits) -> dict:
    """A screen, not a force-field energy or a full PoseBusters validation.

    xyz and ideal_xyz use the SAME molecular graph order. Reference ideal geometry
    is the pinned CCD ideal conformer, not fitted experimental ligand coordinates.
    """
    xyz, ideal = coordinates(xyz), coordinates(ideal_xyz)
    protein = coordinates(protein_xyz)
    if len(xyz) != mol.GetNumAtoms() or xyz.shape != ideal.shape:
        raise ValueError("Geometry check requires complete ligand atom correspondence")
    if len(protein) != len(protein_elements):
        raise ValueError("Protein atom/element correspondence mismatch")
    elements = [a.GetSymbol().upper() for a in mol.GetAtoms()]
    if any(e not in BONDI_RADII_A for e in elements + list(protein_elements)):
        raise ValueError("Unsupported element in the heavy-atom geometry screen")
    # This pilot is deliberately restricted to achiral 25K. Do not pretend this
    # screen validates stereochemistry for arbitrary future chemical series.
    if Chem.FindMolChiralCenters(mol, includeUnassigned=True) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    ):
        raise ValueError("Stereochemical geometry validation is not implemented for this ligand")

    bond_errors = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bond_errors.append(abs(np.linalg.norm(xyz[i] - xyz[j]) - np.linalg.norm(ideal[i] - ideal[j])))

    def angle(points, i, j, k):
        a, b = points[i] - points[j], points[k] - points[j]
        norm = np.linalg.norm(a) * np.linalg.norm(b)
        if norm < 1e-12:
            return None
        return float(np.degrees(np.arccos(np.clip(a @ b / norm, -1, 1))))

    angle_errors = []
    degenerate_angle = False
    for atom in mol.GetAtoms():
        j = atom.GetIdx()
        for i, k in combinations([a.GetIdx() for a in atom.GetNeighbors()], 2):
            actual, expected = angle(xyz, i, j, k), angle(ideal, i, j, k)
            if actual is None or expected is None:
                degenerate_angle = True
            else:
                angle_errors.append(abs(actual - expected))

    ring_deviations = []
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            centered = xyz[list(ring)] - xyz[list(ring)].mean(axis=0)
            _, _, vt = np.linalg.svd(centered)
            ring_deviations.append(float(np.abs(centered @ vt[-1]).max()))

    radii = np.array([BONDI_RADII_A[e] for e in elements])
    topology = Chem.GetDistanceMatrix(mol)
    distances = distance_matrix(xyz, xyz)
    nonbonded = np.triu(topology > 2, k=1)  # exclude bonded and 1-3 neighbours
    ratios = distances[nonbonded] / (radii[:, None] + radii[None, :])[nonbonded]
    overlaps = (
        np.array([BONDI_RADII_A[e] for e in protein_elements])[:, None]
        + radii[None, :] - distance_matrix(protein, xyz)
    )
    max_bond = float(max(bond_errors, default=0))
    max_angle = float(max(angle_errors, default=0))
    max_ring = float(max(ring_deviations, default=0))
    min_ratio = float(ratios.min()) if len(ratios) else None
    severe = int(np.sum(overlaps >= limits["max_protein_ligand_overlap_A"]))
    failures = []
    if max_bond > limits["max_bond_length_error_A"]:
        failures.append("bond_lengths")
    if degenerate_angle or max_angle > limits["max_bond_angle_error_deg"]:
        failures.append("bond_angles")
    if max_ring > limits["max_aromatic_plane_deviation_A"]:
        failures.append("aromatic_planarity")
    if min_ratio is not None and min_ratio < limits["min_nonbonded_vdw_ratio"]:
        failures.append("intraligand_clash")
    if severe:
        failures.append("protein_ligand_clash")
    return {
        "valid": not failures,
        "failed_checks": failures,
        "max_bond_length_error_A": max_bond,
        "max_bond_angle_error_deg": max_angle,
        "degenerate_angle": degenerate_angle,
        "max_aromatic_plane_deviation_A": max_ring,
        "min_nonbonded_vdw_ratio": min_ratio,
        "severe_protein_ligand_overlap_pairs": severe,
        "max_protein_ligand_overlap_A": float(overlaps.max()),
        "stereochemistry": "not_applicable_achiral_graph_without_stereobonds",
        "scope": "heavy_atom_screen_not_full_physical_validation",
    }
