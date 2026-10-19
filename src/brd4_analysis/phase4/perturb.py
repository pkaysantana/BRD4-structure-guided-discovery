"""Deterministic adversarial ligand-coordinate perturbations.

Each transform maps an (N, 3) heavy-atom coordinate array to a new array in
the same atom order and coordinate frame and never mutates its input.
Distortion targets are addressed by CCD atom name resolved through the
molecular graph, so every case is prespecified and reproducible.

Physical framing used by the battery:

* a rigid-body motion of the whole complex changes no internal distance and
  is legitimately invisible to internal-geometry checks;
* an internal distortion changes bonds, angles, torsions or nonbonded
  distances and *should* trip the geometry screen;
* a ligand-only rigid motion is a pose displacement: internal geometry is
  untouched but the ligand's placement relative to the pocket changes.
"""

from __future__ import annotations

from collections import deque

import numpy as np

from brd4_analysis.phase2.metrics import coordinates


def atom_index(mol, name: str) -> int:
    hits = [a.GetIdx() for a in mol.GetAtoms() if a.GetProp("name") == name]
    if len(hits) != 1:
        raise ValueError(f"Expected exactly one atom named {name!r}")
    return hits[0]


def bond_atoms(mol, name_a: str, name_b: str) -> tuple[int, int]:
    i, j = atom_index(mol, name_a), atom_index(mol, name_b)
    if mol.GetBondBetweenAtoms(i, j) is None:
        raise ValueError(f"No bond {name_a}-{name_b}")
    return i, j


def fragment_side(mol, start: int, cut_a: int, cut_b: int) -> frozenset:
    """Atoms reachable from ``start`` once bond (cut_a, cut_b) is removed."""
    blocked = frozenset({cut_a, cut_b})
    seen, queue = {start}, deque([start])
    while queue:
        a = queue.popleft()
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            b = nb.GetIdx()
            if b in seen or {a, b} == blocked:
                continue
            seen.add(b)
            queue.append(b)
    return frozenset(seen)


def _vector3(value, name: str) -> np.ndarray:
    vector = np.asarray(value, dtype=float)
    if vector.shape != (3,) or not np.isfinite(vector).all():
        raise ValueError(f"{name} must be a finite 3-vector")
    return vector


def _finite_scalar(value, name: str) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _cut_partition(mol, i: int, j: int) -> tuple[frozenset, frozenset]:
    """The two fragments produced by cutting a bridge bond.

    Ring bonds and any other non-partitioning cut are invalid: the reachable
    sets overlap, so rotating or translating one "fragment" would move the
    whole molecule and silently fail to apply the requested distortion.
    """
    si, sj = fragment_side(mol, i, i, j), fragment_side(mol, j, i, j)
    if not si.isdisjoint(sj) or si | sj != frozenset(range(mol.GetNumAtoms())):
        raise ValueError("Selected bond does not partition the molecular graph")
    return si, sj


def moved_side(mol, i: int, j: int) -> frozenset:
    """Which fragment a single-bond perturbation displaces: the smaller side
    of a graph-partitioning cut (minimal structural disruption); ties break
    to the side containing the lower minimum atom index."""
    si, sj = _cut_partition(mol, i, j)
    if len(si) != len(sj):
        return si if len(si) < len(sj) else sj
    return si if min(si) < min(sj) else sj


def rotation_matrix(axis, degrees: float) -> np.ndarray:
    """Proper (det=+1) Rodrigues rotation acting on column vectors."""
    axis = _vector3(axis, "Rotation axis")
    degrees = _finite_scalar(degrees, "Rotation angle")
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        raise ValueError("Degenerate rotation axis")
    x, y, z = axis / norm
    c, s = np.cos(np.radians(degrees)), np.sin(np.radians(degrees))
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
    ])


def require_proper_rotation(matrix, *, tolerance: float = 1e-9) -> np.ndarray:
    """Validate an orthonormal rotation matrix with determinant +1."""
    rotation = np.asarray(matrix, dtype=float)
    if rotation.shape != (3, 3) or not np.isfinite(rotation).all():
        raise ValueError("Rotation must be a finite 3x3 matrix")
    if np.linalg.norm(rotation.T @ rotation - np.eye(3)) > tolerance:
        raise ValueError("Rotation matrix is not orthonormal")
    if abs(np.linalg.det(rotation) - 1.0) > tolerance:
        raise ValueError("Rotation matrix must have determinant +1")
    return rotation


def translate(xyz, vector) -> np.ndarray:
    """Rigid translation: no internal distance changes."""
    return coordinates(xyz) + _vector3(vector, "Translation")


def rotate_about(xyz, indices, origin, axis, degrees: float) -> np.ndarray:
    """Proper rotation of ``indices`` about the axis through ``origin``."""
    xyz = coordinates(xyz)
    idx_list = list(indices)
    if not idx_list or not all(isinstance(i, (int, np.integer)) for i in idx_list):
        raise ValueError("Rotation indices must be non-empty integers")
    idx = np.array(sorted(set(idx_list)), dtype=int)
    if len(idx) != len(idx_list) or idx[0] < 0 or idx[-1] >= len(xyz):
        raise ValueError("Rotation indices are duplicated or out of range")
    rotation = rotation_matrix(axis, degrees)
    origin = _vector3(origin, "Rotation origin")
    out = xyz.copy()
    out[idx] = (xyz[idx] - origin) @ rotation.T + origin
    return out


def centroid_scale(xyz, factor: float, centroid=None) -> np.ndarray:
    """Uniform expansion/compression about a point.

    Preserves every angle and molecular handedness and scales all pairwise
    distances by ``factor``. Exact coplanarity is preserved; a nonzero
    absolute plane-deviation measure scales with the molecule.
    """
    xyz = coordinates(xyz)
    factor = _finite_scalar(factor, "Scale factor")
    if factor <= 0:
        raise ValueError("Scale factor must be positive")
    centre = xyz.mean(axis=0) if centroid is None else _vector3(centroid, "Scale centroid")
    return centre + factor * (xyz - centre)


def anisotropic_scale(xyz, axis, factor: float, centroid=None) -> np.ndarray:
    """Scale coordinates along one axis only: distorts lengths and angles."""
    xyz = coordinates(xyz)
    axis = _vector3(axis, "Scaling axis")
    factor = _finite_scalar(factor, "Scale factor")
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        raise ValueError("Degenerate scaling axis")
    n = axis / norm
    centre = xyz.mean(axis=0) if centroid is None else _vector3(centroid, "Scale centroid")
    d = xyz - centre
    return centre + d + (factor - 1.0) * np.outer(d @ n, n)


def reflect_plane(xyz, normal, point=None) -> np.ndarray:
    """Householder mirror (det=-1). Preserves every internal distance and
    angle exactly, so it is invisible to all internal-geometry checks; only
    the protein-frame pose RMSD can respond, and for an achiral ligand even
    that may stay small."""
    xyz = coordinates(xyz)
    normal = _vector3(normal, "Reflection normal")
    norm = np.linalg.norm(normal)
    if norm < 1e-12:
        raise ValueError("Degenerate reflection normal")
    n = normal / norm
    centre = xyz.mean(axis=0) if point is None else _vector3(point, "Reflection point")
    householder = np.eye(3) - 2.0 * np.outer(n, n)
    return centre + (xyz - centre) @ householder


def gaussian_jitter(xyz, sigma: float, seed: int) -> np.ndarray:
    """Fixed-seed isotropic noise; deterministic per (sigma, seed, order)."""
    xyz = coordinates(xyz)
    sigma = _finite_scalar(sigma, "sigma")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    rng = np.random.default_rng(int(seed))
    return xyz + rng.normal(0.0, sigma, xyz.shape)


def bond_stretch(mol, xyz, name_i: str, name_j: str, delta_A: float) -> np.ndarray:
    """Translate the smaller fragment beyond a bridge bond along its axis.

    The i-j covalent length changes by ``delta_A`` and endpoint valence
    angles are preserved, but many nonbonded distances between the two
    fragments necessarily change as well.
    """
    xyz = coordinates(xyz)
    delta_A = _finite_scalar(delta_A, "Bond delta")
    i, j = bond_atoms(mol, name_i, name_j)
    side = moved_side(mol, i, j)
    moved_end, fixed_end = (i, j) if i in side else (j, i)
    direction = xyz[moved_end] - xyz[fixed_end]
    norm = np.linalg.norm(direction)
    if norm < 1e-12:
        raise ValueError("Degenerate bond geometry")
    target_length = norm + delta_A
    if target_length <= 1e-12:
        raise ValueError("Bond stretch would produce a collapsed or reversed bond")
    out = xyz.copy()
    out[list(side)] += delta_A * direction / norm
    return out


def torsion_rotate(mol, xyz, name_i: str, name_j: str, degrees: float) -> np.ndarray:
    """Rotate the smaller fragment beyond a bridge bond about the bond axis.
    Preserves every bond length and every valence angle exactly; only
    dihedral (1-4 and longer) distances change."""
    xyz = coordinates(xyz)
    i, j = bond_atoms(mol, name_i, name_j)
    side = moved_side(mol, i, j)
    on_axis = i if i in side else j
    axis = xyz[j] - xyz[i]
    if np.linalg.norm(axis) < 1e-12:
        raise ValueError("Degenerate bond geometry")
    return rotate_about(xyz, side, xyz[on_axis], axis, degrees)


def angle_bend(mol, xyz, name_i: str, name_j: str, name_k: str, degrees: float) -> np.ndarray:
    """Change the signed valence angle i-j-k by ``degrees``.

    The fragment on the i side of bridge bond i-j is rotated about the axis
    perpendicular to the i-j-k plane through j. Positive values open the
    angle; negative values close it. Every bond length is preserved. The
    cut must partition the graph, the input angle must be nondegenerate,
    and the requested target angle must remain inside the open interval
    (0, 180) degrees.
    """
    xyz = coordinates(xyz)
    degrees = _finite_scalar(degrees, "Angle delta")
    i, j, k = atom_index(mol, name_i), atom_index(mol, name_j), atom_index(mol, name_k)
    if mol.GetBondBetweenAtoms(i, j) is None or mol.GetBondBetweenAtoms(j, k) is None:
        raise ValueError(f"{name_i}-{name_j}-{name_k} is not a valence angle")
    ui = xyz[i] - xyz[j]
    uk = xyz[k] - xyz[j]
    if np.linalg.norm(ui) < 1e-12 or np.linalg.norm(uk) < 1e-12:
        raise ValueError("Degenerate angle geometry")
    axis = np.cross(ui, uk)
    if np.linalg.norm(axis) < 1e-12:
        raise ValueError("Degenerate angle plane")
    current = float(np.degrees(np.arctan2(np.linalg.norm(axis), float(ui @ uk))))
    target = current + degrees
    if not 0.0 < target < 180.0:
        raise ValueError(
            f"Requested angle change would produce target angle {target:.12g} degrees"
        )
    side, other = _cut_partition(mol, i, j)
    if k not in other:
        raise ValueError("Angle arm k must remain on the fixed side")
    return rotate_about(xyz, side, xyz[j], axis, -degrees)


def automorphism_permutation(mol) -> list[int]:
    """The non-identity graph automorphism as a permutation of atom indices.
    Requires exactly one non-identity automorphism (true for 25K)."""
    matches = sorted(set(mol.GetSubstructMatches(mol, uniquify=False, useChirality=True)))
    nontrivial = [m for m in matches if list(m) != list(range(mol.GetNumAtoms()))]
    if len(nontrivial) != 1:
        raise ValueError("Expected exactly one non-identity automorphism")
    return list(nontrivial[0])
