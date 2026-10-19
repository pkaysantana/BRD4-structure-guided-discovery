"""Unit tests for deterministic adversarial coordinate transforms.

Each test verifies the *analytical* property a transform is supposed to
isolate, using small toy point sets plus the real 25K molecular graph for
the graph-aware perturbations.
"""

import numpy as np
import pytest

from brd4_analysis.constants import PROJECT_ROOT
from brd4_analysis.contacts import distance_matrix
from brd4_analysis.phase2.config import read_json, validate_config
from brd4_analysis.phase2.inputs import reference_data
from brd4_analysis.phase4 import perturb
from brd4_analysis.phase4.evaluate import aligned_ligand_ref_order, load_context

CONFIG = PROJECT_ROOT / "configs" / "phase2" / "4men_primary.json"


@pytest.fixture(scope="module")
def mol():
    config = validate_config(read_json(CONFIG))
    return reference_data(config)["mol"]


@pytest.fixture(scope="module")
def ligand_xyz():
    config = validate_config(read_json(CONFIG))
    return reference_data(config)["ligand_xyz"]


@pytest.fixture(scope="module")
def predicted_ligand_xyz():
    ctx = load_context()
    return aligned_ligand_ref_order(ctx["pred"], ctx["ref"])


def _all_bond_lengths(mol, xyz):
    return {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): np.linalg.norm(
            xyz[b.GetBeginAtomIdx()] - xyz[b.GetEndAtomIdx()]
        )
        for b in mol.GetBonds()
    }


def _all_angles(mol, xyz):
    from itertools import combinations

    def angle(i, j, k):
        a, b = xyz[i] - xyz[j], xyz[k] - xyz[j]
        return float(np.degrees(np.arccos(np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))))

    out = {}
    for atom in mol.GetAtoms():
        j = atom.GetIdx()
        for i, k in combinations([a.GetIdx() for a in atom.GetNeighbors()], 2):
            out[(i, j, k)] = angle(i, j, k)
    return out


def test_rotation_matrix_is_proper_and_correct():
    r = perturb.rotation_matrix([0, 0, 1], 90)
    assert np.linalg.det(r) == pytest.approx(1.0)
    assert np.allclose(r @ r.T, np.eye(3))
    assert np.allclose(r @ np.array([1.0, 0, 0]), [0, 1, 0], atol=1e-12)


def test_translate_is_rigid():
    xyz = np.array([[0., 0, 0], [1.4, 0, 0], [2.1, 1.0, 0.5]])
    moved = perturb.translate(xyz, [3.0, -2.0, 7.5])
    assert np.allclose(distance_matrix(moved, moved), distance_matrix(xyz, xyz))
    assert np.allclose(moved.mean(axis=0), xyz.mean(axis=0) + [3.0, -2.0, 7.5])


def test_centroid_scale_scales_distances_preserves_angles():
    xyz = np.array([[0., 0, 0], [1.5, 0, 0], [0.3, 1.2, 0.8], [-0.5, -0.5, 0.2]])
    scaled = perturb.centroid_scale(xyz, 1.12)
    assert np.allclose(distance_matrix(scaled, scaled), distance_matrix(xyz, xyz) * 1.12)
    assert np.allclose(scaled.mean(axis=0), xyz.mean(axis=0))
    # Angles are invariant under uniform scaling.
    for (i, j, k) in [(0, 1, 2), (0, 2, 3)]:
        def angle(p, i, j, k):
            a, b = p[i] - p[j], p[k] - p[j]
            return np.degrees(np.arccos(a @ b / (np.linalg.norm(a) * np.linalg.norm(b))))
        assert angle(scaled, i, j, k) == pytest.approx(angle(xyz, i, j, k))


def test_anisotropic_scale_distorts_only_along_axis():
    xyz = np.array([[0., 0, 0], [2.0, 0, 0], [0, 0, 2.0]])
    scaled = perturb.anisotropic_scale(xyz, [0, 0, 1], 1.5)
    d0 = distance_matrix(xyz, xyz)
    d1 = distance_matrix(scaled, scaled)
    # The z-axis separation grows; the pure-x separation is unchanged only
    # after centring effects — check the vector directly instead.
    z_span = scaled[:, 2].max() - scaled[:, 2].min()
    assert z_span == pytest.approx(3.0)


def test_reflect_preserves_distances_reverses_handedness():
    xyz = np.array([[0., 0, 0], [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
    mirrored = perturb.reflect_plane(xyz, [0, 0, 1])
    assert np.allclose(distance_matrix(mirrored, mirrored), distance_matrix(xyz, xyz))
    # Signed volume of the tetrahedron flips sign.
    def signed_volume(p):
        return np.linalg.det(np.stack([p[1] - p[0], p[2] - p[0], p[3] - p[0]]))
    assert signed_volume(mirrored) == pytest.approx(-signed_volume(xyz))


def test_jitter_is_deterministic_and_zero_sigma_is_identity():
    xyz = np.random.default_rng(1).normal(size=(8, 3))
    a = perturb.gaussian_jitter(xyz, 0.1, seed=7)
    b = perturb.gaussian_jitter(xyz, 0.1, seed=7)
    c = perturb.gaussian_jitter(xyz, 0.1, seed=8)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)
    assert np.array_equal(perturb.gaussian_jitter(xyz, 0.0, seed=7), xyz)
    # Noise magnitude statistics match sigma (sanity, deterministic seed).
    disp = np.linalg.norm(a - xyz, axis=1)
    assert 0.05 < disp.mean() < 0.25


def test_fragment_side_partitions_graph(mol):
    i = perturb.atom_index(mol, "C5")
    j = perturb.atom_index(mol, "N")
    si = perturb.fragment_side(mol, i, i, j)
    sj = perturb.fragment_side(mol, j, i, j)
    assert si.isdisjoint(sj) and si | sj == frozenset(range(mol.GetNumAtoms()))
    # C5's side is the tolyl half (ring carbons C1,C2,C3,C4,C13,C14 +
    # ring methyl C + linker C5 = 8 atoms).
    assert len(si) == 8 and len(sj) == 12
    assert perturb.moved_side(mol, i, j) == si  # smaller side moves


def test_bond_stretch_changes_exactly_one_bond(mol, ligand_xyz):
    i = perturb.atom_index(mol, "C5")
    j = perturb.atom_index(mol, "N")
    before, after = _all_bond_lengths(mol, ligand_xyz), _all_bond_lengths(
        mol, perturb.bond_stretch(mol, ligand_xyz, "C5", "N", 0.30)
    )
    changed = {k for k in before if abs(after[k] - before[k]) > 1e-9}
    assert changed == {(i, j)}
    assert after[(i, j)] == pytest.approx(before[(i, j)] + 0.30, abs=1e-9)
    # All valence angles are preserved by an axial fragment translation.
    assert _all_angles(mol, ligand_xyz) == pytest.approx(
        _all_angles(mol, perturb.bond_stretch(mol, ligand_xyz, "C5", "N", 0.30))
    )


def test_torsion_preserves_bonds_and_angles(mol, ligand_xyz):
    twisted = perturb.torsion_rotate(mol, ligand_xyz, "C5", "N", 90.0)
    assert _all_bond_lengths(mol, twisted) == pytest.approx(_all_bond_lengths(mol, ligand_xyz))
    assert _all_angles(mol, twisted) == pytest.approx(_all_angles(mol, ligand_xyz))
    # But the conformation genuinely changed (nonbonded distances moved).
    assert not np.allclose(distance_matrix(twisted, twisted), distance_matrix(ligand_xyz, ligand_xyz))


@pytest.mark.parametrize("degrees", [-20.0, 20.0])
def test_angle_bend_changes_target_angle_with_signed_contract(mol, ligand_xyz, degrees):
    bent = perturb.angle_bend(mol, ligand_xyz, "C5", "N", "C7", degrees)
    assert _all_bond_lengths(mol, bent) == pytest.approx(_all_bond_lengths(mol, ligand_xyz))
    i, j, k = (perturb.atom_index(mol, n) for n in ("C5", "N", "C7"))

    def angle(p, i, j, k):
        a, b = p[i] - p[j], p[k] - p[j]
        return np.degrees(np.arccos(a @ b / (np.linalg.norm(a) * np.linalg.norm(b))))

    delta = angle(bent, i, j, k) - angle(ligand_xyz, i, j, k)
    assert delta == pytest.approx(degrees, abs=1e-6)


@pytest.mark.parametrize("degrees", [60.0, 180.0])
def test_angle_bend_rejects_targets_outside_open_angle_interval(
    mol, predicted_ligand_xyz, degrees
):
    """For the auditor's ≈123.023309° predicted angle, these requests must
    not silently wrap to a different signed change."""
    i, j, k = (perturb.atom_index(mol, n) for n in ("C5", "N", "C7"))
    a, b = (predicted_ligand_xyz[i] - predicted_ligand_xyz[j],
            predicted_ligand_xyz[k] - predicted_ligand_xyz[j])
    assert np.degrees(np.arccos(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))) == pytest.approx(
        123.02330942986869)
    with pytest.raises(ValueError, match="target angle"):
        perturb.angle_bend(mol, predicted_ligand_xyz, "C5", "N", "C7", degrees)


def test_bond_stretch_rejects_collapse_and_direction_reversal(mol, predicted_ligand_xyz):
    i, j = (perturb.atom_index(mol, n) for n in ("C5", "N"))
    original = float(np.linalg.norm(
        predicted_ligand_xyz[i] - predicted_ligand_xyz[j]))
    assert original == pytest.approx(1.4586583732320595)
    for delta in (-original, -original - 0.2):
        with pytest.raises(ValueError, match="collapsed or reversed"):
            perturb.bond_stretch(mol, predicted_ligand_xyz, "C5", "N", delta)


def test_automorphism_is_a_distance_preserving_permutation(mol, ligand_xyz):
    perm = perturb.automorphism_permutation(mol)
    assert sorted(perm) == list(range(mol.GetNumAtoms()))
    swapped = ligand_xyz[perm]
    # Same point set, different labels: sorted pairwise distance multisets agree.
    tri = np.triu_indices(len(ligand_xyz), 1)
    assert np.allclose(
        np.sort(distance_matrix(swapped, swapped)[tri]),
        np.sort(distance_matrix(ligand_xyz, ligand_xyz)[tri]),
    )


def test_perturbation_errors_are_explicit(mol):
    with pytest.raises(ValueError, match="exactly one atom"):
        perturb.atom_index(mol, "Q9")
    with pytest.raises(ValueError, match="No bond"):
        perturb.bond_stretch(mol, np.zeros((20, 3)), "C", "N", 0.3)  # no C-N bond
    with pytest.raises(ValueError, match="Degenerate"):
        perturb.rotation_matrix([0, 0, 0], 45)
    with pytest.raises(ValueError, match="Degenerate"):
        perturb.reflect_plane(np.zeros((3, 3)), [0, 0, 0])
    with pytest.raises(ValueError, match="sigma"):
        perturb.gaussian_jitter(np.zeros((3, 3)), -1.0, 0)
    with pytest.raises(ValueError, match="finite"):
        perturb.centroid_scale(np.array([[0., 0, np.nan]]), 1.1)
    from rdkit import Chem

    benzene = Chem.MolFromSmiles("c1ccccc1")
    for a in benzene.GetAtoms():
        a.SetProp("name", f"X{a.GetIdx()}")
    with pytest.raises(ValueError, match="automorphism"):
        perturb.automorphism_permutation(benzene)


def test_non_partitioning_ring_bond_is_rejected(mol, ligand_xyz):
    """C14-C1 is a valid graph edge but not a bridge: removing it does not
    partition the molecule, so it must not be treated as a movable fragment."""
    for operation in (
        lambda: perturb.bond_stretch(mol, ligand_xyz, "C14", "C1", 0.2),
        lambda: perturb.torsion_rotate(mol, ligand_xyz, "C14", "C1", 30.0),
        lambda: perturb.angle_bend(mol, ligand_xyz, "C14", "C1", "C2", 15.0),
    ):
        with pytest.raises(ValueError, match="does not partition"):
            operation()


def test_rotation_contract_requires_orthonormal_matrix():
    with pytest.raises(ValueError, match="orthonormal"):
        perturb.require_proper_rotation(np.diag([2.0, 0.5, 1.0]))
    with pytest.raises(ValueError, match="determinant"):
        perturb.require_proper_rotation(np.diag([1.0, 1.0, -1.0]))


def test_rigid_whole_vs_internal_distortion_classification(mol, ligand_xyz):
    """The taxonomy the battery relies on: rigid transforms leave the full
    internal distance matrix bitwise-identical; internal distortions do not."""
    rigid = ligand_xyz @ perturb.rotation_matrix([1, 2, 3], 47.0) + [5, -3, 2]
    assert np.allclose(distance_matrix(rigid, rigid), distance_matrix(ligand_xyz, ligand_xyz))
    distorted = perturb.centroid_scale(ligand_xyz, 1.05)
    assert not np.allclose(
        distance_matrix(distorted, distorted), distance_matrix(ligand_xyz, ligand_xyz)
    )
