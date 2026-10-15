"""Contact-geometry unit tests on synthetic coordinates."""

import numpy as np
import pandas as pd

from brd4_analysis.contacts import distance_matrix, threshold_counts


def test_distance_matrix_known_values():
    a = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]])
    b = np.array([[0.0, 0.0, 5.0]])
    d = distance_matrix(a, b)
    assert d.shape == (2, 1)
    assert d[0, 0] == 5.0
    assert d[1, 0] == np.sqrt(9 + 16 + 25)


def test_distance_matrix_rejects_bad_shape():
    try:
        distance_matrix(np.zeros((2, 2)), np.zeros((2, 3)))
    except ValueError:
        return
    raise AssertionError("expected ValueError")


def test_threshold_counts_boundary_inclusive():
    """A contact at exactly 3.5 Angstrom counts within 3.5A."""
    df = pd.DataFrame(
        {
            "within_3.5A": [True, False, True],
            "within_4.0A": [True, True, True],
            "within_5.0A": [True, True, True],
        }
    )
    counts = threshold_counts(df).set_index("threshold_A")["n_residues"]
    assert counts[3.5] == 2
    assert counts[4.0] == 3
    assert counts[5.0] == 3


def test_threshold_counts_on_real_data(data_dir):
    from brd4_analysis.parse import load_entry
    from brd4_analysis.contacts import contacts_by_residue

    entry = load_entry(data_dir / "raw" / "4MEN.pdb")
    contacts = contacts_by_residue(entry)
    # monotonicity: larger thresholds must capture >= residues
    n35 = contacts["within_3.5A"].sum()
    n40 = contacts["within_4.0A"].sum()
    n50 = contacts["within_5.0A"].sum()
    assert 0 < n35 <= n40 <= n50
    # BRD4 acetyl-lysine pocket: Asn140 must be within 5 A of the ligand
    asn = contacts[contacts["residue"].str.contains("ASN140")]
    assert len(asn) == 1 and bool(asn.iloc[0]["within_5.0A"])
